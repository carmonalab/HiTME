# HiT class definition ##############################################################################

## Build S4 objects to store data
methods::setClass(Class = "HiT",
                  slots = list(
                    metadata = "data.frame",
                    predictions = "list",
                    composition = "list",
                    aggregated_profile = "list",
                    version = "list"
                  )
)



# ProjecTILs on multiple reference maps #############################################################

ProjecTILs.classifier.multi <- function(object,
                                        ref.maps,
                                        bparam = NULL,
                                        layer1_link = "CellOntology_ID") {

  # count how many cells are found by scGate for each map reference
  if (layer1_link %in% names(object@meta.data)) {
    filter.cells <- F
    map.celltypes <- lapply(ref.maps,
                            function(x) {
                              ct <- x@misc$layer1_link
                              nrow(object@meta.data[object@meta.data[[layer1_link]] %in% ct,])
                            })

    present <- names(ref.maps[map.celltypes > 1])
    no.present <- names(ref.maps)[!names(ref.maps) %in% present]

    if (length(present) == 0) {
      message(paste("No cells linked to reference maps found by", layer1_link,
                    ".\nNot running Projectils"))
      run <- FALSE
    } else {
      run <- TRUE
      message("Performing Projectils classification for ",
              paste(present, collapse = ", "))

      if (length(ref.maps) != length(present)) {
        message("Not doing mapping for ",
                paste(no.present, collapse = ", ") )
        ref.maps <- ref.maps[present]
      }
    }

  } else {
    filter.cells <- TRUE

    message("Not scGate classification found in this object\n",
            "Running default scGate (or layer1 classification) filtering within Projectils)")
    run <- TRUE
  }

  if (run) {
    suppressWarnings(require(ProjecTILs))

    functional.clusters <-
      BiocParallel::bplapply(
        X = names(ref.maps),
        BPPARAM = bparam,
        function(m) {
          if (filter.cells) {
            # don't subset if no scGate gating is found in metadata, instead run scGate within projectils
            subset.object <- object
          } else {
            map.celltype <- ref.maps[[m]]@misc$layer1_link
            subset.object <- object[,object@meta.data[[layer1_link]] %in% map.celltype]

          }

          if (ncol(subset.object)>0) {
            message("\nRunning Projectils for ", m, " reference")
            # it is mandatory to make run in serial (ncores = 1 and BPPARAM SerialParam)
            subset.object.pred <- ProjecTILs:::classifier.singleobject(subset.object,
                                                                       ref = ref.maps[[m]],
                                                                       filter.cells = filter.cells,
                                                                       ncores = 1)

            return(subset.object.pred)

          } else {
            message("Not Running Projectils classifier for reference map",
                    paste(map.celltype, collapse = ", "), ". No cells found")
          }
        }
      )


    if (!all(unlist(lapply(functional.clusters, function(x) {is.null(x)})))) {
      functional.clusters <- data.table::rbindlist(lapply(functional.clusters,
                                                          data.table::setDT,
                                                          keep.rownames = TRUE)) %>%
        #remove duplicated row.names (classified by different ref.maps)
        dplyr::filter(!duplicated(rn)) %>%
        tibble::column_to_rownames("rn")

      # remove already present functional.cluster columns
      rm <- grep("^functional.cluster",
                 names(object@meta.data))
      if (length(rm) > 0) {
        object@meta.data <- object@meta.data[, -rm, drop = F]
      }
      object@meta.data <- merge(object@meta.data,
                                functional.clusters,
                                by = 0,
                                all.x = T) %>%
        tibble::column_to_rownames("Row.names")

      # save names of reference maps run (will be lost if objects remerged)
      object@misc[["layer2_param"]][["functional.cluster"]][["References_executed"]] <- names(ref.maps)

    } else {
      object@meta.data <- object@meta.data %>%
        dplyr::mutate(functional.cluster = NA,
                      functional.cluster.conf = NA)
    }
  }
  return(object)
}




match_dictionary <- function(cell_type,
                             dictionary = NULL) {

  if (is.null(dictionary)) {
    dict <- c("^T|CTL" = "T cell",
              "^B|B$" = "B cell",
              "^NK" = "NK",
              "Mono" = "Monocyte-like",
              "^Fibro" = "Fibroblast",
              "^Malig|cancer|tumo[u]r" = "Malignant",
              "Myeloid" = "Myeloid",
              "HSC" = "HSC",
              "DC|Dendri" = "DC",
              "^Epi|^EC" = "Epithelial",
              "Mast" = "Mastocyte",
              "Macroph" = "Macrophage",
              "^Endo" = "Endothelial",
              "Neutro" = "Neutrophil",
              "Strom" = "Stromal")
  }

  for (keyword in names(dict)) {
    if (grepl(keyword, cell_type, ignore.case = TRUE)) {
      cell_type <- dict[[keyword]]
      break
    }
  }
  return(cell_type)
}



StandardizeCellNames <- function(cell.names, dictionary = NULL) {
  standarized_cellnames <- sapply(cell.names, match_dictionary, dictionary)
  standarized_cellnames <- unname(standarized_cellnames)

  return(standarized_cellnames)
}



set_parallel_params <- function(ncores,
                                bparam,
                                progressbar)
{
  if (is.null(ncores)) {
    ncores <- 1
  }

  if (ncores > parallelly::availableCores()) {
    ncores <- parallelly::availableCores()
    message("Using more cores available in this computer, reducing number of cores to ", ncores)
  }

  # set parallelization parameters
  if (is.null(bparam)) {
    if (ncores > 1) {
      param <- BiocParallel::MulticoreParam(workers =  ncores,
                                            progressbar = progressbar)
    } else {
      param <- BiocParallel::SerialParam()
    }
  } else {
    param <- bparam
  }
  return(param)
}





compositional_data <- function(data,
                               split.by = NULL,
                               group.by.1 = NULL,
                               useNA = FALSE,
                               clr_zero_impute_perc = 1,
                               only.counts = FALSE) {

  if (all(is.na(data[[group.by.1]]))) {
    if (!only.counts) {
      ctable <- data.frame("celltype" = character(),
                           "cell_counts" = integer(),
                           "freq" = numeric(),
                           "clr" = numeric())
    } else {
      ctable <- data.frame("celltype" = character(),
                           "cell_counts" = integer())
    }
    return(ctable)
  } else {
    # set grouping variables
    gr_vars <- c(split.by, group.by.1)
    gr_vars2 <- c(split.by)

    ctable <- data %>%
      # drop = F keeps all levels of the factor
      dplyr::group_by(dplyr::across(dplyr::all_of(gr_vars)), .drop = F) %>%
      dplyr::summarize(cell_counts = dplyr::n()) %>%
      dplyr::ungroup()

    colnames(ctable)[1] <- "celltype"

    if (!only.counts) {
      ctable <- ctable %>%
        dplyr::filter(if (!useNA) !is.na(.data[["celltype"]])
                      else rep(TRUE, n())) %>%
        dplyr::group_by(across(all_of(gr_vars2))) %>%
        dplyr::mutate(freq = cell_counts/sum(cell_counts) * 100,
                      !!"celltype" := coalesce(.data[["celltype"]], "NA")) %>%
        as.data.frame() %>%
        na.omit()

      if (nrow(ctable) > 0) {
        # compute clr
        clr.df <- ctable %>%
          dplyr::select(-cell_counts) %>%
          # add pseudocount
          dplyr::mutate(freq = freq + clr_zero_impute_perc) %>%
          tidyr::pivot_wider(names_from = "celltype",
                             values_from = "freq")

        # accommodate df for clr transformation
        ## Remove character columns
        num_cols_bool_idx <- sapply(clr.df, is.numeric)
        num_cols <- names(clr.df)[num_cols_bool_idx]
        chr_cols <- names(clr.df)[!num_cols_bool_idx]
        clr.df.ref <- clr.df %>% dplyr::select(all_of(num_cols))

        clr <- Hotelling::clr(clr.df.ref)

        # add extra cols (if any)
        clr <- cbind(clr.df[,chr_cols],clr)  %>%
          tidyr::pivot_longer(-chr_cols,
                              names_to = "celltype",
                              values_to = "clr")
        # join clr df to main dataframe
        ctable <- dplyr::left_join(ctable,
                                   clr,
                                   by = c(chr_cols, "celltype"))
      }
    }

    return(ctable)
  }
}


# get.cluster.score helper functions ##############################################################################

### Pre-process pseudobulk count data

preproc_pseudobulk <- function (matrix,
                                metadata,
                                cluster.by,
                                nVarGenes = 500,
                                gene.filter = "HVG",
                                black.list = NULL) {

  suppressMessages({
    suppressWarnings({
      matrix <- DESeq2.normalize(matrix = matrix,
                                 metadata = metadata,
                                 cluster.by = cluster.by,
                                 nVarGenes = nVarGenes,
                                 gene.filter = gene.filter,
                                 black.list = black.list)
    })
  })

  return(matrix)
}



### Just a helper for preproc_pseudobulk, additional pre-processing steps are needed

DESeq2.normalize <- function(matrix,
                             metadata,
                             cluster.by,
                             nVarGenes = 500,
                             gene.filter = "HVG",
                             black.list = NULL) {

  # Get black list
  if (is.null(black.list)) {
    utils::data("default_black_list")
  }
  black.list <- unlist(black.list)

  # Normalize pseudobulk data using DESeq2
  # do formula for design with the cluster.by elements in order
  dformula <-  stats::formula(paste("~", cluster.by))
  data <- DESeq2::DESeqDataSetFromMatrix(countData = matrix,
                                         colData = metadata,
                                         design = dformula)
  data <- DESeq2::estimateSizeFactors(data)

  nsub <- min(1000, sum(rowMeans(BiocGenerics::counts(data, normalized=TRUE)) > 10 ))

  # transform counts using vst
  data <- DESeq2::vst(data, blind = T, nsub = nsub)
  data <- SummarizedExperiment::assay(data)

  # Remove black listed genes from the matrix
  data <- data[!row.names(data) %in% black.list,]

  # filter genes accordingly
  if (gene.filter == "HVG") {
    # get top variable genes
    rv <- MatrixGenerics::rowVars(data)
    select <- order(rv, decreasing=TRUE)[seq_len(min(nVarGenes, length(rv)))]
    select <- row.names(data)[select]

  } else if (gene.filter == "default_filter") {
    # get predetermined list of genes
    utils::data("GO_accession_default")
    select <- GO_default
  }

  data <- data[select[select %in% row.names(data)],]

  return(data)
}




get.scores <- function(matrix,
                       cluster_labels,
                       scores,
                       modularity.k,
                       dist.method = "euclidean",
                       ntests = 100, # number of shuffling events
                       seed = 22, # seed for random shuffling
                       title = "", # Title for summary
                       # For PCA
                       invisible = c("var", "quali")) {

  matrix <- t(matrix)

  results <- list()

  # Check if there are at least 2 clusters,
  # that there are more than 4 samples and
  # that there are more samples than clusters (not each sample is one separate cluster)
  if (length(unique(cluster_labels)) > 1 &
      length(cluster_labels) > 4 &
      nrow(matrix) > length(unique(cluster_labels))) {

    results[["Feature_matrix"]] <- matrix
    results[["Distance_matrix"]] <- stats::dist(matrix, method = dist.method)


    # Plot PCA ###############################################
    results[["Plots"]][["PCA_feature_space"]] <- plot_PCA(matrix,
                                                          color.cluster.by = cluster_labels,
                                                          label = "var",
                                                          invisible = invisible) +
      ggplot2::ggtitle("PCA on feature matrix")

    results[["Plots"]][["PCA_sample_space"]] <- plot_PCA(as.matrix(results[["Distance_matrix"]]),
                                                         scale. = TRUE,
                                                         color.cluster.by = cluster_labels,
                                                         invisible = c("var", "quali")) +
      ggplot2::ggtitle("PCA on sample distance matrix")

    #
    #     # Clustering ###############################################
    #
    #     ## Find number of clusters with silhouette method ###############################################
    #     n_clust <- factoextra::fviz_nbclust(x = mat,
    #                                         FUNcluster = cluster::pam,
    #                                         method = "silhouette",
    #                                         k.max = 10,
    #                                         print.summary = TRUE) + theme_minimal() + ggtitle("Optimal number of clusters
    #                                                                                            - silhouette method")
    #     print(n_clust)
    #     n_clust <- n_clust$data
    #     max_cluster <- as.numeric(n_clust$clusters[which.max(n_clust$y)])
    #
    #     ## Plot PCA with clustering ###############################################
    #     results[["Plots"]][["PCA_feature_space"]] <- plot_PCA(matrix,
    #                                                           color.cluster.by = cluster_labels,
    #                                                           label = "var",
    #                                                           invisible = invisible) +
    #       ggplot2::ggtitle("PCA on feature matrix")
    #
    #     results[["Plots"]][["PCA_sample_space"]] <- plot_PCA(as.matrix(results[["Distance_matrix"]]),
    #                                                          scale. = TRUE,
    #                                                          color.cluster.by = cluster_labels,
    #                                                          invisible = c("var", "quali")) +
    #       ggplot2::ggtitle("PCA on sample distance matrix")


    # Calculate scores + plots ###############################################

    suppressMessages({
      suppressWarnings({

        for (s in scores) {

          ## Silhouette_isolated (new) ###############################################
          if (s == "Silhouette_isolated") {
            sils <- calc_sil_onelabel(labels = cluster_labels,
                                      dist = results[["Distance_matrix"]],
                                      return_mean_for_permtest = FALSE)

            avg_per_group <- sils %>%
              dplyr::group_by(group) %>%
              dplyr::summarize(avg_sil_width = mean(sil_width))

            avg_sil <- mean(sils[["sil_width"]])

            CI_intervals <- stats::t.test(sils[["sil_width"]])$conf.int

            p_val <- perm_test(fun = calc_sil_onelabel,
                               data = results[["Distance_matrix"]],
                               labels = cluster_labels,
                               obs = avg_sil,
                               ntests = ntests,
                               seed = seed)


            results[["Scores"]][[s]] <- list("samples" = sils,
                                             "avg_per_group" = avg_per_group,
                                             "summary" = avg_sil,
                                             "conf_int" = CI_intervals,
                                             "n" = nrow(sils),
                                             "p_value" = p_val)

            p <- plot_silhouette(sil_scores = results[["Scores"]][[s]],
                                 title = "Silhouette (isolated) plot")

            results[["Plots"]][[s]] <- p
          }

          ## Silhouette (original) ###############################################
          if (s == "Silhouette") {
            sils <- calc_sil(labels = cluster_labels,
                             dist = results[["Distance_matrix"]],
                             return_mean_for_permtest = FALSE)

            avg_per_group <- sils %>%
              dplyr::group_by(group) %>%
              dplyr::summarize(avg_sil_width = mean(sil_width))

            avg_sil <- mean(sils[["sil_width"]])

            CI_intervals <- t.test(sils[["sil_width"]])$conf.int

            p_val <- perm_test(fun = calc_sil,
                               data = results[["Distance_matrix"]],
                               labels = cluster_labels,
                               obs = avg_sil,
                               ntests = ntests,
                               seed = seed)

            results[["Scores"]][[s]] <- list("samples" = sils,
                                             "avg_per_group" = avg_per_group,
                                             "summary" = avg_sil,
                                             "conf_int" = CI_intervals,
                                             "n" = nrow(sils),
                                             "p_value" = p_val)

            p <- plot_silhouette(sil_scores = results[["Scores"]][[s]],
                                 title = "Silhouette plot")

            results[["Plots"]][[s]] <- p
          }

          ## Modularity ###############################################
          if (s == "Modularity") {

            if (length(cluster_labels) >= (modularity.k+1)) {
              g <- scran::buildKNNGraph(matrix,
                                        transposed = TRUE,
                                        k = modularity.k)

              # Calculate modularity score
              modularity_score <- igraph::modularity(g, membership = as.numeric(factor(cluster_labels)))

              # Plotting the graph
              g <- igraph::set_vertex_attr(g, "name", value = cluster_labels)

              p_val <- perm_test(fun = calc_modularity,
                                 data = matrix,
                                 labels = cluster_labels,
                                 obs = modularity_score,
                                 ntests = ntests,
                                 seed = seed)

              results[["Scores"]][[s]] <- list("igraph" = g,
                                               "summary" = modularity_score,
                                               "n" = length(g),
                                               "p_value" = p_val)

              # Using ggraph layout = "kk" instead of default "stress",
              # as "kk" handles disconnected communities much better (showing them separately, instead of fur balling them like "stress")
              p <- ggraph::ggraph(g, layout = 'kk') +
                ggraph::geom_edge_link(color = "grey", edge_width = 0.2) +
                ggraph::geom_node_point(ggplot2::aes(fill = as.factor(names(igraph::V(g)))),
                                        shape = 21,
                                        color = "black",
                                        size = 5) +
                ggplot2::ggtitle(paste("KNN plot with k = ", modularity.k,
                                       "\nModularity score = ", round(modularity_score, 3),
                                       ifelse(!is.null(p_val),
                                              paste("\np-value:",
                                                    format.pval(p_val, digits = 3)), ""))) +
                ggplot2::labs(fill = "Groups") +
                ggplot2::theme(panel.background = element_rect(fill = "white"))

              results[["Plots"]][[s]] <- p
            }
          }
        }
      })
    })

    # Plot dendrogram  (TODO NEEDS REWORK) ###############################################
    # pc2 <- pc$x[,1:ndim] %>%
    #   as.data.frame() %>%
    #   mutate(celltype = metadata[[df.score[x, 1]]])
    #
    # # Calculate the row-wise average grouping by row names
    # pc2 <- aggregate(. ~ celltype, data = pc2, FUN = mean) %>%
    #   tibble::column_to_rownames("celltype") %>%
    #   as.matrix()
    #
    # dist_group <- stats::dist(pc2,
    #                           method = df.score[x, 2])
    # hclust <- stats::hclust(dist_group,
    #                         method = hclust.method)
    # dendo <- ggdendro::ggdendrogram(as.dendrogram(hclust)) +
    #   ggtitle(paste0("Hierarchical clustering dendrogram - ",
    #                  hclust.method))


    # Combine plots ###############################################
    results[["Plots"]][["Summary_plot"]] <- patchwork::wrap_plots(results[["Plots"]],
                                                                  ncol = 2) +
      patchwork::plot_layout(widths = 1) +
      patchwork::plot_annotation(title,
                                 theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                                                           hjust = 0,
                                                                                           size = 20)))
    return(results)

  } else {

    return(NULL)
  }
}





## get.scores helpers ##############################################################################

### Calculate silhouette score of group vs all others (instead of nearest other group)

calc_sil_onelabel <- function(labels,
                              dist,
                              return_mean_for_permtest = TRUE) {

  sils <- lapply(
    unique(labels),
    function(a) {

      clus1 <- ifelse(labels == a, 1, 2)

      sil <- cluster::silhouette(clus1, dist)

      # silhouette score for each sample
      # change names back to character
      sil_clus1 <- as.data.frame(sil) %>%
        dplyr::filter(cluster == 1) %>%
        dplyr::mutate(cluster = a) %>%
        dplyr::arrange(dplyr::desc(sil_width))


      size <- nrow(sil_clus1)

      sil.sumA <- data.frame(group = a,
                             size = size,
                             avg_sil_width = mean(sil_clus1$sil_width))

      return(list("samples" = sil_clus1, # silhouette score for each sample
                  "avg_per_group" = sil.sumA))
    }
  )

  # join results
  sils <- lapply(sils, function(x) {x[["samples"]]}) %>%
    data.table::rbindlist() %>%
    dplyr::rename(group = cluster) %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(dplyr::desc(sil_width), .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rowid = dplyr::row_number())

  if (return_mean_for_permtest) {
    return(mean(sils[["sil_width"]]))
  } else {
    return(sils)
  }
}


### Calculate (classical) silhouette score

calc_sil <- function(labels,
                     dist,
                     return_mean_for_permtest = TRUE) {

  sils <- cluster::silhouette(x = as.numeric(factor(labels)),
                              dist = dist) %>%
    as.data.frame() %>%
    dplyr::rename(group = cluster) %>%
    dplyr::mutate(group = labels) %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(dplyr::desc(sil_width), .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rowid = dplyr::row_number())

  if (return_mean_for_permtest) {
    return(mean(sils[["sil_width"]]))
  } else {
    return(sils)
  }
}



### Calculate modularity score

calc_modularity <- function(labels,
                            matrix,
                            k = 3) {

  g <- scran::buildKNNGraph(matrix,
                            transposed = TRUE,
                            k = k)

  # Calculate modularity score
  modularity_score <- igraph::modularity(g, membership = as.numeric(factor(labels)))

  return(modularity_score)
}


### Permutation test to calculate p-value
# Permutation testing allows for calculating the p-value, e.g. for the silhouette score
# to calculate how likely it is to observe a specific labelling
# For permutation testing labels are randomly shuffled to calculate the random distribution of labels as the null hypothesis

perm_test <- function (fun,
                       data,
                       labels,
                       obs,
                       ntests,
                       seed) {

  if (!(length(ntests) == 1 &&
        is.numeric(ntests) &&
        ntests >= 0)) {
    stop("Please provide a number for ntests permutation testing")
  }

  if (!(length(ntests) == 1 &&
        is.numeric(seed))) {
    stop("Please provide a number for setting seed")
  }

  # perform the shuffling
  if (ntests > 0) {

    random_values <- c()

    for (u in 1:ntests) {

      # seeding for reproducibility
      set.seed(seed + u)
      random_labels <- sample(labels) %>%
        as.factor() %>%
        as.numeric()

      rand_val <- fun(random_labels, data) %>%
        mean()

      random_values <- c(random_values, rand_val)
    }

    # compute p-value
    p_val <- p.val_zscore(obs = obs,
                          random_values = random_values)

    return(p_val)

  } else if (ntests == 0) {
    return(NULL)
  }
}



### Compute p-value based on z score

p.val_zscore <- function(obs,
                         random_values) {

  mean <- mean(random_values)
  sd <- sd(random_values)
  z_score <- (obs - mean) / sd
  p_val <- 1 - stats::pnorm(z_score)
  return(p_val)
}




### get.scores plot helpers

plot_PCA <- function(matrix,
                     scale. = FALSE,
                     label = "all",
                     invisible = c("var", "quali"),
                     geom.var = c("arrow", "text"),
                     col.var = "steelblue",
                     alpha.var = 0.3,
                     repel = TRUE,
                     color.cluster.by = "none") {

  # Remove constant columns with variance = 0
  constant_columns <- apply(matrix, 2, function(col) var(col) == 0)
  constant_columns_indices <- which(constant_columns)
  if (length(constant_columns_indices) > 0) {
    matrix <- matrix[, -constant_columns_indices, drop=FALSE]
  }

  # Otherwise, stats::prcomp would fail with error

  # If there are still columns left
  if (ncol(matrix) > 0) {
    res.pca <- stats::prcomp(matrix, scale. = scale.)
    suppressWarnings(
      suppressMessages(
        p <- factoextra::fviz_pca(res.pca,
                                  habillage = color.cluster.by,
                                  label = label,
                                  pointsize = 3,
                                  invisible = invisible,
                                  geom.var = geom.var,
                                  col.var = col.var,
                                  alpha.var = alpha.var,
                                  repel = repel) +
          # do not rescale x and y-axes, so scale does not get distorted and represents actual distance better
          coord_equal() +
          # Remove shapes added by fviz_pca
          scale_shape_manual(values = c(rep(19, length(unique(color.cluster.by)))))
      )
    )

    return(p)
  }
}


plot_silhouette <- function(sil_scores,
                            title = "title") {

  xend <- length(sil_scores[["samples"]][["sil_width"]])
  ci <- sil_scores[["conf_int"]]
  m <- sil_scores[["summary"]]

  p <- ggplot2::ggplot(sil_scores[["samples"]],
                       ggplot2::aes(x = rowid,
                                    y = sil_width,
                                    fill = group)) +
    ggplot2::geom_bar(position = "dodge", stat = "identity") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::ggtitle(paste0(title,
                            "\nAverage silhouette width = ", round(m, 3),
                            "   95% CI: [", round(ci[1], 3), ", ", round(ci[2], 3), "]",
                            ifelse(!is.null(sil_scores[["p_value"]]),
                                   paste("\np-value: ",
                                         format.pval(sil_scores[["p_value"]],
                                                     digits = 3)), ""))) +
    ggplot2::labs(fill = "Groups") +
    ggplot2::geom_ribbon(aes(x = 1:xend,
                             ymin = ci[1],
                             ymax = ci[2]),
                         alpha = 0.15,
                         inherit.aes = FALSE) +
    ggplot2::annotate("segment",
                      x = 1,
                      xend = xend,
                      y = ci[1],
                      lty = 1,
                      alpha = 0.15) +
    ggplot2::annotate("segment",
                      x = 1,
                      xend = xend,
                      y = ci[2],
                      lty = 1,
                      alpha = 0.15) +
    ggplot2::annotate("segment",
                      x = 1,
                      xend = xend,
                      y = m,
                      lty = 2) +
    ggplot2::theme_bw()

  return(p)
}



### Combine multiple p-values from get.cluster.scores if batching

combine_pvals <- function(list_pvals_n,
                          pval.combine.method) {

  use_fallback_method <- FALSE

  # To prevent error from metap, replace zeros with some small value
  list_pvals_n[["p_value"]][list_pvals_n[["p_value"]] == 0] <- 1e-300

  if (!"p_value" %in% names(list_pvals_n)) {
    return(list_pvals_n)
  } else
    if (pval.combine.method == "weighted_zmethod") {
      if (!length(list_pvals_n[["p_value"]]) == length(list_pvals_n[["n"]])) {
        use_fallback_method <- TRUE
      } else {
        list_pvals_n[["p_value"]] <- metap::sumz(p = list_pvals_n[["p_value"]],
                                                 weights = list_pvals_n[["n"]])[["p"]][1,1]
      }
    } else
      if (pval.combine.method == "fisher" |
          use_fallback_method) {
        list_pvals_n[["p_value"]] <- metap::sumlog(p = list_pvals_n[["p_value"]])[["p"]]
      } else {
        stop("pval.combine.method not recognized.
         Please see documentation for possible methods.")
      }

  list_pvals_n[["n"]] <- sum(list_pvals_n[["n"]])

  return(list_pvals_n)
}
