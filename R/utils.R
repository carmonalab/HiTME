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
      message(paste("No cells linked to reference maps found by", layer1_link, ".\nNot running Projectils"))
      run <- FALSE
    } else {
      run <- TRUE
      message("Performing Projectils classification for ", paste(present, collapse = ", "))

      if (length(ref.maps) != length(present)) {
        message("Not doing mapping for ", paste(no.present, collapse = ", ") )
        ref.maps <- ref.maps[present]
      }
    }

  } else {
    filter.cells <- T

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
        #remove duplicated rownames (classified by different ref.maps)
        dplyr::filter(!duplicated(rn)) %>%
        tibble::column_to_rownames("rn")

      # remove already present functional.cluster columns
      rm <- grep("^functional.cluster", names(object@meta.data))
      if (length(rm) > 0) {
        object@meta.data <- object@meta.data[, -rm, drop = F]
      }
      object@meta.data <- merge(object@meta.data, functional.clusters, by = 0, all.x = T) %>%
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



# Match words with grep ##############################################################################

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


# sapply match_dictionary to all elements in vector ####################################################

StandardizeCellNames <- function(cell.names, dictionary = NULL) {

  standarized_cellnames <- sapply(cell.names,
                                  match_dictionary,
                                  dictionary)
  standarized_cellnames <- unname(standarized_cellnames)

  return(standarized_cellnames)
}



# Calculate silhouette score ##############################################################################

silhoutte_onelabel <- function(labels = NULL, # vector of labels
                               dist = NULL, # distance object
                               ntests = 0, # number of shuffling events
                               seed = 22, # seed for random suffling
                               ncores = parallelly::availableCores() - 2,
                               bparam = NULL,
                               progressbar = TRUE) {

  if (is.null(labels) || !is.vector(labels)) {
    stop("Please provide a vector of the labels")
  }

  if (is.null(dist) || !class(dist) == "dist") {
    stop("Please provide a a dissimilarity object inheriting from class dist or coercible to one")
  }

  if (!is.numeric(ntests)) {
    stop("Please provide a number of shuffling test for the random assignation of the label")
  }

  if (!is.numeric(seed)) {
    stop("Please provide a number for setting seed")
  }

  # set parallelization parameters
  param <- set_parallel_params(ncores = ncores,
                               bparam = bparam,
                               progressbar = progressbar)

  tlabels <- unique(labels)

  len <- length(labels)

  sils <- BiocParallel::bplapply(
    X = tlabels,
    BPPARAM = param,
    function(a) {

      x_one <- ifelse(labels == a, 2, 1) %>%
        as.factor() %>% as.numeric()

      silh <- cluster::silhouette(x_one, dist)

      # silhouette score for each sample
      # change names back to character
      sil.res.cell <- as.data.frame(silh) %>%
        dplyr::filter(cluster == 2) %>%
        dplyr::mutate(cluster = a) %>%
        arrange(desc(sil_width))


      size <- nrow(sil.res.cell)

      sil.sumA <- data.frame(cluster = a,
                             iteration = "NO",
                             size = size,
                             avg_sil_width = mean(sil.res.cell$sil_width))

      # dataframe for results of each iteration of shuffling
      bots.df <- data.frame(matrix(nrow = 0, ncol = 4))
      names(bots.df) <- c("cluster", "iteration", "size", "avg_sil_width")

      if (ntests > 0) {
        # perform the shuffling
        for (u in 1:ntests) {
          # random vector
          vec <- rep(1, length(labels))
          # seeding for reproducibility
          seed <- seed + which(tlabels == a)
          set.seed(seed)
          random_sample <- sample(1:len, size)
          vec[random_sample] <- 2
          vec <- vec %>% as.factor() %>% as.numeric()

          # run silhouette
          silh <- cluster::silhouette(vec, dist)
          sil.res <- as.data.frame(silh) %>%
            dplyr::filter(cluster == 2)

          sil.sumB <- data.frame(cluster = a,
                                 iteration = u,
                                 size = size,
                                 avg_sil_width = mean(sil.res$sil_width))

          bots.df <- rbind(bots.df, sil.sumB)

        }

        # compute p-value
        sil.sumA$p_val <- p.val_zscore(obs = sil.sumA$avg_sil_width,
                                       random.values = bots.df$avg_sil_width)
        sil.sumA$p_val_adj <- p.adjust(sil.sumA$p_val,
                                       method = "fdr",
                                       n = length(tlabels))
        sil.sumA$conf.int.95 <- NA
        # Compute confidence interval
        bots.df <- bots.df %>%
          mutate(conf.int.95 = I(list(quantile(avg_sil_width, c(0.025,0.975)))))
        bots.df$p_val <- NA
        bots.df$p_val_adj <- NA

        sil.sumA <- rbind(sil.sumA, bots.df)

      }

      return(list("cell" = sil.res.cell, # silhouette score for each sample
                  "summary" = sil.sumA))
    }
  )

  # join results
  sil.all <- lapply(sils, function(x) {x[["cell"]]}) %>%
    data.table::rbindlist()

  sil.sum <- lapply(sils, function(x) {x[["summary"]]}) %>%
    data.table::rbindlist()


  # order silhoute width per cell type
  sil.all <- sil.all %>%
    dplyr::mutate(rowid = dplyr::row_number())


  ret.list <- list("cell" = sil.all,
                   "summary" = sil.sum)

  return(ret.list)
}



# Compute p-value based on z score ##############################################################################

p.val_zscore <- function(obs = NULL,
                         random.values = NULL) {

  mean <- mean(random.values)
  sd <- sd(random.values)
  z_score <- (obs - mean) / sd
  p_value <- 1-pnorm(z_score)
  return(p_value)
}



# Plot shuffling interval confidence ######################################################

plot.score <- function(df = NULL,
                       type = "density") {

  if (type == "density") {
    spl <- split(df, df$cluster)
    pl.list <- list()
    for (a in names(spl)) {

      it <- spl[[a]] %>% filter(iteration != "NO") %>%
        mutate(adj_p_val = NA)
      obs <- spl[[a]] %>% filter(iteration == "NO")
      xit <- data.frame(vline_x = unlist(it[1,"conf.int.95"]))
      xobs <- data.frame(vline_x = obs$avg_sil_width)
      pl.list[[a]] <- it %>%
        ggplot(aes(x=avg_sil_width)) +
        geom_density(fill = "grey") +
        geom_vline(data = xit,
                   aes(xintercept = vline_x,
                       color = "confidence_interval95"),
                   linetype = 2) +
        geom_vline(data = xobs,
                   aes(xintercept = vline_x,
                       color = "Observed_value"),
                   linetype = 2) +
        geom_text(x = mean(c(xit$vline_x[2], xobs$vline_x)),
                  y = max(density(it$avg_sil_width)$y)*0.7,
                  label = paste("adj_pval:",round(obs$adj_p_val,2)),
                  vjust = 1, hjust = 1.5) +
        theme_bw() +
        scale_color_manual(name = "linetypes" ,
                           values = c(Observed_value = "blue",
                                      confidence_interval95 = "red")) +
        ggtitle(paste(a, "- Number of samples:", obs$size))

    }
    pl <- ggpubr::ggarrange(plotlist = pl.list, common.legend = T)
  }

  if (type == "barplot") {

    df0 <- df %>%
      filter(iteration == "NO") %>%
      mutate(adj_p_val = ifelse(is.nan(adj_p_val), 1, adj_p_val),
             Sign = ifelse(adj_p_val < 0.05, "Sig", "Not_Sig"))

    pl <- df0 %>%
      ggplot(aes(cluster, avg_sil_width,
                 fill = Sign)) +
      geom_col() +
      theme_bw()
  }
  return(pl)
}






# Set the parallelization parameters using Biocparallel #######################################################

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
    if (ncores>1) {
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



# Plot PCA ##############################################################################

plot_PCA <- function(data = NULL,
                     invisible = c("var", "quali"),
                     color.cluster.by = "none") {
  res.pca <- prcomp(t(data))
  p <- fviz_pca(res.pca,
                habillage = color.cluster.by,
                label = "var",
                pointsize = 3,
                invisible = c("var", "quali"))
  return(p)
}



# Compute compositional data with dplyr ###################################################################

compositional_data <- function(data,
                               split.by = NULL,
                               group.by.1 = NULL,
                               useNA = FALSE,
                               clr_zero_impute_perc = 1,
                               only.counts = FALSE
) {

  # set grouping variables
  gr_vars <- c(split.by, group.by.1)
  gr_vars2 <- c(split.by)


  ctable <- data %>%
    # drop = F keeps all levels of the factor
    dplyr::group_by(across(all_of(gr_vars)), .drop = F) %>%
    dplyr::summarize(cell_counts = dplyr::n()) %>%
    dplyr::ungroup()

  if (!only.counts) {
    ctable <- ctable %>%
      dplyr::filter(if (!useNA) !is.na(.data[[group.by.1]])
                    else rep(TRUE, n())) %>%
      dplyr::group_by(across(all_of(gr_vars2))) %>%
      dplyr::mutate(freq = cell_counts/sum(cell_counts) * 100,
                    !!group.by.1 := coalesce(.data[[group.by.1]], "NA")) %>%
      as.data.frame()

    # compute clr
    clr.df <- ctable %>%
      dplyr::select(-cell_counts) %>%
      # add pseudocount
      dplyr::mutate(freq = freq + clr_zero_impute_perc) %>%
      tidyr::pivot_wider(names_from = group.by.1,
                         values_from = "freq")

    # accomodate df for clr transformation
    ## Remove character columns
    num_cols_bool_idx <- sapply(clr.df, is.numeric)
    num_cols <- names(clr.df)[num_cols_bool_idx]
    chr_cols <- names(clr.df)[!num_cols_bool_idx]
    clr.df.ref <- clr.df %>% dplyr::select(all_of(num_cols))

    clr <- Hotelling::clr(clr.df.ref)

    # add extra cols (if any)
    clr <- cbind(clr.df[,chr_cols],clr)  %>%
      tidyr::pivot_longer(-chr_cols, names_to = group.by.1, values_to = "clr")
    # join clr df to main dataframe
    ctable <- dplyr::left_join(ctable, clr, by = c(chr_cols, group.by.1))
  }

  return(ctable)
}


# Normalize pseudobulk data using DESeq2 ##############################################################################

DESeq2.normalize <- function(matrix,
                             metadata,
                             cluster.by,
                             black.list
) {
  # do formula for design with the cluster.by elements in order
  dformula <- formula(paste("~", paste(cluster.by, collapse =  " + ")))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = matrix,
                                        colData = metadata,
                                        design = dformula)
  dds <- DESeq2::estimateSizeFactors(dds)

  nsub <- min(1000,
              sum(rowMeans(BiocGenerics::counts(dds, normalized=TRUE)) > 5 ))
  # transform counts usign vst
  vsd <- DESeq2::vst(dds, blind = T, nsub = nsub)
  vsd <- SummarizedExperiment::assay(vsd)

  # Remove black listed genes from the matrix
  vsd <- vsd[!rownames(vsd) %in% black.list,]

  return(vsd)

}
