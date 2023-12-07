#############################################################################################################
# HiT class definition
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


#############################################################################################################
# Projectils on multiple reference maps
ProjecTILs.classifier.multi <- function(object,
                                        ref.maps,
                                        bparam = NULL,
                                        layer1_link = "CellOntology_ID"){

  # count how many cells are found by scGate for each map reference
  if(layer1_link %in% names(object@meta.data)){
    filter.cells <- F
    map.celltypes <- lapply(ref.maps,
                            function(x){
                              ct <- x@misc$layer1_link
                              # collpase if multiple cell types, such as MoMac
                              if(length(ct)>1){ct <- paste(ct, collapse = "_")}

                              nrow(object@meta.data[object@meta.data[[layer1_link]] %in% ct,])
                            })
    present <- names(ref.maps[map.celltypes > 0])
    no.present <- names(ref.maps)[!names(ref.maps) %in% present]

    if(length(present) == 0){
      message(paste("No cells linked to reference maps found by", layer1_link, ".\nNot running Projectils"))
      run <- FALSE
    } else {
      run <- TRUE
      message("Performing Projectils classification for ", paste(present, collapse = ", "))

      if(length(ref.maps) != length(present)){
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

  if(run){
  suppressWarnings(require(ProjecTILs))

  functional.clusters <-
    BiocParallel::bplapply(
      X = names(ref.maps),
      BPPARAM = bparam,
      FUN = function(m){
        if(filter.cells){
          # don't subset if no scGate gating is found in metadata, instead run scGate within projectils
          subset.object <- object
        } else {
        map.celltype <- ref.maps[[m]]@misc$layer1_link
        subset.object <- object[,object@meta.data[[layer1_link]] %in% map.celltype]

        }

        if(ncol(subset.object)>0){
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


  if(!all(unlist(lapply(functional.clusters, function(x){is.null(x)})))){
    functional.clusters <- data.table::rbindlist(lapply(functional.clusters,
                                                        data.table::setDT,
                                                        keep.rownames = TRUE)) %>%
                            #remove duplicated rownames (classified by different ref.maps)
                            dplyr::filter(!duplicated(rn)) %>%
                            tibble::column_to_rownames("rn")
    object@meta.data <- merge(object@meta.data, functional.clusters, by = 0, all.x = T) %>%
                        tibble::column_to_rownames("Row.names")

    # save names of reference maps run
    object@misc[["Projectils_param"]] <- names(ref.maps)
  } else {
    object@meta.data <- object@meta.data %>%
                        dplyr::mutate(functional.cluster = NA,
                                      functional.cluster.conf = NA)
  }
  }


  return(object)

}


#############################################################################################################
# function to match words with grep
match_dictionary <- function(cell_type,
                             dictionary = NULL) {

  if(is.null(dictionary)){
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
              "Strom" = "Stromal"
    )
  }

  for (keyword in names(dict)) {
    if (grepl(keyword, cell_type, ignore.case = TRUE)) {
      cell_type <- dict[[keyword]]
      break
    }
  }
  return(cell_type)
}

#############################################################################################################
# complete function applying previous function with sapply to all vector
StandardizeCellNames <- function(cell.names, dictionary = NULL){

  standarized_cellnames <- sapply(cell.names,
                                  match_dictionary,
                                  dictionary)
  standarized_cellnames <- unname(standarized_cellnames)
  return(standarized_cellnames)

}




#############################################################################################################
# Compute the clustering score

get.cluster.score <- function(matrix = NULL,
                              metadata = NULL,
                              cluster.by = c("celltype", "sample"),
                              ndim = NULL,
                              nVarGenes = 500,
                              black.list = "default",
                              ntests = 0,
                              score = c("silhouette"),
                              dist.method = "euclidean"
                              ){

  if(is.null(matrix) || !is.matrix(matrix)){
    stop("Please provide a matrix object.")
  }

  if(is.null(metadata) || !is.data.frame(metadata)){
    stop("Please provide a metadata object as dataframe")
  }

  # Get black list
  if(is.null(black.list) | black.list == "default"){
    black.list <- data("default_black_list")
  }
  black.list <- unlist(black.list)

  # Remove black listed genes from the matrix
  matrix <- matrix[!rownames(matrix) %in% black.list,]



  # dataframe to store score result
  cnames <- c("grouping", "dist_method", "score_method")
  df.score <- expand.grid(cluster.by, dist.method, score)
  # convert to characters, not factors
  df.score <- lapply(df.score, as.character) %>% as.data.frame()
  names(df.score) <- cnames

  # empty list to fill in the loop
  scores <- list()

  # compute common PCA space using DEseq2
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

  # get top variable genes
  rv <- MatrixGenerics::rowVars(vsd)
  select <- order(rv, decreasing=TRUE)[seq_len(min(nVarGenes, length(rv)))]
  vsd <- vsd[select,]

  # remove samples with low variability if needed
  tryCatch({
    pc <- stats::prcomp(t(vsd))
  },
  error = function(e){
    near_zero_var <- caret::nearZeroVar(vsd)
    if(length(near_zero_var) > 0){
      vsd <- vsd[, -near_zero_var]
      metadata <<- metadata[-near_zero_var, ]
      pc <<- stats::prcomp(t(vsd))
    } else if (length(near_zero_var) == ncol(vsd)) {
      message("PCA compute not possible.\n")
    }
  }
  )


  # produce scree plot to know how many dimensions to use
  eigen_val <- pc$sdev^2
  # Filter only eigen values above 1 (Kaiser rule)
  eigen_val <- eigen_val[eigen_val > 1]
  prop_var <- eigen_val / sum(eigen_val)
  # compute the accumulation of the proportional variance
  prop_var_cum <- cumsum(prop_var)

  plot_var <- data.frame(Proportion_variance = prop_var_cum,
                         PC = 1:length(eigen_val)) %>%
    ggplot2::ggplot(ggplot2::aes(PC, Proportion_variance)) +
    ggplot2::geom_col(color = "lightblue")+
    ggplot2::geom_line()+
    ggplot2::geom_vline(xintercept = ndim,
                        color = "red",
                        linetype = 2) +
    ggplot2::geom_label(ggplot2::aes(label = ifelse(PC == ndim,
                                                 round(Proportion_variance,2), NA))) +
    ggplot2::theme_bw() +
    ggtitle(paste0(ndim, " first PC"))

  # compute distance
  for(x in 1:nrow(df.score)){

    # define the grouping by variable
    gr.by <- as.character(df.score[x, 1])

    dist <- stats::dist(pc$x[,1:ndim],
                        method = df.score[x, 2])

    # do not show legend if too many groups
    leg.pos <- ifelse(length(unique(metadata[[gr.by]])) < 30,
                      "right", "none")

    if(df.score[x,3] == "silhouette"){
      silh <- silhoutte_onelabel(metadata[[gr.by]],
                                 dist = dist,
                                 ntests = ntests)

      whole.mean <- mean(silh$cell$sil_width)

      gpl <- silh$cell %>%
        dplyr::rename(!!gr.by := cluster) %>%
        ggplot2::ggplot(ggplot2::aes(rowid, sil_width, fill = .data[[gr.by]])) +
        ggplot2::geom_col() +
        ggplot2::geom_hline(yintercept = whole.mean,
                            color = "black",
                            linetype = 2) +
        ggplot2::labs(y = "Silhouette width",
             title = paste0(gr.by, " - ", df.score[x, 2], "\nAverage Silhoutte width: ", round(whole.mean, 3))) +
        ggplot2::guides(fill=guide_legend(ncol=3))+
        ggplot2::theme(
          panel.background = element_rect(fill = "white"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = leg.pos
        )

    }

    if(df.score[x,3] == "modularity"){
      # transform to network the distance
      graph <- graph.adjacency(
                                as.matrix(as.dist(cor(matrix,
                                                      method="pearson"))),
                                mode="undirected",
                                weighted=TRUE,
                                diag=FALSE
                              )
      # simplify graph
      graph <- simplify(graph, remove.multiple=TRUE, remove.loops=TRUE)

      # Colour negative correlation edges as blue
      E(graph)[which(E(graph)$weight<0)]$color <- "darkblue"
      # Colour positive correlation edges as red
      E(graph)[which(E(graph)$weight>0)]$color <- "darkred"
      # Convert edge weights to absolute values
      E(graph)$weight <- abs(E(graph)$weight)

      # set grouping variable
      V(graph)$group <- as.numeric(as.factor(metadata[[paste0(gr.by,"N")]]))
      # calculate modularity
      mod_score <- igraph::modularity(graph, V(graph)$group)

      V(graph)$label <- NA

      mod.pl <- plot(graph,
                     vertex.color = V(graph)$group,
                     main = "Network with Group Assignments")




    }

    # # plot for PCA
    pc_sum <- summary(pc)
    PC1_varexpl <- pc_sum$importance[2,"PC1"]
    PC2_varexpl <- pc_sum$importance[2,"PC2"]

    # get first 2 PC
    pc.df <- pc$x[,1:2] %>% as.data.frame() %>%
      tibble::rownames_to_column("sample_celltype") %>%
      left_join(., metadata %>% tibble::rownames_to_column("sample_celltype"),
                by = "sample_celltype")

    pc.pl <- pc.df %>%
      ggplot2::ggplot(ggplot2::aes(PC1, PC2, color = .data[[gr.by]])) +
      ggplot2::geom_point() +
      ggplot2::guides(color=guide_legend(ncol=2))+
      labs(title = paste0(gr.by, " - ", " PCA"),
           y = paste0("PC2 (", PC2_varexpl*100, " %)"),
           x = paste0("PC1 (", PC1_varexpl*100, " %)"))+
      ggplot2::theme(
        panel.background = element_rect(fill = "white"),
        legend.position = leg.pos,
        legend.key = element_rect(fill = "white")
      )
  score.type <- as.character(df.score[x,3])
  pl.list <- list( score.type = gpl,
                  "PCA" = pc.pl,
                  "Scree_plot" = plot_var)


  # plot of bootstraping
  conf.pl <- silh$summary %>%
              dplyr::filter(iteration != "NO") %>%
              ggplot2::ggplot(ggplot2::aes(x = avg_sil_width,
                                           y = ..density..,
                                           fill = cluster)) +
              ggplot2::geom_density(alpha = 0.6, show.legend = F) +
              ggplot2::geom_ribbon(aes(ymin = 0, ymax = ..density..), alpha = 0.05)
              ggplot2::facet_wrap(~cluster, ncol = 2) +
              ggplot2::theme_bw()

              ggplot2::geom_vline(aes(xintercept = ifelse(iteration == "NO",
                                                      avg_sil_width, NA)),
                                  color = "Avg Silhoutte width"),
                                  lty = 2, show.legend = T)
              ggplot2::geom_vline(xintercept = quantile(cof_before[,2], c(0.05,0.95))[1],
                             color = "Bootstrap"), lty = 2, show.legend = T)

    # return list
    ret <- list("whole_avgerage" = whole.mean,
                "bygroup_average" = silh$summary,
                "plots" = pl.list)

    scores[[paste(df.score[x,], collapse = "_")]] <- ret
  }


  return(scores)


}

#############################################################################################################
# Function for silhoutte one label

silhoutte_onelabel <- function(labels = NULL, # vector of labels
                               dist = NULL, # distance object
                               ntests = 0, # number of shuffling events
                               seed = 22){# seed for random suffling

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

  tlabels <- unique(labels)

  sil.all <- data.frame(matrix(nrow = 0, ncol = 3))
  names(sil.all) <- c("cluster", "neighbor", "sil_width")


  sil.sum <- data.frame(matrix(nrow = 0, ncol = 4))
  names(sil.sum) <- c("cluster", "iteration", "size", "avg_sil_width")

  len <- length(labels)

  for(a in tlabels){
    x_one <- ifelse(labels == a, 2, 1) %>%
              as.factor() %>% as.numeric()

    silh <- cluster::silhouette(x_one, dist)

    # change names back to character
    sil.res <- as.data.frame(silh) %>%
              dplyr::filter(cluster == 2) %>%
              dplyr::mutate(cluster = a)
    #join for all samples
    sil.all <- rbind(sil.all, sil.res)

    size <- nrow(sil.res)

    sil.sumA <- data.frame(cluster = a,
                           iteration = "NO",
                           size = size,
                           avg_sil_width = mean(sil.res$sil_width))

    bots.df <- data.frame(matrix(nrow = 0, ncol = 4))
    names(bots.df) <- c("cluster", "iteration", "size", "avg_sil_width")

    if(ntests > 0){
    # perform the shuffling
      for(u in 1:ntests){
        # random vector
        vec <- rep(1, length(labels))
        # seeding for reproducibility
        set.seed(seed+u)
        random_sample <- sample(1:len, size)
        vec[random_sample] <- 2
        vec <- vec %>% as.factor() %>% as.numeric()

        # run silhoutte
        silh <- cluster::silhouette(vec, dist)
        sil.res <- as.data.frame(silh) %>%
                    dplyr::filter(cluster == 2)

        sil.sumB <- data.frame(cluster = a,
                               iteration = u,
                               size = size,
                               avg_sil_width = mean(sil.res$sil_width))

        bots.df <- rbind(bots.df, sil.sumB)

      }

      # run wilcox.test in case distribution is not normal (possible)
      wt <-

      bots.df <- bots.df %>%
                  mutate(conf.int.95 = I(list(quantile(avg_sil_width, c(0.025,0.975)))))



    }
  }
  # order silhoute width per cell type
  sil.all <- sil.all %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(sil_width), .by_group = T) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rowid = dplyr::row_number())


  ret.list <- list("cell" = sil.all,
                   "summary" = sil.sum)

  return(ret.list)
}






