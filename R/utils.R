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
                              cluster.by = c("sample", "celltype"),
                              ndim = 10,
                              score = c("silhouette"),
                              dist.method = "euclidean"
                              ){

  if(is.null(matrix) || !is.matrix(matrix)){
    stop("Please provide a matrix object.")
  }

  if(is.null(metadata) || !is.data.frame(metadata)){
    stop("Please provide a metadata object as dataframe")
  }



  # dataframe to store score result
  cnames <- c("grouping", "dist_method", "score_method")
  df.score <- expand.grid(cluster.by, dist.method, score)
  # convert to characters, not factors
  df.score <- lapply(df.score, as.character) %>% as.data.frame()
  names(df.score) <- cnames

  # empty list to fill in the loop

  scores <- list()

  # convert metadata grouping to numeric and factor
  for(n in cluster.by){
    metadata[[paste0(n,"N")]] <- as.numeric(as.factor(metadata[[n]]))
  }


  # compute common PCA space

  mat.scaled <- scale(matrix)


  # remove samples with low variability if needed
  tryCatch({
    pc <- stats::prcomp(t(mat.scaled))
  },
  error = function(e){
    near_zero_var <- caret::nearZeroVar(mat.scaled)
    if(length(near_zero_var) > 0){
      mat.scaled <- mat.scaled[, -near_zero_var]
      metadata <<- metadata[-near_zero_var, ]
      pc <<- stats::prcomp(t(mat.scaled))
    } else if (length(near_zero_var) == ncol(mat.scaled)) {
      message("PCA compute not possible.\n")
    }
  }
  )



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
      silh <- cluster::silhouette(as.numeric(as.factor(metadata[[paste0(gr.by,"N")]])),
                          dist)

      sil.df <- as.data.frame(silh) %>%
                dplyr::rename(!!paste0(gr.by,"N") := cluster) %>%
                left_join(., dplyr::distinct(metadata[,c(paste0(gr.by,"N"), gr.by)]),
                        by = paste0(gr.by,"N")) %>%
                dplyr::group_by_at(dplyr::vars(!!gr.by)) %>%
                dplyr::arrange(desc(sil_width), .by_group = T) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(rowid = dplyr::row_number())

      whole.mean <- mean(sil.df$sil_width)

      g.df.sum <- sil.df %>%
                    dplyr::group_by_at(dplyr::vars(!!gr.by)) %>%
                    dplyr::reframe(
                      size = n(),
                      average.sil.width = mean(sil_width),
                      !!gr.by := .data[[gr.by]]
                      ) %>%
                    dplyr::distinct()

      gpl <- sil.df %>%
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

    # run dimensional reduction on distances
    # mds <- stats::cmdscale(dist) %>%
    #       as.data.frame() %>%
    #       tibble::rownames_to_column("sample_celltype") %>%
    #       left_join(., metadata %>% tibble::rownames_to_column("sample_celltype"),
    #                 by = "sample_celltype")
    # mds.pl <- mds %>%
    #       ggplot2::ggplot(ggplot2::aes(V1, V2, color = .data[[gr.by]])) +
    #       geom_point() +
    #   ggplot2::guides(color=guide_legend(ncol=2))+
    #   labs(title = paste0(gr.by, " - ", df.score[x, 2], " distance. Multidimensional scaling")) +
    #   ggplot2::theme(
    #     panel.background = element_rect(fill = "white"),
    #     axis.text = element_blank(),
    #     axis.title = element_blank(),
    #     axis.ticks = element_blank(),
    #     legend.position = leg.pos,
    #     legend.key = element_rect(fill = "white")
    #   )

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
      geom_point() +
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
                  "PCA" = pc.pl)

    # return list
    ret <- list("whole_avgerage" = whole.mean,
                "bygroup_average" = g.df.sum,
                "plots" = pl.list)

    scores[[paste(df.score[x,], collapse = "_")]] <- ret
  }


  return(scores)


}







