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
          data("Hs2Mm.convert.table")
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
                              group.by = c("sample", "celltype"),
                              score = c("silhouette"),
                              dist.method = "euclidean",
                              bparam = NULL
                              ){

  if(is.null(matrix) || !is.matrix(matrix)){
    stop("Please provide a matrix object.")
  }

  if(is.null(metadata) || !is.data.frame(metadata)){
    stop("Please provide a metadata object as dataframe")
  }


  score <- score[1]



  # dataframe to store score result
  cnames <- c("grouping", "dist_method", "score_method")
  df.score <- expand.grid(group.by, dist.method, score)
  names(df.score) <- cnames

  # convert metadata grouping to numeric and factor
  for(n in group.by){
    metadata[[paste0(n,"N")]] <- as.numeric(as.factor(metadata[[n]]))
  }

  # compute distance
  scores <- BiocParallel::bplapply(
                X = 1:nrow(df.score),
                BPPARAM = bparam,
                function(x){
                  # define the grouping by variable
                  gb <- as.character(df.score[x, 1])

                  dist <- stats::dist(t(matrix),
                               method = df.score[x, 2])

                  # do not show legend if too many groups
                  leg.pos <- ifelse(length(unique(metadata[[gb]])) < 12,
                                    "right", "none")

                  if(df.score[x,3] == "silhouette"){
                    silh <- cluster::silhouette(as.numeric(as.factor(metadata[[paste0(gb,"N")]])),
                                        dist)

                    sum <- summary(silh)

                    sil.df <- as.data.frame(silh) %>%
                              dplyr::rename(!!paste0(gb,"N") := cluster) %>%
                              left_join(., dplyr::distinct(metadata[,c(paste0(gb,"N"), gb)]),
                                      by = paste0(gb,"N")) %>%
                              dplyr::group_by_at(dplyr::vars(!!gb)) %>%
                              dplyr::arrange(desc(sil_width), .by_group = T) %>%
                              dplyr::ungroup() %>%
                              dplyr::mutate(rowid = dplyr::row_number())

                    whole.mean <- mean(sil.df$sil_width)

                    g.df.sum <- sil.df %>%
                                  dplyr::group_by_at(dplyr::vars(!!gb)) %>%
                                  dplyr::summarize(size = n(),
                                                   average.sil.width = mean(sil_width),
                                                   !!gb := .data[[gb]]) %>%
                                  dplyr::distinct()

                    gpl <- sil.df %>%
                      ggplot2::ggplot(ggplot2::aes(rowid, sil_width, fill = .data[[gb]])) +
                      ggplot2::geom_col() +
                      ggplot2::geom_hline(yintercept = whole.mean,
                                          color = "black",
                                          linetype = 2) +
                      ggplot2::labs(y = "Silhouette width",
                           title = paste0(gb, " - Average Silhoutte width: ", round(sil.mean, 3))) +
                      ggplot2::guides(fill=guide_legend(ncol=2))+
                      ggplot2::theme(
                        panel.background = element_rect(fill = "white"),
                        axis.text.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        legend.position = leg.pos
                      )

                  } else if(df.score[x,3] == "modularity"){
                    # adjancy matrix
                    adjacency_matrix <- as.matrix(matrix > 1)
                    ncol(adjacency_matrix) != nrow(adjacency_matrix)

                    # Create a graph object
                    graph <- igraph::graph_from_adjacency_matrix(matrix, mode = "undirected")

                    # Apply a community detection algorithm (e.g., Louvain)
                    membership <- igraph::cluster_louvain(graph)$membership

                    # Calculate modularity
                    modularity_value <- igraph::modularity(graph, membership)
                  }

                # run dimensional reduction on distances
                mds <- stats::cmdscale(dist) %>%
                      as.data.frame() %>%
                      tibble::rownames_to_column("sample_celltype") %>%
                      left_join(., metadata %>% tibble::rownames_to_column("sample_celltype"),
                                by = "sample_celltype")
                mds.pl <- mds %>%
                      ggplot2::ggplot(ggplot2::aes(V1, V2, color = .data[[gb]])) +
                      geom_point() +
                  ggplot2::guides(color=guide_legend(ncol=2))+
                  labs(title = paste0(gb, " - ", df.score[x, 2], " distance. Multidimensional scaling")) +
                  ggplot2::theme(
                    panel.background = element_rect(fill = "white"),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks = element_blank(),
                    legend.position = leg.pos,
                    legend.key = element_rect(fill = "white")
                  )

                # # plot for PCA
                # pc <- stats::prcomp(t(matrix))
                # pc_sum <- summary(pc)
                # PC1_varexpl <- pc_sum$importance[2,"PC1"]
                # PC2_varexpl <- pc_sum$importance[2,"PC2"]
                #
                # # get first 2 PC
                # pc.df <- pc$x[,1:2] %>% as.data.frame() %>%
                #   tibble::rownames_to_column("sample_celltype") %>%
                #   left_join(., metadata %>% tibble::rownames_to_column("sample_celltype"),
                #             by = "sample_celltype")
                #
                # pc.pl <- pc.df %>%
                #   ggplot2::ggplot(ggplot2::aes(PC1, PC2, color = .data[[gb]])) +
                #   geom_point() +
                #   ggplot2::guides(color=guide_legend(ncol=2))+
                #   labs(title = paste0(gb, " - ", " PCA")) +
                #   ggplot2::theme(
                #     panel.background = element_rect(fill = "white"),
                #     axis.text = element_blank(),
                #     axis.title = element_blank(),
                #     axis.ticks = element_blank(),
                #     legend.position = leg.pos,
                #     legend.key = element_rect(fill = "white")
                #   )



                # return list
                ret <- list("whole_avgerage" = whole.mean,
                            "bygroup_average" = g.df.sum,
                            "plots" = list(df.score[x,3] = gpl,
                                           "cmdscale" = mds.p,
                                           "PCA" = pc.pl))
                return(ret)
                }
  )

  names(scores) <- apply(df.score, 1, paste, collapse = "_")

  return(scores)





}







