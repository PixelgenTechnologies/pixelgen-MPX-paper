




convert_pixelgenh5_seuratrds <- 
  function(filename, newfilename) {
    require(tidyverse)
    require(rhdf5)
    require(Seurat)
    require(SeuratDisk)
    
    original_filename <- filename
    if(endsWith(filename, ".dataset.pxl")) {
      filename <- unzip(filename, "adata.h5ad", exdir = tempdir())
      
    }
    
    anndata_hier <- 
      h5read(filename, "/")  
    
    
    if(endsWith(original_filename, ".dataset.pxl")) {
      file.remove(filename)
      
    }
    
    X <- anndata_hier$X
    colnames(X) <-   anndata_hier$obs$component
    rownames(X) <-   anndata_hier$var$marker
    
    if(is.null(rownames(X)[1])) {
      rownames(X) <-   anndata_hier$var$`_index`
    }
    
    seur_obj <- 
      CreateSeuratObject(counts = X)
    
    if(is.null(anndata_hier$layers)) {
      
      # Since 0.6.3 normalized data is stored in obsm:
      for(layr in names(anndata_hier$obsm)) {
        
        if(!layr %in% c("denoised",
                        "normalized_clr", 
                        "normalized_rel")) next
        
        X_assay <- 
          anndata_hier$obsm[[layr]] %>% 
          as_tibble() %>% 
          column_to_rownames("component") %>% 
          t() 
        
        assay <- CreateAssayObject(counts = X_assay)
        
        seur_obj[[layr]] <- assay
      }
      
    } else {
      # In < 0.6.3 normalized data is stored in layers:
      for(layr in names(anndata_hier$layers)) {
        X_assay <- anndata_hier$layers[[layr]]
        colnames(X_assay) <- anndata_hier$obs$component
        rownames(X_assay) <- anndata_hier$var$marker
        
        if(is.null(rownames(X_assay)[1])) {
          rownames(X_assay) <-   anndata_hier$var$`_index`
        }
        
        
        assay <- CreateAssayObject(counts = X_assay)
        
        seur_obj[[layr]] <- assay
      }
    }
    
    
    
    seur_obj@meta.data <- 
      anndata_hier$obs %>% 
      {
        newdata <- .;
        
        
        for(name in names(newdata)) {
          if(length(newdata[[name]]) == dim(anndata_hier$X)[2]) next
          
          indicies <- newdata[[name]]$codes + 1
          categories <- newdata[[name]]$categories
          
          newdata[[name]] <-factor(categories[indicies], categories)
          
          
        }
        
        newdata
      } %>% 
      as.data.frame() %>% 
      select(any_of(c("COMPONENT", 
                      "VERTICES",
                      "EDGES",
                      "ANTIBODIES",
                      "UPIA",
                      "UPIB",
                      "READS_COUNT",
                      "MEAN_READS_PER_EDGE",
                      "MEAN_UPIA_TO_UPIB_EDGE_DEGREE",
                      "MEDIAN_UMI_PER_UPIA",
                      "UMI",
                      "is_filtered")),
             everything()) %>% 
      column_to_rownames("component")
    
    var_meta <- 
      anndata_hier$var %>% 
      as.data.frame()
    rownames(var_meta) <- anndata_hier$var$MARKER
    
    seur_obj@assays$RNA@meta.features <- 
      var_meta
    
    
    saveRDS(seur_obj, newfilename)
    
  }



read_edgelist <- 
  function(filename) {
    require(tidyverse)
    
    message(paste("Reading edgelist:", filename))
    
    unzipped_filename <- 
      unzip(filename, "edgelist.csv.gz", exdir = tempdir())
    
    if(length(unzipped_filename) != 0) {
      
      edgelist <- read_csv(unzipped_filename)
      
    } else {
      require(arrow)
      
      unzipped_filename <- 
        unzip(filename, "edgelist.parquet", exdir = tempdir())
      edgelist <- read_parquet(unzipped_filename)
      
    }
    
    file.remove(unzipped_filename)
    
    return(edgelist)
    
  }

read_pixeldata_item <- 
  function(filename, items = c("colocalization", "polarization", "edgelist")) {
    require(tidyverse)
    
    message(paste("Reading from:", filename))
    read_items <- 
      items %>% 
      set_names(., .) %>% 
      lapply(function(item) {
        message(paste0(item, "..."))
        
        unzipped_filename <- 
          unzip(filename, paste0(item, ".csv.gz"), exdir = tempdir())
        
        
        if(length(unzipped_filename) != 0) {
          
          outdata <- read_csv(unzipped_filename)
          
        } else {
          require(arrow)
          
          unzipped_filename <- 
            unzip(filename,  paste0(item, ".parquet"), exdir = tempdir())
          outdata <- read_parquet(unzipped_filename)
          
        }
        
        file.remove(unzipped_filename)
        
        return(outdata)
      })
    
    
    
    if(length(read_items) == 1) {
      read_items <- read_items[[1]]
    }
    
    return(read_items)
    
  }

expand_frames <- 
  function(x, y, name_sep = "_") {
    expand_grid(x, y, .name_repair = function(x) paste(x, as.numeric(duplicated(x)) + 1,
                                                       sep = name_sep))
  }

cluster_long_data <-  
  function(long_data,
           distance_method = "euclidean",
           clustering_method = "ward.D2", 
           cluster_rows = T,
           cluster_cols = T, 
           fill = NA) {
    suppressMessages(require(tidyverse))
    
    wide_data <- 
      long_data %>% 
      select(1:3) %>% 
      spread(2, 3, fill = fill) %>% 
      column_to_rownames(names(long_data)[1])
    
    order_row <- 
      rownames(wide_data)
    order_col <- 
      colnames(wide_data)
    
    
    if(cluster_rows) {
      order1 <- 
        wide_data %>% 
        dist(method = distance_method) %>% 
        hclust(method = clustering_method) %>% 
        with(labels[order])
    }
    
    if(cluster_cols) {
      order2 <- 
        wide_data %>% 
        t() %>% 
        dist(method = distance_method) %>% 
        hclust(method = clustering_method) %>% 
        with(labels[order])
    }
    
    long_data %>% 
      rename(v1 = 1, 
             v2 = 2,
             val = 3) %>% 
      mutate(v1 = factor(v1, order1),
             v2 = factor(v2, order2)) %>% 
      set_names(names(long_data))
    
  }



edgelist_to_simple_bipart_graph <- 
  function(edgelist) {
    
    # A function to create a simple bipartite graph from an edgelist.
    
    require(igraph)
    require(tidygraph)
    
    # Simplify edgelist 
    edgelist <- 
      edgelist %>% 
      select(upia = any_of(c("upia", "UPIA")),
             upib = any_of(c("upib", "UPIB"))) %>% 
      distinct() 
    
    
    # Create bipart graph
    bipart_graph <- 
      edgelist %>% 
      # Create tidy graph
      as_tbl_graph(directed = F) %>% 
      mutate(node_type = ifelse(name %in% edgelist$upia,
                                "A",
                                "B"))
    
    bipart_graph
    
  }

edgelist_to_simple_Anode_graph <- 
  function(edgelist) {
    
    # A function to create an A node projection graph from an edgelist.
    
    require(igraph)
    require(tidygraph)
    
    # Simplify edgelist 
    edgelist <- 
      edgelist %>% 
      select(upia = any_of(c("upia", "UPIA")),
             upib = any_of(c("upib", "UPIB"))) %>% 
      distinct() 
    
    # Create Anode graph
    anode_graph <- 
      
      edgelist %>% 
      left_join(edgelist,
                by = c("upib"),
                relationship = "many-to-many",
                suffix = c("1", "2")) %>% 
      select(upia1, upia2) %>% 
      distinct() %>% 
      filter(upia1 < upia2) %>% 
      
      # Create tidy graph
      as_tbl_graph(directed = F) 
    
    anode_graph
    
  }

get_component_layout <- 
  function(component_edge_list, layout_method = "fr", dim = 2, 
           normalize_layout = T, k = 0, pivots = 50, weighted = F) {
    
    # Calculates a graph layout for a component's edgelist, and outputs a list with the bipartite graph, layout, and
    # antibody counts per A pixel. 
    
    # layout_method     - The method for calculating the graph layout; Fruchterman-Reingold (fr), Kamada-Kawai (kk), drl.
    # dim               - The number of dimension for the layout, 2 or 3.
    # normalize_layout  - Whether the coordinate system for the layout should be scaled to -1 to 1.
    # k                 - The size of the neighborhood from which to pool counts from in the UPIA antibody count table. 
    #                     Zero is recommended.
    
    require(igraph)
    require(tidygraph)
    
    if(layout_method == "fr") {
      layout_met <- layout_with_fr
    } else if(layout_method == "kk") {
      layout_met <- layout_with_kk
    } else if(layout_method == "drl") {
      layout_met <- layout_with_drl
    } else if(layout_method == "pmds") {
      layout_met <- graphlayouts::layout_with_pmds
    }
    
    
    component_counts <- 
      component_edge_list %>% 
      create_node_markers_counts(k = k)
    
    simple_graph <- 
      component_edge_list %>% 
      select(upia, upib) %>% 
      distinct() %>% 
      edgelist_to_simple_bipart_graph()
    
    if(weighted) {
      
      node_degree <- 
        component_edge_list %>% 
        select(upia, upib) %>% 
        distinct() %>% 
        {list(A = rename(., name = upia),
              B = rename(., name = upib))} %>% 
        map(. %>% 
              group_by(name) %>% 
              count()) %>% 
        bind_rows(.id = "node_type") %>% 
        ungroup() 
      
      node_degree <- 
        node_degree[match(as_tibble(simple_graph)$name, node_degree$name),]
      
      weights <- 
        simple_graph %>%  
        activate(edges) %>% 
        mutate(from_degree = node_degree$n[from],
               to_degree = node_degree$n[to],
               weight = ifelse(from_degree > to_degree, from_degree, to_degree)) %>% 
        as_tibble() %>% 
        pull(weight)
      
      if(layout_method %in% c("fr", "drl")) {
        weights <- 1 / weights
      }
      
    } else {
      weights <- NA
    }
    
    if(layout_method == "pmds") {
      
      graph_layout <- 
        simple_graph %>% 
        layout_met(dim = dim, pivots = pivots, weights = weights) %>% 
        as_tibble(.name_repair = function(x) c("x", "y", "z")[1:length(x)]) %>% 
        bind_cols(as_tibble(simple_graph)) %>% 
        select(upi = name, node_type, any_of(c("x", "y", "z")))
      
    } else {
      
      graph_layout <- 
        simple_graph %>% 
        layout_met(dim = dim, weights = weights) %>% 
        as_tibble(.name_repair = function(x) c("x", "y", "z")[1:length(x)]) %>% 
        bind_cols(as_tibble(simple_graph)) %>% 
        select(upi = name, node_type, any_of(c("x", "y", "z")))
      
    }
    
    
    if(normalize_layout) {
      
      if(dim == 2) {
        scale_fact <- max(abs(c(graph_layout$x, 
                                graph_layout$y)))
        graph_layout$x <- graph_layout$x / scale_fact
        graph_layout$y <- graph_layout$y / scale_fact
      }
      
      if(dim == 3) {
        scale_fact <- max(abs(c(graph_layout$x, 
                                graph_layout$y,
                                graph_layout$z)))
        graph_layout$x <- graph_layout$x / scale_fact
        graph_layout$y <- graph_layout$y / scale_fact
        graph_layout$z <- graph_layout$z / scale_fact
      }
      
      
    }
    
    return(
      list(graph = simple_graph,
           layout = graph_layout,
           upi_counts = component_counts)
    )
  }

create_node_markers_counts <-  
  function(component_edge_list, k = 0) {
    
    
    # Computes and returns a data frame of antibody counts per
    # node (vertex) of the A node graph given a component edge list as input.
    # The parameter k allows to include neighbors (of each node) when computing the
    # counts. K defines the number of levels when searching neighbors. 
    
    
    component_counts <- 
      component_edge_list %>% 
      group_by(upia, marker) %>% 
      count() %>% 
      spread(marker, n, fill = 0) %>% 
      column_to_rownames("upia")
    
    simple_graph <- 
      component_edge_list %>% 
      select(upia, upib) %>% 
      distinct() %>% 
      edgelist_to_simple_Anode_graph()
    
    if(k > 0) {
      adj_mat <- 
        connect(simple_graph, order = k) %>% 
        get.adjacency()
      
      diag(adj_mat) <- 1
      
      return(as.data.frame(as.matrix(adj_mat %*% as.matrix(component_counts))))
      
    } else {
      
      return(as.matrix(component_counts))
      
    }
    
  }



find_com <- 
  function(x, y, z, weight) {
    c(x = weighted.mean(x, w = weight),
      y = weighted.mean(y, w = weight),
      z = weighted.mean(z, w = weight))
  }

procrustes_align <- 
  function(component_list, anchor_markers = NULL, reference = 1, do_scale = F) {
    
    # Function to align multiple components to a reference component (defaults to first in list) 
    # using a center of mass calculation weighted by the marker count. The center of mass is 
    # calculated for all anchor markers, and are aligned using procrustes analysis. By default, 
    # the layouts are only rotated, but they can be scaled by setting do_scale to TRUE.
    
    
    if(is.null(anchor_markers)) {
      component_list_anchors <- component_list
    } else {
      component_list_anchors <- 
        component_list %>% 
        map(. %>% 
              filter(marker %in% anchor_markers))
    }
    
    
    
    component_ids <- 
      component_list %>% 
      map(. %>% 
            pull(component) %>% 
            unique()) %>% 
      unlist()
    
    
    comp_total_center_of_mass <-
      component_list_anchors %>%
      map(. %>% 
            do({
              find_com(.$x,
                       .$y,
                       .$z,
                       .$count) %>% 
                enframe() %>% 
                spread(name, value) %>% 
                {x <- as.data.frame(.); rownames(x) <- "component_center"; x}
            }) %>% 
            ungroup()) %>% 
      set_names(component_ids) 
    
    comp_center_of_mass <-
      component_list_anchors %>%
      imap(~.x %>% 
             group_by(marker) %>% 
             do({
               find_com(.$x,
                        .$y,
                        .$z,
                        .$count) %>% 
                 enframe() %>% 
                 spread(name, value)
             }) %>% 
             ungroup() %>% 
             column_to_rownames("marker") %>% 
             bind_rows(comp_total_center_of_mass[[.y]])) %>% 
      set_names(component_ids)
    
    
    component_list <- 
      component_list %>% 
      set_names(component_ids)
    
    reference_center_of_mass <- comp_center_of_mass[[reference]] 
    
    # Add origo for anchor markers that are missing from some components 
    comp_center_of_mass <- 
      comp_center_of_mass %>% 
      lapply(function(center_of_mass) {
        if(nrow(center_of_mass) == nrow(reference_center_of_mass)) {
          return(center_of_mass)
        } else {
          missing_markers <- rownames(reference_center_of_mass)[!rownames(reference_center_of_mass) %in% 
                                                                  rownames(center_of_mass)]
          add_marker_matrix <- 
            # matrix(c(0,0,0), nrow = length(missing_markers), ncol = 3) 
            rep(center_of_mass["component_center", ], length(missing_markers)) %>% 
            matrix(nrow = 3) %>% 
            t()
          
          rownames(add_marker_matrix) <- missing_markers
          colnames(add_marker_matrix) <- c("x", "y", "z")
          
          new_center_of_mass <- 
            rbind(center_of_mass, add_marker_matrix) %>% 
            apply(MARGIN = 2, as.numeric) 
            
          rownames(new_center_of_mass) <- rownames(rbind(center_of_mass, add_marker_matrix))
            
          return(new_center_of_mass)
        }
      })
    
    procrustes_objs <- 
      lapply(comp_center_of_mass,
             function(component_center_of_mass) {
               vegan::procrustes(reference_center_of_mass,
                                 component_center_of_mass[rownames(reference_center_of_mass),], scale = do_scale)
             })  
    
    names(component_list) %>% 
      set_names(., .) %>% 
      lapply(function(component_name) {
        bind_cols(component_list[[component_name]],
                  predict(procrustes_objs[[component_name]], 
                          newdata = component_list[[component_name]] %>% 
                            select(x, y, z)) %>% 
                    as_tibble()) %>% 
          rename(x_rot = V1, y_rot = V2, z_rot = V3)
      }) 
    
  }




get_density <- 
  function(x, y, h = 0.5, n = 25, lims = c(c(-1, 1), c(-1, 1))) {
    density_2d <- 
      MASS::kde2d(x, y, h = h, n = n, lims = lims)
    
    count <- length(x)
    density_2d$z  %>% 
      as_tibble() %>% 
      mutate(y = row_number()) %>% 
      gather(x, density, -y) %>% 
      mutate(x = str_remove(x, "V") %>% as.integer(),
             x = density_2d$x[x],
             y = density_2d$y[y],
             count = count) 
  }
