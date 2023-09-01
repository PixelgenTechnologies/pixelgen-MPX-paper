




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
