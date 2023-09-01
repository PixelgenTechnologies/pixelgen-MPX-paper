

theme_umap <- 
  function (base_size = 11, base_family = "", base_line_size = base_size/22, 
            base_rect_size = base_size/22) {
    theme_grey(base_size = base_size, base_family = base_family, 
               base_line_size = base_line_size, 
               base_rect_size = base_rect_size) %+replace% 
      theme(panel.background = element_rect(fill = "white", 
                                            colour = NA), 
            panel.border = element_rect(fill = NA, colour = "grey75"), 
            panel.grid = element_blank(), 
            panel.grid.minor = element_line(linewidth = rel(0.5)), 
            strip.background = element_rect(fill = "grey85", 
                                            colour = "grey85"), 
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.key = element_rect(fill = "white", 
                                      colour = NA), 
            complete = TRUE)
  }


theme_pg <- 
  function (base_size = 11, base_family = "", base_line_size = base_size/22, 
            base_rect_size = base_size/22) {
    theme_grey(base_size = base_size, base_family = base_family, 
               base_line_size = base_line_size, 
               base_rect_size = base_rect_size) %+replace% 
      theme(panel.background = element_rect(fill = "white", 
                                            colour = NA), 
            panel.border = element_rect(fill = NA, colour = "grey75"), 
            panel.grid = element_blank(), 
            panel.grid.minor = element_line(linewidth = rel(0.5)), 
            strip.background = element_rect(fill = "grey85", 
                                            colour = "grey85"), 
            legend.key = element_rect(fill = "white", 
                                      colour = NA), 
            complete = TRUE)
  }