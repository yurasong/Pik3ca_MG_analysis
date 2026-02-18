make_mushroom_plot2 <- function(prioritization_table, top_n = 30, show_rankings = FALSE,
                                show_all_datapoints = FALSE, true_color_range = TRUE,
                                use_absolute_rank = FALSE,
                                size = "scaled_avg_exprs", color = "scaled_p_val_adapted",
                                ligand_fill_colors = c("#DEEBF7", "#08306B"),
                                receptor_fill_colors = c("#FEE0D2", "#A50F15"),
                                unranked_ligand_fill_colors = c(alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2)),
                                unranked_receptor_fill_colors = c( alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2)),
                                ...){
  size_ext <-  c("ligand", "receptor"); color_ext <- c("ligand", "receptor")
  if (size == "pct_expressed") size_ext <- c("sender", "receiver")
  if (color == "pct_expressed") color_ext <- c("sender", "receiver")
  
  cols_to_use <- c("sender", "ligand", "receptor", paste0(size, "_", size_ext), paste0(color, "_", color_ext))
  
  if (!all(cols_to_use %in% colnames(prioritization_table))){
    stop(paste(paste0("`", cols_to_use %>% .[!. %in% colnames(prioritization_table)], "`", collapse =", "), "column not in prioritization table"))
  }
  if(!is.logical(show_rankings) | length(show_rankings) != 1)
    stop("show_rankings should be a TRUE or FALSE")
  if(!is.logical(show_all_datapoints) | length(show_all_datapoints) != 1)
    stop("show_all_datapoints should be a TRUE or FALSE")
  if(!is.logical(true_color_range) | length(true_color_range) != 1)
    stop("true_color_range should be a TRUE or FALSE")
  if(!is.logical(use_absolute_rank) | length(use_absolute_rank) != 1)
    stop("use_absolute_rank should be a TRUE or FALSE")
  if(!is.numeric(top_n) | length(top_n) != 1)
    stop("top_n should be a numeric vector of length 1")
  if(length(ligand_fill_colors) != 2)
    stop("ligand_fill_colors should be a vector of length 2")
  if(length(receptor_fill_colors) != 2)
    stop("receptor_fill_colors should be a vector of length 2")
  if(length(unranked_ligand_fill_colors) != 2)
    stop("unranked_ligand_fill_colors should be a vector of length 2")
  if(length(unranked_receptor_fill_colors) != 2)
    stop("unranked_receptor_fill_colors should be a vector of length 2")
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  library("ggnewscale")
  library("ggforce")
  library("shadowtext")
  requireNamespace("cowplot")
  
  if (!"prioritization_rank" %in% colnames(prioritization_table)){
    prioritization_table <- prioritization_table %>% dplyr::mutate(prioritization_rank = rank(desc(prioritization_score)))
  }
  # Add 'relative rank' column which is basically 1:n
  prioritization_table <- prioritization_table %>% dplyr::mutate(relative_rank = rank(desc(prioritization_score)))
  
  # If use_absolute_rank, use 'prioritization_rank' column to filter top_n
  rank_filter_col <- ifelse(use_absolute_rank, "prioritization_rank", "relative_rank")
  
  # Create a new column of ligand-receptor interactions, and filter table to
  # only include LR interactions that appear in the top_n
  filtered_table <- prioritization_table %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = " - "))
  order_interactions <- unique(filtered_table %>% filter(.data[[rank_filter_col]] <= top_n) %>% pull(lr_interaction))
  filtered_table <- filtered_table %>% filter(lr_interaction %in% order_interactions) %>%
    mutate(lr_interaction = factor(lr_interaction, levels = rev(order_interactions)))
  
  # Check if filtered_table is empty
  if (nrow(filtered_table) == 0){
    stop("No ligand-receptor interactions found in the top_n. Please try use_absolute_rank = FALSE or increase top_n.")
  }
  
  # Keep order of senders, if present (if not, sort alphabetically)
  if (!is.factor(filtered_table$sender)){
    filtered_table$sender <- as.factor(filtered_table$sender)
  } else {
    # Drop levels that are not present in the filtered table
    filtered_table$sender <- droplevels(filtered_table$sender)
  }
  
  lr_interaction_vec <- 1:length(order_interactions) %>% setNames(order_interactions)
  
  # Make each ligand and receptor into separate rows (to draw 1 semicircle per row)
  filtered_table <- filtered_table %>% select(c("lr_interaction", all_of(cols_to_use), "prioritization_rank", "relative_rank")) %>%
    pivot_longer(c(ligand, receptor), names_to = "type", values_to = "protein") %>%
    mutate(size = ifelse(type == "ligand", get(paste0(size, "_", size_ext[1])), get(paste0(size, "_", size_ext[2]))),
           color = ifelse(type == "ligand", get(paste0(color, "_", color_ext[1])), get(paste0(color, "_",  color_ext[2])))) %>%
    select(-contains(c("_ligand", "_receptor", "_sender", "_receiver"))) %>%
    mutate(start = rep(c(-pi, 0), nrow(filtered_table))) %>%
    mutate(x = as.numeric(sender), y = lr_interaction_vec[lr_interaction])
  
  # Warning if size column is not scaled between 0 and 1.001
  if (any(filtered_table$size < 0) | any(filtered_table$size > 1.001)){
    stop("Size column is not scaled between 0 and 1. Please use this column as the color instead.")
  }
  
  # Rename size and color columns to be more human-readable
  keywords_adj <- c("LFC", "pval", "", "product", "mean", "adjusted", "expression") %>% setNames(c("lfc", "p", "val", "prod", "avg", "adj", "exprs"))
  size_title <- sapply(stringr::str_split(size, "_")[[1]], function(k) ifelse(is.na(keywords_adj[k]), k, keywords_adj[k])) %>%
    paste0(., collapse = " ") %>%  stringr::str_replace("^\\w{1}", toupper)
  color_title <- sapply(stringr::str_split(color, "_")[[1]], function(k) ifelse(is.na(keywords_adj[k]), k, keywords_adj[k])) %>%
    paste0(., collapse = " ") %>% stringr::str_replace("^\\w{1}", toupper)
  
  color_lims <- c(0,1)
  if (true_color_range) color_lims <- NULL
  
  scale <- 0.5
  
  ncelltypes <- length(unique(filtered_table$sender))
  n_interactions <- length(lr_interaction_vec)
  legend2_df <- data.frame(values = c(0.25, 0.5, 0.75, 1), x=(ncelltypes+2.5):(ncelltypes+5.5), y=rep(floor(n_interactions/3), 4), start=-pi)
  axis_rect <- data.frame(xmin=0, xmax=ncelltypes+1, ymin=0, ymax=n_interactions+1)
  panel_grid_y <- data.frame(x = rep(seq(from = 0.5, to = ncelltypes+0.5, by = 1), each=2),
                             y = c(n_interactions+1, 0), group = rep(1:(ncelltypes+1), each=2))
  panel_grid_x <- data.frame(y = rep(seq(from = 0.5, to = n_interactions+0.5, by = 1), each=2),
                             x = c(ncelltypes+1, 0), group = rep(1:(n_interactions+1), each=2))
  
  theme_args <- list(panel.grid.major = element_blank(),
                     legend.box = "horizontal",
                     panel.background = element_blank())
  
  theme_args[names(list(...))] <- list(...)
  
  # Check if legend.title is in the extra arguments
  # Multiply by ratio of 5/14: see https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  scale_legend_title_size <- ifelse("legend.title" %in% names(theme_args), theme_args$legend.title$size*(5/14), GeomLabel$default_aes$size)
  # Check if legend.text is in the extra arguments
  scale_legend_text_size <- ifelse("legend.text" %in% names(theme_args), theme_args$legend.text$size*(5/14), GeomLabel$default_aes$size)
  
  # Check if legend.justification is in the extra arguments
  if (!"legend.justification" %in% names(theme_args)) {
    theme_args$legend.justification <- c(1, 0.7)
  }
  
  # Check if legend.position is in the extra arguments
  if (!"legend.position" %in% names(theme_args)){
    theme_args$legend.position <- c(1, 0.7)
  }
  
  p1 <- ggplot() +
    # Draw ligand semicircle
    geom_arc_bar(data = filtered_table %>% filter(type=="ligand",  .data[[rank_filter_col]] <= top_n),
                 aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size)*scale,
                     start = start, end = start + pi, fill=color),
                 color = "white") +
    scale_fill_gradient(low = ligand_fill_colors[1] , high=ligand_fill_colors[2] ,
                        limits=color_lims, oob=scales::squish,
                        n.breaks = 3,
                        guide = guide_colorbar(order = 1),
                        name=paste0(color_title, " (", color_ext[1], ")") %>% stringr::str_wrap(width=15)) +
    # Create new fill scale for receptor semicircles
    new_scale_fill() +
    geom_arc_bar(data = filtered_table %>% filter(type=="receptor", .data[[rank_filter_col]] <= top_n),
                 aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size)*scale,
                     start = start, end = start + pi, fill=color),
                 color = "white") +
    # Size legend
    geom_arc_bar(data = legend2_df, aes(x0=x, y0=y, r0=0, r=sqrt(values)*scale, start=start, end=start+pi), fill="black") +
    geom_rect(data = legend2_df, aes(xmin=x-0.5, xmax=x+0.5, ymin=y-0.5, ymax=y+0.5), color="gray90", fill=NA) +
    geom_text(data = legend2_df, aes(label=values, x=x, y=y-0.6), vjust=1, size = scale_legend_text_size) +
    geom_text(data = data.frame(x = (ncelltypes+4), y = floor(n_interactions/3)+1,
                                label = size_title %>% stringr::str_wrap(width=15)),
              aes(x=x, y=y, label=label), size = scale_legend_title_size, vjust=0, lineheight = .75) +
    # Panel grid
    geom_line(data = panel_grid_y, aes(x=x, y=y, group=group), color = "gray90") +
    geom_line(data = panel_grid_x, aes(x=x, y=y, group=group), color = "gray90") +
    # Draw box to represent x and y "axis"
    geom_rect(data = axis_rect, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), color = "black", fill = "transparent") +
    # Other plot information
    scale_fill_gradient(low = receptor_fill_colors[1], high=receptor_fill_colors[2] , limits=color_lims, oob=scales::squish, n.breaks = 3,
                        name=paste0(color_title,  " (", color_ext[2], ")") %>% stringr::str_wrap(width=15),
                        guide = guide_colorbar(order = 2)) +
    scale_y_continuous(breaks=n_interactions:1, labels=names(lr_interaction_vec), expand = expansion(add=c(0,0))) +
    scale_x_continuous(breaks=1:ncelltypes, labels=levels(filtered_table$sender), position="top", expand = expansion(add=c(0,0))) +
    xlab("Sender cell types") + ylab("Ligand-receptor interaction") +
    coord_cartesian() +
    do.call(theme, theme_args)
  
  # Add unranked ligand and receptor semicircles if requested
  if (show_all_datapoints){
    
    # Limits will depend on true_color_range
    unranked_ligand_lims <- c(0,1); unranked_receptor_lims <- c(0,1)
    if (true_color_range){
      # Follow limits of the top_n lr pairs
      unranked_ligand_lims <- filtered_table %>% filter(type=="ligand",  .data[[rank_filter_col]] <= top_n) %>%
        select(color) %>% range
      unranked_receptor_lims <- filtered_table %>% filter(type=="receptor",  .data[[rank_filter_col]] <= top_n) %>%
        select(color) %>% range
    }
    
    p1 <- p1 + new_scale_fill() +
      geom_arc_bar(data = filtered_table %>% filter(type=="ligand", .data[[rank_filter_col]] > top_n),
                   aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size)*scale,
                       start = start, end = start + pi, fill=color),
                   color = "white") +
      scale_fill_gradient(low = unranked_ligand_fill_colors[1], high=unranked_ligand_fill_colors[2],
                          limits=unranked_ligand_lims, oob = scales::oob_squish,
                          guide = "none") +
      new_scale_fill() +
      geom_arc_bar(data = filtered_table %>% filter(type=="receptor", .data[[rank_filter_col]] > top_n),
                   aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size)*scale,
                       start = start, end = start + pi, fill=color),
                   color = "white") +
      scale_fill_gradient(low=unranked_receptor_fill_colors[1], high=unranked_receptor_fill_colors[2],
                          limits=unranked_receptor_lims, oob = scales::oob_squish,
                          guide = "none")
  }
  
  # Add ranking numbers if requested
  if (show_rankings){
    p1 <- p1 + geom_shadowtext(data = filtered_table %>% filter(.data[[rank_filter_col]] <= top_n),
                               aes(x=x, y=y, label=prioritization_rank))
  }
  
  p1
}