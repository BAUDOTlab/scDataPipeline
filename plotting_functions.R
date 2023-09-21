# Function to create a violin plot with threshold lines
create_violin_plot <- function(data, feature, threshold_low	= NULL, threshold_high = NULL, title, subtitle) {
	plot <- VlnPlot(data, features = feature) +
	ggtitle(title, subtitle = subtitle) +
	theme(axis.title.x = element_blank(),
		axis.text = element_text(size = 16),
		axis.text.x = element_text(angle = 0, size = 24, hjust = 0.5, face = "bold"),
		text = element_text(size = 24),
		plot.subtitle = element_text(hjust = 0.5)) +
	NoLegend() +
	scale_y_continuous(labels = scales::percent_format(scale = 1))

	# Add threshold lines if threshold_low is not NULL
	if (!is.null(threshold_low)) {
		plot <- plot + geom_hline(aes(yintercept = threshold_low), color = "#FF0000AA", size = 1.5)
	}

	# Add threshold lines if threshold_high is not NULL
	if (!is.null(threshold_high)) {
		plot <- plot + geom_hline(aes(yintercept = threshold_high), color = "#FF0000AA", size = 1.5)
	}

	return(plot)
}



# Function to generate highlighted outlier plot
generate_outlier_plot <- function(outlier_type) {
  cellsData <- data.frame(umapCoord, metadata[[outlier_type]])
  colnames(cellsData) <- c(colnames(umapCoord), "Outliers")
  
  ggplot(cellsData[c('UMAP_1', 'UMAP_2')], aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(alpha = .5, size = .75, color = "grey") +
    geom_point(data = cellsData[cellsData$Outliers == TRUE, ],
               aes(color = Outliers), alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("TRUE" = "#FF0000AA", "FALSE" = "#44444444")) +
	# NoLegend() +
    theme_nothing() +
    ggtitle(paste0(DATASET, "\n", table(cellsData$Outliers)["TRUE"], " ", outlier_type, " outlier cells")) +
	CenterTitle() +
    theme(plot.title = element_text(size = 30),
          axis.text = element_blank(),
          line = element_blank(),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 18)) +
    guides(color = guide_legend(override.aes = list(size = 4)))
}



# DEFAULT THEME
defaultTheme <- function() {
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 24),
          plot.subtitle = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 20),
          axis.text = element_blank(),
          line = element_blank(),
          panel.background = element_rect(fill="grey95"))
}