######## Function to create a violin plot with threshold lines
create_violin_plot <- function(data, feature, threshold_low	= NULL, threshold_high = NULL, title, subtitle) {
	plot <- VlnPlot(data, features = feature) +
	ggtitle(title, subtitle = subtitle) +
	theme(axis.title.x = element_blank(),
		axis.text = element_text(size = 16),
		axis.text.x = element_text(angle = 0, size = 24, hjust = 0.5, face = "bold"),
		text = element_text(size = 24),
		plot.subtitle = element_text(hjust = 0.5)) +
	NoLegend()

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



######## Function to generate highlighted outlier plot
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

  generate_outlier_table <- function(outlier_type) {
	merged_df <- metadata %>%
	full_join(clusters, by = "cellIDs")
# Calculate the count of outliers and normal cells for each cluster
	cluster_summary <- merged_df %>%
	  group_by(seurat_clusters) %>%
	  summarize(
		outliers = sum( !!sym(outlier_type) == TRUE),
		normal = as.character(sum( !!sym(outlier_type) == FALSE)),
		total_cells = n()
	  ) %>%
	  mutate(percentage_outliers = format(round((outliers / total_cells) * 100,2), nsmall=2),
			outliers = as.character(outliers),
			total_cells = as.character(total_cells),
			seurat_clusters = as.character(seurat_clusters))
	  
	  
	  cluster_summary %>%
		knitr::kable(escape = FALSE) %>%
		  kable_styling(bootstrap_options = c("striped", "hover"),
				  full_width = FALSE)
  }

######## Highlight cells from a specific cluster on a dimreduc plot
highlightClusterPlot <- function(clusterName, seuratObject, reduction = "umap") {
  clusterCells = which( Idents( seuratObject) == clusterName)
  # Create a dimreduc plot with current cluster highlighted (useful for large number of clusters)
  print( DimPlot( seuratObject, 
				  reduction = reduction, 
				  cols="#44444422", 
				  cells.highlight = clusterCells, 
				  cols.highlight = "#FF000088", 
				  sizes.highlight = Seurat:::AutoPointSize(seuratObject)*1.5, 
				  order = clusterCells,  # Plot highlighted cells last
				  group.by=NULL) + 
		   ggtitle( paste(clusterName)) +
		   theme( axis.title.x = element_blank(),
				  axis.title.y = element_blank(),
				  plot.title = element_text( face = "bold",
											 size = rel( 16/14), 
											 hjust = 0.5,
											 vjust = 1, 
											 margin = margin( b = 7)),
				  legend.position = "none"))

}

######## Helper functions for dynamic display, better used in apply functions
set_thresh <- function(thresholds, metadata, metric){
  group_names <- paste0("gp", 1:(length(thresholds)+1), "_", metric)

  # cut() breaks down a data frame column according to the submitted thresholds (infinites are needed to avoid setting minimum and maximum thresholds), 
  # then returns a column with the labels of each groups
  return(cut(metadata, include.lowest=T, breaks=c(-Inf,thresholds, Inf), labels = group_names))
}

threshline <- function(thresh){
  if(thresh != 0){
	geom_hline(aes(yintercept = thresh), color = "#FF0000AA", size = 1.5)
  }
}

threshtext <- function(thresh, i, metric, thresholds){
  geom_text(aes(x=0.5, y=thresh, vjust=-1, label=paste0("gp",i,"_",metric, " - ", thresh)))
}

threshlabel <- function(thresh){
  res = c()
  for(i in 2:length(thresh)){
	if(i==2){
	  res <- append(res,paste0("x < ", thresh[i]))
	} 
	
	if (i==length(thresh)) {
	  res <- append(res, paste0(thresh[i], " < x"))
	} else {
	  res <- append(res, paste0(thresh[i], " < x < ", thresh[i+1]))
	}
  }
  return(res)
}


######## Function to create a dot plot for 40 markers, then will create a new dotPlot
create_dot_plot <- function(SO, geneList, title=NULL) {
  myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
  
  nbMarkers <- dim(geneList)[1]
  if (nbMarkers %% 40 == 0) {
	DotPlot(SO, features = geneList$ENSid) +
	  RotatedAxis() +
	  CenterTitle() +
	  ggtitle(paste0(title, "\n", DATASET)) +
	  scale_colour_gradientn(colours = myPalette(100)) +
	  scale_x_discrete(labels = geneList$GeneName)
  } else {
	nbPlot <- (nbMarkers %/% 40) + 1
	for (i in seq(1, nbPlot)) {
	  print(DotPlot(SO, features = geneList$ENSid[seq(1+(i-1)*(ceiling(nbMarkers/nbPlot)), min(i*ceiling(nbMarkers/nbPlot), nbMarkers))]) +
		RotatedAxis() +
		CenterTitle() +
		ggtitle(paste0(title, " - ", i, "/", (nbMarkers %/% 40) + 1, "\n", DATASET)) +
		scale_colour_gradientn(colours = myPalette(100)) +
		scale_x_discrete(labels = geneList$GeneName[seq(1+(i-1)*(ceiling(nbMarkers/nbPlot)), min(i*ceiling(nbMarkers/nbPlot), nbMarkers))]))
	}
  }
}


######## DEFAULT THEME
defaultTheme <- function() {
	theme_classic() +
	theme(plot.title = element_text(hjust = 0.5, size = 24),
		  plot.subtitle = element_text(hjust = 0.5, size = 20),
		  axis.title = element_text(size = 20),
		  axis.text = element_blank(),
		  line = element_blank())
}

######## BLANK THEME
blankTheme <- function() {
	theme_void() +
	defaultTheme()
}

######## GREY THEME
greyTheme <- function() {
	defaultTheme() +
	theme(panel.background = element_rect(fill="grey95"))
}


######## Create a violinplot showing feature (gene or metadata) values for each cluster (+ number of 'zero' and 'not zero' cells)
violinFeatureByCluster = function(currentFeature, seuratObject, slot = "counts", clustersColor = NULL, yLabel = "Counts", addStats = TRUE, addTitle = TRUE, trimTitle = FALSE) 
{
	if(length( currentFeature)>1) warning( "Several features in 'currentFeature' argument, corresponding values will be averaged...");
	
	# Extract feature values from Seurat object (average if several features given)
	featureValues = apply( FetchData( object = seuratObject, vars = currentFeature, slot = slot), 1, mean);

	# Plot expression values of current annotation for each cluster as violin + jitter
	ggFigure = ggplot( data.frame( Counts = featureValues, Cluster = Idents( seuratObject)), 
					   aes( x = Cluster, y = Counts)) +
		geom_jitter( width = 0.2, height = 0.3, size = 0.4, colour = "#44444444") +
		geom_violin( aes( col = Cluster, fill = Cluster), scale = "width", alpha = 0.4, draw_quantiles = 0.5) +
		theme_minimal() +
		ylab( label = yLabel) +
		theme( legend.position = "None");
	
	# Change ggplot default clusters colors ('clustersColor' should be a named vector of colors)
	if(!is.null(clustersColor))
	{
		ggFigure = ggFigure + 
			scale_color_manual( values = clustersColor) +
			scale_fill_manual( values = clustersColor);
	}
	
	if(addStats)
	{
		# Get the number of cells in which score is zero (or not) for this annotation
		statsZero = tapply( featureValues, 
							Idents( seuratObject), 
							function(x) { paste( table( factor( as.logical( x), levels = c( TRUE, FALSE))), collapse="\n") });
		
		# Create the text layer with stats to be added
		ggFigure = ggFigure + geom_text( data = data.frame( Stat = statsZero, Cluster = names( statsZero)), 
										 aes( x = Cluster, y = max( featureValues)*1.1, label = Stat), 
										 size = 2.5) +
			geom_text( data = data.frame( text = "Pos.\nNull"), 
					   aes( label = text), 
					   x = -0.4, y = max( featureValues)*1.1, hjust = 0, size = 2.5) +
			coord_cartesian( clip = 'off', xlim = c( 0, length(statsZero))); # Allow drawing outside plot region
	} 

    if(addTitle)
    {
		# Remove 'trimTitle' number of ending character(s) (module names have a number suffix added automatically)
		title = if(trimTitle) substr( corr_table[corr_table$ENSid == currentFeature, "GeneName"], 1, nchar( currentFeature)-trimTitle) else corr_table[corr_table$ENSid == currentFeature, "GeneName"];
        		
		# Create the title layer to be added (copy 'FeaturePlot' title style from cowplot)
		ggFigure = ggFigure + ggtitle( label = title) + 
			theme( plot.title =  element_text( face = "bold",
											   size = rel( 16/14), 
											   hjust = 0.5,
											   vjust = 1, 
											   margin = margin( b = 7)));
	}
	
	# Render figure
	#print( ggFigure)
	suppressWarnings(print( ggFigure));
}



