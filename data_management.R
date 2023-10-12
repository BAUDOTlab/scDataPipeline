# Function to export metadata
export_metadata <- function(meta_data, file_name) {
	# Add rownames in the first column & get names of MD to export
	meta_data <- cbind(cellIDs = rownames(meta_data), meta_data)
	addedMD <- names(meta_data)[!names(meta_data) %in% defaultMDnames]

	# Export added MD for the dataset
	write.table(
		meta_data[addedMD],
		file.path(PATH_OUT_ANALYSIS, file_name),
		row.names = FALSE,
		col.names = TRUE,
		sep = ","
		)
}



# Function to export object as RDS
export_rds <- function(SO, file_name) {
	# Reduce to the default MD & save the object
	SO@meta.data <- SO@meta.data[, defaultMDnames]
	saveRDS(
		SO,
		file.path(PATH_RDS_OBJECTS, file_name)
	)
}


# Function to export dimension reduction tables
export_dimred <- function(SO, dr, file_name) {
    coordinates <- as.data.frame(Embeddings(SO, reduction = dr))
    write.table(x = coordinates,
        file = file.path(
            PATH_OUT_ANALYSIS,
			file_name),
        sep = ",",
        row.names = TRUE,
        col.names = TRUE
    )
}


# Function to extract top features
extract_top_features <- function(df, topn = 20) {
	# if "cluster" is in the column names, sort by cluster, then sort by p_val_adj
	if ("cluster" %in% names(markers)) {
		markers <- markers %>%
			dplyr::arrange(cluster, p_val_adj)
		topn <- markers %>%
			group_by(cluster) %>%
			top_n(n = topn, wt = avg_log2FC)
		topnMarkers <- topn %>%
			dplyr::arrange(cluster, p_val_adj)
	} else { # else sort by p_val_adj
		markers <- markers %>%
			dplyr::arrange(p_val_adj)
		topn <- markers %>%
			top_n(n = topn, wt = avg_log2FC)
		topnMarkers <- topn %>%
			dplyr::arrange(p_val_adj)
	}
	return(topnMarkers)
}
