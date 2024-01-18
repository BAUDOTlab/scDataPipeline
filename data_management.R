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
	# if "group1" is in the column names, sort by group1, then sort by padj
	if ("group1" %in% names(df)) {
		markers <- df %>%
			dplyr::arrange(group1, feature_rank)
		topn <- markers %>%
			group_by(group1) %>%
			top_n(n = topn, wt = feature_rank)
		topnMarkers <- topn %>%
			dplyr::arrange(group1, feature_rank)
	} else { # else sort by padj
		markers <- df %>%
			dplyr::arrange(padj)
		topn <- markers %>%
			top_n(n = topn, wt = feature_rank)
		topnMarkers <- topn %>%
			dplyr::arrange(feature_rank)
	}
	return(topnMarkers)
}
