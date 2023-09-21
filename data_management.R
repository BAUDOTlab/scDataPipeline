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


