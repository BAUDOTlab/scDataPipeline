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


# Get sorted groups for manual filter
check_empty_groups <- function(df, groups){
	empty_groups <- setdiff(levels(df), groups)
	if(length(empty_groups) != 0){
		stop(paste0("Error: the following groups are empty: ",paste(empty_groups, collapse=", ", sep=", "), ". Change your group thresholds."))
	}
}

# Function to extract top features
extract_top_features <- function(df, topn = 20) {
	# if "group1" is in the column names, sort by group1, then sort by padj
	if ("cluster" %in% names(df)) {
		markers <- df %>%
			dplyr::arrange(cluster, desc(log_fc))
		topn <- markers %>%
			group_by(cluster) %>%
			top_n(n = topn, wt = log_fc)
		topnMarkers <- topn %>%
			dplyr::arrange(cluster, desc(log_fc))
	} else { # else sort by padj
		markers <- df %>%
			dplyr::arrange(padj)
		topn <- markers %>%
			top_n(n = topn, wt = log_fc)
		topnMarkers <- topn %>%
			dplyr::arrange(log_fc)
	}
	return(topnMarkers)
}

merge_features <- function(dataset_list) {
	library(data.table)

	outfile.path <- file.path(PATH_ROOT, "00_rawdata", paste0("combined_", DATASET, "_features.tsv.gz"))
	if(file.exists(outfile.path)) return()

	features.list <- lapply(dataset_list, function(dataset) {
		config <- parseToml(paste0(PATH_REQUIREMENTS, "globalParameters_", dataset,".toml"))
		features <- as.data.frame(fread(file.path(PATH_ROOT, config$path$input$PATH_INPUT_LABDATA, "features.tsv.gz"), header = FALSE, sep = "\t", col.names = c("ENSid", "GeneName", "Type")))
		
		return(features)
	})

	combined.features <- unique(Reduce(function(x, y) rbind(x, y, all=TRUE), features.list))
	
	outfile <- gzfile(outfile.path, "w")
	write.table(combined.features, outfile, col.names=FALSE, sep = "\t", row.names=FALSE)
	close(outfile)
}