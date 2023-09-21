# Function to check and create directories if they don't exist
checkDirHierarchy <- function() {
    createDirectory(PATH_RDS_OBJECTS, "Create the rds objects directory: ")
    createDirectory(PATH_ROOT, "Create the root directory: ")
    createDirectory(PATH_OUT_HTML, "Create the html output files directory: ")
    createDirectory(PATH_OUT_FIG, "Create the figures directory: ", recursive = TRUE)
}


# Function to create a directory if it doesn't exist
createDirectory <- function(path, message, recursive = FALSE) {
    if (!dir.exists(path)) { 
        print(paste(message, path))
        dir.create(path, recursive = recursive)
    }
}
