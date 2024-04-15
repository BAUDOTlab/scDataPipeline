assign_parameter <- function(variable_name, value){
    assign(variable_name, value, envir = .GlobalEnv)

    .GlobalEnv[["parameters_list"]][[variable_name]] <- paste(value, collapse = ", ")
}

load_parameters <- function(config_file) {
    # Parse the TOML file
    paramsList <- parseToml(config_file)
    
    # Retrieve the values from the TOML file directly accessible
    list2env(paramsList, envir = .GlobalEnv)
    
    # reassign some variables to avoid messing with the lists
    assign_parameter("DATASET",DATASET)
    if(exists("CONDITION")) assign_parameter("CONDITION", CONDITION)
    assign_parameter("ENS_ID_COLUMN", ENS_ID_COLUMN)

    # Parse the options in the order of the model file
    # ================================================
    # parse the paths variables

    # If the environment variable SCDP_PATH_ROOT exists, we are in a docker env and use it as the root
    if (Sys.getenv("SCDP_PATH_ROOT") != "") {
        assign("PATH_ROOT", Sys.getenv("SCDP_PATH_ROOT"), envir=.GlobalEnv)
    } else {
        assign("PATH_ROOT", path$PATH_ROOT, envir=.GlobalEnv)
    }
    normalizePath(PATH_ROOT, mustWork = TRUE)

    # Function to check for the presence of input and output paths, and then
    # assign their full paths to global environment variables
    assignPath <- function(pathName, pathValue, relativeTo = PATH_ROOT) {
        if (pathName %in% names(path$input) || pathName %in% names(path$output)) {
            # Only use the relative path if the base path doesn't exist
            if(!file.exists(pathValue)){
                fullPath <- file.path(relativeTo, pathValue)
            } else {
                fullPath <- pathValue
            }
            assign(pathName, fullPath, envir=.GlobalEnv)
            #cat(paste("Realpath of", pathName, ":", normalizePath(fullPath), "\n"))
            return(TRUE)
        } else {
            cat(paste("The", pathName, "is not provided\n"))
            return(FALSE)
        }
    }

    # Check and assign specific paths
    assignPath("PATH_INPUT_LABDATA", path$input$PATH_INPUT_LABDATA)
    assignPath("PATH_ATLAS", path$input$PATH_ATLAS)
    assignPath("PATH_ATLAS_FILE", path$input$PATH_ATLAS_FILE, PATH_ATLAS)
    assignPath("PATH_RDS_OBJECTS", path$output$PATH_RDS_OBJECTS)
    assignPath("PATH_OUT_HTML", path$output$PATH_OUT_HTML)
    assignPath("PATH_OUT_FIG", file.path(path$output$PATH_OUT_FIG, DATASET))
    assignPath("PATH_GENES_OF_INTEREST", path$input$PATH_GENES_OF_INTEREST)
    if (exists("PATH_MANUAL_ANNOTATION")) {
        assignPath("PATH_MANUAL_ANNOTATION", path$input$PATH_MANUAL_ANNOTATION)
    }

    # parse the qc variables
    if (exists("qc") && is.list(qc)) {
        list2env(qc, envir=.GlobalEnv)
    }

    # parse the process variables
    if (exists("process") && is.list(process)) {
        list2env(process, envir = .GlobalEnv)
    }

    # parse the filters variables
    if (exists("filters") && is.list(filters)) {
        list2env(filters, envir = .GlobalEnv)
    }

    # parse the ctrl variables
    if (exists("ctrl") && is.list(ctrl)) {
        list2env(ctrl, envir = .GlobalEnv)
    }

    #parse the combine variables
    if (exists("combine") && is.list(combine)) {
        list2env(combine, envir = .GlobalEnv)
    }

    rm(paramsList)
}

