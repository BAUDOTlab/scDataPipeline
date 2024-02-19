load_parameters <- function(config_file) {
    # Parse the TOML file
    paramsList <- parseToml(config_file)
    
    # Retrieve the values from the TOML file directly accessible
    list2env(paramsList, envir = .GlobalEnv)
    
    # Parse the options in the order of the model file
    # ================================================
    # parse the paths variables

    # If the environment variable SCDP_PATH_ROOT exists, we are in a docker env and use it as the root
    if (Sys.getenv("SCDP_PATH_ROOT") != "") {
        assign("PATH_ROOT", Sys.getenv("SCDP_PATH_ROOT"), envir = .GlobalEnv)
    } else {
        assign("PATH_ROOT", path$PATH_ROOT, envir = .GlobalEnv)
    }

    # Function to check for the presence of input and output paths, and then
    # assign their full paths to global environment variables
    checkAndAssignPath <- function(pathName, pathValue) {
        if (pathName %in% names(path$input) || pathName %in% names(path$output)) {
            fullPath <- file.path(PATH_ROOT, pathValue)
            assign(pathName, fullPath, envir = .GlobalEnv)
            cat(paste("Realpath of", pathName, ":", normalizePath(fullPath), "\n"))
            return(TRUE)
        } else {
            cat(paste("The", pathName, "is not provided\n"))
            return(FALSE)
        }
    }

    # Check and assign specific paths
    checkAndAssignPath("PATH_INPUT_LABDATA", path$input$PATH_INPUT_LABDATA)
    checkAndAssignPath("PATH_ATLAS", path$input$PATH_ATLAS)
    checkAndAssignPath("PATH_ATLAS_FILE", path$input$PATH_ATLAS_FILE)
    checkAndAssignPath("PATH_RDS_OBJECTS", path$output$PATH_RDS_OBJECTS)
    checkAndAssignPath("PATH_OUT_HTML", path$output$PATH_OUT_HTML)
    checkAndAssignPath("PATH_OUT_FIG", file.path(path$output$PATH_OUT_FIG, DATASET))
    checkAndAssignPath("PATH_GENES_OF_INTEREST", path$input$PATH_GENES_OF_INTEREST)
    if (exists("PATH_MANUAL_ANNOTATION")) {
        checkAndAssignPath("PATH_MANUAL_ANNOTATION", path$input$PATH_MANUAL_ANNOTATION)
    }

    # parse the qc variables
    if (exists("qc") && is.list(qc)) {
        list2env(qc, envir = .GlobalEnv)
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

