load_parameters <- function(config_file) {
    # Parse the TOML file
    paramsList <- parseToml(config_file)
    
    # Retrieve the values from the TOML file directly accessible
    list2env(paramsList, envir = .GlobalEnv)
    
    # Parse the options in the order of the model file
    # ================================================
    # parse the paths variables

    # If the environment variable SCDP_PATH_ROOT exists, we are in a docker env and use it as the root
    if(Sys.getenv("SCDP_PATH_ROOT") != ""){
        assign("PATH_ROOT", Sys.getenv("SCDP_PATH_ROOT"), envir = .GlobalEnv)
    } else {
        assign("PATH_ROOT", path$PATH_ROOT, envir = .GlobalEnv)        
    }

    if ("PATH_INPUT_LABDATA" %in% names(path$input)) {
        assign("PATH_INPUT_LABDATA",file.path(PATH_ROOT, path$input$PATH_INPUT_LABDATA),		envir = .GlobalEnv)
        cat("Realpath of PATH_INPUT_LABDATA:", normalizePath(PATH_INPUT_LABDATA), "\n")
    } else {
        cat("The PATH_INPUT_LABDATA is not provided\n")
    }

    if ("PATH_ATLAS" %in% names(path$input)) {
        assign("PATH_ATLAS",file.path(PATH_ROOT, path$input$PATH_ATLAS), 				envir = .GlobalEnv)
        cat("Realpath of PATH_ATLAS:", normalizePath(PATH_ATLAS), "\n")
    } else {
        cat("The PATH_ATLAS is not provided\n")
    }

    if ("PATH_ATLAS_FILE" %in% names(path$input)) {
        assign("PATH_ATLAS_FILE",file.path(PATH_ATLAS, path$input$PATH_ATLAS_FILE), 		    envir = .GlobalEnv)
        cat("Realpath of PATH_ATLAS_FILE:", normalizePath(PATH_ATLAS_FILE), "\n")
    } else {
        cat("The PATH_ATLAS_FILE is not provided\n")
    }

    if ("PATH_RDS_OBJECTS" %in% names(path$output)) {
        assign("PATH_RDS_OBJECTS",file.path(PATH_ROOT, path$output$PATH_RDS_OBJECTS), 		envir = .GlobalEnv)
        cat("Realpath of PATH_RDS_OBJECTS:", normalizePath(PATH_RDS_OBJECTS), "\n")
    } else {
        cat("The PATH_RDS_OBJECTS is not provided\n")
    }

    if ("PATH_OUT_HTML" %in% names(path$output)) {
        assign("PATH_OUT_HTML",file.path(PATH_ROOT, path$output$PATH_OUT_HTML), 			envir = .GlobalEnv)
        cat("Realpath of PATH_OUT_HTML:", normalizePath(PATH_OUT_HTML), "\n")
     } else {
        cat("The PATH_OUT_HTML is not provided\n")
    }

    if ("PATH_OUT_FIG" %in% names(path$output)) {
        assign("PATH_OUT_FIG",file.path(PATH_ROOT, path$output$PATH_OUT_FIG, DATASET),	envir = .GlobalEnv)
        cat("Realpath of PATH_OUT_FIG:", normalizePath(PATH_OUT_FIG), "\n")
     } else {
        cat("The PATH_OUT_FIG is not provided\n")
    }

    if ("PATH_GENES_OF_INTEREST" %in% names(path$input)) {
        assign("PATH_GENES_OF_INTEREST",file.path(PATH_ROOT, path$input$PATH_GENES_OF_INTEREST),	envir = .GlobalEnv)
        cat("Realpath of PATH_GENES_OF_INTEREST:", normalizePath(PATH_GENES_OF_INTEREST), "\n")
     } else {
        cat("The PATH_GENES_OF_INTEREST is not provided\n")
    }
    
    if (exists("PATH_MANUAL_ANNOTATION")) {
        assign("PATH_MANUAL_ANNOTATION",file.path(PATH_ROOT, path$input$PATH_MANUAL_ANNOTATION),	envir = .GlobalEnv)
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

