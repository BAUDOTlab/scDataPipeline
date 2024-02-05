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

    assign("PATH_INPUT_LABDATA",	file.path(PATH_ROOT, path$input$PATH_INPUT_LABDATA),		envir = .GlobalEnv)
	assign("PATH_ATLAS",         	file.path(PATH_ROOT, path$input$PATH_ATLAS), 				envir = .GlobalEnv)
    assign("PATH_ATLAS_FILE",       file.path(PATH_ATLAS, path$input$PATH_ATLAS_FILE), 		    envir = .GlobalEnv)
    assign("PATH_RDS_OBJECTS",  	file.path(PATH_ROOT, path$output$PATH_RDS_OBJECTS), 		envir = .GlobalEnv)
    assign("PATH_OUT_HTML",     	file.path(PATH_ROOT, path$output$PATH_OUT_HTML), 			envir = .GlobalEnv)
    assign("PATH_OUT_FIG",      	file.path(PATH_ROOT, path$output$PATH_OUT_FIG, DATASET),	envir = .GlobalEnv)
    assign("PATH_GENES_OF_INTEREST",file.path(PATH_ROOT, path$input$PATH_GENES_OF_INTEREST),	envir = .GlobalEnv)
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

