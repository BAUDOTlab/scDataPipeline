load_parameters <- function(config_file) {
    # Parse the TOML file
    paramsList <- parseToml(config_file)
    
    # Retrieve the values from the TOML file directly accessible
    list2env(paramsList, envir = .GlobalEnv)
    
    # Parse the options in the order of the model file
    # ================================================
    
    # parse the paths variables
    assign("PATH_ROOT", path$PATH_ROOT, envir = .GlobalEnv)

    assign("PATH_INPUT_LABDATA",	file.path(PATH_ROOT, path$input$PATH_INPUT_LABDATA),		envir = .GlobalEnv)
	assign("PATH_ATLAS",         	file.path(PATH_ROOT, path$input$PATH_ATLAS), 				envir = .GlobalEnv)
    assign("PATH_RDS_OBJECTS",  	file.path(PATH_ROOT, path$output$PATH_RDS_OBJECTS), 		envir = .GlobalEnv)
    assign("PATH_OUT_HTML",     	file.path(PATH_ROOT, path$output$PATH_OUT_HTML), 			envir = .GlobalEnv)
    assign("PATH_OUT_FIG",      	file.path(PATH_ROOT, path$output$PATH_OUT_FIG, DATASET),	envir = .GlobalEnv)

    # parse the qc variables
    list2env(qc, envir = .GlobalEnv)

    # parse the process variables
    list2env(process, envir = .GlobalEnv)

    # parse the filters variables
    list2env(filters, envir = .GlobalEnv)

    # parse the ctrl++ variables
    list2env(ctrl, envir = .GlobalEnv)


    rm(paramsList)
}
