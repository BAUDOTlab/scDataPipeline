checkDirHierarchy <- function() {
    # Check if mandatory directory already exist. If not, create them
    if (!dir.exists(PATH_RDS_OBJECTS)) { 
        print(paste("Create the rds objects directory: ", PATH_RDS_OBJECTS))
        dir.create(PATH_RDS_OBJECTS)
    }
    
    if (!dir.exists(PATH_ROOT)) { 
        print(paste("Create the root directory: ", PATH_ROOT))
        dir.create(PATH_ROOT)
    }
    
    if (!dir.exists(PATH_OUT_HTML)) { 
        print(paste("Create the html output files directory: ", PATH_OUT_HTML))
        dir.create(PATH_OUT_HTML)
    }
    
    if (!dir.exists(PATH_OUT_FIG)) { 
        print(paste("Create the figures directory: ", PATH_OUT_FIG))
        dir.create(PATH_OUT_FIG, recursive = TRUE)
    }
}
