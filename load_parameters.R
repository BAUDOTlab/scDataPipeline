load_parameters <- function(config_file) {
    print(config_file)
  configs <- readLines(config_file)

  # Get non-empty and not commented lines
  list_filtered <- grep("^[^#\\s]", configs, value = TRUE)

  # seperate every lines according to the = sign
  list_filtered <- strsplit(list_filtered, "\\s*=\\s*")

  for (index in seq(1, length(list_filtered))) {
    # get left and rigth part of the lines
    key_str <- list_filtered[[index]][1]
    var_value <- list_filtered[[index]][2]

    # for composed line, we split it around comma
    if (grepl(" , ", var_value)) {
      splitted <- strsplit(var_value, "\\s*,\\s*")

      for (part in splitted[[1]]) {
        # check if any part of the line is a precreated variable
        # and replace it by its value
        if (exists(part)) {
          pos <- which(part == splitted[[1]])
          splitted[[1]][pos] <- eval(parse(text = part))
        }
      }

      # collapse all the parts to create a single word value
      var_value <- paste(unlist(splitted), collapse = "")
    }

    # create the variable (left part of the line)
    # with the value (right part of the line) with
    # variable interpretation if needed
    assign(key_str, var_value, envir = parent.frame())
  }
}

