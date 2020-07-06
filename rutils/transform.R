to_onehot_matrix <- function(str_vec) {
  as.matrix(mltools::one_hot(data.table::data.table(as.factor(str_vec))))
}

real_is_integer <- function(x) {
  if (length(x) < 1L) {
    return(TRUE)
  }
  if (any(is.infinite(x)) || any(is.nan(x))) {
    return(FALSE)
  }
  all(floor(x) == x)
}

stan_rdump <- function(list, file = "", append = FALSE, envir = parent.frame(),
                       width = options("width")$width, quiet = FALSE) {
  if (is.character(file)) {
    ex <- sapply(list, exists, envir = envir)
    if (!all(ex)) {
      notfound_list <- list[!ex]
      if (!quiet) {
        warning(paste("objects not found: ", paste(notfound_list,
          collapse = ", "
        ), sep = ""))
      }
    }
    list <- list[ex]
    if (!any(ex)) {
      return(invisible(character()))
    }
    if (nzchar(file)) {
      file <- file(file, ifelse(append, "a", "w"))
      on.exit(close(file), add = TRUE)
    }
    else {
      file <- stdout()
    }
  }
  ## for (x in list) {
  ##   if (!is_legal_stan_vname(x) & !quiet) {
  ##     warning(paste("variable name ", x, " is not allowed in Stan",
  ##       sep = ""
  ##     ))
  ##   }
  ## }
  l2 <- NULL
  addnlpat <- stringr::str_c("(.{1,", width, "})(\\s|$)")
  for (v in list) {
    vv <- get(v, envir)
    if (is.data.frame(vv)) {
      vv <- data.matrix(vv)
    }
    else if (is.list(vv)) {
      vv <- data_list2array(vv)
    }
    else if (is.logical(vv)) {
      mode(vv) <- "integer"
    }
    else if (is.factor(vv)) {
      vv <- as.integer(vv)
    }
    if (!is.numeric(vv)) {
      if (!quiet) {
        warning(stringr::str_c("variable ", v, " is not supported for dumping."))
      }
      next
    }
    if (!is.integer(vv) && max(abs(vv)) < .Machine$integer.max &&
      real_is_integer(vv)) {
      storage.mode(vv) <- "integer"
    }
    if (is.vector(vv)) {
      if (length(vv) == 0) {
        cat(v, " <- integer(0)\n", file = file, sep = "")
        next
      }
      if (length(vv) == 1) {
        cat(v, " <- ", as.character(vv), "\n",
          file = file,
          sep = ""
        )
        next
      }
      str <- stringr::str_c(
        v, " <- \nc(", stringr::str_c(vv, collapse = ", "),
        ")"
      )
      str <- gsub(addnlpat, "\\1\n", str)
      cat(str, file = file)
      l2 <- c(l2, v)
      next
    }
    if (is.matrix(vv) || is.array(vv)) {
      l2 <- c(l2, v)
      vvdim <- dim(vv)
      cat(v, " <- \n", file = file, sep = "")
      if (length(vv) == 0) {
        str <- stringr::str_c("structure(integer(0), ")
      }
      else {
        str <- stringr::str_c("structure(c(", stringr::str_c(as.vector(vv),
          collapse = ", "
        ), "),")
      }
      str <- gsub(addnlpat, "\\1\n", str)
      cat(str, ".Dim = c(", stringr::str_c(vvdim, collapse = ", "),
        "))\n",
        file = file, sep = ""
      )
      next
    }
  }
  invisible(l2)
}


quickdump <- function(name, myenv = parent.frame()) {
  stan_rdump(
    c(
      "N", "K", "J", "G", "XCond", "XInd", "S", "Xcg", "B", "P"
    ),
    file = name,
    envir = myenv,
    append = TRUE
  )
}
