createH5 <- function(data, file, dataset_name) {
  if (file.exists(file)) {
    unlink(file, recursive = FALSE)
  }
  h5createFile(file)
  h5write(data, file, dataset_name)
  H5close()
}


createMetaH5 <- function(counts_dir){
  collections <- list.dirs(counts_dir, full.names = FALSE)
  collections <- collections[-1]
  for (collection in collections) {
    destdir <- paste0(counts_dir, '/', collection)

    meta <- data.table()
    h5_files <- list.files(destdir, "\\.h5$", full.names = FALSE)
    if (!length(h5_files)) {
      next
    }
    h5_meta <- fread(file.path(destdir, "meta.txt"), index = "file_name")
    #h5_meta <- h5_meta[h5_meta$file_name %in% h5_files, ]
    filename <- paste0(collection, '.h5')
    createH5(h5_meta, filename, 'meta')
  }
}



createPriorityH5 <- function(counts_dir, force = FALSE, verbose = FALSE){
  if (!dir.exists(counts_dir)) {
    message(paste0('Counts directory ', counts_dir, " does not extist" ))
    return()
  }
  h5_files <- list.files(file.path(counts_dir), "\\.h5$", full.names = TRUE, recursive = TRUE)
  list_dirs <-  list.dirs(counts_dir, full.names = FALSE, recursive = TRUE)
  list_dirs <- c(".", list_dirs)
  priority_file <- file.path(counts_dir, "counts_priority.txt")
  need_create <- TRUE
  if (file.exists(priority_file)) {
    priority <- fread(priority_file)
    if (!(setequal(priority$directory,list_dirs) && length(unique(priority$priority)) == length(priority$priority))) {
      message(paste0("!!! Priority file ", priority_file , " is invalid and will be replaced"))
    } else {
      need_create <- FALSE
    }
  }
  if (need_create) {
    priority <- data.table(directory = list_dirs, priority = seq_along(list_dirs))
    write.table(x = priority, file = priority_file, sep = "\t", eol = "\n", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  createH5(priority, 'priority.h5', 'priority')

}


updateIndexH5 <- function(counts_dir, force = FALSE, verbose = FALSE){
  if (!dir.exists(counts_dir)) {
    message(paste0('Counts directory ', counts_dir, " does not extist" ))
    return()
  }
  meta_name <- file.path(counts_dir, "meta.rda")
  h5_files <- list.files(file.path(counts_dir), "\\.h5$", full.names = TRUE, recursive = TRUE)
  if (!length(h5_files)) {
    return()
  }
  if (!force) {
    meta_time <- as.numeric(file.mtime(meta_name))
    h5_mtime <- max(unlist(lapply(h5_files, file.mtime)))
    dirs_mtime <- lapply(file.path(counts_dir, list_dirs[-1]), file.mtime)
        if (length(dirs_mtime) > 0) {
      dir_mtime <- max(unlist(dirs_mtime))
    } else {
      dir_mtime <- -Inf
    }

    if (file.exists(meta_name) && meta_time > h5_mtime && meta_time > dir_mtime) {
      return()
    }
  }
  if (file.exists(meta_name)) {
    unlink(meta_name)
  }
  DT_counts_meta <- data.table(matrix(ncol = 3, nrow = 0, dimnames = list( NULL, c("accession", "file", "collection_type"))))
  for (cur_dir in list_dirs) {
    dir_path <- file.path(counts_dir, cur_dir)
    dir_path <- sub(pattern = "/./?$", replacement = "", x =  dir_path)
    cur_files <- list.files(path = dir_path, pattern = "\\.h5", recursive = FALSE )
    if (length(cur_files) == 0) {
      next
    }
    if (!file.exists(file.path(dir_path, "meta.txt"))) {
      if (startsWith(x = tolower(basename(dir_path)), prefix =  "archs4")) {
        updateARCHS4meta(archDir = dir_path)
      } else if (startsWith(x = tolower(basename(dir_path)), prefix = "dee2")) {
        updateDEE2meta(destDir = dir_path)
      }
    }
    message(paste0('Populating ', cur_dir , ' counts meta' ))
    if (!phantasus:::validateCountsCollection(collectionDir = dir_path, verbose = verbose)) {
      message(paste0("!! files in ", cur_dir , " are ignored because there is not correct meta file in this directory."))
      next
    }
    DT_part <- getCountsMetaPart(counts_dir = counts_dir, collection_name = cur_dir, verbose = verbose)
    if (length(DT_part)) {
      DT_counts_meta <- rbindlist(l = list(DT_counts_meta, DT_part))
    }
    rm(DT_part)

  }
  DT_counts_meta$chunk <- gsmtochunk(DT_counts_meta$accession)
  DT_counts_meta_split <- split(DT_counts_meta, DT_counts_meta$chunk)
  createIndexH5(DT_counts_meta_split, 'index.h5')
  save(DT_counts_meta, file = meta_name, eval.promises = TRUE)
  rm(DT_counts_meta)
}



getCountsMetaPart <- function(counts_dir, collection_name, verbose){
  destdir <- file.path(counts_dir, collection_name)
  if (!dir.exists(destdir)) {
    return()
  }
  DT_h5_meta <- data.table()
  h5_files <- list.files(destdir, "\\.h5$", full.names = FALSE)
  if (!length(h5_files)) {
    return()
  }
  h5_meta <- fread(file.path(destdir, "meta.txt"), index = "file_name")
  for (input_file in h5_files) {
    if (input_file %in% h5_meta$file_name) {
      full_name <- file.path(destdir, input_file)
      relative_path <- file.path(collection_name, input_file )
      h5f <- H5Fopen(full_name, flags = "H5F_ACC_RDONLY")
      accession <- h5read(h5f, h5_meta[file_name == input_file, ]$sample_id)
      h5_part <- data.table(accession = accession,
                            file = relative_path,
                            collection_type = collection_name, indexes = seq_along(accession))
      H5Fclose(h5f)
      DT_h5_meta <- rbindlist(l = list(DT_h5_meta, h5_part))
    } else {
      if (verbose) {
        message(paste0("!! ", file.path(destdir, input_file), " is ignored"))
      }
    }

  }
  return(DT_h5_meta)
}

createIndexH5 <- function(data, file) {
  h5createFile(file)
  names <- names(DT_counts_meta_new_splited)
  for (i in seq_along(names)) {
    h5write(data[[i]], file, paste0("/",names[i]))
  }
}
