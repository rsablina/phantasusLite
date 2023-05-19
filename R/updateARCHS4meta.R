updateARCHS4meta <- function(file_path) {
  meta <- fread(paste(file_path, "/meta.txt", sep = ""))
  metafilename <- "meta.h5"
  h5createFile(metafilename)
  h5write(meta, metafilename,"meta")
  counts_priority <- fread(paste(file_path, "/counts_priority.txt", sep = ""))
  h5write(counts_priority, metafilename,"priority")
  H5close()
}

updateARCHS4Index <- function(src, file_path) {
  metaindexfilename <- "index1.h5"
  h5createFile(metaindexfilename)
  DT_counts_meta_new_splited <- split(DT_counts_meta_indexes, DT_counts_meta_indexes$chunk)
  names <- names(DT_counts_meta_new_splited)
  for (i in seq_along(names)) {
    h5write(DT_counts_meta_new_splited[[i]], metaindexfilename, paste0("/",names[i]))
  }

}



