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
  load(paste(file_path, "/meta.rda", sep = ""))

}
