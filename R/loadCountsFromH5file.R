
#' @export
#' @import rhdf5client
#'
#'
loadCountsFromH5FileHSDS <- function(es, src, file, sample_id = NULL, gene_id = NULL, gene_id_type = NULL) {
  if (nrow(es) > 0) {
    return(es)
  }
  f <- HSDSFile(src, file)

  if (is.null(sample_id) || is.null(gene_id) || is.null(gene_id_type) || is.null(priority)) {
    metafilepath <- paste(dir,"/meta.h5",sep="")
    metaf <- HSDSFile(src, metafilepath)
    metads <- HSDSDataset(metaf, '/meta')
    metatable <- metads[1:metads@shape]
    if (is.null(sample_id)) {ы
      sample_id <- metatable$sample_id[metatable$file_name == file]
    }
    if (is.null(gene_id)) {
      gene_id <- metatable$gene_id[metatable$file_name == file]
    }
    if (is.null(gene_id_type)) {
      gene_id_type <- metatable$gene_id_type[metatable$file_name == file]
    }
  }

  dg <- HSDSDataset(f, gene_id)
  genes <- dg[1:dg@shape]

  dsamples <- HSDSDataset(f, samples_id)
  cnt <- 1
  samples <- list()
  while (cnt < dsamples@shape) {
    samples <- c(samples, dsamples[cnt:min(c(cnt + 10000, dsamples@shape))])
    cnt <- min(c(cnt + 10001, dsamples@shape))
  }
  datasets <- listDatasets(f)
  if ("info/version" %in% datasets) {
    arch_version <- HSDSDataset(f, '/info/version')
  } else {
    arch_version <- HSDSDataset(f, '/meta/info/version')
  }
  arch_version <- arch_version[1:arch_version@shape[1]]
  if (is.na(arch_version)) {
    arch_version <- 8
  }
  sampleIndexes <- match(es$geo_accession, samples)
  ds <- HSDSDataset(f, '/data/expression')
  smap <- data.frame(sampleIndexes, geo_accession=es$geo_accession)
  smap <- smap[order(smap$sampleIndexes),]
  smap <- smap[!is.na(smap$sampleIndexes),]
  expression <- ds[1:ds@shape[1], smap$sampleIndexes]
  rownames(expression) <- genes
  colnames(expression) <- smap$geo_accession
  es <- es[,es$geo_accession %in% colnames(expression)]
  expression <- expression[,es$geo_accession]
  es2 <- ExpressionSet(assayData = expression,
                           phenoData = phenoData(es[, !is.na(sampleIndexes)]),
                           annotation = annotation(es),
                           experimentData = experimentData(es))

  return(es2)
}

loadCountsFromHSDS <- function(es, src, dir) {
  if (nrow(es) > 0) {
    return(es)
  }

  metafilepath <- paste(dir,"/meta.h5",sep="")
  metaf <- HSDSFile(src, metafilepath)
  priorityds <- HSDSDataset(metaf, '/priority')
  counts_priority <- priorityds[1:priorityds@shape]
  load(paste(countsDir,"meta.rda",sep = "/"))
  sample_amount <- DT_counts_meta[accession %in% es$geo_accession, .(.N), by = list(file, type_fac = factor(x = collection_type, levels = priority))]
  if (nrow(sample_amount) == 0) {
    return(es)
  }
  setorderv(x = sample_amount,cols = c("N","type_fac"),order = c(-1,1))
  destfile <- sample_amount[,.SD[1]]$file

  f <- HSDSFile(src, paste(dir, destfile, sep = "/"))

  metads <- HSDSDataset(metaf, '/meta')
  metatable <- metads[1:metads@shape]
  sample_id <- metatable$sample_id[metatable$file_name == file]
  gene_id <- metatable$gene_id[metatable$file_name == file]
  gene_id_type <- metatable$gene_id_type[metatable$file_name == file]

  dg <- HSDSDataset(f, gene_id)
  genes <- dg[1:dg@shape]

  dsamples <- HSDSDataset(f, samples_id)

  samples <- list()

  datasets <- listDatasets(f)

  if ("info/version" %in% datasets) {
    arch_version <- HSDSDataset(f, '/info/version')
  } else {
    arch_version <- HSDSDataset(f, '/meta/info/version')
  }

  arch_version <- arch_version[1:arch_version@shape[1]]
  if (is.na(arch_version)) {
    arch_version <- 8
  }

  sampleIndexes <- match(es$geo_accession, samples)
  ds <- HSDSDataset(f, '/data/expression')
  smap <- data.frame(sampleIndexes, geo_accession=es$geo_accession)
  smap <- smap[order(smap$sampleIndexes),]
  smap <- smap[!is.na(smap$sampleIndexes),]
  expression <- ds[1:ds@shape[1], smap$sampleIndexes]
  rownames(expression) <- genes
  colnames(expression) <- smap$geo_accession
  es <- es[,es$geo_accession %in% colnames(expression)]
  expression <- expression[,es$geo_accession]

  es2 <- ExpressionSet(assayData = expression,
                       phenoData = phenoData(es[, !is.na(sampleIndexes)]),
                       annotation = annotation(es),
                       experimentData = experimentData(es))

  return(es2)
}






