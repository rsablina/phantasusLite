
#' @export
#' @import rhdf5client
#'
#'
loadCountsHSDS <- function(es, src, file, sample_id = NULL, gene_id = NULL, gene_id_type = NULL) {
  f <- HSDSFile(src, file)

  if (is.null(sample_id) || is.null(gene_id) || is.null(gene_id_type)) {
    metafilepath <- paste(dirname(file),"/newmeta.h5",sep="")
    metaf <- HSDSFile(src, metafilepath)
    metads <- HSDSDataset(metaf, '/meta/metatable')
  }


  dg <- HSDSDataset(f, gene_id)
  genes <- dg[1:dg@shape]

  dsamples <- HSDSDataset(f, samples_id)
  samples <- dsamples[1:dsamples@shape]
  samples <- dsamples[1:200000]

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

loadCountsHSDSPriority <- function(es, src, file) {

}

getGEO <- function(name) {


}


