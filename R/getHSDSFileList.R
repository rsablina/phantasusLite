#' Returns list of collections on HSDS-server
#' @param url, containing url of the server and root domain.
#'
#' @return List of collections on the server
#'
#' @export
#' @import rhdf5client
getHSDSCollectionsList <- function(url='https://ctlab.itmo.ru/hsds/?domain=/counts') {
    src <- httr::parse_url(url)
    dir <- src$query$domain
    src <- paste0(src$scheme,'://',src$hostname,'/',src$path)
    src <- HSDSSource(src)
    domains <- listDomains(src, dir)
    domains <- domains[-grep("*\\.h5$", domains)]
    domains <- gsub(paste0(dir, '/'), '', domains)
    return(domains)
}

#' Returns list of HDF5-files of a collection
#' @param url, containing url of the server and root domain.
#' @param HDF5 collection name
#'
#' @return List of HDF5-files of the collection
#'
#' @export
#' @import rhdf5client
getHSDSCollectionFileList <- function(url='https://ctlab.itmo.ru/hsds/?domain=/counts', collection) {
    src <- httr::parse_url(url)
    dir <- src$query$domain
    dir <- paste0(dir, '/', collection)
    src <- paste0(src$scheme,'://',src$hostname,'/',src$path)
    src <- HSDSSource(src)
    request <- paste0(src@endpoint, "/domains?domain=",
                      dir)
    response <- rhdf5client:::submitRequest(request)
    domains <- response[["domains"]]
    listDomains <- list()
    for (domain in domains) {
        listDomains <- append(listDomains, sub(".*/", "", domain$name))
    }
    listDomains <- unlist(listDomains)
    return(listDomains)
}

#' Returns list of all HDF5-files on HSDS-server
#' @param url, containing url of the server and root domain.
#'
#' @return List of all HDF5-files on the server
#'
#' @export
#' @import rhdf5client
getHSDSFileList <- function(url='https://ctlab.itmo.ru/hsds/?domain=/counts') {
    src <- httr::parse_url(url)
    dir <- src$query$domain
    #dir <- paste0(dir, '/', collection)
    src <- paste0(src$scheme,'://',src$hostname,'/',src$path)
    src <- HSDSSource(src)
    request <- paste0(src@endpoint, "/domains?domain=",
                      dir)
    response <- rhdf5client:::submitRequest(request)
    domains <- response[["domains"]]
    collections <- getHSDSCollectionsList(url)
    hdf5FileList  <- list()
    for (collection in collections) {
        files <- getHSDSCollectionFileList(url, collection)
        hdf5FileList <- append(hdf5FileList, files)
    }
    hdf5FileList <- unlist(hdf5FileList)
    return(hdf5FileList)
}
