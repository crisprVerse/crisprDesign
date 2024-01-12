#' @title getTxDb
#' @description Convenience function for constructing
#'     a \linkS4class{TxDb} object.
#' 
#' @param file File argument for \code{\link{makeTxDbFromGFF}}
#'     (see help page for \code{\link{makeTxDbFromGFF}}).
#'     If \code{NA} (default), function will return a \linkS4class{TxDb}
#'     object from Ensembl using \code{\link{makeTxDbFromEnsembl}}.
#' @param organism String specifying genus and species name
#'     (e.g. "Homo sapiens" for human). Required if \code{file} is not
#'     provided. If \code{file} is provided, this value can be set to \code{NA}
#'     to have organism information as unspecified.
#' @param release Ensembl release version; passed to
#'     \code{\link{makeTxDbFromEnsembl}} when \code{file} is not specified.
#'     See help page for \code{\link{makeTxDbFromEnsembl}}.
#' @param tx_attrib Argument passed to \code{\link{makeTxDbFromEnsembl}}
#'     when \code{file} is not specified. See help page
#'     for \code{\link{makeTxDbFromEnsembl}}.
#' @param ... Additional arguments passed to either
#'     \code{\link{makeTxDbFromGFF}} (if \code{file} is specified) or
#'     \code{\link{makeTxDbFromEnsembl}} if \code{file} is NA.
#' 
#' @return A \linkS4class{TxDb} object.
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @examples
#' if (interactive()){
#'     # To obtain a TxDb for Homo sapiens from Ensembl:
#'     txdb <- getTxDb()
#' 
#'     # To obtain a TxDb from a GFF file:
#'     file='https://www.mirbase.org/ftp/CURRENT/genomes/hsa.gff3'
#'     txdb <- getTxDb(file=file)
#' }
#' @importFrom GenomicFeatures makeTxDbFromEnsembl
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @export
getTxDb <- function(file=NA,
                    organism,
                    release=NA,
                    tx_attrib='gencode_basic',
                    ...
){
    if (is.null(organism)){
        stop("'organism' cannot be NULL")
    }
    if (is.na(file)){
        if (is.na(organism)){
            stop("'organism' must be specified when 'file' is missing")
        }
        txdb <- GenomicFeatures::makeTxDbFromEnsembl(organism=organism,
                                                     release=release,
                                                     tx_attrib=tx_attrib,
                                                     ...)
    } else {
        txdb <- GenomicFeatures::makeTxDbFromGFF(file,
                                                 organism=organism,
                                                 ...)
    }
    return(txdb)
}




#' @title Convert a \linkS4class{TxDb} object into a \linkS4class{GRangesList}
#' @description Convenience function to reformat a \linkS4class{TxDb} object
#'    into a \linkS4class{GRangesList}.
#' 
#' @param txdb A \linkS4class{TxDb} object.
#' @param standardChromOnly Should only standard chromosomes be kept?
#'     TRUE by default. 
#' @param genome Optional string specifying genome. e.g. "hg38", to be
#'     added to \code{genome(txdb)}.
#' @param seqlevelsStyle String specifying which style should be used
#'     for sequence names. "UCSC" by default (including "chr").
#'     "NCBI" will omit "chr" in the sequence names. 
#' @return A named \linkS4class{GRangesList} of length 7 with the
#'     following elements:
#'     \code{transcripts}, \code{exons}, \code{introns}, \code{cds},
#'     \code{fiveUTRs}, \code{threeUTRs} and \code{tss}.
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @seealso \code{\link{getTxDb}} to obtain a \linkS4class{TxDb} object.
#' 
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @importFrom BiocGenerics start<- end<-
#' @importFrom GenomicFeatures tidyTranscripts tidyExons 
#' @importFrom GenomicFeatures transcriptLengths fiveUTRsByTranscript
#' @importFrom GenomicFeatures threeUTRsByTranscript tidyIntrons
#' @importFrom AnnotationDbi select
#' @importFrom methods as
#' 
#' @examples
#' if (interactive()){
#'     # To obtain a TxDb for Homo sapiens from Ensembl:
#'     txdb <- getTxDb()
#'     
#'     # To convert to a GRanges list:
#'     txdb <- TxDb2GRangesList(txdb)
#' }
#' 
#' 
#' @export
#' @importFrom S4Vectors isTRUEorFALSE
TxDb2GRangesList <- function(txdb,
                             standardChromOnly=TRUE,
                             genome=NULL,
                             seqlevelsStyle=c("UCSC", "NCBI")
){
    txdb <- .validateTxDb(txdb)
    seqlevelsStyle <- match.arg(seqlevelsStyle)
    stopifnot("standardChromOnly must be TRUE or FALSE" = {
        S4Vectors::isTRUEorFALSE(standardChromOnly)
    })
    stopifnot("genome must be a character string or NULL" = {
        is.null(genome) ||
            is.vector(genome, mode="character") &&
            length(genome) == 1
    })
    #message("Converting TxDb to GRangesList...")
    gRangesList <- .TxDb2GRangesList(txdb=txdb,
                                     standardChromOnly=standardChromOnly,
                                     genome=genome)
    gRangesList <- .changeSeqlevelsStyle(gRangesList,
                                         seqlevelsStyle=seqlevelsStyle)
    return(gRangesList)
}





#' @importFrom GenomeInfoDb organism keepStandardChromosomes genome<-
#' @importFrom GenomicRanges GRangesList
#' @importFrom S4Vectors metadata metadata<-
.TxDb2GRangesList <- function(txdb,
                              standardChromOnly,
                              genome
){
    organism <- GenomeInfoDb::organism(txdb)
    
    if (requireNamespace("biomaRt") && !is.na(organism)){
        bm <- .getBiomartData(txdb, organism)
    } else {
        bm <- NULL
    }
    
    exonInfo    <- .getExonInfo(txdb)
    exons       <- .getExonsFromTxDB(txdb, bm, exonInfo)
    fiveUTRs    <- .getFiveUtrsFromTxdb(txdb, bm, exonInfo)
    threeUTRs   <- .getThreeUtrsFromTxdb(txdb, bm, exonInfo)
    transcripts <- .getTranscriptsFromTxDB(txdb, bm)
    introns     <- .getIntronsFromTxDb(txdb, bm)
    cds         <- .getCdsFromExons(exons)
    tss         <- .getTssFromTxDb(txdb, bm, exonInfo)
    #chrominfo   <- .getChromInfoFromTxDb(txdb)
    
    ls <- list(transcripts=transcripts,
               exons=exons,
               cds=cds,
               fiveUTRs=fiveUTRs,
               threeUTRs=threeUTRs,
               introns=introns,
               tss=tss)
    ls <- GenomicRanges::GRangesList(ls)
    if (standardChromOnly && !is.na(organism)){
        ls <- GenomeInfoDb::keepStandardChromosomes(ls,
                                                    species=organism,
                                                    pruning.mode="fine")
    }
    S4Vectors::metadata(ls) <- S4Vectors::metadata(txdb)
    
    if (!is.null(genome)){
        GenomeInfoDb::genome(ls) <- genome
    }

    return(ls)
}




# Get comprehensive gene annotation from Biomart
#' @importFrom S4Vectors metadata
.getBiomartData <- function(txdb,
                            organism
){
    ## use Ensembl version, if applicable
    version <- S4Vectors::metadata(txdb)
    version <- version$value[version$name == "Ensembl release"]
    archives <- biomaRt::listEnsemblArchives()
    url <- archives$url[archives$version == version]
    if (length(url) != 0){
        mart <- biomaRt::useMart("ensembl",
                                 host=url)
    } else {
        mart <- biomaRt::useMart("ensembl")
    }
    
    .inferMartDataset <- function(organism){
        organism <- tolower(organism)
        organism <- strsplit(organism, " ")[[1]]
        genus <- organism[-length(organism)]
        genus <- vapply(genus, substr, start=1, stop=1, FUN.VALUE=character(1))
        genus <- paste0(genus, collapse="")
        species <- organism[length(organism)]
        dataset <- paste0(genus, species, "_gene_ensembl")
        return(dataset)
    }

    martDataset <- .inferMartDataset(organism)
    has_bm_dataset <- martDataset %in% biomaRt::listDatasets(mart)$dataset
    if (!has_bm_dataset){
        stop("Organism '", organism, "' not recognized in biomaRt.")
    } else {
        mart <- biomaRt::useDataset(martDataset,
                                    mart=mart)
        attributes <- c("ensembl_transcript_id",
                        "ensembl_gene_id",
                        "external_gene_name")
        filters <- c("ensembl_transcript_id")
        bm <- biomaRt::getBM(attributes=attributes,
                             filters=filters,
                             values="", # obtain all values
                             mart=mart)
    }
    return(bm)
}



# Extract an transcript GRanges from a TxDb object
.getTranscriptsFromTxDB <- function(txdb, bm=NULL){

    transcripts <- GenomicFeatures::tidyTranscripts(txdb)
    tx_select <- AnnotationDbi::select(txdb,
                                       keys=as.character(transcripts$tx_id),
                                       columns=c('TXID', 'CDSNAME', 'TXTYPE'),
                                       keytype='TXID')
    # remove duplicated rows with NA for CDSNAME
    dup <- tx_select$TXID[duplicated(tx_select$TXID)]
    tx_select <- tx_select[!tx_select$TXID %in% dup|!is.na(tx_select$CDSNAME),]
    transcripts$protein_id <- tx_select$CDSNAME[match(transcripts$tx_id,
                                                      tx_select$TXID)]
    transcripts$tx_type <- tx_select$TXTYPE[match(transcripts$tx_id,
                                                  tx_select$TXID)]
    transcripts$tx_id <- transcripts$tx_name
    transcripts$tx_name <- NULL
    if (!is.null(bm)){
        wh <- match(transcripts$tx_id, bm$ensembl_transcript_id)
        transcripts$gene_symbol <- bm$external_gene_name[wh]
    } else {
        transcripts$gene_symbol <- NA
    }
    return(transcripts)
}


# Extract exon annotation from a TxDb object
.getExonInfo <- function(txdb){
    exons <- GenomicFeatures::tidyExons(txdb)
    exon_cols <- c('EXONID', 'CDSSTART', 'CDSEND', 'CDSNAME', 
                   'TXSTART', 'TXEND', 'TXTYPE', 'TXNAME',
                   'EXONNAME', 'GENEID')
    exonInfo <- AnnotationDbi::select(txdb,
                                      keys=as.character(exons$exon_id), 
                                      columns=exon_cols,
                                      keytype='EXONID')
    return(exonInfo)
}


# Extract an exon GRanges from TxDb object
.getExonsFromTxDB <- function(txdb, bm, exonInfo){
    exons <- GenomicFeatures::tidyExons(txdb)
    exons_exonId_txId <- paste0(exons$exon_id, "_", exons$tx_name)
    exonInfo_exonId_txId <- paste0(exonInfo$EXONID, "_", exonInfo$TXNAME)
    select_match <- match(exons_exonId_txId, exonInfo_exonId_txId)
    
    exons$cds_start  <- exonInfo$CDSSTART[select_match]
    exons$cds_end    <- exonInfo$CDSEND[select_match]
    exons$protein_id <- exonInfo$CDSNAME[select_match]
    exons$tx_start   <- exonInfo$TXSTART[select_match]
    exons$tx_end     <- exonInfo$TXEND[select_match]
    exons$tx_type    <- exonInfo$TXTYPE[select_match]
    # for now, make cds_len refer to entire tx
    cds_len <- GenomicFeatures::transcriptLengths(txdb,
                                                  with.cds_len=TRUE)
    exons$cds_len <- cds_len$cds_len[match(exons$tx_id, cds_len$tx_id)]
    # add gene_symbol using biomaRt
    if (!is.null(bm)){
        wh <- match(exons$gene_id, bm$ensembl_gene_id)
        exons$gene_symbol <- bm$external_gene_name[wh]
    } else {
        exons$gene_symbol <- NA
    }
    # rename some columns
    exons$tx_id <- exons$tx_name
    exons$tx_name <- NULL
    exons$exon_id <- exons$exon_name
    exons$exon_name <- NULL
    return(exons)
}


# Extract CDS from exons GRanges object
.getCdsFromExons <- function(exons){
    cds <- exons
    cds <- cds[!is.na(cds$cds_start) & !is.na(cds$cds_end)]
    cds$exon_start <- start(cds)
    cds$exon_end <- end(cds)
    BiocGenerics::start(cds) <- cds$cds_start
    BiocGenerics::end(cds) <- cds$cds_end
    cds$cds_start <- NULL
    cds$cds_end <- NULL
    return(cds)
}  



# Extract chromosome information from TxDb object
#' @importFrom GenomicFeatures as.list
#' @importFrom GenomeInfoDb seqlevels
.getChromInfoFromTxDb <- function(txdb){
    chrominfo <- GenomicFeatures::as.list(txdb)$chrominfo
    standardChrom <- chrominfo$chrom %in% GenomeInfoDb::seqlevels(txdb)
    chrominfo <- chrominfo[standardChrom, , drop=FALSE]
    ## drop?
    if (!any(grepl('^chr', chrominfo$chrom))){
        chrominfo$chrom <- paste0('chr', chrominfo$chrom)
    }
    ## ^drop?
    return(chrominfo)
}


# Extract 5' UTRs from TxDb object
#' @importFrom GenomicFeatures fiveUTRsByTranscript
.getFiveUtrsFromTxdb <- function(txdb,
                                 bm=NULL,
                                 exonInfo
){
    fiveUTRs <- fiveUTRsByTranscript(txdb, use.names=TRUE)
    fiveUTRs <- unlist(methods::as(fiveUTRs, 'GRangesList'))
    names(fiveUTRs)[is.na(names(fiveUTRs))] <- ""
    fiveUTRs$exon_id <- fiveUTRs$exon_name
    fiveUTRs$exon_name <- NULL
    fiveUTRs$tx_id <- names(fiveUTRs)
    wh <- match(fiveUTRs$tx_id, exonInfo$TXNAME)
    fiveUTRs$gene_id <- exonInfo$GENEID[wh]

    #wh <- match(fiveUTRs$exon_id, exonInfo$EXONNAME)
    #fiveUTRs$tx_id   <- exonInfo$TXNAME[wh]
    #fiveUTRs$gene_id <- exonInfo$GENEID[wh]
    if (!is.null(bm)){
        wh <- match(fiveUTRs$tx_id, bm$ensembl_transcript_id)
        fiveUTRs$gene_symbol <- bm$external_gene_name[wh]
    } else {
        fiveUTRs$gene_symbol <- NA
    }
    return(fiveUTRs)
}


# Extract 3' UTRs from TxDb object
#' @importFrom GenomicFeatures threeUTRsByTranscript
.getThreeUtrsFromTxdb <- function(txdb,
                                  bm=NULL,
                                  exonInfo
){
    threeUTRs <- threeUTRsByTranscript(txdb, use.names=TRUE)
    threeUTRs <- unlist(methods::as(threeUTRs, 'GRangesList'))
    names(threeUTRs)[is.na(names(threeUTRs))] <- ""
    threeUTRs$exon_id <- threeUTRs$exon_name
    threeUTRs$exon_name <- NULL
    threeUTRs$tx_id <- names(threeUTRs)
    wh <- match(threeUTRs$tx_id, exonInfo$TXNAME)
    threeUTRs$gene_id <- exonInfo$GENEID[wh]

    #wh <- match(threeUTRs$exon_id, exonInfo$EXONNAME)
    #threeUTRs$tx_id   <- exonInfo$TXNAME[wh]
    #threeUTRs$gene_id <- exonInfo$GENEID[wh]
    if (!is.null(bm)){
        wh <- match(threeUTRs$tx_id, bm$ensembl_transcript_id)
        threeUTRs$gene_symbol <- bm$external_gene_name[wh]
    } else {
        threeUTRs$gene_symbol <- NA
    }
    return(threeUTRs)
}


# Extract Introns TxDb object
.getIntronsFromTxDb <- function(txdb, bm){
    introns <- tidyIntrons(txdb)
    introns$tx_id <- introns$tx_name
    introns$tx_name <- NULL
    if (!is.null(bm)){
        wh <- match(introns$tx_id, bm$ensembl_transcript_id)
        introns$gene_symbol <- bm$external_gene_name[wh]
    } else {
        introns$gene_symbol <- NA
    }
    return(introns)
}



# Extract TSS information from a TxDb object
.getTssFromTxDb <- function(txdb, bm, exonInfo){
    promoters <- promoters(txdb,
                           downstream=1,
                           upstream=0)
    promoters$tx_id <- promoters$tx_name
    promoters$tx_name <- NULL
    wh <- match(promoters$tx_id, exonInfo$TXNAME)
    promoters$gene_id <- exonInfo$GENEID[wh]
    if (!is.null(bm)){
        wh <- match(promoters$tx_id, bm$ensembl_transcript_id)
        promoters$gene_symbol <- bm$external_gene_name[wh]
    } else {
        promoters$gene_symbol <- NA
    }
    return(promoters)
}



.validateTxDb <- function(txdb){
    if (!is(txdb, "TxDb")){
        stop("txdb must be a TxDb object")
    }
    return(txdb)
}
