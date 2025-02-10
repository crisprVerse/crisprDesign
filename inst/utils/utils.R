utils::globalVariables(c("essential_genes_hart2014_human",
                         "essential_genes_hart2014_mouse",
                         "essential_genes_achilles20Q1_human",
                         "essential_genes_achilles20Q1_mouse",
                         "nonessential_genes_hart2014_human",
                         "nonessential_genes_hart2014_mouse",
                         "nonessential_genes_olfactory_human",
                         "nonessential_genes_olfactory_mouse",
                         "igis_human",
                         "igis_mouse"))


#' @importFrom magrittr %>%
#' @importFrom stringr str_extract
#' @importFrom methods is
NULL






#' Provides a list of the available control gene sets
#'
#' @param species character describing species; options are 'human' or 'mouse'
#'
#' @return list storing 'essentials' and 'nonessentials'
#' @export
listAvailableControlGenesets <- function(species=c("human", "mouse")
){
    species <- match.arg(species)
    if (species=="human"){
        essentials <- c("hart2014", "achilles20Q1")
        nonessentials <- c("hart2014", "olfactory")
    } else {
        essentials <- c("hart2014","achilles20Q1")
        nonessentials <- c("hart2014","olfactory")
    }
    return(list(essentials=essentials,
                nonessentials=nonessentials))
}






#' @title Get a list of essential genes from a particular data source
#' @description Get a list of essential genes from a particular data source.
#' 
#' @param source Data source.
#' @param species Either human or mouse.
#' 
#' @return A character vector of gene symbols. 
#' 
#' @examples 
#' getEssentials(source="hart2014", species="mouse")
#' @export
getEssentials <- function(source=c("hart2014", "achilles20Q1"), 
                          species=c("human","mouse")
){
    species <- match.arg(species)
    source <- match.arg(source)
    if (!source %in% listAvailableControlGenesets(species=species)$essentials){
        stop("The indicated set of essential genes is not available.")
    }
    if (source=="hart2014"){
        if (species=="human"){
            x <- essential_genes_hart2014_human
        } else {
            x <- essential_genes_hart2014_mouse
        }
    }
    if (source=="achilles20Q1"){
        if (species=="human"){
            x <- essential_genes_achilles20Q1_human
        } else {
            x <- essential_genes_achilles20Q1_mouse
        }
    }
    return(x)
}





#' @title Get a list of non-essential genes from a particular data source
#' @description Get a list of non-essential genes from a particular data source.
#' 
#' @param source Data source.
#' @param species Either human or mouse.
#' 
#' @return A character vector of gene symbols. 
#' 
#' @examples 
#' getNonessentials(source="hart2014", species="mouse")
#' @export
getNonessentials <- function(source=c("hart2014", "olfactory"),
                             species=c("human","mouse")
){
    species <- match.arg(species)
    source  <- match.arg(source)
    if (!source %in% listAvailableControlGenesets(species=species)$nonessentials){
        stop("The indicated set of non-essential genes is not available.")
    }
    if (source=="hart2014"){
      if (species=="human"){
          x <- nonessential_genes_hart2014_human
      } else {
          x <- nonessential_genes_hart2014_mouse
      }
    }
    if (source=="olfactory"){
        if (species=="human"){
            x <- nonessential_genes_olfactory_human
        } else {
            x <- nonessential_genes_olfactory_mouse
        }
    }
  	return(x)
}




#' @importFrom stringr str_extract
.inferSpecies <- function(se){
    if ("lib" %in% names(metadata(se))){
        species <- str_extract(tolower(metadata(se)$lib),"human|mouse") 
    } else if ("GDBAssayMetadata" %in% names(metadata(se))){
        temp <- metadata(se)[["GDBAssayMetadata"]]
        lib <- temp@annot_id$annot.id
        species <- str_extract(tolower(lib),"human|mouse") 
    } else {
        stop("lib is not part of the metadata. Cannot retrieve species.")
    }
    if (!species %in% c("human", "mouse")){
        stop("inferred species is neither human or mouse")
    }
    return(species)
}







#Transform a string to have only the first letter as uppercase
.simpleCap <- function(x) {
    x <- tolower(x)
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", 
        collapse = " ")
}


#Get gene name from guide name:
.getGeneName <- function(guides, ann) {
    gene.col <- intersect(colnames(ann), c("gene_symbol", "gene.symbol"))
    if(!all(guides %in% ann$id)){
        stop("Not all guides are found in the annotation")
    }
    ann[[gene.col]][match(guides, ann$id)]
}


# returns a named vector of default plotly-themed hex color codes for plotting
# the first 10 elements are plotly colors, remaining are chosen by the author
.getPlotlyColors <- function(var, outline=FALSE){
  var <- unique(var)
  if (outline){       # darker hues
    cols <- c(
      '#083c5a', '#804007', '#165016', '#6b1414', '#4a345f',
      '#462b25', '#723c61', '#404040', '#5e5f11', '#0c5f68',
      '#800080', '#212121', '#337b33', '#76760f', '#80664d', '#646464')
  } else {
    cols <- c(
      '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
      '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
      '#ff00ff', '#424242', '#66f566', '#ebeb1e', '#ffcc99', '#c8c8c8')
  }
  ind <- seq_len(length(var)) %% length(cols)
  ind[ind==0] <- length(cols)
  cols <- cols[ind]
  names(cols) <- var
  return(cols)
}

.makePositive <- function(Y){
    wh0 <- which(Y == 0)
    a <- min(Y)
    if (a < 0) {
        Y <- Y + abs(a)
    }
    Y[wh0] <- 0
    return(Y)
}


# aggregates matrix values factored by row
.aggregateByRow <- function(mat,
                            fact,
                            fun=c('median', 'mean', 'sum')
){
    names <- colnames(mat)
    fun <- match.arg(fun)
    if (fun=='mean'){
        fun <- function(x) colMeans(x)
    } else if (fun=="median"){
        fun <- function(x) colMedians(x)
    } else if (fun=="sum"){
        fun <- function(x) colSums(x)
    }
    # split matrix rows by factor (fact)
    temp <- split(mat, f=fact)
    temp <- lapply(temp, function(x){
        matrix(x, ncol=ncol(mat))
    })
    # apply function, condense list
    temp <- lapply(temp, fun)
    temp <- do.call(rbind, temp)
    colnames(temp) <- names
    return(temp)
}

#' @title Interactive scatter plot using plotly (simple wrapper)
#' 
#' @description Interactive scatter plot using plotly (simple wrapper).
#'     Used mainly by JP for quick plots.
#' 
#' @param x Numeric vector of x values
#' @param y Numeric vector of y values
#' @param labels Data labels
#' @param xlab String for x-axis label
#' @param ylab String for y-axis label
#' @param title String for plot title
#' @param col Numeric vector of length of 1, or equal to the length of x, specifying points color
#' @param pal Character vector specifying color names
#' 
#' @return An interactive plot as a side value.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
#' @importFrom plotly plot_ly layout
#' @importFrom magrittr %>%
iplot <- function(x,
                  y,
                  labels,
                  xlab="",
                  ylab="",
                  title="",
                  col=NULL,
                  pal=NULL
){
    temp <- data.frame(xx=x,
                       yy=y,
                       label=labels)
    zeroline <- TRUE

    if (is.null(col) & is.null(pal)){
        p <- plot_ly(data=temp,
                     x=~xx,
                     y=~yy,
                     text=temp$label)
    }
    if (!is.null(col) & is.null(pal)){
    p <- plot_ly(data=temp,
                 x=~xx,
                 y=~yy, 
                 text=temp$label,
                 color=col)
    }
    if (!is.null(col) & !is.null(pal)){
        p <- plot_ly(data=temp,
                     x=~xx,
                     y=~yy,
                     text=temp$label,
                     color=col,
                     colors=pal)
    }
    p <- p %>%
        layout(title = title,
               xaxis = list(zeroline = zeroline, title=xlab),
               yaxis = list(zeroline = zeroline, title=ylab))
    return(p)
}


getZPrimeFactor <- function(x,y){
    top <- sd(x)+sd(y)
    bottom <- abs(mean(x)-mean(y))
    factor <- 1-3*top/bottom
    factor
}


#' @export
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom magrittr %>%
iplot <- function(x,y,labels,
 xlab="", ylab="",
 title="", col=NULL,
 xlim=c(-10,10),
 ylim=c(0,40)
 ){
  temp <- data.frame(xx=x, yy=y, label=labels)
  zeroline <- TRUE
  type = "scatter"
  mode="markers"
  pal <- c("grey75","firebrick2")

  p <- plot_ly(data = temp, x = ~xx, y = ~yy, 
             text=temp$label,
             type=type,
             mode=mode,
             color=col, colors=pal,
             source="volcano",
             selected=list(marker=list(color='blue'))
             #marker=list(color~col)
  )
  
  p <- p %>%
    plotly::layout(title = title,
           xaxis = list(zeroline = zeroline, title=xlab, range=xlim),
           yaxis = list(zeroline = zeroline, title=ylab, range=ylim)
           )
  p
}

getOverlap <- function(ranks1, ranks2){
    kk <- 1:length(ranks1)
    r <- factor(pmax(ranks1, ranks2), levels=kk)
    tab <- table(r)
    cs <- cumsum(tab)
    cs/kk
}





#' Quick function to extract a subset of reads from NGS fastq files
#' 
#' Quick function to extract a subset of reads from NGS fastq files.
#'
#' @param ngs String specifying NGS project number (eg. "ngs3096").
#'     Reads from the sample specified by \code{sample_index} will be 
#'     read into R. If NULL, \code{file} must be provided instead. 
#' @param file String specifying path to a FASTQ file from which reads will
#'     be read from. Must be specified if \code{ngs} is NULL. 
#' @param samid String specifying the SAMID of the sample to be considered.
#'     Must be part of the NGS project specified by \code{ngs}.
#' @param n_reads Integer specifying number of reads that should be returned.
#'     1000 by default.
#' @param sample_index Integer specifying which sample should be considered
#'     from the list of files available for a given project specified by the 
#'     \code{ngs}. 
#' @param read_mate String specifying if reads from R1 or R2 files should
#'     be returned. Must be either "R1" or "R2". R1" by default. 
#'
#' @return Character vector of sequencing reads.
#'
#' @author Jean-Philippe Fortin
#' @examples
#' \dontrun{
#'     reads <- extractReads(ngs="ngs3565")
#' }
#' @export
#' @importFrom gneDB ngsprojFastq
#' @importFrom readr read_lines
extractReads <- function(ngs=NULL,
                         file=NULL,
                         samid=NULL, 
                         n_reads=1000, 
                         sample_index=1,
                         read_mate=c("R1", "R2")
){

    read_mate <- match.arg(read_mate)
   

    .readsFromFilePartial <- function(file,
                                      n_reads=1000){
        tempFile  <- paste0(tempfile(), ".gz")
        tempFile2 <- gsub(".gz","", tempFile)
        tempFile3 <- gsub(".gz",".txt", tempFile)
        file.copy(file, tempFile, overwrite=TRUE)
        #system(paste0("gunzip -f ", tempFile))
        cmd <- paste0("zcat ",
                      tempFile,
                      " | head -",
                      .makeLongInteger(n_reads*4),
                      " > ",tempFile3)
        system(cmd)
        reads  <- readr::read_lines(tempFile3,
                                    n_max = .makeLongInteger(n_reads*4))
        wh <- seq(2,length(reads),4)
        wh <- wh[wh<length(reads)]
        reads <- reads[wh]
        return(reads)
    }
    
    .readsFromFileComplete <- function(file){
        tempFile <- paste0(tempfile(), ".gz")
        tempFile2 <- gsub(".gz","", tempFile)
        tempFile3 <- gsub(".gz",".txt", tempFile)
        file.copy(file, tempFile, overwrite=TRUE)
        #system(paste0("gunzip -f ", tempFile))
        cmd <- paste0("gunzip -f ", tempFile2)
        system(cmd)
        reads  <- readr::read_lines(tempFile2)
        wh <- seq(2,length(reads),4)
        wh <- wh[wh<=length(reads)]
        reads <- reads[wh]
        return(reads)
    }
    if (is.null(file) & !is.null(ngs)){
        fastq.df = ngsprojFastq(ngs,collapse=FALSE)
        if (read_mate=="R1"){
            fastqs <- fastq.df$READ1_FILE
        } else {
            fastqs <- fastq.df$READ2_FILE
        }
        if (!is.null(samid)){
            fastqs <- fastqs[grepl(samid, fastqs)]
            file <- fastqs[1]
        }  else {
            file <- fastqs[sample_index]  
        }
    } else if (is.null(file) & is.null(ngs)){
        stop("file or ngs must be provided. ")
    }
    if (!file.exists(file)){
        stop("File does not exist.")
    }
    if (!is.null(n_reads)){
        reads <- .readsFromFilePartial(file,
                                       n_reads=n_reads)
    } else {
        reads <- .readsFromFileComplete(file)
    }
    return(reads)
}


#' Quick function to extract barcodes from a list of reads
#' 
#' Quick function to extract barcodes from a list of reads for testing 
#'     purposes. 
#'
#' @param reads Character vector of sequencing reads.
#' @param flank5 String containing the constant sequence on the 5'
#'     flank of the barcode region.
#' @param flank3 String containing the constant sequence on the 3'
#'     flank of the barcode region.
#' @param leftOnly Should only flank5 be considered? 
#'     FALSE by default.
#' @param barcode_len Length of the barcodes to be retrieved.
#' @param ignore_barcode_len Should any barcode length be returned?
#'     FALSE by default. 
#'
#' @return Character vector of barcodes.
#'
#' @author Jean-Philippe Fortin
#' @examples
#' \dontrun{
#'     reads <- extractReads(ngs="ngs3565")
#'     barcodes <- extractBarcodes(reads, flank5="TACCG", flank3="")
#' }
#' @export
#' @importFrom stringr str_extract
extractBarcodes <- function(reads, 
                            flank5="AGTTCG", 
                            flank3="TTCGGACTGT",
                            leftOnly=FALSE, 
                            barcode_len=19,
                            ignore_barcode_len=FALSE
){
    # Only keeping reads that have flanking sequences:
    good1 <- grepl(flank5, reads)
    good2 <- grepl(flank3, reads)
    if (leftOnly){
        good <- which(good1)
    } else {
        good <- which(good1 & good2)
    }
    reads <- reads[good]
    nleft <- nchar(flank5)
    if (!ignore_barcode_len){
        barcodes <- str_extract(reads,
                                paste0(flank5, "[A-Z]+"))
        barcodes <- gsub(flank5,"",barcodes)
        barcodes <- substr(barcodes, 1, barcode_len)
    } else {
        barcodes <- str_extract(reads,
                                paste0(flank5, "[A-Z]+", flank3))
        barcodes <- gsub(flank5,"",barcodes)
        barcodes <- gsub(flank3,"",barcodes)
    }
    return(barcodes)
}





#' @export 
#' @importFrom biomaRt useMart getBM
getEnsemblOrthologs <- function(ids, 
    species.from="mouse", 
    species.to="rat"
){
    species.latin <- list()
    species.latin[["mouse"]] <- "mmusculus"
    species.latin[["human"]] <- "hsapiens"
    species.latin[["rat"]]   <- "rnorvegicus"
    symbols <- list()
    symbols[["mouse"]] <- "mgi_symbol"
    symbols[["human"]] <- "hgnc_symbol"
    symbols[["rat"]]   <- "rgd_symbol"

    marts <- lapply(species.latin, function(x){
        useMart("ensembl", dataset=paste0(x, "_gene_ensembl"))
    })
    names(marts) <- names(species.latin)
    filters <- "ensembl_gene_id"
    # Making attributes:
    x1 <- paste0(species.latin[[species.to]], "_homolog_orthology_type")
    x2 <- paste0(species.latin[[species.to]], "_homolog_ensembl_gene")
    x3 <- "ensembl_gene_id"
    attributes <- c(x1,x2,x3)

    #Getting orthologs:
    df <- getBM(attributes=attributes,
        filters=filters,
        values=ids,
        mart=marts[[species.from]]
    )
    df <- df[df[[x1]]!="",,drop=FALSE]
    col.from <- paste0(species.latin[[species.from]], "_ensembl_gene_id")
    col.to <- paste0(species.latin[[species.to]], "_ensembl_gene_id")
    colnames(df) <- c("type", col.to, col.from)
    df <- df[,c(3,2,1)]
    # Let's add gene symbol:
    df.from <- getBM(attributes=c(symbols[[species.from]], "ensembl_gene_id"),
        filters="ensembl_gene_id",
        values=df[[col.from]],
        mart=marts[[species.from]]
    )
    df.to <- getBM(attributes=c(symbols[[species.to]], "ensembl_gene_id"),
        filters="ensembl_gene_id",
        values=df[[col.to]],
        mart=marts[[species.to]]
    )
    wh.from <- match(df[[col.from]], df.from$ensembl_gene_id)
    wh.to   <- match(df[[col.to]], df.to$ensembl_gene_id)
    df[[paste0(species.latin[[species.from]], "_gene")]] <- df.from[wh.from,1]
    df[[paste0(species.latin[[species.to]], "_gene")]]   <- df.to[wh.to,1]
    return(df)
}


#' @export 
#' @importFrom biomaRt useMart getBM
getSymbolFromEnsembl <- function(ids, 
    species=c("human", "mouse")
){
    species <- match.arg(species)
    species.latin <- list()
    species.latin[["mouse"]] <- "mmusculus"
    species.latin[["human"]] <- "hsapiens"
    species.latin[["rat"]]   <- "rnorvegicus"
    symbols <- list()
    symbols[["mouse"]] <- "mgi_symbol"
    symbols[["human"]] <- "hgnc_symbol"
    symbols[["rat"]]   <- "rgd_symbol"

    marts <- lapply(species.latin, function(x){
        useMart("ensembl", dataset=paste0(x, "_gene_ensembl"))
    })
    names(marts) <- names(species.latin)
    
    df <- getBM(attributes=c(symbols[[species]], "ensembl_gene_id"),
        filters="ensembl_gene_id",
        values=ids,
        mart=marts[[species]]
    )
    return(df)
}


#' @export 
#' @importFrom biomaRt useMart getBM
getEnsemblFromSymbol <- function(ids, 
    species=c("human", "mouse")
){
    species <- match.arg(species)
    species.latin <- list()
    species.latin[["mouse"]] <- "mmusculus"
    species.latin[["human"]] <- "hsapiens"
    species.latin[["rat"]]   <- "rnorvegicus"
    symbols <- list()
    symbols[["mouse"]] <- "mgi_symbol"
    symbols[["human"]] <- "hgnc_symbol"
    symbols[["rat"]]   <- "rgd_symbol"

    marts <- lapply(species.latin, function(x){
        useMart("ensembl", dataset=paste0(x, "_gene_ensembl"))
    })
    names(marts) <- names(species.latin)
    
    df <- getBM(attributes=c(symbols[[species]], "ensembl_gene_id"),
        filters=symbols[[species]],
        values=ids,
        mart=marts[[species]]
    )
    df
}


#' @export
getCellSurvival <- function(r, d){
    .getCellSurvival <- function(r,d){
        if (d!=0){
            result <- ((1+d)-sqrt((1+d)^2-4*d*r))/(2*d)
        } else {
            result <- r
        }
        result
    }
    survivals <- sapply(1:length(r), function(i){
        .getCellSurvival(r[i], d[i])
    })
    survivals
}




