#' @title Add well row and column info to rowData
#' 
#' @description Add well row and column info to rowData
#' 
#' @param se A \linkS4class{SummarizedExperiment} object.
#' 
#' @return A \linkS4class{SummarizedExperiment} object with updated
#'     annotation.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' data(seExample)
#' seExample <- addPlateRowAndColumn(seExample)
#' 
#' @importFrom stringr str_extract
#' @importFrom SummarizedExperiment rowData rowData<-
#' @export
addPlateRowAndColumn <- function(se){
    ann <- rowData(se)
    if (!"Well" %in% colnames(ann)){
        stop("Well must be in colData(se)")
    }
    ann$row    <- str_extract(rowData(se)[["Well"]], "[A-Z]+")
    ann$column <- str_extract(rowData(se)[["Well"]], "[0-9]+")
    ann$column <- as.integer(ann$column)
    rowData(se) <- ann
    return(se)
}






#' @title Add a column in the feature annotation for gene class
#' 
#' @description Add a column in the feature annotation for gene class.
#' 
#' @param se A \linkS4class{SummarizedExperiment} object.
#' @param gene.field String specifying the name of the column
#'     that contains the gene names.
#' @param class.field String specifying the name of the column
#'     that will store gene class. Default is "class".
#' @param default.class What should be the default class for
#'     genes not specified in \code{classes}?
#'     Default is "target". 
#' @param classes Named list specifying the gene classes.
#'     Names of the list correspond to the names of the classes,
#'     each each list element contains the gene names corresponding
#'     the the given class.
#' 
#' @return A \linkS4class{SummarizedExperiment} object with updated
#'     annotation.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @importFrom stringr str_extract
#' @importFrom SummarizedExperiment rowData rowData<-
#' @export
addGeneClassColumn <- function(se,
                               gene.field="Gene",
                               class.field="class",
                               default.class="target",
                               classes=NULL
){
    ann <- rowData(se)
    ann[[class.field]] <- default.class
    if (!gene.field %in% colnames(ann)){
        stop("gene.field must be in colnames(rowData(se))")
    }
    if (!is.null(classes)){
        n_classes <- length(classes)
        for (class in names(classes)){
            wh <- which(ann[[gene.field]] %in% classes[[class]])
            ann[[class.field]][wh] <- class
        }
    }
    rowData(se) <- ann
    return(se)
}




#' @title log2 transformation of assay data
#' 
#' @description log2 transformation of assay data.
#' 
#' @param se A \linkS4class{SummarizedExperiment} object.
#' @param assay Integer vector specifying the indices of the assay
#'     to log2 transform. If NULL (default), all assays
#'     are log2 trasnformed.
#' @param offset Integer specifying offset to use 
#'     in the log2 transformation.
#' 
#' @return A \linkS4class{SummarizedExperiment} object with 
#'     one or many assays log2 transformed.
#' 
#' @examples
#' data(seExample)
#' seExample <- logTransform(seExample, offset=1)
#' 
#' @author Jean-Philippe Fortin
#' 
#' @importFrom stringr str_extract
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom SummarizedExperiment assays assays<-
#' @export
logTransform <- function(se,
                         assay=NULL,
                         offset=1
){
    n_assays <- length(assays(se))
    if (is.null(assay)){
        assaysToTransform <- seq_len(n_assays)
    } else {
        assaysToTransform <- assay
    }
    
    for (assay in assaysToTransform){
        assays(se)[[assay]] <- log2(assays(se)[[assay]]+offset)
    }
    return(se)
}









#' @title Z score transformation using negative controls.
#' 
#' @description Z score transformation using negative controls.
#' 
#' @param se A \linkS4class{SummarizedExperiment} object.
#' @param neg.controls Character vector specifying genes that should
#'     be used for z-score transformation. The genes must be present
#'     in the \code{gene.field} column of \code{rowData(se)}.
#'     If NULL (default), all genes will be considered.
#' @param gene.field String specifying the name of the column
#'     that contains the gene names specified by \code{neg.controls}.
#' @param fun.centering String specifying the method 
#'     to use for centering. Either "median" (default) or "mean".
#' @param fun.scaling String specifying the method
#'     to use for scaling. Either "mad" (default) or "sd".
#' @param scaling.ignore.controls Should the scaling factors be
#'     calculated using all genes? TRUE by default, which
#'     will ignore negative controls.
#' @param class.field String specifying the name of the column
#'     that stored the gene class in \code{rowData(se)}.
#' @param class.levels Gene classes that should be considered
#'     for z-score transformation. If NULL (default), all classes
#'     are considered. If specified, intersection of the specified 
#'     gene classes and the genes specified in \code{neg.controls}
#'     will be considered for z scoring.
#' 
#' @return A \linkS4class{SummarizedExperiment} object with 
#'     assays z-score transformed. 
#' 
#' @author Jean-Philippe Fortin
#' 
#' @importFrom matrixStats colMedians colMads colSds
#' @export
transformToZScores <- function(se,
                               neg.controls=NULL,
                               gene.field="Gene",
                               fun.centering=c("median", "mean"),
                               fun.scaling=c("mad", "sd"),
                               scaling.ignore.controls=TRUE,
                               class.field="class",
                               class.levels="target"
){
    fun.centering <- match.arg(fun.centering)
    fun.scaling   <- match.arg(fun.scaling)
    n_assays <- length(assays(se))
    ann <- rowData(se)
    n_genes <- nrow(se)
    for (assay in seq_len(n_assays)){
        Y <- assays(se)[[assay]]
        wh_centering <- seq_len(n_genes)
        wh_scaling   <- seq_len(n_genes)
        wh_class     <- seq_len(n_genes)
        if (!is.null(neg.controls)){
            wh_centering <- which(ann[[gene.field]] %in% neg.controls)
            if (!scaling.ignore.controls){
                wh_scaling <- which(ann[[gene.field]] %in% neg.controls)
            } 
        }
        if (!is.null(class.levels)){
            wh_class <- which(ann[[class.field]] %in% class.levels) 
        }
        wh_centering <- intersect(wh_centering, wh_class)
        wh_scaling   <- intersect(wh_scaling, wh_class)
            
    
        if (length(wh_centering)==0 | length(wh_scaling)==0 ){
            stop("Combination of neg.controls and class.levels returns
                  0 gene for transformation.")
        }

        #Centering:
        if (fun.centering=="median"){
            factors <- colMedians(Y[wh_centering,,drop=FALSE], na.rm=TRUE)
        } else {
            factors <- colMeans(Y[wh_centering,,drop=FALSE], na.rm=TRUE)
        }
        Y <- sweep(Y, 2 ,factors, "-")
        
        #Scaling:
        if (fun.scaling=="mad"){
            factors <- colMads(Y[wh_scaling,,drop=FALSE], na.rm=TRUE)
        } else {
            factors <- colSds(Y[wh_scaling,,drop=FALSE], na.rm=TRUE)
        }
        Y <- sweep(Y, 2 ,factors, "/")
        assays(se)[[assay]] <- Y
    }
    return(se)
}



#' @title Calculate limma statistics 
#' 
#' @description Calculate limma statistics.
#' 
#' @param se A \linkS4class{SummarizedExperiment} object.
#' @param design Design matrix to pass to limma.
#' @param contrast.matrix Contrast matrix to pass to limma.
#' 
#' @return For each assay in \code{se}, it returns a named
#'     list with the following elements:
#'     - lfc: log-fold changes
#'     - pval: p-lvaues
#'     - fdr: FDR-adjusted p-values. 
#'     
#' 
#' @examples 
#' data(seExample)
#' library(SummarizedExperiment)
#' 
#' classes <- list(neg=c("OR5M9","OR6N1"),
#'     pos=c("KRAS","NRAS"),
#'     rnaimax=c("RNAiMAX"),
#'     empty="empty")
#' seExample <- addGeneClassColumn(seExample, classes=classes)
#' 
#' # Transforming the data:
#' seExample <- logTransform(seExample, offset=1)
#' seExample <- normalizeBetweenPlates(seExample)
#' seExample <- transformToZScores(seExample)
#' 
#' # Creating design matrix:
#' pheno  <- as.factor(colData(seExample)$Condition)
#' pheno  <- relevel(pheno, ref="WT")
#' design <- model.matrix(~pheno-1)
#' colnames(design) <- gsub("pheno", "", colnames(design))
#' contrasts <- c(Mutant = "Mutant",
#'     Wildtype="WT",
#'     Mutant_vs_Wildtype="Mutant - WT")
#' contrast.matrix <- cbind("Wildtype"=c(1,0),
#'     "Mutant"=c(0,1),
#'     "Mutant_vs_Wildtype"=c(-1,1))
#' 
#' # Calculate statistics:
#' stats <- calculateStatistics(seExample,
#'     design=design,
#'     contrast.matrix=contrast.matrix)
#' 
#' @author Jean-Philippe Fortin
#' 
#' @importFrom limma lmFit contrasts.fit eBayes
#' @importFrom stats p.adjust
#' @export
calculateStatistics <- function(se,
                                design,
                                contrast.matrix=NULL
){
    n_assays <- length(assays(se))
    results <- list()
    for (assay in seq_len(n_assays)){
        Y <- assays(se)[[assay]]
        fit <- lmFit(Y, design)
        if (!is.null(contrast.matrix)){
            fit <- contrasts.fit(fit, contrast.matrix)
        }
        fit <- eBayes(fit)
        stats <- list()
        stats[["lfc"]]  <- fit$coefficients
        stats[["pval"]] <- fit$p.value
        stats[["fdr"]]  <- apply(stats[["pval"]], 2, p.adjust, "fdr")
        results[[assay]] <- stats
    }
    names(results) <- names(assays(se))
    return(results)
}



#' @title Between-plate normalization
#' 
#' @description Between-plate normalization.
#' 
#' @param se A \linkS4class{SummarizedExperiment} object.
#' @param neg.controls Character vector specifying genes that should
#'     be used for estimating normalizatin factors.
#'     The genes must be present in the \code{gene.field}
#'     column of \code{rowData(se)}.
#'     If NULL (default), all genes will be considered.
#' @param fun String specifying the method 
#'     to use for normalization. Either "median" (default) or "mean".
#' @param plate.field String specifying the name of the column
#'     that contains the plate ID.
#' @param gene.field String specifying the name of the column
#'     that contains the gene names specified by \code{neg.controls}.
#' @param class.field String specifying the name of the column
#'     that stored the gene class in \code{rowData(se)}.
#' @param class.levels Gene classes that should be considered
#'     for estimating normalization factors. If NULL (default), all classes
#'     are considered. If specified, intersection of the specified 
#'     gene classes and the genes specified in \code{neg.controls}
#'     will be considered for z scoring.
#' 
#' @return A \linkS4class{SummarizedExperiment} object with 
#'     value normalized between plates. 
#' 
#' @examples 
#' data(seExample)
#' seExample <- logTransform(seExample, offset=1)
#' 
#' classes <- list(neg=c("OR5M9","OR6N1"),
#'     pos=c("KRAS","NRAS"),
#'     rnaimax=c("RNAiMAX"),
#'     empty="empty")
#' seExample <- addGeneClassColumn(seExample, classes=classes)
#' 
#' seExample <- normalizeBetweenPlates(seExample)
#' 
#' @author Jean-Philippe Fortin
#' 
#' @importFrom matrixStats colMedians
#' @importFrom SummarizedExperiment rowData assays assays<-
#' @export
normalizeBetweenPlates <- function(se,
                                   neg.controls=NULL,
                                   fun=c("median", "mean"),
                                   plate.field="Plate",
                                   gene.field="Gene",
                                   class.field="class",
                                   class.levels=c("target")
){
    fun <- match.arg(fun)
    n_assays  <- length(assays(se))
    assaysToNormalize <- seq_len(n_assays)
    ann    <- rowData(se)
    plate  <- ann[[plate.field]]
    plates <- unique(plate)
    n_plates <- length(plates)
    n_genes <- nrow(se)
    if (!gene.field %in% colnames(ann)){
        stop("gene.field not found in colnames(rowData(se))")
    }
    if (!class.field %in% colnames(ann)){
        stop("class.field not found in colnames(rowData(se))")
    }
    for (assay in assaysToNormalize){
        Y <- assays(se)[[assay]]
        for (k in seq_len(n_plates)){
            if (is.null(neg.controls)){
                wh_control <- seq_len(n_genes)
            } else {
                wh_control <- which(ann[[gene.field]] %in% neg.controls)
            }
            wh_plate <- which(plate==plates[[k]])
            wh_norm  <- intersect(wh_control, wh_plate)

            # Only normalizing using target genes:
            if (!is.null(class.levels)){
                wh_class <- which(ann[[class.field]] %in% class.levels) 
            } else {
                wh_class <- seq_len(n_genes)
            }
            wh_norm <- intersect(wh_norm, wh_class)
            if (fun=="median"){
                factors <- colMedians(Y[wh_norm,,drop=FALSE], na.rm=TRUE)
            } else {
                factors <- colMeans(Y[wh_norm,,drop=FALSE], na.rm=TRUE)
            }
            Y[wh_plate,] <- sweep(Y[wh_plate,,drop=FALSE], 2 ,factors, "-")
        }
        assays(se)[[assay]] <- Y
    }
    return(se)
}



#' @title Within-plate normalization
#' 
#' @description Within-plate normalization
#' 
#' @param se A \linkS4class{SummarizedExperiment} object.
#' @param fun String specifying the method 
#'     to use for normalization. Either "median" (default) or "mean".
#' @param plate.field String specifying the name of the column
#'     that contains the plate ID.
#' @param what String specifying if plate rows or columns
#'     should be normalized. "row" by default.
#' @param class.field String specifying the name of the column
#'     that stored the gene class in \code{rowData(se)}.
#' @param class.levels Gene classes that should be considered
#'     for estimating normalization factors. If NULL (default), all classes
#'     are considered. If specified, intersection of the specified 
#'     gene classes and the genes specified in \code{neg.controls}
#'     will be considered for z scoring.
#' 
#' @return A \linkS4class{SummarizedExperiment} object with 
#'     value normalized within plates.
#' 
#' #' @examples 
#' data(seExample)
#' seExample <- logTransform(seExample, offset=1)
#' seExample <- normalizeWithinPlate(seExample)
#' 
#' @author Jean-Philippe Fortin
#' @importFrom SummarizedExperiment rowData assays assays<-
#' @importFrom matrixStats colMedians
#' @export
normalizeWithinPlate <- function(se,
                                 fun=c("median", "mean"),
                                 what=c("row", "column"),
                                 plate.field="Plate",
                                 class.field="class",
                                 class.levels=c("target")
){
    fun  <- match.arg(fun)
    what <- match.arg(what)
    n_assays  <- length(assays(se))
    ann <- rowData(se)
    if (!what %in% colnames(ann)){
        se  <- addPlateRowAndColumn(se)
        ann <- rowData(se)
    }
    plate  <- ann[[plate.field]]
    slice  <- paste0(plate, "_", ann[[what]])
    slices <- unique(slice)
    n_slices <- length(slices)
    n_genes <- nrow(se)
    for (assay in seq_len(n_assays)){
        Y <- assays(se)[[assay]]
        for (k in seq_len(n_slices)){
            wh_control <- seq_len(n_genes)
            wh_slice   <- which(slice==slices[[k]])
            wh_norm <- intersect(wh_control, wh_slice)
            # Only normalizing using target genes:
            if (!is.null(class.levels)){
                wh_class <- which(ann[[class.field]] %in% class.levels) 
            } else {
                wh_class <- seq_len(n_genes)
            }
            wh_norm <- intersect(wh_norm, wh_class)
            if (fun=="median"){
                factors <- colMedians(Y[wh_norm,,drop=FALSE], na.rm=TRUE)
            } else {
                factors <- colMeans(Y[wh_norm,,drop=FALSE], na.rm=TRUE)
            }
            Y[wh_slice,] <- sweep(Y[wh_slice,,drop=FALSE], 2 ,factors, "-")
        }
        assays(se)[[assay]] <- Y
    }
    return(se)
}













