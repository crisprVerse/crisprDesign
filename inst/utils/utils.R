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







