% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TxDb2GRangesList.R
\name{getTxDb}
\alias{getTxDb}
\title{getTxDb}
\usage{
getTxDb(file = NA, organism = "Homo sapiens", tx_attrib = "gencode_basic", ...)
}
\arguments{
\item{file}{File argument for \code{\link{makeTxDbFromGFF}}
(see help page for \code{\link{makeTxDbFromGFF}}).
If \code{NA} (default), function will return a \linkS4class{TxDb}
object from Ensembl using \code{\link{makeTxDbFromEnsembl}}.}

\item{organism}{String specifying genus and species name
(e.g. "Homo sapiens" for human).}

\item{tx_attrib}{Argument passed on to \code{\link{makeTxDbFromEnsembl}}
when \code{file} is not specified. See help page
for \code{\link{makeTxDbFromEnsembl}}.}

\item{...}{Additional arguments passed on to either
\code{\link{makeTxDbFromGFF}} (if \code{file} is specified) or
\code{\link{makeTxDbFromEnsembl}} if \code{file} is NA.}
}
\value{
A \linkS4class{TxDb} object.
}
\description{
Convenience function for constructing
    a \linkS4class{TxDb} object.
}
\examples{
if (interactive()){
    # To obtain a TxDb for Homo sapiens from Ensembl:
    txdb <- getTxDb()

    # To obtain a TxDb from a GFF file:
    file='ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3'
    txdb <- getTxDb(file=file)
}
}
\author{
Jean-Philippe Fortin, Luke Hoberecht
}
