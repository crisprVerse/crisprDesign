\name{crisprDesignNews}
\title{crisprDesign News}
\encoding{UTF-8}


\section{Version 1.13.10}{
\itemize{
    \item Added optional argument substitutions to addEditedAlleles"
}
}

\section{Version 1.13.9}{
\itemize{
    \item Fixed condaEnv argument for addOnTargetScores
}
}

\section{Version 1.13.8}{
\itemize{
    \item SNPs can now be annotated from a GRanges object. 
}
}

\section{Version 1.13.6}{
\itemize{
    \item Fixed splicing junction aa annotation when low scores
}
}


\section{Version 1.13.4}{
\itemize{
    \item Fixed the accessor editedAlleles.
}
}

\section{Version 1.13.1}{
\itemize{
    \item Added more annotation columns to the addEditedAlleles function.
}
}


\section{Version 1.7.1}{
\itemize{
    \item Fixed getTssObjectFromTxObject for negative strand.
}
}


\section{Version 1.5.2}{
\itemize{
    \item Added function addReinitiationFlag. 
}
}


\section{Version 1.1.25}{
\itemize{
    \item Added a warning message when having duplicated gRNAs.
}
}

\section{Version 1.1.21}{
\itemize{
    \item Fixed intergenic annotation of alignments when non-standard chromosomes are found.
}
}


\section{Version 1.1.20}{
\itemize{
    \item \code{flattenGuideSet} is now \code{GuideSet2DataFrames}.
}
}

\section{Version 1.1.19}{
\itemize{
    \item Changed argument \code{standard_chr_only} to FALSE by default in the functions \code{getSpacerAlignments} and friends.
}
}


\section{Version 1.1.18}{
\itemize{
    \item Updated the object \code{guideSetExampleWithAlignments} to contain the latest annotations.
}
}

\section{Version 1.1.17}{
\itemize{
    \item Updated the object \code{guideSetExampleFullAnnotation} to contain the latest annotations.
}
}

\section{Version 1.1.12}{
\itemize{
    \item Refactored \code{addGeneAnnotation}.
}
}


\section{Version 1.1.10}{
\itemize{
    \item Added the function \code{getPreMrnaSequence}.
}
}


\section{Version 1.1.8}{
\itemize{
    \item Added features for non-targeting controls (\code{addNtcs}).
}
}

\section{Version 1.0.0}{
\itemize{
    \item New package \pkg{crisprDesign}, for the design of CRISPR gRNA.
}
}
