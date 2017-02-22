\name{BFplot2files}
\alias{BFplot2files}
\title{ Plot multiple DC gene pairs into a file in pdf format}
\description{
The function takes multiple DC gene pairs from a specified file as input, and for each pair, it plot one figure to show the gene expression patterns of the gene pair and merge these figures into a pdf file.
}
\usage{
BFplot2files(dataExp,
       class,
       classlabel,
       edgefilename,
       plotfilename="Plotedges.pdf"
       )
}
\arguments{
   \item{dataExp}{ a data frame or matrix containing expression data. Columns correspond to genes and rows to samples. Such as in \code{\link{Compute_bf}}.}

  \item{class}{ a numeric or character vector contains the corresponding class information of dataExp. By far, package only accepts binary classes. }

  \item{classlabel}{ a numeric or character vector contains only two values. The first used to label class 1 and the second used to label class 2. The order of these two elements must be kept the same throughout the analyses. }


  \item{edgefilename}{an input file giving the DC gene pairs which are required to be plotted. Consists of at least two columns, the first column giving ids (or indexes) for gene1 and the second column giving ids (or indexes) for gene2. gene1 and gene2 are end-nodes which form one DC gene pair. Columns are splited by 'Tab'.}

  \item{plotfilename}{specify a file (in pdf format) to output the multiple plots for all the gene pairs. The default file is "Plotedges.pdf". }
  }
\details{

  The function takes multiple DC gene pairs from a file specified by \code{edgefilename} as input, and for each pair, it plot one figure to show the gene expression patterns of the gene pair. Refer to \link{BFplot} for details.

}

\keyword{ package }

