\name{Compute_bf}
\alias{Compute_bf}
\title{ Calculation of Bayes factors }
\description{
 Estimate the stength of pair-wise differential co-expression by calculating Bayes factors between all genes from a given expression matrix. 
}
\usage{
Compute_bf(dataExp, class,classlabel,bfthr=6,echo=FALSE)
}
\arguments{
  \item{dataExp}{ a data frame or matrix containing expression data. Columns correspond to genes and rows to samples. Must be kept the same throughout the analyses.}

  \item{class}{ a numeric or character vector contains the corresponding class information of dataExp. By far, package only accepts binary classes. }
  \item{classlabel}{ a numeric or character vector contains only two values. The first used to label class 1 and the second used to label class 2. the order of these two elements must be kept the same throughout the analyses. }
  \item{bfthr}{ Bayes factor threshold. The default value is 6. If it's \code{NULL}, all the gene pairs are remained.}
  \item{echo}{ a logical variable, indicates whether print processing information or not. If \code{TRUE}, it will print processing information for every 5000 gene pairs. }
}
\details{

  The function estimate the stength of pair-wise differential co-expression by calculating Bayes factors between all genes from a given expression matrix. Gene pairs with Bayes factors higher than \code{bfthr} will be remained, if \code{bfthr} is not \code{NULL}. 

}
\value{
  A \code{bfmatrix} object. It's a data frame containing the following elements:
  \item{geneid1}{A character vector giving gene ids for gene1.}
  \item{geneid2}{A character vector giving gene ids for gene2.}
  \item{bf.value}{A numerical vector giving the values of Bayes factors for gene pair, gene1-gene2.}
  \item{type}{A character vector indicating which model is selected for the gene pair.}

}
\keyword{ package }
