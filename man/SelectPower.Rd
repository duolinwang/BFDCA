\name{SelectPower}
\alias{SelectPower}
\title{  Choosing the soft-thresholding power by analysing of the network topology }
\description{
 Choose the soft thresholding power based on the criterion of approximate scale-free topology.
}
\usage{
SelectPower(bfsimilarity,powers)
}
  \arguments{
  \item{bfsimilarity}{ an adjacency matrix based on pair-wise Bayes factor generated from \code{\link{BF_similarity}}.}
  
  \item{powers}{ a vector of candidate soft-thresholding powers provided by user, soft-thresholding power will be chosen from it.}
  }
\details{
 This function plots a network topology under a set of soft-thresholding powers (a set of candidate \code{softpower} parameters for \code{\link{BF_WGCNA}}) provided by user. It constructs two panels in the plot, the left panel shows the scale-free fit index (y-axis) as a function of the soft-thresholding power (x-axis). The right panel shows the mean connectivity (degree, y-axis) as a function of the soft-thresholding power (x-axis). 
According to WGCNA [1] (\url{https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf}), they proposed to choose the lowest soft thresholding power and achieve an approximate scale-free topology (scale-free topology fit index reaches 0.9).
}
\value{
  None.
}
\references{ 
[1] B. Zhang and S. Horvath. A general framework for weighted gene co-expression network analysis. Statistical
Applications in Genetics and Molecular Biology, 4(1):Article 17, 2005.
}
\keyword{ package }



