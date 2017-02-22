\name{ALL}
\alias{ALL}
\docType{data}
\title{ALL dataset.}
\description{
Acute Lymphoblastic Leukemia (ALL) dataset. The ALL dataset was obtained from the ALL package [1] in Bioconductor [2], which uses microarray platform HG-U95Av2 with Affymetrix gene identifiers. Our test focused on the B-cell ALL, which consists of 37 samples with BCR/ABL mutation and 42 samples with no cytogenetic abnormalities. Columns correspond to genes and rows to samples.
}
\usage{data(ALL)}
\format{

A list contains class information and expression.
ALL$class is a vector containing class information, with "1" indicates BCR/ABL mutation and "2" indicates no cytogenetic abnormalities.
ALL$data is a data frame containing expression data. It consist of 79 samples and 8638 genes.

}
\examples{
data(ALL)
}
\references{ 
  {[1] Chiaretti S, Li X, Gentleman R, Vitale A, Vignetti M, Mandelli F, et al. Gene expression profile of adult T-cell acute lymphocytic leukemia identifies distinct subsets of patients with different response to therapy and survival. Blood. 2004;103:2771-8.}

 {[2] Gentleman RC, Carey VJ, Bates DM, Bolstad B, Dettling M, Dudoit S, et al. Bioconductor: open software development for computational biology and bioinformatics. Genome Biol. 2004;5:R80.
 }
}
\keyword{datasets}

