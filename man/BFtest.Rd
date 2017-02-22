\name{BFtest}
\alias{BFtest}
\title{ Predict case/control conditions or phenotypes for testing data. }
\description{
 Predict case/control conditions or phenotypes for a user provided testing data using the DC-based prediction model trained previously on a set of training data. 
}
\usage{
BFtest(testdata, 
       model)
}
\arguments{
  \item{testdata}{ a data frame or matrix containing expression of testing data. Columns correspond to genes and rows to samples. The columns must be consistent with the columns in the training data. }

  \item{model}{ DCmodel object represents the prediction model trained by training data. }

}
\details{

 The function predict case/control conditions or phenotypes for a user provided testing data using the DC-based prediction model trained previously on a set of training data. 
 
}
\value{
  
  A vector contains the predicting class labels for the testing data.
 
}
\keyword{ package }
