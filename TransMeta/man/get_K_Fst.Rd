\name{get_K_Fst}
\alias{ get_K_Fst }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct the genetic similarity kernel matrix K 
}
\description{
Construct the genetic similarity kernel matrix K for testing the single variant association in trans-ethnic meta-analysis.}
\usage{
get_K_Fst(FileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FileName}{A vector of characters that list the file names of summary statistics for each GWAS and the file name of a data file that identifies each study's corresponding ancestry group. 
}
}
\details{Please refer to the Vignettes PDF manual for details of the FileName argument.

}
\value{Returns the genetic similarity kernel matrix K for testing the single variant association in trans-ethnic meta-analysis. For example, the returned matrix can be used for the argument K in the Get_TransMeta function.
}

