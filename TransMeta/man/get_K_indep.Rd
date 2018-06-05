\name{get_K_indep}
\alias{ get_K_indep }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct the group-wise independent kernel matrix K 
}
\description{Construct the group-wise independent kernel matrix K for testing the single variant association in trans-ethnic meta-analysis.
}
\usage{
get_K_indep(FileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ FileName }{A vector of characters that list the file names of summary statistics for each GWAS and the file name of a data file that identifies each study's corresponding ancestry group. 
}
}
\details{Please refer to the Vignettes PDF manual for details of the FileName argument.}
\value{Returns the group-wise independent kernel matrix K for testing the single variant association in trans-ethnic meta-analysis. For example, the returned matrix can be used for the argument K in the Get_TransMeta function.
}
\examples{## Please refer to the example in the Vignettes PDF manual 
}

