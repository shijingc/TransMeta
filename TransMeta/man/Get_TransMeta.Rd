\name{Get_TransMeta}
\alias{Get_TransMeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Single variant association test in the GWAS trans-ethnic meta-analysis 
}
\description{
Implementation of the kernel regression based random effect model for the single variant association test in the GWAS trans-ethnic meta-analysis
}
\usage{
Get_TransMeta(FileName, type = "effect size", K.type = "Fst", K = NULL,

r.all = c(0, 0.09, 0.25, 1), method = "davies")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FileName}{A vector of characters that list the file names of summary statistics for each GWAS and the file name of a data file that identifies each study's corresponding ancestry group. Please refer to the Vignettes PDF manual for its details.  }
  \item{type}{The default value is "effect size", which specifies that the input summary statistics are effect size based; otherwise, the Z-score based summary statistics will be used. }
  \item{K.type}{The default value is "Fst", which specifies that the K matrix is of the genetic similarity kernel structure; otherwise, the K matrix would take the group-independent kernel structure.}
  \item{K}{A square matrix for the kernel matrix K, for example, K can be constructed based on the HapMap3 data using the get_K_Fst function. If the user does not provide a matrix for this argument, then the K matrix will be constructed based on the input data (to estimate F.st) with its specified structure from the K.type argument. }
  \item{r.all}{A vector of possible values of \eqn{\rho} to be used, default = c(0, 0.09, 0.25, 1).}
  \item{method}{A method to compute the p-value (default= \sQuote{davies}). \sQuote{davies} represents an exact method that computes the p-value by inverting the characteristic function of the mixture chisq; \sQuote{liu} represents an approximation method that matches the first 3 moments; \sQuote{liu.mod } represents modified \sQuote{liu} method that matches kurtosis instead of skewness to improve tail probability approximation.}
}
\value{
  \item{ SNP }{The SNP id for each variant}
\item{ Pval }{p-value from TransMeta for the association test}
  \item{ Pval_rho }{ p-value of the score test for each \eqn{\rho} in the \sQuote{r.all } argument}
\item{ Opt_rho }{ the optimal \eqn{\rho} in the \sQuote{r.all} argument which results in smallest p-value for the score test }

}
\examples{## Please refer to the example in the Vignettes PDF manual }

