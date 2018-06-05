\name{F_ST}
\alias{F_ST}
\docType{data}
\title{
  A 9 by 9 symmetric matrix of the pair-wise F.st values for the 9 ancestry groups from HMP3 data.
}
\description{
A 9 by 9 symmetric matrix of the pair-wise F.st values for the 9 ancestry groups from HMP3 data. The 9 ancestry groups are (in the following order): CEU - Utah residents with Northern and Western European ancestry from the CEPH collection; CHD - Chinese in Metropolitan Denver, Colorado; GIH - Gujarati Indians in Houston, Texas; JPT_CHB - Japanese in Tokyo, Japan & Han Chinese in Beijing, China; LWK - Luhya in Webuye, Kenya; MEX - Mexican ancestry in Los Angeles, California; MKK - Maasai in Kinyawa, Kenya; TSI - Toscani in Italia; YRI - Yoruba in Ibadan, Nigeria. 
}
\usage{data(F_ST)}
\details{
A 9 \eqn{\times} 9 symmetric matrix of the pair-wise F.st values for the 9 ancestry groups from HMP3 data. The (i,j)-th entry value in the matrix is the pair-wise F.st value between the i-th ancestry group and the j-th ancestry group. Order of the 9 ancestry groups is arranged as follows: CEU, CHD, GIH, JPT+CHB, LWK, MEX, MKK, TSI, YRI.
}
\source{
The pair-wise F.st values are obtained from Table S6 in: 

http://hapmap.ncbi.nlm.nih.gov/downloads/presentations/nature09298-s1.pdf}
\examples{
data(F_ST)
F_ST
}
\keyword{datasets}
