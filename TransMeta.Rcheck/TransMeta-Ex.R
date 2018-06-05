pkgname <- "TransMeta"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('TransMeta')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("F_ST")
### * F_ST

flush(stderr()); flush(stdout())

### Name: F_ST
### Title: A 9 by 9 symmetric matrix of the pair-wise F.st values for the 9
###   ancestry groups from HMP3 data.
### Aliases: F_ST
### Keywords: datasets

### ** Examples

data(F_ST)
F_ST



cleanEx()
nameEx("G_ind")
### * G_ind

flush(stderr()); flush(stdout())

### Name: G_ind
### Title: K matrix of the group-wise independent kernel structure for the
###   9 ancestry groups from HMP3 data.
### Aliases: G_ind
### Keywords: datasets

### ** Examples


data(G_ind)
G_ind



cleanEx()
nameEx("Get_TransMeta")
### * Get_TransMeta

flush(stderr()); flush(stdout())

### Name: Get_TransMeta
### Title: Single variant association test in the GWAS trans-ethnic
###   meta-analysis
### Aliases: Get_TransMeta

### ** Examples
## Please refer to the example in the Vignettes PDF manual 


cleanEx()
nameEx("get_K_indep")
### * get_K_indep

flush(stderr()); flush(stdout())

### Name: get_K_indep
### Title: Construct the group-wise independent kernel matrix K
### Aliases: ' get_K_indep '

### ** Examples
## Please refer to the example in the Vignettes PDF manual 



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
