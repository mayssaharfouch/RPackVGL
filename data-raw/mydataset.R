
# data-raw/mydataset.R
# Data import and processing pipeline
# https://grasshoppermouse.github.io/posts/2017-10-18-put-your-data-in-an-r-package/


library(readxl)


# Data cleaning code here...
# (Do NOT put data analysis code here!)

## input files examples

#plant parameter file
path_plante <- system.file("extdata", "Parametres_plante_exemple.xls", package = "RPackVGL")
plt_par <- read_excel(path_plante, sheet="Fix2")

#...


## output file examples
path_zip0 <- system.file("extdata", "toto_0.zip", package = "RPackVGL")
toto_0 <- read.table(unz(path_zip0, "toto_0.csv"), header=T, sep=";")
outHR_0 <- read.table(unz(path_zip0, "outHR_0.csv"), header=T, sep=";")
BilanN_0 <- read.table(unz(path_zip0, "BilanN_0.csv"), header=T, sep=";")
paramsd_0 <- read.table(unz(path_zip0, "paramsd_0.csv"), header=T, sep=";")



# This should be the last line.
# Note that names are unquoted.
# I like using overwrite = T so everytime I run the script the
# updated objects are saved, but the default is overwrite = F
usethis::use_data(plt_par, toto_0, outHR_0, overwrite = T)


