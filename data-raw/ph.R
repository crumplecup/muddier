library(data.table)
library(magrittr)
setwd('/home/crumplecup/work/')
ph_cdfs <- fread('inher_ph_cdfs.csv')

df_ph_cdf <- ph_cdfs[1,]
ff_wbl_cdf <- ph_cdfs[2,]
fg_wbl_cdf <- ph_cdfs[3,]

df_ph_pmf <- to_pmf(unlist(df_ph_cdf))
ff_wbl_pmf <- to_pmf(unlist(ff_wbl_cdf))
fg_wbl_pmf <- to_pmf(unlist(fg_wbl_cdf))

setwd('/home/crumplecup/work/muddier/')
usethis::use_data(df_ph_cdf, overwrite = T)
usethis::use_data(ff_wbl_cdf, overwrite = T)
usethis::use_data(fg_wbl_cdf, overwrite = T)
usethis::use_data(df_ph_pmf, overwrite = T)
usethis::use_data(ff_wbl_pmf, overwrite = T)
usethis::use_data(fg_wbl_pmf, overwrite = T)







