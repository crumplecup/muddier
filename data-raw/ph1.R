
df_init <- c(0.215108641887659, 0.784891358112342)
df_trans <- matrix(c(-0.000383524637779842, 0,
                     0, -0.00429331795220715), ncol = 2)

ff_init <- c(0.5, 0.5)
ff_trans <- matrix(c(-0.00173661044495578, 0,
                     0, -0.00108379321428908), ncol = 2)

fg_init <- c(0.5, 0.5)
fg_trans <- matrix(c(-0.00157066203705577, 0,
                     0, -0.000396829601339837), ncol = 2)

library(actuar)
df_dist <- pphtype(0:1, df_init, df_trans)

library(data.table)
library(magrittr)
setwd('/home/crumplecup/work/')
ph_cdfs <- fread('ph_cdfs.csv')

df_ph_cdf <- ph_cdfs[1,]
ff_ph_cdf <- ph_cdfs[2,]
fg_ph_cdf <- ph_cdfs[3,]

df_ph_pmf <- to_pmf(unlist(df_ph_cdf))
ff_ph_pmf <- to_pmf(unlist(ff_ph_cdf))
fg_ph_pmf <- to_pmf(unlist(fg_ph_cdf))

index <- char_pmfs %>% rownames %>% as.numeric %>% sort

setwd('/home/crumplecup/work/muddier/')
usethis::use_data(df_ph_cdf)
usethis::use_data(ff_ph_cdf)
usethis::use_data(fg_ph_cdf)
usethis::use_data(df_ph_pmf)
usethis::use_data(ff_ph_pmf)
usethis::use_data(fg_ph_pmf)







