
# correct resident times for inherited age

# subtract inherited age pmf from sample age using convolution
index <- char_pmfs %>% rownames %>% as.numeric
cor_pmfs <- char_pmfs %>% as.matrix

df_ph_pmf <- df_ph_pmf %>% rev
ff_wbl_pmf <- ff_wbl_pmf %>% rev
fg_wbl_pmf <- fg_wbl_pmf %>% rev

# deposit residence time by facies
df_res <- vector(length(index), mode = 'numeric')
ff_res <- vector(length(index), mode = 'numeric')
fg_res <- vector(length(index), mode = 'numeric')


begin <- Sys.time()

for (i in seq_along(charcoal$facies)) {
  if (charcoal$facies[i] == 'DF') {
    cor_pmfs[,i] <- convo(cor_pmfs[,i], df_ph_pmf, index)
    df_res <- df_res + cor_pmfs[ , i]
  }
  if (charcoal$facies[i] == 'FF') {
    cor_pmfs[,i] <- convo(cor_pmfs[,i], ff_wbl_pmf, index)
    ff_res <- ff_res + cor_pmfs[ , i]
  }
  if (charcoal$facies[i] == 'FG') {
    cor_pmfs[,i] <- convo(cor_pmfs[,i], fg_wbl_pmf, index)
    fg_res <- fg_res + cor_pmfs[ , i]
  }
}

df_res <- df_res / length(charcoal$facies[charcoal$facies == 'DF'])
ff_res <- ff_res / length(charcoal$facies[charcoal$facies == 'FF'])
fg_res <- fg_res / length(charcoal$facies[charcoal$facies == 'FG'])

res_times <- matrix(c(df_res, ff_res, fg_res,
                      to_cdf(df_res), to_cdf(ff_res), to_cdf(fg_res),
                      sort(index)), ncol = 7)

end <- Sys.time()
end - begin

setwd('/home/crumplecup/work/muddier/')
usethis::use_data(cor_pmfs, overwrite = T)
usethis::use_data(res_times, overwrite = T)

setwd('/home/crumplecup/work/')
write.csv(res_times, 'res_times.csv')


