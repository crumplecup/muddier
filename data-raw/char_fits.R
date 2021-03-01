library(data.table)
library(plot3D)
setwd('/home/crumplecup/work/')
load('ffia_100.rds')

scatter3D(rec$ti, rec$to, rec$ks, pch = 20, phi = 30, theta = 350,
          xlab = "input", ylab = 'output', zlab = 'fit',
          ticktype = 'detailed')

rec <- rec[rec$ks <= .2, ]
scatter3D(rec$ti, rec$to, rec$ks, pch = 20, phi = 30, theta = 320,
          xlab = "input", ylab = 'output', zlab = 'fit',
          ticktype = 'detailed')


(frame <- focus(rec$ks <= .15, rec))


dfmn <- fread('/home/crumplecup/work/dfmn.csv')
dfmn <- dfmn$x
ffmn <- fread('/home/crumplecup/work/ffmn.csv')
ffmn <- ffmn$x
fgmn <- fread('/home/crumplecup/work/fgmn.csv')
fgmn <- fgmn$x
iat <- fread('/home/crumplecup/work/iat.csv')
dfia <- iat$ia[iat$type == 'DF']
ffia <- iat$ia[iat$type == 'FF']
fgia <- iat$ia[iat$type == 'FG']

pal <- get_palette(c('charcoal', 'crimson', 'ocean', 'gold', 'rose', 'sky', 'hardwood'))
plot(emp_cdf(iat$ia), pch = 20, col = pal[1], log = 'x')
lines(emp_cdf(iat$ia), col = pal[1])
points(emp_cdf(dfia), pch = 20, col = pal[2])
lines(emp_cdf(dfia), col = pal[2])
points(emp_cdf(ffia), pch = 20, col = pal[3])
lines(emp_cdf(ffia), col = pal[3])
points(emp_cdf(fgia), pch = 20, col = pal[4])
lines(emp_cdf(fgia), col = pal[4])

points(emp_cdf(dfmn), pch = 20, col = pal[5])
lines(emp_cdf(dfmn), col = pal[5])
points(emp_cdf(ffmn), pch = 20, col = pal[6])
lines(emp_cdf(ffmn), col = pal[6])
points(emp_cdf(fgmn), pch = 20, col = pal[7])
lines(emp_cdf(fgmn), col = pal[7])

library(data.table)
(ksmin <- rec[rec$ks == min(rec$ks), ])
(kpmin <- rec[rec$kp == min(rec$kp), ])

library(magrittr)
sz <- length(ffmn)
bts <- 100
plot(emp_cdf(ffmn), pch = 20, col = pal[1], log = 'x')
for (i in 1:nrow(ksmin)) {
  lines(emp_cdf(accumulate(sz, ksmin[i,1] %>% unlist, ksmin[i,2] %>% unlist, dfia, bts)),
        lwd = 2, col = get_palette('crimson', .33))
}




