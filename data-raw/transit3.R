library(data.table)
library(magrittr)
library(parallel)
library(plot3D)
options(mc.cores=detectCores())
set.seed(10101)

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

setwd('/home/crumplecup/work/')


# fit input-output rate pairs to debris-flow deposit charcoal ages
rec <- accumulater1(length(ffmn), c(0.0,0.7),
                   c(0.0, 0.7), ffia, ffmn, 50, 100)
dur <- as.difftime(3, units = 'hours')
begin <- Sys.time()
end <- begin + dur
while (Sys.time() < end) {
  rec <- rbind(
    rec, accumulater1(length(ffmn), c(0.0,1.0),
                     c(0.0, 1.0), ffia, ffmn, 50, 100))
  save(rec, file = 'ffia_100.rds')
}

pal <- get_palette(c('charcoal', 'crimson', 'ocean', 'gold', 'rose', 'sky', 'hardwood'))
sz <- length(ffmn)
bts <- 100
plot(emp_cdf(ffmn), pch = 20, col = pal[1], log = 'x')
for (i in 1:nrow(ksmin)) {
  lines(emp_cdf(accumulate1(sz, ksmin[i,1] %>% unlist, ksmin[i,2] %>% unlist, dfia, bts)),
        lwd = 2, col = get_palette('crimson', .33))
}

bit <- recorder1(length(ffmn), ksmin[1,1], ksmin[1,2], ffia, min(ffmn))
plot(emp_cdf(bit))
lines(emp_cdf(ffmn), col = pal[4], lwd = 2)

scatter3D(rec$ti, rec$to, rec$ks, pch = 20, phi = 30, theta = 50,
          xlab = "input", ylab = 'output', zlab = 'fit',
          ticktype = 'detailed')

rec <- fread('/home/crumplecup/projects/reservoir/dfio_fit3.csv')
rec <- fread('/home/crumplecup/projects/reservoir/ffio_fit.csv')
rec <- fread('/home/crumplecup/projects/reservoir/fgio_fit.csv')
rec <- fread('/home/crumplecup/projects/reservoir/dfio_1k.csv')
rec <- fread('/home/crumplecup/projects/reservoir/ffio_1k.csv')
rec <- fread('/home/crumplecup/projects/reservoir/fgio_1k.csv')
rec <- fread('/home/crumplecup/projects/reservoir/dfio_200k.csv')
load('/home/crumplecup/work/ffia_100.rds')

scatter3D(rec$input, rec$output, rec$ks, pch = 20, phi = 50, theta = 20,
          xlab = "input", ylab = 'output', zlab = 'fit',
          ticktype = 'detailed')
(ksmin <- rec[rec$ks == min(rec$ks), ])
(kpmin <- rec[rec$kp == min(rec$kp)])

sub <- rec[rec$ks < .2, ]
png('dfio.png', height = 15, width = 18, units = 'cm', res = 300)
scatter3D(sub$input, sub$output, sub$ks, pch = 20, phi = 50, theta = 230,
          xlab = "input", ylab = 'output', zlab = 'fit',
          ticktype = 'detailed')
dev.off()

sz <- length(fgmn)
bts <- 100
plot(emp_cdf(fgmn), xlab = 'charcoal age', ylab = 'cdf',
     pch = 20, col = get_palette('charcoal', .4), log = 'x')
for (i in 1:nrow(ksmin)) {
  lines(emp_cdf(accumulate(sz, ksmin[i,1] %>% unlist, ksmin[i,2] %>% unlist, dfia, bts)),
        lwd = 2, col = get_palette('crimson', .33))
}
for (i in 1:nrow(kpmin)) {
  lines(emp_cdf(accumulate(sz, kpmin[i,1] %>% unlist, kpmin[i,2] %>% unlist, dfia, bts)),
        lwd = 2, col = get_palette('ocean', .33))
}

legend('bottomright', legend = c('observed', 'K-S fit', 'Kuiper fit'),
       fill = get_palette(c('charcoal', 'crimson', 'ocean'), .7))



rec <- fread('/home/crumplecup/projects/reservoir/ffio_100k.csv')
scatter3D(rec$input, rec$output, rec$ks, pch = 20, phi = 50, theta = 20,
          xlab = "input", ylab = 'output', zlab = 'fit',
          ticktype = 'detailed')
(ksmin <- rec[rec$ks == min(rec$ks)])
(kpmin <- rec[rec$kp == min(rec$kp)])

sz <- length(ffmn)
bts <- 100
plot(emp_cdf(ffmn), xlab = 'charcoal age', ylab = 'cdf',
     pch = 20, col = get_palette('charcoal', .4))
for (i in 1:nrow(ksmin)) {
  lines(emp_cdf(accumulate(sz, ksmin[i,1] %>% unlist, ksmin[i,2] %>% unlist, dfia, bts)),
        lwd = 2, col = get_palette('crimson', .33))
}
for (i in 1:nrow(kpmin)) {
  lines(emp_cdf(accumulate(sz, kpmin[i,1] %>% unlist, kpmin[i,2] %>% unlist, dfia, bts)),
        lwd = 2, col = get_palette('ocean', .33))
}

legend('bottomright', legend = c('observed', 'K-S fit', 'Kuiper fit'),
       fill = get_palette(c('charcoal', 'crimson', 'ocean'), .7))


rec <- fread('/home/crumplecup/projects/reservoir/fgio_100k.csv')
scatter3D(rec$input, rec$output, rec$ks, pch = 20, phi = 50, theta = 20,
          xlab = "input", ylab = 'output', zlab = 'fit',
          ticktype = 'detailed')
(ksmin <- rec[rec$ks == min(rec$ks)])
(kpmin <- rec[rec$kp == min(rec$kp)])

sz <- length(fgmn)
bts <- 100
plot(emp_cdf(fgmn), xlab = 'charcoal age', ylab = 'cdf',
     pch = 20, col = get_palette('charcoal', .4))
for (i in 1:nrow(ksmin)) {
  lines(emp_cdf(accumulate(sz, ksmin[i,1] %>% unlist, ksmin[i,2] %>% unlist, dfia, bts)),
        lwd = 2, col = get_palette('crimson', .33))
}
for (i in 1:nrow(kpmin)) {
  lines(emp_cdf(accumulate(sz, kpmin[i,1] %>% unlist, kpmin[i,2] %>% unlist, dfia, bts)),
        lwd = 2, col = get_palette('ocean', .33))
}

legend('bottomright', legend = c('observed', 'K-S fit', 'Kuiper fit'),
       fill = get_palette(c('charcoal', 'crimson', 'ocean'), .7))



rec <- fread('/home/crumplecup/projects/reservoir/df_boot.csv')
rec <- fread('/home/crumplecup/projects/reservoir/df_boot1.csv')
rec <- fread('/home/crumplecup/projects/reservoir/df_boot2.csv')

med <- apply(rec, 1, sort)
med <- apply(med, 1, median)
plot(emp_cdf(med))
lines(emp_cdf(fgmn))
medgof <- cdf_gof(med, fgmn)[[1]]
gofs <- apply(rec, 1, function(x) cdf_gof(x, med))
ln <- length(gofs)
gof <- 0
for (i in seq_along(gofs)) {
  gof[i] <- gofs[[i]]$ks
}
# K-S distribution
length(gof[gof < medgof]) / ln  # 0.70886 for 100k samples

plot(emp_cdf(dfmn),
     xlab = 'transit time', ylab = 'cdf')
for (i in 1:nrow(ksmin)) {
  lines(emp_cdf(accumulate(sz, 0.22, 0.207, dfia, bts)),
        lwd = 2, col = get_palette('crimson', .33))
}
for (i in 1:nrow(kpmin)) {
  lines(emp_cdf(accumulate(sz, kpmin[i,1], kpmin[i,2], ia_ar, bts)),
        lwd = 2, col = get_palette('ocean', .33))
}
meds <- apply(rec, 2, sort)
meds <- apply(meds, 1, median)
lines(emp_cdf(meds), lwd = 3, col = get_palette('gold'))
apply(rec, 2, function(x) lines(emp_cdf(x)))


df <- data.frame(age = sort(dfmn), facies = 'df')
ff <- data.frame(age = ffmn, facies = 'ff')
fg <- data.frame(age = fgmn, facies = 'ff')
df <- rbind(df, rbind(ff, fg))

d <- dist(df[,1])
c <- hclust(d)
plot(c)

k <- kmeans(df[,1], centers = 11)
k

bit <- lapply(1:1000, function(x) sort(sample(dfmn, length(dfmn), replace = T))) %>% as.data.table
bit <- apply(bit, 1, median)
plot(emp_cdf(bit), log = 'x')
lines(emp_cdf(dfmn), col = get_palette('gold', .7))


library(data.table)
library(plot3D)
rec1 <- fread('/home/crumplecup/work/output/fgio_1kstdy_30k.csv')
rec <- fread('/home/crumplecup/work/output/fgio_10kstdy_30k.csv')
rec3 <- fread('/home/crumplecup/work/output/ff_1kstdy_30k.csv')
rec2 <- fread('/home/crumplecup/work/output/ff_10kstdy_30k.csv')
rec5 <- fread('/home/crumplecup/work/output/df_1kstdy_30k.csv')
rec4 <- fread('/home/crumplecup/work/output/df_10kstdy_30k.csv')
grav <- fread('/home/crumplecup/work/output/fg100k_stdy.csv')
fine <- fread('/home/crumplecup/work/output/ff100k_stdy.csv')
debr <- fread('/home/crumplecup/work/output/df100k_stdy.csv')
debr20 <- fread('/home/crumplecup/work/output/df20k_stdy.csv')
debr1 <- fread('/home/crumplecup/work/output/df1_stdy.csv')
debr5 <- fread('/home/crumplecup/work/output/df5k.csv')
pal <- get_palette(c('sky', 'ocean', 'gold', 'hardwood', 'leaf', 'forest', 'charcoal'))


plot(debr5$input, debr5$ks, col = pal[3], pch = 20)

plot(debr1$input, debr1$ks, col = pal[7], pch = 20)
points(debr20$input, debr20$ks, col = pal[6], pch = 20)
points(debr5$input, debr5$ks, col = pal[2], pch = 20)
debr20[debr20$ks == min(debr20$ks), ]
plot(debr$input, debr$ks, col = pal[7], pch = 20)
points(debr20$input, debr20$ks, col = pal[1], pch = 20)

plot(debr20$input, debr20$ks, col = pal[7], pch = 20)
debr20[debr20$ks == min(debr20$ks), ]

plot(rec1$input, rec1$ks, col = pal[1], pch = 20,
     xlab = 'rate', ylab = 'K-S fit')
points(rec$input, rec$ks, col = pal[2], pch = 20)
points(grav$input, grav$ks, col = pal[7], pch = 20)
legend('bottomright', legend = c('gravels 1k', 'gravels 10k', 'gravels 100k'),
       fill = pal[c(1:2, 7)])
plot(grav$input, grav$ks, col = pal[7], pch = 20)
grav[grav$ks == min(grav$ks), ]

rec3 <- rec3[rec3$input > .3]
plot(rec3$input, rec3$ks, col = pal[3], pch = 20,
     xlab = 'rate', ylab = 'K-S fit')
points(rec2$input, rec2$ks, col = pal[4], pch = 20)
points(fine$input, fine$ks, col = pal[7], pch = 20)
legend('bottomright', legend = c('fines 1k', 'fines 10k', 'fines 100k'),
       fill = pal[c(3:4, 7)])

plot(fine$input, fine$ks, col = pal[7], pch = 20)
fine[fine$ks == min(fine$ks), ]
plot(debr$input, debr$ks, col = pal[7], pch = 20)
debr[debr$ks == min(debr$ks), ]

rec5 <- rec5[rec5$input > .3]
plot(rec5$input, rec5$ks, col = pal[5], pch = 20,
     xlab = 'rate', ylab = 'K-S fit')
points(rec4$input, rec4$ks, col = pal[6], pch = 20)
legend('topright', legend = c('debris flows 1k', 'debris flows 10k'),
       fill = pal[5:6])

rec4[rec4$ks == min(rec4$ks),]

scatter3D(rec$input, rec$output, rec$ks, ticktype = 'detailed', phi = 60, theta = 230)
rec[rec$ks < .18,]

plot(rec1$input, rec1$ks, col = pal[1],
     xlab = 'rate', ylab = 'K-S fit', pch = 20,
     ylim = c(min(rec5$ks), max(rec3$ks)),
     xlim = c(min(rec3$input), max(rec3$input)))
points(rec$input, rec$ks, col = pal[2], pch = 20)
points(rec3$input, rec3$ks, col = pal[3], pch = 20)
points(rec5$input, rec5$ks, col = pal[5], pch = 20)

legend('topright', legend = c('gravels 1k', 'gravels 10k', 'fines 1k', 'fines 10k', 'debris flows 1k', 'debris flows 10k'),
       fill = pal)



