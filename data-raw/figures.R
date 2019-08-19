#png('area_pdel.png')
plot(trans$area,trans$DebrisFlow,col= rgb(red=1, green=0, blue=.2, alpha=.5),pch=20,
     main='Cross-sectional Area vs. Delivery Probability',
     xlab='Cross-sectional Area (m2)',ylab='Delivery Probability' )
abline(lm(trans$DebrisFlow ~ trans$area),col='gray')
#dev.off()
