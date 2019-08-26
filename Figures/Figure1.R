# LT 8/09/2018

# cleaned up 16/08/2019

# analyze Knape dataset

require(TMB)
load("Knape.Rdata")

# bs for full model
bs.fixed = sapply(allmodelfits$full, function(m) summary(m$sd)[1,1])

# bs for random effect model
rnd.par = allmodelfits$randomb$randomb$f$env$parList()
bs.rnd = rnd.par$b + exp(rnd.par$lnu) * rnd.par$bdevs

# use color-blind friendly color palette
cols = c(purp = rgb(.8, .6, .7), black = rgb(0,0,0), blue=rgb(0.35, 0.70, 0.9) , orange=rgb(.8, .4, 0), green=rgb(0, .6, .5))

#postscript(file="Figure1.ps",width=3, height=3)
plot(density(bs.rnd), xlim=c(-1.0, 1.1), ylim=c(0,7),col=cols[3], lwd=2,
     main="",
     xlab = "Estimates of density-dependent parameter (b)", ylab = "density")
lines(density(bs.fixed), lwd=2, col=cols[2])

# add uncertainty on full b from boostraping
points(0.59, 0.1, col=cols[2], pch = 20, cex=1.4)
lines(c(0.5,0.68), rep(0.1,2), col=cols[2], lwd=3)

# add uncertainty on mean b for random b model
sdb = summary(allmodelfits$randomb$randomb$sd)[1,]
points(sdb[1], 0.2, col=cols[3], pch = 20, cex=1.4)
lines(c(sdb[1]-qnorm(0.025)*sdb[2],sdb[1]+qnorm(0.025)*sdb[2]), rep(0.2,2), col=cols[3], lwd=3)

legend("topleft",legend=c("full model", paste("random-", "b")), lwd=2, col=cols[c("black", "blue")])
#dev.off()

