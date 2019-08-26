# LT 8/09/2018

# cleaned 17/08/2019
require(TMB)
load("Hierarchical.Rdata")

# get 95% orders
fit = allmodelfits$hierarchical$hierarchical$fit
sds = allmodelfits$hierarchical$hierarchical$sd
cis = cbind(summary(sds,"random")[1:14,1] + qnorm(0.975) * summary(sds,"random")[1:14,2],
            summary(sds,"random")[1:14,1] - qnorm(0.975) * summary(sds,"random")[1:14,2])
cis= cis[,2:1] * exp(fit$par["lnu0"])
cis = cis + fit$par["b"]

means = summary(sds,"random")[1:14,1] * exp(fit$par["lnu0"]) + fit$par["b"]
cis = cbind(mean=means, cis)

rownames(cis) = levels(order.fact)

# plot figure 2
# use color-blind friendly color palette
cols = c(purp = rgb(.8, .6, .7), black = rgb(0,0,0), blue=rgb(0.35, 0.70, 0.9) , orange=rgb(.8, .4, 0), green=rgb(0, .6, .5))

# random b model
y=0
plot(0,0,xlim=c(-0, 1.5), ylim = c(-20,-0), type="n", ylab="", xlab="density dependence b", yaxt='n', bty='n', ann=FALSE)
points(allmodelfits$randomb$randomb$fit$par["b"], y, col=1, pch = 20)
rndb.sd = summary(allmodelfits$randomb$randomb$sd)
rndb.confint = c(rndb.sd[1,1] + qnorm(0.025)*rndb.sd[1,2], rndb.sd[1,1] + qnorm(0.9725)*rndb.sd[1,2])
lines(rndb.confint, rep(y,2), col=1, lwd=3)
text(1.3, y, "random-b", col=1, cex=0.9)

# hierarchical model
y = y - 1
points(allmodelfits$hierarchical$hierarchical$fit$par["b"], y, col=1, pch = 20)
hib.sd = summary(allmodelfits$hierarchical$hierarchical$sd)
hib.confint = c(hib.sd[1,1] + qnorm(0.025)*hib.sd[1,2], hib.sd[1,1] + qnorm(0.9725)*hib.sd[1,2])
lines(hib.confint, rep(y,2), col=1, lwd=3)
text(1.3, y, "hierarchical-b", col=1, cex=0.9)

# get the insects
insects = levels(droplevels(subset(order.fact, class.fact=="Insecta")))
y = y - 1
text(1.3, y, "Insecta", col=cols["green"], cex=0.9)
for (insect in insects) {
  y = y - 1
  points(cis[insect,1], y, col=cols["green"], pch = 20)
  lines(cis[insect, 2:3], rep(y,2), col=cols["green"], lwd=3)
  text(1.3, y, insect, col=cols["green"], cex=0.7)
}

# get the fishes
fishes = levels(droplevels(subset(order.fact, class.fact=="Osteichthyes")))
y = y - 1
text(1.3, y, "Osteichthyes", col=cols["blue"], cex=0.9)
for (fish in fishes) {
  y = y - 1
  points(cis[fish,1], y, col=cols["blue"], pch = 20)
  lines(cis[fish, 2:3], rep(y,2), col=cols["blue"], lwd=3)
  text(1.3, y, fish, col=cols["blue"], cex=0.7)
}

# get the birds
birds = levels(droplevels(subset(order.fact, class.fact=="Aves")))
y = y - 1
text(1.3, y, "Aves", col=cols["orange"], cex=0.9)
for (bird in birds) {
  y = y - 1
  points(cis[bird,1], y, col=cols["orange"], pch = 20)
  lines(cis[bird, 2:3], rep(y,2), col=cols["orange"], lwd=3)
  text(1.3, y, bird, col=cols["orange"], cex=0.7)
}

# get the mammals
mammals = levels(droplevels(subset(order.fact, class.fact=="Mammalia")))
y = y - 1
text(1.3, y, "Mammalia", col=cols["purp"], cex=0.9)
for (mam in mammals) {
  y = y - 1
  points(cis[mam,1], y, col=cols["purp"], pch = 20)
  lines(cis[mam, 2:3], rep(y,2), col=cols["purp"], lwd=3)
  text(1.3, y, mam, col=cols["purp"], cex=0.7)
}