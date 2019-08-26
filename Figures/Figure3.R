# LT 8/09/2018

# cleaned 17/08/2019
require(TMB)
load("Hierarchical.Rdata")

# use color-blind friendly color palette
cols = c(purp = rgb(.8, .6, .7), black = rgb(0,0,0), blue=rgb(0.35, 0.70, 0.9) , orange=rgb(.8, .4, 0), green=rgb(0, .6, .5))

g = allmodelfits$hierarchical$hierarchical$f
# rebuilds bs for hierarchical models
bs= g$env$parList()$b + exp(g$env$parList()$lnu0) * g$env$parList()$bdevs0[order.fact] +
  exp(g$env$parList()$lnu1) * g$env$parList()$bdevs1

# plot density of bs per class
#insects
plot(density(subset(bs, class.fact=="Insecta")), xlim=c(0.5,1),
     col=cols["green"], lwd=2,
     xlab="strength of density dependence (b)",
     main = "Density plots of estimates of b")
#fishes
lines(density(subset(bs, class.fact=="Osteichthyes")), col=cols["blue"],lwd=2)
lines(density(subset(bs, class.fact=="Aves")), col=cols["orange"],lwd=2)
lines(density(subset(bs, class.fact=="Mammalia")), col=cols["purp"],lwd=2)

legend("topright",legend=c("Insecta", "Osteichthyes", "Aves", "Mammalia"), lwd=2, col=cols[c("green", "blue", "orange", "purp")])

