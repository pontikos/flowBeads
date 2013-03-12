library(flowCore)
library(flowBeads)
b <- BeadFlowFrame(fcs.filename='~nikolas/CAD54_2008FEB11_TREG_6BEADS_110208.FCS')
gb <- gateBeads(b)


plot(gb)

flow.data <- read.FCS('~nikolas/CAD54_2008FEB11_TREG_6BEADS_110208.FCS')
flow.data.2 <- normaliseBeads(bead.data=gb, flow.data=flow.data)


beads.day1 <- BeadFlowFrame('beads-day1.fcs')
g.beads.day1 <- gateBeads(beads.day1, verbose=T)
beads.day2 <- BeadFlowFrame('beads-day2.fcs')
g.beads.day2 <- gateBeads(beads.day2, verbose=T)
person.day1 <- read.FCS('person-day1.fcs')
person.day2 <- read.FCS('person-day2.fcs')

#
lim <- 1.5
trans <- g.beads.day1@trans
x.day1 <- person.day1@exprs[,5]
x.day1 <- x.day1[trans(x.day1) > lim]
day1.trans <- g.beads.day1@mef.transform$APC$fun
x.day1.mef <- day1.trans(x.day1)
#
x.day2 <- person.day2@exprs[,5]
x.day2 <- x.day2[trans(x.day2) > lim]
day2.trans <- g.beads.day2@mef.transform$APC$fun
x.day2.mef <- day2.trans(x.day2)

#
par(mfrow=c(2,1))
plot( density(trans(x.day1)), main=abs(mean(trans(x.day1))-mean(trans(x.day2))) )
abline(v=mean(trans(x.day1)), col='black')
lines(density(trans(x.day2)), col='red')
abline(v=mean(trans(x.day2)), col='red')
#
plot( density(trans(x.day1.mef)), col='black', lty=2, main=abs(mean(trans(x.day1.mef))-mean(trans(x.day2.mef))) )
abline(v=mean(trans(x.day1.mef)), col='black', lty=2)
lines(density(trans(x.day2.mef)), col='red', lty=2)
abline(v=mean(trans(x.day2.mef)), col='red', lty=2)

par(mfrow=c(2,1))
d <- rbind( data.frame(x=trans(x.day1), y=1), data.frame(x=trans(x.day2), y=2) )
boxplot( x ~ y, data=d, horizontal=T, main=abs(mean(trans(x.day1))-mean(trans(x.day2))) )
d <- rbind( data.frame(x=trans(x.day1.mef), y=1), data.frame(x=trans(x.day2.mef), y=2) )
boxplot( x ~ y, data=d, horizontal=T, main=abs(mean(trans(x.day1.mef))-mean(trans(x.day2.mef))) )
par(mfrow=c(1,1))

plot(trans(x.day1), trans(x.day1.mef), pch='.')
l <- line(trans(x.day1), trans(x.day1.mef))
abline(coef(l), col='blue')
print(l)

par(mfrow=c(2,1))
plot(g.beads.day1, 'APC')
plot(g.beads.day2, 'APC')


data(gbeads1)
data(gbeads2)
l <- relTransform(gbeads1, gbeads2)
a <- l$ALEXA.700
plot(a$mfi1, a$mfi2)
abline(b=1,a=0)
plot(a$fun(a$mfi1), a$mfi2)
abline(b=1,a=0)
points(a$beta*a$mfi1+a$alpha, a$mfi2, col='red')


