
### Polynomial Wavelet Analysis
###
### (C) Leonard Mada
### draft 0.2
###
### Cardano-type Polynomials:
### theory developed during 2018-2019


### Related / Useful Materials:

# 1. Ecological Forecasting: Statistical Methods Series:
#    wsyn: Wavelet approaches to Studies of Synchrony in Ecology
#    https://www.youtube.com/watch?v=sxVHFwAkhdA
# 2. Package wsyn: Wavelet Approaches to Studies of Synchrony
#    in Ecology and Other Fields
#    https://cran.r-project.org/web/packages/wsyn/index.html
# 3. Sheppard LW, Defriez EJ, Reid PC, Reuman DC.
#    Synchrony is more than its top-down and climatic parts:
#    interacting Moran effects on phytoplankton in British seas.
#    PLoS Comput Biol 2019; 15(3): e1006744.
#    https://doi.org/10.1371/journal.pcbi.1006744


############

### History:

# draft 0.2: 2020-05-27
# - generating code;
# - disclosed solution (brief);
# draft 0.1: 2020-03-18
# - R code to explore some properties;
# pre-0.1: pre 2020-02
# - initial thoughts about Wavelets;
# 2018-2019:
# - development of consistent Theory about the Polynomials;
#   [was disseminated, but was never published]

### TODO:
# - prepare some brief material for:
#   Generalization of the Cardano formula;
# - disseminate material to enable proper presentations
#   e.g. Two Brilliant Ways to Solve This Non Factorable Cubic Equation
#   https://www.youtube.com/watch?v=Jz-Z0jfs2V4

#############

# Constants for Wavelets:
c = 1
d = 0

### "Quadratic" / Cardan-type Polynomials
# Odd-Order
f3  = function(x) x^3 - 3*c*x - 2*d
f5  = function(x) x^5 - 5*c*x^3 + 5*c^2*x - 2*d
f7  = function(x) x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d
f9  = function(x) x^9 - 9*c*x^7 + 27*c^2*x^5 - 30*c^3*x^3 + 9*c^4*x - 2*d
f11 = function(x) x^11 - 11*c*x^9 + 44*c^2*x^7 - 77*c^3*x^5 + 55*c^4*x^3 - 11*c^5*x - 2*d
f13 = function(x) x^13 - 13*c*x^11 + 65*c^2*x^9 - 156*c^3*x^7 + 182*c^4*x^5 - 91*c^5*x^3 + 13*c^6*x - 2*d
f15 = function(x) x^15 - 15*c*x^13 + 90*c^2*x^11 - 275*c^3*x^9 + 450*c^4*x^7 - 378*c^5*x^5 + 140*c^6*x^3 - 15*c^7*x - 2*d
# Even-Order
f2 = function(x) x^2 - 2*c - 2*d
f4 = function(x) x^4 - 4*c*x^2 + 2*c^2 - 2*d
f6 = function(x) x^6 - 6*c*x^4 + 9*c^2*x^2 - 2*c^3 - 2*d
f8 = function(x) x^8 - 8*c*x^6 + 20*c^2*x^4 - 16*c^3*x^2 + 2*c^4 - 2*d
f10 = function(x) x^10 - 10*c*x^8 + 35*c^2*x^6 - 50*c^3*x^4 + 25*c^4*x^2 -2*c^6 - 2*d
f12 = function(x) x^12 - 12*c*x^10 + 54*c^2*x^8 - 112*c^3*x^6 + 105*c^4*x^4 - 36*c^5*x^2 + 2*c^6 - 2*d
f14 = function(x) x^14 - 14*c*x^12 + 77*c^2*x^10 - 210*c^3*x^8 + 294*c^4*x^6 - 196*c^5*x^4 + 49*c^6*x^2 -2*c^7 - 2*d

# All / All Odds
f.all = c(f3, f5, f7, f9, f11, f13, f15)
cols.all = c("green", "springgreen2", "lightseagreen", "red", "blue", "purple", "pink", "yellow")

### Test Polynomial Solutions
n = 14
# Parameters
c = 1
d = 2
# Roots
det=(d^2 - c^n)^(1/2)
det
x = (d+det)^(1/n) + (d-det)^(1/n)
x

#############
### Generator
# Even-powered Polynomials
poly.cardan.gen = function(n) {
	if(n %% 2 == 0) {
		coeff = list("c0"=1, "c2"=c(1, -2))
	}
	for(i in seq(4, n, by=2)) {
		coeff = expand(i, coeff)
	}
	return(coeff)
}
expand = function(n, coeff) {
	last = n/2
	ch = choose(n, seq(n/2, n-1, by=1))
	tmp = 0
	for(id in last:1) {
		tmp = tmp + c(rep(0,last-id), coeff[[id]]) * ch[id]
	}
	ch = c(1, -tmp)
	coeff[[paste("c", n, sep="")]] = ch
	return(coeff)
}

###
coeff = poly.cardan.gen(16)
coeff

# generate polynomial
n = 14
paste(coeff[[n/2 + 1]], "*c^", seq(0, n/2), "*x^", seq(n, 0, by=-2), " + ", sep="", collapse="")


# not needed
fh = list()
for(id in 1:length(f.all)) {
	f = f.all[[id]]
	fh = c(fh, function(x) {
		f(x) / 2
	})
}
lapply(f.all, function(x) x(-2))

##################

##################
### Properties ###

# These polynomials have very interesting properties:
# 1.) P(0) = -2*d;
# 2.) for c = 1, d = 0:
# all values of P(x) are between:[-2, +2] for x in [-2, 2];
# P(-2) = -2; P(2) = 2;

### Visualization:
from = -2.005
to = 2.005
ylim = c(-3, 3) # NULL

curve(f15(x), from=from, to=to, ylim=ylim)
curve(f13(x), from=from, to=to, ylim=ylim, add=T, col="springgreen2", n=256)
curve(f11(x), from=from, to=to, ylim=ylim, add=T, col="lightseagreen", n=256)
curve(f9(x), from=from, to=to, add=T, col="green", ylim=ylim, n=256)
curve(f7(x), from=from, to=to, add=T, col="red", ylim=ylim, n=256)
curve(f5(x), from=from, to=to, add=T, col="blue", ylim=ylim, n=256)
curve(f3(x), from=from, to=to, add=T, col="purple", ylim=ylim, n=256)


curve(f15(x), from=from, to=to, ylim=ylim)
curve(f3(x), from=from, to=to, add=T, col="purple", ylim=ylim, n=256)

### Composition: + normalization


# It is possible to compose these polynomials:
# - this enables the construction of polynomial "wavelets";

f = function(x, n) (f3(x)/3 + f5(x)*(2/3 - 1/n) + f7(x)/n + 2*d)
curve(f(x, 2), from=from, to=to, ylim=ylim, col="springgreen2")
for(n in 3:15) curve(f(x, n), from=from, to=to, ylim=ylim, col="red", add=T)

f = function(x, n) (f5(x)/3 + f7(x)*(2/3 - 1/n) + f9(x)/n + 2*d)
for(n in 2:15) curve(f(x, n), from=from, to=to, ylim=ylim, col="lightblue", add=T)


#####################
### Zero Baseline ###

f = function(x, n) ((f3(x) - f5(x))/3 + (f7(x) - f3(x))*(2/3 - 1/n) + (f9(x) - f11(x))/n)/2
curve(f(x, 2), from=from, to=to, ylim=ylim, col="springgreen2")
for(n in 3:15) curve(f(x, n), from=from, to=to, ylim=ylim, col="lightblue", add=T)

f = function(x, n) ((f3(x) - f7(x))/3 + (f5(x) - f7(x))*(2/3 - 1/n) + (f9(x) - f11(x))/n)/2
for(n in 2:15) curve(f(x, n), from=from, to=to, ylim=ylim, col="red", add=T)


###
id.all = 2:length(f.all)
id2 = 3
curve(f3(x), from=from, to=to, ylim=ylim, col="green")
curve(f3(f.all[[1]](x)), from=from, to=to, ylim=ylim, col=cols.all[[1]], add=T)
for(id in id2) curve(f3(f.all[[id]](x)), from=from, to=to, ylim=ylim, col=cols.all[[id]], add=T)
lines(c(-2,2), c(-2,2), col="red")


id.all = 2:length(f.all)
id2 = 3
curve(f5(x), from=from, to=to, ylim=ylim, col="green")
curve(f5(f.all[[2]](x)), from=from, to=to, ylim=ylim, col=cols.all[[1]], add=T)
for(id in id2) curve(f5(f.all[[id]](x)), from=from, to=to, ylim=ylim, col=cols.all[[id]], add=T)
lines(c(-2,2), c(-2,2), col="red")


### Recursive
id = 2
curve(f.all[[id]](x), from=from, to=to, ylim=ylim, col="green")
curve(f.all[[id]](f.all[[id]](x)), from=from, to=to, ylim=ylim, col=cols.all[[4]], add=T)
curve(f.all[[id]](f.all[[id]](f.all[[id]](x))), from=from, to=to, ylim=ylim, col=cols.all[[5]], add=T)
lines(c(-2,2), c(-2,2), col="red")

id = 2
x = seq(-2, 2, by=0.0001)
plot.init(x.lim=c(-2, 2), y.lim=c(-3, 3), c("x", "f"))
lines(x, f.all[[id]](x), col="green")
lines(x, f.all[[id]](f.all[[id]](x)), col=cols.all[[4]])
lines(x, f.all[[id]](f.all[[id]](f.all[[id]](x))), col=cols.all[[5]])
lines(c(-2,2), c(-2,2), col="red")

#####################
#####################

# install.packages("wavelets")

library(wavelets)

data(ecg)

x.dwt = dwt(ecg)
plot(x.dwt)

###
start = c(0)
end = c()
seq.seq = c(start, seq(-2, 0, by=0.1), end)
f5.flt = - f5(seq.seq) / 4
f7.flt = - f7(seq.seq) / 4
f9.flt = - f9(seq.seq) / 4
f11.flt = - f11(seq.seq) / 4
f13.flt = - f13(seq.seq) / 4
f15.flt = - f15(seq.seq) / 4
f.flt = list(f5.flt, f7.flt, f9.flt, f11.flt, f13.flt, f15.flt)

i = 5
flt = f.flt[[i]] - f.flt[[i-1]]
x.dwt = dwt(ecg, filter=flt)
plot(x.dwt)


