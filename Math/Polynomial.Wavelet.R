
### Polynomial Wavelet Analysis
###
### (C) Leonard Mada
### draft 0.1: 2020-03-18
### pre-01: 2020-02;
###
### Cardan-type Polynomials:
### theory developed during 2018-2019

c = 1
d = 0

### "Quadratic" / Cardan-type Polynomials
f3  = function(x) x^3 - 3*c*x - 2*d
f5  = function(x) x^5 - 5*c*x^3 + 5*c^2*x - 2*d
f7  = function(x) x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d
f9  = function(x) x^9 - 9*c*x^7 + 27*c^2*x^5 - 30*c^3*x^3 + 9*c^4*x - 2*d
f11 = function(x) x^11 - 11*c*x^9 + 44*c^2*x^7 - 77*c^3*x^5 + 55*c^4*x^3 - 11*c^5*x - 2*d
f13 = function(x) x^13 - 13*c*x^11 + 65*c^2*x^9 - 156*c^3*x^7 + 182*c^4*x^5 - 91*c^5*x^3 + 13*c^6*x - 2*d
f15 = function(x) x^15 - 15*c*x^13 + 90*c^2*x^11 - 275*c^3*x^9 + 450*c^4*x^7 - 378*c^5*x^5 + 140*c^6*x^3 - 15*c^7*x - 2*d


### Properies

# These polynomials have very interesting properties:
# 1.) P(0) = -2*d;
# 2.) for c = 1:
# all values of P(x) are between:[-2*d - 2, -2*d + 2] for x in [-2, 2];

### Visualization:
from = -2.005
to = 2.005
ylim = NULL # c(-15, 10)

curve(f15(x), from=from, to=to, ylim=ylim)
curve(f13(x), from=from, to=to, ylim=ylim, add=T, col="springgreen2")
curve(f11(x), from=from, to=to, ylim=ylim, add=T, col="lightseagreen")
curve(f9(x), from=from, to=to, add=T, col="green", ylim=ylim)
curve(f7(x), from=from, to=to, add=T, col="red", ylim=ylim)
curve(f5(x), from=from, to=to, add=T, col="blue", ylim=ylim)
curve(f3(x), from=from, to=to, add=T, col="purple", ylim=ylim)

### Composition: + normalization

# It is possible to compose these polynomials:
# - this enables the construction of polynomial "wavelets";

f = function(x, n) (f3(x)/3 + f5(x)*(2/3 - 1/n) + f7(x)/n + 2 + 2*d)
curve(f(x, 2), from=from, to=to, ylim=ylim, col="springgreen2")
for(n in 3:15) curve(f(x, n), from=from, to=to, ylim=ylim, col="red", add=T)

f = function(x, n) (f5(x)/3 + f7(x)*(2/3 - 1/n) + f9(x)/n + 2 + 2*d)
for(n in 2:15) curve(f(x, n), from=from, to=to, ylim=ylim, col="lightblue", add=T)


