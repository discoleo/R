########################################
##
## Mathematical Models of Tumour Growth
##
## Leonard Mada
##
## draft v.0.1a


### Models:
# 1. Simple Models: NLS, ODEs;
# 2. Linear Combinations of Base-types;
# 3. Non-Linear Combinations of Base-types;
# Note:
# - does NOT cover PDEs;


####################

### Helper Functions

library(demodelr)
library(dplyr)
library(ggplot2)


eval.Fun = function(x, nm, FUN, params) {
	sapply(x, function(x) {
		params[[nm]] = x;
		do.call(FUN, params);
	})
}
eval.Exp = function(x, nm, e, params) {
	sapply(x, function(x) {
		params[[nm]] = x;
		eval(e, params);
	})
}
plot.curve = function()

####################
####################

### Model Types

### NLS:

### Michaelis-Menten type:
# V = Vmax * t / (1 + b*t)
### Saturated Exponential:
# V = Vmax * (1 - exp(-b*t))
### Exponential Fractions:
# V = Vmax * (1/(k + exp(- b*t)) - 1/(k + 1)) * k*(k+1);


### ODE:

# dVdt = k * V * (Vmax - V)
# dVdt = a * t^p / (t^n + 1)^k

### Integrals

# V = a * I( x^p / (x^n + 1)^k ), x on [0, t]


#################

#################
### NLS Types ###
#################

### Michaelis-Menten type:
# V = Vmax * t^n / (1 + b*t^n)

xlim = c(0, 15); ylim = c(0, 2.1);
#
MM.eq = expression(Vmax * t^n / (b + t^n))
col = c("#0000FFA0", "#2496F2A0", "#9664F8A0")
params = list(Vmax = 2, b = 1, n = 1)
curve(eval.Exp(x, "t", MM.eq, params), col = col[1], xlim=xlim, ylim=ylim)
params = list(Vmax = 2, b = 2/5, n = 1)
curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[2])
params = list(Vmax = 2, b = 3, n = 1)
curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[3])


### Saturated Exponential type:
# V = Vmax * (1 - exp(-b*t^n)^k)

exps.eq = expression(Vmax * (1 - exp(-b*t^n))^k)
col = c("#FF0000A0", "#FC2490A0", "#FC9064A0", "#FE8090A0")
#
params = list(Vmax = 2, b = 1, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[1])
params = list(Vmax = 2, b = 1, n = 1/2, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[2])
params = list(Vmax = 2, b = 1/2, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[3])
params = list(Vmax = 2, b = 1/5, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[4])


#################
#################

### Integrals

### I( x^p / (x^n + 1)^k )
frI = function(up, Vmax, p, n) {
	FUN = \(x) x^p / (x^n + 1);
	tmp = integrate(FUN, lower = 0, upper = up, rel.tol = 1E-8)$value;
	Vmax * sin(pi*(p+1)/n) * n / pi * tmp;
}
col = c("#0000FFA0", "#2496F2A0", "#9664F8A0")
# Ref:
params = list(Vmax = 2, b = 1, n = 1)
curve(eval.Exp(x, "t", MM.eq, params), col = "#A0A0A0B0",
	xlim=xlim, ylim=ylim, ylab = "Growth")
#
params = list(Vmax = 2, p = 0, n = 2)
curve(eval.Fun(x, "up", frI, params), col = col[1], add = T)
params = list(Vmax = 2, p = 1/2, n = 2)
curve(eval.Fun(x, "up", frI, params), col = col[2], add = T)
params = list(Vmax = 2, p = 1/3, n = 3)
curve(eval.Fun(x, "up", frI, params), col = col[3], add = T)
params = list(Vmax = 2, p = 4/3, n = 3)
curve(eval.Fun(x, "up", frI, params), col = col[3], add = T)

