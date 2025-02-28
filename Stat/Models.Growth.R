########################################
##
## Mathematical Models of Tumour Growth
##
## Leonard Mada
##
## draft v.0.1b


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


### Eval Formulas:
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

### Plot:
curve.ref = function(col = "#F02424E2", lwd = 2, params = NULL, Vmax = 2,
		col.h = "#FF3224B0") {
	ylim = c(0, Vmax + 0.1);
	if(is.null(params)) params = list(Vmax = Vmax, b = 1, n = 1);
	curve(eval.Exp(x, "t", MM.eq, params), col = col, lwd=lwd,
		xlim=xlim, ylim=ylim, ylab = "Growth");
	abline(h = Vmax, col = col.h, lwd=lwd);
}
plot.curve = function() {
	# TODO
}

####################
####################

### Introduction

# Animal models are often used to study malignant processes.
# Experimental models are often based on mice or rats.
# Models include spontaneous or induced tumour models,
# transgenic tumours as well as transplanted tumours.

# Measuring tumour sizes in living animals is often plagued
# by variability and measurement errors. However, the rate of growth
# of the tumour may impact the response to therapeutic interventions.

# 1) Y. Zhou et al. Experimental mouse models for translational
#    human cancer research. Front Immunol. 2023 Mar 10;14:1095388.
#    https://doi.org/10.3389/fimmu.2023.1095388. PMID: 36969176;


### Effects of Growth Rate
# T50, T75: time to reach 50% or 75% of final volume;


### Model Types

### NLS:

### Michaelis-Menten type:
# V = Vmax * t^n / (b + t^n)
### Saturated Exponential:
# V = Vmax * (1 - exp(-k*t))
### Exponential Fractions:
# V = Vmax * (1/(k + exp(- b*t)) - 1/(k + 1)) * k*(k+1);
### Atan type:
# V = Vmax * atan(k*x^p) * 2/pi;


### ODE:

# dVdt = k * V * (Vmax - V)
# dVdt = a * t^p / (t^n + 1)^k

### Integrals

# V = a * I( x^p / (x^n + 1)^k ), x on [0, t]


#################

#################
### NLS Types ###
#################

col.blue    = c("#0000FFC0", "#72A2FAA0", "#A888FAA0", "#88A8FAA0");
col.magenta = c("#F000FFC0", "#F8A2A2A0", "#F864B4A0", "#F8B4FFA0");
col.green   = c("#00FF00A0", "#32FC90A0", "#64FC90A0", "#90FE64A0");
col.brown   = c("#F8B464A0");


### Michaelis-Menten type:
# V = Vmax * t^n / (b + t^n)

MM.eq = expression(Vmax * t^n / (b + t^n));

curve.MM = function(col, lwd = 2, legend = TRUE, Vmax = 2) {
	params = list(Vmax = Vmax, b = 1, n = 1)
	# curve(eval.Exp(x, "t", MM.eq, params), col = col[1], xlim=xlim, ylim=ylim)
	curve.ref(params = params, col = col[1], lwd=lwd);
	params = list(Vmax = Vmax, b = 2/5, n = 1)
	curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[2], lwd=lwd)
	params = list(Vmax = Vmax, b = 2, n = 1)
	curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[3], lwd=lwd)
	params = list(Vmax = Vmax, b = 3, n = 1)
	curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[4], lwd=lwd)
	if(legend) {
		legend(12, 0.5, legend = paste0("MM", 1:4), fill = col);
	}
}

xlim = c(0, 15); ylim = c(0, 2.1);
# xlim = c(0, 100);  ylim = c(0, 2.1);
# xlim = c(0, 1000); ylim = c(0, 2.1);
#
lwd = 2;

curve.MM(col.blue, lwd=lwd)
#
col = col.magenta
params = list(Vmax = 2, b = 1/2, n = 1/2)
curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[2], lwd=lwd)
params = list(Vmax = 2, b = 1, n = 1/2)
curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[3], lwd=lwd)
params = list(Vmax = 2, b = 2, n = 1/2)
curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[4], lwd=lwd)


### Saturated Exponential type:
# V = Vmax * (1 - exp(-b*t^n)^k)

exps.eq = expression(Vmax * (1 - exp(-b*t^n))^k)

curve.MM(col.blue, lwd=lwd)
#
col = col.green;
params = list(Vmax = 2, b = 1, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[1], lwd=lwd)
params = list(Vmax = 2, b = 1, n = 1/2, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[2], lwd=lwd)
# Param b:
params = list(Vmax = 2, b = 1/2, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[3], lwd=lwd)
params = list(Vmax = 2, b = 1/5, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[4], lwd=lwd)


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


###################
###################

### ODE

### Logistic Growth

log.eq = V ~ k * V^p * (Vmax-V)

# Ref:
params = list(Vmax = 2, b = 1, n = 1)
curve(eval.Exp(x, "t", MM.eq, params), col = "#909090B0",
	xlim=xlim, ylim=ylim, ylab = "Growth")
#
params = list(Vmax = 2, k = 1, p = 1)
sol = rk4(log.eq, c(V=0.1), param = params, deltaT=0.1, n_steps = 160)
lines(sol, col = col[2])
params = list(Vmax = 2, k = 1, p = 1/3)
sol = rk4(log.eq, c(V=0.1), param = params, deltaT=0.1, n_steps = 160)
lines(sol, col = col[3])
params = list(Vmax = 2, k = 1, p = 1.5)
sol = rk4(log.eq, c(V=0.1), param = params, deltaT=0.1, n_steps = 160)
lines(sol, col = col[3])


log2.eq = V ~ k * V^p * (Vmax^(1/2) - V^(1/2))
#
params = list(Vmax = 2, k = 1, p = 1)
sol = rk4(log2.eq, c(V=0.1), param = params, deltaT=0.1, n_steps = 160)
lines(sol, col = "blue")

# TODO

