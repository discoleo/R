#########################
##
## Mathematical Models
##  of Tumour Growth
##
## Leonard Mada
##
## draft v.0.1f


### Introduction

# Animal models are often used to study malignant processes.
# Experimental models are often based on mice or rats.
# Models include spontaneous or induced tumour models,
# transgenic tumours as well as transplanted tumours. (Ref 1)

# Measuring tumour sizes in living animals is often plagued
# by variability and measurement errors. However, the rate of growth
# of the tumour may impact the response to therapeutic interventions.

# Measurements in living animals usually comprise some area
# (of the tumour), while the volume is proportional to Area^1.5,
# justifying the use of fractional exponents in various formulas.

# 1) Y. Zhou et al. Experimental mouse models for translational
#    human cancer research. Front Immunol. 2023 Mar 10;14:1095388.
#    https://doi.org/10.3389/fimmu.2023.1095388. PMID: 36969176;


### Evaluation of Growth Rate
# - T25, T50, T75: time to reach 25%, 50% or 75% of final volume;
# - Ratios of the various Txx parameters;


### Model Types

# 1. Simple Models: NLS, ODEs;
# 2. Linear Combinations of Base-types;
# 3. Non-Linear Combinations of Base-types;
# Note:
# - does NOT cover PDEs;
# - Combinations: fare badly when evaluated with information criteria;

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


####################
####################

### Helper Functions

library(demodelr)
library(dplyr)
library(ggplot2)

### Palette
col.blue    = c("#0000FFC0", "#72A2FAA0", "#A888FAA0", "#88A8FAA0");
col.magenta = c("#F000FFC0", "#F8A2A2A0", "#F864B4A0", "#F8B4FFA0");
col.green   = c("#00FF00A0", "#90FE64A0", "#32FC90A0", "#64FC90A0");
col.brown   = c("#F8B464A0");
# c("#FF0000A0", "#FC2490A0", "#FC9064A0", "#FE8090A0")


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

solve.eq = function(eq, params, Vx = Vmax / 2, Vmax = 2, ..., lim = c(0, 100)) {
	params = c(as.list(params), Vmax=Vmax);
	eqf = function(x)
		eval(eq, c(params, t=x)) - Vx;
	uniroot(eqf, lim, ...);
}

### Plot:
curve.ref = function(col = "#F02424E2", lwd = 2, params = NULL, Vmax = 2,
		col.h = "#FF3224B0", xlab = "Time", ylab = "Growth") {
	ylim = c(0, Vmax + 0.1);
	if(is.null(params)) params = list(Vmax = Vmax, b = 1, n = 1);
	curve(eval.Exp(x, "t", MM.eq, params), col=col, lwd=lwd,
		xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab);
	abline(h = Vmax, col = col.h, lwd=lwd);
}
plot.curve = function() {
	# TODO
}


### Michaelis-Menten Models

MM.eq = expression(Vmax * t^n / (b + t^n));

# Vx = solve V - Vx == 0;
curve.MM = function(col, b = c(1, 2/5, 2, 3), lwd = 2, legend = TRUE,
		Vmax = 2, Vx = Vmax / 2, xy.txt = c(8, 0.5)) {
	doS = (Vx != 0);
	params = list(Vmax = Vmax, b = b[1], n = 1);
	if(doS) sol1 = solve.eq(MM.eq, params, Vx=Vx);
	# curve(eval.Exp(x, "t", MM.eq, params), col = col[1], xlim=xlim, ylim=ylim)
	curve.ref(params = params, col = col[1], lwd=lwd);
	if(length(b) == 1) return(); # Hack: Reference curve;
	params$b = b[2];
	if(doS) sol2 = solve.eq(MM.eq, params, Vx=Vx);
	curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[2], lwd=lwd)
	params$b = b[3];
	if(doS) sol3 = solve.eq(MM.eq, params, Vx=Vx);
	curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[3], lwd=lwd)
	params$b = b[4];
	if(doS) sol4 = solve.eq(MM.eq, params, Vx=Vx);
	curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[4], lwd=lwd)
	if(legend) {
		legend(xy.txt[1], xy.txt[2], legend = paste0("MM", seq_along(b)), fill = col);
	}
	if(doS) {
		sol = rbind(sol1, sol2, sol3, sol4);
		return(sol);
	}
}
curve.MMextn = function(b, n, col, Vx = Vmax / 2, Vmax = 2, lwd = 2, xy.legend) {
	params  = list(Vmax = Vmax, b = b[1], n = n[1]);
	doSolve = (Vx != 0);
	len = length(b);
	if(len == 0) return();
	sol = list();
	for(id in seq(len)) {
		params$b = b[id]; params$n = n[id];
		if(doSolve) sol[[id]] = solve.eq(MM.eq, params, Vmax=Vmax, Vx=Vx);
		curve(eval.Exp(x, "t", MM.eq, params), add = TRUE,
			col = col[id], lwd=lwd);
	}
	xy = xy.legend;
	legend(xy[1], xy[2], legend = paste0("MMv", seq(len)), fill = col);
	if(doSolve) {
		sol = do.call(rbind, sol);
		return(sol);
	}
}

#################

#################
### NLS Types ###
#################

### Michaelis-Menten type:
# V = Vmax * t^n / (b + t^n)

MM.eq = expression(Vmax * t^n / (b + t^n));

### Init:
xlim = c(0, 15); ylim = c(0, 2.1);
# xlim = c(0, 100);  ylim = c(0, 2.1);
# xlim = c(0, 1000); ylim = c(0, 2.1);
#
lwd = 2;


### Michaelis-Menten: Basic
Vx = 1; Vmax = 2;
curve.MM(Vx=Vx, Vmax=Vmax, col=col.blue, lwd=lwd)

### Michaelis Menten: Power-Variants
xy = c(11.5, 0.5); col = col.magenta[2:4];
b  = c(1/2, 1, 2); n = rep(1/2, 3);
sol = curve.MMextn(b=b, n=n, Vx=Vx, Vmax=Vmax, col=col, lwd=lwd, xy.legend = xy);
# Vmax / 2
abline(h = 1, col = "green", lty = 2)
abline(v = 4, col = "green", lty = 2)
print(sol)


### Note:
# Ratio (T50 - T25) / (T75 - T50):
# - Independent on b;
# - Depends only on exponent n;
# Values of Ratio:
# - n = 1:   R =  4/3;
# - n = 1/2: R = 10/9;
# - n = 3/2: R = 1.48;


### Saturated Exponential type:
# V = Vmax * (1 - exp(-b*t^n)^k)

exps.eq = expression(Vmax * (1 - exp(-b*t^n))^k)

curve.ExpS = function(b, n = 1, k = 1, col, Vmax = 2, lwd = 2,
		labels = NULL, xy.labels = NULL) {
	len = length(b);
	if(len == 0) return();
	if(length(n) == 1 && len > 1) n = rep(n, len);
	if(length(k) == 1 && len > 1) k = rep(k, len);
	params = list(Vmax = Vmax, b = b[1], n = n[1], k = k[1]);
	# Curves:
	for(id in seq(len)) {
		params$b = b[id]; params$n = n[id]; params$k = k[id];
		curve(eval.Exp(x, "t", exps.eq, params), add = TRUE,
			col = col[id], lwd=lwd);
	}
	# Legend:
	if(! is.null(labels)) {
		xy = xy.labels;
		if(is.null(xy)) xy = c(11.5, Vmax / 4); # Hardcoded!
		legend(xy[1], xy[2], labels, fill=col);
	}
}
curve.ExpSMix = function(b = c(1, 1, 1/2, 1/5), n = c(1, 1/2, 1, 1), k = 1,
		Vmax = 2, col, lwd = 2, xy.labels = c(11.5, Vmax / 4)) {
	doLegend = ! is.null(xy.labels);
	lbls = if(doLegend) {
		paste0("ES", rep(c("n","b"), each=2), 1:2);
	} else NULL;
	curve.ExpS(b=b, n=n, k=k, Vmax=Vmax, col=col, lwd=lwd,
		labels = lbls, xy.labels = xy.labels);
}

###
curve.MM(col.blue, lwd=lwd)
#
xy = c(11.5, 0.5); col = col.green;
curve.ExpSMix(Vmax = 2, col=col, xy.labels = xy)


### Variation of n:
curve.MM(col.blue, lwd=lwd)
#
col = col.green;
params = list(Vmax = 2, b = 1, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[1], lwd=lwd)
params = list(Vmax = 2, b = 1, n = 1/2, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[2], lwd=lwd)
params = list(Vmax = 2, b = 1, n = 1/3, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[3], lwd=lwd)
params = list(Vmax = 2, b = 1, n = 1/4, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[4], lwd=lwd)


### Variation of b:
curve.MM(col.blue, lwd=lwd)
#
col = col.green;
params = list(Vmax = 2, b = 1, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[1], lwd=lwd)
params = list(Vmax = 2, b = 1/2, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[2], lwd=lwd)
params = list(Vmax = 2, b = 1/3.5, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[3], lwd=lwd)
params = list(Vmax = 2, b = 1/5.5, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[4], lwd=lwd)


# Simple vs Extended
curve.ref(col.blue[1])
col = col.green;
params = list(Vmax = 2, b = 1, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[1], lwd=lwd)
params = list(Vmax = 2, b = 1, n = 1/2, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[2], lwd=lwd)
params = list(Vmax = 2, b = 1, n = 1/3, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[3], lwd=lwd)
params = list(Vmax = 2, b = 1, n = 1/4, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[4], lwd=lwd)
#
col = col.magenta;
params = list(Vmax = 2, b = 1, n = 1, k = 1)
# curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[1], lwd=lwd)
params = list(Vmax = 2, b = 1/2, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[2], lwd=lwd)
params = list(Vmax = 2, b = 1/3.5, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[3], lwd=lwd)
params = list(Vmax = 2, b = 1/5.5, n = 1, k = 1)
curve(eval.Exp(x, "t", exps.eq, params), add = T, col = col[4], lwd=lwd)
#
col.tmp = c(col.blue[1], col.magenta[-1]);
legend(10, 0.75, legend = c("MM1", paste0("ExpExt", 1:3)), fill = col.tmp);


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

