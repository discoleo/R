#########################
##
## Mathematical Models
##  of Tumour Growth
##
## Leonard Mada
##
## draft v.0.1t


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
# - Ratios of differences of Txx, e.g.:
#   R = (T75 - T50) / (T50 - T25);

# Note:
# - Estimating Tnn requires accurate knowledge of Vmax;
# - Estimated Vmax may be inaccurate, shifting Tn1 / Tn2
#   significantly off the correct value;
# => R = (T75/T50 - 1) / (1 - T25/T50) should be more stable;
# - if Vmax is underestimated, T75/T50 is lower than the correct value,
#   but T25/T50 is larger, which makes (1 - T25/T50) smaller,
#   ultimately balancing the ratio;


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
# dVdt = a * t^p * exp(-k * t^n)


### Integrals

# V = a * I( x^p / (x^n + 1)^k ), x on [0, t]
# - can be reformulated as an ODE (see above);
# - where a = scaled Vmax;


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
col.grey    = c("#909090B0");
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
# xlim = defined globally;
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
curve.MM = function(col, b = c(1, 1/2, 2, 3), n = 1, lwd = 2, legend = TRUE,
		Vmax = 2, Vx = Vmax / 2, xy.txt = c(8, 0.5)) {
	doS = (Vx != 0);
	params = list(Vmax = Vmax, b = b[1], n = n[1]);
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
		curve(eval.Exp(x, "t", MM.eq, params), add = TRUE,
			col = col[id], lwd=lwd);
		if(doSolve) sol[[id]] = solve.eq(MM.eq, params, Vmax=Vmax, Vx=Vx);
	}
	xy = xy.legend;
	legend(xy[1], xy[2], legend = paste0("MMv", seq(len)), fill = col);
	if(doSolve) {
		sol = do.call(rbind, sol);
		return(sol);
	}
}


### Saturated Exponential Models

exps.eq = expression(Vmax * (1 - exp(-b*t^n))^k)

curve.ExpS = function(b = 1, n = 1, k = 1,
		Vx = Vmax / 2, Vmax = 2, col = 1, lwd = 2,
		labels = NULL, xy.labels = NULL) {
	len = c(length(b), length(n), length(k));
	if(any(len == 0)) return();
	len = max(len);
	# Params:
	if(length(b) == 1 && len > 1) b = rep(b, len);
	if(length(n) == 1 && len > 1) n = rep(n, len);
	if(length(k) == 1 && len > 1) k = rep(k, len);
	params = list(Vmax = Vmax, b = b[1], n = n[1], k = k[1]);
	sol    = list();
	if(length(col) == 1 && len > 1) col = rep(col, len);
	doSolve = (Vx != 0);
	# Curves:
	for(id in seq(len)) {
		params$b = b[id]; params$n = n[id]; params$k = k[id];
		curve(eval.Exp(x, "t", exps.eq, params), add = TRUE,
			col = col[id], lwd=lwd);
		if(doSolve) sol[[id]] = solve.eq(exps.eq, params, Vmax=Vmax, Vx=Vx);
	}
	# Legend:
	if(! is.null(labels)) {
		xy = xy.labels;
		if(is.null(xy)) xy = c(11.5, Vmax / 4); # Hardcoded!
		legend(xy[1], xy[2], labels, fill=col);
	}
	if(doSolve) {
		sol = do.call(rbind, sol);
		return(sol);
	}
}

# Mix of various Saturated Exponentials
curve.ExpSMix = function(b = c(1, 1, 1/2, 1/5), n = c(1, 1/2, 1, 1), k = 1,
		Vx = Vmax / 2, Vmax = 2, col = 1, lwd = 2,
		xy.labels = c(11.5, Vmax / 4)) {
	doLegend = ! is.null(xy.labels);
	lbls = if(doLegend) {
		paste0("ES", rep(c("n","b"), each=2), 1:2);
	} else NULL;
	curve.ExpS(b=b, n=n, k=k, Vx=Vx, Vmax=Vmax, col=col, lwd=lwd,
		labels = lbls, xy.labels = xy.labels);
}


### ATAN Model

atan.eq = expression(Vmax * (atan(b * t^n) * 2/pi)^k)

curve.Atan = function(b = 1, n = 1, k = 1,
		Vx = Vmax / 2, Vmax = 2,
		col = 1, lwd = 2, labels = NULL, xy.labels = NULL) {
	len = c(length(b), length(n), length(k));
	if(any(len == 0)) return();
	len = max(len);
	# Params:
	if(length(b) == 1 && len > 1) b = rep(b, len);
	if(length(n) == 1 && len > 1) n = rep(n, len);
	if(length(k) == 1 && len > 1) k = rep(k, len);
	params = list(Vmax=Vmax, b = b[1], n = n[1], k = k[1]);
	sol    = list();
	if(length(col) == 1 && len > 1) col = rep(col, len);
	doSolve = (Vx != 0);
	for(id in seq(len)) {
		params$b = b[id]; params$n = n[id]; params$k = k[id];
		curve(eval.Exp(x, "t", atan.eq, params), add = TRUE,
			col = col[id], lwd=lwd);
		if(doSolve) sol[[id]] = solve.eq(atan.eq, params, Vmax=Vmax, Vx=Vx);
	}
	# Legend:
	if(! is.null(xy.labels)) {
		if(is.null(labels)) {
			labels = paste0("AT", seq(len));
		}
		xy = xy.labels;
		legend(xy[1], xy[2], labels, fill=col);
	}
	if(doSolve) {
		sol = do.call(rbind, sol);
		return(sol);
	}
}

### Integrals

### I( x^p / (x^n + 1)^k )
frI = function(up, Vmax, p, n, k=1) {
	FUN = \(x) x^p / (x^n + 1)^k;
	tmp = integrate(FUN, lower = 0, upper = up, rel.tol = 1E-8)$value;
	# Vmax * sin(pi*(p+1)/n) * n / pi * tmp; # k == 1;
	scale = gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n;
	Vmax * tmp /scale;
}

curve.IntFr = function(p = c(0, 1/3, 1/2, 4/3), n = 3, k = 1,
		col = 1, lwd = 2, Vx = Vmax / 2, Vmax = 2,
		labels = NULL, xy.labels = NULL) {
	len = c(length(p), length(n), length(k));
	if(any(len == 0)) return();
	if(Vx != 0) doSolve = TRUE;
	len = max(len);
	sol = list();
	eq  = expression(frI(up=t, Vmax=Vmax, p=p, n=n, k=k));
	# Params:
	if(length(p) == 1 && len > 1) b = rep(b, len);
	if(length(n) == 1 && len > 1) n = rep(n, len);
	if(length(k) == 1 && len > 1) k = rep(k, len);
	params = list(Vmax=Vmax, p = p[1], n = n[1], k = k[1]);
	for(id in seq(len)) {
		params$p = p[id]; params$n = n[id];
		curve(eval.Fun(x, "up", frI, params), col = col[id], lwd=lwd, add = TRUE);
		if(doSolve) sol[[id]] = solve.eq(eq, params, Vmax=Vmax, Vx=Vx);
	}
	# Legend:
	if(! is.null(labels)) {
		xy = xy.labels;
		if(is.null(xy)) xy = c(11.5, Vmax / 4); # Hardcoded!
		legend(xy[1], xy[2], labels, fill=col);
	}
	if(doSolve) {
		sol = do.call(rbind, sol);
		return(sol);
	}
}

### ODE

### Logistic Growth

log.eq = V ~ k * V^p * (Vmax^n - V^n)^m

curve.odeLog = function(k = 1, p = 1, m = 1, n = 1,
		Vmax = 2, V0 = 0.1, col = 1, lwd = 2,
		labels = NULL, xy.labels = NULL,
		n_steps = 160, dt = 0.1) {
	len = c(length(k), length(p), length(n));
	if(any(len == 0)) return(data.frame(t = numeric(0), V = numeric(0)));
	len = max(len);
	# Parameters:
	if(length(k) == 1 && len > 1) k = rep(k, len);
	if(length(p) == 1 && len > 1) p = rep(p, len);
	if(length(n) == 1 && len > 1) n = rep(n, len);
	params = list(Vmax = Vmax, k=k[1], p=p[1], m=m[1], n=n[1]);
	V0 = c(V = V0);
	for(id in seq(len)) {
		params$k = k[id]; params$p = p[id]; params$n = n[id];
		sol = rk4(log.eq, V0, parameters = params, deltaT=dt, n_steps = n_steps);
		lines(sol, col = col[id], lwd=lwd);
	}
	# Legend:
	if(! is.null(xy.labels)) {
		if(is.null(labels)) {
			labels = paste0("OL", seq(len));
		} else if(length(labels) == 1 && len > 1) {
			labels = paste0(labels, seq(len));
		}
		xy = xy.labels;
		legend(xy[1], xy[2], labels, fill=col);
	}
}
# ODE: Logistic Eq with Initial conditions
# - Effect of Initial condition;
curve.odeLogI0 = function(p = c(1, 0.6, 1.25), V0 = c(0.1, 0.005),
		k = 1, n = 1, lwd = 2,
		xy1 = c(8, 0.5), xy2 = c(11.5, 0.5), adj.xy2 = c(0.25, 0.05)) {
	xy1 = xy1 + V0[1]; xy2 = xy2 + V0[1];
	curve.ref(col.grey[1], lwd=lwd);
	#
	col = col.magenta; lbls = "OL";
	curve.odeLog(k = k, p = p, n = n, V0=V0[1], col=col, lwd=lwd,
		xy=xy1, labels = lbls);
	abline(h = V0, lty = 2, col = "#B06464B2")
	# Initial condition
	col = col.blue; lbls = "OL0";
	curve.odeLog(k = k, p = p, n = n, V0=V0[2], col=col, lwd=lwd,
		xy=xy2, labels = lbls);
	xy2 = xy2 + adj.xy2;
	lblV0 = paste0("V0 = ", V0[2]);
	text(xy2[1], xy2[2], lblV0, col = col[1], adj = c(0, 0));
}

# Variation of Power n
curve.odeLogPn = function(n = c(1, 1/2), p = c(1,1/3,1.5), k = 1,
		V0 = 0.1, Vmax = 2,
		lwd = 2, xy1 = c(8, 0.5), xy2 = c(11.5, 0.5), xy.adj = c(0.5, 0.05)) {
	curve.ref(col.grey[1], lwd=lwd, Vmax=Vmax); # Ref
	xy1 = xy1 + V0; xy2 = xy2 + V0;
	#
	col = col.green;
	curve.odeLog(k = k, p = p, n = n[1], V0=V0, Vmax=Vmax,
		col=col, lwd=lwd, xy=xy1, labels = "OL");
	# Variant: n
	col = col.magenta;
	curve.odeLog(k = k, p = p, n = n[2], V0=V0, Vmax=Vmax,
		col=col, lwd=lwd, xy=xy2, labels = "OLn");
	xy2  = xy2 + xy.adj;
	lbln = paste0("n = ", n[2]);
	text(xy2[1], xy2[2], lbln, col = col[1], adj = c(0, 0))
	#
	abline(h = V0, lty = 2, col = "#B06464B2")
}


#################

#################
### NLS Types ###
#################

### Init:
xlim = c(0, 15); ylim = c(0, 2.1);
lwd  = 2;


### Michaelis-Menten type:
# V = Vmax * t^n / (b + t^n)

# Fully Generalized Equation:
# V = Vmax * t^(n*k) / (b + t^n)^k

MM.eq = expression(Vmax * t^n / (b + t^n));

### Alternative Init:
xlim = c(0, 15); ylim = c(0, 2.1);
# xlim = c(0, 100);  ylim = c(0, 2.1);
# xlim = c(0, 1000); ylim = c(0, 2.1);


### Michaelis-Menten: Basic
Vx = 1; Vmax = 2;
curve.MM(Vx=Vx, Vmax=Vmax, col=col.blue, lwd=lwd)

### Michaelis Menten: Power-Variants
xy = c(11.5, 0.5); col = col.magenta[2:4];
b  = c(1/2, 1, 2); n = rep(1/2, 3);
sol = curve.MMextn(b=b, n=n, Vx=Vx, Vmax=Vmax, col=col, lwd=lwd, xy.legend = xy);
# h = Vmax / 2
abline(h = 1, col = "green", lty = 2)
abline(v = 4, col = "green", lty = 2)
print(sol)


### Note:
# Ratios:
#   R  = (T75 - T50) / (T50 - T25) or
#   Rw = (T75 - T25) / (T75 - T50):
# - Independent on b;
# - Depend only on exponent n;
# Values of Ratio:
# - n = 1:   R = 3; Rw =  4/3;
# - n = 1/2: R = 9; Rw = 10/9;
# - n = 3/2: R = 3^(2/3); Rw = 1.48;
# - Conjecture: R = 3^(1/n);


#########################

### Saturated Exponential
# V = Vmax * (1 - exp(-b*t^n)^k)

exps.eq = expression(Vmax * (1 - exp(-b*t^n))^k)


### Basic: Mixed Variants
curve.MM(col=col.blue, lwd=lwd)
#
xy = c(11.5, 0.5); col = col.green;
curve.ExpSMix(Vmax = 2, col=col, xy.labels = xy)


### Exploration: Exponent k
# Magenta curve approaches green curve as time progresses;
curve.MM(b = 1, col.blue, lwd=lwd)
xy1 = c(8, 0.5); xy2 = c(11.5, 0.5);
curve.ExpSMix(Vmax = 2, col=col.green, xy.labels = xy1)
curve.ExpSMix(Vmax = 2, col=col.magenta, xy.labels = xy2, k = 4/3)
text(xy2[1] + 0.5, xy[2] + 0.05, "k = 4/3", col = col.magenta[1], adj = c(0, 0))


### Variation of n:
curve.MM(col=col.blue, lwd=lwd)
#
xy = c(11.5, 0.5); col = col.green;
n = c(1, 1/2, 1/3, 1/4); b = rep(1, 4);
lbls = paste0("ESn", 1:4);
curve.ExpS(n=n, b=b, k = 1, Vmax = 2, col=col, lwd=lwd, labels = lbls, xy.labels = xy)


### Variation of b:
curve.MM(col=col.blue, lwd=lwd)
#
xy = c(11.5, 0.5); col = col.green;
n = rep(1, 4); b = c(1, 1/2, 1/3.5, 1/5.5);
lbls = paste0("ESb", 1:4);
curve.ExpS(n=n, b=b, k = 1, Vmax = 2, col=col, lwd=lwd, labels = lbls, xy.labels = xy)


### Simple vs Power
# Green & Magenta curves separate much later;
curve.ref(col=col.blue[1], lwd=lwd)
xy1 = c(8, 0.5); xy2 = c(11.5, 0.5);
col = col.green; lbls = paste0("ESb", 1:4);
b = 1 / c(1.5, 2, 3.5, 5.5); n = rep(1, 4);
curve.ExpS(n=n, b=b, k = 1, Vmax = 2, col=col, lwd=lwd, labels = lbls, xy.labels = xy1)
col = col.magenta; lbls = paste0("ESnb", 1:4);
b = 1 / c(1.5, 2, 3.5, 5.5); n = 10/11;
curve.ExpS(n=n, b=b, k = 1, Vmax = 2, col=col, lwd=lwd, labels = lbls, xy.labels = xy2)
text(xy2[1] + 0.5, xy[2] + 0.05, "n = 0.91", col = col.magenta[1], adj = c(0, 0))


### Power vs Extended (2)
# Note: curves are not matched;
curve.ref(col=col.blue[1])
xy1 = c(8, 0.5); xy2 = c(11.5, 0.5);
col = col.green;
n = 1 / c(1, 2, 3, 4); b = rep(1, 4);
lbls = paste0("ESn", 1:4);
curve.ExpS(n=n, b=b, k = 1, Vmax = 2, col=col, lwd=lwd, labels = lbls, xy.labels = xy1)
#
col = col.magenta;
b = 1 - 1 / c(12, 3, 2, 1.5); n = 2/3;
# b = 1 - 1 / c(15, 4, 2.5, 1.75); n = 2/3;
# n = 1 / c(1, 2, 3, 4); b = rep(1, 4);
lbls = paste0("ESnb", 1:4);
curve.ExpS(n=n, b=b, k = 1, Vmax = 2, col=col, lwd=lwd, labels = lbls, xy.labels = xy2)


### Note:
# Ratios:
#   R  = (T75 - T50) / (T50 - T25) or
#   Rw = (T75 - T25) / (T75 - T50):
# - Independent on b;
# - Depend only on exponent n;
# Values of Ratio:
# - n = 1:   R = log(2) / (log(3) - log(2)); # 1.70
# - n = 1/2: R = 3*log(2)^2 / (log(2)^2 - log(4/3)^2); # 3.62
# - n = 3/2: # R = 1.32;
#   R = log(2)^(2/3) * (2^(2/3) - 1) / (log(2)^(2/3) - log(4/3)^(2/3));


#################
### ATAN Type ###

# V = Vmax * (atan(b * t^n) * 2/pi)^k

atan.eq = expression(Vmax * (atan(b * t^n) * 2/pi)^k)


### Basic:
# xlim = c(0, 50)
Vx = 1; n = 1;
curve.MM(col.blue, lwd=lwd, n=n, Vx=Vx)
#
xy = c(11.5, 0.5); col = col.green;
b  = c(2, 1, 1/2, 1/3);
curve.Atan(b=b, n=n, Vx=Vx, Vmax = 2, col=col, lwd=lwd, xy.labels = xy)

# Note:
# - Part 1: Atan resembles MM for a significant period of time;
# - Part 2: Atan grows faster for a little bit longer:
#   => MM flattens (saturates) slightly earlier than Atan;
# - Part 3: during the later saturation phase,
#   MM grows again a little bit faster than Atan (but only during phase 3,
#   when Atan suffers a much "stronger" saturation); 


### Variation of n:
curve.MM(col.blue, lwd=lwd)
#
xy = c(11.5, 0.5); col = col.green;
n  = c(2, 1, 1/2, 1/3); b = 1;
curve.Atan(b=b, n=n, Vmax = 2, col=col, xy.labels = xy)
# all curves pass through (1,1) (when b = 1);


### Simple vs Power
# Green & Magenta curves separate after some delay;
curve.ref(col.blue[1], lwd=lwd)
xy1 = c(8, 0.5); xy2 = c(11.5, 0.5);
col = col.green; lbls = paste0("ATb", 1:4);
b = 1 / c(1, 0.75, 2, 3); n = 1;
curve.Atan(n=n, b=b, k = 1, Vmax = 2, col=col, lwd=lwd, labels = lbls, xy.labels = xy1)
col = col.magenta; lbls = paste0("ATnb", 1:4);
n = 5/7; # n = 4/5; # overlap with MM;
curve.Atan(n=n, b=b, k = 1, Vmax = 2, col=col, lwd=lwd, labels = lbls, xy.labels = xy2)
text(xy2[1] + 0.5, xy[2] + 0.05, "n = 0.7", col = col.magenta[1], adj = c(0, 0))


### Note:
# Ratios:
#   R  = (T75 - T50) / (T50 - T25) or
#   Rw = (T75 - T25) / (T75 - T50):
# - Independent on b;
# - Depend only on exponent n;
# Values of Ratio:
# - n = 1:   R = sqrt(2) + 1;
#            (tan(3*pi/8) - 1) / (1 - tan(pi/8));
# - n = 1/2: R = tan(3*pi/8)^2;
#            (tan(3*pi/8)^2 - 1) / (1 - tan(pi/8)^2);
# - n = 3/2: R = tan(3*pi/8)^(2/3);
#            (tan(3*pi/8)^(2/3) - 1) / (1 - tan(pi/8)^(2/3))


#################
#################

### Integrals

### I( x^p / (x^n + 1)^k )

# Note:
# - can be reformulated as an ODE;
# - but the integral is directly computable in R;

### Basic: Mixed Variants
curve.MM(col.blue, lwd=lwd, Vx=0)
#
p = c(0, 0.3, 2/3, 1.1); n = c(2,2, 3,3);
lbls = paste0("IF", rep(c("p","n"), each=2), 1:2);
curve.IntFr(p=p, n=n, col = col.green, lwd=lwd, Vmax=Vmax, labels=lbls)


### Simple vs Power
curve.ref(col.blue[1], lwd=lwd)
xy1 = c(8, 0.5); xy2 = c(11.5, 0.5);
#
Vx = 1;
p = c(0,1/2,1,4/3); n = 2.75;
lbls = paste0("IFp", 1:4);
curve.IntFr(p=p, n=n, col = col.green, lwd=lwd, Vx=Vx, Vmax=Vmax, labels=lbls, xy=xy1)
p = c(0,1/2,1,4/3); n = 2.9;
lbls = paste0("IFnp", 1:4);
curve.IntFr(p=p, n=n, col = col.magenta, lwd=lwd, Vx=Vx, Vmax=Vmax, labels=lbls, xy=xy2)
text(xy2[1] + 0.5, xy2[2] + 0.05, "n = 2.9", col = col.magenta[1], adj = c(0, 0))


###################
###################

### ODE

### Logistic Growth

log.eq = V ~ k * V^p * (Vmax^n - V^n)^m

# Specified Variant:
log2.eq = V ~ k * V^p * (Vmax^(1/2) - V^(1/2))

# Simple:
# log.eq = V ~ k * V^p * (Vmax - V);

# Fully Generalized:
# V ~ k * (P1(V) - P1(0))^m1 * (P2(Vmax) - P2(V))^m2;
# where P1, P2 = 2 functions;


### Variation of Power n

curve.odeLogPn(lwd=lwd)


### Initial Conditions:
# Note: highly dependent on V0!

curve.odeLogI0(n = 0.75, lwd=lwd)



# TODO: more;

