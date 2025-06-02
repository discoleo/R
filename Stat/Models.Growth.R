#########################
##
## Mathematical Models
##  of Tumour Growth
##
## Leonard Mada
##
## draft v.0.2c


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
# V = Vmax * t^(n*m) / (b + t^n)^m
### Saturated Exponential:
# V = Vmax * (1 - exp(-k*t^n))^m
### Exponential Fractions:
# V = Vmax * (1 - exp(- b1*t^n1))^m1 / (1 + 1/k * exp(- b2*t^n2))^m2
# [old] V = Vmax * (1/(k + exp(- b*t)) - 1/(k + 1)) * k*(k+1);
### Atan type:
# V = Vmax * atan(k*t^p) * 2/pi;
### Gompertz type:
# V = Vmax * exp(- k1 * exp(- k2*t^n))
# V = Vmax * exp(- k1 * exp(- k2*t^n)) * atan(k2 * t^m)*2/pi


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

# Note: moved to file:
# Models.Growth.Helper.R

source("Models.Growth.Helper.R")


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


### Michaelis-Menten: Basic w. Delta(Tnn)
curve.MM(Vx=Vx, Vmax=Vmax, col=col.blue, lwd=lwd)
### Tnn for curve with b = 1;
# T50: h = Vmax / 2
for(val in c(0.5,1,1.5)) abline(h = val, col = "green", lty = 2);
for(val in c(1/3,1,3)) abline(v = val, col = "green", lty = 2);
#
lines(c(1/3, 1), c(0.5,0.5), col = "green", lwd=lwd)
lines(c(3, 1), c(1.5,1.5), col = "green", lwd=lwd)


### Michaelis Menten: Power-Variants
curve.MM(Vx=Vx, Vmax=Vmax, col=col.blue, lwd=lwd)
xy = c(11.5, 0.5); col = col.magenta[2:4];
b  = c(1/2, 1, 2); n = rep(1/2, 3);
sol = curve.MMextn(b=b, n=n, Vx=Vx, Vmax=Vmax, col=col, lwd=lwd, xy.legend = xy);
xy2 = xy + c(0.5, 0.05);
text(xy2[1], xy2[2], "n = 0.5", col = col.magenta[1], adj = c(0, 0))
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
text(xy2[1] + 0.5, xy2[2] + 0.05, "k = 4/3", col = col.magenta[1], adj = c(0, 0))


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
b = 1 / c(1.5, 2, 3.5, 5.5); n = 1;
curve.ExpS(n=n, b=b, k = 1, Vmax = 2, col=col, lwd=lwd, labels = lbls, xy.labels = xy1)
col = col.magenta; lbls = paste0("ESnb", 1:4);
n = 10/11; # same b;
curve.ExpS(n=n, b=b, k = 1, Vmax = 2, col=col, lwd=lwd, labels = lbls, xy.labels = xy2)
text(xy2[1] + 0.5, xy2[2] + 0.05, "n = 0.91", col = col.magenta[1], adj = c(0, 0))


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
b  = c(1, 2, 1/2, 1/3);
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
n  = c(1, 2, 1/2, 1/3); b = 1;
curve.Atan(b=b, n=n, Vmax = 2, col=col, xy.labels = xy, labels = "ATn")
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


###################

### Gompertz Models

gompertz.eq = expression(Vmax * exp(- k1 * exp(- k2*t^n)));
# Note: V(0) > 0;
gompertzAt.eq = expression(Vmax * exp(- k1 * exp(- k2*t^n)) * atan(k2 * t^m)*2/pi);

### Extension for negative t:
# V = Vmax * exp(- k1 * exp(- sign(t) * k2 * abs(t)^n));
### Extension: V(0) = 0
# V = Vmax * (exp(- k1 * exp(- k2*t^n)) - exp(-k1)) / (1 - exp(-k1));
# Note: the atan-extension seems more versatile;

# TODO


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

####################

# use.power = extra power;
fit.atan = function(data, init = NULL, use.power = FALSE) {
	if(use.power) {
		atan.frm = V ~ Vmax * (atan(b * t^n) * 2/pi)^k;
	} else {
		atan.frm = V ~ Vmax * (atan(b * t^n) * 2/pi);
	}
	if(is.null(init)) {
		init = list(Vmax = max(data$V), b=1, n=1);
		if(use.power) init$k = 1;
	}
	sol = nls(atan.frm, data, start=init);
	return(sol);
}
plot.pred = function(sol, data, add.zero = TRUE, ...) {
	tt = data$t; y = data$V;
	if(add.zero) { tt = c(0, tt); y = c(0, y); }
	yp = predict(sol, list(t = tt));
	matplot(tt, cbind(y, yp), type = "l", ...);
}

###
Vmax = 2
params = list(Vmax=Vmax, b = 1, n = 1)

### Uniform Time:
n  = 8; div = 2;
# n = 12; div = 3;
# n = 15; div = 3;
# n = 100; div = 10;
tt = seq(n) / div; params$t = tt;
y  = eval(MM.eq, params)
mm.df = data.frame(t = params$t, V = y)

###
sol = fit.atan(mm.df)
plot.pred(sol, mm.df)
print(sol)

# Note: very good fit if NO error-term;


### Uniform Time + Error:
n  = 8; div = 2;
# n = 12; div = 2;
# n = 12; div = 3;
# n = 15; div = 3;
# n = 25; div = 3;
# n = 25; div = 5;
tt = seq(n) / div; params$t = tt;
y  = eval(MM.eq, params);
mm.df = data.frame(t = params$t, V = y);
mm.df.err = mm.df;
# Error:
# Note: highly susceptible to random errors;
err = rnorm(n, 0, 0.1);
mm.df.err$V = mm.df$V + err;

###
sol = fit.atan(mm.df.err)
plot.pred(sol, mm.df)
points(tt, mm.df.err$V, col = "red")
print(sol)

# TODO: fit also the MM to the data;

