
### Leonard Mada
###
### Double Integrals
### draft 0.2

# Numerical Integration in R
# including:
# - Gamma(1/n);
# - Int[0, Inf] (e^(-x^n));
# - double integrals:
#   generalization to various powers;
# - I I ... d phi d theta = Gamma(1/3)^3

####################

### helper functions

# fixed inner Limits
Inner.Integral = function(y, inner.f, inner.limits=c(0, 1), ...) {
	sapply(y, 
    function(z) {
		integrate(
			function(x) inner.f(x, z), inner.limits[1], inner.limits[2], ...
			)$value
	})
}
# variable inner Limits
Inner.L.Integral = function(y, inner.f, inner.limits.f, f.lower=0, ...) {
	sapply(y, 
    function(z) {
		integrate(
			function(x) inner.f(x, z), f.lower, inner.limits.f(z), ...
			)$value
	})
}

### simple Integral
int.f = function(i.f, lower=0, upper=Inf, print=TRUE, ...) {
	rez.int = integrate(i.f, lower=lower, upper=upper, ...)
	if(print) {
		cat("Integral = ")
		print(rez.int)
	}
	return(rez.int$value)
}

### Concrete Functions

### Exponential functions

e.f = function(x, pow=4) {
	exp(-x^pow)
}
e.2D.f = function(x, y, pow=4) {
	exp( -x^pow - y^pow)
}
eC.2D.f = function(c, d, pow=4) {
	exp(-2*d) / 2 / sqrt(d^2 - c^pow)
}

# generate Inner Integrals
Inner.gen = function(pow=4) {
	ePow.2D.f = function(x, y) {
		exp( -x^pow - y^pow)
	}
	eInner.f = function(x) {
		Inner.Integral(x, ePow.2D.f, inner.limits=c(0, Inf))
	}
	return(eInner.f)
}
# I(n) = Int[0, Inf] Int[0, x^(2/n)] (e^(-2*x) / (sqrt(x^2 - y^n) * 2*n) ) dy dx
Inner.L.gen = function(pow=4, subdivisions=512) {
	eCPow.2D.f = function(c, d) {
		exp(-2*d) / sqrt(d^2 - c^pow) * 2 / pow
	}
	eInner.f = function(d) {
		Inner.L.Integral(d, eCPow.2D.f, inner.limits.f=function(lim) lim^(2/pow),
			subdivisions=subdivisions, rel.tol=1E-10, stop.on.error=F)
	}
	return(eInner.f)
}


### Cases: Pow = 4
e4Inner.f = Inner.gen(4)

# e4C.In.f = Inner.L.gen(4)
e4C.In.f = function(d) {
	Inner.L.Integral(d, e4C.2D.f, inner.limits.f=sqrt,
		subdivisions=512, rel.tol=1E-10, stop.on.error=F)
}

################

### n = 4
n = 4

# direct Integral
g4 = int.f(e.f, pow=n)
# check result
g4 * n / gamma(1/n)

# a double integral
# int.f(e4Inner.f)
int.f(Inner.gen(n))
# check result
(gamma(1/n) / n)^2

### change of variables
# int.f(e4C.In.f, subdivisions=512)
int.f(Inner.L.gen(n), subdivisions=512)




int.f(function(x) exp(-2*x)/x^(1-2/n))
# integrate e^(-2*x)/x^(1-2/5) dx from x=0 to x=Inf

int.f(function(x) 1/sqrt(1-x^2)/x^(1-2/n), upper=1)
# integrate 1/sqrt(1-x^2)/x^(1-2/5) dx from x=0 to x=1

gamma(2/5) / 2^(2/5) *
sqrt(pi) * gamma(1/5) / gamma(7/10) / 5^2 * 2
#
(gamma(1/5) / 5)^2


####################

### Test:
# Int[0, pi/2] Int[0, pi/2] sin(f)/(sin(f)^3 * (sin(t)+cos(t)) * (1 - sin(2*t)/2) + cos(f)^3) dt df
# Result:
# gamma(1/3)^3 / 3^2


spherical.f = function(t, f) {
	sin(f)/(sin(f)^3 * (sin(t)+cos(t)) * (1 - sin(2*t)/2) + cos(f)^3)
}

# Spherical
InnerSp3.Integral = function(y, inner.f=spherical.f, inner.limits=c(0, pi/2), ...) {
	Inner.Integral(y, inner.f=spherical.f, inner.limits=inner.limits, ...)
}

# evaluate
x = int.f(InnerSp.Integral, upper=pi/2)

x^(1/3) * 3^(2/3)
gamma(1/3)

#################


