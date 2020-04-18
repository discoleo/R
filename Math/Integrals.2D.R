
### Leonard Mada
###
### Double Integrals
### draft 0.1


# including:
# - Gamma(1/n);
# - Int[0, Inf] (e^(-x^n));
# - double integrals:
#   genralization to various powers;


### helper functions

Inner.Integral = function(y, inner.f, inner.limits=c(0, 1), ...) {
	sapply(y, 
    function(z) {
		integrate(
			function(x) inner.f(x, z), inner.limits[1], inner.limits[2], ...
			)$value
	})
}
Inner.L.Integral = function(y, inner.f, inner.limits.f, f.lower=0, ...) {
	sapply(y, 
    function(z) {
		integrate(
			function(x) inner.f(x, z), f.lower, inner.limits.f(z), ...
			)$value
	})
}

int.f = function(i.f, lower=0, upper=Inf, print=TRUE, ...) {
	rez.int = integrate(i.f, lower=lower, upper=upper, ...)
	if(print) {
		cat("Integral = ")
		print(rez.int)
	}
	return(rez.int$value)
}

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
Inner.L.gen = function(pow=4) {
	eCPow.2D.f = function(c, d) {
		exp(-2*d) / sqrt(d^2 - c^pow) * 2 / pow
	}
	eInner.f = function(d) {
		Inner.L.Integral(d, eCPow.2D.f, inner.limits.f=function(lim) lim^(2/pow),
			subdivisions=512, rel.tol=1E-10, stop.on.error=F)
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
g4 * n / gamma(1/n)

# a double integral
# int.f(e4Inner.f)
int.f(Inner.gen(n))
# check result
(gamma(1/n) / n)^2

### change of variables
# int.f(e4C.In.f, subdivisions=512)
int.f(Inner.L.gen(n), subdivisions=512)
