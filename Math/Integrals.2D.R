
### Leonard Mada
###
### Double Integrals


# including:
# - Gamma(1/n);
# - Int[0, Inf] (e^(-x^n));
# - double integrals;


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

e4.f = function(x) {
	exp(-x^4)
}
e4.2D.f = function(x, y) {
	exp(-x^4-y^4)
}
e4C.2D.f = function(c, d) {
	exp(-2*d) / 2 / sqrt(d^2 - c^4)
}
e4Inner.f = function(x) {
	Inner.Integral(x, e4.2D.f, inner.limits=c(0, Inf))
}
e4C.In.f = function(d) {
	Inner.L.Integral(d, e4C.2D.f, inner.limits.f=sqrt,
		subdivisions=512, rel.tol=1E-10, stop.on.error=F)
}

################

# direct Integral
g4 = int.f(e4.f)
g4 * 4 / gamma(1/4)

# a double integral
int.f(e4Inner.f)
(gamma(1/4) / 4)^2

# change of variables
int.f(e4C.In.f, subdivisions=512)
