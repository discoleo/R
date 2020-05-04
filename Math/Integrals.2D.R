
### Leonard Mada
###
### Double Integrals
### draft 0.3b

# Numerical Integration in R
# - with emphasis of double integrals;
#   [and mostly focused on Gamma-derivatives]
# including:
# - Gamma(1/n);
# - Int[0, Inf] (e^(-x^n));
# - double integrals:
#   generalization to various powers;
# - I I ... d phi d theta = Gamma(1/3)^3
#   and workout to a simple integral;

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

sphericalTan.f = function(x, t) {
	x / (x^3 * (sin(t)+cos(t)) * (1 - sin(2*t)/2) + 1)
}

InnerSp3.Integral = function(y, inner.f=sphericalTan.f, inner.limits=c(0, Inf), ...) {
	Inner.Integral(y, inner.f=inner.f, inner.limits=c(0, Inf), subdivisions=4096, ...)
}

# evaluate
x = int.f(InnerSp3.Integral, upper=pi/2, subdivisions=4096)

# Test
x^(1/3) * 3^(2/3)
gamma(1/3)


# Graph Test
# (but df has to be included as well)
t = pi/10
curve(spherical.f(t, x)*cos(x)^2, from=0, to=pi/2)
curve(sphericalTan.f(tan(x), t), from=0, to=pi/2, add=T, col="red", lty=3)


###################

sphericalTan2.f = function(x, t) {
	x / (x^3 + 1) / ((sin(t)+cos(t)) * (1 - sin(2*t)/2))^(2/3)
}

InnerSp32.Integral = function(y, inner.f=sphericalTan2.f, inner.limits=c(0, Inf), ...) {
	Inner.Integral(y, inner.f=inner.f, inner.limits=c(0, Inf), subdivisions=4096, ...)
}

# evaluate
x = int.f(InnerSp32.Integral, upper=pi/2, subdivisions=4096)

# Test
x^(1/3) * 3^(2/3)
gamma(1/3)


####################
####################

sphericalTrig.f = function(t) {
	1 / ((sin(t)+cos(t)) * (1 - sin(t)*cos(t)))^(2/3)
}

x = int.f(function(x) x /(x^3+1), lower=0, upper=Inf) * int.f(sphericalTrig.f, lower=0, upper=pi/2)
# Test
x^(1/3) * 3^(2/3)
gamma(1/3)

################

x = 2*pi / 3^(3/2) * int.f(sphericalTrig.f, lower=0, upper=pi/2)
# Test
x^(1/3) * 3^(2/3)
gamma(1/3)

################

x = int.f(sphericalTrig.f, lower=0, upper=pi/2)

# Test
x^(1/3) * 3^(1/6) * (2*pi)^(1/3)
gamma(1/3)

################

### x = sin(t) + cos(t)
sphericalFr.f = function(x) {
	1 / (x * (3 - x^2))^(2/3) / sqrt(2 - x^2)
}

x = int.f(sphericalFr.f, lower=1, upper=sqrt(2))

# Test
x^(1/3) * 3^(1/6) * 2^(8/9) * pi^(1/3)
gamma(1/3)


### alternative:
### x = (sin(t) + cos(t))^(1/3)
sphericalFr.f = function(x) {
	1 / ((3 - x^6)^(2/3) * sqrt(2 - x^6))
}

x = int.f(sphericalFr.f, lower=1, upper=2^(1/6))
# Test
x
gamma(1/3)^3 / 4 / 2^(2/3) / 3^(3/2) / pi

# Test: Gamma(1/3)
x^(1/3) * 2^(8/9) * 3^(1/2) * pi^(1/3)
gamma(1/3)


### Lower limit: 0 / different!
x = int.f(sphericalFr.f, lower=0, upper=2^(1/6))

# Test
x
gamma(1/3)^3 * 2^(1/3) / 3^(3/2) / pi * sin(pi/3)^2 * sin(pi/6)
gamma(1/3)^3 / 4 / 2^(2/3) / 3^(1/2) / pi


### Lower limit: 0, upper = 1 / different!
x = int.f(sphericalFr.f, lower=0, upper=1)

# Test
x
gamma(1/3)^3 / 4 / 2^(2/3) / 3^(1/2) / pi * 2/3

