####################
###
### Leonard Mada
###
### Integrals: Complex Analysis
### draft 0.1b


### I( (1 - cos(x^n)) / x^(n+1) )

# Ref:
# Complex Analysis: Integral of (1-cos(x))/x^2 using Contour Integration
# https://www.youtube.com/watch?v=PmQiZ5ipAAA


integrate.cosFr = function(n, iter=2000, tol=1E-8, print=TRUE) {
	r = integrate(function(x) (1-cos(x^n))/x^(n+1), lower=-Inf, upper=Inf, subdivisions=iter, abs.tol=tol);
	if(print) print(r$abs.error);
	return(r$value);
}
integrate.cosFrAbs = function(n, iter=2000, tol=1E-8, print=TRUE) {
	r = integrate(function(x) (1-cos(abs(x)^n))/abs(x)^(n+1), lower=-Inf, upper=Inf, subdivisions=iter, abs.tol=tol);
	if(print) print(r$abs.error);
	return(r$value);
}
###
# r = pi / n;  # for n = odd;

###
n = 3
integrate.cosFr(n) * n

###
n = 4
integrate.cosFrAbs(n) * n

###
n = 5
integrate.cosFr(n) * n

###
n = 1/2
integrate.cosFrAbs(n) * n

