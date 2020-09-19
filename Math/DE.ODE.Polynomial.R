
########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs
###
### draft v.0.1a-plot


### History

### draft v.0.1a-plot:
# - added diagnostic plots (+ tangents);
### draft v.01a:
# - initial draft:
#   some ODEs based on Cardan polynomials;


###################

###################
### Terminology ###

# let f(x), h(x) be 2 functions which are differentiable;
# p(x), q(x) = functions to be found;
# dp, dq:
# dp = d(p(x)); dq = d(q(x));
# Note: q(x) is an intermediary function used for various derivations;


####################

### helper functions
line.tan = function(x, col="red", dx=5) {
	slope = dp(x)
	x.max = ifelse( (abs(x) >= 1), dx*x, 10);
	isInf = abs(slope) == Inf
	x.max[isInf] = x[isInf]
	lines(c(x, x.max), c(p(x), p(x) + (x.max-x)*slope), col=col)
	return(slope)
}


##########################
### Cardan-Polynomials ###
##########################

### System:
# p^n + q^n = 2*f(x)
# p*q = h(x)
# => q = h / p
# =>
# n*p^(n-1)*dp + n*q^(n-1)*dq = 2 * df
# p^(n-1)*dp + q^(n-1)*dq = 2/n * df
# p^(n-1)*dp + (h/p)^(n-1) * (dh/p - h*dp/p^2) = 2/n * df
# p^(2*n)*dp + h^(n-1) * (p*dh - h*dp) - 2/n * p^(n+1) * df = 0
#
# p^(2*n)*dp - h^n*dp + h^(n-1)*p*dh - 2/n * p^(n+1) * df = 0

### Solutions
# p = (f + sqrt(f^2 - h^n))^(1/n)
# p = (f - sqrt(f^2 - h^n))^(1/n)
# Note:
# - 2 basic solutions are possible;
# - it is possible to rotate these solutions using the roots of unity;



############
### Examples

#########
### n = 2
p^4*dp - h^2*dp + h*p*dh - p^3 * df = 0

###
# h(x) = x
# f(x) = 1
p^4*dp - x^2*dp + x*p = 0
# p = sqrt(f + sqrt(f^2 - h^2))
p = function(x) {
	sqrt(1 + sqrt(1 - x^2))
}
curve(p, from=-1, to=1)


###
# h(x) = x
# f(x) = 2*x
p^4*dp - x^2*dp - 2*p^3 + x*p = 0
# p = sqrt(f + sqrt(f^2 - h^2))
p = function(x) {
	sqrt(2*x + x*sqrt(3))
}
curve(p, from=0, to=3)

#########

#########
### n = 3
p^6*dp - h^3*dp + h^2*p*dh - 2/3 * p^4 * df = 0

###
# h(x) = x
# f(x) = 1
p^6*dp - x^3*dp + x^2*p = 0
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (1 + sqrt(1 - x^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
curve(p, from=-3, to=1)

###
# h(x) = x
# f(x) = 3*x
p^6*dp - x^3*dp - 2*p^4 + x^2*p = 0
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (3*x + sqrt(9*x^2 - x^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	div = (p(x)^6 - x^3)
	dp = 2*p(x)^4 - x^2*p(x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=-3, to=9)
sapply(c(0.001, 2*(1:3), -1), line.tan, dx=3)

###
# h(x) = x
# f(x) = x^3 + 3*x
p^6*dp - x^3*dp - 2*(x^2+1)*p^4 + x^2*p = 0
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (x^3 + 3*x + sqrt((x^3 + 3*x)^2 - x^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	div = (p(x)^6 - x^3)
	dp = 2*(x^2+1)*p(x)^4 - x^2*p(x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=-3, to=9)
sapply(c(0.001, 2*(1:3), -1), line.tan, dx=2)


