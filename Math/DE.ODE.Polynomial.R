
########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs
###
### draft v.0.1b


### History

### draft v.0.1b:
# - added classic/full Cardan Polynomials (P3);
### draft v.0.1a-plot - v.0.1a-px:
# - added diagnostic plots (+ tangent lines);
# - added more examples (v.0.1a-px);
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
line.tan = function(x, col="red", dx=5, p=p, dp=dp) {
	slope = dp(x)
	x.max = ifelse( (abs(x) >= 1), dx*x, 10);
	isInf = abs(slope) == Inf
	x.max[isInf] = x[isInf]
	lines(c(x, x.max), c(p(x), p(x) + (x.max-x)*slope), col=col)
	return(slope)
}
rootn = function(r, n) {
	ifelse( (Im(r) == 0 & Re(r) >= 0), r^(1/n), - (-r)^(1/n) )
}
### round()
round0 = function(m, tol=1E-7) {
	m[abs(Re(m)) < tol & abs(Im(m)) < tol] = 0
	isNotNA =  ! is.na(m)
	isZero = (Re(m) != 0) & (abs(Re(m)) < tol)
	if(sum(isZero[isNotNA]) > 0) {
		m[isZero] = complex(re=0, im=Im(m[isZero]))
	}
	isZero = (Im(m) != 0) & (abs(Im(m)) < tol)
	if(sum(isZero[isNotNA]) > 0) {
		m[isZero] = Re(m[isZero])
	}
	return(m)
}

##########################

########################
###     Order 1      ###
###   Non-Linear     ###
########################


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



################
### Examples ###

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
p^6*dp - h^3*dp - 2/3 * p^4 * df + h^2*p*dh = 0
# where df, dh = given;

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


###
# h(x) = x
# f(x) = 3/2 * log(x)
p^6*dp - h^3*dp - 2/3 * p^4 * df + h^2*p*dh = 0
p^6*dp - x^3*dp - 1/x*p^4 + x^2*p = 0
x*p^6*dp - x^4*dp - p^4 + x^3*p = 0
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (3/2 * log(x) + sqrt((3/2 * log(x))^2 - x^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	div = x*(p(x)^6 - x^3)
	dp = p(x)^4 - x^3*p(x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=0, to=0.8)
sapply(c(0.001, (1:3)/6), line.tan, dx=2)


###
# h(x) = 1/(x^2 + 1)
# f(x) = 3/2 * log(x^2 + 1)
p^6*dp - h^3*dp - 2/3 * p^4 * df + h^2*p*dh = 0
(x^2 + 1)^4 * p^6*dp - (x^2 + 1)*dp - 2*x*(x^2 + 1)^3 * p^4 - 2*x*p = 0
# p = (f + sqrt(f^2 - h^3))^(1/3)
p = function(x, n=3) {
	r = (3/2 * log(x^2+1) + sqrt(9/4 * log(x^2+1)^2 - 1/(x^2+1)^n))
	ifelse( (r >= 0), r^(1/n), - (-r)^(1/n) )
}
dp = function(x) {
	x.mult = (x^2+1)^3
	p.x = p(x)
	div = (x^2+1)*(x.mult*p.x^6 - 1)
	dp = 2*x*(x.mult*p.x^4 - p.x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(p, from=-5, to=5)
sapply(c(-(1:4), 1:4), line.tan, dx=3)


##########################
### Cardan-Polynomials ###
### Full Root          ###
##########################


### System:
# y = p + q;
# where:
# p^n + q^n = 2*f(x)
# p*q = h(x)

### Solutions
# y = p + q, where:
# p = (f + sqrt(f^2 - h^n))^(1/n)
# p = (f - sqrt(f^2 - h^n))^(1/n)
# Note:
# - it is possible to rotate these solutions using the roots of unity;

################
### Examples ###

#########
### n = 3
# y^3 - 3*h*y - 2*f = 0
# 3*y^2*dy - 3*h*dy - 3*y*dh - 2*df = 0
# y^2*dy - h*dy - y*dh - 2/3*df = 0

###
# h(x) = x
# f(x) = x^3
y^2*dy - h*dy - y*dh - 2/3*df = 0
y^2*dy - x*dy - y - 2*x^2 = 0
# y = (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, n=3) {
	r1 = (x^3 + sqrt(x^6 - x^n + 0i))
	r2 = (x^3 - sqrt(x^6 - x^n + 0i))
	r = round0(rootn(r1, n=n) + rootn(r2, n=n))
	return(r)
}
dy = function(x) {
	y.x = y(x)
	div = (y.x^2 - x)
	dp = (y.x + 2*x^2)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-2, to=2)
sapply(c(-(4:1)/5, (1:4)/5, 1.01), line.tan, dx=3, p=y, dp=dy)

###
# h(x) = x^2
# f(x) = 3*log(x)
y^2*dy - h*dy - y*dh - 2/3*df = 0
y^2*dy - x^2*dy - 2*x*y - 2/x = 0
x*y^2*dy - x^3*dy - 2*x^2*y - 2 = 0
# y = (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, n=3) {
	r1 = (3*log(x + 0i) + sqrt(9*log(x + 0i)^2 - x^(2*n) + 0i))
	### imaginary parts do NOT cancel for: x < 0;
	# sign.x = sign(x)
	# sign.x[x >= 0] = 1
	# r2 = (sign.x * 3*log(x + 0i) - sqrt(9*log(x + 0i)^2 - x^(2*n) + 0i))
	r1 = round0(rootn(r1, n=n))
	r2 = x^2 / r1
	return( round0(r1 + r2) )
}
dy = function(x) {
	y.x = y(x)
	div = (x*y.x^2 - x^3)
	dp = (2*x^2*y.x + 2)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=0.01, to=2, ylim=c(-5, 0.5)) # xlim=c(-2, 2)
sapply((1:4)/5, line.tan, dx=3, p=y, dp=dy)
# log(x) is complex for x < 0;
# curve(y, from=-2, to=-0.01, add=T, ylim=c(-5, 0.5))
# sapply(c(-(4:1)/5), line.tan, dx=3, p=y, dp=dy)


#########
### n = 5
# y^5 - 5*h*y^3 + 5*h^2*y - 2*f = 0
# 5*y^4*dy - 15*h*y^2*dy - 5*y^3*dh + 5*h^2*dy + 10*h*y*dh - 2*df = 0
# y^4*dy - 3*h*y^2*dy + h^2*dy - y^3*dh + 2*h*y*dh - 2/5 * df = 0

### TODO:
# - concrete examples;


