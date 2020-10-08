
########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs
###
### draft v.0.1c-tr2


### History

### Order 1 Non-Liniar
### draft v.0.1c - v.0.1c-tr2:
# - added symmetrically shifted, eg:
#   y^2*dy - 2*y*dy - (x-1)*dy - y = x^2 - 2;
#   y^2*dy + (x^2 - x)*dy + (2*x-1)*y = 9*x^2 - 4*x; (v.0.1c-sh2)
# - added transformed base-polynomials:
#   y^3*dy + (x+1)*dy - 3*x*y^2 - y = 0; (v.0.1c-tr)
#   x^2*y*dy + (x+1)*dy - x*y^2 - 1/3*y = 0; (v.0.1c-tr2)
### draft v.0.1b-sh:
# - added classic/full Cardan Polynomials (P3);
# - added shifted version (P3) (v.0.1b-sh);
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
# q = (f - sqrt(f^2 - h^n))^(1/n)
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


###############
###############

### Transformed Polynomial

#########
### n = 3
# y^3 - 3*h*y - 2*f = 0
# y^2 - 3*h - 2*f/y = 0
# 2*y*dy - 3*dh - 2*df/y + 2*f*dy/y^2 = 0
# y^3*dy - 3/2*dh*y^2 - df*y + f*dy = 0

### Examples
# h(x) = x^2
# f(x) = x + 1
2*y^3*dy - 3*dh*y^2 - 2*df*y + 2*f*dy = 0
2*y^3*dy - 2*3*x*y^2 - 2*y + 2*(x+1)*dy = 0
y^3*dy + (x+1)*dy - 3*x*y^2 - y = 0
# y = (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, n=3) {
	det = sqrt((x+1)^2 - x^6 + 0i)
	r1 = round0(rootn(x+1 + det, n=3))
	r2 = round0(rootn(x+1 - det, n=3))
	return( round0(r1 + r2) )
}
dy = function(x) {
	y.x = y(x)
	div = (y.x^3 + x + 1)
	dp = (3*x*y.x^2 + y.x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-3, to=3)
sapply(c(-2, c(-(4:1), 1:4)/5, 2), line.tan, dx=3, p=y, dp=dy)

### Variant
# y^3*dy - 3/2*dh*y^2 - df*y + f*dy = 0
# (3*h*y + 2*f)*dy - 3/2*dh*y^2 - df*y + f*dy = 0
# 3*h*y*dy + 3*f*dy - 3/2*dh*y^2 - df*y = 0

### Examples
# h(x) = x^2
# f(x) = x + 1
3*h*y*dy + 3*f*dy - 3/2*dh*y^2 - df*y = 0
3*x^2*y*dy + 3*(x+1)*dy - 3*x*y^2 - y = 0
x^2*y*dy + (x+1)*dy - x*y^2 - 1/3*y = 0
# y = (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, n=3) {
	det = sqrt((x+1)^2 - x^6 + 0i)
	r1 = round0(rootn(x+1 + det, n=3))
	r2 = round0(rootn(x+1 - det, n=3))
	return( round0(r1 + r2) )
}
dy = function(x) {
	y.x = y(x)
	div = (x^2*y.x + x + 1)
	dp = (x*y.x^2 + 1/3 * y.x)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-3, to=3)
sapply(c(-2, c(-(4:1), 1:4)/5, 2), line.tan, dx=3, p=y, dp=dy)

###############
### Shifted ###

#########
### n = 3
### Shifted
# y^3 + 3*s*y^2 - 3*h*y - 2*f = 0
# 3*y^2*dy + 6*s*y*dy + 3*y^2*ds - 3*h*dy - 3*y*dh - 2*df = 0
# =>
# y^2*dy + 2*s*y*dy - h*dy + y^2*ds - y*dh - 2/3*df = 0
# (y^2 + 2*s*y - h)*dy + y^2*ds - y*dh - 2/3*df = 0

### Solution: reduce shift
# y => y - s =>
# y^3 - 3*(h + s^2)*y - 2*f + 2*s^3 + 3*h*s = 0

###
# h(x) = x
# f(x) = x^3
# s(x) = 1/2
y^2*dy + 2*s*y*dy - h*dy + y^2*ds - y*dh - 2/3*df = 0
y^2*dy + y*dy - x*dy - y - 2*x^2 = 0
# y = - s + (f + sqrt(f^2 - h^3))^(1/3) + (f - sqrt(f^2 - h^3))^(1/3)
y = function(x, n=3) {
	d = x^3 - 3/4*x - 1/8
	det = sqrt(d^2 - (x + 1/4)^3 + 0i)
	r1 = (d + det); r2 = (d - det)
	r = round0(rootn(r1, n=n) + rootn(r2, n=n))
	# shift back
	r = r - 1/2
	return(r)
}
dy = function(x) {
	y.x = y(x)
	div = (y.x^2 + y.x - x)
	dp = (y.x + 2*x^2)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-2, to=2)
sapply(c(-(4:1)/5, (1:6)/5), line.tan, dx=3, p=y, dp=dy)


#########
### n = 3
### y-Shifted
# y^3 + 3*s*y^2 - 3*h*y - 2*f = 0
# 3*y^2*dy + 6*s*y*dy + 3*y^2*ds - 3*h*dy - 3*y*dh - 2*df = 0
# =>
# y^2*dy + 2*s*y*dy - h*dy + y^2*ds - y*dh - 2/3*df = 0
# (y^2 + 2*s*y - h)*dy + y^2*ds - y*dh - 2/3*df = 0

### y-Shift: u = y*m
# (m^2*y^2 + 2*s*m*y - h)*(m*dy + y*dm) + m^2*y^2*ds - m*y*dh - 2/3*df = 0
# (m^3*y^2 + 2*s*m^2*y - h*m)*dy + m^2*y^3*dm + (m^2*ds + 2*m*s*dm)*y^2 - (m*dh + h*dm)*y - 2/3*df = 0
# or
### y-Shift: u = y + n
# ((y+n)^2 + 2*s*(y+n) - h)*(dy + dn) + (y + n)^2*ds - (y + n)*dh - 2/3*df = 0
# (y^2 + 2*(s+n)*y + n^2 + 2*s*n - h)*(dy + dn) + ds*y^2 - (dh - 2*n*ds)*y - 2/3*df + n^2*ds - n*dh = 0
# (y^2 + 2*(s+n)*y + n^2 + 2*s*n - h)*dy + (ds + dn)*y^2 - (dh - 2*n*ds - 2*(s+n))*y - 2/3*df + n^2*ds - n*dh +(n^2 + 2*s*n - h)*dn = 0

### Solution: reduce shift
# y => y - s =>
# y^3 - 3*(h + s^2)*y - 2*f + 2*s^3 + 3*h*s = 0

###
# h(x) = x - 1/4
# f(x) = x^3
# s(x) = 1/2
# n(x) = -1/2
(y^2 + 2*(s+n)*y + n^2 + 2*s*n - h)*dy + (ds + dn)*y^2 - (dh - 2*n*ds - 2*(s+n))*y - 2/3*df + n^2*ds - n*dh +(n^2 + 2*s*n - h)*dn = 0
y^2*dy - x*dy - y - 2*x^2 + 1/2 = 0
# simple: d(y^3 - 3*x*y) = 6*x^2 - 3/2;
# [but useful to test formulas]
# y = n - s + (f.sh + sqrt(f.sh^2 - h.sh^3))^(1/3) + (f.sh - sqrt(f.sh^2 - h.sh^3))^(1/3)
y = function(x, n=3) {
	d = x^3 - 3/4*x + 3/4*1/4 - 1/8
	det = sqrt(d^2 - (x - 1/4 + 1/4)^3 + 0i)
	r1 = (d + det); r2 = (d - det)
	r = round0(rootn(r1, n=n) + rootn(r2, n=n))
	# shift back
	r = r - 1/2 + 1/2
	return(r)
}
dy = function(x) {
	y.x = y(x)
	div = (y.x^2 - x)
	dp = (y.x + 2*x^2 - 1/2)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-2, to=2)
sapply(c(-(4:1)/5, (1:6)/5), line.tan, dx=3, p=y, dp=dy)


# df = from free term;
# y*dy, dy, y^2, y: 4 parameters needed;
# but currently only 3 available: h(x), s(x), n(x);


#############

#############
### n = 5 ###

# y^5 - 5*h*y^3 + 5*h^2*y - 2*f = 0
# 5*y^4*dy - 15*h*y^2*dy - 5*y^3*dh + 5*h^2*dy + 10*h*y*dh - 2*df = 0
# y^4*dy - 3*h*y^2*dy + h^2*dy - y^3*dh + 2*h*y*dh - 2/5 * df = 0

### TODO:
# - concrete examples;



######################
######################

### Shifted Symmetrically

### Base System:
# (p-s)^3 + (q-s)^3 = r
# p*q = c
# y = p + q

y^2*dy - 2*s*y*dy + (s^2 - c)*dy - y^2*ds + (2*s*ds - dc)*y = 1/3*dr - 2*c*ds - 2*s*dc + 2*s^2*ds

### Examples

### Example 1:
# s = 1
# c = x
# r = x^3
y^2*dy - 2*y*dy - (x-1)*dy - y = x^2 - 2
###
y = function(x, n=3) {
	s = 1
	d = 1/2 * (x^3 - 3*x + 1)
	det = sqrt(d^2 - x^3 + 0i)
	r1 = (d + det); r2 = (d - det)
	r = round0(rootn(r1, n=n) + rootn(r2, n=n))
	# shift back
	r = r + s
	return(r)
}
dy = function(x) {
	y.x = y(x)
	div = (y.x^2 - 2*y.x - x + 1)
	dp = (y.x + x^2 - 2)
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-2, to=2)
sapply(c(-(4:1)/5, (1:6)/5), line.tan, dx=3, p=y, dp=dy)

###########

### Shifted Symmetrically

### Base System:
# (p-s1)^3 + (q-s1)^3 + 3*b1*(p+q) = r
# (p-s2)*(q-s2) = c
# y = p + q

### Solution to Polynomial:
# (p-s2)*(q-s2) = c # =>
# p*q - s2*(p+q) + s2^2 - c # = 0
# p*q = s2*y + c - s2^2

p^3 + q^3 - 3*s1*(p^2+q^2) + 3*s1^2*(p+q) - 2*s1^3 + 3*b1*(p+q) = r
y^3 - 3*p*q*y - 3*s1*(y^2 - 2*p*q) + 3*(s1^2+b1)*y - 2*s1^3 - r = 0
y^3 - 3*(s2*y + c - s2^2)*y - 3*s1*(y^2 - 2*(s2*y + c - s2^2)) + 3*(s1^2+b1)*y - 2*s1^3 - r = 0
y^3 - 3*(s1+s2)*y^2 + 3*(s1^2 + s2^2 + 2*s1*s2 + b1 - c)*y + 6*s1*c - 2*s1^3 - 6*s1*s2^2 - r # = 0
# Shift: y => y + (s1+s2)
# Solution is based on this polynomial:
y^3 + (3*b1 - 3*c)*y + (3*b1*s1 + 3*b1*s2 + 3*c*s1 - 3*c*s2 - r - 3*s1*s2^2 + 3*s1^2*s2 - s1^3 + s2^3)

### ODE
y^2*dy - 2*(s1+s2)*y*dy - (ds1+ds2)*y^2 + (s1^2 + s2^2 + 2*s1*s2 + b1 - c)*dy +
 + (2*s1*ds1 + 2*s2*ds2 + 2*ds1*s2 + 2*s1*ds2 + db1 - dc)*y # =
# = - (2*ds1*c + 2*s1*dc - 2*s1^2*ds1 - 2*ds1*s2^2 - 4*s1*s2*ds2 - 1/3 * dr)

### Example 1:
# s1 = x
# s2 = -x
# b1 = x^2
# c = x
# r = x^3
y^2*dy + (x^2 - x)*dy + (2*x-1)*y = 9*x^2 - 4*x
###
y = function(x, n=3) {
	s1 = x; s2 = -x;
	b1 = x^2; c = x; r = x^3;
	d = -1/2 * (3*b1*s1 + 3*b1*s2 + 3*c*s1 - 3*c*s2 - r - 3*s1*s2^2 + 3*s1^2*s2 - s1^3 + s2^3)
	det = sqrt(d^2 - (c - b1)^3 + 0i)
	r1 = (d + det); r2 = (d - det)
	r.sol = round0(rootn(r1, n=n) + rootn(r2, n=n))
	# shift back
	r.sol = r.sol + s1 + s2
	return(r.sol)
}
dy = function(x) {
	y.x = y(x)
	div = (y.x^2 + x^2 - x)
	dp =  - y.x * (2*x-1) + 9*x^2 - 4*x
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-3, to=3)
# a nice local minimum
sapply(c(-(4:1)/5, (1:6)/5), line.tan, dx=3, p=y, dp=dy)


