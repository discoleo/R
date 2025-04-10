########################
###
### Leonard Mada
### [the one and only]
###
### Exact Integration of Trigonometric Radicals
### - Higher Powers
###   Integral( 1 / tan(x)^(1/p) ) dx
###   Integral( tan(x)^(1/p) ) dx
###
### draft v.0.2a


###############
### History ###


### draft v.0.1b - v.0.1f:
# - added: Integral( tan(x)^(1/p) ) dx;
# - examples with decomposition into fractions of unity:
#  -- decomposition for p == 5;
#  -- more decompositions (p = 7, p = 9, generic); [v.0.1e]
#  -- more decompositions: tan(x)^(1/p); [v.0.1f]
# - simplified some of the decompositions; [v.0.1d]
### draft v.0.1a:
# - initial draft:
#   Integral( 1 / tan(x)^(1/p) ) dx;


##############
### Theory ###

### Prerequisites:
# - see Fractions of Unity:
#   Integrals.Fractions.Unity.R;
# - see derived Integrals:
#   Integrals.Fractions.Unity.Derived.R;


########################
### 1 / tan(y)^(1/p) ###

### I(p) = I( 1 / tan(y)^(1/p) ) dy
### x^p = tan(y)^(p-1) =>
# I(p) = p/(p-1) * I( 1 / (x^(2*p/(p-1)) + 1) ) dx
### x = t^(p-1) => dx = (p-1)*t^(p-2)
# I(p) = p * I( t^(p-2) / (t^(2*p) + 1) ) dt
# = ...;
# - can be decomposed into a sum of Fractions of Unity;


####################
### tan(y)^(1/p) ###

### I(p) = I( tan(y)^(1/p) ) dy
### x^p = tan(y)^(p-1) =>
# I(p) = p/(p-1) * I( x^(2/(p-1)) / (x^(2*p/(p-1)) + 1) ) dx
### x = t^(p-1) => dx = (p-1)*t^(p-2)
# I(p) = p * I( t^p / (t^(2*p) + 1) ) dt
# = ...;
# - can be decomposed into a sum of Fractions of Unity;


### Examples

### n = 5
# I = I( 1 / tan(y)^(1/5) ) dy
# x^5 = tan(y)^4 =>
# I = 5/4 * I( 1 / (x^(5/2) + 1) ) dx # Rationalizing =>
#   = 5/4 * I( (x^(5/2) - 1) / (x^5 - 1) ) dx
#   = 5/4 * I( x^(5/2) / (x^5 - 1) ) dx - 5/4 * I( 1 / (x^5 - 1) ) dx
# x = t^2 => dx = 2*t*dt
# I = 5/2 * I( t^6 / (t^10 - 1) ) dt - 5/4 * I( 1 / (x^5 - 1) ) dx;
### Alternative: instead of Rationalizing
# x = t^2 => dx = 2*t*dt
# I = 5/2 * I( t / (t^5 + 1) ) dt;
# for simple fraction Decomposition, see:
# Integrals.Fractions.Unity.R;


### n = 7
# I = I( 1 / tan(y)^(1/7) ) dy
# x^7 = tan(y)^6 =>
# I = 7/6 * I( 1 / (x^(7/3) + 1) ) dx
# x = t^3 => dx = 3*t^2*dt
# I = 7/2 * I( t^2 / (t^7 + 1) ) dt;


### n = 9
# I = I( 1 / tan(y)^(1/9) ) dy
# x^9 = tan(y)^8 =>
# I = 9/8 * I( 1 / (x^(9/4) + 1) ) dx
# x = t^4 => dx = 4*t^3*dt
# I = 9/2 * I( t^3 / (t^9 + 1) ) dt;


### n = n # ODD
# I = I( 1 / tan(y)^(1/n) ) dy
# x^n = tan(y)^(n-1) =>
# I = n/(n-1) * I( 1 / (x^(2*n/(n-1)) + 1) ) dx
# x = t^((n-1)/2) => dx = (n-1)/2*t^((n-3)/2)*dt
# I = n/2 * I( t^((n-3)/2) / (t^n + 1) ) dt;


####################

### helper functions

# see DE.ODE.Helper.R
# TODO: move to specific Helper file;

### Other

tanp = function(x, p=5, inv=TRUE) {
	if(inv) 1 / rootn(tan(x), p)
	else rootn(tan(x), p)
}
sin.rp = function(x, p=5) {
	rootn(sin(x), p) / cos(x)
}
unity.rp = function(x, x.pow, n, p=2, b0=1) {
	if(missing(x.pow)) {
		1 / (rootn(x^n, p) + b0)
	} else {
		# TODO: multiple powers (but not critical)
		x.rn = rootn(x^n, p)
		x^x.pow / (x.rn + b0)
	}
}
unity.conj.rp = function(x, n, p=2, b0=1) {
	r = rootn(x^n, p);
	### Powers
	pow = 0:(p-1);
	inv.pow = p - 1 - pow;
	sign = rep(c(-1, 1), p %/% 2);
	isOdd = (p %% 2 == 1)
	if(isOdd) sign = c(1, sign);
	b0.pow = b0^inv.pow * sign;
	num = sapply(r, function(x) sum(x^pow * b0.pow) )
	num / (x^n + if(isOdd) b0^p else -b0^p )
}
convert.range = function(x, n, p=n-1, FUN=tan) {
	rootn(FUN(x)^p, n)
}

####################
####################

################
### Examples ###

lower = 1 + 1E-3
upper = 3; # 7/3

##########
### n == 5
p = 5;
rg = convert.range(c(lower, upper), n=p)
#
integrate(tanp, lower = lower, upper = upper, p=p, rel.tol=1E-10)
p/(p-1) * integrate(unity.rp, lower = rg[1], upper = rg[2], n=p, p=(p-1)/2)$value
2*p/(p-1) * integrate(unity.rp, lower = sqrt(rg[1]), upper = sqrt(rg[2]), x.pow=1, n=p, p=1)$value
# Rationalizing not needed!
p/(p-1) * integrate(unity.conj.rp, lower = rg[1], upper = rg[2], n=p, p=(p-1)/2)$value

# convert.range(c(lower, upper)): must NOT include 1;
# [because of the numerical integration (of the Rationalized fractions)]
lower = 1 + 1E-3; upper = 2; rg = convert.range(c(lower, upper), n=p);
# numerical integration:
integrate(tanp, lower = lower, upper = upper, p=p, rel.tol=1E-10)
5/2 * integrate(unity.rp, lower = sqrt(rg[1]), upper = sqrt(rg[2]), x.pow=1, n=5, p=1)$value # b0 == 1!
# [old] based on Rationalization:
# I = 5/2 * I( t^6 / (t^10 - 1) ) dt - 5/4 * I( 1 / (x^5 - 1) ) dx;
5/2 * integrate(unity.rp, lower = sqrt(rg[1]), upper = sqrt(rg[2]), x.pow=6, n=10, p=1, b0=-1)$value -
	5/4 * integrate(unity.rp, lower = rg[1], upper = rg[2], n=5, p=1, b0=-1)$value


##################
### (tan(x))^(1/n)
rg = convert.range(c(lower, upper), n=p)
#
integrate(tanp, lower = lower, upper = upper, p=p, inv=F, rel.tol=1E-10)
p/(p-1) * integrate(unity.rp, lower = rg[1], upper = rg[2], x.pow=1/2, n=p, p=(p-1)/2)$value
# p * I( t^5 / (t^10 + 1) ) dt;
p * integrate(unity.rp, lower = rootn(rg[1], p-1), upper = rootn(rg[2], p-1), x.pow=p, n=2*p, p=1)$value
# p/2 * I( u^2 / (u^5 + 1) ) du;
p/2 * integrate(unity.rp, lower = rootn(rg[1], (p-1)/2), upper = rootn(rg[2], (p-1)/2), x.pow=(p-1)/2, n=p, p=1)$value


##########
### n == 7
p = 7;
rg = convert.range(c(lower, upper), n=p)
#
integrate(tanp, lower = lower, upper = upper, p=p, rel.tol=1E-10)
p/(p-1) * integrate(unity.rp, lower = rg[1], upper = rg[2], n=p, p=(p-1)/2)$value
3*p/(p-1) * integrate(unity.rp, lower = rootn(rg[1], 3), upper = rootn(rg[2], 3), x.pow=3 - 1, n=p, p=1)$value
# Rationalizing not needed!
p/(p-1) * integrate(unity.conj.rp, lower = rg[1], upper = rg[2], n=p, p=(p-1)/2)$value

### (tan(x))^(1/n)
integrate(tanp, lower = lower, upper = upper, p=p, inv=F)
p/(p-1) * integrate(unity.rp, lower = rg[1], upper = rg[2], x.pow=1/3, n=p, p=(p-1)/2)$value


##########
### n == 9
p = 9;
rg = convert.range(c(lower, upper), n=p)
#
integrate(tanp, lower = lower, upper = upper, p=p, rel.tol=1E-10)
p/(p-1) * integrate(unity.rp, lower = rg[1], upper = rg[2], n=p, p=(p-1)/2)$value
4*p/(p-1) * integrate(unity.rp, lower = rootn(rg[1], 4), upper = rootn(rg[2], 4), x.pow=4 - 1, n=p, p=1)$value
# Rationalizing not needed!
p/(p-1) * integrate(unity.conj.rp, lower = rg[1], upper = rg[2], n=p, p=(p-1)/2)$value

### (tan(x))^(1/n)
integrate(tanp, lower = lower, upper = upper, p=p, inv=F)
p/(p-1) * integrate(unity.rp, lower = rg[1], upper = rg[2], x.pow=1/4, n=p, p=(p-1)/2)$value


##########
### n == 8
p = 8;
rg = convert.range(c(lower, upper), n=p)
#
integrate(tanp, lower = lower, upper = upper, p=p, rel.tol=1E-10)
p/(p-1) * integrate(unity.rp, lower = abs(rg[1]), upper = abs(rg[2]), n=2*p, p=(p-1))$value # TODO: check abs();
p * integrate(unity.rp, lower = abs(rootn(rg[1], p-1)), upper = abs(rootn(rg[2], p-1)), x.pow=p - 2, n=2*p, p=1)$value



###########################
###########################


### I(p) = I( sin(y)^(1/p) / cos(y) ) dy
### sin(y) = x^p => cos(y) dy = p*x^(p-1) dx =>
# I(p) = p * I( x^p / (1 - x^(2*p)) ) dx;

lower = 0; upper = pi/3;

p = 9;
rg = convert.range(c(lower, upper), n=p, p=1, FUN=sin)
#
integrate(sin.rp, lower = lower, upper = upper, p=p, rel.tol=1E-10)
-p * integrate(unity.rp, lower = rg[1], upper = rg[2], x.pow=p, n=2*p, p=1, b0=-1)$value

# [the classical way is simpler]
# I(p) = I( (tan(x))^(1/p) ) dx
### tan(x) = i*sin(y) =>
# (1 + tan(x)^2) dx = i*cos(y) dy =>
# dx = i/cos(y) dy =>
# I(p) =
# = 1i^(1 + 1/p) * I( sin(y)^(1/p) / cos(y) ) dy;


####################
####################

### on [0, pi/2]

###
integrate(\(x) sqrt(sin(x)), 0, pi/2)
sqrt(2/pi) * gamma(3/4)^2


### I( sin(x)^p * cos(x)^q )
p = 1/3; q = - 1/5;
integrate(\(x) sin(x)^p * cos(x)^q, 0, pi/2)
gamma((p+1)/2) * gamma((q+1)/2) / gamma((p+q)/2 + 1) / 2


####################
### on [0, pi/4] ###

### I( sin(x)^(5/3) / cos(x)^(1/3) )
integrate(\(x) sin(x)^(5/3) / cos(x)^(1/3), 0, pi/4)
(beta(1/3, 1/2) - beta(1/3, 1)) / 2^(2 + 2/3)

###
integrate(\(x) cos(x)^(5/3) / sin(x)^(1/3), 0, pi/4)
(beta(1/3, 1/2) + beta(1/3, 1)) / 2^(2 + 2/3)

### I( 1 / (cos(x) * sin(x)^(1/3)) )
integrate(\(x) 1 / (cos(x) * sin(x)^(1/3)), 0, pi/4)
- (digamma(1/3) + Euler)/2 + 3/4*log((2^(2/3) + 2^(1/3) + 1)/3) +
	- sqrt(3)/2 * atan((2^(1/3) + 1/2) * 2/sqrt(3)) +
	+ sqrt(3)/2 * atan((1 + 1/2) * 2/sqrt(3));


### I( sin(x)^(4/3) / cos(x)^(2/3) )
integrate(\(x) sin(x)^(4/3) / cos(x)^(2/3), 0, pi/4)
(beta(1/6, 1/2) - beta(1/6, 1)) / 2^(2 + 1/3)

###
integrate(\(x) cos(x)^(4/3) / sin(x)^(2/3), 0, pi/4)
(beta(1/6, 1/2) + beta(1/6, 1)) / 2^(2 + 1/3)


### I( sin(x)^(2/3) / cos(x)^(4/3) )
integrate(\(x) sin(x)^(2/3) / cos(x)^(4/3), 0, pi/4, rel.tol=1E-8)
(gamma(-1/6)*gamma(1/2)/gamma(1/3) - gamma(-1/6)*gamma(1)/gamma(5/6)) / 2^(2 - 1/3)
-6 * (beta(5/6, 1/2)/3 - beta(5/6, 1)*5/6) / 2^(2 - 1/3)
# Note: beta(5/6, 1)*5/6 == 1;
# (beta(-1/6, 1/2) - beta(-1/6, 1)) / 2^(2 - 1/3)

###
integrate(\(x) cos(x)^(2/3) / sin(x)^(4/3) - 1/x^(4/3), 0, pi/4, rel.tol=1E-8)
-6 * (beta(5/6, 1/2)/3 + beta(5/6, 1)*5/6) / 2^(2 - 1/3) + 3*(pi/4)^(-1/3)
# divergent: hypothetical result;
integrate(\(x) cos(x)^(2/3) / sin(x)^(4/3), 0, pi/4, rel.tol=1E-8)
-6 * (beta(5/6, 1/2)/3 + beta(5/6, 1)*5/6) / 2^(2 - 1/3);
# (beta(-1/6, 1/2) + beta(-1/6, 1)) / 2^(2 - 1/3)


### I( sin(x)^(1/3) / cos(x)^(5/3) )
integrate(\(x) sin(x)^(1/3) / cos(x)^(5/3), 0, pi/4, rel.tol=1E-8)
(gamma(-1/3)*gamma(1/2)/gamma(1/6) - gamma(-1/3)*gamma(1)/gamma(2/3)) / 2^(2 - 2/3)
- (beta(2/3, 1/2)/2 - 3) / 2^(2 - 2/3)
# (beta(-1/3, 1/2) - beta(-1/3, 1)) / 2^(2 - 2/3)

###
integrate(\(x) cos(x)^(1/3) / sin(x)^(5/3) - 1/x^(5/3), 0, pi/4, rel.tol=1E-8)
- (beta(2/3, 1/2)/2 + 3) / 2^(2 - 2/3) + 3/2*(pi/4)^(-2/3)
# divergent: hypothetical result;
integrate(\(x) cos(x)^(1/3) / sin(x)^(5/3), 0, pi/4, rel.tol=1E-8)
- (beta(2/3, 1/2)/2 + 3) / 2^(2 - 2/3)
# (beta(-1/3, 1/2) + beta(-1/3, 1)) / 2^(2 - 2/3)


# Derivation:
integrate(\(x) 2^(4/3) * sin(x)^(4/3) / cos(x)^(2/3), 0, pi/4)
integrate(\(x) cos(x)^(4/3) * (1 - sin(x))/cos(x)^2, 0, pi/2)
integrate(\(x) cos(x)^(-2/3) - sin(x)*cos(x)^(-2/3), 0, pi/2)
(beta(1/6, 1/2) - beta(1/6, 1))/2
# - alternative method: sum & diff;


###################
###################

### I( 1 / (1 + cos(x)^(1/n)) )
# Maths 505: One of my favorite integrals (so far)
# https://www.youtube.com/watch?v=FRfpBNhfHFc
# Note: Rationalization => 3 integrals;
# - can be extended to all integer fractions;
# - see also: Fractions.Rationalisation.R;


### I( 1 / (1 + cos(x)^(1/n)) )
integrate(\(x) 1 / (1 + sin(x)^(1/3)), 0, pi/2)
integrate(\(x) 1 / (1 + cos(x)^(1/3)), 0, pi/2)
1/2 * gamma(-1/2) * (gamma(5/6) / gamma(1/3) + gamma(7/6) / gamma(2/3) +
	- gamma(4/3) / gamma(5/6) - gamma(2/3) / gamma(1/6)) + 1;

### Helper:
integrate(\(x) cos(x)^(1/3) / (1 + cos(x)), 0, pi/2)
integrate(\(x) (cos(x)^(1/3) - cos(x)^(4/3)) / sin(x)^2, 0, pi/2)
# 1/2 * beta(2/3, -1/2) - 1/2 * beta(7/6, -1/2)
1/2 * gamma(-1/2) * (gamma(2/3) / gamma(1/6) - gamma(7/6) / gamma(2/3))

#
integrate(\(x) cos(x)^(2/3) / (1 + cos(x)), 0, pi/2)
integrate(\(x) (cos(x)^(2/3) - cos(x)^(5/3)) / sin(x)^2, 0, pi/2)
# 1/2 * beta(5/6, -1/2) - 1/2 * beta(4/3, -1/2)
1/2 * gamma(-1/2) * (gamma(5/6) / gamma(1/3) - gamma(4/3) / gamma(5/6))

#
integrate(\(x) 1 / (1 + cos(x)), 0, pi/2)
1;


### I( sin(x)^p / (1 - cos(x)^(1/3)) )
p = 4/3
integrate(\(x) sin(x)^p / (1 - cos(x)^(1/3)), 0, pi/2)
integrate(\(x) sin(x)^(p - 2) * (1 + cos(x)) *
	(cos(x)^(2/3) + cos(x)^(1/3) + 1), 0, pi/2);
1/2 * gamma((p-1)/2) * (
	gamma(5/6) / gamma(p/2 + 1/3) + gamma(4/3) / gamma(p/2 + 5/6) +
	+ gamma(2/3) / gamma(p/2 + 1/6) + gamma(7/6) / gamma(p/2 + 2/3) +
	+ gamma(1/2) / gamma(p/2) + 1 / gamma((p+1)/2)
)

# Note: same value;
integrate(\(x) 2^(4/3 - 1) * sin(x/2)^(4/3 - 2) * cos(x/2)^(4/3) *
	(cos(x)^(2/3) + cos(x)^(1/3) + 1), 0, pi/2)
# TODO: try to isolate components;


### Fraction Rationalization

### Ex 1:
p = 1/3
integrate(\(x) sin(x)^p / (cos(x)^(2/3) - 2*cos(x)^(1/3) + 2), 0, pi/2)
integrate(\(x) sin(x)^p *
	(cos(x)^(4/3) + 2*cos(x) + 2*cos(x)^(2/3) + 4*cos(x)^(1/3) + 4) /
	(cos(x)^2 + 4*cos(x) + 8), 0, pi/2)
# can be split now into 5 separate integrals;
# TODO


# Derivation:

# source("Polynomials.Helper.R")
pp = mult.pm(as.pm("x^2 - 2*x + 2"),
	as.pm("x^4 + b3*x^3 + b2*x^2 + b1*x + b0"));
tmp = lapply(as.coeff.pm(pp, "x"), print.pm)

# Test:
b = c(4, 4, 2, 2);
pp = mult.pm(as.pm("x^2 - 2*x + 2"),
	as.pm("x^4 + b[4]*x^3 + b[3]*x^2 + b[2]*x + b[1]"))
print.pm(pp, lead = "x")
# x^6 + 4*x^3 + 8


### Ex 2:
p = 1/3
integrate(\(x) sin(x)^p / (cos(x) - cos(x)^(2/3) - cos(x)^(1/3) + 2), 0, pi/2)
integrate(\(x) sin(x)^p *
	(cos(x)^2 + cos(x)^(5/3) + 2*cos(x)^(4/3) + 3*cos(x) +
		+ 3*cos(x)^(2/3) + 2*cos(x)^(1/3) + 4) /
	(cos(x)^3 + 2*cos(x)^2 + 5*cos(x) + 8), 0, pi/2)
# can be split now into 7 separate integrals;
# TODO


# Derivation:

# source("Polynomials.Helper.R")
pp = mult.pm(as.pm("x^3 - x^2 - x + 2"),
	as.pm("x^6 + b5*x^5 + b4*x^4 + b3*x^3 + b2*x^2 + b1*x + b0"));
tmp = lapply(as.coeff.pm(pp, "x"), print.pm)

# Test:
b = c(4,2,3,3,2,1)
pp = mult.pm(as.pm("x^3 - x^2 - x + 2"),
	as.pm("x^6 + b[6]*x^5 + b[5]*x^4 + b[4]*x^3 + b[3]*x^2 + b[2]*x + b[1]"))
print.pm(pp, lead = "x")
# x^9 + 2*x^6 + 5*x^3 + 8


####################
####################

### Other

### I( (1/tan(x)^2 - tan(x)^2)^(1/4) )
# Maths 505: An awesome calculus result I cooked up
# https://www.youtube.com/watch?v=WoxHJyOOma4
# - see also similar section in file:
#   Integrals.Fractions.Unity.Radicals.R;
# - that section is much more extensive;

### Pow = 1
integrate(\(x) sqrt(1/tan(x) - tan(x)) / 2, 0, pi/4)
integrate(\(x) sqrt(1 - x^4) / (x^4 + 1), 0, 1) # varia
pi/4

### Pow = 2
integrate(\(x) (1/tan(x)^2 - tan(x)^2)^(1/4), 0, pi/4)
integrate(\(x) 1/sqrt(2) * (cos(x) / sin(x)^2)^(1/4), 0, pi/2)
sqrt(2)/4 * beta(5/8, 1/4)

### Pow = 1/2
integrate(\(x) sqrt(1/tan(x)) - sqrt(tan(x)), 0, pi/4)
integrate(\(x) sqrt(1 - sin(x)) / sqrt(sin(x)) * sqrt(2)/2, 0, pi/2)
integrate(\(x) sqrt(1 - cos(x)) / sqrt(cos(x)) * sqrt(2)/2, 0, pi/2)
integrate(\(x) sin(x/2) / sqrt(cos(x)), 0, pi/2)
log(3 + 2*sqrt(2)) * sqrt(2)/2

