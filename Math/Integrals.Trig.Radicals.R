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
### draft v.0.1e


###############
### History ###


### draft v.0.1b - v.0.1e:
# - added: Integral( tan(x)^(1/p) ) dx;
# - examples with decomposition into fractions of unity:
#  -- decomposition for p == 5;
#  -- more decompositions (p = 7, p = 9, generic); [v.0.1e]
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


### I(p) = Integral( 1 / tan(y)^(1/p) dy )
#
# x^p = tan(y)^(p-1) =>
# I(p) = p/(p-1) * I( 1 / (x^(2*p/(p-1)) + 1) )
# = p/(p-1) * I( Rationalized(1 / (x^(2*p/(p-1)) + 1)) )
# = ...;
# - can be decomposed into a sum of Fractions of Unity;


### I(p) = Integral( tan(y)^(1/p) dy )
#
# x^p = tan(y)^(p-1) =>
# I(p) = p/(p-1) * I( x^(2/(p-1)) / (x^(2*p/(p-1)) + 1) )
# = p/(p-1) * I( Rationalized(x^(2/(p-1)) / (x^(2*p/(p-1)) + 1)) )
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
convert.range = function(x, n, p=n-1) {
	rootn(tan(x)^p, n)
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
lower = 1 + 1E-3; upper = 2; rg = convert.range(c(lower, upper), n=p);
# numerical integration:
integrate(tanp, lower = lower, upper = upper, p=p, rel.tol=1E-10)
# I = 5/2 * I( t^6 / (t^10 - 1) ) dt - 5/4 * I( 1 / (x^5 - 1) ) dx;
5/2 * integrate(unity.rp, lower = sqrt(rg[1]), upper = sqrt(rg[2]), x.pow=6, n=10, p=1, b0=-1)$value -
	5/4 * integrate(unity.rp, lower = rg[1], upper = rg[2], n=5, p=1, b0=-1)$value
5/2 * integrate(unity.rp, lower = sqrt(rg[1]), upper = sqrt(rg[2]), x.pow=1, n=5, p=1)$value # b0 == 1!


### (tan(x))^(1/n)
integrate(tanp, lower = lower, upper = upper, p=p, inv=F)
p/(p-1) * integrate(unity.rp, lower = rg[1], upper = rg[2], x.pow=1/2, n=p, p=(p-1)/2)$value


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



