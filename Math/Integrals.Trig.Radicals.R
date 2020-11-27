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
### draft v.0.1b


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
	# TODO: generalize
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


lower = 1 + 1E-3
upper = 3; # 7/3

###
p = 5;
rg = convert.range(c(lower, upper), n=p)
#
integrate(tanp, lower = lower, upper = upper, p=p)
p/(p-1) * integrate(unity.rp, lower = rg[1], upper = rg[2], n=p, p=(p-1)/2)$value
p/(p-1) * integrate(unity.conj.rp, lower = rg[1], upper = rg[2], n=p, p=(p-1)/2)$value

### (tan(x))^(1/n)
integrate(tanp, lower = lower, upper = upper, p=p, inv=F)
p/(p-1) * integrate(unity.rp, lower = rg[1], upper = rg[2], x.pow=1/2, n=p, p=(p-1)/2)$value


### n == 7
p = 7;
rg = convert.range(c(lower, upper), n=p)
#
integrate(tanp, lower = lower, upper = upper, p=p)
p/(p-1) * integrate(unity.rp, lower = rg[1], upper = rg[2], n=p, p=(p-1)/2)$value
p/(p-1) * integrate(unity.conj.rp, lower = rg[1], upper = rg[2], n=p, p=(p-1)/2)$value

### (tan(x))^(1/n)
integrate(tanp, lower = lower, upper = upper, p=p, inv=F)
p/(p-1) * integrate(unity.rp, lower = rg[1], upper = rg[2], x.pow=1/3, n=p, p=(p-1)/2)$value


### n == 9
p = 9;
rg = convert.range(c(lower, upper), n=p)
#
integrate(tanp, lower = lower, upper = upper, p=p)
p/(p-1) * integrate(unity.rp, lower = rg[1], upper = rg[2], n=p, p=(p-1)/2)$value
p/(p-1) * integrate(unity.conj.rp, lower = rg[1], upper = rg[2], n=p, p=(p-1)/2)$value

### (tan(x))^(1/n)
integrate(tanp, lower = lower, upper = upper, p=p, inv=F)
p/(p-1) * integrate(unity.rp, lower = rg[1], upper = rg[2], x.pow=1/4, n=p, p=(p-1)/2)$value



