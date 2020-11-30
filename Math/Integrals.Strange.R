
### Leonard Mada
###
### Integrals: Strange
###
### draft v.0.1c


### TODO:
# - find & fix bugs;


### I( 1 / (x^sqrt(2) + 1) ) dx

# y = 1/2 * (x^sqrt(2) - x^(-sqrt(2)))
# dy = sqrt(2) * sqrt(y^2 + 1) / x * dx;
# I = I( x * 1/ (x * (x^sqrt(2) + 1)) ) dx
# = 1/sqrt(2) * x * I( 1 / (sqrt(y^2 + 1) * (y + 1 + sqrt(y^2 + 1))) ) dy -
#   1/sqrt(2) * I( I( 1 / (sqrt(y^2 + 1) * (y + 1 + sqrt(y^2 + 1))) )) dy;

# y = tan(z) =>
# I2 = I( 1 / (sin(z) + cos(z) + 1) ) dz

# conversion function
convert.powsq = function(x, b=-1, pow=sqrt(2)) {
	1/2 * (x^pow + b * x^(-pow))
}
# actual function
pow.fr = function(x, pow=sqrt(2)) {
	1 / (x^pow + 1)
}
pow.subst = function(x, pow=sqrt(2)) {
	x.sq = sqrt(x^2 + 1)
	div = pow * x.sq * (x + 1 + x.sq)
	1/div
}
integrate.pow = function(range, FUN=pow.fr, rel.tol=1E-10, ...) {
	integrate(FUN, lower=range[1], upper=range[2], rel.tol=rel.tol, ...)
}
integrate.xpow = function(range) {
	FUN = function(y) {
		x = (y + sqrt(y^2 + 1))^(1/sqrt(2))
		x * pow.subst(y)
	}
	integrate(FUN, lower=range[1], upper=range[2], subdivisions=1000)
}
integrate.2Dpow = function(xrg, yrg, FUN=pow.subst, subdivisions=1000, rel.tol=1E-7) {
	# TODO: ???
	I1 = xrg[2] * sincos.exact(atan(yrg[2])) - xrg[1] * sincos.exact(atan(yrg[1]))
	# TODO: replace double integral
	I2 = integrate(function(x)
		sapply(x, function(x) integrate.pow(c(yrg[1], x), FUN=pow.subst,
			subdivisions=subdivisions, rel.tol=rel.tol)$value),
		lower=yrg[1], upper=yrg[2])$value
	print(c(I1, I2))
	return(sqrt(2)*I1 - I2) # check if corect ???
}
### 1/(sin(x) + a * cos(x) + b)
sincos.exact = function(x, a=1, b=1) {
	c2 = 1/sqrt(a^2 + 1)
	# TODO: correct issues when a < 0!
	alfa = acos(c2) # atan2(a*c2, c2) # sign(a) * acos(c2)
	c1 = b * c2
	a.sq = sqrt(complex(re = 1 - c1^2, im=0))
	# print(a.sq)
	#
	# c2/2 * log((1 - cos(x + alfa))/(1 + cos(x + alfa)))
	# with fractions; alternative: atan;
	r1 = c2/2 * 1/a.sq * log((a.sq - cos(x + alfa))/(a.sq + cos(x + alfa)))
	# use result from 1/(sin(x)^2 - a^2);
	r2 = c2/2 * 1/a.sq * log((a.sq * tan(x + alfa) - c1)/(a.sq * tan(x + alfa) + c1))
	r1 - r2
	#
	r = c2/2 * 1/a.sq * log(
		(a.sq - cos(x + alfa)) / (a.sq + cos(x + alfa)) /
		(a.sq * tan(x + alfa) - c1)* (a.sq * tan(x + alfa) + c1))
}


### Test

rg = c(1, 2)
rg2 = convert.powsq(rg)
rg2

integrate.pow(rg)
# TODO: fix integral
integrate.2Dpow(rg, rg2)
integrate.xpow(rg2)


###########

### Test: Solution of I2

rg = c(0, pi/5)
rg2 = tan(rg)

### OK
integrate(function(z) 1 / (sin(z) + cos(z) + 1), lower=rg[1], upper=rg[2])
integrate.pow(rg2, FUN=pow.subst, pow=1)
diff(sincos.exact(rg))

