
########################
###
### Leonard Mada
### [the one and only]
###
### Exact Integration
### of Polynomial Fractions
### derived from Roots of Unity
###
### draft v.0.1b



### Pre-Requisits

### Roots of Unity:
# - Base-Integrals:
#   x^k / (x^n - 1) dx;
#   x^k / (x^n + 1) dx;
#   see file: Integrals.Fractions.Unity.R;
# - Powers of Roots of Unity:
#   x^k / (x^n - 1)^p dx, x^k / (x^n + 1)^p dx
#   see file: Integrals.Fractions.Unity.Powers.R;


###############
### History ###

# draft v.0.1b:
# - added more details to the examples with fractinal powers;
# draft v.0.1a:
# - moved examples from Integrals.Fractions.Unity.R
#   to this file;
# - added examples with fractional powers;


################
### Introduction


### A.) Fractional Powers
### B.) Trigonometric Derivations
### C.) Derived Polynomials

### Terminology

### I[k, n]
# I[k, n] = Integral x^k / (x^n + 1) dx;
# or
# I[k, n] = Integral x^k / (x^n - 1) dx;
# (depending on context)


#########################
### A.) Fractional Powers

### I[k/m, j/n]
# Integral x^(k/m) / (x^(j/n) - 1) dx
# =>
# y = x^(1/(m*n))
# x = y^(m*n)
# dx = m*n*y^(m*n - 1) * dy
# =>
# I[k/m, j/n] = I  m*n * y^(k*n + m*n - 1) / (x^(m*j) - 1) dy
# = m*n * I[k*n + m*n - 1, m*j];
# Note: limits of integration also change!


#################################
### B.) Trigonometric Derivations
# - TODO: document + expand;
# - see also file:
#   Integrals.Exercises.R;


###########################
### C.) Derived Polynomials
# - TODO


######################
######################

### helper function
I.f = function(f, lim) {
	integrate(f, lower=lim[1], upper=lim[2])$value
}
I.pf = function(b, n, lim) {
	if(length(b) > 1 && length(n) > 1 && length(b) != length(n)) {
		stop("Differing lengths!")
	}
	pow = n + 1
	coeff = b / pow
	sum(coeff * lim[2]^pow) - sum(coeff * lim[1]^pow)
}

#####################

################
### Examples ###
################

### A.) Fractional Powers

lim = c(1.2, 1.5)
#
integrate(function(x) x^(3/5)/(x^(5/4) + 1), lower=lim[1], upper=lim[2]) # ==
integrate(function(x) 20 * x^31/(x^25 + 1), lower=lim[1]^(1/20), upper=lim[2]^(1/20)) # ==
integrate(function(x) 20*x^6 - 20 * x^6/(x^25 + 1), lower=lim[1]^(1/20), upper=lim[2]^(1/20)) # ==
I.pf(20, 6, lim^(1/20)) - 20 * I.f(function(x) x^6/(x^25 + 1), lim^(1/20))
# I[6] = Integral x^6/(x^25 + 1) can be computed based on I[0] = 1/(x^25 + 1);
# I[0] => [1/x] => I[23] => [half] => I[11] => [half] => I[5]
#      => [half] => I2 => [1/x] => I[21] => [half] => I[10]
#      => [1/x] => I[13] => [half] => I[6];
# I[6] = f(I[0]) using the reverse sequence as per:
#   Integrals.Fractions.Unity.R;
# I[0] can be computed using the fraction decomposition as per:
#   Integrals.Fractions.Unity.R;


lim = c(1.2, 1.5)
#
integrate(function(x) (x^(3/5) + x^(1/3)) / (x^(5/4) + 1), lower=lim[1], upper=lim[2]) # ==
I.f(function(x) 20 * x^31/(x^25 + 1), lim=lim^(1/20)) +
	I.f(function(x) 12 * x^15/(x^15 + 1), lim=lim^(1/12)) # ==
I.f(function(x) 20*x^6 - 20 * x^6/(x^25 + 1), lim=lim^(1/20)) +
	I.f(function(x) 12 - 12/(x^15 + 1), lim=lim^(1/12)) # ==
I.pf(20, 6, lim=lim^(1/20)) - 20*I.f(function(x) x^6/(x^25 + 1), lim=lim^(1/20)) +
	I.pf(12, 0, lim=lim^(1/12)) - 12*I.f(function(x) 1/(x^15 + 1), lim=lim^(1/12))
# I[6, 25] & I[0, 15] can be computed exactly as per:
# Integrals.Fractions.Unity.R;


######################

### Other

### n = 5
# (and some n = 10 equivalent to n = 5)

lim = c(2, 3)
# trivially computable: see Integrals.Fractions.Unity.R;
integrate(function(x) 2 * x^7/(x^10 + 1), lower=lim[1], upper=lim[2])
integrate(function(x) x^3/(x^5 + 1), lower=lim[1]^2, upper=lim[2]^2)
### derived Integrals: *ALL* are equal!
# tan(y) = x => dx = 1/cos(y)^2 dy;
integrate(function(x) 2 * sin(x)^7*cos(x)/(sin(x)^10 + cos(x)^10), lower=atan(lim[1]), upper=atan(lim[2]))
# y = sin(x)^2 => dy = 2*sin(x)*cos(x)*dx;
integrate(function(x) x^3/(x^5 + (1-x)^5), lower=sin(atan(lim[1]))^2, upper=sin(atan(lim[2]))^2)
integrate(function(x) x^3/(x^5 + (1-x)^5), lower=1-1/(lim[1]^2+1), upper=1-1/(lim[2]^2+1))
# y = cos(x)^2 => dy = -2*sin(x)*cos(x)*dx
integrate(function(x) (x - 1)^3/(x^5 + (1-x)^5), lower=cos(atan(lim[1]))^2, upper=cos(atan(lim[2]))^2)
integrate(function(x) (x - 1)^3/(x^5 + (1-x)^5), lower=1/(lim[1]^2+1), upper=1/(lim[2]^2+1))
# from x^3 / Q(x): x = (y + 1)/2 => dx = 1/2 * dy
integrate(function(x) 2*(1+x)^3/((1+x)^5 + (1-x)^5), lower=1-2/(lim[1]^2+1), upper=1-2/(lim[2]^2+1))
integrate(function(x) (1+x)^3/(5*x^4 + 10*x^2 + 1), lower=1-2/(lim[1]^2+1), upper=1-2/(lim[2]^2+1))

# derived: (x^3 - (x-1)^3) / Q(x) = (3*x^2 - 3*x + 1) / Q(x)
# the 2nd Integral can be easily computed exactly;
I.f(function(x) (3*x^2 - 3*x + 1)/(x^5 + (1-x)^5), lim) # ==
I.f(function(x) x^3/(x^5 + 1), lim = c(1/(1-lim[1]) - 1, 1/(1-lim[2]) - 1)) -
	I.f(function(x) x^3/(x^5 + 1), lim = c(1/lim[1] - 1, 1/lim[2] - 1))



### n == 6 & 12
# paramater: can be changed;
n = 12
#
lim = c(2, 3)
#
n.half = n/2
# trivially computable: see Integrals.Fractions.Unity.R;
integrate(function(x) 1/2 * x^(n.half-2)/(x^n.half + 1), lower=lim[1]^2, upper=lim[2]^2)
integrate(function(x) x^(n-3)/(x^n + 1), lower=lim[1], upper=lim[2])
### derived Integrals: *ALL* are equal!
# tan(y) = x => dx = 1/cos(y)^2 dy;
integrate(function(x) sin(x)^(n-3)*cos(x)/(sin(x)^n + cos(x)^n), lower=atan(lim[1]), upper=atan(lim[2]))
# y = sin(x)^2 => dy = 2*sin(x)*cos(x)*dx;
integrate(function(x) 1/2 * x^(n.half - 2)/(x^n.half + (1-x)^n.half), lower=sin(atan(lim[1]))^2, upper=sin(atan(lim[2]))^2)
integrate(function(x) 1/2 * x^(n.half - 2)/(x^n.half + (1-x)^n.half), lower=1-1/(lim[1]^2+1), upper=1-1/(lim[2]^2+1))
# y = cos(x)^2 => dy = -2*sin(x)*cos(x)*dx
integrate(function(x) -1/2 * (1 - x)^(n.half - 2)/(x^n.half + (1-x)^n.half), lower=cos(atan(lim[1]))^2, upper=cos(atan(lim[2]))^2)
integrate(function(x) -1/2 * (1 - x)^(n.half - 2)/(x^n.half + (1-x)^n.half), lower=1/(lim[1]^2+1), upper=1/(lim[2]^2+1))
# from x^(n.half - 2) / Q(x): x = (y + 1)/2 => dx = 1/2 * dy
integrate(function(x) (1+x)^(n.half - 2)/((1+x)^n.half + (1-x)^n.half), lower=1-2/(lim[1]^2+1), upper=1-2/(lim[2]^2+1))



### n = 7
# (and some n = 14 equivalent to n = 7)
# see above for parametric version;

lim = c(2, 3)
# trivially computable: see Integrals.Fractions.Unity.R;
integrate(function(x) 2 * x^11/(x^14 + 1), lower=lim[1], upper=lim[2])
integrate(function(x) x^5/(x^7 + 1), lower=lim[1]^2, upper=lim[2]^2)
### derived Integrals: *ALL* are equal!
# tan(y) = x => dx = 1/cos(y)^2 dy;
integrate(function(x) 2 * sin(x)^11*cos(x)/(sin(x)^14 + cos(x)^14), lower=atan(lim[1]), upper=atan(lim[2]))
# y = sin(x)^2 => dy = 2*sin(x)*cos(x)*dx;
integrate(function(x) x^5/(x^7 + (1-x)^7), lower=sin(atan(lim[1]))^2, upper=sin(atan(lim[2]))^2)
integrate(function(x) x^5/(x^7 + (1-x)^7), lower=1-1/(lim[1]^2+1), upper=1-1/(lim[2]^2+1))
# y = cos(x)^2 => dy = -2*sin(x)*cos(x)*dx
integrate(function(x) (x - 1)^5/(x^7 + (1-x)^7), lower=cos(atan(lim[1]))^2, upper=cos(atan(lim[2]))^2)
integrate(function(x) (x - 1)^5/(x^7 + (1-x)^7), lower=1/(lim[1]^2+1), upper=1/(lim[2]^2+1))
# from x^5 / Q(x): x = (y + 1)/2 => dx = 1/2 * dy
integrate(function(x) 2*(1+x)^5/((1+x)^7 + (1-x)^7), lower=1-2/(lim[1]^2+1), upper=1-2/(lim[2]^2+1))
integrate(function(x) (1+x)^5/(7*x^6 + 35*x^4 + 21*x^2 + 1), lower=1-2/(lim[1]^2+1), upper=1-2/(lim[2]^2+1))

# derived: (x^5 - (x-1)^5) / Q(x) = (5*x^4 - 10*x^3 + 10*x^2 - 5*x + 1) / Q(x)
# the 2nd Integral can be easily computed exactly;
I.f(function(x) (5*x^4 - 10*x^3 + 10*x^2 - 5*x + 1)/(x^7 + (1-x)^7), lim) # ==
I.f(function(x) x^5/(x^7 + 1), lim = c(1/(1-lim[1]) - 1, 1/(1-lim[2]) - 1)) -
	I.f(function(x) x^5/(x^7 + 1), lim = c(1/lim[1] - 1, 1/lim[2] - 1))


###########

### [old]
# TODO: cleanup

n = 11/2
#
n.half = n/2
integrate(function(x) x^(n-3)/(x^n + 1), lower=lower, upper=upper)
integrate(function(x) sin(x)^(n-3)*cos(x)/(sin(x)^n + cos(x)^n), lower=atan(lower), upper=atan(upper))
integrate(function(x) 1/2 * x^(n.half - 2)/(x^n.half + (1-x)^n.half), lower=sin(atan(lower))^2, upper=sin(atan(upper))^2)

###########
###########

lower = 0.2 # 0.6
upper = 0.8
integrate(function(x) 1/2 * x^3/(x^5 + 1), lower=2/(1-lower)-1, upper=2/(1-upper)-1 )
integrate(function(x) (1+x)^3/((1+x)^5 + (1-x)^5), lower=lower, upper=upper)

### x^3 + x
log( ((1+upper)^5 + (1-upper)^5) / ((1+lower)^5 + (1-lower)^5) ) * 1/40
integrate(function(x) (x^3 + x)/((1+x)^5 + (1-x)^5), lower=lower, upper=upper)

### only x
(atan((upper^2+1)*sqrt(5)/2i) - atan((lower^2+1)*sqrt(5)/2i)) * sqrt(5)/2i / 20
integrate(function(x) x/((1+x)^5 + (1-x)^5), lower=lower, upper=upper)
# => x^3, 3*x^2 + 1

### 5*x^4 + 10*x^2 + 1
(upper - lower)/2
integrate(function(x) (5*x^4 + 10*x^2 + 1)/((1+x)^5 + (1-x)^5), lower=lower, upper=upper)


############
############

n = 5
lower = 0.2
upper = 0.8
integrate(function(x) 1/(x^(n/2) - 1), lower=lower, upper=upper)
integrate(function(x) 2 * x/(x^n - 1), lower=sqrt(lower), upper=sqrt(upper))
integrate(function(x) (x^(n/2) + 1)/(x^n - 1), lower=lower, upper=upper)


integrate(function(x) 1/(x^(n/3) - 1), lower=lower, upper=upper)
integrate(function(x) 3 * x^2/(x^n - 1), lower=lower^(1/3), upper=upper^(1/3))


integrate(function(x) x^(1/2)/(x^n - 1), lower=lower, upper=upper)
integrate(function(x) 2 * x^2/(x^(2*n) - 1), lower=sqrt(lower), upper=sqrt(upper))


integrate(function(x) x^(1/2)/(x^(n/3) - 1), lower=lower, upper=upper)
integrate(function(x) 6 * x^8/(x^(2*n) - 1), lower=lower^(1/6), upper=upper^(1/6))


