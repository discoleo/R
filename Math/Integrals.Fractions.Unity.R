
########################
###
### Leonard Mada
### [the one and only]
###
### Exact Integration of Polynomial Fractions
### - Roots of unity:
###   Integral( 1 / (x^n - 1) )dx
### - Roots of minus unity:
###   Integral( 1 / (x^n + 1) )dx
### - TODO: Polynomial fractions:
###   Integral( P(x) / (x^n + 1) )dx
###
### version pre-1.alfa [draft]
### 2020-02-26
### - all powers: odd + even;
### - both roots of unity & minus unity;
### - TODO: fix "sporadic" bug for interval (-0.n, 0.n)
### 2020-02-24
### - all odd powers;
### 2020-02-23
### - all prime powers;
### [< 2020-02-20] the many calculations;


### helper functions ###

# Fraction Decomposition: 1/(x^n - 1)
decompose.fr = function(n, type=c("U+", "U-")) {
	### v 2.0
	### works both with n = odd & n = even!
	type = match.arg(type)
	
	### generate functions:
	# simple eval
	if(type == "U-") {
		diff = -1;
	} else {
		diff = 1;
	}
	eval.fr = function(x) {
		1 / (x^n - diff)
	}
	
	if(n %% 2 == 1) {
		# roots of unity
		m = roots1.conj(n)
		# coefficients
		b0 = 1/n
		a.sol = b0 * m$m.sum
		b.sol = rep(-2*b0, (n-1)/2)
		#
		m.shift = m$m.sum/2
		D = b.sol + a.sol*m.shift
		m.sq = sqrt(1 - m.shift^2)
		rez.f = function(x) {
			b0*log(x - 1) +
			sum(a.sol/2*log((x - m$m.conj[,1])*(x - m$m.conj[,2]))) +
			sum(D / m.sq * atan((x - m.shift)/m.sq))
		}
		if(type == "U-") {
			u.minus = -1; # complex(re=cos(pi/n), im=sin(pi/n))
			u.mult = 1; # - u.minus;
		} else {
			u.minus = 1;
			u.mult = 1;
		}
		# Exact integration
		integrate.exact = function(low, upper) {
			low = low / u.minus;
			upper = upper / u.minus;
			if(Re(low) < 1 | Re(upper) < 1) {
				low = as.complex(low)
				upper = as.complex(upper)
			}
			rez = rez.f(upper) - rez.f(low)
			return(rez * u.mult)
		}
	} else {
		if(type == "U-") {
			# TODO: bug with root of unity
			u.minus = complex(re=cos(pi/n), im=sin(pi/n))
			u.mult = - u.minus;
		} else {
			u.minus = 1;
			u.mult = 1;
		}
		f1 = decompose.fr(n/2)$integrate
		f2 = decompose.fr(n/2, type="U-")$integrate
		integrate.exact = function(low, upper) {
			low = low / u.minus;
			upper = upper / u.minus;
			(f1(low, upper) - f2(low, upper))/2 * u.mult;
		}
	}
	
	integrate.numeric = function(low, upper, subdivisions=4*1024, rel.tol=1E-10) {
		integrate(eval.fr, lower=low, upper=upper, subdivisions=subdivisions, rel.tol=rel.tol)
	}
	integrate.all = function(low, upper, subdivisions=4*1024, rel.tol=1E-10) {
		r.exact = integrate.exact(low, upper)
		r.num = try(
			integrate.numeric(low, upper, subdivisions, rel.tol)
		)
		r = list(exact=r.exact, num=r.num)
		return(r)
	}
	#
	return (list(
		n=n,
		eval=eval.fr,
		integrate = integrate.exact,
		integrate.num = integrate.numeric,
		integrate.all = integrate.all))
}


# generate Roots of Unity
roots1.conj = function(n, computeRotation=FALSE) {
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	n_1 = n - 1
	n.half = n_1/2
	i = 1:(n.half)
	m.all = m^(1:n_1)
	m.m = matrix(c(m.all[i], m.all[n-i]), ncol=2)
	m.sum = m.m[,1] + m.m[,2]
	#
	r = list(m=m, m.all=m^(1:n_1), m.conj=m.m, m.sum=m.sum)
	#
	if(computeRotation) {
		id.m = matrix( c(rep(i, n.half) * rep(i, rep(n.half, n.half))), ncol=n.half )
		id.m = id.m %% n
		id.m[id.m > n.half] = n - id.m[id.m > n.half]
		rot.m = matrix(m.sum[id.m], ncol=n.half)
		r$rot = rot.m
	}
	return(r)
}

### Separate functions
### Exact Integral
rez.f = function(x, a, b, b0, m.conj) {
	m.sum = m.conj[,1] + m.conj[,2]
	m.shift = m.sum/2
	D = b + a*m.shift
	m.sq = sqrt(1 - m.shift^2)
	#
	b0*log(x - 1) +
	sum(a/2*log((x - m.conj[,1])*(x - m.conj[,2]))) +
	sum(D / m.sq * atan((x - m.shift)/m.sq))
}
eval.fr = function(x, n) {
	1 / (x^n - 1)
}

######################
######################

# Note:
# - great care has been taken to exclude any errors;
# - if any errors still swept in:
#   it is easy to redo the calculations and correct the errors!

###############
### Case n = 5:
# 1 / (x^5 - 1)

f = decompose.fr(5)

lower = 1.001
upper = 1024
f$integrate.all(lower, upper)
# 1.137622+0i
# 1.137622 with absolute error < 2e-14


lower = 1 + 1E-8
upper = 1024
f$integrate.all(lower, upper)
# 3.439807+0i
# 3.439807 with absolute error < 1.4e-11


####################

###############
### Case n = 7:
# 1 / (x^7 - 1)

f = decompose.fr(7)

lower = 1.001
upper = 1024
f$integrate.all(lower, upper)
# 0.7468605+0i
# 0.7468605 with absolute error < 1.4e-14


lower = 1 + 1E-8
upper = 1024
f$integrate.all(lower, upper)
# 2.391136+0i
# 2.391136 with absolute error < 1.2e-11


####################

################
### Case n = 11:
# 1 / (x^11 - 1)

f = decompose.fr(11)

lower = 1.001
upper = 1024
f$integrate.all(lower, upper)
# 0.42502+0i
# 0.42502 with absolute error < 8e-15


lower = 1 + 1E-8
upper = 1024
f$integrate.all(lower, upper)
# 1.471195+0i
# 1.471195 with absolute error < 5.1e-12


################
### Case n = 13:
# 1 / (x^13 - 1)

# TODO


################
### Case n = 17:
# 1 / (x^17 - 1)

f = decompose.fr(17)

lower = 1.001
upper = 1024
f$integrate.all(lower, upper)
# 0.246099+0i
# 0.246099 with absolute error < 4.9e-15


lower = 1 + 1E-8
upper = 1024
f$integrate.all(lower, upper)
# 0.9228601+0i
# 0.9228601 with absolute error < 3.5e-12


################
### Case n = 23:
# 1 / (x^23 - 1)

f = decompose.fr(23)

lower = 1.001
upper = 1024
f$integrate.all(lower, upper)
# 0.167701+0i
# ERROR


lower = 1 + 1E-8
upper = 1024
f$integrate.all(lower, upper)
# 0.6677857+0i
# 0.6677857 with absolute error < 2.1e-12


#################
### Case n = 101:
# 1 / (x^101 - 1)

f = decompose.fr(101)

lower = 1.001
upper = 1024
f$integrate.all(lower, upper)
# 0.02335263+0i
# ERROR: wide discrepancy:
# 4.879488e-33 with absolute error < 9.7e-33


lower = 1 + 1E-8
upper = 1024
f$integrate.all(lower, upper)
# 0.1368511+0i
# ERROR: wide discrepancy:
# 5.118336e-33 with absolute error < 1e-32




########################
########################
######## OLD ###########
########################

### TODO:
# clean up;


### Case n = 5:
# 1 / (x^5 - 1)

### Explicit calculations
b0 - b1 - b2 # = 1
(b1*(m^2 + m^3) + b2*(m + m^4) -(a1 + a2) + b0 + b1 + b2)*x # = 0
(a1*(m^2 + m^3) + a2*(m + m^4) - b1*(m^2 + m^3) - b2*(m + m^4) + a1 + a2 + b0 - b1 - b2)*x^2 # = 0
(-a1*(m^2 + m^3) - a2*(m + m^4) - (a1 + a2) + b0 + b1 + b2)*x^3 # = 0
(a1 + a2 + b0)*x^4 # = 0
###
b0 = b1 + b2 + 1
a1 + a2 + b0 # = 0
- (a1 + a2) + 2*b0 + b1*(m^2 + m^3) + b2*(m + m^4) # = 1
a1*(m^2 + m^3) + a2*(m + m^4) - b1*(m^2 + m^3) - b2*(m + m^4) + a1 + a2 + 1 # = 0
# =>
a1*(m^2 + m^3) + a2*(m + m^4) + 2*b0 # = 0
-a1*(m^2 + m^3) - a2*(m + m^4) - (a1 + a2) + b0 + b1 + b2 # = 0
# =>
5*b0 - 1 # = 0
# =>
a1*(m^2 + m^3) + a2*(m + m^4) # = -2*b0 = -2/5
a1 + a2 # = -b0 = -1/5
# =>
b1*(m^2 + m^3) + b2*(m + m^4) # = 1 + (a1 + a2) - 2*b0 = 1 - 3*b0
b1 + b2 # = b0 - 1
###
b0 - b # = 1
bS2 + b0 + b - a # = 0
aS2 - bS2 + b0 - b + a # = 0
- aS2 + b0 + b - a # = 0
b0 + a # = 0


### Test
n = 5
# n_1 = n - 1
# m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
# m.all = m^(1:n_1)

n = 5
m.l = roots1.conj(n)
m = m.l$m
m.sum = m.l$m.sum
m.conj = m.l$m.conj
###

###
b0 = 1/5
#
b = b0 - 1
a = -b0
aS2 =  3*b0 - 1
bS2 = -3*b0 + 1
#
#
a.sol = b0 * m.sum
b.sol = rep(-2*b0, (n-1)/2)
#
# a1 = (aS2 - a*(m^1 + m^4)) / m.diff
# a2 = a - a1
# old
# a1 = b0 * (-2 + (m + m^4)  ) /  m.diff
# a1 = b0 / ((m^2 + m^3)^2 - 1)
# a1 = b0 / (m^1 + m^4 + 1)
# a2 = b0 * (-2 + (m^2 + m^3)) / -m.diff
# a2 = b0 / (m^2 + m^3 + 1)
#
# b1 = (1 - 3*b0 - (b0-1)*(m + m^4))   / m.diff
# b2 = (1 - 3*b0 - (b0-1)*(m^2 + m^3)) / -m.diff

# Test
x = 1.5
div = 1/(x^5 - 1)
#
b0/(x - 1) + sum( (a.sol*x + b.sol) / (x - m.conj[,1]) / (x - m.conj[,2]) )
div

##################
### Exact Integral
rez.f = function(x, a, b, b0, m.m) {
	m.sum = m.m[,1] + m.m[,2]
	m.shift = m.sum/2
	D = b + a*m.shift
	m.sq = sqrt(1 - m.shift^2)
	#
	b0*log(x - 1) +
	sum(a/2*log((x - m.m[,1])*(x - m.m[,2]))) +
	sum(D / m.sq * atan((x - m.shift)/m.sq))
}
fr.f = function(x, n) {
	1 / (x^n - 1)
}

# n_1 = n - 1
# i = 1:(n_1/2)
# m.all = m^(1:n_1)
# m.m = matrix(c(m.all[i], m.all[n-i]), ncol=2)

# Test Integral
low = 1.00 + 8E-11 # may need setting rel.tol larger!
top = 1024
# subdivisions = 4*256 # 4*1024
integrate(fr.f, lower=low, upper=top, subdivisions=4*1024, rel.tol=1E-10, n=n)
rez.f(top, a.sol, b.sol, b0, m.conj) - rez.f(low, a.sol, b.sol, b0, m.conj)


#####################
#####################

### n = 7
(b0 - b1 - b2 - b3) +
(- a1 - a2 - a3 + b0 - b1*m - b1*m^6 - b2*m^2 - b2*m^5 - b3*m^3 - b3*m^4)*x^1 +
(- a1*m - a1*m^6 - a2*m^2 - a2*m^5 - a3*m^3 - a3*m^4 + b0 - b1 - b1*m^2 - b1*m^5 - b2 - b2*m^3 - b2*m^4 - b3 - b3*m - b3*m^6)*x^2 +
(- a1 - a1*m^2 - a1*m^5 - a2 - a2*m^3 - a2*m^4 - a3 - a3*m - a3*m^6 + b0 + b1 + b1*m^2 + b1*m^5 + b2 + b2*m^3 + b2*m^4 + b3 + b3*m + b3*m^6)*x^3 +
(a1 + a1*m^2 + a1*m^5 + a2 + a2*m^3 + a2*m^4 + a3 + a3*m + a3*m^6 + b0 + b1*m + b1*m^6 + b2*m^2 + b2*m^5 + b3*m^3 + b3*m^4)*x^4 +
(a1*m + a1*m^6 + a2*m^2 + a2*m^5 + a3*m^3 + a3*m^4 + b0 + b1 + b2 + b3)*x^5 +
(a1 + a2 + a3 + b0)*x^6

#
a = a1 + a2 + a3
b = b1 + b2 + b3
aS1 = a1*m^1 + a1*m^6 + a2*m^2 + a2*m^5 + a3*m^3 + a3*m^4
aS2 = a1*m^2 + a1*m^5 + a2*m^3 + a2*m^4 + a3*m^1 + a3*m^6
bS1 = b1*m^1 + b1*m^6 + b2*m^2 + b2*m^5 + b3*m^3 + b3*m^4
bS2 = b1*m^2 + b1*m^5 + b2*m^3 + b2*m^4 + b3*m^1 + b3*m^6
# =>
        b0 - b # = 1
  - a + b0 - bS1 # = 0
- aS1 + b0 - bS2 - b # = 0
- aS2 + b0 + bS2 + b - a # = 0
  aS2 + b0 + bS1     + a # = 0
  aS1 + b0 + b # = 0
    a + b0 # = 0
# => sum =>

### Solution
n = 7
#
m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
#
b0 = 1/7
a = - b0
b = b0 - 1
aS1 = -b - b0
aS2 = -2*b0
bS1 = 2*b0
bS2 = 2*b0
###
a1 = b0 * (m^1+m^6)
a2 = -b0 / (m^1+m^6+1)
a3 = -b0 / (m^2+m^5+1)
b1 = (b0 - 1) * 2 / (n-1)
b2 = (b0 - 1) * 2 / (n-1)
b3 = (b0 - 1) * 2 / (n-1)
#
a = c(a1, a2, a3)
b = c(b1, b2, b3)
a
b

# [not needed]
# aS3 = -a - aS1 - aS2
# bS3 = -b - bS1 - bS2
# [not needed]
# aS13 = aS1 - (m^3 + m^4)*a
# aS23 = aS2 - (m^1 + m^6)*a
# (m^1 + m^6 - m^3 - m^4)*a1 + (m^2 + m^5 - m^3 - m^4)*a2
# (m^2 + m^5 - m^1 - m^6)*a1 + (m^3 + m^4 - m^1 - m^6)*a2
# a1 = ((m^3 + m^4 - m^1 - m^6)*aS13 - (m^2 + m^5 - m^3 - m^4)*aS23) / -7
# a2 = ((m^2 + m^5 - m^1 - m^6)*aS13 - (m^1 + m^6 - m^3 - m^4)*aS23) /  7
# a3 = a - a1 - a2
#
# bS13 = bS1 - (m^3 + m^4)*b
# bS23 = bS2 - (m^1 + m^6)*b
# b1 = ((m^3 + m^4 - m^1 - m^6)*bS13 - (m^2 + m^5 - m^3 - m^4)*bS23) / -7
# b2 = ((m^2 + m^5 - m^1 - m^6)*bS13 - (m^1 + m^6 - m^3 - m^4)*bS23) /  7
# b3 = b - b1 - b2

n = 7
m.l = roots1.conj(n)
m = m.l$m
m.sum = m.l$m.sum
m.conj = m.l$m.conj
#
b0 = 1/n
a.sol = b0 * m.sum
b.sol = rep(-2*b0, (n-1)/2)

### Test
x = 2
#
div = 1/(x^n - 1)
b0/(x - 1) + sum( (a.sol*x + b.sol) / (x - m.conj[,1]) / (x - m.conj[,2]))
div

##################
### Exact Integral

# Test Integral
low = 1.00 + 8E-11 # may need setting rel.tol larger!
top = 1024
# subdivisions = 4*256 # 4*1024
integrate(fr.f, lower=low, upper=top, subdivisions=4*1024, rel.tol=1E-10, n=n)
rez.f(top, a.sol, b.sol, b0, m.conj) - rez.f(low, a.sol, b.sol, b0, m.conj)
# curve(fr.f, from=low, to=top)


#####################

### n = 9
# TODO:
# m^3 + m^6 = -1
(- b + b0 - 1) +
(- a + b0 - bS1)*x^1 +
(- aS1 - b + b0 - bS2)*x^2 +
(- a - aS2 + b0 + b - 3*b3 - bS1)*x^3 +
(a - 3*a3 - aS1 + b0 - b + 3*b3 + bS1)*x^4 +
(- a + 3*a3 + aS1 + b + b0 + bS2)*x^5 +
(a + aS2 + b0 + bS1)*x^6 +
(aS1 + b + b0)*x^7 +
(a + b0)*x^8


(- b + b0 - 1) +
(- a + b0 - bS1)*x^1 +
(- aS1 - b + b0 - bS2)*x^2 +
(- a - aS2 + b0 + b1 + b2 - 2*b3 + b4 - bS1)*x^3 +
(a1 + a2 - 2*a3 + a4 - aS1 + b0 - b1 - b2 + 2*b3 - b4 + bS1)*x^4 +
(- a1 - a2 + 2*a3 - a4 + aS1 + b + b0 + bS2)*x^5 +
(a + aS2 + b0 + bS1)*x^6 +
(aS1 + b + b0)*x^7 +
(a + b0)*x^8


### TODO !!!

###
b0 - b = 1
b0 - bS1 - a = 0
b0 - aS1 - bS2 - b = 0
b0 - aS2 - bS1 + b - a = 0
b0 - aS1 + a = 0

### Solution
n = 9
#
m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
#
b0 = 1/n # TODO: ???
###
a = ...
b1 = (b0 - 1) * 2 / (n-1)
b2 = (b0 - 1) * 2 / (n-1)
b3 = (b0 - 1) * 2 / (n-1)
b4 = (b0 - 1) * 2 / (n-1)
#
a = c(a1, a2, a3, a4)
b = rep(b1, (n-1)/2) # c(b1, b2, b3, b4)
a
b

### Test
x = 2
n = 9
#
i = 1:((n-1)/2)
m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
#
div = 1/(x^n - 1)
b0/(x - 1) + sum( (a*x + b) / (x - m^i) / (x - m^(n - i)))
div

##################
### Exact Integral
rez.f = function(x, a, b, b0, m.m) {
	m.sum = m.m[,1] + m.m[,2]
	m.shift = m.sum/2
	D = b + a*m.shift
	m.sq = sqrt(1 - m.shift^2)
	#
	b0*log(x - 1) +
	sum(a/2*log((x - m.m[,1])*(x - m.m[,2]))) +
	sum(D / m.sq * atan((x - m.shift)/m.sq))
}
fr.f = function(x, n) {
	1 / (x^n - 1)
}

n_1 = n - 1
i = 1:(n_1/2)
m.all = m^(1:n_1)
m.m = matrix(c(m.all[i], m.all[n-i]), ncol=2)

# Test Integral
low = 1.00 + 8E-11 # may need setting rel.tol larger!
top = 1024
# subdivisions = 4*256 # 4*1024
integrate(fr.f, lower=low, upper=top, subdivisions=4*1024, rel.tol=1E-10, n=n)
rez.f(top, a, b, b0, m.m) - rez.f(low, a, b, b0, m.m)
# curve(fr.f, from=low, to=top)


#####################

### n = 11
(- b + b0) +
(- a + b0 - bS1)*x^1 +
(- aS1 - b + b0 - bS2)*x^2 +
(- a - aS2 + b + b0 + bS2 + bS4 + bS5)*x^3 +
(a + aS2 + aS4 + aS5 - b + b0 - bS2 - bS4)*x^4 +
(- a - aS2 - aS4 + b + b0 + bS2 + bS4)*x^5 +
(a + aS2 + aS4 - b + b0 - bS2 - bS4 - bS5)*x^6 +
(- a - aS2 - aS4 - aS5 + b + b0 + bS2)*x^7 +
(a + aS2 + b0 + bS1)*x^8 +
(aS1 + b + b0)*x^9 +
(a + b0)*x^10 

#
n  = 11
b0 = 1/n
b = b0 - 1
a = - b0
aS1 = - 2*b0 + 1
bS1 = 2*b0
bS2 = 2*b0
aS2 = - 2*b0
aS4 = - 2*b0
aS5 = - 2*b0
bS4 = 2*b0
bS5 = - 9*b0 + 1
#
aS3 = -2*b0
aS = c(aS1, aS2, aS3, aS4, aS5)
sum(aS) + a
# all b are the same!
b.sol = -2*b0   # (b0 - 1) * 2 / (n-1)



### Solution
n = 11
#
m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
n_1 = n - 1
n.half = n_1/2
i = 1:(n.half)
m.all = m^(1:n_1)
m.m = matrix(c(m.all[i], m.all[n-i]), ncol=2)
m.sum = m.m[,1] + m.m[,2]
#
id.m = matrix( c(rep(i, n.half) * rep(i, rep(n.half, n.half))), ncol=n.half )
id.m = id.m %% n
id.m[id.m > n.half] = n - id.m[id.m > n.half]
rot.m = matrix(m.sum[id.m], ncol=n.half)
solve(rot.m ,aS)

a.sol = m.sum * b0


### Test
x = 2
#
div = 1/(x^n - 1)
b0/(x - 1) + sum( (a.sol*x + b.sol) / (x - m^i) / (x - m^(n - i)))
div

##################
### Exact Integral
rez.f = function(x, a, b, b0, m.m) {
	m.sum = m.m[,1] + m.m[,2]
	m.shift = m.sum/2
	D = b + a*m.shift
	m.sq = sqrt(1 - m.shift^2)
	#
	b0*log(x - 1) +
	sum(a/2*log((x - m.m[,1])*(x - m.m[,2]))) +
	sum(D / m.sq * atan((x - m.shift)/m.sq))
}
fr.f = function(x, n) {
	1 / (x^n - 1)
}


# Test Integral
low = 1.00 + 8E-11 # may need setting rel.tol larger!
top = 1024
# subdivisions = 4*256 # 4*1024
integrate(fr.f, lower=low, upper=top, subdivisions=4*1024, rel.tol=1E-10, n=n)
rez.f(top, a.sol, b.sol, b0, m.m) - rez.f(low, a.sol, b.sol, b0, m.m)
# curve(fr.f, from=low, to=top)



#####################
#####################

### n = 13
n = 13
m.l = roots1.conj(n)
m = m.l$m
m.sum = m.l$m.sum
m.conj = m.l$m.conj
#
b0 = 1/n
a.sol = b0 * m.sum
b.sol = rep(-2*b0, (n-1)/2)

### Test
x = 2
#
div = 1/(x^n - 1)
b0/(x - 1) + sum( (a.sol*x + b.sol) / (x - m.conj[,1]) / (x - m.conj[,2]))
div

##################
### Exact Integral

# Test Integral
low = 1.00 + 8E-10 # may need setting rel.tol larger!
top = 1024
# subdivisions = 4*256 # 4*1024
integrate(fr.f, lower=low, upper=top, subdivisions=4*1024, rel.tol=1E-10, n=n)
rez.f(top, a.sol, b.sol, b0, m.conj) - rez.f(low, a.sol, b.sol, b0, m.conj)
# curve(fr.f, from=low, to=top)



#####################
#####################

### n = 17
n = 17
m.l = roots1.conj(n)
m = m.l$m
m.sum = m.l$m.sum
m.conj = m.l$m.conj
#
b0 = 1/n
a.sol = b0 * m.sum
b.sol = rep(-2*b0, (n-1)/2)

### Test
x = 2
#
div = 1/(x^n - 1)
b0/(x - 1) + sum( (a.sol*x + b.sol) / (x - m.conj[,1]) / (x - m.conj[,2]))
div

##################
### Exact Integral

# Test Integral
low = 1.00 + 8E-10 # may need setting rel.tol larger!
top = 1024
# subdivisions = 4*256 # 4*1024
integrate(fr.f, lower=low, upper=top, subdivisions=4*1024, rel.tol=1E-10, n=n)
rez.f(top, a.sol, b.sol, b0, m.conj) - rez.f(low, a.sol, b.sol, b0, m.conj)
# curve(fr.f, from=low, to=top)



#####################
#####################

### n = 19
n = 19
m.l = roots1.conj(n)
m = m.l$m
m.sum = m.l$m.sum
m.conj = m.l$m.conj
#
b0 = 1/n
a.sol = b0 * m.sum
b.sol = rep(-2*b0, (n-1)/2)

### Test
x = 2
#
div = 1/(x^n - 1)
b0/(x - 1) + sum( (a.sol*x + b.sol) / (x - m.conj[,1]) / (x - m.conj[,2]))
div

##################
### Exact Integral

# Test Integral
low = 1.00 + 8E-10 # may need setting rel.tol larger!
top = 1024
# subdivisions = 4*256 # 4*1024
integrate(fr.f, lower=low, upper=top, subdivisions=4*1024, rel.tol=1E-10, n=n)
rez.f(top, a.sol, b.sol, b0, m.conj) - rez.f(low, a.sol, b.sol, b0, m.conj)
# curve(fr.f, from=low, to=top)



##################
##################

### Case n = 15

(- b + b0 - 1) +
(- a + b0 - bS1)*x^1 +
(- aS1 - b + b0 - bS2)*x^2 +
(- a - aS2 + b0 + b1 + b1*m^6 + b1*m^9 + b2 + b2*m^3 + b2*m^12 + b3 + b3*m^3 + b3*m^12 + b4 + b4*m^6 + b4*m^9 - b5 + b5*m^3 + b5*m^6 + b5*m^9 + b5*m^12 + b6 + b6*m^6 + b6*m^9 + b7 + b7*m^3 + b7*m^12 - bS1)*x^3 +
(a1 + a1*m^6 + a1*m^9 + a2 + a2*m^3 + a2*m^12 + a3 + a3*m^3 + a3*m^12 + a4 + a4*m^6 + a4*m^9 - a5 + a5*m^3 + a5*m^6 + a5*m^9 + a5*m^12 + a6 + a6*m^6 + a6*m^9 + a7 + a7*m^3 + a7*m^12 - aS1 + b0 + b0*m + b0*m^2 + b0*m^4 + b0*m^5 + b0*m^7 + b0*m^8 + b0*m^10 + b0*m^11 + b0*m^13 + b0*m^14 - 2*b1 - b1*m^2 - b1*m^3 - b1*m^4 - b1*m^6 - b1*m^9 - b1*m^11 - b1*m^12 - b1*m^13 - 2*b2 - b2*m^3 - b2*m^4 - b2*m^6 - b2*m^7 - b2*m^8 - b2*m^9 - b2*m^11 - b2*m^12 - b3*m - b3*m^2 - b3*m^4 - b3*m^5 - b3*m^7 - b3*m^8 - b3*m^10 - b3*m^11 - b3*m^13 - b3*m^14 - 2*b4 - b4*m - b4*m^3 - b4*m^6 - b4*m^7 - b4*m^8 - b4*m^9 - b4*m^12 - b4*m^14 - 2*b5 - b5*m^3 - 2*b5*m^5 - b5*m^6 - b5*m^9 - 2*b5*m^10 - b5*m^12 - b6*m - b6*m^2 - b6*m^4 - b6*m^5 - b6*m^7 - b6*m^8 - b6*m^10 - b6*m^11 - b6*m^13 - b6*m^14 - 2*b7 - b7*m - b7*m^2 - b7*m^3 - b7*m^6 - b7*m^9 - b7*m^12 - b7*m^13 - b7*m^14)*x^4 +
(- 2*a1 - a1*m^2 - a1*m^3 - a1*m^4 - a1*m^6 - a1*m^9 - a1*m^11 - a1*m^12 - a1*m^13 - 2*a2 - a2*m^3 - a2*m^4 - a2*m^6 - a2*m^7 - a2*m^8 - a2*m^9 - a2*m^11 - a2*m^12 - a3*m - a3*m^2 - a3*m^4 - a3*m^5 - a3*m^7 - a3*m^8 - a3*m^10 - a3*m^11 - a3*m^13 - a3*m^14 - 2*a4 - a4*m - a4*m^3 - a4*m^6 - a4*m^7 - a4*m^8 - a4*m^9 - a4*m^12 - a4*m^14 - 2*a5 - a5*m^3 - 2*a5*m^5 - a5*m^6 - a5*m^9 - 2*a5*m^10 - a5*m^12 - a6*m - a6*m^2 - a6*m^4 - a6*m^5 - a6*m^7 - a6*m^8 - a6*m^10 - a6*m^11 - a6*m^13 - a6*m^14 - 2*a7 - a7*m - a7*m^2 - a7*m^3 - a7*m^6 - a7*m^9 - a7*m^12 - a7*m^13 - a7*m^14 - b0 - b0*m^3 - b0*m^5 - b0*m^6 - b0*m^9 - b0*m^10 - b0*m^12 + 2*b1 + b1*m^2 + b1*m^4 + b1*m^5 + b1*m^6 + b1*m^7 + b1*m^8 + b1*m^9 + b1*m^10 + b1*m^11 + b1*m^13 + 2*b2 + b2*m + b2*m^3 + b2*m^4 + b2*m^5 + b2*m^7 + b2*m^8 + b2*m^10 + b2*m^11 + b2*m^12 + b2*m^14 + 2*b3 + 2*b3*m^3 + b3*m^5 + 2*b3*m^6 + 2*b3*m^9 + b3*m^10 + 2*b3*m^12 + 2*b4 + b4*m + b4*m^2 + b4*m^5 + b4*m^6 + b4*m^7 + b4*m^8 + b4*m^9 + b4*m^10 + b4*m^13 + b4*m^14 - b5 - b5*m^5 - b5*m^10 + 2*b6 + 2*b6*m^3 + b6*m^5 + 2*b6*m^6 + 2*b6*m^9 + b6*m^10 + 2*b6*m^12 + 2*b7 + b7*m + b7*m^2 + b7*m^3 + b7*m^4 + b7*m^5 + b7*m^10 + b7*m^11 + b7*m^12 + b7*m^13 + b7*m^14)*x^5 +
(2*a1 + a1*m^2 + a1*m^4 + a1*m^5 + a1*m^6 + a1*m^7 + a1*m^8 + a1*m^9 + a1*m^10 + a1*m^11 + a1*m^13 + 2*a2 + a2*m + a2*m^3 + a2*m^4 + a2*m^5 + a2*m^7 + a2*m^8 + a2*m^10 + a2*m^11 + a2*m^12 + a2*m^14 + 2*a3 + 2*a3*m^3 + a3*m^5 + 2*a3*m^6 + 2*a3*m^9 + a3*m^10 + 2*a3*m^12 + 2*a4 + a4*m + a4*m^2 + a4*m^5 + a4*m^6 + a4*m^7 + a4*m^8 + a4*m^9 + a4*m^10 + a4*m^13 + a4*m^14 - a5 - a5*m^5 - a5*m^10 + 2*a6 + 2*a6*m^3 + a6*m^5 + 2*a6*m^6 + 2*a6*m^9 + a6*m^10 + 2*a6*m^12 + 2*a7 + a7*m + a7*m^2 + a7*m^3 + a7*m^4 + a7*m^5 + a7*m^10 + a7*m^11 + a7*m^12 + a7*m^13 + a7*m^14 + b0 + b0*m^3 - b0*m^5 + b0*m^6 + b0*m^9 - b0*m^10 + b0*m^12 - 2*b1 + b1*m - b1*m^2 - b1*m^3 - b1*m^6 - b1*m^9 - b1*m^12 - b1*m^13 + b1*m^14 - 2*b2 + b2*m^2 - b2*m^3 - b2*m^4 - b2*m^6 - b2*m^9 - b2*m^11 - b2*m^12 + b2*m^13 - b3*m - b3*m^4 - b3*m^5 - b3*m^10 - b3*m^11 - b3*m^14 - 2*b4 - b4*m^3 + b4*m^4 - b4*m^6 - b4*m^7 - b4*m^8 - b4*m^9 + b4*m^11 - b4*m^12 - 2*b5 - b5*m^3 - b5*m^6 - b5*m^9 - b5*m^12 - b6*m^2 - b6*m^5 - b6*m^7 - b6*m^8 - b6*m^10 - b6*m^13 - 2*b7 - b7*m - b7*m^3 - b7*m^6 + b7*m^7 + b7*m^8 - b7*m^9 - b7*m^12 - b7*m^14)*x^6 +
(- 2*a1 + a1*m - a1*m^2 - a1*m^3 - a1*m^6 - a1*m^9 - a1*m^12 - a1*m^13 + a1*m^14 - 2*a2 + a2*m^2 - a2*m^3 - a2*m^4 - a2*m^6 - a2*m^9 - a2*m^11 - a2*m^12 + a2*m^13 - a3*m - a3*m^4 - a3*m^5 - a3*m^10 - a3*m^11 - a3*m^14 - 2*a4 - a4*m^3 + a4*m^4 - a4*m^6 - a4*m^7 - a4*m^8 - a4*m^9 + a4*m^11 - a4*m^12 - 2*a5 - a5*m^3 - a5*m^6 - a5*m^9 - a5*m^12 - a6*m^2 - a6*m^5 - a6*m^7 - a6*m^8 - a6*m^10 - a6*m^13 - 2*a7 - a7*m - a7*m^3 - a7*m^6 + a7*m^7 + a7*m^8 - a7*m^9 - a7*m^12 - a7*m^14 - b0*m - b0*m^2 - b0*m^4 - 2*b0*m^5 - b0*m^7 - b0*m^8 - 2*b0*m^10 - b0*m^11 - b0*m^13 - b0*m^14 + 2*b1 - b1*m + b1*m^2 + b1*m^3 + b1*m^6 + b1*m^9 + b1*m^12 + b1*m^13 - b1*m^14 + 2*b2 - b2*m^2 + b2*m^3 + b2*m^4 + b2*m^6 + b2*m^9 + b2*m^11 + b2*m^12 - b2*m^13 + b3*m + b3*m^4 + b3*m^5 + b3*m^10 + b3*m^11 + b3*m^14 + 2*b4 + b4*m^3 - b4*m^4 + b4*m^6 + b4*m^7 + b4*m^8 + b4*m^9 - b4*m^11 + b4*m^12 + 2*b5 + b5*m^3 + b5*m^6 + b5*m^9 + b5*m^12 + b6*m^2 + b6*m^5 + b6*m^7 + b6*m^8 + b6*m^10 + b6*m^13 + 2*b7 + b7*m + b7*m^3 + b7*m^6 - b7*m^7 - b7*m^8 + b7*m^9 + b7*m^12 + b7*m^14)*x^7 +
(2*a1 - a1*m + a1*m^2 + a1*m^3 + a1*m^6 + a1*m^9 + a1*m^12 + a1*m^13 - a1*m^14 + 2*a2 - a2*m^2 + a2*m^3 + a2*m^4 + a2*m^6 + a2*m^9 + a2*m^11 + a2*m^12 - a2*m^13 + a3*m + a3*m^4 + a3*m^5 + a3*m^10 + a3*m^11 + a3*m^14 + 2*a4 + a4*m^3 - a4*m^4 + a4*m^6 + a4*m^7 + a4*m^8 + a4*m^9 - a4*m^11 + a4*m^12 + 2*a5 + a5*m^3 + a5*m^6 + a5*m^9 + a5*m^12 + a6*m^2 + a6*m^5 + a6*m^7 + a6*m^8 + a6*m^10 + a6*m^13 + 2*a7 + a7*m + a7*m^3 + a7*m^6 - a7*m^7 - a7*m^8 + a7*m^9 + a7*m^12 + a7*m^14 + b0 + b0*m^3 - b0*m^5 + b0*m^6 + b0*m^9 - b0*m^10 + b0*m^12 - 2*b1 - b1*m^2 - b1*m^4 - b1*m^5 - b1*m^6 - b1*m^7 - b1*m^8 - b1*m^9 - b1*m^10 - b1*m^11 - b1*m^13 - 2*b2 - b2*m - b2*m^3 - b2*m^4 - b2*m^5 - b2*m^7 - b2*m^8 - b2*m^10 - b2*m^11 - b2*m^12 - b2*m^14 - 2*b3 - 2*b3*m^3 - b3*m^5 - 2*b3*m^6 - 2*b3*m^9 - b3*m^10 - 2*b3*m^12 - 2*b4 - b4*m - b4*m^2 - b4*m^5 - b4*m^6 - b4*m^7 - b4*m^8 - b4*m^9 - b4*m^10 - b4*m^13 - b4*m^14 + b5 + b5*m^5 + b5*m^10 - 2*b6 - 2*b6*m^3 - b6*m^5 - 2*b6*m^6 - 2*b6*m^9 - b6*m^10 - 2*b6*m^12 - 2*b7 - b7*m - b7*m^2 - b7*m^3 - b7*m^4 - b7*m^5 - b7*m^10 - b7*m^11 - b7*m^12 - b7*m^13 - b7*m^14)*x^8 +
(- 2*a1 - a1*m^2 - a1*m^4 - a1*m^5 - a1*m^6 - a1*m^7 - a1*m^8 - a1*m^9 - a1*m^10 - a1*m^11 - a1*m^13 - 2*a2 - a2*m - a2*m^3 - a2*m^4 - a2*m^5 - a2*m^7 - a2*m^8 - a2*m^10 - a2*m^11 - a2*m^12 - a2*m^14 - 2*a3 - 2*a3*m^3 - a3*m^5 - 2*a3*m^6 - 2*a3*m^9 - a3*m^10 - 2*a3*m^12 - 2*a4 - a4*m - a4*m^2 - a4*m^5 - a4*m^6 - a4*m^7 - a4*m^8 - a4*m^9 - a4*m^10 - a4*m^13 - a4*m^14 + a5 + a5*m^5 + a5*m^10 - 2*a6 - 2*a6*m^3 - a6*m^5 - 2*a6*m^6 - 2*a6*m^9 - a6*m^10 - 2*a6*m^12 - 2*a7 - a7*m - a7*m^2 - a7*m^3 - a7*m^4 - a7*m^5 - a7*m^10 - a7*m^11 - a7*m^12 - a7*m^13 - a7*m^14 - b0 - b0*m^3 - b0*m^5 - b0*m^6 - b0*m^9 - b0*m^10 - b0*m^12 + 2*b1 + b1*m^2 + b1*m^3 + b1*m^4 + b1*m^6 + b1*m^9 + b1*m^11 + b1*m^12 + b1*m^13 + 2*b2 + b2*m^3 + b2*m^4 + b2*m^6 + b2*m^7 + b2*m^8 + b2*m^9 + b2*m^11 + b2*m^12 + b3*m + b3*m^2 + b3*m^4 + b3*m^5 + b3*m^7 + b3*m^8 + b3*m^10 + b3*m^11 + b3*m^13 + b3*m^14 + 2*b4 + b4*m + b4*m^3 + b4*m^6 + b4*m^7 + b4*m^8 + b4*m^9 + b4*m^12 + b4*m^14 + 2*b5 + b5*m^3 + 2*b5*m^5 + b5*m^6 + b5*m^9 + 2*b5*m^10 + b5*m^12 + b6*m + b6*m^2 + b6*m^4 + b6*m^5 + b6*m^7 + b6*m^8 + b6*m^10 + b6*m^11 + b6*m^13 + b6*m^14 + 2*b7 + b7*m + b7*m^2 + b7*m^3 + b7*m^6 + b7*m^9 + b7*m^12 + b7*m^13 + b7*m^14)*x^9 +
(2*a1 + a1*m^2 + a1*m^3 + a1*m^4 + a1*m^6 + a1*m^9 + a1*m^11 + a1*m^12 + a1*m^13 + 2*a2 + a2*m^3 + a2*m^4 + a2*m^6 + a2*m^7 + a2*m^8 + a2*m^9 + a2*m^11 + a2*m^12 + a3*m + a3*m^2 + a3*m^4 + a3*m^5 + a3*m^7 + a3*m^8 + a3*m^10 + a3*m^11 + a3*m^13 + a3*m^14 + 2*a4 + a4*m + a4*m^3 + a4*m^6 + a4*m^7 + a4*m^8 + a4*m^9 + a4*m^12 + a4*m^14 + 2*a5 + a5*m^3 + 2*a5*m^5 + a5*m^6 + a5*m^9 + 2*a5*m^10 + a5*m^12 + a6*m + a6*m^2 + a6*m^4 + a6*m^5 + a6*m^7 + a6*m^8 + a6*m^10 + a6*m^11 + a6*m^13 + a6*m^14 + 2*a7 + a7*m + a7*m^2 + a7*m^3 + a7*m^6 + a7*m^9 + a7*m^12 + a7*m^13 + a7*m^14 + b0 + b0*m + b0*m^2 + b0*m^4 + b0*m^5 + b0*m^7 + b0*m^8 + b0*m^10 + b0*m^11 + b0*m^13 + b0*m^14 - b1 + b1*m - b1*m^6 - b1*m^9 + b1*m^14 - b2 + b2*m^2 - b2*m^3 - b2*m^12 + b2*m^13 - b3 - b4 + b4*m^4 - b4*m^6 - b4*m^9 + b4*m^11 + b5 - b5*m^3 + b5*m^5 - b5*m^6 - b5*m^9 + b5*m^10 - b5*m^12 - b6 - b7 - b7*m^3 + b7*m^7 + b7*m^8 - b7*m^12)*x^10 +
(- a1 + a1*m - a1*m^6 - a1*m^9 + a1*m^14 - a2 + a2*m^2 - a2*m^3 - a2*m^12 + a2*m^13 - a3 - a4 + a4*m^4 - a4*m^6 - a4*m^9 + a4*m^11 + a5 - a5*m^3 + a5*m^5 - a5*m^6 - a5*m^9 + a5*m^10 - a5*m^12 - a6 - a7 - a7*m^3 + a7*m^7 + a7*m^8 - a7*m^12 + b - b0*m^3 - b0*m^6 - b0*m^9 - b0*m^12 + bS2)*x^11 +
(a + aS2 + b0 + bS1)*x^12 +
(aS1 + b + b0)*x^13 +
(a + b0)*x^14 

