
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
### - Polynomial fractions:
###   Integral( P(x) / (x^n - 1) )dx
### - Polynomial fractions:
###   Integral( P(x) / (x^n - 1)^2 )dx
###
### version 1 [RC1] [draft]
###
### 2020-03-01
### - polynomial fractions: P(x) / (x^n - 1)^2
###   Cases: n=3, n=5, n=7;
### 2020-02-29
### 2020-02-28
### - Cases: n = 9, n = 15;
### 2020-02-27
### - polynomial fractions: P(x) / (x^n - 1)
###  -- n = all primes (straightforward);
###     TODO: generating function for any n (prime);
###  -- TODO: n = 9; [SOLVED]
### - removed old derivations;
### 2020-02-26
### - all powers: odd + even;
### - both roots of unity & minus unity;
### - TODO:
###   fix "sporadic" bug for interval (-0.n, 0.n)
### 2020-02-24
### - all odd powers;
### 2020-02-23
### - all prime powers;
### [< 2020-02-20] the many calculations;


### Part A:
### Base fraction: 1 / (x^n - 1)
### Part B:
### Polynomial fraction: P(x) / (x^n - 1)
### Part C:
### Polynomial fraction: P(x) / (x^n - 1)^2


### Publishing:
### - it is easy nowadays to find Plagiarism;
### - it is easy to involve me in specific projects,
###   including collaborations on articles/books/other projects;

########################

### Part A

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

# TODO: add results


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


#########################
#########################
#########################

### Part B

#########################
###
### Polynomial Fractions
###
#########################

# generic functions
I.num.gen = function(n, base.pow=1) {
	f = function(low, upper, pow) {
		integrate(function(x) x^pow/(x^n - 1)^base.pow, low, upper)
	}
	return(f)
}
I.n1.gen = function(n) {
	f = function(low, upper) {
		r = log((upper^n - 1)/(low^n - 1))
		r = r / n
		return(r)
	}
	return(f)
}
I.gen = function(n) {
	int.f = decompose.fr(n);
	I0 = int.f$integrate
	
	In_1 = I.n1.gen(n) # n - 1
	In_2 = Iinv.gen(I0) # n - 2
	Ihalf = Ihalf.gen(I0)
	Ihalf_inv = Iinv.gen(Ihalf)
	
	r = c(I0, In_1, In_2, Ihalf, Ihalf_inv)
	return(r)
}
Iinv.gen = function(int.f, sign.inverse=FALSE, limits.inverse=FALSE) {
	mult = ifelse(sign.inverse, -1, 1);
	f = function(low, upper) {
		if(limits.inverse) {
			low = -low; upper = -upper; # WHY is needed?
		}
		int.f(1/low, 1/upper) * mult
	}
	return(f)
}
Ihalf.gen = function(int.f, sign.inverse=FALSE, asSum=TRUE, div=1) {
	# div is used for: 4*x^n = 1/(x^n-1)^2 - 1/(x^n+1)^2
	Ihalf = function(low, upper) {
		if(sign.inverse) {
			low = -low;
			upper = -upper;
		}
		# valid only for odd base powers & even exponents
		if(Re(low) < 0 | Re(upper) < 0) {
			low = as.complex(low)
			upper = as.complex(upper)
		}
		low = sqrt(low)
		upper = sqrt(upper)
		if(asSum) {
			r = int.f(low, upper) + int.f(-low, -upper)
		} else {
			r = int.f(low, upper) - int.f(-low, -upper)
		}
		return(r / div)
	}
}
Idx.gen = function(n, dx.pow, div=dx.pow) {
	I = decompose.fr(n);
	I0 = I$integrate
	Id = function(low, upper) {
		low = low^dx.pow
		upper = upper^dx.pow
		r = I0(low, upper) / div
		return(r)
	}
	return(Id)
}
# Composite Fractions
Idxfr.gen = function(n) {
	I.base = I.gen(n);
	I0.base = I.base[[1]]
	# I0.base = (x^n - 1) / (x^n - 1)^2
	#
	# d (1 / (x^n - 1)) = -n*x^(n-1) / (x^n - 1)^2
	# d (x / (x^n - 1)) = -((n-1)*x^n + 1) / (x^n - 1)^2
	# d (x^j / (x^n - 1)) = -((n-j)*x^(n + j - 1) + j*x^(j - 1)) / (x^n - 1)^2
	Ipow = function(x, pow, base.pow=1) {
		if(pow == 0) {
			r = 1 / (x^n - 1)^base.pow
		} else {
			r = x^pow / (x^n - 1)^base.pow
		}
		return(r)
	}
	#
	I0 = function(low, upper) {
		r = -((n-1)*I0.base(low, upper) + Ipow(upper, 1) - Ipow(low, 1)) / n;
		return(r)
	}
	In = function(low, upper) {
		r = I0.base(low, upper) + I0(low, upper)
		return(r)
	}
	In_1 = function(low, upper) {
		r = -(Ipow(upper, 0) - Ipow(low, 0)) / n
		return(r)
	}
	#  d ((x^n - 1)^2) = 2*n * x^(n-1)*(x^n - 1)
	I2n_1 = function(low, upper) {
		r = log( (upper^n - 1) / (low^n - 1) ) / n +
			In_1(low, upper)
		return(r)
	}
	# d (x^j / (x^n - 1)) = -((n-j)*x^(n + j - 1) + j*x^(j - 1)) / (x^n - 1)^2
	Ifrhigh.gen = function(pow, Ipow_1) {
		# I.pow = pow - 1
		Ijpow = function(low, upper) {
			r = - Ipow(upper, pow) + Ipow(low, pow) - (n - pow) * Ipow_1(low, upper)
			r = r / pow
			return(r) # r.pow = j - 1
		}
		return(Ijpow)
	}
	Ifrlow.gen = function(pow, Ipow_1) {
		# I.pow = pow - 1
		Ijpow = function(low, upper) {
			r = - Ipow(upper, pow) + Ipow(low, pow) - pow * Ipow_1(low, upper)
			r = r / (n - pow)
			return(r) # r.pow = n + j - 1
		}
		return(Ijpow)
	}
	r.f = c(I0, I2n_1, In, In_1, Ifrhigh.gen, Ifrlow.gen);
	#
	return(r.f)
}

##############

### Test

### Case n = 5
n = 5

I.num = I.num.gen(n)
# EXACT INTEGRALS
I = I.gen(n)
#
I0 = I[[1]]
I4 = I[[2]]
I3 = I[[3]]
I2 = I[[4]]
I1 = I[[5]]

### Test
low = 2
upper = 4

I4(low, upper)
I.num(low, upper, 4)

I3(low, upper)
I.num(low, upper, 3)

I2(low, upper)
I.num(low, upper, 2)

I1(low, upper)
I.num(low, upper, 1)

I0(low, upper)
I.num(low, upper, 0)

#####################

### Case n = 7
n = 7

I.num = I.num.gen(n)
# EXACT INTEGRALS
I = I.gen(n)
#
I0 = I[[1]]
I6 = I[[2]]
I5 = I[[3]]
I3 = I[[4]]
I2 = I[[5]]
I4 = Ihalf.gen(I2)
I1 = Iinv.gen(I4)

### Test
low = 2
upper = 4

I6(low, upper)
I.num(low, upper, 6)

I5(low, upper)
I.num(low, upper, 5)

I4(low, upper)
I.num(low, upper, 4)

I3(low, upper)
I.num(low, upper, 3)

I2(low, upper)
I.num(low, upper, 2)

I1(low, upper)
I.num(low, upper, 1)

I0(low, upper)
I.num(low, upper, 0)


#####################

### Case n = 9
n = 9

I.num = I.num.gen(n)
# EXACT INTEGRALS
I = I.gen(n)
#
I0 = I[[1]]
I8 = I[[2]]
I7 = I[[3]]
I4 = I[[4]]
I3 = I[[5]]
I6 = Ihalf.gen(I4)
I1 = Iinv.gen(I6)
# Different cycle: x^2 / ((x^3)^3 - 1)
I2 = Idx.gen(3, 3)
I5 = Ihalf.gen(I2)
# 6*x^8 - 3*(x^5 + x^2)
# I52 = 2*I8 - Ideriv.gen(n, pow)

### Test
low = 2
upper = 4

I8(low, upper)
I.num(low, upper, 8)

I7(low, upper)
I.num(low, upper, 7)

I6(low, upper)
I.num(low, upper, 6)

I5(low, upper)
I.num(low, upper, 5)

I4(low, upper)
I.num(low, upper, 4)

I3(low, upper)
I.num(low, upper, 3)

I2(low, upper)
I.num(low, upper, 2)

I1(low, upper)
I.num(low, upper, 1)

I0(low, upper)
I.num(low, upper, 0)


#####################

### Case n = 11
n = 11

I.num = I.num.gen(n)
# EXACT INTEGRALS
I = I.gen(n)
#
I0  = I[[1]]
I10 = I[[2]]
I9  = I[[3]]
I5  = I[[4]]
I4  = I[[5]]
I7 = Ihalf.gen(I4)
I2 = Iinv.gen(I7)
I6 = Ihalf.gen(I2)
I3 = Iinv.gen(I6)
I8 = Ihalf.gen(I6, sign.inverse=TRUE)
I1 = Iinv.gen(I8, sign.inverse=TRUE)
I.all = c(I10, I9, I8, I7, I6, I5, I4, I3, I2, I1, I0)

### Test
low = 2
upper = 4

for(i in 1:n) {
	cat("\nPow = "); cat(n - i); cat("\n")
	cat(I.all[[i]](low, upper))
	cat("\n")
	print(I.num(low, upper, n - i))
}

# TODO: loss of accuracy for I8 and I1 !


#####################

### Case n = 15
n = 15

I.num = I.num.gen(n)
# EXACT INTEGRALS
I = I.gen(n)
#
I0  = I[[1]]
I14 = I[[2]]
I13 = I[[3]]
I7  = I[[4]]
I6  = I[[5]]
I10 = Ihalf.gen(I6)
I3  = Iinv.gen(I10)
I12 = Ihalf.gen(I10)
I1  = Iinv.gen(I12, sign.inverse=TRUE, limits.inverse=TRUE)
# Different cycle: x^2/((x^3)^5 - 1)
I2  = Idx.gen(5, 3)
I11 = Iinv.gen(I2)
I8  = Ihalf.gen(I2)
I5  = Iinv.gen(I8)
# Different cycle: x^4/((x^5)^3 - 1)
I4  = Idx.gen(3, 5)
I9  = Iinv.gen(I4)
#
I.all = c(I14, I13, I12, I11, I10, I9, I8, I7, I6, I5, I4, I3, I2, I1, I0)

### Test
low = 2
upper = 4

for(i in 1:n) {
	cat("\nPow = "); cat(n - i); cat("\n")
	cat(I.all[[i]](low, upper))
	cat("\n")
	print(I.num(low, upper, n - i))
}

# TODO: loss of accuracy for I1 !


#########################
#########################
#########################

### Part C

#########################
###
### Polynomial Fractions
###   Composit Fractions
###
#########################

### P(x) / (x^n - 1)^2


#####################

### Case n = 3
n = 3

I.num = I.num.gen(n, 2)
# EXACT INTEGRALS
I = Idxfr.gen(n)
#
I0 = I[[1]]
I5 = I[[2]]
I3 = I[[3]]
I2 = I[[4]]
#
I1 = Iinv.gen(I3, sign.inverse=TRUE)
I4 = Iinv.gen(I0, sign.inverse=TRUE)


### Test
low = 2
upper = 4

I5(low, upper)
I.num(low, upper, 5)

I4(low, upper)
I.num(low, upper, 4)

I3(low, upper)
I.num(low, upper, 3)

I2(low, upper)
I.num(low, upper, 2)

I1(low, upper)
I.num(low, upper, 1)

I0(low, upper)
I.num(low, upper, 0)


#####################

### Case n = 5
n = 5

I.num = I.num.gen(n, 2)
# EXACT INTEGRALS
I = Idxfr.gen(n)
#
I0 = I[[1]]
I9 = I[[2]]
I5 = I[[3]]
I4 = I[[4]]
#
I8 = Iinv.gen(I0, sign.inverse=TRUE)
I3 = Iinv.gen(I5, sign.inverse=TRUE)
I2 = Ihalf.gen(I0, div=2)
I6 = Iinv.gen(I2, sign.inverse=TRUE)
I7 = I[[6]](3, I2)
I1 = Iinv.gen(I7, sign.inverse=TRUE)
#
I.all = c(I9, I8, I7, I6, I5, I4, I3, I2, I1, I0)


### Test
low = 2
upper = 4

pow = 2
for(i in 1:(pow*n)) {
	cat("\nPow = "); cat(pow*n - i); cat("\n")
	cat(I.all[[i]](low, upper))
	cat("\n")
	print(I.num(low, upper, pow*n - i))
}


#####################

### Case n = 7
n = 7

I.num = I.num.gen(n, 2)
# EXACT INTEGRALS
I = Idxfr.gen(n)
#
I0  = I[[1]]
I13 = I[[2]]
I7  = I[[3]]
I6  = I[[4]]
#
I12 = Iinv.gen(I0, sign.inverse=TRUE)
I5  = Iinv.gen(I7, sign.inverse=TRUE)
I3  = Ihalf.gen(I0, div=2)
I9  = Iinv.gen(I3, sign.inverse=TRUE)
I10 = I[[6]](4, I3)
I2  = Iinv.gen(I10, sign.inverse=TRUE)
I4  = Ihalf.gen(I2, div=2)
I8  = Iinv.gen(I4, sign.inverse=TRUE)
I11 = I[[6]](5, I4)
I1  = Iinv.gen(I11, sign.inverse=TRUE)
#
I.all = c(I13, I12, I11, I10, I9, I8, I7, I6, I5, I4, I3, I2, I1, I0)


### Test
low = 2
upper = 4

pow = 2
for(i in 1:(pow*n)) {
	cat("\nPow = "); cat(pow*n - i); cat("\n")
	cat(I.all[[i]](low, upper))
	cat("\n")
	print(I.num(low, upper, pow*n - i))
}




