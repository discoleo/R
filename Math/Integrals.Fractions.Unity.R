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


###########################

###############
### History ###
###############

### version 1 [RC1] [draft]
#
### 2023-01-14:
# - shortcut formulas moved to new file:
#   Integrals.Fractions.Unity.Definite.R;
### 2023-01-09:
# - shortcut formulas for I on [0, Inf];
### 2021-12-08:
# - cleanup: simplification of showcases;
### 2020-09-13:
# - cleanup: moved derived examples to:
#   Integrals.Fractions.Unity.Derived.R;
### 2020-05-20, 2020-06-03, 2020-06-04:
# - a new type of fractions: TODO: move to separate file;
# - shortcut for even powers (example);
# - only minor edit (unfortunately no time for more work);
### 2020-03-01
# - polynomial fractions: P(x) / (x^n - 1)^2
#   Cases: n=3, n=5, n=7;
### 2020-02-29
### 2020-02-28
### - polynomial fractions: P(x) / (x^n - 1)
###  -- Cases: n = 9, n = 15;
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
### TODO: move to: Integrals.Fractions.Unity.Powers.R;
### [are solvable using a different approach!]


### Publishing:
### - I am open to various projects,
###   including collaborations on articles/books/other projects;
### - Academic: it is easy nowadays to find Plagiarism,
###   but I am happy to get involved in projects related to these topics;


### Examples

### Fraction Decompositions

### ODD Powers:
n = 7 # e.g. 7, 9, 11;
### Roots of unity
# m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
len = (n-1)/2;
c1 = 2*cos(2*seq(len)*pi/n);
### Coefficients
b0 = 1/n;
b  = 2*b0; a = b0 * c1;
### Tests
x = 3 # e.g. 2, 3, pi, 4 # some arbitrary value for testing

### Partial Fractions
1/(x^n - 1) # ==
b0/(x - 1) + sum( (a*x - b) / (x^2 - c1*x + 1) )
1/n * (1/(x - 1) + sum( 2*(cs*x - 1) / (x^2 - 2*cs*x + 1) ))

###
1/(x^n + 1) # ==
b0/(x + 1) + sum( (a*x + b) / (x^2 + c1*x + 1) )
1/n * (1/ (x + 1) + sum( 2*(cs*x + 1) / (x^2 + 2*cs*x + 1) ))


### EVEN Powers:
# shortcut for even powers
# [was computed previously using the difference of the 2 lower powers]
# TODO: improve also code;
n = 10
n.2 = n/2 - 1;
# m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
c1 = 2*cos(2*seq(n.2)*pi/n)
b = -2/n; a = 1/n;
### Test
x = 2 # any value;
# used to test the fraction decomposition;
1/(x^n - 1) # ==
sum((a*c1*x + b) / (x^2 - c1*x + 1)) - b/(x^2-1)

###
x = 3^(1/5)
n = 3; # but for (4*n)!
c1 = 2*cos(seq(1, 4*n, by=2)*pi/(4*n));
1/(x^(4*n) + 1) # ==
sum((c1*x + 2) / (x^2 + c1*x + 1)) / (4*n)


### Integrals

### I( x^p / (x^n + 1) )
n = 7; # ODD Integer
p = 3; # Integer in [0, n-1]
lim = 4/5; # Arbitrary Interval: can be > 1;
integrate(\(x) x^p / (x^n + 1), 0, lim)
id = seq(2, n-1, by=2); x = lim; sg = if(p %% 2 == 0) 1 else -1;
cs = cos(id*pi/n); csp = cos(id*(p+1)*pi/n);
sn = sin(id*pi/n); snp = sin(id*(p+1)*pi/n);
sg/n * (log(x+1) +
	+ sum(csp*log(x^2 + 2*cs*x + 1) +
	+ 2*snp * (atan((x + cs)/sn) - atan(cs/sn))));

###
n = 8; # EVEN Integer: both 0 (mod 4) & 2 (mod 4);
p = 3; # Integer in [0, n-1]
lim = 4/5; # Arbitrary Interval: can be > 1;
integrate(\(x) x^p / (x^n + 1), 0, lim)
sg = if(p %% 2 == 0) 1 else -1;
id = seq(1, n, by=2); x = lim;
cs = cos(id*pi/n); csp = cos(id*(p+1)*pi/n);
sn = sin(id*pi/n); snp = sin(id*(p+1)*pi/n);
sg/n * sum(csp*log(x^2 + 2*cs*x + 1) +
	+ 2*snp * (atan((x + cs)/sn) - atan(cs/sn)));


### Infinite Integrals
# - shortcut formulas are available;

n = 5
integrate(function(x) 1/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^3/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi/n) / n
#
integrate(function(x) x/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^2/(x^n + 1), lower=0, upper=Inf)
pi / sin(2*pi/n) / n

###
n = 7
integrate(function(x) 1/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^5/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi/n) / n
#
integrate(function(x) x/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^4/(x^n + 1), lower=0, upper=Inf)
pi / sin(2*pi/n) / n
#
integrate(function(x) x^2/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^3/(x^n + 1), lower=0, upper=Inf)
pi / sin(3*pi/n) / n

### Note:
# - the remaining shortcut formulas to specific Definite Integrals
#   have been moved to file:
#   Integrals.Fractions.Unity.Definite.R;


########################
########################

##############
### Part A ###
##############

### Helper Functions ###

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
	} else if(n == 2) {
		if(type == "U-") {
			integrate.exact = function(lower, upper) {
				atan(upper) - atan(lower)
			}
		} else {
			integrate.exact = function(lower, upper) {
				log((upper-1)*(lower+1)/(upper+1)/(lower-1)) / 2
			}
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
# TODO: refactor & use directly trigonometric functions;
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
intNumeric = function(n, range, type="positive", poles=NULL, eps=1E-5) {
	type = pmatch(type, c("positive", "negative"));
	if(is.na(type)) stop("Type not supported!");
	k0 = if(type == 1) 1 else -1;
	if(is.null(poles) ||
			range[1] > poles[1] || range[2] < poles[1]) {
		r = integrate(function(x) 1/(x^n + k0), lower=range[1], upper=range[2]);
	} else {
		# TODO: multiple poles;
		r = integrate(function(x) 1/(x^n + k0), lower=range[1], upper=poles[1] - eps)$value;
		r = r + integrate(function(x) 1/(x^n + k0), lower=poles[1] + eps, upper=range[2])$value;
	}
	return(r)
}
intExact = function(n, range, type="positive") {
	type = pmatch(type, c("positive", "negative"));
	if(is.na(type)) stop("Type not supported!");
	len = if(n %% 2 == 1) (n-1)/2 else (n/2) - 1;
	c1  = 2*cos(2*seq(len)*pi/n);
	b0  = 1/n; b = -2*b0; a = b0 * c1;
	sh = c1/2; Dsq = sqrt(1 - sh^2);
	#
	if(n %% 2 == 1) {
		if(type == 1) {
			# 1 / (x^n + 1)
			lfr = (range[2]+1) / (range[1]+1);
			if(lfr < 0) lfr = lfr + 0i;
			r = b0 * log(lfr);
			r = sum(r, a/2*log((range[2]^2 + c1*range[2] + 1) /
				(range[1]^2 + c1*range[1] + 1)) );
			r = sum(r, -(a/2*c1 + b) / Dsq * atan((range[2] + sh)/Dsq),
					(a/2*c1 + b) / Dsq * atan((range[1] + sh)/Dsq));
		} else {
			# 1 / (x^n - 1)
			lfr = (range[2]-1) / (range[1]-1);
			if(lfr < 0) lfr = lfr + 0i;
			r = b0 * log(lfr);
			r = sum(r, a/2*log((range[2]^2 - c1*range[2] + 1) /
				(range[1]^2 - c1*range[1] + 1)) );
			r = sum(r, (a/2*c1 + b) / Dsq * atan((range[2] - sh)/Dsq),
					-(a/2*c1 + b) / Dsq * atan((range[1] - sh)/Dsq));
		}
	} else {
		# TODO
	}
	return(r)
}
### Exact Integral
# [old]
rez.f = function(x, a, b, b0, m.conj) {
	# [old]
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


### [new]
n = 5
rng = c(1.5, 3) # Range
intNumeric(n, rng)
intExact(n, rng)

n = 5
rng = c(1.5, 3)
intNumeric(n, rng, type="neg")
intExact(n, rng, type="neg")

### Poles
n = 5
rng = c(-1.5, 3)
intNumeric(n, rng, poles=-1)
intExact(n, rng) - 1i*pi/n # half Residue?

n = 5
rng = c(0, 3)
intNumeric(n, rng, type="neg", poles=1)
intExact(n, rng, type="neg") - 1i*pi/n # half Residue?


#########
### [old]
lower = 1.001
upper = 1024
f$integrate.all(lower, upper)
# 1.137622+0i
# 1.137622 with absolute error < 2e-14

###
lower = 1 + 1E-8
upper = 1024
f$integrate.all(lower, upper)
# 3.439807+0i
# 3.439807 with absolute error < 1.4e-11


### Roots of minus Unity
f = decompose.fr(5, type="U-")

lower = 1.001
upper = 1024
f$integrate.all(lower, upper)
# 0.1801464+0i
# 0.1801464 with absolute error < 2e-14

####################

###############
### Case n = 6:
# 1 / (x^6 - 1)

# uses a different type of decomposition!
# 1/(x^3 - 1) - 1/(x^3 + 1)
f = decompose.fr(6)

lower = 1.001
upper = 1024
f$integrate.all(lower, upper)
# 0.9053347+0i
# 0.9053347 with absolute error < 1.6e-14


lower = 1 + 1E-8
upper = 1024
f$integrate.all(lower, upper)
# 2.823739+0i
# 2.823739 with absolute error < 2.6e-12

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

f = decompose.fr(13)

lower = 1.001
upper = 1024
f$integrate.all(lower, upper)
# 0.3448444+0i
# 0.3448444 with absolute error < 6.5e-15


lower = 1 + 1E-8
upper = 1024
f$integrate.all(lower, upper)
# 1.229993+0i
# 1.229993 with absolute error < 4.5e-12


################
### Case n = 16:
# 1 / (x^16 - 1)

# uses a different type of decomposition
f = decompose.fr(16)

lower = 1.001
upper = 1024
f$integrate.all(lower, upper)
# 0.2656526+0i
# 0.2656526 with absolute error < 5.2e-15


lower = 1 + 1E-8
upper = 1024
f$integrate.all(lower, upper)
# 0.9847423+0i
# 0.9847423 with absolute error < 2.8e-12


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
# Numeric: ERROR
# possible with lower upper limit;
upper = 4


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
upper = 4; # OK - same result: 0.02335263;

lower = 1 + 1E-8
upper = 1024
f$integrate.all(lower, upper)
# 0.1368511+0i
# ERROR: wide discrepancy:
# 5.118336e-33 with absolute error < 1e-32


################

### Even Powers: Direct Formula


### Case: 2*(2*n+1)
n = 6
x = sqrt(3/5);
cs = 2*cos(seq(1,n, by=2)*pi/n);
div = n;
1/(x^n + 1)
r = 2/(x^2 + 1) + (cs[1]*x + 2)/(x^2+cs[1]*x+1) + (cs[3]*x + 2)/(x^2 + cs[3]*x+1)
r / div


###
n = 10
x = sqrt(3/5);
cs = 2*cos(seq(1,n, by=2)*pi/n);
div = n;
1/(x^n + 1)
r = 2/(x^2 + 1) + (cs[1]*x + 2)/(x^2+cs[1]*x+1) + (cs[2]*x + 2)/(x^2 + cs[2]*x+1) +
	+ (cs[4]*x + 2)/(x^2+cs[4]*x+1) + (cs[5]*x + 2)/(x^2 + cs[5]*x+1)
r / div


### Case: 4*n
n = 4
x = sqrt(3/5);
cs = 2*cos(seq(1,n, by=2)*pi/n);
1/(x^n + 1)
r = (cs[1]*x + 2)/(x^2+cs[1]*x+1) + (cs[2]*x + 2)/(x^2 + cs[2]*x+1)
r / n

###
n = 8
x = sqrt(3/5);
cs = 2*cos(seq(1,n, by=2)*pi/n);
1/(x^n + 1)
r = (cs[1]*x + 2)/(x^2+cs[1]*x+1) + (cs[2]*x + 2)/(x^2+cs[2]*x+1) +
	+ (cs[3]*x + 2)/(x^2 + cs[3]*x+1) + (cs[4]*x + 2)/(x^2 + cs[4]*x+1)
r / n

###
n = 16
x = sqrt(3/5);
cs = 2*cos(seq(1,n, by=2)*pi/n);
1/(x^n + 1)
r = (cs[1]*x + 2)/(x^2+cs[1]*x+1) + (cs[2]*x + 2)/(x^2+cs[2]*x+1) +
	+ (cs[3]*x + 2)/(x^2 + cs[3]*x+1) + (cs[4]*x + 2)/(x^2 + cs[4]*x+1) +
	+ (cs[5]*x + 2)/(x^2 + cs[5]*x+1) + (cs[6]*x + 2)/(x^2 + cs[6]*x+1) +
	+ (cs[7]*x + 2)/(x^2 + cs[7]*x+1) + (cs[8]*x + 2)/(x^2 + cs[8]*x+1)
r / n


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
		if(Re(low) < 0 || Re(upper) < 0) {
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
	return(Ihalf)
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

### Case n = 13
n = 13

# 0 => 11  6 => 5  9 => 2 => 7 => 4 => 8 => 3 10 => 1

I.num = I.num.gen(n)
# EXACT INTEGRALS
I = I.gen(n)
#
I0  = I[[1]]
I12 = I[[2]]
I11 = I[[3]]
I6  = I[[4]]
I5  = I[[5]]
I9  = Ihalf.gen(I6)
I2  = Iinv.gen(I9)
I7  = Ihalf.gen(I2)
I4  = Iinv.gen(I7)
I8  = Ihalf.gen(I4, sign.inverse=TRUE)
I3  = Iinv.gen(I8, sign.inverse=TRUE)
I10 = Iinv.gen(I8, sign.inverse=TRUE)
I1  = Iinv.gen(I10)
I.all = c(I12, I11, I10, I9, I8, I7, I6, I5, I4, I3, I2, I1, I0)

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


##################################
##################################

### Sequence Generation:
n = 101

seq.df = data.frame(id=1:(n-2), prev=NA, type=NA)
levels(seq.df$type) = c("half", "inv")

seq.start = 0
pow.seq = c(seq.start, rep(NA, n - 2))
pos = 1

while(pos <= length(pow.seq)) {
	current = pow.seq[pos]
	pos = pos + 1
	# Inverse
	inv = (n - current - 2)
	cat(inv); cat(", ");
	isAdded = FALSE;
	if(inv > 0 && is.na(seq.df$prev[inv])) {
		seq.df$prev[inv] = current;
		seq.df$type[inv] = "inv";
		# if(inv %% 2 == 0) {
		if(inv > 0) {
			pow.seq[pos] = inv
			isAdded = TRUE;
		}
	}
	# Half
	if(current %% 2 == 1) { next; }
	half = (n + current - 1) / 2
	if(half > 0 && is.na(seq.df$prev[half])) {
		seq.df$prev[half] = current;
		seq.df$type[half] = "half";
		if(half > 0) {
			pow.seq[ifelse(isAdded, pos+1, pos)] = half
		}
	}
}

table(is.na(seq.df$prev))
table(seq.df$type)
length(unique(pow.seq))

# [1]  0 99 50 49 75 24 62 37 81 18 59 40 70 29 85 14 57 42 71 28 64 35 82 17 91  8 54 45
# [29] 77 22 61 38 69 30 65 34 67 32 66 33 83 16 58 41 79 20 60 39 80 19 90  9 95  4 52 47
# [57] 76 23 88 11 94  5 97  2 51 48 74 25 87 12 56 43 78 21 89 10 55 44 72 27 86 13 93  6
# [85] 53 46 73 26 63 36 68 31 84 15 92  7 96  3 98  1


### Sequence Generation:
n = 10403
# outside of main cycle: 202 numbers ;-)
# [smiling Diffie]



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


#####################

### Case n = 11
n = 11

I.num = I.num.gen(n, 2)
# EXACT INTEGRALS
I = Idxfr.gen(n)
#
I0  = I[[1]]
I21 = I[[2]]
I11 = I[[3]]
I10 = I[[4]]
#
I20 = Iinv.gen(I0, sign.inverse=TRUE)
I9  = Iinv.gen(I11, sign.inverse=TRUE)
I5  = Ihalf.gen(I0, div=2)
I15 = Iinv.gen(I5, sign.inverse=TRUE)
I16 = I[[6]](6, I5) # n + pow - 1
I4  = Iinv.gen(I16, sign.inverse=TRUE)
I7  = Ihalf.gen(I4, div=2)
I13 = Iinv.gen(I7, sign.inverse=TRUE)
I18 = I[[6]](8, I7) # n + pow - 1
I2  = Iinv.gen(I18, sign.inverse=TRUE)
I6  = Ihalf.gen(I2, div=2)
I14 = Iinv.gen(I6, sign.inverse=TRUE)
I8  = Ihalf.gen(I6, sign.inverse=TRUE, div=-2)
I12 = Iinv.gen(I8, sign.inverse=TRUE)
I19 = I[[6]](9, I8) # n + pow - 1
I1  = Iinv.gen(I19, sign.inverse=TRUE)
I17 = I[[6]](7, I6) # n + pow - 1
I3  = Iinv.gen(I17, sign.inverse=TRUE)
#
I.all = c(I21, I20, I19, I18, I17, I16, I15, I14, I13, I12, I11, I10,
	I9, I8, I7, I6, I5, I4, I3, I2, I1, I0)


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

# TODO: ERRORS:
# I8 => I19; needs sign.inverse;
# I1 ???



##########################
##########################

### Derived Examples:
# - moved to:
#   Integrals.Fractions.Unity.Derived.R;
# - TODO: remove from here;

lower = 2
upper = 3
integrate(function(x) 1/2 * x^3/(x^5 + 1), lower=lower^2, upper=upper^2)
integrate(function(x) x^7/(x^10 + 1), lower=lower, upper=upper)
integrate(function(x) sin(x)^7*cos(x)/(sin(x)^10 + cos(x)^10), lower=atan(lower), upper=atan(upper))
integrate(function(x) 1/2 * x^3/(x^5 + (1-x)^5), lower=sin(atan(lower))^2, upper=sin(atan(upper))^2)
integrate(function(x) 1/2 * x^3/(x^5 + (1-x)^5), lower=1-1/(lower^2+1), upper=1-1/(upper^2+1))
integrate(function(x) (1+x)^3/((1+x)^5 + (1-x)^5), lower=1-2/(lower^2+1), upper=1-2/(upper^2+1))

###
n = 12
#
lower = 2
upper = 3
#
n.half = n/2
integrate(function(x) 1/2 * x^(n.half-2)/(x^n.half + 1), lower=lower^2, upper=upper^2)
integrate(function(x) x^(n-3)/(x^n + 1), lower=lower, upper=upper)
integrate(function(x) sin(x)^(n-3)*cos(x)/(sin(x)^n + cos(x)^n), lower=atan(lower), upper=atan(upper))
integrate(function(x) 1/2 * x^(n.half - 2)/(x^n.half + (1-x)^n.half), lower=sin(atan(lower))^2, upper=sin(atan(upper))^2)
integrate(function(x) 1/2 * x^(n.half - 2)/(x^n.half + (1-x)^n.half), lower=1-1/(lower^2+1), upper=1-1/(upper^2+1))
integrate(function(x) (1+x)^(n.half - 2)/((1+x)^n.half + (1-x)^n.half), lower=1-2/(lower^2+1), upper=1-2/(upper^2+1))
