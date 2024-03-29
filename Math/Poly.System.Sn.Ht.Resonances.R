########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Heterogeneous Symmetric:
### Resonances
###
### draft v.0.1L


### Resonances in Polynomial Systems
# System: 3..n Variables
# Type: Heterogeneous Symmetric

### Description:
# - this work started as an exploration of resonances
#   in polynomial systems;
# - however, it is more widely applicable to Matrix Theory;
# - Class 2 polynomials play an important role
#   in Matrix Theory;


# Note:
# - moved to this new file from file:
#   Poly.System.Hetero.Symmetric.S3.Mixed.R;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")
source("Polynomials.Helper.Matrix.R")


### Class 2 Poly:
# - calculate the Class 2 polynomial;
# - unity = only a placeholder;
# TODO:
# - proper name: expand.u vs poly.unity;
# - Optimization: replace.pm is very slow;
expand.u.pm = function(x, n, unity="u") {
	# TODO: n = non-prime;
	if(n < 2) stop("Invalid n!");
	if(n == 2) return(x);
	#
	pR  = x;
	tmp = x; tmp[, unity] = (tmp[, unity] * (n - 1)) %% n;
	pR  = mult.pm(pR, tmp);
	pR[, unity] = pR[, unity] %% n;
	# quasi-optimization;
	isEven = (n %% 2) == 0;
	if(isEven) {
		tmp = replace.pm(x, -1, unity);
		pR  = mult.pm(pR, tmp);
	}
	n2b = (n - 1) %/% 2;
	if(n > 4) for(i in seq(2, n2b)) {
		tmp1 = x;
		tmp1[, unity] = (tmp1[, unity] * i) %% n;
		tmp2 = x;
		tmp2[, unity] = (tmp2[, unity] * (n - i)) %% n;
		tmp = mult.pm(tmp1, tmp2);
		tmp[, unity] = tmp[, unity] %% n;
		pR = mult.pm(pR, tmp);
		pR[, unity] = pR[, unity] %% n;
	}
	# Simplify: SLOW !!!
	if(isEven) {
		pR = replace.pm(pR, -1, unity, pow = n2b + 1);
		hasUnity = ! is.na(match(unity, names(pR)));
		if(hasUnity) {
			# tricky: alternating;
			if(n2b %% 2 == 1) {
				pInv = data.frame(u = seq(0, n2b - 2), coeff=c(-1, 1));
				pInv = rbind(pInv, c(n2b - 1, -1));
			} else
				pInv = data.frame(u = seq(0, n2b - 1), coeff=c(-1, 1));
			names(pInv)[1] = unity;
			pR = replace.pm(pR, pInv, unity, pow=n2b);
		}
	} else {
		pInv = data.frame(u = seq(0, n-2), coeff=-1);
		names(pInv)[1] = unity;
		pR = replace.pm(pR, pInv, unity, pow=n-1);
	}
	# Test:
	hasUnity = ! is.na(match(unity, names(pR)));
	if(hasUnity && any(pR[, unity] != 0)) {
		warning("Something went wrong!");
	}
	pR = drop.pm(pR);
	if(! inherits(pR, "pm")) class(pR) = c("pm", class(pR));
	return(pR);
}

### Resonances
# n = number of variables;
resonance = function(p, n=3) {
	# "Resonance" with Roots of unity;
	# currently only for 2-variable terms: x^p1*y^p2;
	sg = if(n %% 2 == 1) 1 else -1;
	p.all = p[1]^n + sg*p[2]^n;
	r = list(p=p.all, f=factors(p.all), p.trivial = sum(p));
	return(r);
}

diag.T3 = function(x, n) {
	if(length(x) != 3) {
		stop("Number of powers is not 3!");
	}
	m = diag(x[1], n);
	for(nc in seq(n-1)) m[nc, nc + 1] = x[2];
	for(nc in seq(n-2)) m[nc, nc + 2] = x[3];
	m[n-1, 1] = x[3];
	m[n, 1] = x[2]; m[n, 2] = x[3];
	return(m);
}
diag.T4 = function(x, n) {
	if(length(x) != 4) {
		stop("Number of powers is not 4!");
	}
	m = diag(x[1], n);
	for(nc in seq(n-1)) m[nc, nc + 1] = x[2];
	for(nc in seq(n-2)) m[nc, nc + 2] = x[3];
	for(nc in seq(n-3)) m[nc, nc + 3] = x[4];
	m[n-2, 1] = x[4];
	m[n-1, 1] = x[3]; m[n-1, 2] = x[4];
	m[n, 1] = x[2]; m[n, 2] = x[3]; m[n, 3] = x[4];
	return(m);
}

### Generators
# n = number of variables;
polyRes3 = function(n) {
	p = det.mpm(diag.lpm(c("k", "n", "p"), n=n));
	p = sort.pm(p, c("k", "n"));
	if( ! inherits(p, "pm")) class(p) = c("pm", class(p));
	return(p);
}
polyRes4 = function(n) {
	p = det.mpm(diag.lpm(c("k1", "n1", "n2", "k2"), n=n));
	p = sort.pm(p, c("k1", "n1", "n2", "k2"));
	if( ! inherits(p, "pm")) class(p) = c("pm", class(p));
	return(p);
}

### Test

# n = number of variables;
test.res.T4 = function(k, pow=c(2,1,1,1), n=length(k)) {
	sapply(seq(0, n-1), function(i) {
		id = c(i, i+1, i+2, i+3);
		id = (id %% n) + 1;
		sum(k[id] * pow);
	});
}
test.res.T3 = function(k, pow=c(3,1,1), n=length(k)) {
	sapply(seq(0, n-1), function(i) {
		id = c(i, i+1, i+2);
		id = (id %% n) + 1;
		sum(k[id] * pow);
	});
}
test.res.T2 = function(k, n, pow) {
	sapply(seq(0, n-1), function(i) {
		id = c(i, i+1);
		id = (id %% n) + 1;
		sum(k[id] * pow);
	});
}


#####################
#####################

##################
### Resonances ###
##################

######################
### Roots of Unity ###
######################

### Sys[3]: 3 Variables
# - if (x, y, z) is a solution, then:
#   (x, y, z) * m^(k1, k2, k3) is also a solution;


### Derivation:

### Order: k + n
# - Base-Eq with terms: x^k*y^n, y^k*z^n, z^k*x^n
# - roots of unity: m^p = 1;
# - reusing: "x" = power of m in x;

k*x + n*y # = 0 (mod p);
k*y + n*z # = 0 (mod p);
k*z + n*x # = 0 (mod p);

### Sum =>
(k+n)*(x+y+z) # = 0 (mod p);

### Reduction =>
(n^2 + k^2 - n*k)*x # = 0 (mod p);

### Trivial Powers
# p = (k+n) # or one of its Divisors;

### Non-Trivial Powers
# p = (n^2 + k^2 - n*k) # or one of its Divisors;

### Combinations:
# p = Divisors of (k+n)*(n^2 + k^2 - n*k)

#############
### Examples:

### Note:
# - entire system has to be compatible with these "rotations",
#   if they should represent valid transformations;

#######
### Ex:
# k = 2; n = 1;
# Base Monomials: (x^2*y, y^2*z, z^2*x);
### Trivial:
p = 3;
### Non-Trivial & Combinations:
p = c(3, 9); # p[non-trivial] is the same!
# ex: p = 9 =>
# new solution: (x,y,z) * (m^1, m^7, m^4);
# new solution: (x,y,z) * (m^2, m^5, m^8);


#######
### Ex:
# k = 3; n = 1;
### Trivial:
p = 4;
### Non-Trivial & Combinations:
p = c(7, 14, 28);
# ex: p = 7 =>
# new solution: (x,y,z) * (m^1, m^4, m^2);
# ex: p = 14 =>
# new solution: (x,y,z) * (m^1, m^11, m^9);
# ex: p = 28 =>
# new solution: (x,y,z) * (m^1, m^25, m^9);


#######
### Ex:
# k = 3; n = 2;
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(7, 35);
# ex: p = 7 =>
# new solution: (x,y,z) * (m^1, m^2, m^4);
# ex: p = 35 =>
# new solution: (x,y,z) * (m^1, m^16, m^11);


#####################
### Dual Resonances
### Conjugate Powers:
# (n^2 + k^2 - n*k) is the same!

###
# (n, 1) => p | (n^2 - n + 1);
# (n, n-1) => p | (n^2 - n + 1);

###
# (n, i) => p | (n^2 + i^2 - n*i);
# (n, n - i) => p | (n^2 + i^2 - n*i);

### Examples:

### Ex 1:
p = 7;
# n = c(3, 2) & c(1, 3);
# New solutions:
# (x,y,z) * (m^1, m^2, m^4);
# (x,y,z) * (m^3, m^6, m^5);
# - see a solved example in file:
#   Poly.System.Hetero.Symmetric.S3.Mixed.Dual.R;


### Ex 2:
p = 13;
# n = c(4, 3) & c(1, 4);
# New solutions:
# (x,y,z) * (m^1, m^3, m^9);
# (x,y,z) * (m^2, m^6, m^5);
# (x,y,z) * (m^4, m^12, m^10);
# (x,y,z) * (m^7, m^8, m^11);


######################
######################

######################
### Roots of Unity ###
### Generalization ###
######################

### Sys[i]:

### System with [i] Variables
# - for i = any odd;
### Order: k + n
# - Base-Eq with terms:
#   x[1]^k*x[2]^n, x[2]^k*x[3]^n, ..., x[i-1]^k*x[i]^n, x[i]^k*x[1]^n;

### All Powers
# p = Divisors of (k^i + n^i)


#######
### Ex:
# i = 5; k = 2; n = 1;
# p = Divisors(33);
# Base Monomials: (x1^2*x2, ..., x5^2*x1);
### Trivial:
p = 3;
### Non-Trivial & Combinations:
p = c(11, 33);
# ex: p = 11 =>
# new solutions:
# (x1,x2,x3,x4,x5) * (m^1, m^9, m^4, m^3, m^5);
# (x1,x2,x3,x4,x5) * (m^2, m^7, m^8, m^6, m^10);


#######
### Ex:
# i = 5; k = 3; n = 1;
# p = Divisors(244);
### Trivial:
p = 4;
### Non-Trivial & Combinations:
p = c(61, 122, 244);
# ex: p = 61 =>
# new solution:
# (x1,x2,x3,x4,x5) * (m^1, m^58, m^9, m^34, m^20);


#######
### Ex:
# i = 5; k = 3; n = 2;
# p = Divisors(275);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(11, 25, 55, 275);
# ex: p = 11 =>
# new solution:
# (x1,x2,x3,x4,x5) * (m^1, m^4, m^5, m^9, m^3);
# ex: p = 25 =>
# new solution:
# (x1,x2,x3,x4,x5) * (m^1, m^11, m^21, m^6, m^16);


###################

### System with [i] Variables
# - for i = even;
### Order: k + n
# - Base-Eq with terms:
#   x[1]^k*x[2]^n, x[2]^k*x[3]^n, ..., x[i-1]^k*x[i]^n, x[i]^k*x[1]^n;

### All Powers
# p = Divisors of (k^i - n^i)


#######
### Ex:
# i = 4; k = 2; n = 1;
# p = Divisors(15);
### Trivial:
p = 3;
### Non-Trivial & Combinations:
p = c(5, 15);
# ex: p = 5 =>
# new solution:
# (x1,x2,x3,x4) * (m^1, m^3, m^4, m^2);


#######
### Ex:
# i = 4; k = 3; n = 1;
# p = Divisors(80);
### Trivial:
p = 4;
### Non-Trivial & Combinations:
p = c(5, 8, 10); # and higher
# ex: p = 5 =>
# new solution:
# (x1,x2,x3,x4) * (m^1, m^2, m^4, m^3);
# ex: p = 8 =>
# new solution:
# (x1,x2,x3,x4) * (m^1, m^5, m^1, m^5);


#######
### Ex:
# i = 4; k = 3; n = 2;
# p = Divisors(65);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(13, 65);
# ex: p = 13 =>
# new solution:
# (x1,x2,x3,x4) * (m^1, m^5, m^12, m^8);


#######
### Ex:
# i = 4; k = 4; n = 1;
# p = Divisors(255);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(3, 17); # and higher
# ex: p = 17 =>
# new solution:
# (x1,x2,x3,x4) * (m^1, m^13, m^16, m^4);


#######
### Ex:
# i = 6; k = 2; n = 1;
resonance(c(2,1), n=6)
# p = Divisors(63);
### Trivial:
p = 3;
### Non-Trivial & Combinations:
p = c(7, 9, 21, 63);
# ex: p = 7 =>
# new solution:
# (x1,x2,...,x6) * (m^1, m^5, m^4, m^6, m^2, m^3);
# ex: p = 9 =>
# new solution:
# (x1,x2,...,x6) * (m^1, m^7, m^4, m^1, m^7, m^4);


####################

################
### 3 Powers ###
################

### Sys: 4 Variables

### Order (k,n,p)
# p = Divisors(k^4 + p^4 - n^4 + 4*k*p*n^2 - 2*k^2*p^2)
#   = Divisors((k+n+p) * (...));

divisors.S4T3 = function(n) {
	# Note: (k+n+p) is omitted;
	k = n[1]; p = n[3]; n = n[2];
	d = k^3 + p^3 - n^3 - k*p*(k + p) - n*(k - p)^2 + n^2*(k + p);
	return(d);
}

### Special Cases:
### k = p;
# d = (2*k + n)*(2*k - n)*n^2;
### n = k + p;
# d = 0;
# - equivalent to: (k,k,0) + (0,p,p)
# => trivial Solution: (m^j, m^(pp-j), m^j, m^(pp-j));


### Order: 2+1+1
### 4 Variables
# i = 4; p = c(2,1,1)
# p = Divisors(16);
### Trivial:
p = 4;
### Non-Trivial & Combinations:
p = c(8, 16); # also = 4;
# - there are various quasi-non-trivial solutions;
### Ex: p = 8 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^3, m^7, m^3, m^7);
# (3,7,3), (7,3,7), (3,7,3), (7,3,7)


### Order: 3+1+1
### 4 Variables
# i = 4; p = c(3,1,1)
# p = Divisors(75);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(3, 15, 25, 75, 125, 375);
# - there are also quasi-non-trivial solutions for p = 5;
# - p = 3 reduces to an E2-type;
### Ex: p = 3 =>
# "non-trivial" solution:
# (x1,x2,x3,x4) * (m^1, m^2, m^1, m^2);
# (1,2,1), (2,1,2), (1,2,1), (2,1,2)
### Ex: p = 5 =>
# "non-trivial" solutions:
# (x1,x2,x3,x4) * (m^1, m^2, m^0, m^4);
# (1,2,0), (2,0,4), (0,4,1), (4,1,2)
# (x1,x2,x3,x4) * (m^3, m^2, m^4, m^0);
# (3,2,4), (2,4,0), (4,0,3), (0,3,2)
### Ex: p = 15 =>
# "non-trivial" solution:
# (x1,x2,x3,x4) * (m^1, m^11, m^1, m^11);
# (1,11,1), (11,1,11), (1,11,1), (11,1,11)


### Order: 4+1+1
### 4 Variables
# i = 4; p = c(4,1,1)
# p = Divisors(6*40);
### Trivial:
p = 6;
### Non-Trivial & Combinations:
p = c(4,5,8); # and many more
# - there are various quasi-non-trivial solutions;
### Ex: p = 5 =>
# "non-trivial" solution:
# (x1,x2,x3,x4) * (m^1, m^2, m^4, m^3);
# (1,2,4), (2,4,3), (4,3,1), (3,1,2)
### Ex: p = 8 =>
# "non-trivial" solution:
# (x1,x2,x3,x4) * (m^5, m^7, m^5, m^7);
# (5,7,5), (7,5,7), (5,7,5), (7,5,7)


### Special Cases:

### Order: 2+3+1
# d = 0
p = 5;
# (x1,x2,x3,x4) * (m^1, m^2, m^2, m^0);
# (1,2,2), (2,2,0), (2,0,1), (0,1,2)
# (x1,x2,x3,x4) * (m^2, m^4, m^4, m^0);
# (2,4,4), (4,4,0), (4,0,2), (0,2,4)
# [*Trivial*] (x1,x2,x3,x4) * (m^1, m^4, m^1, m^4);
# (1,4,1), (4,1,4), (1,4,1), (4,1,4)
p = 7;
# (x1,x2,x3,x4) * (m^1, m^6, m^1, m^6);
# (1,6,1), (6,1,6), (1,6,1), (6,1,6)
# (x1,x2,x3,x4) * (m^2, m^5, m^2, m^5);
# (2,5,2), (5,2,5), (2,5,2), (5,2,5)


####################
####################

### Sys: 5 Variables

### Order: (k,n,p)
# p = Divisors(k^5 + n^5 - 5*k*n^3*p + 5*k^2*n*p^2 + p^5)
# Note: unfortunate repetition of "p";

divisors.S5T3 = function(n) {
	# Note: (k+n+p) is omitted;
	k = n[1]; p = n[3]; n = n[2];
	k^4 + n^4 + p^4 - k^3*n - k^3*p - k*n^3 - n^3*p - k*p^3 - n*p^3 +
		+ k^2*n^2 + k^2*p^2 + n^2*p^2 + k*n*p*(2*(k + p) - 3*n);
}


### Special Cases:

### Order: (p1,1,1)
# p = Divisors(p1^5 + 5*p1^2 - 5*p1 + 2)
#   = (p1 + 2) * (p1^4 - 2*p1^3 + 4*p1^2 - 3*p1 + 1)
# - Note: finding the explicit/individual powers is a headache;

### Order: 2+1+1
### 5 Variables
# i = 5; p = c(2,1,1)
# p = Divisors(88);
### Trivial:
p = 4;
### Non-Trivial & Combinations:
p = c(11, 22, 44);
# ex: p = 11 =>
# new solutions:
# (x1,x2,...,x5) * (m^1, m^4, m^5, m^9, m^3);
# (1,4,5), (4,5,9), (5,9,3), (9,3,1), (3,1,4)
# (x1,x2,...,x5) * (m^2, m^8, m^10, m^7, m^6);
# ex: p = 22 =>
# new solution:
# (x1,x2,...,x5) * (m^1, m^15, m^5, m^9, m^3);
# (1,15,5), (15,5,9), (5,9,3), (9,3,1), (3,1,15)
k = c(1,4,5,9,3);
k = c(2,8,10,7,6);
test.res.T3(k, pow=c(2,1,1))


### Order: 3+1+1
# i = 5; p = c(3,1,1)
# p = Divisors(25*11)
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(11); # & combinations;
# ex: p = 11 =>
# new solutions:
# (x1,x2,...,x6) * (m^1, m^5, m^3, m^4, m^9);
# (1,5,3), (5,3,4), (3,4,9), (4,9,1), (9,1,5);
# (x1,x2,...,x6) * (m^2, m^10, m^6, m^8, m^7);
k = c(1,5,3,4,9);
k = c(2,10,6,8,7);
test.res.T3(k, pow=c(3,1,1))


### 6 Variables
# i = 6; p = c(3,1,1)
# p = Divisors(720);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(9); # possible others
# ex: p = 9 =>
# new solution:
# (x1,x2,...,x6) * (m^1, m^2, m^4, m^8, m^7, m^5);
# (1,2,4), (2,4,8), (4,8,7), (8,7,5), (7,5,1), (5,1,2)


####################

### Sys: 5 Variables

### Order: c(p1, p2, p1)
# p = Divisors((2*p1 + p2) * (p1^2 + p1*p2 - p2^2)^2);

### Special Case: c(p1, 1, p1)
# p = Divisors((2*p1 + 1) * (p1^2 + p1 - 1)^2);

### Ex 1: c(2,1,2)
# p = Divisors(125);


### Ex 2: c(3,1,3)
# p = Divisors(7 * 11^2);


### Ex 3: c(4,5,4)
# p = Divisors(13 * 11^2);


##########
### Order: c(p1, p2, 1)
# p = Divisors(p1^5 + p2^5 + 5*p1^2*p2 - 5*p1*p2^3 + 1)
#   = (p1 + p2 + 1) * (...);
divisors.S5T3 = function(p) {
	# assumes p[3] == 1;
	p1 = p[1]; p2 = p[2];
	d = p1^4 + p2^4 - p1^3*p2 - p1*p2^3 + p1^2*p2^2 - p1^3 - p2^3 + 2*p1^2*p2 - 3*p1*p2^2 +
		+ p1^2 + p2^2 + 2*p1*p2 - p1 - p2 + 1;
	return(d);
}


######################
######################

### Note:
# - the polynomials are always divisible by (k + n + p),
#   but the number of monomials in the result rises rapidly;
# - number monomials = 28 (for v = 7);
# - number monomilas = 55 (for v = 10);
# - number monomials = 91 (for v = 13);


### Sys: 7 variables

### Order: c(k,n,p)
# p = Divisors(k^7 + n^7 + p^7 - 7*k*n^5*p - 7*k^3*n*p^3 + 14*k^2*n^3*p^2)


x = 2*cos(pi*c(1,3,5)/7) + 2;
x^3 - 7*x^2 + 14*x - 7 # = 0


######################
######################

### Sys: 9 variables

### Order: c(k,n,p)
# p = Divisors(k^9 + n^9 + p^9 - 9*k*n^7*p + 9*k^4*n*p^4 + 27*k^2*n^5*p^2 - 30*k^3*n^3*p^3)


x = 2*cos(pi*c(1,3,5,7)/9) + 2;
x^4 - 9*x^3  + 27*x^2 - 30*x + 9 # = 0


######################
######################

### Sys: 11 variables

### Order: c(k,n,p)
# p = Divisors(k^11 + n^11 + p^11 - 11*k*n^9*p - 11*k^5*n*p^5 + 44*k^2*n^7*p^2 - 77*k^3*n^5*p^3 + 55*k^4*n^3*p^4)


x = 2*cos(pi*c(1,3,5,7,9)/11) + 2;
x^5 - 11*x^4 + 44*x^3 - 77*x^2 + 55*x - 11 # = 0


######################
######################

### Sys: 13 variables

### Order: c(k,n,p)
# p = Divisors(
#		k^13 + n^13 + p^13 - 13*k*n^11*p + 65*k^2*n^9*p^2 - 156*k^3*n^7*p^3 + 182*k^4*n^5*p^4 +
#		- 91*k^5*n^3*p^5 + 13*k^6*n*p^6)

# Note: takes some seconds to compute!
# det.mpm(diag.lpm(c("k", "n", "p"), n=13))

x = 2*cos(pi*c(seq(1,11, by=2))/13) + 2;
x^6 - 13*x^5 + 65*x^4 - 156*x^3 + 182*x^2 - 91*x + 13 # = 0


######################
######################

### Even Number of Vars:

### Note:
# - formula is still applicable;
# - only the "n"-leading term has a changed sign:
#   k^v - n^v + p^v, or n^v - k^v - p^v;
# - the polynomials are always divisible by (k + n + p),
#   but the number of monomials in the result rises rapidly;


### Sys: 6 variables

### Order: c(k,n,p)
# p = Divisors(k^6 - n^6 + p^6 + 6*k*n^4*p - 9*k^2*n^2*p^2  + 2*k^3*p^3)

# x is relative to exponent variable n;
x = 2*cos(pi*c(seq(1,5, by=2))/6) + 2;
x^3 - 6*x^2 + 9*x - 2 # = 0


#######################
#######################

### Sys: 8 variables

### Order: c(k,n,p)
# p = Divisors(k^8 - n^8 + p^8 + 8*k*n^6*p - 20*k^2*n^4*p^2 + 16*k^3*n^2*p^3 - 2*k^4*p^4)

# x is relative to exponent variable n;
x = 2*cos(pi*c(seq(1,7, by=2))/8) + 2;
x^4 - 8*x^3 + 20*x^2 - 16*x + 2 # = 0
x^3 - 6*x^2 + 9*x - 2 # = 0


######################
######################

### Sys: 10 variables

### Order: c(k,n,p)
# p = Divisors(
#		k^10 - n^10 + p^10 + 10*k*n^8*p - 35*k^2*n^6*p^2 + 50*k^3*n^4*p^3 - 25*k^4*n^2*p^4 + 2*k^5*p^5)

x = 2*cos(pi*c(seq(1,9, by=2))/10) + 2;
x^5 - 10*x^4 + 35*x^3 - 50*x^2 + 25*x - 2 # = 0


######################
######################

### Sys: 12 variables

### Order: c(k,n,p)
# p = Divisors(
#		k^12 - n^12 + p^12 + 12*k*n^10*p - 54*k^2*n^8*p^2 + 112*k^3*n^6*p^3 - 105*k^4*n^4*p^4 +
#			+ 36*k^5*n^2*p^5 - 2*k^6*p^6)

x = 2*cos(pi*c(seq(1,11, by=2))/12) + 2;
x^6 - 12*x^5 + 54*x^4 - 112*x^3 + 105*x^2 - 36*x + 2 # = 0

### Div:
p = toPoly.pm("k^12 - n^12 + p^12 + 12*k*n^10*p - 54*k^2*n^8*p^2 + 112*k^3*n^6*p^3 - 105*k^4*n^4*p^4 +
		+ 36*k^5*n^2*p^5 - 2*k^6*p^6");
pR = replace.pm(p, toPoly.pm("- k - p"), "n")
nrow(pR)

# pR = div.pm(p, "k + n + p", "n")


####################
####################
####################

################
### 4 Powers ###
################

### 5 Variables:

polyRes4(5)

k1^5 + n1^5 + n2^5 + k2^5 - 5*k1^3*k2*n2 - 5*k1*k2^3*n1 - 5*k1*n1^3*n2 - 5*k2*n1*n2^3 +
	+ 5*k1^2*k2*n1^2 + 5*k1*k2^2*n2^2 + 5*k1^2*n1*n2^2 + 5*k2^2*n1^2*n2


################
### 7 Variables:

polyRes4(7)

k1^7 + n1^7 + n2^7 + k2^7 - 7*k1*n1^5*n2 - 7*k2*n1*n2^5 +
	+ 7*k1^4*k2^2*n1 + 7*k1^2*k2^4*n2 + 7*k1^4*k2*n2^2 + 7*k1*k2^4*n1^2 + 7*k1^2*k2*n1^4 + 7*k1*k2^2*n2^4 +
	- 7*k1^3*n1*n2^3 - 7*k2^3*n1^3*n2 + 14*k1^2*n1^3*n2^2 + 14*k2^2*n1^2*n2^3 +
	- 21*k1^3*k2*n1^2*n2 - 21*k1*k2^3*n1*n2^2

k1^7 + n1^7 + n2^7 + k2^7 +
	+ 7*k1^2*k2*n1^4 - 21*k1^3*k2*n1^2*n2 + 7*k1^4*k2*n2^2 +
	+ 7*k1*k2^2*n2^4 - 21*k1*k2^3*n1*n2^2 + 7*k1*k2^4*n1^2 +
	- 7*k1*n1^5*n2 + 14*k1^2*n1^3*n2^2 - 7*k1^3*n1*n2^3 +
	- 7*k2*n1*n2^5 + 14*k2^2*n1^2*n2^3 - 7*k2^3*n1^3*n2 +
	# ???
	+ 7*k1^4*k2^2*n1 + 7*k1^2*k2^4*n2;


# Class 2 Poly:
p = expand.u.pm(toPoly.pm("k1 + n1*u + n2*u^2 + k2*u^3"), n=7)
pR = mult.pm(p, toPoly.pm("k1+k2+n1+n2"))
# print.pm(pR)
# p0 = polyClip();
# diff.pm(p0, pR); # SUCCESS !

x = 2*cos(2*(1:2)*pi/5) + 2;
7*x^2 - 21*x + 7 # = 0
x = 2 - 2*cos(2*(1:3)*pi/7);
x^3 - 7*x^2 + 14*x - 7 # = 0

# [old]
x = 2*cos(pi*c(seq(1,14, by=2))/7) + 1;
x^7 - 7*x^6 + 14*x^5 - 21*x^3 + 7*x^2 + 7*x + 1 # = 0


### Special Cases:

### Case: k2 = k1; n2 = n1;
2*(n1 + k1) * (n1^3 - 4*n1^2*k1 + 3*n1*k1^2 + k1^3)^2

# x = n1/k1;
x = - 2*cos(2*seq(1,3)*pi/7) + 1;
x^3 - 4*x^2 + 3*x + 1 # = 0


### Case: k2 = -k1; n2 = - n1;
# Note: (k1 + k2) == 0;
7*(k1 + k2) * (n1^3 + 2*n1^2*k1 - n1*k1^2 - k1^3)^2

x = - 2*cos(2*seq(1,3)*pi/7) - 1;
x^3 + 2*x^2 - x - 1 # = 0


### Case: n2 = 0
k1^7 + k2^7 + n1^7 + 7*k1^4*k2^2*n1 + 7*k1^2*k2*n1^4 + 7*k1*k2^4*n1^2

# Class 2 Poly:
p = expand.u.pm(toPoly.pm("k1 + n1*u + k2*u^3"), n=7)
pR = mult.pm(p, toPoly.pm("k1+k2+n1"))
print.pm(pR)

# same as above!
k1^7 + k2^7 + n1^7 + 7*k1^4*k2^2*n1 + 7*k1^2*k2*n1^4 + 7*k1*k2^4*n1^2


# n1 = - k1; n2 = 0;
k2*(k2^6 + 7*k2^3*k1^3 - 7*k2*k1^5 + 7*k1^6)


u = unity(7, all=T); u = u[-1];
x = u^3 - u;
x^6 + 7*x^3 - 7*x + 7 # = 0


################
### 9 Variables:

polyRes4(9)

k1^9 + n1^9 + n2^9 + k2^9 + 3*k1^6*k2^3 + 3*k1^3*k2^6 +
	- 9*k1^5*k2*n2^3 - 9*k1*k2^5*n1^3 + 18*k1^4*k2^2*n1^3 + 18*k1^2*k2^4*n2^3 +
		+ 9*k1^2*k2*n1^6 + 9*k1*k2^2*n2^6 +
	+ 9*k1^4*n1*n2^4 + 9*k2^4*n1^4*n2 - 30*k1^3*n1^3*n2^3 - 30*k2^3*n1^3*n2^3 +
		+ 27*k1^2*n1^5*n2^2 + 27*k2^2*n1^2*n2^5 - 9*k1*n1^7*n2 - 9*k2*n1*n2^7 +
	- 27*k1^5*k2^2*n1*n2 - 27*k1^2*k2^5*n1*n2 + 54*k1^4*k2*n1^2*n2^2 + 54*k1*k2^4*n1^2*n2^2 +
	- 45*k1^3*k2*n1^4*n2 - 45*k1*k2^3*n1*n2^4

k1^9 + n1^9 + n2^9 + k2^9 +
	# x = 2*cos(2*(1:3)*pi/7) + 2;
	+ 9*k1^2*k2*n1^6 - 45*k1^3*k2*n1^4*n2 + 54*k1^4*k2*n1^2*n2^2 - 9*k1^5*k2*n2^3 +
	+ 9*k1*k2^2*n2^6 - 45*k1*k2^3*n1*n2^4 + 54*k1*k2^4*n1^2*n2^2 - 9*k1*k2^5*n1^3 +
	# x = 2 - 2*cos(2*(1:4)*pi/9);
	- 9*k1*n1^7*n2 + 27*k1^2*n1^5*n2^2 - 30*k1^3*n1^3*n2^3 + 9*k1^4*n1*n2^4 +
	- 9*k2*n1*n2^7 + 27*k2^2*n1^2*n2^5 - 30*k2^3*n1^3*n2^3 + 9*k2^4*n1^4*n2 +
	# ...
	- 27*k1^5*k2^2*n1*n2 + 18*k1^4*k2^2*n1^3 +
	- 27*k1^2*k2^5*n1*n2 + 18*k1^2*k2^4*n2^3 +
	+ 3*k1^6*k2^3 + 3*k1^3*k2^6;


# Class 2 Poly:
p = expand.u.pm(toPoly.pm("k1 + n1*u + n2*u^2 + k2*u^3"), n=9)
# problem with non-prime powers;
p = replace.pm(p, toPoly.pm("-1 - u^3"), "u", pow=6)
pR = mult.pm(p, toPoly.pm("k1+k2+n1+n2"))
# print.pm(pR)
# p0 = polyClip();
# diff.pm(p0, pR); # SUCCESS !


# [old]
x = 2*cos(2*(1:3)*pi/7) + 2;
9*x^3 - 45*x^2 + 54*x - 9 # = 0
x = 2 - 2*cos(2*(1:4)*pi/9);
x^4 - 9*x^3 + 27*x^2 - 30*x + 9 # = 0
x = 2*cos(pi*c(seq(1,18, by=2))/18) + 1;
x^9 - 9*x^8 + 27*x^7 - 21*x^6 - 36*x^5 + 54*x^4 + 9*x^3 - 27*x^2 + 2 # = 0
# - possible the k1^j * n1^(9 - 2*j) * n2^j
#   and the equivalent k2^j * n2^(9 - 2*j) * n1^j parts;
x = 2*cos(pi*c(seq(1,18, by=2))/9) + 1;
x^9 - 9*x^8 + 27*x^7 - 21*x^6 - 36*x^5 + 54*x^4 + 9*x^3 - 27*x^2 + 0*x + 4 # = 0
x = 2*cos(pi*c(seq(1,9, by=2))/9) + 2;
x^5 - 9*x^4 + 27*x^3 - 30*x^2 + 9*x # = 0


### Special Cases:

### Case: k1 = k2; n1 = n2;
2*(k1 + n1) * (2*k1 - n1)^2 * (n1^3 - 3*n1^2*k1 + k1^3)^2

# x = n1/k1;
x = 1 - 2*cos(2*seq(1,4)*pi/9);
(x - 2)*(x^3 - 3*x^2 + 1) # = 0


### Case: k2 = -k1; n2 = - n1;
# Note: (k1 + k2) == 0;
9*n1^2 * (k1 + k2) * (n1^6 + 6*n1^5*k1 + 9*n1^4*k1^2 - 6*n1^3*k1^3 - 18*n1^2*k1^4 + 9*k1^6)
9*n1^2 * (k1 + k2) * (n1^3 + 3*n1^2*k1 - 3*k1^3)^2

x = 2*cos(c(1,3,5,7)*pi/9) - 1;
x*(x^3 + 3*x^2 - 3) # = 0


### Case: n2 = 0;
k1^9 + k2^9 + 3*k1^3*k2^6 - 9*k1*k2^5*n1^3 +
	+ n1^9 + 9*k1^2*k2*n1^6 + 18*k1^4*k2^2*n1^3 + 3*k1^6*k2^3
# n1 = n2 = 0;
k1^9 + k2^9 + 3*k1^6*k2^3 + 3*k1^3*k2^6
(k1^3 + k2^3)^3

# n1 = - k1; n2 = 0;
k2*(k2^8 + 3*k2^5*k1^3 + 9*k2^4*k1^4 + 3*k2^2*k1^6 - 18*k2*k1^7 + 9*k1^8)
k2*(k2^2 + 3*k2*k1 + 3*k1^2)*(k2^6 - 3*k1*k2^5 + 6*k1^2*k2^4 - 6*k1^3*k2^3 + 9*k1^4*k2^2 - 9*k1^5*k2 + 3*k1^6)


u = unity(9, all=T); u = u[-1];
x = u - u^6;
x^8 + 3*x^5 + 9*x^4 + 3*x^2 - 18*x + 9 # = 0
#
(x^2 + 3*x + 3)*(x^6 - 3*x^5 + 6*x^4 - 6*x^3 + 9*x^2 - 9*x + 3)


#################
### 11 Variables:

polyRes4(11)


k1^11 + n1^11 + n2^11 + k2^11 +
	+ 11*k1^2*k2*n1^8 - 77*k1^3*k2*n1^6*n2 + 165*k1^4*k2*n1^4*n2^2 - 110*k1^5*k2*n1^2*n2^3 + 11*k1^6*k2*n2^4 +
	+ 11*k1*k2^2*n2^8 - 77*k1*k2^3*n1*n2^6 + 165*k1*k2^4*n1^2*n2^4 - 110*k1*k2^5*n1^3*n2^2 + 11*k1*k2^6*n1^4 +
	#
	- 11*k1*n1^9*n2 + 44*k1^2*n1^7*n2^2 - 77*k1^3*n1^5*n2^3 + 55*k1^4*n1^3*n2^4 - 11*k1^5*n1*n2^5 +
	- 11*n1*k2*n2^9 + 44*k2^2*n1^2*n2^7 - 77*k2^3*n1^3*n2^5 + 55*k2^4*n1^4*n2^3 - 11*k2^5*n1^5*n2 +
	# ???
 - 11*k1^7*k2^3*n2 + 66*k1^6*k2^2*n1*n2^2 + 22*k1^3*k2^6*n2^2 +
 - 11*k1^3*k2^7*n1 + 66*k1^2*k2^6*n1^2*n2 + 22*k1^6*k2^3*n1^2 +
 - 110*k1^5*k2^2*n1^3*n2 - 110*k1^2*k2^5*n1*n2^3 +
 + 33*k1^4*k2^2*n1^5 + 33*k1^2*k2^4*n2^5


# Class 2 Poly:
p = expand.u.pm(toPoly.pm("k1 + n1*u + n2*u^2 + k2*u^3"), n=11)
pR = mult.pm(p, toPoly.pm("k1+k2+n1+n2"))
# print.pm(pR)
# p0 = polyClip();
# diff.pm(p0, pR); # SUCCESS !


# [old]
x = 2*cos(2*(1:4)*pi/9) + 2;
11*x^4 - 77*x^3 + 165*x^2 - 110*x + 11 # = 0
x = 2 - 2*cos(2*(1:5)*pi/11); # k1^j*n1^(11-2*j)*n2^j; j = 0:5;
x^5 - 11*x^4 + 44*x^3 - 77*x^2 + 55*x - 11 # = 0


### Special Cases:

### Case: k1 = k2; n1 = n2;
2*(k1 + n1) * (n1^5 - 6*n1^4*k1 + 10*n1^3*k1^2 - n1^2*k1^3 - 6*n1*k1^4 + k1^5)^2

# x = n1/k1;
x = 1 - 2*cos(2*seq(1,5)*pi/11);
x^5 - 6*x^4 + 10*x^3 - x^2 - 6*x + 1 # = 0


### Case: k2 = -k1; n2 = - n1;
# Note: (k1 + k2) == 0;
11*(k1 + k2) * (n1^5 + 4*n1^4*k1 + 2*n1^3*k1^2 - 5*n1^2*k1^3 - 2*n1*k1^4 + k1^5)^2

x = - 2*cos(2*seq(1,5)*pi/11) - 1;
x^5 + 4*x^4 + 2*x^3 - 5*x^2 - 2*x + 1 # = 0


### Case: n2 = 0
# n1 = - k1; n2 = 0;
k2*(k2^10 + 11*k2^6*k1^4 + 11*k2^5*k1^5 + 22*k2^2*k1^8 - 33*k2*k1^9 + 11*k1^10)

u = unity(11, all=T); u = u[-1];
x = u^8 - u;
x^10 + 11*x^6 + 11*x^5 + 22*x^2 - 33*x + 11 # = 0


#################
### 13 Variables:

polyRes4(13)


k1^13 + n1^13 + n2^13 + k2^13 +
	+ 13*k1^2*k2*n1^10 - 117*k1^3*k2*n1^8*n2 + 364*k1^4*k2*n1^6*n2^2 - 455*k1^5*k2*n1^4*n2^3 +
		+ 195*k1^6*k2*n1^2*n2^4 - 13*k1^7*k2*n2^5 +
	+ 13*k1*k2^2*n2^10 - 117*k1*k2^3*n1*n2^8 + 364*k1*k2^4*n1^2*n2^6 - 455*k1*k2^5*n1^3*n2^4 +
		+ 195*k1*k2^6*n1^4*n2^2 - 13*k1*k2^7*n1^5 +
	- 13*k1*n1^11*n2 + 65*k1^2*n1^9*n2^2 - 156*k1^3*n1^7*n2^3 + 182*k1^4*n1^5*n2^4 +
		- 91*k1^5*n1^3*n2^5 + 13*k1^6*n1*n2^6 +
	- 13*k2*n1*n2^11 + 65*k2^2*n1^2*n2^9 - 156*k2^3*n1^3*n2^7 + 182*k2^4*n1^4*n2^5 +
		- 91*k2^5*n1^5*n2^3 + 13*k2^6*n1^6*n2 +
	# ???
 + 13*k1^8*k2^4*n1 + 26*k1^8*k2^3*n2^2 +
 + 13*k1^4*k2^8*n2 + 26*k1^3*k2^8*n1^2 +
 - 130*k1^7*n1^2*k2^3*n2 - 130*k1^7*n1*k2^2*n2^3 + 52*k1^4*n1^7*k2^2 +
 - 130*k1^3*n1*k2^7*n2^2 - 130*k1^2*n1^3*k2^7*n2 + 52*k1^2*k2^4*n2^7 +
 + 390*k1^6*n1^3*k2^2*n2^2 + 390*k1^2*n1^2*k2^6*n2^3 +
 + 65*k1^6*n1^4*k2^3 + 65*k1^3*k2^6*n2^4 - 273*k1^5*n1^5*k2^2*n2 - 273*k1^2*n1*k2^5*n2^5;


# Class 2 Poly:
p = expand.u.pm(toPoly.pm("k1 + n1*u + n2*u^2 + k2*u^3"), n=13)
pR = mult.pm(p, toPoly.pm("k1+k2+n1+n2"))
# print.pm(pR)
# p0 = polyClip();
# diff.pm(p0, pR); # SUCCESS !


# [old]
x = 2*cos(2*(1:5)*pi/11) + 2;
13*x^5 - 117*x^4 + 364*x^3 - 455*x^2 + 195*x - 13 # = 0
x = 2 - 2*cos(2*(1:6)*pi/13); # k1^j*n1^(11-2*j)*n2^j; j = 0:6;
x^6 - 13*x^5 + 65*x^4 - 156*x^3 + 182*x^2 - 91*x + 13 # = 0


### Special Cases:

### Case: k1 = k2; n1 = n2;
2*(k1 + n1) * (n1^6 - 7*n1^5*k1 + 15*n1^4*k1^2 - 6*n1^3*k1^3 - 11*n1^2*k1^4 + 6*n1*k1^5 + k1^6)^2


x = - 2*cos(2*seq(1,6)*pi/13) + 1;
x^6 - 7*x^5 + 15*x^4 - 6*x^3 - 11*x^2 + 6*x + 1 # = 0


### Case: k2 = -k1; n2 = - n1;
# Note: (k1 + k2) == 0;
13*(k1 + n1) * (n1^6 - 5*n1^5*k1 + 5*n1^4*k1^2 + 6*n1^3*k1^3 - 7*n1^2*k1^4 - 2*n1*k1^5 + k1^6)^2

x = 2*cos(2*seq(1,6)*pi/13) + 1;
x^6 - 5*x^5 + 5*x^4 + 6*x^3 - 7*x^2 - 2*x + 1 # = 0


### Case: n2 = 0
# n1 = - k1; n2 = 0;
k2*(k2^12 + 26*k1^5*k2^7 + 13*k1^6*k2^6 - 13*k1^9*k2^3 + 65*k1^10*k2^2 - 52*k1^11*k2 + 13*k1^12)

u = unity(13, all=T); u = u[-1];
x = u - u^8;
x^12 + 26*x^7 + 13*x^6 - 13*x^3 + 65*x^2 - 52*x + 13 # = 0


################
################

### Note:
# - All formulas can be computed using Class 2 polynomials;

# [old]
# - formula seems dependent on: Order (mod 6);
### Order: 6*v
(k1^(2*v) - k2^(2*v))^3 + ...
### Order: 6*v + 2
(k1^(3*v+1) - n2^(3*v+1))^2 - (k2^(3*v+1) - n1^(3*v+1))^2 + ...
### Order: 6*v + 4
(k1^(3*v+2) + n2^(3*v+2))^2 - (k2^(3*v+2) + n1^(3*v+2))^2 + ...


################
### 6 Variables:

polyRes4(6)

k1^6 - n1^6 + n2^6 - k2^6 - 3*k1^4*k2^2 + 3*k1^2*k2^4 + 2*k1^3*n2^3 - 2*k2^3*n1^3 +
	- 6*k1^2*k2*n1^3 + 6*k1*k2^2*n2^3 + 6*k1*n1^4*n2 - 6*k2*n1*n2^4 +
	- 9*k1^2*n1^2*n2^2 + 9*k2^2*n1^2*n2^2 + 12*k1^3*k2*n1*n2 - 12*k1*k2^3*n1*n2

k1^6 - n1^6 + n2^6 - k2^6 +
	- 3*k1^4*k2^2 + 3*k1^2*k2^4 +
	+ 2*k1^3*n2^3 - 9*k1^2*n1^2*n2^2 + 6*k1*n1^4*n2 +
	- 2*k2^3*n1^3 + 9*k2^2*n1^2*n2^2 - 6*k2*n1*n2^4 +
	+ 12*k1^3*k2*n1*n2 - 6*k1^2*k2*n1^3 +
	- 12*k1*k2^3*n1*n2 + 6*k1*k2^2*n2^3;

# Class 2 Poly:
p = expand.u.pm(toPoly.pm("k1 + n1*u + n2*u^2 + k2*u^3"), n=6)
pR = mult.pm(p, toPoly.pm("k1+k2+n1+n2"))
# print.pm(pR)
# p0 = polyClip();
diff.pm(p0, pR); # SUCCESS !


# [old]
x = 2*cos(seq(1,4, by=2)*pi/4) + 2;
3*x^2 - 12*x + 6 # = 0
x = 2 - 2*cos(seq(1,6,by=2)*pi/6); # k1^j*n1^(6-2*j)*n2^j; j = 0:3;
x^3 - 6*x^2 + 9*x - 2 # = 0
x = 2*cos(pi*c(seq(1,12, by=2))/6) + 1;
x^6 - 6*x^5 + 9*x^4 + 4*x^3 - 12*x^2 + 4 # = 0
x = 2*cos(pi*c(seq(1,6, by=2))/6) + 1;
x^3 - 3*x^2 + 2 # = 0


################
### 8 Variables:

polyRes4(8)

k1^8 - n1^8 + n2^8 - k2^8 - 2*k1^4*n2^4 + 2*k2^4*n1^4 +
 + 8*k1^5*k2^2*n2 - 8*k1^2*k2^5*n1 - 12*k1^4*n1^2*k2^2 + 12*k1^2*k2^4*n2^2 +
 + 16*k1^3*n1^2*n2^3 - 16*k2^3*n1^3*n2^2 +
 - 8*k1^2*k2*n1^5 + 8*k1*k2^2*n2^5 +
 + 8*k1*n1^6*n2 - 8*k2*n1*n2^6 +
 - 24*k1^4*k2*n1*n2^2 + 24*k1*k2^4*n1^2*n2 - 20*k1^2*n1^4*n2^2 + 20*n1^2*k2^2*n2^4 +
 + 32*k1^3*k2*n1^3*n2 - 32*k1*k2^3*n1*n2^3

# Class 2 Poly:
p = expand.u.pm(toPoly.pm("k1 + n1*u + n2*u^2 + k2*u^3"), n=8)
pR = mult.pm(p, toPoly.pm("k1+k2+n1+n2"))
# print.pm(pR)
# p0 = polyClip();
diff.pm(p0, pR); # SUCCESS !


# [old]
x = 2*cos(pi*c(seq(1,16, by=2))/8) + 1;
x^8 - 8*x^7 + 20*x^6 - 8*x^5 - 30*x^4 + 24*x^3 + 12*x^2 - 8*x + 1 # = 0
x = 2*cos(pi*c(seq(1,8, by=2))/8) + 1;
x^4 - 4*x^3 + 2*x^2 + 4*x - 1 # = 0
x = 2*cos(pi*c(seq(1,4, by=2))/4) + 1;
x^2 - 2*x - 1 # = 0 => k1^8 - 2*k1^4*n2^4 - n2^4;


#################
### 10 Variables:

polyRes4(10)

k1^10 - n1^10 + n2^10 - k2^10 + 2*k1^5*n2^5 - 2*k2^5*n1^5 +
 + 10*k1*n1^8*n2 - 10*n1*k2*n2^8 - 10*k1^2*n1^7*k2 + 10*k1*k2^2*n2^7 +
 - 10*k1^6*n1*k2^3 - 15*k1^6*k2^2*n2^2 + 10*k1^3*k2^6*n2 - 35*k1^2*n1^6*n2^2 + 15*k1^2*n1^2*k2^6 +
 + 35*n1^2*k2^2*n2^6 + 60*k1^5*n1^2*k2^2*n2 + 40*k1^5*n1*k2*n2^3 + 60*k1^3*n1^5*k2*n2 +
 - 40*k1*n1^3*k2^5*n2 - 25*k1^4*k2^2*n1^4 + 25*k1^2*k2^4*n2^4 +
 - 25*k1^4*n1^2*n2^4 + 25*k2^4*n1^4*n2^2 + 50*k1^3*n1^4*n2^3 - 50*n1^3*k2^3*n2^4 +
 - 100*k1^4*k2*n1^3*n2^2 + 100*k1*k2^4*n1^2*n2^3 - 60*k1^2*n1*k2^5*n2^2 - 60*k1*n1*k2^3*n2^5


# Class 2 Poly:
p = expand.u.pm(toPoly.pm("k1 + n1*u + n2*u^2 + k2*u^3"), n=10)
pR = mult.pm(p, toPoly.pm("k1+k2+n1+n2"))
# print.pm(pR)
# p0 = polyClip();
diff.pm(p0, pR); # SUCCESS !


# [old]
x = 2*cos(pi*c(seq(1,20, by=2))/10) + 1;
x^10 - 10*x^9 + 35*x^8 - 40*x^7 - 35*x^6 + 98*x^5 - 15*x^4 - 60*x^3 + 15*x^2 + 10*x + 1 # = 0
x = 2*cos(pi*c(seq(2,20, by=2))/10) + 1;
x^10 - 10*x^9 + 35*x^8 - 40*x^7 - 35*x^6 + 98*x^5 - 15*x^4 - 60*x^3 + 15*x^2 + 10*x - 3 # = 0


#################
### 12 Variables:

polyRes4(12)


k1^12 - n1^12 + n2^12 - k2^12 - 3*k1^8*k2^4 + 3*k1^4*k2^8 +
 + 12*k1*n1^10*n2 - 12*n1*k2*n2^10 - 12*k1^2*n1^9*k2 + 12*k1*k2^2*n2^9 +
 - 54*k1^2*n1^8*n2^2 + 54*n1^2*k2^2*n2^8 + 48*k1^7*n1*k2^3*n2 +
 + 24*k1^7*k2^2*n2^3 + 96*k1^3*n1^7*k2*n2 - 48*k1^3*n1*k2^7*n2 - 24*k1^2*n1^3*k2^7 - 96*k1*n1*k2^3*n2^7 +
 - 40*k1^6*n1^3*k2^3 - 180*k1^6*n1^2*k2^2*n2^2 - 60*k1^6*n1*k2*n2^4 - 2*k1^6*n2^6 - 42*k1^4*n1^6*k2^2 +
 + 112*k1^3*n1^6*n2^3 + 40*k1^3*k2^6*n2^3 + 180*k1^2*n1^2*k2^6*n2^2 + 42*k1^2*k2^4*n2^6 + 60*k1*n1^4*k2^6*n2 +
 + 2*n1^6*k2^6 - 112*n1^3*k2^3*n2^6 + 180*k1^5*n1^4*k2^2*n2 + 240*k1^5*n1^3*k2*n2^3 +
 + 36*k1^5*n1^2*n2^5 - 252*k1^4*n1^5*k2*n2^2 - 180*k1^2*n1*k2^5*n2^4 - 240*k1*n1^3*k2^5*n2^3 + 252*k1*n1^2*k2^4*n2^5 +
 - 36*n1^5*k2^5*n2^2 - 105*k1^4*n1^4*n2^4 + 105*n1^4*k2^4*n2^4


# Class 2 Poly:
p = expand.u.pm(toPoly.pm("k1 + n1*u + n2*u^2 + k2*u^3"), n=12)
p = replace.pm(p, toPoly.pm("u^2 - 1"), "u", pow=4);
pR = mult.pm(p, toPoly.pm("k1+k2+n1+n2"))
# print.pm(pR)
# p0 = polyClip();
diff.pm(p0, pR); # SUCCESS !


