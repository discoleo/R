########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Heterogeneous Symmetric:
### Resonances
###
### draft v.0.1h


### Resonances in Polynomial Systems
# System: 3..n Variables
# Type: Heterogeneous Symmetric


# Note:
# - moved to this new file from file:
#   Poly.System.Hetero.Symmetric.S3.Mixed.R;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")
source("Polynomials.Helper.Matrix.R")


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

### Generators
# n = number of variables;
polyRes3 = function(n) {
	p = det.mpm(diag.lpm(c("k", "n", "p"), n=n));
	p = sort.pm(p, c("k", "n"));
	if( ! inherits(p, "pm")) class(p) = c("pm", class(p));
	return(p);
}

### Test

# n = number of variables;
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

