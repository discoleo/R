########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Heterogeneous Symmetric:
### Resonances
###
### draft v.0.1a


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


### Resonances
resonance = function(p, n=3) {
	# "Resonance" with Roots of unity;
	# currently only for 2-variable terms: x^p1*y^p2;
	sg = if(n %% 2 == 1) 1 else -1;
	p.all = p[1]^n + sg*p[2]^n;
	r = list(p=p.all, f=factors(p.all), p.trivial = sum(p));
	return(r);
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


######################
######################

######################
### Roots of Unity ###
### Generalization ###
######################

### Sys[i]:

### System with [i] Variables
# - for i = all odd;
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
# new solution:
# (x1,x2,x3,x4,x5) * (m^1, m^9, m^4, m^3, m^5);


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

## Order (a,b,c)
# a^2*c^4 - 2*a^4*c^2 + 4*b^2*a^3*c - b^4*a^2 + a^6


### Order: 2+1+1
### 4 Variables
# i = 4; p = c(2,1,1)
# p = Divisors(2^6);
### Trivial:
p = 4;
### Non-Trivial & Combinations:
p = c();
# - there are various quasi-non-trivial solutions;
### Ex: p = 8 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^3, m^7, m^3, m^7);
# (3,7,3), (7,3,7), (3,7,3), (7,3,7)


### Order: 3+1+1
### 4 Variables
# i = 4; p = c(3,1,1)
# p = Divisors(3^3*5^2);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(15);
# - the computed ones do NOT work;
# - but there are quasi-non-trivial solutions for p = 5;
### Ex: p = 5 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^1, m^2, m^0, m^4);
# (1,2,0), (2,0,4), (0,4,1), (4,1,2)
# (3,2,4), (2,4,0), (4,0,3), (0,3,2)
### Ex: p = 15 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^1, m^11, m^1, m^11);
# (1,11,1), (11,1,11), (1,11,1), (11,1,11)


### Order: 4+1+1
### 4 Variables
# i = 4; p = c(4,1,1)
# p = Divisors(2^8*3*5);
### Trivial:
p = 6;
### Non-Trivial & Combinations:
p = c();
# - some of the computed ones do NOT work;
# - but there are various quasi-non-trivial solutions;
### Ex: p = 8 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^5, m^7, m^5, m^7);
# (5,7,5), (7,5,7), (5,7,5), (7,5,7)


### Order: (p1,1,1)

### 5 Variables:
# p = Divisors(p1^3*(p1-1)*(2*p1^2 + 3))
# - but condition is often insufficient;

### Order: 2+1+1
### 5 Variables
# i = 5; p = c(2,1,1)
# p = Divisors(88);
### Trivial:
p = 4;
### Non-Trivial & Combinations:
p = c(11, 22); # & higher
# ex: p = 11 =>
# new solution:
# (x1,x2,...,x5) * (m^1, m^4, m^5, m^9, m^3);
# (1,4,5), (4,5,9), (5,9,3), (9,3,1), (3,1,4)
# ex: p = 22 =>
# new solution:
# (x1,x2,...,x5) * (m^1, m^15, m^5, m^9, m^3);
# (1,15,5), (15,5,9), (5,9,3), (9,3,1), (3,1,15)


### Order: 3+1+1
### 6 Variables
# i = 6; p = c(3,1,1)
# p = Divisors(???);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(9); # possible others
# ex: p = 9 =>
# new solution:
# (x1,x2,...,x6) * (m^1, m^2, m^4, m^8, m^7, m^5);
# (1,2,4), (2,4,8), (4,8,7), (8,7,5), (7,5,1), (5,1,2)

