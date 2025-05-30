
## Polynomials: Class 3
## Roots: All Real


### Various Examples

# Based on:
# => cos(2*pi / (2*n+1));
# => subset of 2 * cos(2*pi / higher);
# => cos(pi / (n+1));


####################

### Helper Functions

source("Polynomials.Helper.R")


####################
####################

###############
### Order 4 ###

### Based on cos(2*pi/9)
r = 2*cos(2*pi*seq(4) / 9)
# Note: 1 trivial integer root;

x = r; poly.calc0(x, digits = 5)
x^4 + x^3 - 3*x^2 - 2*x + 1

x = r^2 + 2*r; poly.calc0(x, digits = 5)
x^4 - 5*x^3 - 3*x^2 + 4*x + 1


### Based on cos(2*pi/15)
r = 2*cos(2*pi*c(1,2,4,7)/15)

x = r; poly.calc0(x, digits = 5)
x^4 - x^3 - 4*x^2 + 4*x + 1

x = r^2 - r; poly.calc0(x, digits = 5)
x^4 - 8*x^3 + 14*x^2 - 7*x + 1

x = r^3 - r - 1; poly.calc0(x, digits = 5)
x^4 + 4*x^3 - 19*x^2 - 16*x + 1


### Based on cos(pi/5)
r = 2*cos(pi*seq(4)/5)

# Trivial Poly:
x = r; poly.calc0(x, digits = 5)
x^4 - 3*x^2 + 1

# Trivial Roots:
x = r^2 + r; poly.calc0(x, digits = 5)
x^4 - 6*x^3 + 8*x^2 - 2*x - 1

x = r*3 + r^2 - 2*r; poly.calc0(x, digits = 5)
x^4 - 6*x^3 + 8*x^2 - 2*x - 1

# Less-Trivial:
x = r^2 + 2*r + 1; poly.calc0(x, digits = 5)
x^4 - 10*x^3 + 23*x^2 - 10*x + 1
(x^2 - 5*x + 1)^2 - 4*x^2


###################
###################

###############
### Order 5 ###
###############

### Based on cos(2*pi/11)
r = 2*cos(2*pi*seq(5)/11)

#
x = r; poly.calc0(x, digits = 5)
x^5 + x^4 - 4*x^3 - 3*x^2 + 3*x + 1

x = r^2 - r; poly.calc0(x, digits = 5)
x^5 - 10*x^4 + 29*x^3 - 25*x^2 + 3*x + 1

x = r^2 + r; poly.calc0(x, digits = 5)
x^5 - 8*x^4 + 19*x^3 - 15*x^2 + x + 1

x = r^2 + 2*r; poly.calc0(x, digits = 5)
x^5 - 7*x^4 + 2*x^3 + 17*x^2 + 9*x + 1

x = r^4 - r^2 + 5*r; poly.calc0(x, digits = 5)
x^5 - 11*x^4 - 55*x^3 + 2*121*x^2 + 3*121*x - 121


### Gen:

###
s = c(1,-1,2,-2,-1);
x = sapply(seq(5), \(id) {
	sum(s * (2*cos(2*id*pi/11))^(4:0));
})
# Note: NOT very practical;
s4 = s[1]; s3 = s[2]; s2 = s[3]; s1 = s[4]; s0 = s[5];
poly.calc0(x);
x^5 - 44*x^4 + 429*x^3 - 638*x^2 - 2156*x - 473 # = 0
x^5 - (25*s4 - 4*s3 + 9*s2 - 1*s1 + 5*s0)*x^4 +
	+ (184*s4^2 - 36*s4*s3 - 31*s3^2 + 147*s4*s2 - 20*s3*s2 + 28*s2^2 - 9*s4*s1 - 21*s3*s1 +
		- 5*s2*s1 - 4*s1^2 + 100*s4*s0 - 16*s3*s0 + 36*s2*s0 - 4*s1*s0 + 10*s0^2)*x^3 +
	+ ...; # see below for free term;

# TODO: ugly;


### Derivation:

p = lapply(seq(5), \(id) {
	r = 2*cos(2*pi*id / 11);
	s = diag(5);
	colnames(s) = paste0("s", 4:0);
	p = as.data.frame(s);
	p$coeff = r^(4:0);
	as.pm(p);
})
print.pm(sum.lpm(p))
print.pm(sum.lpm(list(
	prod.pm(p[[1]], sum.lpm(p[ -1])),
	prod.pm(p[[2]], sum.lpm(p[3:5])),
	prod.pm(p[[3]], sum.lpm(p[4:5])),
	prod.pm(p[[4]], p[[5]]) )))

print.pm(prod.pm(p), digits = 3, width = 80) # Monster:

s4^5 - 3*s4^4*s3 - 3*s4^3*s3^2 + 4*s4^2*s3^3 + s4*s3^4 - s3^5 + 15*s4^4*s2 +
 - 3*s4^3*s3*s2 - 16*s4^2*s3^2*s2 + 2*s4*s3^3*s2 + 3*s3^4*s2 + 35*s4^3*s2^2 +
 - 8*s4^2*s3*s2^2 - 15*s4*s3^2*s2^2 + 3*s3^3*s2^2 + 28*s4^2*s2^3 - 5*s4*s3*s2^3 +
 - 4*s3^2*s2^3 + 9*s4*s2^4 - s3*s2^4 + s2^5 - 42*s4^4*s1 - 29*s4^3*s3*s1 + 58*s4^2*s3^2*s1 +
 + 12*s4*s3^3*s1 - 15*s3^4*s1 - 97*s4^3*s2*s1 - 30*s4^2*s3*s2*s1 +
 + 51*s4*s3^2*s2*s1 + 3*s3^3*s2*s1 - 79*s4^2*s2^2*s1 - 13*s4*s3*s2^2*s1 + 16*s3^2*s2^2*s1 +
 - 26*s4*s2^3*s1 - 2*s3*s2^3*s1 - 3*s2^4*s1 - 60*s4^3*s1^2 + 86*s4^2*s3*s1^2 +
 + 32*s4*s3^2*s1^2 - 35*s3^3*s1^2 - 67*s4^2*s2*s1^2 + 65*s4*s3*s2*s1^2 + 8*s3^2*s2*s1^2 +
 - 25*s4*s2^2*s1^2 + 15*s3*s2^2*s1^2 - 3*s2^3*s1^2 + 31*s4^2*s1^3 + 20*s4*s3*s1^3 +
 - 28*s3^2*s1^3 + 21*s4*s2*s1^3 + 5*s3*s2*s1^3 + 4*s2^2*s1^3 + 4*s4*s1^4 - 9*s3*s1^4 +
 + s2*s1^4 - s1^5 + 155*s4^4*s0 + 68*s4^3*s3*s0 - 180*s4^2*s3^2*s0 +
 - 27*s4*s3^3*s0 + 42*s3^4*s0 + 441*s4^3*s2*s0 - 26*s4^2*s3*s2*s0 - 200*s4*s3^2*s2*s0 +
 + 29*s3^3*s2*s0 + 384*s4^2*s2^2*s0 - 45*s4*s3*s2^2*s0 - 58*s3^2*s2^2*s0 + 130*s4*s2^3*s0 +
 - 12*s3*s2^3*s0 + 15*s2^4*s0 + 94*s4^3*s1*s0 - 255*s4^2*s3*s1*s0 +
 - 56*s4*s3^2*s1*s0 + 97*s3^3*s1*s0 + 43*s4^2*s2*s1*s0 - 221*s4*s3*s2*s1*s0 + 30*s3^2*s2*s1*s0 +
 - 3*s4*s2^2*s1*s0 - 51*s3*s2^2*s1*s0 - 3*s2^3*s1*s0 - 113*s4^2*s1^2*s0 +
 - 25*s4*s3*s1^2*s0 + 79*s3^2*s1^2*s0 - 84*s4*s2*s1^2*s0 + 13*s3*s2*s1^2*s0 +
 - 16*s2^2*s1^2*s0 - 3*s4*s1^3*s0 + 26*s3*s1^3*s0 + 2*s2*s1^3*s0 + 3*s1^4*s0 + 403*s4^3*s0^2 +
 - 149*s4^2*s3*s0^2 - 157*s4*s3^2*s0^2 + 60*s3^3*s0^2 + 580*s4^2*s2*s0^2 +
 - 124*s4*s3*s2*s0^2 - 86*s3^2*s2*s0^2 + 255*s4*s2^2*s0^2 - 32*s3*s2^2*s0^2 + 35*s2^3*s0^2 +
 - 40*s4^2*s1*s0^2 - 139*s4*s3*s1*s0^2 + 67*s3^2*s1*s0^2 - 31*s4*s2*s1*s0^2 +
 - 65*s3*s2*s1*s0^2 - 8*s2^2*s1*s0^2 - 38*s4*s1^2*s0^2 + 25*s3*s1^2*s0^2 +
 - 15*s2*s1^2*s0^2 + 3*s1^3*s0^2 + 184*s4^2*s0^3 - 36*s4*s3*s0^3 - 31*s3^2*s0^3 +
 + 147*s4*s2*s0^3 - 20*s3*s2*s0^3 + 28*s2^2*s0^3 - 9*s4*s1*s0^3 - 21*s3*s1*s0^3 - 5*s2*s1*s0^3 +
 - 4*s1^2*s0^3 + 25*s4*s0^4 - 4*s3*s0^4 + 9*s2*s0^4 - s1*s0^4 + s0^5;


###################
###################

### Order 6

# Variants based on:
# Type sq-free: cos(2*pi/13), cos(2*pi/21);
# Type sq: cos(pi/7), cos(pi/8), cos(pi/9), cos(pi/12);
# Other: from cos(pi/10);

### Based on cos(2*pi/13)
r = 2*cos(2*pi*seq(6)/13)

#
x = r; poly.calc0(x)
x^6 + x^5 - 5*x^4 - 4*x^3 + 6*x^2 + 3*x - 1

#
x = r^2 - r; poly.calc0(x)
x^6 - 12*x^5 + 47*x^4 - 69*x^3 + 32*x^2 + 3*x - 1

#
x = r^3 - r^2 - r; poly.calc0(x)
x^6 + 14*x^5 + 47*x^4 + 22*x^3 - 20*x^2 - 10*x - 1

#
x = r^3 + r^2 - r; poly.calc0(x)
x^6 - 8*x^5 + 5*x^4 + 20*x^3 - 20*x^2 + 2*x + 1

#
x = r^4 - 2*r^2 + r; poly.calc0(x, digits = 5)
x^6 - 8*x^5 + 5*x^4 + 46*x^3 + 19*x^2 - 11*x + 1

#
x = r^4 + r^3 - 2*r^2 - 2*r; poly.calc0(x, digits = 5)
x^6 - 7*x^5 + 2*x^4 + 33*x^3 + 2*x^2 - 7*x + 1

#
x = r^4 + r^3 - 2*r^2 - r; poly.calc0(x, digits = 5)
x^6 - 6*x^5 - 11*x^4 + 6*x^3 + 15*x^2 + 7*x + 1


### Based on cos(2*pi/21)
id = c(1,2,4,5,8,10)
r = 2*cos(2*pi*id/21)


x = r; poly.calc0(x)
x^6 - x^5 - 6*x^4 + 6*x^3 + 8*x^2 - 8*x + 1

x = r^2 - r; poly.calc0(x)
x^6 - 12*x^5 + 46*x^4 - 62*x^3 + 16*x^2 + 11*x + 1

x = r^3 - 2*r; poly.calc0(x)
x^6 + x^5 - 13*x^4 - 6*x^3 + 15*x^2 + 8*x + 1

x = r^3 + r^2 - r; poly.calc0(x)
x^6 - 13*x^5 + 29*x^4 + 78*x^3 - 48*x^2 + x + 1

x = r^4 - 2*r^2 + r; poly.calc0(x)
x^6 - 16*x^5 + 74*x^4 - 57*x^3 - 180*x^2 + 10*x + 1

x = r^4 - r^2 + r - 1; poly.calc0(x)
x^6 - 23*x^5 + 158*x^4 - 225*x^3 - 579*x^2 - 214*x + 1


### Based on cos(pi/7)
r = 2*cos(pi*seq(6)/7)

# Trivial Poly:
x = r; poly.calc0(x)
x^6 - 5*x^4 + 6*x^2 - 1

# Non-Trivial:
x = r^2 + r; poly.calc0(x, digits = 5)
x^6 - 10*x^5 + 32*x^4 - 38*x^3 + 13*x^2 + 2*x - 1

x = r^2 + 2*r + 1; poly.calc0(x, digits = 5)
x^6 - 16*x^5 + 82*x^4 - 154*x^3 + 101*x^2 - 22*x + 1

x = (r + 1)^2 - (r^2 + r)^2; poly.calc0(x, digits = 5)
x^6 + 20*x^5 + 38*x^4 - 70*x^3 - 107*x^2 + 26*x + 1


### Based on cos(pi/8)
id = seq(7)[-4]
r = 2*cos(pi*id/8)

# Trivial Poly
x = r; poly.calc0(x)
x^6 - 6*x^4 + 10*x^2 - 4

# Less-Trivial
x = r^2 + r - 1; poly.calc0(x)
x^6 - 6*x^5 + 5*x^4 + 12*x^3 - 7*x^2 - 2*x + 1
(x^2 - 2*x - 1) * (x^4 - 4*x^3 - 2*x^2 + 4*x - 1)


### Based on cos(pi/9)
id = seq(8)[ - c(3,6)]
r = 2*cos(pi*id/9)

# Trivial Poly
x = r; poly.calc0(x)
x^6 - 6*x^4 + 9*x^2 - 1

# Non/Less-Trivial
x = (r + 1)^2; poly.calc0(x, digits = 5)
x^6 - 18*x^5 + 105*x^4 - 226*x^3 + 198*x^2 - 72*x + 9
# Note: "Squared" Poly;
x = sqrt(x); # P(x^2) is a shifted trivial poly;
(x^6 + 9*x^4 - 12*x^2 + 3)^2 - (6*x^5 - 4*x^3)^2

#
x = - r^3 + r^2 + 3*r - 1; poly.calc0(x, digits = 5)
x^6 - 6*x^5 + 6*x^4 + 18*x^3 - 33*x^2 + 12*x - 1

#
x = (r^2 - r - 1)^2; poly.calc0(x, digits = 5)
x^6 - 30*x^5 + 225*x^4 - 398*x^3 + 210*x^2 - 36*x + 1


### Derived from cos(pi/10)
id = seq(9)[-5]
r0 = 2*cos(pi*id/10)
r = r0[ - c(2,5)]

x = r^2 - r; poly.calc0(x, digits = 5)
x^6 - 14*x^5 + 69*x^4 - 140*x^3 + 95*x^2 + 10*x - 5

x = r^3 - r^2 - r + 1; poly.calc0(x, digits = 5)
x^6 + 10*x^5 + 9*x^4 - 76*x^3 - 41*x^2 + 18*x - 1


### Based on cos(pi/12)
id = c(1,2,5,7,10,11)
r = 2*cos(pi*id/12)

# Trivial Poly
x = r; poly.calc0(x)
x^6 - 7*x^4 + 13*x^2 - 3

# Non-Trivial
x = r^2 + r - 1; poly.calc0(x)
x^6 - 8*x^5 + 13*x^4 + 16*x^3 - 19*x^2 + 1


####################
####################

### Order 7

# TODO

### Based on cos(2*pi/21)
# Trivial Root: see also Order 6;
id = seq(10)[ - c(3,6,9)]
r = 2*cos(2*pi*id/21)

x = r; poly.calc0(x, digits = 5)
x^7 - 7*x^5 + 14*x^3 - 7*x + 1


####################
####################

### Order 8

### Based on cos(2*pi/17)
r = 2*cos(2*pi*seq(8)/17)

x = r; poly.calc0(x, digits = 5)
x^8 + x^7 - 7*x^6 - 6*x^5 + 15*x^4 + 10*x^3 - 10*x^2 - 4*x + 1

x = r^2 - r; poly.calc0(x, digits = 5)
x^8 - 16*x^7 + 95*x^6 - 261*x^5 + 338*x^4 - 177*x^3 + 7*x^2 + 13*x + 1

# TODO: other types;


### Based on cos(pi/10)
id = seq(9)[-5]
r = 2*cos(pi*id/10)

# Trivial Poly
x = r; poly.calc0(x)
x^8 - 8*x^6 + 21*x^4 - 20*x^2 + 5

x = (r+1)^2; poly.calc0(x, digits = 5)
x^8 - 24*x^7 + 214*x^6 - 876*x^5 + 1655*x^4 - 1316*x^3 + 414*x^2 - 44*x + 1


### Based on cos(pi/15)
id = seq(14)[ - c(3,5,6,9,10,12)]
r = 2*cos(pi*id/15)

# Trivial Poly
x = r; poly.calc0(x, digits = 5)
x^8 - 9*x^6 + 26*x^4 - 24*x^2 + 1

# Non-Trivial:
x = r^2 - 2*r; poly.calc0(x, digits = 5)
x^8 - 18*x^7 + 97*x^6 - 100*x^5 - 274*x^4 - 50*x^3 + 48*x^2 + 16*x + 1


####################
####################

### Order 9

### Based on cos(2*pi/19)
r = 2*cos(2*pi * seq(9) / 19)

#
x = r; poly.calc0(x, digits = 5)
x^9 + x^8 - 8*x^7 - 7*x^6 + 21*x^5 + 15*x^4 - 20*x^3 - 10*x^2 + 5*x + 1

x = r^2 - r; poly.calc0(x, digits = 5)
x^9 - 18*x^8 + 125*x^7 - 425*x^6 + 743*x^5 - 631*x^4 + 189*x^3 + 28*x^2 - 14*x + 1

x = r^3 + r^2 - r; poly.calc0(x, digits = 5)
x^9 - 14*x^8 + 47*x^7 + 18*x^6 - 175*x^5 + 116*x^4 + 50*x^3 - 46*x^2 + x + 1

# TODO: remaining types;


### Based on cos(2*pi/27)
id = seq(13)[ - c(3,6,9,12)]
r = 2*cos(2*pi*id/27)

#
x = r; poly.calc0(x, digits = 5)
x^9 - 9*x^7 + 27*x^5 - 30*x^3 + 9*x + 1

x = r^2 - r; poly.calc0(x, digits = 5)
x^9 - 18*x^8 + 126*x^7 - 438*x^6 + 801*x^5 - 747*x^4 + 300*x^3 - 18*x^2 - 9*x + 1

