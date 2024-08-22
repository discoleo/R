
## Polynomials: Class 3
## Roots: All Real


### Various Examples

# Based on:
# => cos(2*pi / (2*n+1));
# => subset of 2 * cos(2*pi / higher);
# => cos(pi / (n+1));


###################

### Helper Functions

source("Polynomials.Helper.R")


####################
####################

### Order 4

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

### Order 6
# Examples

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

