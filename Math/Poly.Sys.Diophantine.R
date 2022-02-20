########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Diophantine Equations
###
### draft v.0.1b-Sol1


####################
####################

### Helper Functions

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;

source("Polynomials.Helper.R")


##############################


### Elliptic Curves

# History of elliptic curves rank records
# https://web.math.pmf.unizg.hr/~duje/tors/rankhist.html

library(gmp)

### Examples

x0 = c(-8, -1, 49/4)
y0 = c(12, 9, 231/8)

### Base:
x = x0; y = y0;
y^2 - x^3 + 82*x

### 1st Square:
x = x0^2; y = y0^2;
y^2 - x^3 + 2*82*x^2 - 82^2*x

###
x = x0^2 - 164/3; y = y0^2;
y^2 - x^3 + 1/12*164^2*x - 1/4*(164/3)^3

###
x = 9*x0^2 - 3*164; y = 27*y0^2;
y^2 - x^3 + 27*82^2*x - 54*82^3

### 2nd Square:
x = as.bigq(9*x0^2 - 3*164)^2; y = as.bigq(27*y0^2)^2;
y^2 - 108*82^3*y - x^3 + 54*82^2*x^2 - 27^2*82^4*x + 54^2*82^6

###
x = as.bigq(9*x0^2 - 3*164)^2; y = as.bigq(27*y0^2)^2 - 54*82^3;
y^2 - x^3 + 54*82^2*x^2 - 27^2*82^4*x

###
x = (as.bigq(9*x0^2 - 3*164)/82)^2; y = (as.bigq(27*y0^2)^2 - 54*82^3)/82^3;
y^2 - x^3 + 54*x^2 - 27^2*x
y^2 - x*(x - 27)^2


### Sol
p = 82;
x = (as.bigq(9*x0^2)/p - 6)^2; y = (as.bigq(27*y0^2)^2/p^3 - 54);
y^2 - x*(x - 27)^2

### p = p0^2 + 1;
p = 37
# TODO: Factor *9;
k = p*9 - 6;
x = sqrt((k + 6)*p) / 3;
y = sqrt((k-3)*p*sqrt((k+6)*p)/3) / 3
y^2 - x^3 + p*x


####################

### x^3 + y^3 + z^3 = 3*x*y*z

### Examples

x = 16; y = 9; z = -25;
x^3 + y^3 + z^3 - 3*x*y*z

