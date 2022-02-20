########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Diophantine Equations
###
### draft v.0.1b


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


### Examples

x0 = c(-8, -1, 49/4)
y0 = c(12, 9, 231/8)

### Base:
x = x0; y = y0;
y^2 - x^3 + 82*x

###
x = x0^2; y = y0^2;
y^2 - x^3 + 2*82*x^2 - 82^2*x

###
x = x0^2 - 164/3; y = y0^2;
y^2 - x^3 + 1/12*164^2*x - 1/4*(164/3)^3

###
x = 9*x0^2 - 3*164; y = 27*y0^2;
y^2 - x^3 + 27*82^2*x - 54*82^3



####################

### x^3 + y^3 + z^3 = 3*x*y*z

### Examples

x = 16; y = 9; z = -25;
x^3 + y^3 + z^3 - 3*x*y*z

