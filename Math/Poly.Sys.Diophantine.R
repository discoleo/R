########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Diophantine Equations
###
### draft v.0.1k


####################
####################

### Helper Functions

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;

source("Polynomials.Helper.R")


### Elliptic Curves: Roots
roots.elliptic = function(xy, b, scale=1) {
	# b = c(b3, b2, b1, b0);
	if(length(b) == 3) b = c(b, 0);
	if(length(xy) > 2) {
		slope = (xy[4] - xy[3]) / (xy[2] - xy[1]);
		y0 = xy[3] - slope*xy[1];
		isNotOrigin = TRUE;
	} else {
		slope = xy[2] / xy[1];
		isNotOrigin = TRUE;
	}
	bc = b; bc[2] = bc[2] - slope^2;
	if(isNotOrigin) {
		bc[3] = bc[3] - 2*slope*y0;
		bc[4] = bc[4] - y0^2;
	}
	r = roots(bc);
	y = sapply(r, function(r) {
		sqrt(sum(b * r^seq(3, 0)));
	})
	sol = cbind(x=r, y=y);
	sol * scale;
}

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

#################
### Special Cases

### Case 1:
### x^3 - (p0^2 + 1)*x
p0 = 6
p = p0^2 + 1;
x = p; y = p*p0;
y^2 - x^3 + p*x


### Case 2:
### x^3 + (p0^2 - 1)*x
p0 = 6
p = p0^2 - 1;
x = p; y = p*p0;
y^2 - x^3 - p*x


### Case 3: Squares
p = 11
x0 = 6; y0 = sqrt(x^3 - p*x);
x = x0; y = y0;
y^2 - x^3 + p*x
# Squaring =>
x = x0^2; y = y0^2;
print(c(x, y))
y^2 - x^3 + 2*p*x^2 - p^2*x


### Particular Examples:

### Ex 1:
### x^3 - 7*x
x = c(4, 0, -7/4); y = c(6, 0, -21/8);
y^2 - x^3 + 7*x
# Self-Shift =>
x = 0; y = 0; # trivial solution
x = -4; y = -6; # initially trivial, now shifted
x = c(-4, -23/4); y = c(-6, -69/8);
(y + 6)^2 - (x + 4)^3 + 7*(x + 4)
y^2 + 12*y - x^3 - 12*x^2 - 41*x


### Ex 2:
x0 = c(12, 1/4); y0 = c(30, 5/8);
x = x0; y = y0;
y^2 - x^3 + 6*x^2 - 3*x
# Variant: Shift =>
x = x0 - 2; y = y0;
y^2 - x^3 + 9*x + 10


roots.elliptic(c(12, 30), c(1,-6,3))


### Generator Techniques

### Swap with Pythagorean Triple

# Base
x = c(0, 4, -15/4); y = c(1, 3, 7/8);
y^2 - x^3 + 14*x - 1
# Swap: y^2 => z^2 - x^2
x = c(0, 4); y = c(1, 5)
y^2 - x^3 - x^2 + 14*x - 1


### Multiplicative Generator

### Base
# Note: x[1:2] are the same;
x1 = c(0, 3, 1); y1 = c(0, 6, 2);
x = x1; y = y1;
y^2 - x*(x^2 + 3)
x2 = c(0, 3, 13); y2 = c(0, 12, 52);
x = x2; y = y2;
y^2 - x*(x^2 + 39)
# Prod =>
x = x1[1:2]; y = (y1*y2)[1:2];
y^2 - x^2*(x^2 + 3)*(x^2 + 39)
x = x^2; x = c(x, 13); y = c(y, 104);
y^2 - x*(x + 3)*(x + 39)
y^2 - x*(x^2 + 42*x + 117)
### Variants:
x0 = x; y0 = y;
### V1: Shift
x = x0 + 14; y = y0;
y^2 - x^3 + 3*14^2*x - 117*x + 14*117 - 2*14^3
y^2 - x^3 + 3*157*x - 14*275


### Ex 2:
# Note: xb & yb are the same;
xb = 3; yb = 7;
x = c(xb, 0, -2); y = c(yb, 4, 2);
y^2 - x*(x^2 + 2) - 16
#
x = c(xb, -1, 1/4); y = c(yb, 1, 23/8);
y^2 - x*(x^2 + 5) - 7
# Prod =>
(y^2 - 16)*(y^2 - 7) - x^2*(x^2 + 2)*(x^2 + 5)
x = c(xb^2, 0, -23/9); y = c(yb^2, 16, 179/27);
(y - 16)*(y - 7) - x*(x + 2)*(x + 5)
y^2 - 23*y - x*(x^2 + 7*x + 10) + 7*16
y^2 - 23*y - x^3 - 7*x^2 - 10*x + 112
### Variants:
x0 = x; y0 = y;
# V1: Shifts
x = x0; y = y0 - 23/2;
(y + 23/2)^2 - 23*(y + 23/2) - x^3 - 7*x^2 - 10*x + 112
y^2 - x^3 - 7*x^2 - 10*x - 81/4
# V2: Self-Shift
x0s = x0[1]; y0s = y0[1];
x = x0 - x0s; y = y0 - y0s;
(y+y0s)^2 - 23*(y+y0s) - (x+x0s)^3 - 7*(x+x0s)^2 - 10*(x+x0s) + 112
y^2 + 2*49*y - 23*y - x^3 - 27*x^2 - 3*81*x - 7*x^2 - 7*18*x - 10*x
y^2 + 75*y - x^3 - 34*x^2 - 379*x
cbind(x, y)


roots.elliptic(c(0, 3, 4, 7), c(1, 0, 2, 16))
roots.elliptic(c(-1, 3, 1, 7), c(1, 0, 5, 7))
roots.elliptic(c(9,0, 37.5,4.5), c(1, 7, 10, 81/4))


### Ex 3:
k = c(1, 2)
pE = toPoly.pm("(y - b^2 + a^3 + k[1]*a^2)*(y - b^2 + a^3 + k[2]*a^2) - x*(x + k[1]*a)*(x + k[2]*a)")
pE = sort.pm(pE, c("y", "x"))
print.pm(pE, lead="y")

a = 5
b = 2
x = a^2; y = b^2;
y^2 + (2*a^3 + 3*a^2 - 2*b^2)*y - x^3 - 3*a*x^2 - 2*a^2*x + a^6 + 3*a^5 - 2*a^3*b^2 + 2*a^4 + b^4 - 3*a^2*b^2
# the Example:
y^2 + 317*y - x^3 - 15*x^2 - 50*x + 24966
# TODO: shifts;


#######################
#######################

#######################
### Diophantine Eqs ###
#######################

### x^3 + y^3 + z^3 = 3*x*y*z

### Examples

x = 16; y = 9; z = -25;
x^3 + y^3 + z^3 - 3*x*y*z

