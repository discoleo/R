########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S6: C3-Hetero-Symmetric
### Useful Formulas
###
### draft v.0.1a


### Formulas:

# - Formulas & derivations;
# - Useful for Ht & S6/C3 Systems;


####################

### Helper Functions


source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


### Debug
x = sqrt(c(2,3,5, 7,11,13));
x[3] = - x[3]; x[5] = -x[5];
x1 = x[1]; x2 = x[2]; x3 = x[3];
y1 = x[4]; y2 = x[5]; y3 = x[6];

### Notation:
sx = x1 + x2 + x3; sy = y1 + y2 + y3;
px = x1 * x2 * x3; py = y1 * y2 * y3;
E2x = (x1 + x2)*x3 + x1*x2;
E2y = (y1 + y2)*y3 + y1*y2;


###########################

### Basic

A1 = x1*y1 + x2*y2 + x3*y3;
B1 = x1*y2 + x2*y3 + x3*y1;
C1 = x1*y3 + x2*y1 + x3*y2;
E2ABC = (A1 + B1)*C1 + A1*B1;

### Sum:
A1 + B1 + C1 - sx*sy # = 0

### E2:
E2ABC - E2y*(x1^2 + x2^2 + x3^2) - E2x*(y1^2 + y2^2 + y3^2) - E2x*E2y # = 0
E2ABC - E2y*sx^2 - E2x*sy^2 + 3*E2x*E2y # = 0

### E3:
A1*B1*C1 - py*(x1^3 + x2^3 + x3^3) - px*(y1^3 + y2^3 + y3^3) - 3*px*py +
	- (x1^2*x2 + x2^2*x3 + x3^2*x1)*(y1^2*y2 + y2^2*y3 + y3^2*y1) +
	- (x1*x2^2 + x2*x3^2 + x3*x1^2)*(y1*y2^2 + y2*y3^2 + y3*y1^2) # = 0

# [intermediary]
A2a = x1^2*x2 + x2^2*x3 + x3^2*x1;
A2b = x1*x2^2 + x2*x3^2 + x3*x1^2;
B2a = y1^2*y2 + y2^2*y3 + y3^2*y1;
B2b = y1*y2^2 + y2*y3^2 + y3*y1^2;

### Sum:
A2a + A2b - E2x*sx + 3*px # = 0
B2a + B2b - E2y*sy + 3*py # = 0

### Prod:
A2a * A2b - (x1^3*x2^3 + x2^3*x3^3 + x3^3*x1^3) +
	- px*(x1^3 + x2^3 + x3^3) - 3*px^2 # = 0
A2a * A2b - (- 3*E2x*px*sx + E2x^3 + 3*px^2) +
	- px*(sx^3 - 3*E2x*sx + 3*px) - 3*px^2 # = 0
A2a * A2b + 6*E2x*px*sx - E2x^3 - 9*px^2 - px*sx^3 # = 0

# similarly:
B2a * B2b + 6*E2y*py*sy - E2y^3 - 9*py^2 - py*sy^3 # = 0

# =>
A1*B1*C1 - py*(sx^3 - 3*E2x*sx + 3*px) - px*(sy^3 - 3*E2y*sy + 3*py) - 3*px*py +
	- A2a*B2a - (E2x*sx - 3*px - A2a)*(E2y*sy - 3*py - B2a) # = 0
A1*B1*C1 - py*(sx^3 - 3*E2x*sx) - px*(sy^3 - 3*E2y*sy) - 18*px*py +
	- 2*A2a*B2a + (E2x*sx - 3*px)*B2a + (E2y*sy - 3*py)*A2a +
	- E2x*sx*E2y*sy + 3*E2x*py*sx + 3*E2y*px*sy # = 0


# TODO

