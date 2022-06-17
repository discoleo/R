########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### DE Systems: S3 Lambert-type
###
### draft v.0.1a


#############
### Types ###
#############

### Symmetric:
# - TODO
### Hetero-Symmetric:
# - Simple, derived from:
#   exp(y1) = b1*y2 + R1;
# - 2 Exponentials, derived from:
#   exp(y1) + exp(y2) = b1*y3 + R1;


####################

### Helper functions

source("Polynomials.Helper.R")


# include: DE.ODE.Helper.R;
source("DE.ODE.Helper.R")


#########################
#########################

### Hetero:
### Symmetric & non-Symmetric

### Base-System:
# exp(y1) + exp(y2) = b1*y3 + R1
# exp(y2) + exp(y3) = b2*y1 + R2
# exp(y3) + exp(y1) = b3*y2 + R3

# Note:
# - "slightly easier" to solve when:
#   b1 = b2 = b3, and R1 = R2 = R3;
#   => trivial solution: y1 = y2 = y3;
#   => Non-trivial solution: ???

### D =>
exp(y1)*dy1 + exp(y2)*dy2 - b1*dy3 - db1*y3 - dR1 # = 0
exp(y2)*dy2 + exp(y3)*dy3 - b2*dy1 - db2*y1 - dR2 # = 0
exp(y3)*dy3 + exp(y1)*dy1 - b3*dy2 - db3*y2 - dR3 # = 0

# assuming (b1, b2, b3) = ct;
exp(y1)*dy1 + exp(y2)*dy2 - b1*dy3 - dR1 # = 0
exp(y2)*dy2 + exp(y3)*dy3 - b2*dy1 - dR2 # = 0
exp(y3)*dy3 + exp(y1)*dy1 - b3*dy2 - dR3 # = 0

### Solve Linear System
# - "linear" in: exp(y1), exp(y2), exp(y3);
# Det = 2*dy1*dy2*dy3;

### exp(y1) =
(- b2*dy1 + b3*dy2 + b1*dy3 + dR1 - dR2 + dR3) / (2*dy1);
### exp(y2) =
(+ b2*dy1 - b3*dy2 + b1*dy3 + dR1 + dR2 - dR3) / (2*dy2);
### exp(y3) =
(+ b2*dy1 + b3*dy2 - b1*dy3 - dR1 + dR2 + dR3) / (2*dy3);

### Substitute in Base-system:
exp(y1) + exp(y2) - b1*y3 - R1 # = 0
-2*(b1*y3 + R1)*dy1*dy2 +
	(- b2*dy1 + b3*dy2 + b1*dy3 + dR1 - dR2 + dR3)*dy2 +
	(b2*dy1 - b3*dy2 + b1*dy3 + dR1 + dR2 - dR3)*dy1 # = 0
# =>
# ...


### Special Case:
# b1 = b2 = b3 = b = const,
# R1 = R2 = R3 = R
2*b*y3*dy1*dy2 - b*dy1^2 - b*dy2^2 +
	+ 2*R*dy1*dy2 + 2*b*dy1*dy2 - b*dy1*dy3 - b*dy2*dy3 +
	- dR*dy1 - dR*dy2 # = 0

# TODO: check;

