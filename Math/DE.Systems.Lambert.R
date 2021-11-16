########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### DE Systems: Lambert
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


####################

### Helper functions

library(pracma)
# needed for Lambert W;


# include: DE.ODE.Helper.R;
source("DE.ODE.Helper.R")


#########################
#########################

### Hetero:
### Symmetric & non-Symmetric

# exp(y1) = b1*y2 + R1
# exp(y2) = b2*y1 + R2

# Note:
# - "slightly easier" to solve when:
#   b1 = b2, R1 = R2;
#   => trivial solution: y1 = y2;

### D =>
exp(y1)*dy1 - b1*dy2 - db1*y2 - dR1 # = 0
exp(y2)*dy2 - b2*dy1 - db2*y1 - dR2 # = 0

# DE System:
(b1*y2 + R1)*dy1 - b1*dy2 - db1*y2 - dR1 # = 0
(b2*y1 + R2)*dy2 - b2*dy1 - db2*y1 - dR2 # = 0

### Special Cases:
# b1 = b2 = ct;
(b*y2 + R1)*dy1 - b*dy2 - dR1 # = 0
(b*y1 + R2)*dy2 - b*dy1 - dR2 # = 0

# TODO: check & solve;

