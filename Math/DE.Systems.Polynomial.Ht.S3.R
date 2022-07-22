########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### DE Systems: Polynomial
### Base: Hetero-Symmetric
###
### draft v.0.1a


####################

### Helper Functions

# include: DE.ODE.Helper.R;
source("DE.ODE.Helper.R")
source("DE.ODE.Helper.Math.R")


########################
########################
########################

########################
### Hetero-Symmetric ###
########################

### S3 ###

### Base-System
# y1^n = n*b1*y2 + n*R1
# y2^n = n*b2*y3 + n*R2
# y3^n = n*b3*y1 + n*R3

### Note:
# - the original polynomial system is much easier to solve when:
#   R1 = R2 = R3 = R(x), and
#   b1 = b2 = b3 = b(x);
#  [see Poly.System.Hetero.Symmetric.S3.R]


### D =>
y1^(n-1)*dy1 - b1*dy2 - db1*y2 - dR1 # = 0
# * y1 =>
(n*b1*y2 + n*R1)*dy1 - b1*y1*dy2 - db1*y1*y2 - dR1*y1 # = 0

### ODE System:
n*b1*y2*dy1 - b1*y1*dy2 - db1*y1*y2 + n*R1*dy1 - dR1*y1 # = 0
n*b2*y3*dy2 - b2*y2*dy3 - db2*y2*y3 + n*R2*dy2 - dR2*y2 # = 0
n*b3*y1*dy3 - b3*y3*dy1 - db3*y3*y1 + n*R3*dy3 - dR3*y3 # = 0

# TODO: check;

