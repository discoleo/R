########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### DE Systems: Non-Polynomial
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

### Basic ###

### Base-System:
# x*atan(y1) = b1*y2 + R1
# x*atan(y2) = b2*y1 + R2

### Note:
# - the original polynomial system is easier to solve when:
#   R1 = R2 = R(x), and
#   b1 = b2 = b(x);


### D =>
atan(y1) + x/(y1^2 + 1)*dy1 - b1*dy2 - db1*y2 - dR1 # = 0
# *x*(y1^2 + 1) =>
x^2*dy1 - b1*x*(y1^2 + 1)*dy2 +
	+ (b1*y2 + R1)*(y1^2 + 1) - db1*x*y2*(y1^2 + 1) - x*dR1*(y1^2 + 1) # = 0
b1*x*y1^2*dy2 - x^2*dy1 + b1*x*dy2 +
	- (b1 - db1*x)*y1^2*y2 - (R1 - x*dR1)*y1^2 +
	- (b1 - db1*x)*y2 - (R1 - x*dR1) # = 0

### ODE System:
b1*x*y1^2*dy2 - x^2*dy1 + b1*x*dy2 +
	- (b1 - db1*x)*y1^2*y2 - (R1 - x*dR1)*y1^2 +
	- (b1 - db1*x)*y2 - (R1 - x*dR1) # = 0
b2*x*y2^2*dy1 - x^2*dy2 + b2*x*dy1 +
	- (b2 - db2*x)*y1*y2^2 - (R2 - x*dR2)*y2^2 +
	- (b2 - db2*x)*y1 - (R2 - x*dR2) # = 0

# TODO: check;


################
################

### Compound ###

### Base-System:
# y2*atan(y1) = b1*y1 + R1
# y1*atan(y2) = b2*y2 + R2

### Note:
# - the original polynomial system is easier to solve when:
#   R1 = R2 = R(x), and
#   b1 = b2 = b(x);

# TODO

