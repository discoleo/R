########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Hetero-Symmetric S4: Mixed
### E2-Type: Asymmetric
###
### draft v.0.1a


### E2-Type:
# x1^k*x2^n + x2^k*x3^n + x3^k*x4^n + x4^k*x1^n = R2;
# where: n != k (Asymmetric);


### History

### v.0.1a
# - moved E2 types (Asymmetric) from file:
#   Poly.System.S4.HtMixed.Basic.R;
# - for basic Theory: see previous file;


####################
####################

### Helper Functions

source("Poly.System.S4.HtMixed.Basic.Helper.R")


###################


###################
### Type: E21a  ###
###################

###############
### Order 1 ###
###############

x1 + x2 + x3 + x4 - R1 # = 0
x1^2*x2 + x2^2*x3 + x3^2*x4 + x4^2*x1 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

# TODO: complicated;


### Derivation:

E21b + E12b - (E2 - E2a)*S + E3 # = 0
# ???
# Epoly.distinct(c(2,1), v=4)
E21 - S*E2 + 3*E3 = 0
 
