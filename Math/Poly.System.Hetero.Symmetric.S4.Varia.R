########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric: Various
###
### draft v.0.1b


### Various Hetero-Symmetric Systems


####################

### Helper Functions

source("Polynomials.Helper.R")


####################
####################

### S4 based on:
# James H. Davenport. SSc2017: Example for Groebner Basis
# https://people.bath.ac.uk/masjhd/Slides/SC2School2017/SC2-G4.pdf

x1 + x2 + x3 + x4 # = 0
x1*x2 + x2*x3 + x3*x4 + x4*x1 # = 0
x1*x2*x3 + x2*x3*x4 + x3*x4*x1 + x4*x1*x2 # = 0
x1*x2*x3*x4 # = R4

### Solution:

### Trivial Cases:

### T.1) x1 = x2 = x3 = x4
# - NO solutions;

### T.2) x1 = x3; x2 = x4;
# - NO solutions;

### T.3) x3 = - x1; x4 = - x2;
# - system is undetermined;
# - infinitely many solutions;
(x1*x2)^2 - R4 # = 0


### Examples:

R = 3
x1 = sqrt(2)
#
x2 = sqrt(R/x1^2);
x2 = c(x2, - x2); x1 = c(x1, x1);
x3 = - x1; x4 = - x2;

### Test:
x1 + x2 + x3 + x4 # = 0
x1*x2 + x2*x3 + x3*x4 + x4*x1 # = 0
x1*x2*x3 + x2*x3*x4 + x3*x4*x1 + x4*x1*x2 # = 0
x1*x2*x3*x4 # = R4

### Derivation:

### Eq 2 =>
(x1 + x3)*(x2 + x4) # = 0
# =>
# x3 = - x1;
# => Eq 1 =>
# x2 + x4 = 0;
# x4 = - x2;


###################
### Generalization:

### G.1.) x1 + x2 + x3 + x4 = R1
# x3 = - x1; x2 + x4 = R1;
# - NO solutions;

### G.2.) (x1 + x3)*(x2 + x4) = R2
# (x1 + x3)^2 = - R2;


### Derivation:

### G.1.) x1 + x2 + x3 + x4 = R1

### Eq 2 =>
x3 = - x1;
# Eq 1 =>
x2 + x4 - R1 # = 0

### Eq 3:
x2*x3*x4 + x4*x1*x2 # = 0 =>
x1*x2*x3 + x3*x4*x1 # = 0 =>
x2 + x4 # = 0 # Contradiction!
# => NO solutions!


### G.2.) (x1 + x3)*(x2 + x4) = R2
# (x1 + x3)^2 = - R2;

# TODO

