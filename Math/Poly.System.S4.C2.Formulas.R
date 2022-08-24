########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: Hetero-Symmetric
### Useful Formulas
###
### draft v.0.1a


### Formulas:

# - Formulas & derivations;
# - Useful for Ht & C2 Systems;

# - Applicable for systems described in:
#   Poly.System.S4.C2.Symmetric.R;
#   Poly.System.S4.C2.R;
#   Poly.System.Hetero.Symmetric.S4.R;


####################

### Helper Functions


source("Poly.System.S4.C2.Helper.R")


####################
####################

### Notation:
s1 = x1 + x3; s2 = x2 + x4;
p1 = x1 * x3; p2 = x2 * x4;
sp = p1 + p2; ps = s1 * s2;
S  = s1 + s2;
E4 = p1 * p2;


#################
#################

#################
### Composite ###
#################

### Formula for:
### x1^2*x2 + x2^2*x3 + x3^2*x4 + x4^2*x1 = R

### Derivation:

### Part Eq 1:
A1 = x1^2*x2 + x3^2*x4;
B1 = x1^2*x4 + x3^2*x2;

### A1 + B1 =>
A1 + B1 - (x1^2 + x3^2)*(x2 + x4) # = 0;
A1 + B1 - (s1^2 - 2*p1)*s2 # = 0;

### A1 * B1 =>
A1 * B1 - x2*x4*(x1^4 + x3^4) - (x1*x3)^2*(x2^2 + x4^2) # = 0
A1 * B1 - p2*(s1^4 - 4*p1*s1^2 + 2*p1^2) - p1^2*(s2^2 - 2*p2) # = 0

### Eq 1 =>
A1*((s1^2 - 2*p1)*s2 - A1) - p2*(s1^4 - 4*p1*s1^2 + 2*p1^2) - p1^2*(s2^2 - 2*p2) # = 0
A1^2 - A1*(s1^2 - 2*p1)*s2 + p2*(s1^4 - 4*p1*s1^2 + 2*p1^2) + p1^2*(s2^2 - 2*p2) # = 0
A1^2 - A1*(s1^2 - 2*p1)*s2 + p2*s1^4 - 4*p1*p2*s1^2 + p1^2*s2^2 # = 0
# =>
A1^2 - A1*(s1^2 - 2*p1)*s2 + p2*s1^4 - 4*E4*s1^2 + p1^2*s2^2 # = 0


### Part Eq 2:
# - similarly:
A2 = x2^2*x3 + x4^2*x1;
B2 = x2^2*x3 + x4^2*x1;

# =>
A2^2 - A2*(s2^2 - 2*p2)*s1 + p1*(s2^4 - 4*p2*s2^2 + 2*p2^2) + p2^2*(s1^2 - 2*p1) # = 0
# =>
A2^2 - A2*(s2^2 - 2*p2)*s1 + p1*s2^4 - 4*E4*s2^2 + p2^2*s1^2 # = 0


### Initial Eq:
# A1 + A2 = R;

A1 + A2 - R # = 0
A1^2 - A1*(s1^2 - 2*p1)*s2 + p2*s1^4 - 4*E4*s1^2 + p1^2*s2^2 # = 0
A2^2 - A2*(s2^2 - 2*p2)*s1 + p1*s2^4 - 4*E4*s2^2 + p2^2*s1^2 # = 0

# TODO:
# - check & derive final equation;

###
p1 = toPoly.pm("A1 + A2 - R")
p2 = toPoly.pm("A1^2 - A1*(s1^2 - 2*p1)*s2 + p2*s1^4 - 4*E4*s1^2 + p1^2*s2^2")
p3 = toPoly.pm("A2^2 - A2*(s2^2 - 2*p2)*s1 + p1*s2^4 - 4*E4*s2^2 + p2^2*s1^2")

pR = solve.lpm(p1, p2, p3, xn=c("A1", "A2"))
str(pR)
# 99 Monomials: seems a lot!

