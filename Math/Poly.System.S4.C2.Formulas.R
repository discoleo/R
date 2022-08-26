########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: Hetero-Symmetric
### Useful Formulas
###
### draft v.0.1c


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
# x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];
s1 = x1 + x3; s2 = x2 + x4;
p1 = x1 * x3; p2 = x2 * x4;
sp = p1 + p2; ps = s1 * s2;
S  = s1 + s2;
E4 = p1 * p2;

E2 = x1*x2 + x1*x3 + x1*x4 + x2*x3 + x2*x4 + x3*x4;
E3 = x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4;


###################
###################

###################
### Elementary  ###
### Polynomials ###
###################

##########
### E2 ###
##########

### Eq:
sp + ps - E2 # = 0

### Derivation:
E2 = x1*x2 + x1*x3 + x1*x4 + x2*x3 + x2*x4 + x3*x4;
p1 + p2 + x1*x2 + x1*x4 + x2*x3 + x3*x4 - E2 # = 0
p1 + p2 + s1*s2 - E2 # = 0
# =>
sp + ps - E2 # = 0


##########
### E3 ###
##########

### Eq:
E3^2 - sp*S*E3 + ps*sp^2 + E4*S^2 - 4*ps*E4 # = 0

### Derivation:
E3 = x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4;
x1*x3*(x2 + x4) + x2*x4*(x1 + x3) - E3 # = 0
p1*s2 + p2*s1 - E3 # = 0

A = p1*s2 + p2*s1;
B = p1*s1 + p2*s2;

### A + B =>
A + B - (p1 + p2)*(s1 + s2) # = 0
A + B - sp*S # = 0

### A * B =>
A * B - s1*s2*(p1^2 + p2^2) - p1*p2*(s1^2 + s2^2) # = 0
A * B - ps*(sp^2 - 2*E4) - E4*(S^2 - 2*ps) # = 0

### =>
A*(sp*S - A) - ps*(sp^2 - 2*E4) - E4*(S^2 - 2*ps) # = 0
A^2 - sp*S*A + ps*sp^2 + E4*S^2 - 4*ps*E4 # = 0


##############

##############
### Helper ###
##############

### p1*s2 + p2*s1
A = p1*s2 + p2*s1;
# =>
A^2 - sp*S*A + ps*sp^2 + E4*S^2 - 4*ps*E4 # = 0

#####################

### p1*s2^2 + p2*s1^2
A = p1*s2^2 + p2*s1^2;
B = p1*s1^2 + p2*s2^2;
# =>
A^2 - sp*(S^2 - 2*ps)*A + E4*(S^4 - 4*ps*S^2) + ps^2*sp^2 # = 0

### A + B =>
A + B - (p1 + p2)*(s1^2 + s2^2) # = 0
A + B - sp*(S^2 - 2*ps) # = 0

### A * B =>
A * B - E4*(s1^4 + s2^4) - ps^2*(p1^2 + p2^2) # = 0
A * B - E4*(S^4 - 4*ps*S^2 + 2*ps^2) - ps^2*(sp^2 - 2*E4) # = 0
A * B - E4*(S^4 - 4*ps*S^2) - ps^2*sp^2 # = 0

### =>
A*(A - sp*(S^2 - 2*ps)) + E4*(S^4 - 4*ps*S^2) + ps^2*sp^2 # = 0
A^2 - sp*(S^2 - 2*ps)*A + E4*(S^4 - 4*ps*S^2) + ps^2*sp^2 # = 0


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


#################
#################

################
### E3-Types ###
################

#############
### E211a ###
#############

### Formula:
E211a = x1^2*x2*x3 + x2^2*x3*x4 + x3^2*x4*x1 + x4^2*x1*x2;
# TODO

### Derivation:
p1*(x1*x2 + x3*x4) + p2*(x3*x2 + x1*x4) - E211a # = 0

### Method 1:
A = x1*x2 + x3*x4;
B = x1*x4 + x3*x2;

### A + B =>
A + B - (x1+x3)*(x2+x4) # = 0
A + B - s1*s2 # = 0
A + B - ps # = 0

### A * B =>
A * B - p1*(x2^2 + x4^2) - p2*(x1^2 + x3^2) # = 0
A * B - p1*(s2^2 - 2*p2) - p2*(s1^2 - 2*p1) # = 0
A * B - p1*s2^2 - p2*s1^2 + 4*E4 # = 0

A1 = p1*s2^2 + p2*s1^2;
# System =>
A * B - A1 + 4*E4 # = 0
A1^2 - sp*(S^2 - 2*ps)*A1 + E4*(S^4 - 4*ps*S^2) + ps^2*sp^2 # = 0
# =>
A^2*B^2 - sp*S^2*A*B + 2*sp*ps*A*B + 8*E4*A*B + E4*S^4 + sp^2*ps^2 +
	- 4*E4*sp*S^2 - 4*E4*ps*S^2 + 8*E4*sp*ps + 16*E4^2 # = 0
# =>
A^2*(ps - A)^2 + (2*sp*ps + 8*E4 - sp*S^2)*A*(ps - A) + E4*S^4 + sp^2*ps^2 +
	- 4*E4*sp*S^2 - 4*E4*ps*S^2 + 8*E4*sp*ps + 16*E4^2 # = 0

A^4 - 2*ps*A^3 - (2*sp*ps + 8*E4 - sp*S^2 - ps^2)*A^2 +
	+ (2*sp*ps^2 + 8*E4*ps - ps*sp*S^2)*A +
	+ E4*S^4 + sp^2*ps^2 +
	- 4*E4*sp*S^2 - 4*E4*ps*S^2 + 8*E4*sp*ps + 16*E4^2 # = 0
# similarly:
B^4 - 2*ps*B^3 - (2*sp*ps + 8*E4 - sp*S^2 - ps^2)*B^2 +
	+ (2*sp*ps^2 + 8*E4*ps - ps*sp*S^2)*B +
	+ E4*S^4 + sp^2*ps^2 +
	- 4*E4*sp*S^2 - 4*E4*ps*S^2 + 8*E4*sp*ps + 16*E4^2 # = 0

### System:
p1*A + p2*B - E211a # = 0
A^4 - 2*ps*A^3 - (2*sp*ps + 8*E4 - sp*S^2 - ps^2)*A^2 +
	+ (2*sp*ps^2 + 8*E4*ps - ps*sp*S^2)*A +
	+ E4*S^4 + sp^2*ps^2 +
	- 4*E4*sp*S^2 - 4*E4*ps*S^2 + 8*E4*sp*ps + 16*E4^2 # = 0
B^4 - 2*ps*B^3 - (2*sp*ps + 8*E4 - sp*S^2 - ps^2)*B^2 +
	+ (2*sp*ps^2 + 8*E4*ps - ps*sp*S^2)*B +
	+ E4*S^4 + sp^2*ps^2 +
	- 4*E4*sp*S^2 - 4*E4*ps*S^2 + 8*E4*sp*ps + 16*E4^2 # = 0

# TODO: solve or find simpler alternative;

