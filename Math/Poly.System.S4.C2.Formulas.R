########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: Hetero-Symmetric
### Useful Formulas
###
### draft v.0.2i


### Formulas:

# - Formulas & derivations;
# - Useful for Ht & C2 Systems;

# - Applicable for systems described in:
#   Poly.System.S4.C2.Symmetric.R;
#   Poly.System.S4.C2.R;
#   Poly.System.Hetero.Symmetric.S4.R;


### Sections

### Basic:
# A.) Elementary Polynomials: E2, E3
# B.) Helper: p1^k*s2^m + p2^k*s1^m
### Composite:
# C.) Type E2a
# D.) Type E3a
# E.) Symmetric Types


####################

### Helper Functions


source("Poly.System.S4.C2.Helper.R")


####################
####################

### Motivation:

# - Characteristic polynomial of Polynomial systems
#   with special types of Symmetries (2xC2 or Ht-4)
#   can be transformed to a polynomial of lower Order;
# - Order poly(S, E4, sp, ps) = 1/4 * Order(x1, x2, x3, x4);
# - Order poly(s1, s2, p1, p2) = 1/2 * Order(x1, x2, x3, x4);

### Aim:
# - compute poly(S, E4, sp, ps);


### Debug
x = sqrt(c(2,3,5,7));
x[3] = - x[3];
x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];

### Notation:
s1 = x1 + x3; s2 = x2 + x4;
p1 = x1 * x3; p2 = x2 * x4;
sp = p1 + p2; ps = s1 * s2;
S  = s1 + s2;
E4 = p1 * p2;
# used for Reductions:
p1s1 = p1*s1 + p2*s2;

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

### Alternatives:
# E3 = sp*S - p1s1
E3 - sp*S + p1s1;

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

### Note:
# - Eqs: A & B always satisfy the same equation;


##################

### Ea: p*s
### p1*s2 + p2*s1
A = p1*s2 + p2*s1;
B = p1*s1 + p2*s2;
# =>
A^2 - sp*S*A + ps*sp^2 + E4*S^2 - 4*ps*E4 # = 0
B^2 - sp*S*B + ps*sp^2 + E4*S^2 - 4*ps*E4 # = 0
p1s1^2 - sp*S*p1s1 + ps*sp^2 + E4*S^2 - 4*ps*E4 # = 0

### Alternatives:
# A = sp*S - p1s1;
A - sp*S + p1s1;
# B = p1s1;


#####################

### Ea: p*s^2
### p1*s2^2 + p2*s1^2
A = p1*s2^2 + p2*s1^2;
B = p1*s1^2 + p2*s2^2;
# =>
A^2 - sp*(S^2 - 2*ps)*A + E4*(S^4 - 4*ps*S^2) + ps^2*sp^2 # = 0

### Alternatives:
# - may be useful for reductions;
#   A = S*(p1*s2 + p2*s1) - sp*ps;
#   B = S*p1s1 - sp*ps;
#   A = - S*p1s1 + sp*S^2 - sp*ps;
A - S*(p1*s2 + p2*s1) + sp*ps # = 0
B - S*(p1*s1 + p2*s2) + sp*ps # = 0

# other:
sp*B - (p1*s1 + p2*s2)^2 - E4*S^2 + 4*ps*E4 # = 0
# - can be reduced to the alternative;


### Derivation:

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


#####################

### Ea: p^2*s^2
### p1^2*s2^2 + p2^2*s1^2
A = p1^2*s2^2 + p2^2*s1^2;
B = p1^2*s1^2 + p2^2*s2^2;

### Alternatives:
# A = sp*S*(p1*s2 + p2*s1) - ps*sp^2 - E4*S^2 + 2*ps*E4;
# A = - sp*S*p1s1 + sp^2*S^2 - ps*sp^2 - E4*S^2 + 2*ps*E4;
# B = sp*S*p1s1 - ps*sp^2 - E4*S^2 + 2*ps*E4;
A - (p1*s2 + p2*s1)^2 + 2*E4*ps # = 0
A - sp*S*(p1*s2 + p2*s1) + ps*sp^2 + E4*S^2 - 2*ps*E4 # = 0
A + sp*S*p1s1 - sp^2*S^2 + ps*sp^2 + E4*S^2 - 2*ps*E4 # = 0
# B:
B - sp*S*p1s1 + ps*sp^2 + E4*S^2 - 2*ps*E4 # = 0


### Derivation:

### A + B =>
A + B - (p1^2 + p2^2)*(s1^2 + s2^2) # = 0
A + B - (sp^2 - 2*E4)*(S^2 - 2*ps) # = 0

### A * B =>
A * B - ps^2*(p1^4 + p2^4) - E4^2*(s1^4 + s2^4) # = 0
A * B - ps^2*(sp^4 - 4*E4*sp^2 + 2*E4^2) - E4^2*(S^4 - 4*ps*S^2 + 2*ps^2) # = 0
A * B - ps^2*(sp^4 - 4*E4*sp^2) - E4^2*(S^4 - 4*ps*S^2) - 4*E4^2*ps^2 # = 0

# =>
A^2 - A*(sp^2 - 2*E4)*(S^2 - 2*ps) +
	+ ps^2*(sp^4 - 4*E4*sp^2) + E4^2*(S^4 - 4*ps*S^2) + 4*E4^2*ps^2 # = 0


#####################

### Ea: p*s^3
### p1*s2^3 + p2*s1^3
A = p1*s2^3 + p2*s1^3;
B = p1*s1^3 + p2*s2^3;

### Formula:
A^2 - sp*(S^3 - 3*ps*S)*A +
	+ E4*(S^6 - 6*ps*S^4 + 9*ps^2*S^2 - 4*ps^3) + ps^3*sp^2 # = 0

### Alternatives:
# A = (p1*s2 + p2*s1)*(S^2 - ps) - ps*sp*S;
# A = (ps - S^2)*p1s1 + sp*S^3 - 2*ps*sp*S;
# B = (S^2 - ps)*p1s1 - ps*sp*S;
A - (p1*s2 + p2*s1)*(S^2 - ps) + ps*sp*S # = 0
B + (ps - S^2)*p1s1 + ps*sp*S # = 0


### Derivation:

### A + B
A + B - (p1 + p2)*(s1^3 + s2^3) # = 0
A + B - sp*(S^3 - 3*ps*S) # = 0

### A * B
A * B - E4*(s1^6 + s2^6) - ps^3*(p1^2 + p2^2) # = 0
A * B - E4*(S^6 - 6*ps*S^4 + 9*ps^2*S^2 - 2*ps^3) - ps^3*(sp^2 - 2*E4) # = 0
A * B - E4*(S^6 - 6*ps*S^4 + 9*ps^2*S^2 - 4*ps^3) - ps^3*sp^2 # = 0

### Alternatives:
A - (p1*s2 + p2*s1)*(s1^2 + s2^2) + ps*(p1*s1 + p2*s2) # = 0
A - (p1*s2 + p2*s1)*(S^2 - 2*ps) + ps*p1s1 # = 0
# based on: p1*s2 + p2*s1;
A - (p1*s2 + p2*s1)*(S^2 - ps) + ps*sp*S # = 0
# based on: p1*s1 + p2*s2;
A - (S^2 - ps)*(S*sp - p1s1) + ps*sp*S # = 0
A + (S^2 - ps)*p1s1 - sp*S^3 + 2*ps*sp*S # = 0


#####################

### Ea: p*s^4
### p1*s2^4 + p2*s1^4
A = p1*s2^4 + p2*s1^4;
B = p1*s1^4 + p2*s2^4;

### Formula:
A^2 - sp*(S^4 - 4*ps*S^2 + 2*ps^2)*A +
	+ ps^4*sp^2 + E4*(S^8 - 8*ps*S^6 + 20*ps^2*S^4 - 16*ps^3*S^2) # = 0

### Alternatives:
A # =
(2*ps*S - S^3)*p1s1 + sp*S^4 - 3*sp*ps*S^2 + sp*ps^2 # = 0
B # =
- (2*ps*S - S^3)*p1s1 - sp*ps*S^2 + sp*ps^2 # = 0

### Derivation:

### A + B =>
A + B - (p1 + p2)*(s1^4 + s2^4) # = 0
A + B - sp*(S^4 - 4*ps*S^2 + 2*ps^2) # = 0

### A * B =>
A * B - ps^4*(p1^2 + p2^2) - E4*(s1^8 + s2^8) # = 0
A * B - ps^4*(sp^2 - 2*E4) - E4*(S^8 - 8*ps*S^6 + 20*ps^2*S^4 - 16*ps^3*S^2 + 2*ps^4) # = 0
A * B - ps^4*sp^2 - E4*(S^8 - 8*ps*S^6 + 20*ps^2*S^4 - 16*ps^3*S^2) # = 0


### Alternatives:
A - (p1*s2 + p2*s1)*(s1^3 + s2^3) + ps*(p1*s1^2 + p2*s2^2) # = 0
A - (sp*S - p1s1)*(S^3 - 3*ps*S) + ps*(p1*s1^2 + p2*s2^2) # = 0
A + (S^3 - 3*ps*S)*p1s1 - sp*S*(S^3 - 3*ps*S) + ps*(p1*s1^2 + p2*s2^2) # = 0
A + (S^3 - 3*ps*S)*p1s1 - sp*S*(S^3 - 3*ps*S) +
	+ ps*(S*(p1*s1 + p2*s2) - sp*ps) # = 0
A + (S^3 - 2*ps*S)*p1s1 - sp*S^4 + 3*sp*ps*S^2 - sp*ps^2 # = 0


#########################

### Ea: p^2*s^4
### p1^2*s2^4 + p2^2*s1^4
A = p1^2*s2^4 + p2^2*s1^4;
B = p1^2*s1^4 + p2^2*s2^4;

### Formula:
# TODO

### Alternatives:
A # =
(2*sp*ps*S - sp*S^3)*p1s1 +
	+ sp^2*S^4 + sp^2*ps^2 - 3*ps*sp^2*S^2 - 2*ps^2*E4 - E4*S^4 + 4*ps*E4*S^2;
B # =
(sp*S^3 - 2*sp*ps*S)*p1s1 +
	+ ps^2*sp^2 - ps*sp^2*S^2 - 2*ps^2*E4 - E4*S^4 + 4*ps*E4*S^2;


### Derivation:
# - non-standard derivation;

A - (p1*s2^2 + p2*s1^2)^2 + 2*ps^2*E4 # = 0
# Reduction =>
A - (S*p1s1 - sp*S^2 + sp*ps)^2 + 2*ps^2*E4 # = 0
A - S^2*p1s1^2 + (2*sp*S^3 - 2*sp*ps*S)*p1s1 +
	 - sp^2*S^4 - sp^2*ps^2 + 2*sp^2*ps*S^2 + 2*ps^2*E4 # = 0
# Reduction =>
A + (sp*S^3 - 2*sp*ps*S)*p1s1 +
	- sp^2*S^4 - sp^2*ps^2 + 3*ps*sp^2*S^2 + 2*ps^2*E4 + E4*S^4 - 4*ps*E4*S^2 # = 0

### A + B
A + B - (p1^2 + p2^2)*(s1^4 + s2^4) # = 0
A + B - (sp^2 - 2*E4)*(S^4 - 4*ps*S^2 + 2*ps^2) # = 0
A + B - sp^2*S^4 + 4*ps*sp^2*S^2 - 2*ps^2*sp^2 + 2*E4*(S^4 - 4*ps*S^2 + 2*ps^2) # = 0


#########################

### Ea: p^3*s^2
### p1^3*s2^2 + p2^3*s1^2
A = p1^3*s2^2 + p2^3*s1^2;
B = p1^3*s1^2 + p2^3*s2^2;

### Alternatives:
A # =
(sp^2*S - E4*S)*(p1*s2 + p2*s1) - ps*sp^3 - sp*E4*S^2 + 3*sp*ps*E4;
# based on: (p1*s1 + p2*s2):
- (sp^2*S - E4*S)*p1s1 + sp^3*S^2 - ps*sp^3 - 2*sp*E4*S^2 + 3*sp*ps*E4;
# B = (sp^2*S - E4*S)*p1s1 - sp*E4*S^2 - ps*sp^3 + 3*sp*ps*E4;


### Derivation:

### A + B =>
A + B - (p1^3 + p2^3)*(s1^2 + s2^2) # = 0
A + B - (sp^3 - 3*sp*E4)*(S^2 - 2*ps) # = 0
A + B - (sp^3*S^2 - 3*sp*E4*S^2 - 2*ps*sp^3 + 6*sp*ps*E4) # = 0

### A * B =>
A * B - E4^3*(s1^4 + s2^4) - ps^2*(p1^6 + p2^6) # = 0
A * B - E4^3*(S^4 - 4*ps*S^2 + 2*ps^2) +
	- ps^2*(sp^6 - 6*E4*sp^4 + 9*E4^2*sp^2 - 2*E4^3) # = 0

### Alternatives:

### Derivation:
A - sp*(p1^2*s2^2 + p2^2*s1^2) + E4*(p1*s2^2 + p2*s1^2) # = 0
A - sp*(p1^2*s2^2 + p2^2*s1^2) + E4*S*(p1*s2 + p2*s1) - sp*ps*E4 # = 0
A - sp*(sp*S*(p1*s2 + p2*s1) - ps*sp^2 - E4*S^2 + 2*ps*E4) +
	+ E4*S*(p1*s2 + p2*s1) - sp*ps*E4 # = 0
A + (E4*S - sp^2*S)*(p1*s2 + p2*s1) + ps*sp^3 + sp*E4*S^2 - 3*sp*ps*E4 # = 0

### B:
B - (sp^2*S - E4*S)*p1s1 + sp*E4*S^2 + ps*sp^3 - 3*sp*ps*E4 # = 0


#################
#################

#################
### Composite ###
#################

#################
### Type: E2a ###
#################

##############
###  E21a  ###
##############

### Formula for:
E21a = x1^2*x2 + x2^2*x3 + x3^2*x4 + x4^2*x1

### Derivation:

A1 = x1^2*x2 + x3^2*x4;
B1 = x1^2*x4 + x3^2*x2;
#
A2 = x2^2*x3 + x4^2*x1;
B2 = x2^2*x3 + x4^2*x1;


### s2*A1 + s1*A2

### s2*A1 =>
s2*A1 - (x1^2*x2^2 + x3^2*x4^2) - p2*(x1^2 + x3^2) # = 0
s2*A1 - (x1^2*x2^2 + x3^2*x4^2) - p2*(s1^2 - 2*p1) # = 0
s2*A1 - (x1^2*x2^2 + x3^2*x4^2) - p2*s1^2 + 2*E4 # = 0

### s1*A2 =>
s1*A2 - (x2^2*x3^2 + x4^2*x1^2) - p1*(x2^2 + x4^2) # = 0
s1*A2 - (x2^2*x3^2 + x4^2*x1^2) - p1*(s2^2 - 2*p2) # = 0
s1*A2 - (x2^2*x3^2 + x4^2*x1^2) - p1*s2^2 + 2*E4 # = 0

# =>
s2*A1 + s1*A2 - (x1^2+x3^2)*(x2^2+x4^2) - (p1*s2^2 + p2*s1^2) + 4*E4 # = 0
s2*A1 + s1*A2 - (s1^2 - 2*p1)*(s2^2 - 2*p2) - (p1*s2^2 + p2*s1^2) + 4*E4 # = 0
s2*A1 + s1*A2 + (p1*s2^2 + p2*s1^2) - ps^2 # = 0
# Reduction =>
s2*A1 + s1*A2 - S*p1s1 + sp*S^2 - sp*ps - ps^2 # = 0


### Part Eq 1:

### A1 + B1 =>
A1 + B1 - (s1^2 - 2*p1)*s2 # = 0;
# alternative:
A1 + B1 + 2*p1*s2 + ps*s2 - ps*S # = 0;

### A1 * B1 =>
A1 * B1 - p2*(s1^4 - 4*p1*s1^2) - p1^2*s2^2 # = 0
# alternative:
A1 * B1 - p1*s2*sp*S - p2*s1*(S^3 - 2*ps*S) +
	+ p2*(ps*S^2 - ps^2 - sp*ps) + 3*s1*E4*S + sp^2*ps + E4*S^2 - 5*ps*E4 # = 0

# =>
A1^2 + A1*(2*p1*s2 + ps*s2 - ps*S) + p1*s2*sp*S + p2*s1*(S^3 - 2*ps*S) +
	- p2*(ps*S^2 - ps^2 - sp*ps) - 3*s1*E4*S - sp^2*ps - E4*S^2 + 5*ps*E4 # = 0
# similarly:
A2^2 + A2*(2*p2*s1 + ps*s1 - ps*S) + p2*s1*sp*S + p1*s2*(S^3 - 2*ps*S) +
	- p1*(ps*S^2 - ps^2 - sp*ps) - 3*s2*E4*S - sp^2*ps - E4*S^2 + 5*ps*E4 # = 0

### E21a:
A1 + A2 - E21a # = 0
# ? OR ?
E21a^2 - 2*E21a*A1 + A1^2 - A2^2 # = 0

###
pE  = toPoly.pm("A1 + A2 - E21a")
pA1 = toPoly.pm("A1^2 + A1*(2*p1*s2 + ps*s2 - ps*S) + p1*s2*sp*S + p2*s1*(S^3 - 2*ps*S) +
	- p2*(ps*S^2 - ps^2 - sp*ps) - 3*s1*E4*S - sp^2*ps - E4*S^2 + 5*ps*E4")
pA2 = toPoly.pm("A2^2 + A2*(2*p2*s1 + ps*s1 - ps*S) + p2*s1*sp*S + p1*s2*(S^3 - 2*ps*S) +
	- p1*(ps*S^2 - ps^2 - sp*ps) - 3*s2*E4*S - sp^2*ps - E4*S^2 + 5*ps*E4")

pR = solve.pm(pE, pA1, "A1")
pR = solve.pm(pR$Rez, pA2, "A2")
pR = pR$Rez;
pR = orderVars.pm(pR, c("s1","s2","p1","p2","coeff"));
str(pR)


replaceSym = function(p, pow=1, xn=c("s1", "s2"), xr=c("S", "ps"), xflt=c("p1", "p2")) {
	isNot = TRUE;
	if( ! is.null(xflt)) {
		isNot = (p[, xflt[1]] == 0) & (p[, xflt[2]] == 0);
	}
	isS1 = (p[, xn[1]] == pow) & isNot;
	isS2 = (p[, xn[2]] == pow) & isNot;
	pTmp = p[isS2, ];
	# Check existence:
	if(nrow(pTmp) == 0) return(p);
	# Check symmetry:
	pTmp[, xn[1]] = pTmp[, xn[2]];
	pTmp[, xn[2]] = 0;
	pTmp0 = diff.pm(pTmp, p[isS1, ]);
	if(nrow(pTmp0) > 0) {
		warning("Error: ",  xn[1], " & ", xn[2], " do NOT cancel out!");
	} else if(pow == 1) {
		p[isS1, xr[1]] = p[isS1, xr[1]] + p[isS1, xn[1]];
		p[isS1, xn[1]] = 0;
		#
		p = p[ ! isS2, ];
		p = aggregate0.pm(p);
	} else {
		if(pow == 2) {
			pRepl = data.frame(S=c(2,0), p=c(0,1), coeff=c(1,-2));
		} else if(pow == 3) {
			pRepl = data.frame(S=c(3,1), p=c(0,1), coeff=c(1,-3));
		} else {
			warning("Not yet implemented!");
			return(p);
		}
		names(pRepl) = c(xr, "coeff");
		p = p[ ! (isS1 | isS2), ];
		pTmp = replace.pm(pTmp, pRepl, xn[1], pow=pow);
		p = sum.pm(p, pTmp);
	}
	return(p);
}
simplifyPS0 = function(p) {
	p = replace.pm.m(p, c("s1", "s2"), "ps");
	p = replace.pm.m(p, c("p1", "p2"), "E4");
	#
	p = replaceSym(p, pow=1);
	p = replaceSym(p, pow=2);
	p = replaceSym(p, pow=3);
	#
	replaceSymP = function(p, pow) {
		replaceSym(p, pow=pow, xn=c("p1", "p2"), xr=c("sp", "E4"), xflt=c("s1", "s2"));
	}
	p = replaceSymP(p, pow=1);
	p = replaceSymP(p, pow=2);
	p = replaceSymP(p, pow=3);
	return(invisible(p));
}
simplifyPS = function(p, iter=1) {
	p = simplifyPS0(p);
	for(i in seq(iter)) {
		p = replace.pm(p, pS1, "s1", pow=2)
		p = replace.pm(p, pS2, "s2", pow=2)
		p = replace.pm(p, pP1, "p1", pow=2)
		p = replace.pm(p, pP2, "p2", pow=2)
		p = simplifyPS0(p)
	}
	return(invisible(p));
}

pS1 = toPoly.pm("S*s1 - ps")
pS2 = toPoly.pm("S*s2 - ps")
pP1 = toPoly.pm("sp*p1 - E4")
pP2 = toPoly.pm("sp*p2 - E4")

pT = simplifyPS(pR, iter=2)
str(pT)

eval.pm(pT, list(E21a=E21a, S=S, E4=E4, s1=s1, s2=s2, p1=p1, p2=p2, sp=sp, ps=ps))


### pC11 * (s1*p1 + s2*p2)
pC11 = pT[pT$s1 == 1 & pT$p1 == 1, ];
pC11$s1 = 0; pC11$p1 = 0;
pC11 = drop.pm(pC11);
print.pm(pC11, lead=NA)

### pC12 * (s1*p2 + s2*p1)
pC12 = pT[pT$s1 == 1 & pT$p2 == 1, ];
pC12$s1 = 0; pC12$p2 = 0;
pC12 = drop.pm(pC12);
print.pm(pC12, lead=NA)

### pC00
pC00 = pT[pT$s1 == 0 & pT$s2 == 0 & pT$p1 == 0 & pT$p2 == 0, ];
pC00 = drop.pm(pC00);
# print.pm(pC00, lead=NA)

### [old]
pSP12Pow2 = toPoly.pm("(p1*s2 + p2*s1)^2");
pSP12Pow2 = simplifyPS(pSP12Pow2, iter=1);
pSP12Pow2 = orderVars.pm(pSP12Pow2, c("s1","s2","p1","p2","coeff"));

### Formula for E21a:
# - computable, but result is probably too big (364 monomials);
pC11 = toPoly.pm("4*E21a*ps^3 - ps^4*S - 2*E21a*ps^2*S^2 - 6*E4*ps*S^3 + ps^3*S^3 +
	+ 2*E21a*ps^2*sp + 6*E4*ps*S*sp");
pC12 = toPoly.pm("4*E21a^3 + 48*E21a*E4*ps + 6*E21a*ps^3 - 10*E21a^2*ps*S - 40*E4*ps^2*S +
	- 2*ps^4*S - 20*E21a*E4*S^2 + 4*E21a*ps^2*S^2 + 2*E21a^2*S^3 + 48*E4*ps*S^3 +
	+ 2*ps^3*S^3 - 2*E21a*ps*S^4 - 10*E4*S^5 + 10*E21a*ps^2*sp + 6*E21a^2*S*sp +
	+ 30*E4*ps*S*sp - 3*ps^3*S*sp - 18*E21a*ps*S^2*sp - 10*E4*S^3*sp + 15*ps^2*S^3*sp +
	+ 4*E21a*S^4*sp - 8*ps*S^5*sp + S^7*sp - 4*E21a*ps*sp^2 + 10*ps^2*S*sp^2 +
	+ 4*E21a*S^2*sp^2 - 10*ps*S^3*sp^2 + 2*S^5*sp^2 - 2*ps*S*sp^3 + S^3*sp^3");
#
pSP11 = toPoly.pm("sp11^2 - sp*S*sp11 + ps*sp^2 + E4*S^2 - 4*ps*E4");
pSP12 = toPoly.pm("sp12^2 - sp*S*sp12 + ps*sp^2 + E4*S^2 - 4*ps*E4");

###
# pRR = toPoly.pm("C11*sp11 + C12*sp12 + C00");
# alternative: better;
pRR = toPoly.pm("C11*sp11 + C12*(sp*S - sp11) + C00");
#
pRR = solve.pm(pRR, pSP11, "sp11")
pRR = pRR$Rez;

pRR = replace.pm(pRR, pC12, "C12")
table(pRR$C00)
table(pRR$C11)
pRR = replace.pm(pRR, pC11, "C11")
pRR = replace.pm(pRR, pC00, "C00") # requires pC00!

str(pRR)
# but 364 monomials;

# acceptable precision with the smaller E21a;
eval.pm(pRR, list(S=S, E4=E4, ps=ps, sp=sp, E21a=E21a))


### [old]
# - use of pSP12Pow2: sill huge final result;
pSP12Pow2 = toPoly.pm("sp*S*sp12 + 4*E4*ps - E4*S^2 - sp^2*ps");
pRR = replace.pm(pRR, pSP12Pow2, "sp12", pow=2);
pRR = solve.pm(pRR, pSP12, "sp12");
pRR = pRR$Rez;
table(pRR$C12)


#################
#################

################
### E3-Types ###
################

#####################
### Simple Types  ###
### Ea[n+k, n, k] ###
#####################

#############
### E211a ###
#############

### Formula for:
E211a = x1^2*x2*x3 + x2^2*x3*x4 + x3^2*x4*x1 + x4^2*x1*x2;

### Eq E211a: (the Monster)
E211a^4 - 2*ps*sp*E211a^3 +
	+ (sp^3*S^2 - 2*ps*sp^3 + ps^2*sp^2 - 4*sp*E4*S^2 + 2*ps^2*E4 + 8*ps*sp*E4 +
		- 8*sp^2*E4 + 32*E4^2)*E211a^2 +
	+ (2*ps^2*sp^4 + 4*ps*sp^2*E4*S^2- ps*sp^4*S^2 - 2*ps^3*sp*E4 + 8*ps*sp^3*E4 +
		- 8*ps^2*sp^2*E4 - 32*ps*sp*E4^2)*E211a +
	+ sp^4*E4*S^4 - 8*sp^2*E4^2*S^4 + ps^2*sp^6 - 4*sp^5*E4*S^2 - 4*ps*sp^4*E4*S^2 +
		+ ps^2*sp^3*E4*S^2 + 8*ps*sp^5*E4 - 8*ps^2*sp^4*E4 + 16*E4^3*S^4 - 2*ps^3*sp^3*E4 +
		+ 32*sp^3*E4^2*S^2 - 4*ps^2*sp*E4^2*S^2 + 32*ps*sp^2*E4^2*S^2 + ps^4*E4^2 +
		+ 16*sp^4*E4^2 + 8*ps^3*sp*E4^2 - 64*ps*sp^3*E4^2 - 64*ps*E4^3*S^2 - 64*sp*E4^3*S^2 +
		+ 8*ps^2*sp^2*E4^2 + 32*ps^2*E4^3 + 128*ps*sp*E4^3 - 128*sp^2*E4^3 + 256*E4^4;


### Derivation:
p1*(x1*x2 + x3*x4) + p2*(x3*x2 + x1*x4) - E211a # = 0

A = x1*x2 + x3*x4;
B = x1*x4 + x3*x2;
#
p1*A + p2*B - E211a # = 0

### Alternative:
# - used with Method 3 (see "Old Formulas"), but not needed;
p1*A + p2*(ps - A) - E211a # = 0


###########
### Step 1:

### A + B =>
A + B - (x1+x3)*(x2+x4) # = 0
A + B - s1*s2 # = 0
A + B - ps # = 0

### A * B =>
A * B - p1*(x2^2 + x4^2) - p2*(x1^2 + x3^2) # = 0
A * B - p1*(s2^2 - 2*p2) - p2*(s1^2 - 2*p1) # = 0
A * B - (p1*s2^2 + p2*s1^2) + 4*E4 # = 0
# reduction =>
A * B + S*p1s1 - sp*S^2 + sp*ps + 4*E4 # = 0


###########
### Step 2:

A2 = p1*A + p2*B;
B2 = p2*A + p1*B;

### A2 + B2
A2 + B2 - sp*ps # = 0

### A2 * B2
A2 * B2 - E4*(A^2 + B^2) - A*B*(p1^2 + p2^2) # = 0
A2 * B2 - E4*(ps^2 - 2*A*B) - A*B*(sp^2 - 2*E4) # = 0
A2 * B2 - E4*(ps^2 + 2*S*p1s1 - 2*sp*S^2 + 2*sp*ps + 8*E4) +
	+ (S*p1s1 - sp*S^2 + sp*ps + 4*E4)*(sp^2 - 2*E4) # = 0
A2 * B2 + (sp^2*S - 4*E4*S)*p1s1 - 16*E4^2 - E4*(ps^2 - 4*sp^2 - 4*sp*S^2 + 4*sp*ps) +
	- (S^2 - ps)*sp^3 # = 0

# =>
A2^2 - sp*ps*A2 - (sp^2*S - 4*E4*S)*p1s1 +
	+ 16*E4^2 + E4*(ps^2 - 4*sp^2 - 4*sp*S^2 + 4*sp*ps) + (S^2 - ps)*sp^3 # = 0


### Step 3:

pE  = toPoly.pm("A2 - E211a")
pA2 = toPoly.pm("A2^2 - sp*ps*A2 - (sp^2*S - 4*E4*S)*p1s1 +
	+ 16*E4^2 + E4*(ps^2 - 4*sp^2 - 4*sp*S^2 + 4*sp*ps) + (S^2 - ps)*sp^3")
pP1S1 = toPoly.pm("p1s1^2 - sp*S*p1s1 + ps*sp^2 + E4*S^2 - 4*ps*E4")

pR = solve.pm(pE, pA2, "A2")
pR = pR$Rez
pR = solve.pm(pR, pP1S1, "p1s1")
pR = pR$Rez
str(pR)
# 41 Monomials

pR = sort.pm(pR, "E4", xn2="E211a")
invisible(toCoeff(pR, "E211a", print=TRUE))


#####################
#####################

#############
### E321a ###
#############

### Formula for:
E321a = x1^3*x2^2*x3 + x2^3*x3^2*x4 + x3^3*x4^2*x1 + x4^3*x1^2*x2;


### Derivation:

E321a - p1*(x1^2*x2^2 + x3^2*x4^2) - p2*(x2^2*x3^2 + x4^2*x1^2) # = 0

###
A1 = x1*x2 + x3*x4;
A2 = x2*x3 + x4*x1;

E321a - p1*(A1^2 - 2*E4) - p2*(A2^2 - 2*E4) # = 0
E321a - (p1*A1^2 + p2*A2^2) + 2*sp*E4 # = 0

### Sum:
A1 + A2 - ps # = 0

### Prod:
A1 * A2 - p1*(x2^2 + x4^2) - p2*(x1^2 + x3^2) # = 0
A1 * A2 - p1*(s2^2 - 2*p2) - p2*(s1^2 - 2*p1) # = 0
A1 * A2 - (p1*s2^2 + p2*s1^2) + 4*E4 # = 0
# Reduction =>
A1 * A2 + S*p1s1 - sp*S^2 + sp*ps + 4*E4 # = 0

### Eqs:
A1^2 - ps*A1 - S*p1s1 + sp*S^2 - sp*ps - 4*E4 # = 0
A2^2 - ps*A2 - S*p1s1 + sp*S^2 - sp*ps - 4*E4 # = 0

# =>
E321a - p1*(ps*A1 + S*p1s1 - sp*S^2 + sp*ps + 4*E4) +
	- p2*(ps*A2 + S*p1s1 - sp*S^2 + sp*ps + 4*E4) + 2*sp*E4 # = 0
E321a - ps*(p1*A1 + p2*A2) - p1*(S*p1s1 - sp*S^2 + sp*ps + 4*E4) +
	- p2*(S*p1s1 - sp*S^2 + sp*ps + 4*E4) + 2*sp*E4 # = 0
E321a - ps*(p1*A1 + p2*A2) - sp*S*p1s1 + sp^2*S^2 - sp^2*ps - 2*sp*E4 # = 0

############
### Part. 2:
B1 = p1*A1 + p2*A2;
B2 = p2*A1 + p1*A2;

### Sum:
B1 + B2 - sp*ps # = 0

### Prod:
B1 * B2 - E4*(A1^2 + A2^2) - A1*A2*(p1^2 + p2^2) # = 0
B1 * B2 - E4*(ps^2 - 2*A1*A2) - A1*A2*(sp^2 - 2*E4) # = 0
B1 * B2 + (sp^2*S - 4*E4*S)*p1s1 - E4*(ps^2 - 4*sp^2 - 4*sp*S^2 + 4*sp*ps + 16*E4) +
	+ (ps - S^2)*sp^3 # = 0

### Eq:
B1^2 - sp*ps*B1 - (sp^2*S - 4*E4*S)*p1s1 +
	+ E4*(ps^2 - 4*sp^2 - 4*sp*S^2 + 4*sp*ps + 16*E4) + (S^2 - ps)*sp^3 # = 0


#
pE = toPoly.pm("E321a - ps*B1 - sp*S*p1s1 + sp^2*S^2 - sp^2*ps - 2*sp*E4")
pB1 = toPoly.pm("B1^2 - sp*ps*B1 - (sp^2*S - 4*E4*S)*p1s1 +
	+ E4*(ps^2 - 4*sp^2 - 4*sp*S^2 + 4*sp*ps + 16*E4) + (S^2 - ps)*sp^3")
pP1S1 = toPoly.pm("p1s1^2 - sp*S*p1s1 + ps*sp^2 + E4*S^2 - 4*ps*E4")

pR = solve.pm(pE, pB1, "B1")
pR = pR$Rez
pR = solve.pm(pR, pP1S1, "p1s1")
pR = pR$Rez
str(pR)
# 97 Monomials

pR = sort.pm(pR, "E4", xn2="E321a")
invisible(toCoeff(pR, "E321a", print=TRUE))

E321a^4 + (- 8*E4*sp + 2*sp^2*S^2 - 4*sp^2*ps - 2*sp*ps^2)*E321a^3 +
	(24*E4^2*sp^2 + 32*E4^2*ps^2 + 2*E4*sp^2*S^4 - 12*E4*sp^3*S^2 - 8*E4*sp^2*ps*S^2 +
		- 4*E4*sp*ps^2*S^2 + 2*E4*ps^4 + 24*E4*sp^3*ps + 8*E4*sp*ps^3 + 4*E4*sp^2*ps^2 +
		+ sp^4*S^4 - 4*sp^4*ps*S^2 - 2*sp^3*ps^2*S^2 + 6*sp^4*ps^2 + sp^2*ps^4 + 4*sp^3*ps^3)*E321a^2 +
	(- 32*E4^3*sp^3 - 128*E4^3*sp*ps^2 - 8*E4^2*sp^3*S^4 - 16*E4^2*sp*ps^2*S^4 + 24*E4^2*sp^4*S^2 +
		+ 32*E4^2*sp^3*ps*S^2 + 64*E4^2*sp*ps^3*S^2 + 48*E4^2*sp^2*ps^2*S^2 - 48*E4^2*sp^4*ps +
		- 40*E4^2*sp*ps^4 + 8*E4^2*sp^3*ps^2 - 96*E4^2*sp^2*ps^3 + 2*E4*sp^4*S^6 - 4*E4*sp^5*S^4 +
		- 12*E4*sp^4*ps*S^4 + 2*E4*sp^3*ps^2*S^4 + 16*E4*sp^5*ps*S^2 + 16*E4*sp^4*ps^2*S^2 +
		+ 6*E4*sp^2*ps^4*S^2 - 8*E4*sp^3*ps^3*S^2 - 2*E4*sp*ps^6 - 24*E4*sp^5*ps^2 +
		- 12*E4*sp^2*ps^5 - 12*E4*sp^3*ps^4 - sp^5*ps^2*S^4 + 2*sp^6*ps^2*S^2 + 4*sp^5*ps^3*S^2 +
		- 4*sp^6*ps^3 - 2*sp^5*ps^4)*E321a +
	+ 16*E4^4*sp^4 + 128*E4^4*sp^2*ps^2 + 256*E4^4*ps^4 + 8*E4^3*sp^4*S^4 + 16*E4^3*ps^4*S^4 +
	- 16*E4^3*sp^5*S^2 - 64*E4^3*ps^5*S^2 - 32*E4^3*sp^4*ps*S^2 - 64*E4^3*sp*ps^4*S^2 +
	- 80*E4^3*sp^3*ps^2*S^2 + 32*E4^3*ps^6 + 32*E4^3*sp^5*ps + 128*E4^3*sp*ps^5 +
	- 16*E4^3*sp^4*ps^2 - 56*E4^3*sp^2*ps^4 + 160*E4^3*sp^3*ps^3 + E4^2*sp^4*S^8 +
	- 4*E4^2*sp^5*S^6 - 8*E4^2*sp^4*ps*S^6 - 4*E4^2*sp^3*ps^2*S^6 + 4*E4^2*sp^6*S^4 +
	+ 24*E4^2*sp^5*ps*S^4 + 36*E4^2*sp^4*ps^2*S^4 + 24*E4^2*sp^3*ps^3*S^4 +
	- 2*E4^2*sp^2*ps^4*S^4 - 16*E4^2*sp^6*ps*S^2 - 4*E4^2*sp*ps^6*S^2 - 24*E4^2*sp^5*ps^2*S^2 +
	+ 8*E4^2*sp^2*ps^5*S^2 - 80*E4^2*sp^4*ps^3*S^2 - 28*E4^2*sp^3*ps^4*S^2 + E4^2*ps^8 +
	+ 8*E4^2*sp*ps^7 + 24*E4^2*sp^6*ps^2 + 12*E4^2*sp^2*ps^6 - 16*E4^2*sp^5*ps^3 +
	- 8*E4^2*sp^3*ps^5 + 68*E4^2*sp^4*ps^4 + E4*sp^4*ps^4*S^4 - 4*E4*sp^7*ps^2*S^2 +
	- 4*E4*sp^5*ps^4*S^2 - 4*E4*sp^4*ps^5*S^2 + 8*E4*sp^7*ps^3 - 4*E4*sp^6*ps^4 +
	+ 2*E4*sp^4*ps^6 + 8*E4*sp^5*ps^5 + sp^8*ps^4;


### [old]
# A2 = ps - A1;
E321a - ps*(p1*A1 - p2*A1 + ps*p2) - sp*S*p1s1 + sp^2*S^2 - sp^2*ps - 2*sp*E4 # = 0
#
dp = p1 - p2;
E321a - ps*dp*A1 - ps^2*p2 - sp*S*p1s1 + sp^2*S^2 - sp^2*ps - 2*sp*E4 # = 0


#####################

#############
### E312a ###
#############

### Formula for:
E312a = x1^3*x2*x3^2 + x2^3*x3*x4^2 + x3^3*x4*x1^2 + x4^3*x1*x2^2;


### Derivation:

E312a - p1^2*(x1*x2 + x3*x4) - p2^2*(x2*x3 + x4*x1) # = 0

###
A1 = x1*x2 + x3*x4;
A2 = x2*x3 + x4*x1;

E312a - p1^2*A1 - p2^2*A2 # = 0

### Sum:
A1 + A2 - ps # = 0

### Prod:
A1 * A2 - p1*(x2^2 + x4^2) - p2*(x1^2 + x3^2) # = 0
A1 * A2 - p1*(s2^2 - 2*p2) - p2*(s1^2 - 2*p1) # = 0
A1 * A2 - (p1*s2^2 + p2*s1^2) + 4*E4 # = 0
# Reduction =>
A1 * A2 + S*p1s1 - sp*S^2 + sp*ps + 4*E4 # = 0

### Eqs:
A1^2 - ps*A1 - S*p1s1 + sp*S^2 - sp*ps - 4*E4 # = 0
A2^2 - ps*A2 - S*p1s1 + sp*S^2 - sp*ps - 4*E4 # = 0

############
### Part. 2:
B1 = p1^2*A1 + p2^2*A2;
B2 = p2^2*A1 + p1^2*A2;

### Sum:
B1 + B2 - ps*(p1^2 + p2^2) # = 0
B1 + B2 - ps*(sp^2 - 2*E4) # = 0

### Prod:
B1 * B2 - E4^2*(A1^2 + A2^2) - A1*A2*(p1^4 + p2^4) # = 0
B1 * B2 - E4^2*(ps^2 - 2*A1*A2) - A1*A2*(sp^4 - 4*sp^2*E4 + 2*E4^2) # = 0
B1 * B2 - E4^2*(ps^2 + 2*(S*p1s1 - sp*S^2 + sp*ps + 4*E4)) +
	+ (S*p1s1 - sp*S^2 + sp*ps + 4*E4)*(sp^4 - 4*sp^2*E4 + 2*E4^2) # = 0
B1 * B2 + (sp^4*S - 4*sp^2*E4*S)*p1s1 - E4^2*(ps^2 + 16*sp^2) +
	+ 4*sp^2*E4*(sp*S^2 + sp^2 - sp*ps) - S^2*sp^5 + ps*sp^5 # = 0

### Eq:
B1^2 - ps*(sp^2 - 2*E4)*B1 - (sp^4*S - 4*sp^2*E4*S)*p1s1 +
	+ E4^2*(ps^2 + 16*sp^2) +
	- 4*sp^2*E4*(sp*S^2 + sp^2 - sp*ps) + S^2*sp^5 - ps*sp^5 # = 0


### Step 3:
pE = toPoly.pm("E312a - B1")
pB1 = toPoly.pm("B1^2 - ps*(sp^2 - 2*E4)*B1 - (sp^4*S - 4*sp^2*E4*S)*p1s1 +
	+ E4^2*(ps^2 + 16*sp^2) +
	- 4*sp^2*E4*(sp*S^2 + sp^2 - sp*ps) + S^2*sp^5 - ps*sp^5")
pP1S1 = toPoly.pm("p1s1^2 - sp*S*p1s1 + ps*sp^2 + E4*S^2 - 4*ps*E4")

pR = solve.pm(pE, pB1, "B1")
pR = pR$Rez
pR = solve.pm(pR, pP1S1, "p1s1")
pR = pR$Rez
str(pR)
# 48 Monomials

pR = sort.pm(pR, "E4", xn2="E312a")
invisible(toCoeff(pR, "E312a", print=TRUE))

### Eq:
E312a^4 + (4*ps*E4 - 2*sp^2*ps)*E312a^3 +
	(32*sp^2*E4^2 + 6*ps^2*E4^2 - 4*sp^3*S^2*E4 - 8*sp^4*E4 + 8*sp^3*ps*E4 +
		- 4*sp^2*ps^2*E4 + sp^5*S^2 - 2*sp^5*ps + sp^4*ps^2)*E312a^2 +
	(64*sp^2*ps*E4^3 + 4*ps^3*E4^3 - 8*sp^3*ps*S^2*E4^2 - 48*sp^4*ps*E4^2 +
		+ 16*sp^3*ps^2*E4^2 - 2*sp^2*ps^3*E4^2 + 6*sp^5*ps*S^2*E4 + 8*sp^6*ps*E4 +
		- 12*sp^5*ps^2*E4 - sp^7*ps*S^2 + 2*sp^7*ps^2)*E312a +
	256*sp^4*E4^4 + 32*sp^2*ps^2*E4^4 + ps^4*E4^4 + 16*sp^4*S^4*E4^3 - 64*sp^5*S^2*E4^3 +
		- 64*sp^4*ps*S^2*E4^3 - 4*sp^3*ps^2*S^2*E4^3 - 128*sp^6*E4^3 + 128*sp^5*ps*E4^3 +
		- 8*sp^4*ps^2*E4^3 + 8*sp^3*ps^3*E4^3 - 8*sp^6*S^4*E4^2 + 32*sp^7*S^2*E4^2 +
		+ 32*sp^6*ps*S^2*E4^2 + sp^5*ps^2*S^2*E4^2 + 16*sp^8*E4^2 - 64*sp^7*ps*E4^2 +
		+ 16*sp^6*ps^2*E4^2 - 2*sp^5*ps^3*E4^2 + sp^8*S^4*E4 - 4*sp^9*S^2*E4 +
		- 4*sp^8*ps*S^2*E4 + 8*sp^9*ps*E4 - 8*sp^8*ps^2*E4 + sp^10*ps^2 # = 0


#####################
#####################

#############
### E422a ###
#############

### Formula for:
E422a = x1^4*x2^2*x3^2 + x2^4*x3^2*x4^2 + x3^4*x4^2*x1^2 + x4^4*x1^2*x2^2;

### Eq E422a:



### Derivation:
p1^2*((x1*x2)^2 + (x3*x4)^2) + p2^2*((x3*x2)^2 + (x1*x4)^2) - E422a # = 0

### Method 1:
# - without initial reduction;

###########
### Step 1:
A1 = (x1*x2)^2 + (x3*x4)^2;
B1 = (x1*x4)^2 + (x3*x2)^2;
#
p1^2*A1 + p2^2*B1 - E422a # = 0

### A1 + B1
A1 + B1 - (x1^2 + x3^2)*(x2^2 + x4^2) # = 0
A1 + B1 - (s1^2 - 2*p1)*(s2^2 - 2*p2) # = 0
A1 + B1 + 2*(p1*s2^2 + p2*s1^2) - ps^2 - 4*E4 # = 0
# Reduction =>
A1 + B1 - 2*S*p1s1 + 2*sp*S^2 - 2*sp*ps - ps^2 - 4*E4 # = 0

### A1 * B1
A1 * B1 - p1^2*(x2^4 + x4^4) - p2^2*(x1^4 + x3^4) # = 0
A1 * B1 - p1^2*(s2^4 - 4*p2*s2^2 + 2*p2^2) - p2^2*(s1^4 - 4*p1*s1^2 + 2*p1^2) # = 0
A1 * B1 - (p1^2*s2^4 + p2^2*s1^4) + 4*E4*(p1*s2^2 + p2*s1^2) - 4*E4^2 # = 0
# Reduction =>
A1 * B1 - (p1^2*s2^4 + p2^2*s1^4) - 4*E4*S*p1s1 + 4*E4*(sp*S^2 - sp*ps) - 4*E4^2 # = 0
A1 * B1 - (4*E4*S + 2*sp*ps*S - sp*S^3)*p1s1 +
	+ E4*(S^4 + 4*sp*S^2 - 4*ps*S^2 + 2*ps^2 - 4*sp*ps) - 4*E4^2 +
	- sp^2*S^4 - sp^2*ps^2 + 3*ps*sp^2*S^2 # = 0


###########
### Step 2:
A2 = p1^2*A1 + p2^2*B1;
B2 = p2^2*A1 + p1^2*B1;

### A2 + B2 =>
A2 + B2 - (p1^2 + p2^2)*(A1 + B1) # = 0
A2 + B2 - (sp^2 - 2*E4)*(2*S*p1s1 - 2*sp*S^2 + 2*sp*ps + ps^2 + 4*E4) # = 0
A2 + B2 - 2*(sp^2*S - 2*E4*S)*p1s1 +
	- (sp^2 - 2*E4)*(- 2*sp*S^2 + 2*sp*ps + ps^2 + 4*E4) # = 0

### A2 * B2 =>
A2 * B2 - E4^2*(A1^2 + B1^2) - A1*B1*(p1^4 + p2^4) # = 0
A2 * B2 - E4^2*((A1 + B1)^2 - 2*A1*B1) - A1*B1*(sp^4 - 4*E4*sp^2 + 2*E4^2) # = 0
# =>
A2 * B2 - E4^2*((2*S*p1s1 - 2*sp*S^2 + 2*sp*ps + ps^2 + 4*E4)^2 - 2*A1*B1) +
	- A1*B1*(sp^4 - 4*E4*sp^2 + 2*E4^2) # = 0

###########
### Step 3:
A2 - E422a # = 0
# =>
E422a^2 - (A2 + B2)*E422a + A2*B2 # = 0

###
pE = toPoly.pm("E422a^2 - SA2B2*E422a + A2B2")
pSA2B2 = toPoly.pm("- 2*(sp^2*S - 2*E4*S)*p1s1 +
	- (sp^2 - 2*E4)*(- 2*sp*S^2 + 2*sp*ps + ps^2 + 4*E4)"); # actually: - (A2 + B2);
pSA2B2$coeff = - pSA2B2$coeff;
pA2B2 = toPoly.pm("- E4^2*((2*S*p1s1 - 2*sp*S^2 + 2*sp*ps + ps^2 + 4*E4)^2 - 2*A1B1) +
	- A1B1*(sp^4 - 4*E4*sp^2 + 2*E4^2)"); # actually: - (A2 * B2);
pA2B2$coeff = - pA2B2$coeff;
pA1B1 = toPoly.pm("(4*E4*S + 2*sp*ps*S - sp*S^3)*p1s1 +
	- E4*(S^4 + 4*sp*S^2 - 4*ps*S^2 + 2*ps^2 - 4*sp*ps) + 4*E4^2 +
	+ sp^2*S^4 + sp^2*ps^2 - 3*ps*sp^2*S^2");
pR = replace.pm(pA2B2, pA1B1, "A1B1")
#
pE = replace.pm(pE, pSA2B2, "SA2B2");
pE = replace.pm(pE, pR, "A2B2");

pP1S1 = toPoly.pm("p1s1^2 - sp*S*p1s1 + ps*sp^2 + E4*S^2 - 4*ps*E4")

pR = solve.pm(pE, pP1S1, "p1s1")
pR = pR$Rez
str(pR)
# 213 Monomials: is probably the Minimal polynomial;

pR = sort.pm(pR, "E422a", xn2="E4")
# invisible(toCoeff(pR, "E422a", print=TRUE))

eval.pm(pR, list(E4=E4, sp=sp, ps=ps, S=S, E422a=E422a))


#####################
#####################

##############
###  Sums  ###
##############

###################
### E21a + E12a ###
###################

E21a = x1^2*x2 + x2^2*x3 + x3^2*x4 + x4^2*x1;
E12a = x1*x2^2 + x2*x3^2 + x3*x4^2 + x4*x1^2;

### Formula for:
### E21a + E12a
(E21a + E12a)^2 + 2*S*(sp - ps)*(E21a + E12a) +
	- 2*sp*ps*S^2 + ps^2*S^2 + 4*sp^2*ps + 4*S^2*E4 - 16*ps*E4 # = 0

### Alternatives:
E21a + E12a - 2*p1s1 - ps*S + 2*sp*S # = 0

### Derivation:
E21a + E12a - s1*(x2^2 + x4^2) - s2*(x1^2 + x3^2) # = 0
E21a + E12a - s1*(s2^2 - 2*p2) - s2*(s1^2 - 2*p1) # = 0
E21a + E12a - ps*S + 2*(p1*s2 + p2*s1) # = 0
# =>
E21a + E12a - 2*p1s1 - ps*S + 2*sp*S # = 0

#
pEs = toPoly.pm("E21a + E12a - 2*p1s1 - ps*S + 2*sp*S")
pP1S1 = toPoly.pm("p1s1^2 - sp*S*p1s1 + ps*sp^2 + E4*S^2 - 4*ps*E4")
pR = solve.pm(pEs, pP1S1, "p1s1")
pR = pR$Rez;
pR = sort.pm(pR, c("E21a", "E12a"))
print.pm(pR, lead=NA)


#####################

#####################
### E211a + E112a ###
#####################

E211a = x1^2*x2*x3 + x2^2*x3*x4 + x3^2*x4*x1 + x4^2*x1*x2;
E112a = x1*x2*x3^2 + x2*x3*x4^2 + x3*x4*x1^2 + x4*x1*x2^2;

### Sum
E211a + E112a - sp*ps # = 0

### Diff
E211a - E112a - (p1 - p2)*(x1 - x3)*(x2 - x4) # = 0
(E211a - E112a)^2 - (sp^2 - 4*E4)*(s1^2 - 4*p1)*(s2^2 - 4*p2) # = 0
# =>
(E211a - E112a)^2 - (sp^2 - 4*E4)*(ps^2 - 4*(p1*s2^2 + p2*s1^2) + 16*E4) # = 0

### Prod
A = x1*x2 + x3*x4;
B = x1*x4 + x3*x2;
pAB = A*B;
pAB - p1*s2^2 - p2*s1^2 + 4*E4 # = 0

E211a * E112a - E4*ps^2 - A*B*sp^2 + 4*A*B*E4 # = 0
E211a * E112a - E4*ps^2 - pAB*(sp^2 - 4*E4) # = 0
# =>
(E211a - E112a)^2 - (sp^2 - 4*E4)*(ps^2 - 4*pAB) # = 0
(E211a - E112a)^2 - (sp^2 - 4*E4)*ps^2 + 4*(E211a * E112a - E4*ps^2) # = 0
# redundant
(E211a - E112a)^2 - sp^2*ps^2 + 4*E211a * E112a # = 0


### Prod: E211a * E112a

# Helper:
E22a = (x1*x2)^2 + (x2*x3)^2 + (x3*x4)^2 + (x4*x1)^2;
E301a = x1^3*x3 + x2^3*x4 + x3^3*x1 + x4^3*x2;
E323a = x1^3*x2^2*x3^3 + x2^3*x3^2*x4^3 + x3^3*x4^2*x1^3 + x4^3*x1^2*x2^3;
p1s1  = p1*s1 + p2*s2;

### Prod:
E211a * E112a - E323a - 4*E4^2 - E301a*E4 - E22a*E4 # = 0
E211a * E112a - (4*E4*S - sp^2*S)*p1s1 +
	- 4*sp*ps*E4 + 4*sp^2*E4 - ps^2*E4 + 4*sp*E4*S^2 + ps*sp^3 - sp^3*S^2 - 16*E4^2 # = 0

# from E211a:
E211a^2 - ps*sp*E211a + (4*E4*S - sp^2*S)*p1s1 +
	+ sp^3*S^2 + ps^2*E4 - 4*sp^2*E4 - 4*sp*E4*S^2 + 16*E4^2 + 4*sp*ps*E4 - ps*sp^3 # = 0
# => unfortunately redundant!
E211a * E112a + E211a^2 - ps*sp*E211a # = 0


#####################
#####################

#####################
### E311a + E113a ###
#####################

E311a = x1^3*x2*x3 + x2^3*x3*x4 + x3^3*x4*x1 + x4^3*x1*x2;
E113a = x1*x2*x3^3 + x2*x3*x4^3 + x3*x4*x1^3 + x4*x1*x2^3;

E311a + E113a - p1*s2*(s1^2 - 2*p1) - p2*s1*(s2^2 - 2*p2) # = 0
E311a + E113a - p1*s2*s1^2 - p2*s1*s2^2 + 2*p1^2*s2 + 2*p2^2*s1 # = 0
E311a + E113a - ps*(p1*s1 + p2*s2) + 2*p1^2*s2 + 2*p2^2*s1 # = 0

# TODO


#################
#################

#################
### Symmetric ###
#################

#################
### Type E2a  ###
### E[n,n]    ###
#################

############
### E22a ###
############

### Formula for:
E22a = (x1*x2)^2 + (x2*x3)^2 + (x3*x4)^2 + (x4*x1)^2;

### Derivation:
E22a - (x1^2 + x3^2)*(x2^2 + x4^2) # = 0
E22a - (s1^2 - 2*p1)*(s2^2 - 2*p2) # = 0
E22a + 2*(p2*s1^2 + p1*s2^2) - ps^2 - 4*E4 # = 0

A = p2*s1^2 + p1*s2^2;
# =>
E22a + 2*A - ps^2 - 4*E4 # = 0
A^2 - sp*(S^2 - 2*ps)*A + E4*(S^4 - 4*ps*S^2) + ps^2*sp^2 # = 0
# =>
E22a^2 - (8*E4 - 2*sp*S^2 + 2*ps^2 + 4*ps*sp)*E22a +
	+ 16*E4^2 + 4*E4*(S^4 - 4*ps*S^2 - 2*sp*S^2 + 2*ps^2 + 4*ps*sp) +
	- 2*ps^2*sp*S^2 + ps^4 + 4*ps^3*sp + 4*ps^2*sp^2 # = 0

# alternative: reduction;
E22a + 2*(p2*s1^2 + p1*s2^2) - ps^2 - 4*E4 # = 0
# =>
E22a + 2*S*(p1*s2 + p2*s1) - 2*sp*ps - ps^2 - 4*E4 # = 0
E22a - 2*S*(p1*s1 + p2*s2) + 2*sp*S^2 - 2*sp*ps - ps^2 - 4*E4 # = 0


#######################
#######################

#######################
### Type Quasi-E2a  ###
### E[k,0,n,0]      ###
#######################

#############
### E301a ###
#############

### Formula for:
E301a = x1^3*x3 + x2^3*x4 + x3^3*x1 + x4^3*x2;

### Derivation:
E301a - p1*(x1^2 + x3^2) - p2*(x2^2 + x4^2) # = 0
E301a - p1*(s1^2 - 2*p1) - p2*(s2^2 - 2*p2) # = 0
E301a - (p1*s1^2 + p2*s2^2) + 2*(p1^2 + p2^2) # = 0
E301a - (p1*s1^2 + p2*s2^2) + 2*(sp^2 - 2*E4) # = 0
# reduction =>
E301a - S*(p1*s1 + p2*s2) + sp*ps + 2*sp^2 - 4*E4 # = 0

# TODO


#####################
#####################

#####################
### Type E3a sym  ###
### E[k,n,k]      ###
#####################

#############
### E121a ###
#############

### Formula for:
E121a = x1*x2^2*x3 + x2*x3^2*x4 + x3*x4^2*x1 + x4*x1^2*x2;

### Alternatives:
E121a # =
- S*p1s1 + sp*S^2 - sp*ps - 4*E4;


### Derivation:
E121a - p1*(x2^2 + x4^2) - p2*(x1^2 + x3^2) # = 0
E121a - p1*(s2^2 - 2*p2) - p2*(s1^2 - 2*p1) # = 0
E121a - (p1*s2^2 + p2*s1^2) + 4*E4 # = 0
# Reduction =>
E121a + S*p1s1 - sp*S^2 + sp*ps + 4*E4 # = 0


#############

#############
### E323a ###
#############

### Formula for:
E323a = x1^3*x2^2*x3^3 + x2^3*x3^2*x4^3 + x3^3*x4^2*x1^3 + x4^3*x1^2*x2^3;

### Alternatives:
E323a # =
(E4*S - sp^2*S)*p1s1 +
	+ 3*sp*ps*E4 - 2*E4*sp^2 + 4*E4^2 - 2*sp*E4*S^2 - ps*sp^3 + sp^3*S^2;

### Derivation:
E323a - p1^3*(x2^2 + x4^2) - p2^3*(x1^2 + x3^2) # = 0
E323a - p1^3*(s2^2 - 2*p2) - p2^3*(s1^2 - 2*p1) # = 0
E323a - p1^3*s2^2 - p2^3*s1^2 + 2*E4*(p1^2 + p2^2) # = 0
E323a - p1^3*s2^2 - p2^3*s1^2 + 2*E4*(sp^2 - 2*E4) # = 0
# Reduction =>
E323a - p1*(sp*p1 - E4)*s2^2 - p2*(sp*p2 - E4)*s1^2 + 2*E4*sp^2 - 4*E4^2 # = 0
E323a - sp*p1^2*s2^2 - sp*p2^2*s1^2 + E4*(p1*s2^2 + p2*s1^2) + 2*E4*sp^2 - 4*E4^2 # = 0
E323a - sp*(p1^2*s2^2 + p2^2*s1^2) + E4*(S*(p1*s2 + p2*s1) - sp*ps) + 2*E4*sp^2 - 4*E4^2 # = 0
E323a - sp*(sp*S*(p1*s2 + p2*s1) - ps*sp^2 - E4*S^2 + 2*ps*E4) +
	+ E4*S*(sp*S - p1s1) - sp*ps*E4 + 2*E4*sp^2 - 4*E4^2 # = 0
E323a - sp^2*S*(p1*s2 + p2*s1) - E4*S*p1s1 +
	- 3*sp*ps*E4 + 2*E4*sp^2 - 4*E4^2 + 2*sp*E4*S^2 + ps*sp^3 # = 0
E323a + (sp^2*S - E4*S)*p1s1 +
	- 3*sp*ps*E4 + 2*E4*sp^2 - 4*E4^2 + 2*sp*E4*S^2 + ps*sp^3 - sp^3*S^2 # = 0


