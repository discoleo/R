########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: Hetero-Symmetric
### Useful Formulas
###
### draft v.0.1n


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

### Motivation:

# - Characteristic polynomial of Polynomial systems
#   with special types of Symmetries (2xC2 or Ht-4)
#   can be transformed to a polynomial of lower Order;
# - Order poly(S, E4, sp, ps) = 1/4 * Order(x1, x2, x3, x4);
# - Order poly(s1, s2, p1, p2) = 1/2 * Order(x1, x2, x3, x4);

### Aim:
# - compute poly(S, E4, sp, ps);


### Debug
x = sqrt(c(2,3,5,7))
x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];

### Notation:
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

### Note:
# - Eqs: A & B always satisfy the same equation;


##################

### Ea: p*s
### p1*s2 + p2*s1
A = p1*s2 + p2*s1;
B = p1*s1 + p2*s2;
# =>
A^2 - sp*S*A + ps*sp^2 + E4*S^2 - 4*ps*E4 # = 0


#####################

### Ea: p*s^2
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

### Alternatives:
# - may be useful for reductions;
A - S*(p1*s2 + p2*s1) + sp*ps # = 0
B - S*(p1*s1 + p2*s2) + sp*ps # = 0


#####################

### Ea: p^2*s^2
### p1^2*s2^2 + p2^2*s1^2
A = p1^2*s2^2 + p2^2*s1^2;
B = p1^2*s1^2 + p2^2*s2^2;
# =>

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

### Ea: p*s^4
### p1*s2^4 + p2*s1^4
A = p1*s2^4 + p2*s1^4;
B = p1*s1^4 + p2*s2^4;
# =>

### A + B =>
A + B - (p1 + p2)*(s1^4 + s2^4) # = 0
A + B - sp*(S^4 - 4*ps*S^2 + 2*ps^2) # = 0

### A * B =>
A * B - ps^4*(p1^2 + p2^2) - E4*(s1^8 + s2^8) # = 0
A * B - ps^4*(sp^2 - 2*E4) - E4*(S^8 - 8*ps*S^6 + 20*ps^2*S^4 - 16*ps^3*S^2 + 2*ps^4) # = 0
A * B - ps^4*sp^2 - E4*(S^8 - 8*ps*S^6 + 20*ps^2*S^4 - 16*ps^3*S^2) # = 0

# =>
A^2 - A*sp*(S^4 - 4*ps*S^2 + 2*ps^2) +
	+ ps^4*sp^2 + E4*(S^8 - 8*ps*S^6 + 20*ps^2*S^4 - 16*ps^3*S^2) # = 0


#################
#################

#################
### Composite ###
#################

### Formula for:
E21a = x1^2*x2 + x2^2*x3 + x3^2*x4 + x4^2*x1

### Derivation:

A1 = x1^2*x2 + x3^2*x4;
B1 = x1^2*x4 + x3^2*x2;
#
A2 = x2^2*x3 + x4^2*x1;
B2 = x2^2*x3 + x4^2*x1;

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
simplifyPS = function(p) {
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

pT = simplifyPS(pR)
str(pT)

eval.pm(pT, list(E21a=E21a, S=S, E4=E4, s1=s1, s2=s2, p1=p1, p2=p2, sp=sp, ps=ps))

# TODO;


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
### Eq E211a: (the Monster)
# - see Method 2;
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
p1*A + p2*B - E211a # = 0

### Method 1:

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
# similarly: (same formula)
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

# see Method 2: simpler alternative;


#############
### Method 2:
A1 = p1*A + p2*B; # = E211a
B1 = p1*B + p2*A;

### A1 + B1 =>
A1 + B1 - (p1 + p2)*(A + B) # = 0
A1 + B1 - sp*ps # = 0

### A1 * B1
A1 * B1 - p1*p2*(A^2 + B^2) - A*B*(p1^2 + p2^2) # = 0
A1 * B1 - E4*(ps^2 - 2*A*B) - A*B*(sp^2 - 2*E4) # = 0
A1 * B1 - E4*ps^2 - A*B*sp^2 + 4*A*B*E4 # = 0
# =>
A1 * B1 - E4*(ps^2 - 2*p1*s2^2 - 2*p2*s1^2 + 8*E4) +
	- (p1*s2^2 + p2*s1^2 - 4*E4)*(sp^2 - 2*E4) # = 0
A1 * B1 - E4*(ps^2 - 2*p1*s2^2 - 2*p2*s1^2 - 4*sp^2) +
	- (p1*s2^2 + p2*s1^2)*(sp^2 - 2*E4) - 16*E4^2 # = 0
A1 * B1 - E4*(ps^2 - 4*sp^2) +
	- (p1*s2^2 + p2*s1^2)*(sp^2 - 4*E4) - 16*E4^2 # = 0

### Transform 2:
A2 = p1*s2^2 + p2*s1^2
# System: =>
A2^2 - sp*(S^2 - 2*ps)*A2 + E4*(S^4 - 4*ps*S^2) + ps^2*sp^2 # = 0
A1 * B1 - E4*(ps^2 - 4*sp^2) - A2*(sp^2 - 4*E4) - 16*E4^2 # = 0

### Eq E211a:
# - Formula: see at the beginning of the section;


### Additional Formulas
# - necessary for robust computations linked to E211a;

p1 = toPoly.pm("p1*A - p2*A + ps*p2 - E211a");
p2 = toPoly.pm("A*(A - ps) + p1*s2^2 + p2*s1^2 - 4*E4");
pR = solve.pm(p1, p2, "A");
pR = pR$Rez;

pR = replace.fr.pm(pR, data.frame(E4=1, coeff=1), data.frame(p1=1, coeff=1), xn="p2")
pDiv = toPoly.pm("p1^2 - sp*p1 + E4");
pR = div.pm(pR, pDiv, "p1");
pR = pR$Rem; # Remainder
pR = sort.pm(pR, "p1", xn2 = c("ps"))
tmp = toCoeff(pR, "p1", print=TRUE)


#####################
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

