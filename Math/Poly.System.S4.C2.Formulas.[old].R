########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: Hetero-Symmetric
### Formulas: Old Derivations
###
### draft v.0.1a


### Formulas:
# - Old derivations;


####################

### Helper Functions


source("Poly.System.S4.C2.Helper.R")


####################
####################

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
A1 + B1 - (x1^2 + x3^2)*(x2 + x4) # = 0;
A1 + B1 - (s1^2 - 2*p1)*s2 # = 0;
# alternative:
A1 + B1 - (S*s1 - ps - 2*p1)*s2 # = 0;
A1 + B1 + 2*p1*s2 + ps*s2 - ps*S # = 0;


### A1 * B1 =>
A1 * B1 - x2*x4*(x1^4 + x3^4) - (x1*x3)^2*(x2^2 + x4^2) # = 0
A1 * B1 - p2*(s1^4 - 4*p1*s1^2 + 2*p1^2) - p1^2*(s2^2 - 2*p2) # = 0
A1 * B1 - p2*(s1^4 - 4*p1*s1^2) - p1^2*s2^2 # = 0
# alternative:
A1 * B1 - p2*((S*s1 - ps)^2 - 4*p1*s1^2) - (sp*p1 - E4)*(S*s2 - ps) # = 0
A1 * B1 - p2*(S*s1 - ps)^2 + 4*E4*s1^2 - (sp*p1 - E4)*(S*s2 - ps) # = 0
A1 * B1 - p2*(S^2*(S*s1 - ps) + ps^2 - 2*ps*S*s1) + 4*E4*s1^2 - (sp*p1 - E4)*(S*s2 - ps) # = 0
A1 * B1 - p2*(s1*S^3 - 2*s1*ps*S - ps*S^2 + ps^2) + 4*E4*(S*s1 - ps) +
	- p1*s2*sp*S + p1*sp*ps + s2*E4*S - ps*E4 # = 0
A1 * B1 - p1*s2*sp*S - p2*s1*(S^3 - 2*ps*S) +
	+ p2*(ps*S^2 - ps^2) + p1*sp*ps + 4*s1*E4*S + s2*E4*S - 5*ps*E4 # = 0
A1 * B1 - p1*s2*sp*S - p2*s1*(S^3 - 2*ps*S) +
	+ p2*(ps*S^2 - ps^2 - sp*ps) + 3*s1*E4*S + sp^2*ps + E4*S^2 - 5*ps*E4 # = 0


### [old]
### Eq 1 =>
A1*((s1^2 - 2*p1)*s2 - A1) - p2*(s1^4 - 4*p1*s1^2) - p1^2*s2^2 # = 0
A1^2 - A1*(s1^2 - 2*p1)*s2 + p2*(s1^4 - 4*p1*s1^2) + p1^2*s2^2 # = 0
A1^2 - A1*(s1^2 - 2*p1)*s2 + p2*s1^4 - 4*p1*p2*s1^2 + p1^2*s2^2 # = 0
# =>
A1^2 - A1*(s1^2 - 2*p1)*s2 + p2*s1^4 - 4*E4*s1^2 + p1^2*s2^2 # = 0


### Part Eq 2:
# - similarly:
A2 = x2^2*x3 + x4^2*x1;
B2 = x2^2*x3 + x4^2*x1;

# =>
A2^2 - A2*(s2^2 - 2*p2)*s1 + p1*(s2^4 - 4*p2*s2^2) + p2^2*s1^2 # = 0
# =>
A2^2 - A2*(s2^2 - 2*p2)*s1 + p1*s2^4 - 4*E4*s2^2 + p2^2*s1^2 # = 0


### Initial Eq:
# A1 + A2 = E21a;

A1 + A2 - E21a # = 0
A1^2 - A1*(s1^2 - 2*p1)*s2 + p2*s1^4 - 4*E4*s1^2 + p1^2*s2^2 # = 0
A2^2 - A2*(s2^2 - 2*p2)*s1 + p1*s2^4 - 4*E4*s2^2 + p2^2*s1^2 # = 0

# TODO:
# - derive final equation;

###
p1 = toPoly.pm("A1 + A2 - E21a")
p2 = toPoly.pm("A1^2 - A1*(s1^2 - 2*p1)*s2 + p2*s1^4 - 4*E4*s1^2 + p1^2*s2^2")
p3 = toPoly.pm("A2^2 - A2*(s2^2 - 2*p2)*s1 + p1*s2^4 - 4*E4*s2^2 + p2^2*s1^2")

pR = solve.lpm(p1, p2, p3, xn=c("A1", "A2"))
str(pR)
# 99 Monomials: seems a lot!

pR = pR[[2]]$Rez;
pR$ps = 0;
pR = pR[ , c("E21a", "s1", "s2", "ps", "p1", "p2", "E4", "coeff")];

### E4:
powE4 = sapply(seq(nrow(pR)), function(id) {
	min(pR[id, c("p1", "p2")]);
})
# isE4 = powE4 > 0;
pR$E4 = pR$E4 + powE4;
pR$p1 = pR$p1 - powE4;
pR$p2 = pR$p2 - powE4;

### ps:
powPs = sapply(seq(nrow(pR)), function(id) {
	min(pR[id, c("s1", "s2")]);
})
pR$ps = pR$ps + powPs;
pR$s1 = pR$s1 - powPs;
pR$s2 = pR$s2 - powPs;

# TODO:
# - absolute monster;

### Method 2:

pA1 = toPoly.pm("A1^2 - A1*Q1 - R1")
pA2 = toPoly.pm("A2^2 - A2*Q2 - R2")
pE = toPoly.pm("(E21a^2 - A1*Q1 - A2*Q2 - SR)^2 - 4*A1*A2*Q1*Q2 - 4*A1*Q1*R2 - 4*A2*Q2*R1 - 4*PR")

pR = pE;
pR = solve.pm(pR, pA1, "A1")$Rez
pR = solve.pm(pR, pA2, "A2")$Rez
str(pR)
# 1009 monomials!

pTmpR = pR[pR$R1 > 0 & pR$R2 == 0, ];
nrow(pTmpR)
# 2*242 monomials;
# but non-trivially R-symmetric;

### R:
R1 = p2*s1^4 - 4*E4*s1^2 + p1^2*s2^2;
R2 = p1*s2^4 - 4*E4*s2^2 + p2^2*s1^2;
SR = R1 + R2;

A14 = p2*s1^4 + p1*s2^4;
A22 = p1^2*s2^2 + p2^2*s1^2;
SR + 4*E4*(S^2 - 2*ps) - A14 - A22

# A14:
A14^2 - A14*sp*(S^4 - 4*ps*S^2 + 2*ps^2) +
	+ ps^4*sp^2 + E4*(S^8 - 8*ps*S^6 + 20*ps^2*S^4 - 16*ps^3*S^2) # = 0
