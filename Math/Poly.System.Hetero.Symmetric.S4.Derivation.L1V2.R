########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
###  == Derivation ==
###  Type: L1 V2
###
### draft v.0.1b-reduced


####################
####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")

E2af = function(x, n=2) {
	if(is.null(dim(x))) {
		x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];
	} else {
		stop("TODO");
	}
	if(length(n) == 1) {
		E2a = x1^n*x2 + x2^n*x3 + x3^n*x4 + x4^n*x1;
	} else {
		p = n[2];
		n = n[1];
		E2a = x1^n*x2^p + x2^n*x3^p + x3^n*x4^p + x4^n*x1^p;
	}
	return(E2a);
}
E3af = function(x, n=2) {
	if(is.null(dim(x))) {
		x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];
	} else {
		stop("TODO");
	}
	if(length(n) == 1) {
		E3a = x1^n*x2*x3 + x2^n*x3*x4 + x3^n*x4*x1 + x4^n*x1*x2;
	} else {
		p1 = n[2]; p2 = n[3];
		n  = n[1];
		E3a = x1^n*x2^p1*x3^p2 + x2^n*x3^p1*x4^p2 + x3^n*x4^p1*x1^p2 + x4^n*x1^p1*x2^p2;
	}
	return(E3a);
}


####################
####################

################
### Type V2b ###
################

###############
### Order 2 ###
###############

### V2b:
### x1^2 + b*x2*x3 = R
x1^2 + b*x2*x3 # - R
x2^2 + b*x3*x4 # - R
x3^2 + b*x4*x1 # - R
x4^2 + b*x1*x2 # - R

### Solution:

### Eq 1: Sum(x1*...) =>
S^3 - 3*E2*S + (b + 3)*E3 - R*S # = 0

### Eq 2: Eq H2 + Eq H4 =>
E21 - (E21a + E12a) - (sp*S - E3) # = 0
S*E2 - 3*E3 - (E21a + E12a) - ((E2 - E11a)*S - E3) # = 0
b*(E21a + E12a) + S^3 - 2*E2*S - 4*R*S + 2*b*E3 # = 0
b*E21a + S^3 - 2*E2*S + (b - 4)*R*S - (b^2 - 2*b)*E3 # = 0
(b+1)*S^3 - 2*(b+1)*E2*S +
	+ (b^2 - 2*b - 4)*R*S - (b+1)*(b^2 - 2*b)*E3 # = 0
# Reduction =>
(b + 1)*E2*S - (b + 1)*(b^2 - b + 3)*E3 + (b^2 - b - 3)*R*S # = 0

### Eq 3: Eq H5 + Eq H6
E211a + E112a - ps*sp # = 0
E211a + E112a - E11a*(E2 - E11a) # = 0
b^3*E211a + b^3*E112a + b*(S^2 - 2*E2 - 4*R)*(S^2 + (b - 2)*E2 - 4*R) # = 0
2*b^2*E3*S + (5*b^2 - 4*b - 4)*E2^2 +
	- (5*b^2 - 4*b - 4)*E2*S^2 + 2*(4*b^2 - 8*b - 8)*R*E2 +
	(b^2 - b - 1)*S^4 - (2*b^2 - 8*b - 8)*R*S^2 - 16*b*R^2 - 16*R^2 # = 0
# Reduction =>
(b^3 - 4*b - 3)*E3*S - (5*b^2 - 4*b - 4)*E2^2 +
	+ (2*b^2 - b - 1)*E2*S^2 - 2*(4*b^2 - 8*b - 8)*R*E2 +
	+ (b^2 - 7*b - 7)*R*S^2 + 16*b*R^2 + 16*R^2 # = 0
6*(b + 1)*E3*S + (5*b^2 - 4*b - 4)*E2^2 +
	- 2*b^2*E2*S^2 + 8*(b^2 - 2*b - 2)*R*E2 +
	- (2*b^2 - 8*b - 10)*R*S^2 - 16*b*R^2 - 16*R^2 # = 0


### Helper Eqs:

### Eq H1: E11a
# Sum =>
S^2 - 2*E2 + b*E11a - 4*R # = 0
# b*E11a = - S^2 + 2*E2 + 4*R;

### Eq H2: E21a
# Sum(x2*...) =>
(b+1)*E21a - R*S # = 0
# (b+1)*E21a = R*S;

### Eq H3: E31a
# Sum(x2^2*...) =>
E22a + b*E31a - R*(S^2 - 2*E2) # = 0
E22 - E22b + b*E31a - R*(S^2 - 2*E2) # = 0
E2^2 - 2*S*E3 + 4*E4 - E11a^2 + b*E31a - R*S^2 + 2*R*E2 # = 0
b^3*E31a + b^2*E2^2 + 2*R*b^2*E2 - 2*b^2*S*E3 + 4*b^2*E4 +
	- (S^2 - 2*E2 - 4*R)^2 - b^2*R*S^2 # = 0
b^3*E31a + (b^2 - 4)*E2^2 + 2*R*(b^2 - 8)*E2 + 4*E2*S^2 - 2*b^2*S*E3 + 4*b^2*E4 +
	- S^4 - (b^2 - 8)*R*S^2 - 16*R^2 # = 0

### Eq H4: E12a
# Sum(x4*...) =>
E12a + b*E3 - R*S # = 0
# E12a = - b*E3 + R*S;

### Eq H5: E211a
# Sum(x1^2*...) =>
S4 + b*E211a - R*(S^2 - 2*E2) # = 0
b*E211a + S^4 - R*S^2 + 2*E2^2 - 4*E2*S^2 + 2*R*E2 + 4*E3*S - 4*E4 # = 0

### Eq H6: E112a
# Sum(x4^2*...) =>
E22a + b*E112a - R*(S^2 - 2*E2) # = 0
E22 - (E11a^2 - 2*E4) + b*E112a - R*(S^2 - 2*E2) # = 0
E2^2 - 2*E3*S - E11a^2 + b*E112a - R*S^2 + 2*R*E2 + 4*E4 # = 0
b^3*E112a + b^2*E2^2 - 2*b^2*E3*S +
	- (S^2 - 2*E2 - 4*R)^2 - b^2*R*S^2 + 2*b^2*R*E2 + 4*b^2*E4 # = 0
b^3*E112a + (b^2 - 4)*E2^2 + 4*E2*S^2 - 2*b^2*E3*S - S^4 +
	- (b^2 - 8)*R*S^2 + 2*(b^2 - 8)*R*E2 + 4*b^2*E4 - 16*R^2 # = 0


### Alternatives:

# - Diffs usually do NOT work with low powers;


### Debug:
R = 3; b = -2;
x1 =  0.6180339888 - 0.5257311121i;
x2 = -1.6180339887 + 0.8506508083i;
x3 =  0.6180339888 + 0.5257311121i;
x4 = -1.6180339887 - 0.8506508083i;

x = c(x1,x2,x3,x4)
s1 = x1 + x3; s2 = x2 + x4;
p1 = x1 * x3; p2 = x2 * x4;
sp = p1 + p2; ps = s1 * s2;
S = s1 + s2; E4 = p1 * p2;
E2 = sp + ps;
E3 = p1*s2 + p2*s1;
#
E11a = ps;
E21  = S*E2 - 3*E3;
E21a = E2af(x, n=2);
E12a = E2af(x, n=c(1,2));
E31a = E2af(x, n=3);
#
E22  = E2^2 - 2*S*E3 + 2*E4;
E22a = E2af(x, n=c(2,2));
E22b = (x1*x3)^2 + (x2*x4)^2;
E33a = E2af(x, n=c(3,3));
#
S4 = sum(x^4);
E211a = E3af(x, n=2);
E112a = E3af(x, n=c(1,1,2));

