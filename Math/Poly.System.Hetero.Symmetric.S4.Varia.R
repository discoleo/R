########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric: Various
###
### draft v.0.1c


### Various Hetero-Symmetric Systems


####################

### Helper Functions

source("Polynomials.Helper.R")

test.S4Ht.Simple = function(sol, R) {
	x1 = sol[,1]; x2 = sol[,2];
	x3 = sol[,3]; x4 = sol[,4];
	#
	err1 = x1 + x2 + x3 + x4;
	err2 = x1*x2 + x2*x3 + x3*x4 + x4*x1;
	err3 = x1*x2*x3 + x2*x3*x4 + x3*x4*x1 + x4*x1*x2;
	err4 = x1*x2*x3*x4;
	err = rbind(err1, err2, err3, err4);
	err = round0(err);
	return(err);
}

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
x1*x3*(x2 + x4) # = 0
# R1, R4 != 0 => Contradiction!
# => NO solutions!


### G.2.) (x1 + x3)*(x2 + x4) = R2
# (x1 + x3)^2 = - R2;

### Eq 1 =>
# x2 + x4 = - (x1 + x3);

### Eq 3 =>
x1*x3*(x2 + x4) + x2*x4*(x1 + x3) # = 0
# & Eq 1 =>
(x1 + x3)*(x1*x3 - x2*x4) # = 0
# Eq 2: (x1 + x3) != 0 =>
# x1*x3 = x2*x4

### Eq 4 =>
(x1*x3)^2 - R4 # = 0

### Derived System:
# (x1*x3)^2 = R4
# (x1 + x3)^2 = -R2

### Solver:
solve.S4Ht.P1Simple = function(R) {
	R2 = R[1]; R4 = R[2];
	x13 = rootn(R4, 2);
	xs  = rootn(-R2, 2);
	# all variants:
	x13 = c(x13, - x13); xs = c(xs, -xs);
	x13 = c(x13, x13); xs = rep(xs, each=2);
	xd  = rootn(xs^2 - 4*x13, 2);
	# Sol: x1 & x3
	x1 = (xs + xd)/2; x3 = xs - x1;
	# x2 & x4:
	xs  = R2 / xs;
	x24 = R4 / x13;
	xd  = rootn(xs^2 - 4*x24, 2);
	# Sol: x2 & x4
	x2 = (xs + xd)/2; x4 = xs - x2;
	sol = cbind(x1=x1, x2=x2, x3=x3, x4=x4);
	# TODO: all valid permutations;
	return(sol);
}

### Examples:

R = c(-2, 3)
sol = solve.S4Ht.P1Simple(R);
test.S4Ht.Simple(sol);

###
R = c(-2, -1)
sol = solve.S4Ht.P1Simple(R);
test.S4Ht.Simple(sol);

