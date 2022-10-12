########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S2
### Asymmetric Derived from Ht/Symmetric
### Transform: Fractions
###
### draft v.0.1a


#######################

### Helper Functions

source("Polynomials.Helper.R")


#######################


########################
### Hetero-Symmetric ###
########################

### Division-Type

### Base System: P[2]
# x^2 + b2*x*y + b1*y = R
# y^2 + b2*x*y + b1*x = R

### Transform:
# x => x;
# y => y/x;

### Eq 1:
x^3 + b2*x*y - R*x + b1*y # = 0

### Eq 2:
b1*x^3 - R*x^2 + y^2 + b2*x^2*y # = 0
y^2 + b2*x^2*y - R*x^2 - b1*b2*x*y + b1*R*x - b1^2*y # = 0

### System:
x^3 + b2*x*y - R*x + b1*y # = 0
y^2 + b2*x^2*y - R*x^2 - b1*b2*x*y + b1*R*x - b1^2*y # = 0


### Solver:

solve.S2AsymFr.P2 = function(R, b, debug=TRUE) {
	b1 = b[1]; b2 = b[2];
	S  = b1;
	if(debug) print(S);
	xy = - (S^2 + b1*S - 2*R) / (2*b2 - 2);
	d = sqrt(S^2 - 4*xy + 0i);
	x = (S + d)/2; y = S - x;
	x = c(x, y); y = x[c(2, 1)];
	# Tr:
	y = x*y;
	#
	sol = cbind(x, y);
	return(sol);
}
test.S2AsymFr.P2 = function(sol, R, b) {
	x = sol[,1]; y = sol[,2];
	b1 = b[1]; b2 = b[2];
	err1 = x^3 + b2*x*y - R*x + b1*y;
	err2 = y^2 + b2*x^2*y - R*x^2 - b1*b2*x*y + b1*R*x - b1^2*y;
	err = rbind(err1, err2);
	err = round0(err);
	return(err);
}

### Example:

R = -2; b = c(3,-2);
sol = solve.S2AsymFr.P2(R, b)

test.S2AsymFr.P2(sol, R=R, b=b)


### Ex 2:
R = 3; b = c(-1, 4);
sol = solve.S2AsymFr.P2(R, b)

test.S2AsymFr.P2(sol, R=R, b=b)


##############

### Variant:
### Transform:
# x => x / (y + f1);
# y => y / (x + f2);

# TODO

