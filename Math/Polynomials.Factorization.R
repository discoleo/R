########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials:
### Factorization
###
### draft v.0.1a


### Factorization of Polynomials
# - some experimental approaches;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;


##########################


###############
### Order 4 ###
###############

# - alternative to classic solution for the P4;

### x^4 + B2*x^2 + B1*x + B0
(x^2 + b1*x + c1)*(x^2 - b1*x + c2)
### =>
x^4 + (c1 + c2 - b1^2)*x^2 + b1*(c1 - c2)*x + c1*c2

### Solution:
c1 + c2 - b1^2 - B2 # = 0
b1*(c1 - c2) - B1 # = 0
c1*c2 - B0 # = 0

### Eq2^2 =>
b1^2*(c1 - c2)^2 - B1^2 # = 0
### S = c1 + c2 =>
b1^2*(S^2 - 4*c1*c2) - B1^2 # = 0
b1^2*(S^2 - 4*B0) - B1^2 # = 0

### Eq 1 =>
(c1 + c2)*(S^2 - 4*B0) - b1^2*(S^2 - 4*B0) - B2*(S^2 - 4*B0) # = 0
S*(S^2 - 4*B0) - B1^2 - B2*(S^2 - 4*B0) # = 0
S^3 - B2*S^2 - 4*B0*S + 4*B0*B2 - B1^2 # = 0

### Solver:
factorize.P4 = function(b, invert=FALSE) {
	if(length(b) > 3) {
		if(b[4] != 0) stop("P4 must be in reduced form! Not yet implemented.")
	}
	coeff = c(1, - b[3], - 4*b[1], 4*b[1]*b[3] - b[2]^2)
	S = roots(coeff)
	#
	c.d = sqrt(S^2 - 4*b[1] + 0i)
	#
	c.d = c(c.d, -c.d)
	S = c(S, S)
	isZero = round0(c.d) == 0;
	b1 = b[2] / c.d;
	if(any(isZero)) print("Warning: Diff(c) == 0;")
	b1[isZero] = sqrt(S - b[3] + 0i); # TODO: "+/-" vs correct variant;
	c1 = (S + c.d)/2;
	c2 = (S - c.d)/2;
	if(invert) b1 = - b1;
	# TODO: select correct variant;
	# - it seems all are correct;
	p1 = sapply(seq_along(S), function(id) c(1, -b1[id], c1[id]))
	p2 = sapply(seq_along(S), function(id) c(1,  b1[id], c2[id]))
	return(list(p1=p1, p2=p2))
}

### Examples:

b = c(1, 1, 3)
p = factorize.P4(b)
p

polymul(p$p1[,1], p$p2[,1])


### Ex 2:
b = c(1, -1, 0)
p = factorize.P4(b)
p

polymul(p$p1[,1], p$p2[,1])

