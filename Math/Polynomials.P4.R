########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: P[4]
### v.0.1a


### Solve Method for P[4]

### Experimental Approaches:
# - evaluating solutions for P[4];


### Solver:
solve.P4 = function(b, debug=TRUE) {
	# - assumes b3 = 0; b2 = 0;
	# - assumes all 4 roots are complex;
	a.root = roots(c(64, 0, -16*b[1], -b[2]^2));
	a.root = rootn(a.root, 2)
	if(debug) {
		cat("All intermediary (P[3]) roots:\n")
		print(a.root);
	}
	isReal = (Im(round0(a.root)) == 0)
	a.root = a.root[isReal];
	bc = b[2] / (4*a.root);
	z1sq = a.root^2 + bc;
	z2sq = a.root^2 - bc;
	z1 = rootn(z1sq, 2); z2 = rootn(z2sq, 2);
	return(list(a=a.root, z1=z1, z2=z2));
}

### Case 1:
# - all roots complex;

b = c(1, 1)
x = roots(c(1, 0, 0, rev(b)))
a = Re(x)
x; a

a.root = solve.P4(b)
a.root


### Derivation:
64*a^6 - 16*b[1]*a^2 - b[2]^2


#######################
#######################

####################
### Polynomials: ###
###  of Class 1  ###
####################

### Examples:

### Ex 1:
K = 2
#
k = K^(1/4)
x = k^3 + 2*k^2 - 2*k # the root
#
x^4 - 8*K*(K+4)*x - 16*K + 56*K^2 - K^3


### Ex 2:
K = 3
s = 2
#
k = K^(1/4)
x = s^2*k^3 + 2*s*k^2 - 2*k
#
x^4 - 8*s*K*(K*s^4 + 4)*x - 16*K + 56*s^4*K^2 - s^8*K^3

y = K*s^4 + 4
x^4 - 8*s*K*y*x - 256*K + 64*K*y - K*y^2

