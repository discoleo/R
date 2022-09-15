########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: P[4]
### v.0.1b


### Solve Method for P[4]

# - proper way to solve a P[4]:
#   based on a C2-decomposition;


####################

### Helper Functions

test.P4 = function(sol, R=NULL) {
	err1 = apply(sol, 1, sum);
	err4 = apply(sol, 1, prod);
	s1 = sol[,1] + sol[,3];
	s2 = sol[,2] + sol[,4];
	p1 = sol[,1] * sol[,3];
	p2 = sol[,2] * sol[,4];
	err2 = p1 + p2 + s1*s2;
	err3 = p1*s2 + p2*s1;
	err = cbind(err1, err2, err3, err4);
	if( ! is.null(R)) {
		err = err - rep(R, each=nrow(err));
	}
	err = round0(err);
	return(err);
}


### Debug:
R = c(-1,2,-3,2)
R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
b = c(1, -R1, R2, -R3, R4)
sol = roots(b)

x1 = sol[1]; x2 = sol[2]; x3 = sol[3]; x4 = sol[4];
S = x1 + x2 + x3 + x4;
s1 = x1 + x3; s2 = x2 + x4;
p1 = x1 * x3; p2 = x2 * x4;
ps = s1 * s2; sp = p1 + p2;
#
E2 = sp + ps;
E3 = p1*s2 + p2*s1;
E4 = p1 * p2;


###########################
###########################

### Method 1
### C2-Symmetry

# - the proper way to solve a P[4];

### Transformed System:
# E2 = R2; E3 = R3;
S - R1 # = 0
sp + ps - E2 # = 0
E3^2 - sp*S*E3 + ps*sp^2 + E4*S^2 - 4*ps*E4 # = 0
E4 - R4 # = 0

### Eq:
ps^3 - 2*E2*ps^2 + (E2^2 + E3*S - 4*E4)*ps + E4*S^2 - E2*E3*S + E3^2 # = 0


### Solver:
solve.P4 = function(R, debug=TRUE, all=FALSE) {
	S = R[1]; E2 = R[2]; E3 = R[3]; E4 = R[4];
	coeff = c(1, - 2*E2, (E2^2 + E3*S - 4*E4), E4*S^2 - E2*E3*S + E3^2);
	ps = roots(coeff);
	sp = E2 - ps;
	if(debug) print(ps);
	FUN = function(id, s, p) roots(c(1, -s[id], p[id]));
	# Step 2:
	len = length(ps);
	s1 = sapply(seq(len), FUN, s = rep(S, len), p = ps);
	s1 = as.vector(s1);
	s2 = S - s1;
	sp = rep(sp, each=2);
	# Step 3: robust;
	p1 = (s1*sp - E3) / (s1 - s2);
	p2 = sp - p1;
	# Step 4:
	len = length(s1);
	x13 = sapply(seq(len), FUN, s=s1, p=p1);
	x13 = t(x13);
	x1 = x13[,1]; x3 = x13[,2];
	#
	x24 = sapply(seq(len), FUN, s=s2, p=p2);
	x24 = t(x24);
	x2 = x24[,1]; x4 = x24[,2];
	#
	sol = cbind(x1, x2, x3, x4);
	if(all) {
		# TODO: all proper permutations;
	}
	invisible(sol);
}

### Examples:

### Ex 1:
R = c(-1,2,-3,2)
sol = solve.P4(R)

test.P4(sol)


### Ex 2:
R = c(3,-1,-3,2)
sol = solve.P4(R)

test.P4(sol)


######################
######################

### [old]

### Experimental Approaches:
# - evaluating solutions for P[4];

### Method 2:
### Old Experiments

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


#######################
#######################

x = roots(c(1,0,0,3,3))
y = 1/x

x^4 + 3*x + 3 # = 0
3*y^4 + 3*y^3 + 1 # = 0

x^4 + 6*x + 6 + 9*(y^4 + 2*y^3 + y^2) # = 0
x + 3*y^3 + 3*y^2 # = 0

x^3 + 3 + 3*y # = 0

p1 = toPoly.pm("x^4 + 3*x + 3")
p2 = toPoly.pm("3*x^4 + 3*x^3 + 1")
p = toPoly.pm("-9*x - (-x^3-3)^3 - 3*(-x^3-3)^2")
p
div.pm(p, p1)

