########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S4:
### Mixed Type with Resonances
###
### draft v.0.1b-more


### Heterogeneous Symmetric
### Polynomial Systems: 4 Variables
### Mixed: Hetero + Symmetric

### Example:
# x1^p*x2^n*x3^n + x2^p*x3^n*x4^n + x3^p*x4^n*x1^n + x4^p*x1^n*x2^n = R1
# x1^k + x2^k + x^3^k + x^4^k = R2
# ...


###############

###############
### History ###
###############


### draft v.0.1b:
# - [started] classic roots;
# - more entangled roots; [v.0.1b-more]
### draft v.0.1a:
# - Proof of concept:
#   Rotations Entangled with Roots of Unity;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R

### other functions
# ...


#####################
#####################

##################
### Rotations: ###
### x1^n*x2*x3 ###
##################

### Order: 3+1+1
# x1^3*x2*x3 + x2^3*x3*x4 + x3^3*x4*x1 + x4^3*x1*x2 = R1
# x1^5 + x2^5 + x3^5 + x4^5 = R2
# x1^10 + x2^10 + x3^10 + x4^10 = R3
# (x1*x2)^5 + (x2*x3)^5 + (x3*x4)^5 + (x4*x1)^5 = R4


### Classic Solution:

# a toy model for computing robust roots;
solve.byx1.S4M5.classic = function(x1, R) {
	A = roots(c(1, -R[2], R[4]));
	A = rep(A, each=length(x1));
	B = R[2] - A;
	x1 = rep(x1, each=2);
	x3p5  = A - x1^5;
	x24p5 = (R[2]^2 - R[3] - 2*R[4]) / 2 - x1^5*x3p5;
	# non-robust
	# x2^10 - B*x2^5 + x24p5 = 0
	len = length(x3p5);
	x2p5 = sapply(seq(len), function(id) roots(c(1, - B[id], x24p5[id])));
	x2p5 = as.vector(x2p5);
	# repeat x1
	x1 = rep(x1, each=2);
	B = rep(B, each=2);
	x3p5 = rep(x3p5, each=2);
	x4p5 = B - x2p5;
	print(cbind(x1, x2p5, x3p5, x4p5))
	# TODO: ???
	m = unity(5, all=TRUE);
	x = lapply(seq_along(x1), function(id) {
		x2 = rootn(x2p5[id], 5) * m;
		x3 = rootn(x3p5[id], 5) * m;
		x4 = rootn(x4p5[id], 5) * m;
		expand.grid(x1[id], x2, x3, x4);
	})
	sol = do.call(rbind, x);
	names(sol) = paste0("x", seq(4));
	return(sol);
}
test.R1 = function(x) {
	sum(x[1]^3*x[2]*x[3], x[2]^3*x[3]*x[4], x[3]^3*x[4]*x[1], x[4]^3*x[1]*x[2]);
}

### Test:
# x1 & R: see below;
x = solve.byx1.S4M5.classic(x1, R);
R1r = apply(x, 1, test.R1);
R1r = round(R1r, 2);
# only 5 root-tuples are real!
x[R1r == R[1], ]


### Derivation:

### let:
A = x1^5 + x3^5;
B = x2^5 + x4^5;
# Eq 2 =>
A + B - R2 # = 0
# Eq 4 =>
A*B - R4 # = 0
# =>
A^2 - R2*A + R4 # = 0

### (Eq 2)^2 - Eq 3 - 2*Eq 4 =>
2*(x1*x3)^5 + 2*(x2*x4)^5 - R2^2 + R3 + 2*R4 # = 0
# (x1*x3)^5 + (x2*x4)^5 = (R2^2 - R3 - 2*R4) / 2;
# if (x1*x3)^5 is known =>
x2^10 - B*x2^5 + (R2^2 - R3 - 2*R4) / 2 - (x1*x3)^5 # = 0

### Eq 2 =>
# x1^5 + x3^5 = R2 - (x2^5 + x4^5) # sq =>
x1^10 + x3^10 + 2*(x1*x3)^5 - (x2^10 + x4^10 + 2*(x2*x4)^5 - 2*R2*(x2^5 + x4^5) + R2^2) # = 0
x2^10 + x4^10 + 2*(x2*x4)^5 - (x1^10 + x3^10 + 2*(x1*x3)^5 - 2*R2*(x1^5 + x3^5) + R2^2) # = 0
### Sum => [redundant]


### Test
x1^3*x2*x3 + x2^3*x3*x4 + x3^3*x4*x1 + x4^3*x1*x2 # - R1
x1^5 + x2^5 + x3^5 + x4^5 # - R2
x1^10 + x2^10 + x3^10 + x4^10 # - R3
(x1*x2)^5 + (x2*x3)^5 + (x3*x4)^5 + (x4*x1)^5 # - R4


### Debug:
R = c(1, -2, -1, 3);
R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
x1 =  0.0451568676 + 0.1440038706i;
x2 =  0.9916018740 + 0.5052387172i;
x3 = -0.1370521599 - 1.1076829770i;
x4 =  0.7725898004 + 0.1223347600i;
x = c(x1, x2, x3, x4);

m = unity(5, all=TRUE);
x = sapply(seq_along(m), function(id) x*m[id]);
x = t(x)
x1 = x[,1]; x2 = x[,2]; x3 = x[,3]; x4 = x[,4];

x0 = x[1,];
sol0 = x0 * m[2]^c(1,0,2,3);
sol1 = x0 * m[2]^c(1,2,0,4);
sol2 = x0 * m[2]^c(3,2,4,0);
sol3 = x0 * m[2]^c(1,4,3,0);
sol4 = x0 * m[2]^c(1,3,4,2);
x = rbind(x, sol0, sol1, sol2, sol3, sol4);
x = rbind(x, x[,c(2,3,4,1)], x[,c(3,4,1,2)], x[,c(4,1,2,3)]);
rownames(x) = NULL
x1 = x[,1]; x2 = x[,2]; x3 = x[,3]; x4 = x[,4];



##################
##################
##################

##################
### Resonances ###
##################

### Order: 3+1+1
### 4 Variables
# i = 4; p = c(3,1,1)
# p = Divisors(3^3*5^2);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(15);
# - the computed ones do NOT work;
# - but there are quasi-non-trivial solutions for p = 5;
### Ex: p = 5 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^1, m^2, m^0, m^4);
# (1,2,0), (2,0,4), (0,4,1), (4,1,2)
# (1,3,4), (3,4,2), (4,2,1), (2,1,3)
# (1,4,3), (4,3,0), (3,0,1), (0,1,4)
# (1,0,2), (0,2,3), (2,3,1), (3,1,0)
# (3,2,4), (2,4,0), (4,0,3), (0,3,2)
### Ex: p = 15 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^1, m^11, m^1, m^11);
# (1,11,1), (11,1,11), (1,11,1), (11,1,11)
