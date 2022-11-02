########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S3:
### Mixed Type: Dual / Multiple E2a Eqs
###
### draft v.0.1b-sol


### Heterogeneous Symmetric
### Polynomial Systems: 3 Variables
### Mixed: Hetero + Symmetric
### with 2 or 3 Rotations

### Example:
# x^n1*y^p1 + y^n1*z^p1 + z^n1*x^p1 = R1
# x^n2*y^p2 + y^n2*z^p2 + z^n2*x^p2 = R2

### Eq 3:
# S = R3 or E2 = R3 or E3 = R3;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


# library(polynom)
# library(pracma)

# the functions are in the file:
# Polynomials.Helper.R


### Other Functions

test.Ht3Dual = function(x, y, z, R, n=2, p=1, b=0, type) {
	# Types Eq 3: E3, E3*S, E3/S;
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	x.sum = x + y + z;
	e3f = if(missing(type)) 1 else if(match("M3", type) > 0) x.sum else if(match("D3", type) > 0) 1/x.sum;
	if( ! missing(type)) print(paste0("Type: ", type))
	### Test
	err1 = x^p*y^n + y^p*z^n + z^p*x^n + b[1]*x.sum
	err2 = x^p*z^n + y^p*x^n + z^p*y^n + if(length(b) < 2) 0 else b[2]*x.sum
	err3 = x*y*z * e3f + if(length(b) < 3) 0 else b[3]*x.sum
	err = rbind(err1, err2, err3)
	if( ! missing(R)) {
		err = err - rep(R, each=length(x))
	}
	err = round0(err)
	return(err)
}

###########################
###########################

###########################
### Mixed Systems       ###
### Type 2: 2 Rotations ###
###########################

### Dual Eq: p = 1
### x*y^n + y*z^n + z*x^n = R1
### x*z^n + y*x^n + z*y^n = R2

###############
### Order 2 ###
###############

x*y^2 + y*z^2 + z*x^2 - R1 # = 0
x*z^2 + y*x^2 + z*y^2 - R2 # = 0
x*y*z - R3 # = 0

### Extensions:
### Extension M3:
# x*y*z*(x+y+z) = R3
### Extension D3:
# x*y*z / (x+y+z) = R3

### Solution:
### Eq 1 + Eq 2 =>
E2*S - 3*E3 - R1 - R2 # = 0

### see Section "Simple System":
R3*S^6 - (R1*R2 + 6*R1*R3 + 6*R2*R3 + 9*R3^2)*S^3 + (R1 + R2 + 3*R3)^3 # = 0

### Extension M3:
R3*S^8 - R1*R2*S^6 - 6*R3*(R1 + R2)*S^5 - 9*R3^2*S^4 + (R1 + R2)^3*S^3 +
	+ 9*R3*(R1 + R2)^2*S^2 + 27*R3^2*(R1 + R2)*S + 27*R3^3 # = 0

### Solver:
solve.Ht3Dual = function(R, b=0, type, sort=FALSE, debug=TRUE) {
	# TODO: R1 == R2
	if(missing(type)) {
		type = 1;
		RS3 = (R[1] + R[2] + 3*R[3]);
		coeff = c(R[3], 0, 0, - (R[1]*R[2] + 6*(R[1]+R[2])*R[3] + 9*R[3]^2),
			0, 0, RS3^3)
		if(b[1] != 0) {
			coeff = coeff + c(0, 0, b[1]*(R[2] + 6*R[3]),
				-b[1]^3, 3*b[1]^2*RS3, -3*b[1]*RS3^2, 0);
		}
		if(length(b) > 1) {
			if(b[2] != 0) {
				coeff = coeff + c(0, 0, b[2]*(R[1] + 6*R[3]),
					-b[2]^3, 3*b[2]^2*RS3, -3*b[2]*RS3^2, 0);
			}
			if(b[1] != 0 && b[2] != 0) {
				b12 = b[1]*b[2]; b12s = b[1]+b[2];
				coeff = coeff + c(0, -b12, 0,
					-3*b12*b12s, 6*b12*RS3, 0, 0);
			}
		}
	} else if(match("M3", type) > 0) {
		type = 3;
		R12 = R[1] + R[2]
		coeff = c(R[3], 0, - R[1]*R[2], - 6*R[3]*(R12), - 9*R[3]^2, (R12)^3,
			9*R[3]*(R12)^2, 27*R[3]^2*(R12), 27*R[3]^3)
	}
	S = roots(coeff);
	if(debug) print(S);
	len = length(S);
	b2 = if(length(b) > 1) b[2] else 0;
	b3 = if(length(b) > 2) b[3] else 0;
	R1 = R[1] - b[1]*S; R2 = R[2] - b2*S; R3 = R[3] - b3*S;
	E3 = if(type == 1) R3 else if(type == 3) R3 / S;
	E2 = (R1 + R2 + 3*E3);
	if(any(S == 0)) {
		print("Warning: Div by 0!")
		E2.zero = - rootn(R1^2 + 9*R3^2 + 3*R1*R3, 3);
		E2 = ifelse(S == 0, E2.zero, E2 / S);
		# E3 cannot be 0!
		E3 = ifelse(S == 0, R3, E3); # TODO: check if correct (by Type)!
		x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	} else {
		E2 = E2 / S;
		x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	}
	S  = matrix(S,  ncol=len, nrow=3, byrow=T);
	E3 = matrix(E3, ncol=len, nrow=3, byrow=T);
	R1 = rep(R1, each=3); R2 = rep(R2, each=3);
	yz = E3/x;
	yz.s = S - x;
	### robust:
	yz.d = (R2 - R1) / (x^2 + yz - x*yz.s);
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z));
	if(sort) sol = sort.sol(sol);
	return(sol);
}

### Examples:

R = c(1, 3, -1);
sol = solve.Ht3Dual(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.Ht3Dual(x, y, z)

round0.p(poly.calc(x))

err = 1 - 37*x^6 + 76*x^9 - 37*x^12 + x^18 # trivial & symmetric;
round0(err)

### Classic Polynomial:
# - degenerate generalized-symmetric P18;
R1 = R[1]; R2 = R[2]; R3 = R[3];
R3*x^18 - (R1*R2 - 3*R3^2)*x^15 +
	+ (R1^3 + R2^3 - 5*R1*R2*R3 + 6*R3^3)*x^12 +
	- (R1^2*R2^2 - 2*R1^3*R3 - 2*R2^3*R3 + 6*R1*R2*R3^2 - 7*R3^4)*x^9 +
	+ R3^2*(R1^3 + R2^3 - 5*R1*R2*R3 + 6*R3^3)*x^6 +
	- R3^4*(R1*R2 - 3*R3^2)*x^3 + R3^7


#########
### Ex 2:
R = c(-1, 2, 1);
sol = solve.Ht3Dual(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.Ht3Dual(x, y, z, n=2)

### degenerate!
round0.p(poly.calc(x[1:6]^3))
x = x^3
err = 1 + 5*x + 23*x^2 + 29*x^3 + 23*x^4 + 5*x^5 + x^6 # symmetric;
round0(err)


#########
### Ex 3: Ext A1
R = c(0, -1, 1);
b = 1;
sol = solve.Ht3Dual(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.Ht3Dual(x, y, z, b=b, n=2)

###
round0.p(poly.calc(x))
err = 1 + 3*x + 3*x^2 + 4*x^3 + 9*x^4 + 10*x^5 + 10*x^6 + 13*x^7 + 11*x^8 + 12*x^9 +
	+ 12*x^10 + 7*x^11 + 7*x^12 + 5*x^13 + 4*x^14 + 2*x^15 + 2*x^16 + x^18;
round0(err)


#########
### Ex 4: Ext A2
R = c(0, 0, 1);
b = c(1,-1);
sol = solve.Ht3Dual(R, b=b, sort=TRUE)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.Ht3Dual(x, y, z, b=b, n=2)

###
round0.p(poly.calc(x))
err = 1 - x^2 + 3*x^3 + 8*x^5 + 6*x^6 - 9*x^7 + 17*x^8 + 7*x^9 + 19*x^11 + 6*x^12 +
	+ 10*x^14 + 3*x^15 + x^17 + x^18;
round0(err)


#################
### Extension M3:

### Ex 1:
R = c(1, 3, -1);
sol = solve.ht3Dual(R, type="M3")
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.Ht3Dual(x, y, z, n=2, type="M3")

round0.p(poly.calc(x)) * 27

err = 1 + 20*x + 166*x^2 + 736*x^3 + 1862*x^4 + 2596*x^5 + 1305*x^6 - 2488*x^7 - 7403*x^8 +
	- 8428*x^9 - 2887*x^10 + 3848*x^11 + 9086*x^12 + 7032*x^13 - 81*x^14 - 1160*x^15 - 3402*x^16 +
	- 108*x^17 - 180*x^18 - 432*x^19 - 27*x^20 + 81*x^22 + 27*x^24;
round0(err)


### Ex 2:
R = c(0, 3, -3);
sol = solve.ht3Dual(R, type="M3")
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3Dual(x, y, z, n=2, type="M3")

round0.p(poly.calc(x))

err = 27 + 135*x + 270*x^2 + 270*x^3 + 108*x^4 - 135*x^5 - 405*x^6 - 522*x^7 - 315*x^8 - 18*x^9 +
	+ 135*x^10 + 231*x^11 + 222*x^12 + 72*x^13 - 8*x^15 - 27*x^16 - 9*x^17 + 3*x^18 - 9*x^19 +
	- 3*x^20 + x^24
round0(err)


### TODO: R1 == R2;


### Classic Polynomial:
# TODO


#######################
#######################

################
### Variant: ###
### E2-Type  ###
################

############
### Order 2: n = 2
x*y^2 + y*z^2 + z*x^2 - R1 # = 0
x*z^2 + y*x^2 + z*y^2 - R2 # = 0
x*y + x*z + y*z - R3 # = 0

### Extensions:
### Extension M3:
# (x*y + x*z + y*z)*(x+y+z) = R3
### Extension D3:
# (x*y + x*z + y*z) / (x+y+z) = R3

### Solution:
### Eq 1 + Eq 2 =.
x*y^2 + y*z^2 + z*x^2 + x*z^2 + y*x^2 + z*y^2 - R1 - R2 # = 0
E2*S - 3*E3 - R1 - R2 # = 0
# 3*E3 = E2*S - R1 - R2 # = 0
# 3*E3 = R3*S - R1 - R2 # = 0

### see Section "Simple System":
E3*S^3 - (R1+6*E3)*E2*S + R1^2 + E2^3 + 9*E3^2 + 3*R1*E3 # = 0
E3*S^3 - (R1+6*E3)*R3*S + R1^2 + R3^3 + 9*E3^2 + 3*R1*E3 # = 0
E3*S^3 - (R1+2*(R3*S - R1 - R2))*R3*S + R1^2 + R3^3 + (R3*S - R1 - R2)^2 + R1*(R3*S - R1 - R2) # = 0
E3*S^3 - (2*R3*S - 2*R1 - 2*R2)*R3*S + R3^3 + (R3*S - R1 - R2)^2 - R1*R2 # = 0
3*E3*S^3 - 3*R3^2*S^2 + 3*R1^2 + 3*R2^2 + 3*R3^3 + 3*R1*R2 # = 0
(R3*S - R1 - R2)*S^3 - 3*R3^2*S^2 + 3*R1^2 + 3*R2^2 + 3*R3^3 + 3*R1*R2 # = 0
R3*S^4 - (R1 + R2)*S^3 - 3*R3^2*S^2 + 3*R1^2 + 3*R2^2 + 3*R3^3 + 3*R1*R2 # = 0

### Extension A2: power 1;
# (x*y + x*z + y*z) + b1*(x+y+z) = R3;


### Solver:
solve.Ht3Dual.E2 = function(R, b=0, type, tol=1E-5, debug=TRUE) {
	# TODO: analyse special cases: R3 == 0, R1 == R2;
	if(missing(type) || type == 1) {
		type = 1;
		coeff = c(R[3], - (R[1] + R[2]), - 3*R[3]^2, 0, 3*(R[1]^2 + R[2]^2 + R[3]^3 + R[1]*R[2]));
		if(all(b == 0)) {
		} else {
			if(b[1] != 0) {
				coeff = coeff + c(b[1], 0, 3*b[1]^2, -6*b[1]*R[1] - 3*b[1]*R[2], 0);
			}
			if(length(b) > 1 && b[2] != 0) {
				coeff = coeff + c(b[2], 0, 3*b[2]^2, -6*b[2]*R[2] - 3*b[2]*R[1], 0);
			}
			if(length(b) > 1 && b[1] != 0 && b[2] != 0) {
				coeff = coeff + c(0, 0, 3*b[1]*b[2], 0, 0);
			}
			if(length(b) > 2 && b[3] != 0) {
				coeff = c(0, coeff) +
					c(- b[3], -3*b[3]^2, 6*b[3]*R[3] - 3*b[3]^3, 9*b[3]^2*R[3], -9*b[3]*R[3]^2, 0)
			}
			# TODO: power 2;
		}
	} else if(match("M3", type) > 0) {
		type = 3;
		R12 = R[1] + R[2]
		coeff = c() # TODO
	}
	S = roots(coeff)
	if(debug) print(S);
	b2  = if(length(b) > 1) b[2] else 0;
	b31 = if(length(b) > 2) b[3] else 0; # Ext A[E2]: power 1;
	b32 = if(length(b) > 3) b[4] else 0; # Ext A[E2]: power 2;
	R1 = R[1] - b[1]*S;
	R2 = R[2] - b2*S; R3 = R[3] - b31*S - b32*S^2;
	E2 = R3; # Ext A;
	E3 = (E2*S - R1 - R2) / 3;
	len = length(S);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	S  = matrix(S,  ncol=len, nrow=3, byrow=T)
	E2 = matrix(E2, ncol=len, nrow=3, byrow=T);
	R1 = rep(R1, each=3); R2 = rep(R2, each=3);
	# E3 = matrix(E3, ncol=len, nrow=3, byrow=T) # not used;
	yz.s = S - x
	yz = E2 - x*yz.s
	### robust: R[1] == R[2]
	y = ifelse(R1 == R2, {
		isZero = round0(x*yz.s - x^2 - yz, tol=tol) == 0
		y = ifelse(isZero, NA, yz.s/2)
		if(any(isZero)) {
			print("Special case: == 0")
			print(paste0("Zeros: ", sum(isZero)))
			# TODO: still NOT perfect;
			# - only half of (y, z) are correct!
			# - but even these are numerically unstable!
			# - also fails for R[3] != 0;
			# - remaining: z = x; y = yz.s - x;
			y[isZero] = - yz.s[isZero]
		}
		y;
		}, {
		yz.d = (R2 - R1) / (x^2 + yz - x*yz.s);
		y = (yz.s + yz.d) / 2;
		y;
	});
	z = yz.s - y;
	cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
}

### Examples:

R = c(0, 3, -3);
sol = solve.Ht3Dual.E2(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x*y^2 + y*z^2 + z*x^2 # - R[1] # = 0
x*z^2 + y*x^2 + z*y^2 # - R[2] # = 0
x*y + x*z + y*z # - R[3] # = 0

round0.p(poly.calc(x))

err = 27 - 57*x + 19*x^2 + 37*x^3 - 9*x^4 - 33*x^5 + 3*x^6 + 18*x^7 + 3*x^8 +
	- 6*x^9 - 3*x^10 + x^11 + x^12
round0(err)


#########
### Ex 2:
R = c(1, 2, 0);
sol = solve.Ht3Dual.E2(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x*y^2 + y*z^2 + z*x^2 # - R[1] # = 0
x*z^2 + y*x^2 + z*y^2 # - R[2] # = 0
x*y + x*z + y*z # - R[3] # = 0

round0.p(poly.calc(x)) # trivial polynomial for R3 == 0;

err = 1 + 3*x^3 - 4*x^6 + x^9
round0(err)


#########
### Ex 3:
# special Test!
R = c(1, 1, 0);
sol = solve.Ht3Dual.E2(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x*y^2 + y*z^2 + z*x^2 # - R[1] # = 0
x*z^2 + y*x^2 + z*y^2 # - R[2] # = 0
x*y + x*z + y*z # - R[3] # = 0

# trivial polynomial for R3 == 0;
round0.p(poly.calc(x)) * 27

err = 1 + 3*x^3 - 4*x^6 + x^9
round0(err)


### Extensions:

### A Eq 1&2:
R = c(1, -1, 2);
b = c(1,1,0)
sol = solve.Ht3Dual.E2(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
S = (x + y + z);
x*y^2 + y*z^2 + z*x^2 + b[1]*S # - R[1] # = 0
x*z^2 + y*x^2 + z*y^2 + b[2]*S # - R[2] # = 0
x*y + x*z + y*z + b[3]*S # - R[3] # = 0

### Classic Polynomial: P[12] => P[6]
round0.p(poly.calc(x) * 6)
err = 128 + 352*x^2 + 448*x^4 + 334*x^6 + 154.5*x^8 + 43.5*x^10 + 6*x^12
round0(err)


########
### AE2: power 1;
R = c(1, 2, 3);
b = c(0,0,1)
sol = solve.Ht3Dual.E2(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x*y^2 + y*z^2 + z*x^2 # - R[1] # = 0
x*z^2 + y*x^2 + z*y^2 # - R[2] # = 0
x*y + x*z + y*z + b[3]*(x + y + z) # - R[3] # = 0

### Classic Polynomial:
round0.p(poly.calc(x)) * 3

err = 31/9 - 5*x + 40*x^2 + 167*x^3 - 256*x^4 + 1346*x^5 - 1136*x^6 + 46*x^7 - 557*x^8 +  
	+ 333*x^9 - 18*x^10 + 81*x^11 - 33*x^12 + 9*x^13 + 3*x^15
round0(err)


########################
########################

###############
### Order 3 ###
###############

### n = 3
x*y^3 + y*z^3 + z*x^3 - R1 # = 0
x*z^3 + y*x^3 + z*y^3 - R2 # = 0
x*y*z - R3 # = 0

### Solution:
### Eq 1 + Eq 2 =>
x*y^3 + y*z^3 + z*x^3 + x*z^3 + y*x^3 + z*y^3 - R1 - R2 # = 0
E2*(S^2 - 2*E2) - R3*S - R1 - R2 # = 0
2*E2^2 - E2*S^2 + R3*S + R1 + R2 # = 0
# E2 = (S^2 +/- sqrt(Det))/4

### see Section "Simple System":
R3*S^5 - 5*E2*R3*S^3 + (7*R3^2 - R1*E2)*S^2 + (E2^2*R3 + R1*R3)*S +
	+ R1^2 + 2*R1*E2^2 + E2^4 # = 0
- 128*R1*R2 + 64*R1^2 + 64*R2^2 - 32*Dsq*R1*S^2 - 32*Dsq*R2*S^2 + 1728*R3^2*S^2 - 320*Dsq*R3*S^3 +
	- 64*R1*S^4 - 64*R2*S^4 - 96*R3*S^5 + 8*Dsq*S^6 + 8*S^8 # = 0
4*Dsq*R1*S^2 + 4*Dsq*R2*S^2 + 40*Dsq*R3*S^3 - Dsq*S^6 =
	S^8 - 12*R3*S^5 - 8*R1*S^4 - 8*R2*S^4 + 216*R3^2*S^2 - 16*R1*R2 + 8*R1^2 + 8*R2^2
R3*S^13 - 26*R3^2*S^10 - (13*R1*R3 + 13*R2*R3)*S^9 - R1*R2*S^8 + 119*R3^3*S^7 +
	+ (186*R1*R3^2 + 186*R2*R3^2)*S^6 + (90*R1*R2*R3 + 39*R1^2*R3 + 39*R2^2*R3)*S^5 +
	+ (8*R1*R2^2 + 8*R1^2*R2 + 729*R3^4)*S^4 +
	(- 108*R1*R2*R3^2 + 54*R1^2*R3^2 + 54*R2^2*R3^2)*S^2 +
	(- 4*R1*R2^3 + 6*R1^2*R2^2 - 4*R1^3*R2 + R1^4 + R2^4)
### 2*E2^2 = E2*S^2 - R3*S - R1 - R2
R3*S^5 - 5*E2*R3*S^3 + (7*R3^2 - R1*E2)*S^2 + (E2^2*R3 + R1*R3)*S +
	+ R1^2 + 2*R1*E2^2 + E2^4 # = 0
R3*S^5 - 5*E2*R3*S^3 + (7*R3^2 - R1*E2)*S^2 + (1/2*(E2*S^2 - R3*S - R1 - R2)*R3 + R1*R3)*S +
	+ R1^2 + R1*(E2*S^2 - R3*S - R1 - R2) + 1/4 * (E2*S^2 - R3*S - R1 - R2)^2 # = 0
E2*(S^6 - 40*R3*S^3 - 4*S^2*(R1 + R2)) +
	7*R3*S^5 - (R1 + R2)*S^4 + 54*R3^2*S^2 + 2*(R1 - R2)^2 # = 0


### Solution
solve.ht3Dual = function(R, b=0) {
	if(length(b) == 1 && b[1] == 0) {
		coeff = c(
			R[3], 0, 0, - 26*R[3]^2, - (13*R[1]*R[3] + 13*R[2]*R[3]), - R[1]*R[2], 119*R[3]^3,
			(186*R[1]*R[3]^2 + 186*R[2]*R[3]^2), (90*R[1]*R[2]*R[3] + 39*R[1]^2*R[3] + 39*R[2]^2*R[3]),
			(8*R[1]*R[2]^2 + 8*R[1]^2*R[2] + 729*R[3]^4), 0,
			(- 108*R[1]*R[2]*R[3]^2 + 54*R[1]^2*R[3]^2 + 54*R[2]^2*R[3]^2), 0,
			(- 4*R[1]*R[2]^3 + 6*R[1]^2*R[2]^2 - 4*R[1]^3*R[2] + R[1]^4 + R[2]^4))
	} else {
		# TODO
	}
	if(length(b) > 1) {
		# Ext 2:
		# TODO
	}
	S = roots(coeff)
	len = length(S)
	print(S)
	b2 = if(length(b) > 1) b[2] else 0; # TODO: Ext 2;
	# TODO: Ext 2;
	# E2 = (S^2 + sign*sqrt(S^4 - 8*R[3]*S - 8*R[1] - 8*R[2]))/4
	# robust
	div = - (S^6 - 40*R[3]*S^3 - 4*S^2*(R[1] + R[2]));
	E2 = (7*R[3]*S^5 - (R[1] + R[2])*S^4 + 54*R[3]^2*S^2 + 2*(R[1] - R[2])^2) / div
	#
	x = sapply(seq(length(S)), function(id) roots(c(1, -S[id], E2[id], -R[3])))
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	### TODO: Div by 0
	yz = R[3]/x
	yz.s = S - x
	### robust:
	yz.d = (R[2] - R[1]) / (x^3 + yz*yz.s - x*(yz.s^2 - yz))
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

### Examples:

R = c(1, 3, -1);
sol = solve.ht3Dual(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3Dual(x, y, z, n=3)

# order 39
round0.p(poly.calc(x))

# order 39
err = 1 + 5*x^3 - 6*x^4 + 10*x^6 - 9*x^7 + 9*x^8 + 11*x^9 - 41*x^10 + 4*x^11 - 52*x^12 + 168*x^13 +
	+ 100*x^14 + 20*x^15 - 65*x^16 - 41*x^17 + 28*x^18 + 250*x^19 + 85*x^20 - 27*x^21 - 25*x^22 +
	+ 3*x^23 + 33*x^24 + 65*x^25 - 54*x^26 - 30*x^27 + 3*x^28 - x^30 - 16*x^32 - 4*x^33 + 3*x^34 + x^39
round0(err)


#########
### Ex 2:
R = c(0, -1, -1);
sol = solve.ht3Dual(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3Dual(x, y, z, n=3)

# order 39
round0.p(poly.calc(x))


######################
######################

####################
### Orders 2 & 3 ###
####################

### System:
x^3*y + y^3*z + z^3*x - R1 # = 0
x^2*y + y^2*z + z^2*x - R2 # = 0
x*y*z - R3 # = 0

### Solution:

### E21a + E21b:
E21a + E21b - E2*S + 3*E3 # = 0

### E21a * E21b:
E21a * E21b - E3*S^3 + 6*E2*E3*S - E2^3 - 9*E3^2 # = 0

### E31a + E31b:
E31a + E31b - E2*(S^2 - 2*E2) + E3*S # = 0

### E31a * E31b:
E31a * E31b - E3*S^5 + 5*E2*E3*S^3 - 7*E3^2*S^2 - E2^2*E3*S - E2^4 # = 0

### Alternative:
E31a - E21a*S - 2*E3*S + E2^2 + E3*S # = 0

### Eq S:
E3^2*S^6 + (2*E21a^2*E3 - 19*E21a*E3^2 - 7*E3^3)*S^3 - 5*(2*E21a*E3*E31a - 3*E3^2*E31a)*S^2 +
	- (E21a*E31a^2 - 9*E3*E31a^2)*S +
	+ E21a^4 + 6*E21a^3*E3 + 27*E21a^2*E3^2 + 54*E21a*E3^3 + 81*E3^4 + E31a^3 # = 0

### Solver:

solve.S3HtMix.D32 = function(R, debug=TRUE, all=FALSE) {
	coeff = coeff.S3HtMix.D32(R);
	S = roots(coeff);
	if(debug) print(S);
	E3 = R[3]; E31a = R[1]; E21a = R[2];
	E2x0 = E3*S^3 + E21a^2 + 3*E21a*E3 + 9*E3^2;
	E2 = E2x0 / (5*E3*S + E31a);
	#
	len = length(S);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3)));
	x = as.vector(x);
	S = rep(S, each=3); E2 = rep(E2, each=3);
	# Robust:
	s = S - x; e2 = E2 - s*x;
	# (x^2 + e2)*y + s*x*z = E21a + e2*x;
	y = (s^2*x - E21a - e2*x) / (s*x - x^2 - e2);
	z = s - y;
	sol = cbind(x, y, z);
	return(sol);
}
coeff.S3HtMix.D32 = function(R) {
	E31a = R[1]; E21a = R[2]; E3 = R[3];
	coeff = c(E3^2, 0, 0, (2*E21a^2*E3 - 19*E21a*E3^2 - 7*E3^3),
		- 5*(2*E21a*E3*E31a - 3*E3^2*E31a), - (E21a*E31a^2 - 9*E3*E31a^2),
		E21a^4 + 6*E21a^3*E3 + 27*E21a^2*E3^2 + 54*E21a*E3^3 + 81*E3^4 + E31a^3);
	return(coeff);
}

### Examples:

###
R = c(-1,3,2)
sol = solve.S3HtMix.D32(R)


### Test:
x = sol[,1]; y = sol[,2]; z = sol[,3];
x^3*y + y^3*z + z^3*x # - R[1]
x^2*y + y^2*z + z^2*x # - R[2]
x*y*z # - R[3]


### Debug:
R = c(3, -5, 2);
x =  0.7431564372 + 1.5834791598i;
y =  0.3769770644 + 0.7079940406i;
z = -0.8543963617 - 1.1410556005i;
S = x+y+z; E2 = (x+y)*z + x*y; E3 = x*y*z;
E21a = x^2*y + y^2*z + z^2*x;
E21b = x*y^2 + y*z^2 + z*x^2;
E31a = x^3*y + y^3*z + z^3*x;
E31b = x*y^3 + y*z^3 + z*x^3;


### Eqs:
E21a^2 + 3*E3*E21a - E2*S*E21a + 9*E3^2 + E3*S^3 - 6*E2*E3*S + E2^3 # = 0

E31a^2 + E3*S*E31a - E2*S^2*E31a + 2*E2^2*E31a + E3*S^5 + 7*E3^2*S^2 - 5*E2*E3*S^3 + E2^2*E3*S + E2^4 # = 0


### Derivation:

prod.S3E2ab(2,1)

###
pE21 = toPoly.pm("E21a * E21b - E3*S^3 + 6*E2*E3*S - E2^3 - 9*E3^2")
pE21b = toPoly.pm("E21a + E21b - E2*S + 3*E3")
pE21 = solve.pm(pE21, pE21b, "E21b")$Rez

pE31 = toPoly.pm("E31a * E31b - E3*S^5 + 5*E2*E3*S^3 - 7*E3^2*S^2 - E2^2*E3*S - E2^4")
pE31b = toPoly.pm("E31a + E31b - E2*(S^2 - 2*E2) + E3*S")
pE31 = solve.pm(pE31, pE31b, "E31b")$Rez

### [old]
# pR = solve.pm(pE21, pE31, "E2")
# str(pR) # 172 Monomials


p0 = toPoly.pm("E31a - E21a*S - 2*E3*S + E2^2 + E3*S")

pR = solve.pm(p0, pE21, "E2")

