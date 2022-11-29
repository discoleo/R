########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S3:
### Mixed Type
###
### draft v.0.4a


### Heterogeneous Symmetric
### Polynomial Systems: 3 Variables
### Mixed: Hetero + Symmetric

### Example:
# x^p*y^n + y^p*z^n + z^p*x^n = R1
# x*y + x*z + y*z = R2
# x*y*z = R3


###############

###############
### History ###
###############


### draft v.0.4a:
# - [refactor] moved Dual-Systems to file:
#   Poly.System.Hetero.Symmetric.S3.Mixed.Dual.R;
### draft v.0.3c:
# - [started] Resonances for 3 Powers:
#  -- formula for 5 variables;
#  -- some quasi-non-trivial examples with 4 variables;
### draft v.0.3b - v.0.3b-spCase:
# - classic Polynomial: P[15] for S3P31;
# - [started] classic Polynomial: P[77] for Mixed S3P31 + Symmetric P7; [v.0.3b-clPP77]
# - solved Special Case: S == 0; [v.0.3b-spCase]
### draft v.0.3a:
# - solved: Mixed Ht S3P31 + Symmetric P7;
#   x^7 + y^7 + z^7 = R2;
### draft v.0.2p - v.0.2p-full:
# - Resonances with roots of unity;
### draft v.0.2o-clean:
# - more cleanup: moved to Derivation;
### draft v.0.2n - v.0.2n-clPoly:
# - solved: Mixed Ht S3P21 + Symmetric P3;
#   x^3 + y^3 + z^3 = R2;
# - Classic polynomial: degenerate P[27];
### draft v.0.2m - v.0.2m-ext:
# - extensions of type A to the HtDual E3-variant system;
# - extensions of type A to the HtDual E2-variant system:
#   P[12] => P[6]: degeneracy possible;
### draft v.0.2l:
# - moved Derivations to new file:
#   Poly.System.Hetero.Symmetric.S3.Mixt.Derivation.R;
# - extension type A1 for the Order 3 System;
### draft v.0.2k:
# - parametric P[9] for simple case;
### draft v.0.2j:
# - variant: E2*S = R2;
### draft v.0.2i-exp - v.0.2i-ex:
# - minor experiments with redundancy;
# - more examples;
### draft v.0.2h - v.0.2h-ext1:
# - Dual system with E2 = R3;
# - Extension: E2 + b1*(x+y+z) = R3; [v.0.2h-ext1]
### draft v.0.2g:
# - classic Polynomial for the simple Dual system:
#   degenerate P18: pseudo-P6;
### draft v.0.2f:
# - linear extension (type A3) to the Order 2 system;
### draft v.0.2e:
# - combined variant:
#   (x*y^2 + y*z^2 + z*x^2) - a*(x*z^2 + y*x^2 + z*y^2) = R1;
### draft v.0.2d:
# - Extensions to the Order 2 system:
#   M3 extension; [v.0.2d]
### draft v.0.2c - v.0.2c-robust:
# - added system with 2 hetero-symmetric equations:
#   Order 2: trivial polynomials;
#   Order 3: intermediary polynomial of order 13; [v.0.2c-ord3]
# - [DONE]: fixed robust roots! [v.0.2c-robust]
### draft v.0.2a - v.0.2b-sp: [08-12-2020]
# - Generalization:
#   x^p*y^n + y^p*z^n + z^p*x^n = R1;
# - solved: n = 3, p = 2; [full in v.0.2b]
# - [TODO][DONE]: robust special case; [v.0.2b-sp]
### draft v.0.1c: [08-12-2020]
# - full robust solution: using dS / dR1;
### draft v.0.1b: [07-12-2020]
# - solved + robust:
#   x*y^3 + y*z^3 + z*x^3 = R1;
### draft v.0.1a-ext1 - v.0.1a-ext2:
# - Extensions:
#   x*y^2 + y*z^2 + z*x^2 + b1*(x+y+z) = R1;
#   x*y + x*z + y*z + b2*(x+y+z) = R2;
### draft v.0.1a:
# - moved to new file
#   from Poly.System.Hetero.Symmetric.S3.R;
#
#########################
### former branch v.0.1c:
### in Poly.System.Hetero.Symmetric.S3.R
### draft v.0.1c-pre-alpha - v.0.1c-exact:
# - first look & solved + exact/robust solution: [v.0.1c-exact]
#   x*y^2 + y*z^2 + z*x^2 = R1;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


# library(polynom)
# library(pracma)

# the functions are in the file:
# Polynomials.Helper.R

### other functions

resonance = function(p, n=3) {
	# "Resonance" with Roots of unity;
	# currently only for 2-variable terms: x^p1*y^p2;
	sg = if(n %% 2 == 1) 1 else -1;
	p.all = p[1]^n + sg*p[2]^n;
	r = list(p=p.all, f=factors(p.all), p.trivial = sum(p));
	return(r);
}
test.S3HtM = function(sol, b.ext=0, R=NULL, n) {
	test.ht3(sol, b=b.ext, R=R, n=n[1], p=n[2]);
}
test.ht3 = function(x, y, z, R=NULL, n=2, p=1, b=0) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	### Test
	x.sum = x + y + z
	err1 = x^p*y^n + y^p*z^n + z^p*x^n + b[1]*x.sum
	err2 = x*y + y*z + z*x + if(length(b) < 2) 0 else b[2]*x.sum
	err3 = x*y*z + if(length(b) < 3) 0 else b[3]*x.sum
	err = rbind(err1, err2, err3)
	if( ! is.null(R)) {
		err = err - rep(R, each=length(x))
	}
	err = round0(err)
	return(err)
}
test.Ht3Dual = function(x, y, z, R, n=2, p=1, b=0, type) {
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

############################

##########################
### Mixed Systems      ###
### Type 1: 1 Rotation ###
##########################

### Generalization:
### x^p*y^n + y^p*z^n + z^p*x^n
### x*y + x*z + y*z
### x*y*z

################################

##############
### Simple ###
### p = 1  ###
##############

### x*y^n + y*z^n + z*x^n = R1

###############
### Order 2 ###
###############

### n = 2
x*y^2 + y*z^2 + z*x^2 - R1 # = 0
x*y + y*z + z*x - R2 # = 0
x*y*z - R3 # = 0

### Solution:

### Eq:
E3*S^3 - (R1+6*E3)*E2*S + R1^2 + E2^3 + 9*E3^2 + 3*R1*E3 # = 0


### Solver:
solve.Ht3 = function(R, b=0, debug=TRUE) {
	if(all(b == 0)) {
		coeff = c(R[3], 0, - (R[1]+6*R[3])*R[2], R[1]^2 + R[2]^3 + 9*R[3]^2 + 3*R[1]*R[3])
	} else if(length(b) < 3) {
		coeff = c(R[3], (b[1]^2 + b[1]*R[2]), - (R[1]*R[2] + 3*b[1]*R[3] + 6*R[2]*R[3] + 2*b[1]*R[1]),
			R[1]^2 + R[2]^3 + 9*R[3]^2 + 3*R[1]*R[3])
		if(length(b) > 1) {
			# Ext 2:
			coeff = coeff - c(b[2]^3 + b[1]*b[2], -b[2]*(R[1] + 3*b[2]*R[2] + 6*R[3]), 3*b[2]*R[2]^2, 0)
		}
	} else {
		# Ext 3:
		coeff = c(-b[3], R[3], (b[1]^2 + 3*b[1]*b[3] + 9*b[3]^2 + b[1]*R[2] + 6*b[3]*R[2]),
			- (R[1]*R[2] + 3*b[1]*R[3] + 6*R[2]*R[3] + 2*b[1]*R[1] + 3*b[3]*R[1] + 18*b[3]*R[3]),
			R[1]^2 + R[2]^3 + 9*R[3]^2 + 3*R[1]*R[3])
		coeff = coeff - c(0, b[2]^3 + b[1]*b[2] + 6*b[2]*b[3],
			- b[2]*(R[1] + 3*b[2]*R[2] + 6*R[3]), 3*b[2]*R[2]^2, 0) # Ext 2
	}
	S = roots(coeff);
	if(debug) print(S);
	len = length(S)
	b2 = if(length(b) > 1) b[2] else 0; # Ext 2;
	b3 = if(length(b) > 2) b[3] else 0; # Ext 3;
	x = sapply(S, function(S) roots(c(1, -S, R[2] - b2*S, - R[3] + b3*S)))
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	yz = (R[3] - b3*S) / x
	yz.s = S - x
	### robust:
	# x*y^2 - (x^2+yz)*y + (x^2+yz)*yz.s + b1*S - R1
	# x*y^2 + x*yz - x*y*yz.s = 0 # x*y*(y+z - yz.sum) = 0
	# (x^2+yz - x*yz.s)*y + x*yz - (x^2+yz)*yz.s - b1*S + R1 = 0
	y = - (x*yz - (x^2+yz)*yz.s - b[1]*S + R[1]) / (x^2+yz - x*yz.s)
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

### Examples

### Ex 1:
R = c(1, 1, 1)
b = 0
sol = solve.Ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

round0.p(poly.calc(x))

err = -1 + 3*x - 3*x^2 + 4*x^3 + x^4 - 4*x^5 + 11*x^6 - 4*x^7 + x^9
round0(err)

### P[9]
R1 = R[1]; R2 = R[2]; R3 = R[3];
R3*x^9 - R2*(R1 + 3*R3)*x^7 + (R1^2 + R2^3 + 3*R1*R3 + 6*R3^2)*x^6 +
	- R2^2*(R1 + 3*R3)*x^5 + R1*R2*R3*x^4 + R3*(R2^3 + 3*R3^2)*x^3 +
	- 3*R2^2*R3^2*x^2 + 3*R2*R3^3*x - R3^4


#########
### Ex 2:
R = c(-3*6, -6, 6)
b = 0
sol = solve.Ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

round0.p(poly.calc(x))
err = -216 - 648*x - 648*x^2 - 108*x^3 + 108*x^4 + x^9
round0(err)
# also R = c(-3*3, -3, 3)
# -27 - 81*x - 81*x^2 + 27*x^4 + 9*x^6 + x^9


#########
### Ex 3:
k = 1 # trivial
R = c(3*k^3, 3*k^2, k^3)
b = 1 # not trivial anymore
sol = solve.Ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

round0.p(poly.calc(x))


#########
### Ex 3:
R = c(0, 1, 1)
sol = solve.Ht3(R)

### Test
test.ht3(sol)

round0.p(poly.calc(sol[,1]))
x = sol[,1]
err = -1 + 3*x - 3*x^2 + 4*x^3 - 3*x^5 + 7*x^6 - 3*x^7 + x^9
round0(err)


###############
### Extensions:

### Extension A1:
### x*y^2 + y*z^2 + z*x^2 + b1*(x+y+z) = R1
### Extension A2:
### x*y + x*z + y*z + b2*(x+y+z) = R2;
### Extension A3:
### x*y*z + b3*(x+y+z) = R3;

### Sum - R1*initial Eq
R3*S^3 + (b1^2 + b1*R2)*S^2 - (R1*R2 + 3*b1*R3 + 6*R2*R3 + 2*b1*R1)*S +
	+ R1^2 + R2^3 + 9*R3^2 + 3*R1*R3 # = 0

### Extension A3:
# - includes A1, but A2 was missed;
#   [see function solve.Ht3() for complete variant]
# E3 = R3 - b3*S
-b3*S^4 + R3*S^3 + (b1^2 + 3*b1*b3 + 9*b3^2 + b1*R2 + 6*b3*R2)*S^2 +
	- (R1*R2 + 3*b1*R3 + 6*R2*R3 + 2*b1*R1 + 3*b3*R1 + 18*b3*R3)*S +
	+ R1^2 + R2^3 + 9*R3^2 + 3*R1*R3 # = 0


### Examples

###
R = c(1, 1, 1); b = 1;
sol = solve.Ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

round0.p(poly.calc(x))
err = -1 + 3*x - x^2 + 8*x^4 - 13*x^5 + 15*x^6 - 9*x^7 + 2*x^8 + x^9
round0(err)


### Ex 2:
R = c(0, 0, 1);
b = 1; # b = -2;
sol = solve.Ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

round0.p(poly.calc(x[4:9]))
# works only with b = 1, b = -2
err = 1 + 2*x^2 - 2*x^3 + 3*x^4 - 2*x^5 + x^6
# err = 1 - x^2 - 2*x^3 + 3*x^4 + x^5 + x^6
round0(err)


### Ex 3:
R = c(0, 1, -1); b = 1;
sol = solve.Ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

round0.p(poly.calc(x))

err = 1 + 3*x + x^2 - 5*x^4 - 10*x^5 - 11*x^6 - 6*x^7 - 2*x^8 + x^9
round0(err)


###
# R3 = (b[2]^2 + b[1]*b[2]) +/- 1;
R = c(0, 1, 9); b = c(1, 2);
sol = solve.Ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

poly.calc(x)

# numeric instability with 1st root!
err = -729 - 19521*x - 8649*x^2 - 6428*x^3 - 4979*x^4 - 2234*x^5 - 653*x^6 - 154*x^7 - 122*x^8 + x^9
round0(err)


###
# R3 = (b[2]^2 + b[1]*b[2]) +/- 1;
R = c(0, 1, -2); b = c(2, -1);
sol = solve.Ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

poly.calc(x)

err = 8 - 72*x + 60*x^2 - 42*x^3 + 60*x^4 - 12*x^5 + 31*x^6 + 9*x^7 + 21*x^8 + x^9
round0(err)


### Ext 3:
R = c(1, 1, 0); b = c(-7, 1, 1);
sol = solve.Ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

round0.p(poly.calc(x))

err = -2 + 15*x - 65*x^2 + 71*x^3 - 15*x^4 + 96*x^5 + 14*x^6 - 95*x^7 - 57*x^8 +
	- 87*x^9 - 36*x^10 + x^12
round0(err)


### Ext 3, ex 2:
R = c(0, 1, 0); b = c(-7, 1, 1);
sol = solve.Ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

round0.p(poly.calc(x))

err = -1 + x - 32*x^2 + 67*x^3 - 39*x^4 + 86*x^5 - 4*x^6 - 79*x^7 - 25*x^8 +
	- 75*x^9 - 35*x^10 + x^12
round0(err)


#########################
#########################

### x*y^n + y*z^n + z*x^n = R1

###############
### Order 3 ###
###############

x*y^3 + y*z^3 + z*x^3 - R1 # = 0
x*y + y*z + z*x - R2 # = 0
x*y*z - R3 # = 0

### Solution:

### Eq:
E3*S^5 - 5*E2*E3*S^3 + (7*E3^2 - R1*E2)*S^2 + (E2^2*E3 + R1*E3)*S +
	+ R1^2 + 2*R1*E2^2 + E2^4 # = 0
R3*S^5 - 5*R2*R3*S^3 + (7*R3^2 - R1*R2)*S^2 + (R2^2*R3 + R1*R3)*S +
	+ R1^2 + 2*R1*R2^2 + R2^4 # = 0


### Solver:
solve.Ht3.S3P3 = function(R, b=0, debug=TRUE) {
	coeff = c(R[3], 0, - 5*R[2]*R[3], (7*R[3]^2 - R[1]*R[2]),
			(R[2]^2*R[3] + R[1]*R[3]), R[1]^2 + 2*R[1]*R[2]^2 + R[2]^4)
	if(any(b != 0)) {
		if(b[1] != 0) {
			coeff = coeff + c(0, 0, b[1]*R[2], -b[1]*R[3] + b[1]^2,
				-2*b[1]*R[1] - 2*b[1]*R[2]^2, 0);
		}
	}
	S = roots(coeff);
	if(debug) print(S);
	len = length(S);
	b2 = if(length(b) > 1) b[2] else 0;
	b3 = if(length(b) > 2) b[3] else 0;
	R1 = R[1] - b[1]*S;
	R2 = R[2] - b2*S;
	R3 = R[3] - b3*S;
	x = sapply(seq(len), function(id) roots(c(1, -S[id], R2[id], -R3[id])))
	S = matrix(S, ncol=len, nrow=3, byrow=T);
	R1 = rep(R1, each=3); R2 = rep(R2, each=3); R3 = rep(R3, each=3);
	yz = R3/x;
	yz.s = S - x;
	### robust:
	x3 = ifelse(R1 == 0, {
		# with chain rule!
		# x3 = (5*R[3]*S^4 - 15*R[2]*R[3]*S^2 + 14*R[3]^2*S + R[2]^2*R[3]) # * dS/dR
		dS = - R2*S^2 + R3*S + 2*R2^2;
		x3 = - dS;
		x3
	}, { # else {
		x3 = (R3*S^5 - 5*R2*R3*S^3 + 7*R3^2*S^2 + R2^2*R3*S + R2^4) / R1;
		x3
	})
	yz.d = (x3 - R1) / (x^3 + yz*yz.s - x*(yz.s^2 - yz));
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
}

### Examples:

R = c(1, 1, -2);
sol = solve.Ht3.S3P3(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x*y^3 + y*z^3 + z*x^3 # - R[1] # = 0
x*y + y*z + z*x # - R[2] # = 0
x*y*z # - R[3] # = 0


round0.p(poly.calc(x) * R[3])

err = 32 + 80*x + 80*x^2 + 120*x^3 + 130*x^4 + 61*x^5 + 36*x^6 + 6*x^7 - 9.5*x^8 - 17*x^9 +
	- 19*x^10 - 3*x^11 - 3.5*x^12 + x^15
round0(err)


### Classic Polynomial:
R1 = R[1]; R2 = R[2]; R3 = R[3];
R3*x^15 - (R1*R2 - 2*R3^2)*x^12 + R3*(R1 - 4*R2^2)*x^11 + (R1^2 + R2^4 + 9*R2*R3^2)*x^10 +
	+ R3*(3*R1*R2 - 4*R2^3 - 4*R3^2)*x^9 - (R1*R2^3 + R1*R3^2 - 6*R2^2*R3^2)*x^8 +
	+ R2*R3*(2*R1*R2 + R3^2)*x^7 - R3^2*(R1*R2 + 5*R2^3 + 3*R3^2)*x^6 +
	+ R2^2*R3*(R2^3 + 15*R3^2)*x^5 - 5*R2*R3^2*(R2^3 + 3*R3^2)*x^4 +
	+ 5*R3^3*(2*R2^3 + R3^2)*x^3 - 10*R2^2*R3^4*x^2 + 5*R2*R3^5*x - R3^6


#########
### Ex 2:
R = c(0, 2, -2);
sol = solve.Ht3.S3P3(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x*y^3 + y*z^3 + z*x^3 # - R[1] # = 0
x*y + y*z + z*x # - R[2] # = 0
x*y*z # - R[3] # = 0

round0.p(poly.calc(x))

err = 32 + 160*x + 320*x^2 + 400*x^3 + 400*x^4 + 272*x^5 + 104*x^6 + 8*x^7 - 48*x^8 - 48*x^9 +
	- 44*x^10 - 16*x^11 - 4*x^12 + x^15
round0(err)


#########
### Ex 3:
R = c(-1,0,1);
b = 1
sol = solve.Ht3.S3P3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x*y^3 + y*z^3 + z*x^3 + b[1]*(x+y+z) # - R[1] # = 0
x*y + y*z + z*x # - R[2] # = 0
x*y*z # - R[3] # = 0

round0.p(poly.calc(x))

err = -1 + 5*x^3 - 3*x^6 - x^8 - 4*x^9 + x^10 + x^11 + 2*x^12 + x^15
round0(err)


#############################
#############################

### Generalization:
### x^p*y^n + y^p*z^n + z^p*x^n

#####################
### Higher Powers ###
### p > 1         ###
#####################

### p = 2
### x^2*y^n + y^2*z^n + z^2*x^n = R1

### Order 3: n = 3, p = 2
x^2*y^3 + y^2*z^3 + z^2*x^3 - R1 # = 0
x*y + y*z + z*x - R2 # = 0
x*y*z - R3 # = 0

### Solution:

# Sum() - R1*(initial Eq) + R1^2 =>
R3^2*S^4 + (2*R1*R3 + R2*R3^2)*S^2 - (R1*R2^2 + 5*R2^3*R3)*S +
	+ R1^2 + R2^5 + 7*R2^2*R3^2 + R1*R2*R3 # = 0

### Solver:
solve.Ht3.S3P32 = function(R, b=0) {
	if(length(b) == 1 && b[1] == 0) {
		coeff = c(R[3]^2, 0, (2*R[1]*R[3] + R[2]*R[3]^2), - (R[1]*R[2]^2 + 5*R[2]^3*R[3]),
			R[1]^2 + R[2]^5 + 7*R[2]^2*R[3]^2 + R[1]*R[2]*R[3])
	} else {
		# TODO
	}
	if(length(b) > 1) {
		# Ext 2:
		# TODO
	}
	S = roots(coeff)
	print(S)
	b2 = if(length(b) > 1) b[2] else 0; # TODO: Ext 2;
	x = sapply(S, function(x) roots(c(1,-x, R[2] - b2*x, -R[3])))
	S = matrix(S, ncol=4, nrow=3, byrow=T)
	yz = R[3]/x
	yz.s = S - x
	### robust:
	if(R[1] == 0) {
		# with chain rule!
		dS = 2*R[3]*S^2 - R[2]^2*S + R[2]*R[3]
		x3 = - dS
	} else {
		x3 = (R[3]^2*S^4 + R[2]*R[3]^2*S^2 - 5*R[2]^3*R[3]*S + R[2]^5 + 7*R[2]^2*R[3]^2) / R[1]
	}
	yz.d = (x3 - R[1]) / (x^3*yz.s + yz^2 - x^2*(yz.s^2 - yz))
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

### Examples:

R = c(1, 1, -1);
sol = solve.Ht3.S3P32(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2*y^3 + y^2*z^3 + z^2*x^3 # - R[1] # = 0
x*y + y*z + z*x # - R[2] # = 0
x*y*z # - R[3] # = 0


poly.calc(x)

err = 1 + 4*x + 6*x^2 + 8*x^3 + 12*x^4 + 10*x^5 + 13*x^6 + 14*x^7 + 12*x^8 + 8*x^9 + 3*x^10 + x^12
round0(err)


#########
### Ex 2:
R = c(0, 2, -1);
sol = solve.Ht3.S3P32(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2*y^3 + y^2*z^3 + z^2*x^3 # - R[1] # = 0
x*y + y*z + z*x # - R[2] # = 0
x*y*z # - R[3] # = 0


poly.calc(x)

err = 1 + 8*x + 24*x^2 + 36*x^3 + 42*x^4 + 56*x^5 + 86*x^6 + 108*x^7 + 92*x^8 + 44*x^9 + 10*x^10 + x^12
round0(err)


###############
### Extensions:

### TODO: A1, A2, A3;



#############################
#############################

################
### Variants ###
################

### Variant: E2*S = R2

### Order 2: n = 2
x*y^2 + y*z^2 + z*x^2 - R1 # = 0
(x*y + y*z + z*x)*S - R2 # = 0
x*y*z - R3 # = 0

### Solution:

### Eq:
E3*S^3 - (R1+6*E3)*E2*S + R1^2 + E2^3 + 9*E3^2 + 3*R1*E3 # = 0
# =>
R3*S^6 + (R1^2 - R1*R2)*S^3 + 3*(3*R3^2 + R1*R3 - 2*R2*R3)*S^3 + R2^3 # = 0


### Solver:
solve.ht3v.S3P21 = function(R, b=0, debug=TRUE) {
	if(all(b == 0)) {
		coeff = c(R[3], 0, 0, (R[1]^2 - R[1]*R[2]) + 3*(3*R[3]^2 + R[1]*R[3] - 2*R[2]*R[3]),
			0, 0, R[2]^3);
	} else {
		# TODO
	}
	S = roots(coeff);
	if(debug) print(S);
	S = S[S != 0];
	len = length(S)
	b2 = if(length(b) > 1) b[2] else 0; # Ext 2;
	b3 = if(length(b) > 2) b[3] else 0; # Ext 3;
	R2 = R[2]/S - b2;
	R3 = R[3] - b3*S;
	x = sapply(seq(len), function(id) roots(c(1, -S[id], R2[id], - R3[id])))
	# S = matrix(S, ncol=len, nrow=3, byrow=T);
	S = rep(S, each=3);
	R3 = rep(R3, each=3);
	yz = R3 / x;
	yz.s = S - x;
	### robust:
	# x*y^2 - (x^2+yz)*y + (x^2+yz)*yz.s + b1*S - R1
	# x*y^2 + x*yz - x*y*yz.s = 0 # x*y*(y+z - yz.sum) = 0
	# (x^2+yz - x*yz.s)*y + x*yz - (x^2+yz)*yz.s - b1*S + R1 = 0
	y = - (x*yz - (x^2+yz)*yz.s - b[1]*S + R[1]) / (x^2+yz - x*yz.s)
	z = yz.s - y;
	cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
}

### Examples:

R = c(0, -1, 1)
sol = solve.ht3v.S3P21(R);
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test = test.ht3(sol);
test[2,] = test[2,] * (x+y+z);
test

###
round0.p(poly.calc(x))
# reducible to symmetric P[6];


#############################
#############################

################
### Variants ###
################

### Mixed with Symmetric
### x^3 + y^3 + z^3 = R2

### Ht Order 2+1:
# x*y^2 + y*z^2 + z*x^2 = R1
# x^3 + y^3 + z^3 = R2
# x*y*z = R3

### Solution:

### Eq 1:
E3*S^3 - (R1+6*E3)*E2*S + R1^2 + E2^3 + 9*E3^2 + 3*R1*E3 # = 0

### Eq 2:
S^3 - 3*E2*S + 3*E3 - R2 # = 0

### Auxiliary Eqs:
# 3*E2*S = S^3 + 3*E3 - R2
# Case: S = 0
# E2^3 + (R1^2 + 9*E3^2 + 3*R1*E3) = 0; # all 3 solutions

### Eq S:
S^9 - (9*R1 + 3*R2 + 18*E3)*S^6 + (3*R2^2 + 36*R2*E3 + 108*E3^2 + 9*R1*R2 + 54*R1*E3 + 27*R1^2)*S^3 +
	- R2^3 + 27*E3^3 + 9*R2^2*E3 - 27*R2*E3^2

### Solver:
solve.S3Mixed.P21P3 = function(R, b=0, debug=TRUE) {
	R1 = R[1]; R2 = R[2]; E3 = R[3];
	coeff = c(1, 0, 0, - (9*R1 + 3*R2 + 18*E3), 0, 0,
		(3*R2^2 + 36*R2*E3 + 108*E3^2 + 9*R1*R2 + 54*R1*E3 + 27*R1^2), 0, 0,
		- R2^3 + 27*E3^3 + 9*R2^2*E3 - 27*R2*E3^2);
	S = roots(coeff);
	if(debug) print(S);
	len = length(S);
	E2 = (S^3 + 3*E3 - R2) / (3*S);
	if(any(S == 0)) E2[S == 0] = roots(c(1,0,0, 9*E3^2 + R1^2 + R1*R2));
	E3 = rep(E3, len);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	S = rep(S, each=3); E3 = rep(E3, each=3);
	yz.s = S - x; yz = E3 / x;
	# robust
	y = - (x*yz - (x^2+yz)*yz.s - b[1]*S + R[1]) / (x^2 + yz - x*yz.s)
	z = yz.s - y
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z));
	return(sol);
}

### Examples:
R = c(-1,-3,1)
sol = solve.S3Mixed.P21P3(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test:
x*y^2 + y*z^2 + z*x^2 # - R1
x^3 + y^3 + z^3 # - R2
x*y*z # - R3

### Classic Polynomial:
# - degenerate P27;
round0.p(poly.calc(x))

R1 = R[1]; R2 = R[2]; R3 = R[3];
x^27 - 3*R2*x^24 + 3*(R2^2 - R2*R3 - 3*R3^2)*x^21 +
	+ (6*R2^2*R3 + 18*R2*R3^2 - R2^3 - 3*R3^3)*x^18 +
	- (R1^3*R2 + 6*R1^3*R3 + 3*R2^3*R3 + 6*R2^2*R3^2 - 24*R2*R3^3 - 27*R3^4)*x^15 +
	+ (R1^3*R2^2 + 6*R1^3*R2*R3 - 3*R2^3*R3^2 - 21*R2^2*R3^3 - 21*R2*R3^4 + 18*R3^5)*x^12 +
	- (R1^6 - 6*R1^3*R2*R3^2 - 9*R1^3*R3^3 + R2^3*R3^3 + 15*R2^2*R3^4 + 45*R2*R3^5 + 24*R3^6)*x^9 +
	+ (R1^3*R2*R3^3 + 6*R1^3*R3^4 - 3*R2^2*R3^5 - 21*R2*R3^6 - 27*R3^7)*x^6 +
	- (3*R2*R3^7 + 9*R3^8)*x^3 - R3^9


###########################
###########################

### Mixed with Symmetric P7
### x^7 + y^7 + z^7 = R2

### Ht Order 3+1:
# x*y^3 + y*z^3 + z*x^3 = R1
# x^7 + y^7 + z^7 = R2
# x*y*z = R3


### Solution:

### Note:
# - Eq S is P[28];
# - but classic polynomial is a degenerate P[84];
# TODO:
# - evaluate if the P[12] obtained from the P[84] could be used to back-solve the P[28];

### Eq 1:
R3*S^5 - 5*E2*R3*S^3 + (7*R3^2 - R1*E2)*S^2 + (E2^2*R3 + R1*R3)*S +
	+ R1^2 + 2*R1*E2^2 + E2^4 # = 0

### Eq 2:
S^7 - 7*E2*S^5 + 7*E3*S^4 + 14*E2^2*S^3 - 21*E3*E2*S^2 - 7*E2^3*S + 7*E3^2*S + 7*E3*E2^2 - R2 # = 0

### Eq S: P[28]
# see Derivation;


### Solver:
solve.S3P31SymmP7 = function(R, debug=TRUE) {
	coeff = coeff.S3P31SymmP7(R);
	S = roots(coeff);
	S = round0(S);
	if(debug) print(S);
	len = length(S);
	E2 = E2.S3P31SymmP7(S, R);
	E3 = R[3]; E3 = rep(E3, len);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	S = rep(S, each=3); E2 = rep(E2, each=3); E3 = rep(E3, each=3);
	yz.s = S - x; yz = E3 / x;
	# robust
	R1 = R[1]; R2 = E2; R3 = E3;
	x3 = if(R1 == 0) {
		# with chain rule!
		# x3 = (5*R[3]*S^4 - 15*R[2]*R[3]*S^2 + 14*R[3]^2*S + R[2]^2*R[3]) # * dS/dR
		dS = - R2*S^2 + R3*S + 2*R2^2;
		x3 = - dS;
		x3
	} else {
		x3 = (R3*S^5 - 5*R2*R3*S^3 + 7*R3^2*S^2 + R2^2*R3*S + R2^4) / R1;
		x3
	};
	xdiv = (x^3 + yz*yz.s - x*(yz.s^2 - yz));
	yz.d = (x3 - R1) / xdiv;
	if(any(S == 0)) {
		isZero = (S == 0);
		half = sum(isZero) / 2;
		yz.d.half = sqrt(yz.s[isZero][seq(half)]^2 - 4*yz[isZero][seq(half)] + 0i);
		yz.d[isZero][seq(half)] = yz.d.half;
		# - sqrt();
		yz.d[isZero][seq(half+1, 2*half)] = - yz.d.half;
	}
	y = (yz.s + yz.d) / 2;
	z = yz.s - y;
	cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
}
coeff.S3P31SymmP7 = function(R) {
	R1 = R[1]; R2 = R[2]; E3 = R[3];
	coeff = c(1, 0, 0, - 56*E3, - 28*R1, 0, 1190*E3^2, - 4*R2 + 1148*E3*R1, 294*R1^2,
		- 11564*E3^3, 560*R2*E3 - 15386*E3^2*R1, - 210*R2*R1 - 9114*E3*R1^2,
		47383*E3^4 - 1274*R1^3, - 13370*R2*E3^2 + 62818*E3^3*R1,
		6*R2^2 + 476*R2*E3*R1 + 81928*E3^2*R1^2, - 58996*E3^5 + 784*R2*R1^2 + 19208*E3*R1^3,
		91728*R2*E3^3 + 93982*E3^4*R1 + 6517*R1^4,
		1792*R2^2*E3 + 41062*R2*E3^2*R1 - 173558*E3^3*R1^2,
		148862*E3^6 + 161*R2^2*R1 + 10682*R2*E3*R1^2 - 5831*E3^2*R1^3,
		- 48118*R2*E3^4 - 461678*E3^5*R1 + 5978*R2*R1^3 - 25382*E3*R1^4,
		12194*R2^2*E3^2 - 72128*R2*E3^3*R1 + 124509*E3^4*R1^2 + 4802*R1^5,
		- 4*R2^3 - 67228*E3^7 + 3150*R2^2*E3*R1 + 9016*R2*E3^2*R1^2 - 92610*E3^3*R1^3,
		- 104272*R2*E3^5 + 177674*E3^6*R1 + 1323*R2^2*R1^2 - 4802*R2*E3*R1^3 + 50421*E3^2*R1^4,
		6272*R2^2*E3^3 + 26068*R2*E3^4*R1 - 57624*E3^5*R1^2 + 686*R2*R1^4 - 4802*E3*R1^5,
		105*R2^3*E3 + 117649*E3^8 - 1666*R2^2*E3^2*R1 - 1715*R2*E3^3*R1^2 + 72030*E3^4*R1^3 + 2401*R1^6,
		- 4802*R2*E3^6 + 77*R2^3*R1 - 33614*E3^7*R1 + 833*R2^2*E3*R1^2 + 1029*R2*E3^2*R1^3 - 7203*E3^3*R1^4,
		# 49*(15*R2^2*E3^4 + 210*R2*E3^5*R1 + 735*E3^6*R1^2 + 2*R2^2*R1^3 + 28*R2*E3*R1^4 + 98*E3^2*R1^5),
		49*(R2 + 7*R1*E3)^2 * (2*R1^3 + 15*E3^4),
		- 14*(R2 + 7*R1*E3)^3 * E3^2,
		(R2 + 7*R1*E3)^4);
	return(coeff)
}
E2.S3P31SymmP7 = function(S, R, digits=4) {
	R1 = R[1]; R2 = R[2]; E3 = R[3];
	px0 = 5*S^17 - 57*E3*S^14 + 42*R1*S^13 - 343*E3^2*S^11 - 3*R2*S^10 - 49*E3*R1*S^10 - 175*R1^2*S^9 +
		+ 2541*E3^3*S^8 - 89*R2*E3*S^7 - 574*E3^2*R1*S^7 - 42*R2*R1*S^6 + 245*E3*R1^2*S^6 - 245*E3^4*S^5 +
		- 147*R1^3*S^5 - 637*E3^3*R1*S^4 - 2*R2^2*S^3 - 343*E3^5*S^2 - 21*R2*R1^2*S^2 - 98*E3*R1^3*S^2 +
		+ 7*R2*E3^3*S + 49*E3^4*R1*S - R2^2*E3 - 14*R2*E3^2*R1 - 49*E3^3*R1^2;
	pDiv = 22*S^15 - 448*E3*S^12 + 84*R1*S^11 + 2030*E3^2*S^9 + 26*R2*S^8 + 98*E3*R1*S^8 - 98*R1^2*S^7 +
		- 392*E3^3*S^6 + 154*R2*E3*S^5 - 1176*E3^2*R1*S^5 + 14*R2*R1*S^4 + 294*E3*R1^2*S^4 - 686*E3^4*S^3 +
		- 98*R1^3*S^3 + 28*R2*E3^2*S^2 + 196*E3^3*R1*S^2 + R2^2*S - 49*E3^2*R1^2*S;
	isZero = (round(pDiv, digits) == 0);
	if(any(isZero)) {
		isSZero = (S == 0);
		isZero  = isZero & (! isSZero);
		# TODO: check & improve! R = c(1,-1,-1)
		px0[isZero]  = roots(c(1,-2,2));
		pDiv[isZero] = 1;
		# S == 0
		px0[isSZero]  = roots(c(7*E3, 0, -R2));
		pDiv[isSZero] = 1;
		print("Gargamel!")
	}
	return(px0 / pDiv);
}

### Examples:

R = c(-1,2,-2)
sol = solve.S3P31SymmP7(R);
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x*y^3 + y*z^3 + z*x^3 # - R1
x^7 + y^7 + z^7 # - R2
x*y*z # - R3 # = 0

round0.p(poly.calc(x), tol=0.5)


### Ex 2:
# S = 0
R = c(1,-7,1)
sol = solve.S3P31SymmP7(R);
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x*y^3 + y*z^3 + z*x^3 # - R1
x^7 + y^7 + z^7 # - R2
x*y*z # - R3 # = 0

round0.p(poly.calc(x), tol=0.5)


### Classic Polynomial:
R1 = R[1]; R2 = R[2]; R3 = R[3];
x^84 - 4*R2*x^77 + (6*R2^2 - 7*R1*R2*R3 - 21*R1^2*R3^2) * x^70 +
	+ (63*R1^2*R2*R3^2 + 21*R1*R2^2*R3 - 4*R2^3 - 4*R3^7) * x^63 +
	+ (2*R1^7 + R2^4 - 21*R1*R2^3*R3 - 49*R1^2*R2^2*R3^2 + 84*R1^3*R2*R3^3 +
		+ 168*R1^4*R3^4 - R2*R3^7 + 7*R1*R3^8) * x^56 +
	+ (- 4*R1^7*R2 + 7*R1*R2^4*R3 - 7*R1^2*R2^3*R3^2 - 168*R1^3*R2^2*R3^3 - 336*R1^4*R2*R3^4 +
		+ 14*R2^2*R3^7 + 7*R1*R2*R3^8 + 63*R1^2*R3^9) * x^49 +
	+ (3*R1^7*R2^2 + 7*R1^8*R2*R3 + 14*(2*R1^9 + R1^2*R2^4)*R3^2 + 77*R1^3*R2^3*R3^3 + 105*R1^4*R2^2*R3^4 +
		- 294*R1^5*R2*R3^5 - 392*R1^6*R3^6 - 11*R2^3*R3^7 - 14*R1*R2^2*R3^8 - 63*R1^2*R2*R3^9 +
		+ 98*R1^3*R3^10 + 6*R3^14) * x^42 +
	+ 0; # TODO


###########################
###########################
###########################

########################
### Combined Variant ###
########################

############
### Order 2: n = 2
(x*y^2 + y*z^2 + z*x^2) - b*(x*z^2 + y*x^2 + z*y^2) - R1 # = 0
x*y + x*z + y*z - R2 # = 0
x*y*z - R3 # = 0

### Solution:
(b+1)*(x*y^2 + y*z^2 + z*x^2) - b*(E2*S - 3*E3) - R1 # = 0
(b+1)*(x*z^2 + y*x^2 + z*y^2) - (E2*S - 3*E3) + R1 # = 0
# =>
R3*(x^3+y^3+z^3) + (x^3*y^3+x^3*z^3+y^3*z^3) +
	- 1/(b+1) * (b*(E2*S - 3*E3) + R1)*(x^2*y + y^2*z + x*z^2) + 3*R3^2 # = 0
R3*(S^3 - 3*E2*S + 3*E3) + (R2^3 - 3*R3*(R2*S - R3)) +
	- 1/(b+1) * (b*(E2*S - 3*E3) + R1)*(x^2*y + y^2*z + x*z^2) + 3*R3^2 # = 0
R3*S^3 - 6*R2*R3*S + R2^3 + 9*R3^2 +
	- 1/(b+1) * (b*(E2*S - 3*E3) + R1)*(x^2*y + y^2*z + x*z^2) # = 0
R3*S^3 - 6*R2*R3*S + R2^3 + 9*R3^2 +
	- 1/(b+1)^2 * (b*(E2*S - 3*E3) + R1)*(E2*S - 3*E3 - R1) # = 0
(b+1)^2*R3*S^3 - (b+1)^2*6*R2*R3*S + (b+1)^2*R2^3 + 9*(b+1)^2*R3^2 +
	- (b*(E2*S - 3*E3) + R1)*(E2*S - 3*E3 - R1) # = 0
(b+1)^2*R3*S^3 - (b+1)^2*6*R2*R3*S + (b+1)^2*R2^3 + 9*(b+1)^2*R3^2 +
	- b*(R2*S - 3*R3)^2 + (b-1)*R1*(R2*S - 3*R3) + R1^2 # = 0
(b+1)^2*R3*S^3 - b*R2^2*S^2 - 6*(b+1)^2*R2*R3*S + 6*b*R2*R3*S + (b-1)*R1*R2*S +
	+ R1^2 + (b+1)^2*R2^3 + 9*((b+1)^2-b)*R3^2 - 3*(b-1)*R1*R3 # = 0


### Solution
solve.S3HtM.DualSum21 = function(R, a, b=0) {
	if(all(b == 0)) {
		coeff = c((a+1)^2*R[3], - a*R[2]^2, - 6*(a+1)^2*R[2]*R[3] + 6*a*R[2]*R[3] + (a-1)*R[1]*R[2],
			R[1]^2 + (a+1)^2*R[2]^3 + 9*((a+1)^2-a)*R[3]^2 - 3*(a-1)*R[1]*R[3])
	} else {
		# TODO
		print("Not yet implemented!")
	}
	S = roots(coeff)
	len = length(S)
	print(S)
	b2 = if(length(b) > 1) b[2] else 0; # Ext 2;
	x = sapply(S, function(x) roots(c(1, -x, R[2] - b2*x, -R[3])))
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	yz = R[3]/x
	yz.s = S - x
	### robust:
	R1A = (a*(R[2]*S - 3*R[3]) + R[1]) / (a+1) # TODO: a = -1
	R1B = ((R[2]*S - 3*R[3]) - R[1]) / (a+1) # TODO: a = -1
	if(any(round0(R1A - R1B) == 0)) print("Numerical instability!")
	yz.d = - (R1A - R1B) / (x^2 + yz - x*yz.s)
	y = (yz.s + yz.d)/2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}
test.S3HtM.P21 = function(sol, a, b.ext=0, R=NULL, tol=1E-8) {
	x = sol[,1]; y = sol[,2]; z = sol[,3];
	err1 = (x*y^2 + y*z^2 + z*x^2) - a*(x*z^2 + y*x^2 + z*y^2);
	err2 = x*y + x*z + y*z;
	err3 = x*y*z;
	err = rbind(err1, err2, err3);
	if( ! is.null(R)) {
		err = err - R;
	}
	err = round0(err, tol=tol);
	return(err);
}

### Example:
a = 1
R = c(1,1,1)
sol = solve.S3HtM.DualSum21(R, a=a)

test.S3HtM.P21(sol, a=a)

poly.calc0(x)
err = -1 + 3*x - 3.25*x^2 + 4.5*x^3 - 1.75*x^4 - x^5 + 4.5*x^6 - 1.5*x^7 - 0.25*x^8 + x^9
round0(err)

### Ex 2:
a = 2
R = c(1,3,-2)
sol = solve.S3HtM.DualSum21(R, a=a)

test.S3HtM.P21(sol, a=a)


### Special Case:
# Non-Oriented Ht:
# => all permutations are enabled!

# TODO:

R = c(0, 2, 3)
a = 1
sol = solve.S3HtM.DualSum21(R, a=a)

test.S3HtM.P21(sol, a=a)

### Debug
R = c(0, 2, 3)
a = 1

x =  1.089990536315 - 1.250695049316i;
y = -0.148968812161 + 1.079762780452i;
z =  1.089990536315 - 1.250695049316i;
sol = cbind(x, y, z)
# every permutation is valid;
sol = rbind(sol, sol[c(1,3,2)])
test.S3HtM.P21(sol, a=a)


###
source("Polynomials.Helper.Solvers.Num.R")

solve.S3HtM.Num = function(x, R, a=1) {
	x = matrix(x, ncol=3);
	xc = x[2,]; x = x[1,] + 1i*xc;
	x = matrix(x, nrow=1);
	y = test.S3HtM.P21(matrix(x, nrow=1), R=R, a=a, tol=1E-15);
	y = rbind(Re(y), Im(y));
	y = as.vector(y);
	return(y);
}

x0 = c(1.09-1.2507i, -0.149+1.0798i, 1.09-1.2507i);
x = solve.all(solve.S3HtM.Num, x0, R=R)


### Test
(x*y^2 + y*z^2 + z*x^2) - a*(x*z^2 + y*x^2 + z*y^2) # - R[1] # = 0
x*y + x*z + y*z # - R[2] # = 0
x*y*z # - R[3] # = 0


############################
############################

####################
### HtComponents ###
####################

x^2*y + b*z # = R
y^2*z + b*x # = R
z^2*x + b*y # = R

### Solution:

# - see file:
#   Poly.System.Hetero.Symmetric.S3.Leading.R;
# - section: Mixt-Order: 2+1;

### [main technique]
### Sum =>
(x^2*y + y^2*z + z^2*x) + b*S - 3*R # = 0
# "Inversion/Rotation" =>
E3*S^3 - (3*R - b*S + 6*E3)*E2*S + (3*R - b*S)^2 + E2^3 + 9*E3^2 + 3*(3*R - b*S)*E3 # = 0
E3*S^3 + 9*E3^2 - 3*b*E3*S + 9*R*E3 - 6*E2*E3*S + E2^3 + b*E2*S^2 - 3*R*E2*S +
	 + b^2*S^2 - 6*b*R*S + 9*R^2 # = 0


##################
##################
##################

##################
### Resonances ###
##################

######################
### Roots of Unity ###
######################

### Order: k + n
# - Base-Eq with terms: x^k*y^n, y^k*z^n, z^k*x^n
# - roots of unity: m^p = 1;
# - reusing: "x" = power of m in x;

k*x + n*y # = 0 (mod p);
k*y + n*z # = 0 (mod p);
k*z + n*x # = 0 (mod p);

### Sum =>
(k+n)*(x+y+z) # = 0 (mod p);

### Reduction =>
(n^2 + k^2 - n*k)*x # = 0 (mod p);

### Trivial Powers
# p = (k+n) # or one of its Divisors;

### Non-Trivial Powers
# p = (n^2 + k^2 - n*k) # or one of its Divisors;

### Combinations:
# p = Divisors of (k+n)*(n^2 + k^2 - n*k)

#############
### Examples:

### Note:
# - entire system has to be compatible with these "rotations",
#   if they should represent valid transformations;

#######
### Ex:
# k = 2; n = 1;
### Trivial:
p = 3;
### Non-Trivial & Combinations:
p = c(3, 9); # p[non-trivial] is the same!
# ex: p = 9 =>
# new solution: (x,y,z) * (m^1, m^7, m^4);


#######
### Ex:
# k = 3; n = 1;
### Trivial:
p = 4;
### Non-Trivial & Combinations:
p = c(7, 14, 28);
# ex: p = 7 =>
# new solution: (x,y,z) * (m^1, m^4, m^2);
# ex: p = 14 =>
# new solution: (x,y,z) * (m^1, m^11, m^9);
# ex: p = 28 =>
# new solution: (x,y,z) * (m^1, m^25, m^9);


#######
### Ex:
# k = 3; n = 2;
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(7, 35);
# ex: p = 7 =>
# new solution: (x,y,z) * (m^1, m^2, m^4);
# ex: p = 35 =>
# new solution: (x,y,z) * (m^1, m^16, m^11);


######################
### Roots of Unity ###
### Generalization ###
######################

### System with [i] Variables
# - for i = all odd;
### Order: k + n
# - Base-Eq with terms:
#   x[1]^k*x[2]^n, x[2]^k*x[3]^n, ..., x[i-1]^k*x[i]^n, x[i]^k*x[1]^n;

### All Powers
# p = Divisors of (k^i + n^i)


#######
### Ex:
# i = 5; k = 2; n = 1;
# p = Divisors(33);
### Trivial:
p = 3;
### Non-Trivial & Combinations:
p = c(11, 33);
# ex: p = 11 =>
# new solution:
# (x1,x2,x3,x4,x5) * (m^1, m^9, m^4, m^3, m^5);


#######
### Ex:
# i = 5; k = 3; n = 1;
# p = Divisors(244);
### Trivial:
p = 4;
### Non-Trivial & Combinations:
p = c(61, 122, 244);
# ex: p = 61 =>
# new solution:
# (x1,x2,x3,x4,x5) * (m^1, m^58, m^9, m^34, m^20);


#######
### Ex:
# i = 5; k = 3; n = 2;
# p = Divisors(275);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(11, 25, 55, 275);
# ex: p = 11 =>
# new solution:
# (x1,x2,x3,x4,x5) * (m^1, m^4, m^5, m^9, m^3);
# ex: p = 25 =>
# new solution:
# (x1,x2,x3,x4,x5) * (m^1, m^11, m^21, m^6, m^16);


###################

### System with [i] Variables
# - for i = even;
### Order: k + n
# - Base-Eq with terms:
#   x[1]^k*x[2]^n, x[2]^k*x[3]^n, ..., x[i-1]^k*x[i]^n, x[i]^k*x[1]^n;

### All Powers
# p = Divisors of (k^i - n^i)


#######
### Ex:
# i = 4; k = 2; n = 1;
# p = Divisors(15);
### Trivial:
p = 3;
### Non-Trivial & Combinations:
p = c(5, 15);
# ex: p = 5 =>
# new solution:
# (x1,x2,x3,x4) * (m^1, m^3, m^4, m^2);


#######
### Ex:
# i = 4; k = 3; n = 1;
# p = Divisors(80);
### Trivial:
p = 4;
### Non-Trivial & Combinations:
p = c(5, 8, 10); # and higher
# ex: p = 5 =>
# new solution:
# (x1,x2,x3,x4) * (m^1, m^2, m^4, m^3);
# ex: p = 8 =>
# new solution:
# (x1,x2,x3,x4) * (m^1, m^5, m^1, m^5);


#######
### Ex:
# i = 4; k = 3; n = 2;
# p = Divisors(65);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(13, 65);
# ex: p = 13 =>
# new solution:
# (x1,x2,x3,x4) * (m^1, m^5, m^12, m^8);


#######
### Ex:
# i = 4; k = 4; n = 1;
# p = Divisors(255);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(3, 17); # and higher
# ex: p = 17 =>
# new solution:
# (x1,x2,x3,x4) * (m^1, m^13, m^16, m^4);


#######
### Ex:
# i = 6; k = 2; n = 1;
resonance(c(2,1), n=6)
# p = Divisors(63);
### Trivial:
p = 3;
### Non-Trivial & Combinations:
p = c(7, 9, 21, 63);
# ex: p = 7 =>
# new solution:
# (x1,x2,...,x6) * (m^1, m^5, m^4, m^6, m^2, m^3);
# ex: p = 9 =>
# new solution:
# (x1,x2,...,x6) * (m^1, m^7, m^4, m^1, m^7, m^4);


####################

################
### 3 Powers ###
################

## Order (a,b,c)
# a^2*c^4 - 2*a^4*c^2 + 4*b^2*a^3*c - b^4*a^2 + a^6


### Order: 2+1+1
### 4 Variables
# i = 4; p = c(2,1,1)
# p = Divisors(2^6);
### Trivial:
p = 4;
### Non-Trivial & Combinations:
p = c();
# - there are various quasi-non-trivial solutions;
### Ex: p = 8 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^3, m^7, m^3, m^7);
# (3,7,3), (7,3,7), (3,7,3), (7,3,7)


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
# (3,2,4), (2,4,0), (4,0,3), (0,3,2)
### Ex: p = 15 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^1, m^11, m^1, m^11);
# (1,11,1), (11,1,11), (1,11,1), (11,1,11)


### Order: 4+1+1
### 4 Variables
# i = 4; p = c(4,1,1)
# p = Divisors(2^8*3*5);
### Trivial:
p = 6;
### Non-Trivial & Combinations:
p = c();
# - some of the computed ones do NOT work;
# - but there are various quasi-non-trivial solutions;
### Ex: p = 8 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^5, m^7, m^5, m^7);
# (5,7,5), (7,5,7), (5,7,5), (7,5,7)


### Order: (p1,1,1)

### 5 Variables:
# p = Divisors(p1^3*(p1-1)*(2*p1^2 + 3))
# - but condition is often insufficient;

### Order: 2+1+1
### 5 Variables
# i = 5; p = c(2,1,1)
# p = Divisors(88);
### Trivial:
p = 4;
### Non-Trivial & Combinations:
p = c(11, 22); # & higher
# ex: p = 11 =>
# new solution:
# (x1,x2,...,x5) * (m^1, m^4, m^5, m^9, m^3);
# (1,4,5), (4,5,9), (5,9,3), (9,3,1), (3,1,4)
# ex: p = 22 =>
# new solution:
# (x1,x2,...,x5) * (m^1, m^15, m^5, m^9, m^3);
# (1,15,5), (15,5,9), (5,9,3), (9,3,1), (3,1,15)


### Order: 3+1+1
### 6 Variables
# i = 6; p = c(3,1,1)
# p = Divisors(???);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(9); # possible others
# ex: p = 9 =>
# new solution:
# (x1,x2,...,x6) * (m^1, m^2, m^4, m^8, m^7, m^5);
# (1,2,4), (2,4,8), (4,8,7), (8,7,5), (7,5,1), (5,1,2)

