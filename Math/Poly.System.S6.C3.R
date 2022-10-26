########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S6: C3-Hetero-Symmetric
### Basic Types
###
### draft v.0.1b


####################

### Helper Functions


source("Polynomials.Helper.R")


#####################
#####################

##############
### Theory ###
##############

### Basic Systems:

### Eqs 1 & 2:
# x1^n1 + x2^n1 + x3^n1 = R1
# y1^n2 + y2^n2 + y3^n2 = R2

### A.) Type 2-Ht:
# x1*y1 + x2*y2 + x3*y3 = R3
# x1*y2 + x2*y3 + x3*y1 = R4
# x1*x2*x3 = R5
# y1*y2*y3 = R6

### Symmetries:
# - if {(x1, x2, x3), (y1, y2, y3)} is a solution,
#   so are also all concordant C3-cyclic permutations:
#   e.g. {(x2, x1, x3), (y2, y1, y3)};


### B.) Type 1-Ht & E2:
# - has additional symmetries;
# x1*y1 + x2*y2 + x3*y3 = R3
# E2x + E2y = R4
# x1*x2*x3 = R5
# y1*y2*y3 = R6
### Variant Eq 4:
# E2x * E2y = R4


### C.) Type 1-Ht & 2 x E2:
# - has additional symmetries;
# x1*y1 + x2*y2 + x3*y3 = R3
# E2x = R4
# E2y = R5
# x1*x2*x3*y1*y2*y3 = R6


###################
###################

#################
### Type 2-Ht ###
#################

### Basic System: Eqs 3 & 4
# x1*y1 + x2*y2 + x3*y3 = R3
# x1*y2 + x2*y3 + x3*y1 = R4


### Solution:

### Sum:
A1 + B1 + C1 - sx*sy # = 0
# C1 = sx*sy - A1 - B1;

### E2:
E2ABC + 3*E2x*E2y - E2y*sx^2 - E2x*sy^2 # = 0

### E3:
A1*B1*C1 - py*(sx^3 - 3*E2x*sx) - px*(sy^3 - 3*E2y*sy) - 18*px*py +
	- 2*A2a*B2a + (E2x*sx - 3*px)*B2a + (E2y*sy - 3*py)*A2a +
	- E2x*sx*E2y*sy + 3*E2x*py*sx + 3*E2y*px*sy # = 0
# with:
A2a^2 - (E2x*sx - 3*px)*A2a - (6*E2x*px*sx - E2x^3 - 9*px^2 - px*sx^3) # = 0
B2a^2 - (E2y*sy - 3*py)*B2a - (6*E2y*py*sy - E2y^3 - 9*py^2 - py*sy^3) # = 0


### Solver:

solver.S6C3.P1Ht2 = function(R, debug=TRUE, all=FALSE) {
	coeff = coeff.S6C3.P1Ht2(R);
	E2x = roots(coeff);
	if(debug) print(E2x);
	#
	sx = R[1]; sy = R[2]; px = R[5]; py = R[6];
	A1 = R[3]; B1 = R[4];
	C1 = sx*sy - A1 - B1;
	E2ABC = (A1 + B1)*C1 + A1*B1;
	E2y = (sy^2*E2x - E2ABC) / (3*E2x - sx^2);
	#
	len = length(E2x);
	x1 = sapply(seq(len), function(id) {
		roots(c(1, -sx, E2x[id], -px));
	})
	# Robust: x2 & x3 are exactly determined;
	# Note: all 12 solutions valid when A1 == B1;
	if(round0(A1 - B1) == 0) {
		tmp = x1[c(1,3,2), ];
		x1  = cbind(x1, tmp);
		E2x = rep(E2x, 2); E2y = rep(E2y, 2);
		x2  = as.vector(x1[2,]); x3 = as.vector(x1[3,]); x1 = as.vector(x1[1,]);
	}  else {
		x1 = as.vector(x1);
		E2x = rep(E2x, each=3); E2y = rep(E2y, each=3);
		x23 = sx - x1; # px23 = E2x - x1*x23;
		len = length(x1);
		x2 = sapply(seq(len), function(id) {
			lst = list(x1=x1[id], x23=x23[id], E2x=E2x[id], E2y=E2y[id],
				A1=A1, B1=B1, C1=C1, sx=sx, sy=sy, px=px, py=py);
			x2 = x2.S6C3.P1Ht2(lst);
		})
		x3 = x23 - x2;
	}
	y1  = (x1*A1 + x3*B1 + x2*C1 - E2x*sy) / (sx^2 - 3*E2x);
	y23 = sy - y1;
	y2  = (x3*y23 - A1 + x1*y1) / (x3 - x2);
	y3  = y23 - y2;
	sol = cbind(x1, x2, x3, y1, y2, y3);
	if(all) {
		# TODO
	}
	return(sol);
}
x2.S6C3.P1Ht2 = function(vars) {
	x1 = vars$x1; x23 = vars$x23; E2x = vars$E2x; E2y = vars$E2y;
	A1 = vars$A1; B1 = vars$B1; C1 = vars$C1;
	sx = vars$sx; sy = vars$sy; px = vars$px; py = vars$py;
	cc = c((C1 - B1)^3, (C1 - B1)^2*(3*(A1*x1 + B1*x23) - sy*sx^2),
		3*sy^2*E2x^2*(B1 - C1) + 2*sx^2*sy^2*E2x*(C1 - B1) +
			+ 2*(B1 - C1)*sx^2*sy*(A1*x1 + B1*x23) +
			+ 3*(C1 - B1)*(A1*x1 + B1*x23)^2 +
			+ (C1 - B1)*(3*E2x - sx^2)^2*E2y,
		2*sy^3*E2x^3 - sx^2*sy^3*E2x^2 + sx^2*sy^2*E2x*(A1*x1 + B1*x23) +
			- sy^2*(3*E2x - sx^2)*(A1*x1 + B1*x23)*E2x +
			+ (A1*x1 + B1*x23)^3 - sx^2*sy*(A1*x1 + B1*x23)^2 +
			+ (3*E2x - sx^2)^3*py +
			+ (3*E2x - sx^2)^2*(A1*x1 + B1*x23 - sy*E2x)*E2y
	);
	cc = rev(cc);
	if(round0(B1 - C1) == 0) {
		warning("Special Case!");
		return(NA); # TODO
	} else {
		cc = cc / cc[4];
		cc = cc[-4];
		cc = cc - c(-px, E2x, -sx);
	}
	if(round0(cc[3]) == 0) return( - cc[1] / cc[2]);
	cc = cc / cc[3];
	dd = c(px, - E2x + cc[1], sx + cc[2]);
	if(round0(dd[3]) == 0) return( - dd[1] / dd[2]);
	dd = dd / dd[3];
	cc = cc - dd;
	if(round0(cc[2]) == 0) {
		warning("Something went wrong in the robust computation of x2!");
		return(NA);
	}
	return( - cc[1] / cc[2]);
}
coeff.S6C3.P1Ht2 = function(R) {
	sx = R[1]; sy = R[2]; px = R[5]; py = R[6];
	A1 = R[3]; B1 = R[4];
	C1 = sx*sy - A1 - B1;
	E2ABC = (A1 + B1)*C1 + A1*B1;
	coeff = c(sy^6 - 54*sy^3*py + 729*py^2,
		- 6*sy^4*E2ABC + 54*sy^3*sx^2*py + 162*sy*E2ABC*py - 1458*sx^2*py^2,
		sy^4*sx^2*E2ABC + 9*sy^2*E2ABC^2 - 18*sy^3*sx^4*py - 189*sy*sx^2*E2ABC*py + 1215*sx^4*py^2 +
			+ 9*sy^3*sx*C1*B1*A1 - 243*sx*py*C1*B1*A1,
		- 2*sy^2*sx^2*E2ABC^2 - 4*E2ABC^3 - 2*sy^6*sx^3*px + 9*sy^4*sx*E2ABC*px + 2*sy^3*sx^6*py +
			+ 81*sy*sx^4*E2ABC*py + 54*sy^3*sx^3*px*py - 243*sy*sx*E2ABC*px*py - 540*sx^6*py^2 +
			- 6*sy^3*sx^3*C1*B1*A1 - 9*sy*sx*E2ABC*C1*B1*A1 - 27*sy^3*px*C1*B1*A1 + 297*sx^3*py*C1*B1*A1 +
			+ 729*px*py*C1*B1*A1 - 27*C1^2*B1^2*A1^2,
		sx^2*E2ABC^3 + 6*sy^4*sx^3*E2ABC*px - 27*sy^2*sx*E2ABC^2*px - 15*sy*sx^6*E2ABC*py - 54*sy^3*sx^5*px*py +
			+ 243*sy*sx^3*E2ABC*px*py + 135*sx^8*py^2 + sy^3*sx^5*C1*B1*A1 + 6*sy*sx^3*E2ABC*C1*B1*A1 +
			+ 81*sy*E2ABC*px*C1*B1*A1 - 135*sx^5*py*C1*B1*A1 - 729*sx^2*px*py*C1*B1*A1 + 27*sx^2*C1^2*B1^2*A1^2,
		- sy^4*sx^5*E2ABC*px + 18*sx*E2ABC^3*px + sy*sx^8*E2ABC*py + 18*sy^3*sx^7*px*py - 81*sy*sx^5*E2ABC*px*py +
			- 18*sx^10*py^2 - sy*sx^5*E2ABC*C1*B1*A1 + 9*sy^3*sx^4*px*C1*B1*A1 - 54*sy*sx^2*E2ABC*px*C1*B1*A1 +
			+ 27*sx^7*py*C1*B1*A1 + 243*sx^4*px*py*C1*B1*A1 - 9*sx^4*C1^2*B1^2*A1^2,
		sy^2*sx^5*E2ABC^2*px - 4*sx^3*E2ABC^3*px + sy^6*sx^6*px^2 - 9*sy^4*sx^4*E2ABC*px^2 +
			+ 27*sy^2*sx^2*E2ABC^2*px^2 - 27*E2ABC^3*px^2 - 2*sy^3*sx^9*px*py + 9*sy*sx^7*E2ABC*px*py +
			+ sx^12*py^2 - 2*sy^3*sx^6*px*C1*B1*A1 + 9*sy*sx^4*E2ABC*px*C1*B1*A1 - 2*sx^9*py*C1*B1*A1 +
			- 27*sx^6*px*py*C1*B1*A1 + sx^6*C1^2*B1^2*A1^2
	);
	return(coeff);
}
### Test:
test.S6C3.P1Ht2 = function(sol, R=NULL) {
	x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3];
	y1 = sol[,4]; y2 = sol[,5]; y3 = sol[,6];
	
	err1 = x1 + x2 + x3;
	err2 = y1 + y2 + y3;
	err3 = x1*y1 + x2*y2 + x3*y3;
	err4 = x1*y2 + x2*y3 + x3*y1;
	err5 = x1*x2*x3;
	err6 = y1*y2*y3;
	err = rbind(err1, err2, err3, err4, err5, err6);
	if( ! is.null(R)) err = err - R;
	err = round0(err);
	return(err);
}

### Examples

R = c(2,3,-1,4,5,6)
sol = solver.S6C3.P1Ht2(R)

test.S6C3.P1Ht2(sol)


### Ex 2: Special Case
R = c(2,3,-2,-2,5,6)
sol = solver.S6C3.P1Ht2(R)

test.S6C3.P1Ht2(sol)


### Ex 3: Special Case
R = c(-1,3,-1,-1,5,6)
sol = solver.S6C3.P1Ht2(R)

test.S6C3.P1Ht2(sol)


### Ex 4:
R = c(-1,5,-1,2,3,-1)
sol = solver.S6C3.P1Ht2(R)

test.S6C3.P1Ht2(sol)


### Test:
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3];
y1 = sol[,4]; y2 = sol[,5]; y3 = sol[,6];

x1 + x2 + x3 # = R1
y1 + y2 + y3 # = R2
x1*y1 + x2*y2 + x3*y3 # = R3
x1*y2 + x2*y3 + x3*y1 # = R4
x1*x2*x3 # = R5
y1*y2*y3 # = R6


### Debug

R = c(2,3,-1,4,5,6)
x1 = -0.0751764371 - 1.2207329994i;
x2 = -0.3636407258 + 1.6144160902i;
x3 =  2.4388171629 - 0.3936830909i;
y1 =  0.1828061717 + 1.2261038046i;
y2 =  2.8373554636 + 0.4578248858i;
y3 = -0.0201616352 - 1.6839286904i;


sx = x1 + x2 + x3; sy = y1 + y2 + y3;
px = x1 * x2 * x3; py = y1 * y2 * y3;
E2x = (x1 + x2)*x3 + x1*x2;
E2y = (y1 + y2)*y3 + y1*y2;
#
A1 = x1*y1 + x2*y2 + x3*y3;
B1 = x1*y2 + x2*y3 + x3*y1;
C1 = x1*y3 + x2*y1 + x3*y2;
E2ABC = (A1 + B1)*C1 + A1*B1;
# [intermediary]
A2a = x1^2*x2 + x2^2*x3 + x3^2*x1;
A2b = x1*x2^2 + x2*x3^2 + x3*x1^2;
B2a = y1^2*y2 + y2^2*y3 + y3^2*y1;
B2b = y1*y2^2 + y2*y3^2 + y3*y1^2;


### Derivation:
p2 = toPoly.pm(...); # see formula for E3 above;
p1 = toPoly.pm("E2ABC + 3*E2x*E2y - E2y*sx^2 - E2x*sy^2");
pA2a = toPoly.pm("A2a^2 - (E2x*sx - 3*px)*A2a - (6*E2x*px*sx - E2x^3 - 9*px^2 - px*sx^3)");
pB2a = toPoly.pm("B2a^2 - (E2y*sy - 3*py)*B2a - (6*E2y*py*sy - E2y^3 - 9*py^2 - py*sy^3)");

# Note: B2a cancels out:
pR = solve.lpm(p2, pA2a, pB2a, xn=c("A2a", "B2a"))
pR = drop.pm(pR[[2]]$Rez)

pR2 = solve.pm(p1, pR, "E2y")
str(pR2)
toCoeff(pR2$Rez, "E2x")


###
x1*A1 + x3*B1 + x2*C1 - y1*(sx^2 - 2*E2x) - E2x*(sy - y1) # = 0
x1*A1 + x3*B1 + x2*C1 - y1*(sx^2 - 3*E2x) - E2x*sy # = 0

y1*(sx - 2*x1) + x2*y3 + x3*y2 + x1*sy - B1 - C1 # = 0

### x3*A1 - x1*B1 =>
(x2*x3 - x1^2)*y2 + (x3^2 - x1*x2)*y3 - A1*x3 + B1*x1 # = 0
(x1^2 - x2*x3)*y1 - (x2^2 - x1*x3)*y3 - A1*x1 + B1*x2 # = 0

###
pYA = toPoly.pm("x1*y1 + x2*y2 + (x23 - x2)*(sy - y1 - y2) - A1")
pYB = toPoly.pm("x1*y2 + x2*(sy - y1 - y2) + (x23 - x2)*y1 - B1")
pY1 = toPoly.pm("x1*A1 + (x23 - x2)*B1 + x2*C1 - y1*(sx^2 - 3*E2x) - E2x*sy");
pY2 = toPoly.pm("y1^3 - sy*y1^2 + E2y*y1 - py");
pY = solve.pm(pYA, pYB, "y2")$Rez; # redundant;
pR = solve.pm(pY1, pY2, "y1")
pR = pR$Rez;
pX = toPoly.pm("x2^3 - sx*x2^2 + E2x*x2 - px")
pR2 = solve.pm(pR, pX, "x2", stop.at=1)
table(pR2[[2]]$x2)


#########################
#########################

### Variant: E2x + E2y

### Basic System: Eqs 5 & 6
# E2x + b*E2y = R5
# x1*x2*x3*y1*y2*y3 = R6

### Solution:

### Sum:
A1 + B1 + C1 - sx*sy # = 0
# C1 = sx*sy - A1 - B1;

### E2:
E2ABC + 3*E2x*E2y - E2y*sx^2 - E2x*sy^2 # = 0
# =>
b*E2ABC + 3*E2x*(E2s - E2x) - (E2s - E2x)*sx^2 - b*E2x*sy^2 # = 0
3*E2x^2 - (3*E2s + sx^2 - b*sy^2)*E2x - b*E2ABC + E2s*sx^2 # = 0


### Solver:

solve.S6C3.P1Ht2VE2 = function(R, b, debug=TRUE, all=FALSE) {
	sx = R[1]; sy = R[2]; E2s = R[5]; E3 = R[6];
	A1 = R[3]; B1 = R[4];
	C1 = sx*sy - A1 - B1;
	E2ABC = (A1 + B1)*C1 + A1*B1;
	coeff = c(3, - (3*E2s + sx^2 - b*sy^2), - b*E2ABC + E2s*sx^2);
	E2x = roots(coeff);
	E2y = (E2s - E2x) / b;
	if(debug) print(E2x);
	# px, py:
	len = length(E2x);
	px = sapply(seq(len), function(id) {
		coeff = coeff.S6C3.P1Ht2VE2(R, list(E2x=E2x[id], E2y=E2y[id]));
		roots(coeff);
	});
	px  = as.vector(px);
	E2x = rep(E2x, each=4);
	E2y = rep(E2y, each=4);
	py  = E3 / px;
	if(debug) print(px);
	#
	len = length(px);
	x1 = sapply(seq(len), function(id) {
		roots(c(1, -sx, E2x[id], -px[id]));
	})
	# Robust: x2 & x3 are exactly determined;
	# Note: all 24 solutions valid when A1 == B1;
	if(round0(A1 - B1) == 0) {
		if(debug) print("Special Case:")
		tmp = x1[c(1,3,2), ];
		x1  = cbind(x1, tmp);
		E2x = rep(E2x, 2); E2y = rep(E2y, 2);
		px  = rep(px, 2); py = rep(py, 2);
		x2  = as.vector(x1[2,]); x3 = as.vector(x1[3,]); x1 = as.vector(x1[1,]);
	}  else {
		x1 = as.vector(x1);
		E2x = rep(E2x, each=3); E2y = rep(E2y, each=3);
		px  = rep(px, each=3); py = rep(py, each=3);
		x23 = sx - x1; # px23 = E2x - x1*x23;
		len = length(x1);
		x2 = sapply(seq(len), function(id) {
			lst = list(x1=x1[id], x23=x23[id], E2x=E2x[id], E2y=E2y[id],
				A1=A1, B1=B1, C1=C1, sx=sx, sy=sy, px=px[id], py=py[id]);
			x2 = x2.S6C3.P1Ht2VE2(lst);
		})
		x3 = x23 - x2;
	}
	y1  = (x1*A1 + x3*B1 + x2*C1 - E2x*sy) / (sx^2 - 3*E2x);
	y23 = sy - y1;
	y2  = (x3*y23 - A1 + x1*y1) / (x3 - x2);
	y3  = y23 - y2;
	sol = cbind(x1, x2, x3, y1, y2, y3);
	if(all) {
		# TODO
	}
	return(sol);
}
x2.S6C3.P1Ht2VE2 = function(vars) {
	# same as: x2.S6C3.P1Ht2();
	x1 = vars$x1; x23 = vars$x23; E2x = vars$E2x; E2y = vars$E2y;
	A1 = vars$A1; B1 = vars$B1; C1 = vars$C1;
	sx = vars$sx; sy = vars$sy; px = vars$px; py = vars$py;
	cc = c((C1 - B1)^3, (C1 - B1)^2*(3*(A1*x1 + B1*x23) - sy*sx^2),
		3*sy^2*E2x^2*(B1 - C1) + 2*sx^2*sy^2*E2x*(C1 - B1) +
			+ 2*(B1 - C1)*sx^2*sy*(A1*x1 + B1*x23) +
			+ 3*(C1 - B1)*(A1*x1 + B1*x23)^2 +
			+ (C1 - B1)*(3*E2x - sx^2)^2*E2y,
		2*sy^3*E2x^3 - sx^2*sy^3*E2x^2 + sx^2*sy^2*E2x*(A1*x1 + B1*x23) +
			- sy^2*(3*E2x - sx^2)*(A1*x1 + B1*x23)*E2x +
			+ (A1*x1 + B1*x23)^3 - sx^2*sy*(A1*x1 + B1*x23)^2 +
			+ (3*E2x - sx^2)^3*py +
			+ (3*E2x - sx^2)^2*(A1*x1 + B1*x23 - sy*E2x)*E2y
	);
	cc = rev(cc);
	if(round0(B1 - C1) == 0) {
		warning("Special Case!");
		return(NA); # TODO
	} else {
		cc = cc / cc[4];
		cc = cc[-4];
		cc = cc - c(-px, E2x, -sx);
	}
	if(round0(cc[3]) == 0) return( - cc[1] / cc[2]);
	cc = cc / cc[3];
	dd = c(px, - E2x + cc[1], sx + cc[2]);
	if(round0(dd[3]) == 0) return( - dd[1] / dd[2]);
	dd = dd / dd[3];
	cc = cc - dd;
	if(round0(cc[2]) == 0) {
		warning("Something went wrong in the robust computation of x2!");
		return(NA);
	}
	return( - cc[1] / cc[2]);
}
coeff.S6C3.P1Ht2VE2 = function(R, E2) {
	sx = R[1]; sy = R[2]; E3 = R[6];
	A1 = R[3]; B1 = R[4]; E2x = E2$E2x; E2y = E2$E2y;
	C1 = sx*sy - A1 - B1; ABC = A1*B1*C1;
	coeff = c((sy^2 - 3*E2y)^3,
		sx^3*sy^2*E2y^2 - 4*sx^3*E2y^3 + sx*sy^4*E2x*E2y - 9*sx*sy^2*E2x*E2y^2 + 18*sx*E2x*E2y^3 +
			+ ABC*sy*(9*E2y - 2*sy^2),
		- 2*E3*sx^3*sy^3 + 9*E3*sx^3*sy*E2y + 9*E3*sx*sy^3*E2x - 27*E3*sx*sy*E2x*E2y + sx^2*E2x^2*E2y^3 +
			+ sy^2*E2x^3*E2y^2 - 4*E2x^3*E2y^3 - 27*ABC*E3 - ABC*sx*sy*E2x*E2y + ABC^2,
		(sx^4*sy*E2x*E2y + sx^2*sy^3*E2x^2 - 9*sx^2*sy*E2x^2*E2y - 4*sy^3*E2x^3 +
			+ 18*sy*E2x^3*E2y - 2*ABC*sx^3 + 9*ABC*sx*E2x)*E3,
		E3^2*(sx^2 - 3*E2x)^3
	);
	return(coeff);
}
### Test
test.S6C3.P1Ht2VE2 = function(sol, b, R=NULL) {
	x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3];
	y1 = sol[,4]; y2 = sol[,5]; y3 = sol[,6];
	#
	err1 = x1 + x2 + x3;
	err2 = y1 + y2 + y3;
	err3 = x1*y1 + x2*y2 + x3*y3;
	err4 = x1*y2 + x2*y3 + x3*y1;
	err5 = (x1+x2)*x3 + x1*x2 + b*(y1+y2)*y3 + b*y1*y2;
	err6 = x1*x2*x3*y1*y2*y3;
	err = rbind(err1, err2, err3, err4, err5, err6);
	if( ! is.null(R)) err = err - R;
	err = round0(err);
	return(err);
}

### Examples:

### Ex 1:
R = c(2,-1,3,4,5,-3)
b = 2
sol = solve.S6C3.P1Ht2VE2(R, b=b);

test.S6C3.P1Ht2VE2(sol, b=b)


### Ex 2:
R = c(2,4,-5,4,-1,-2)
b = -1
sol = solve.S6C3.P1Ht2VE2(R, b=b);

test.S6C3.P1Ht2VE2(sol, b=b)


### Ex 3: Special Case
R = c(2,-1,3,3,5,-3)
b = -2
sol = solve.S6C3.P1Ht2VE2(R, b=b);

test.S6C3.P1Ht2VE2(sol, b=b)


#########
### Test:
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3];
y1 = sol[,4]; y2 = sol[,5]; y3 = sol[,6];

x1 + x2 + x3 # = R1
y1 + y2 + y3 # = R2
x1*y1 + x2*y2 + x3*y3 # = R3
x1*y2 + x2*y3 + x3*y1 # = R4
(x1 + x2)*x3 + x1*x2 + b*(y1 + y2)*y3 + b*y1*y2 # = R5
x1*x2*x3*y1*y2*y3 # = R6


### Debug

R = c(2,3,-1,4,5,6); b = 1;
x1 = -0.2816739496 + 0.6599722415i;
x2 =  2.1005351736 + 0.9357918602i;
x3 =  0.1811387760 - 1.5957641017i;
y1 = -0.0595357591 - 0.5027349019i;
y2 =  0.6616831651 + 1.5701963341i;
y3 =  2.3978525940 - 1.0674614321i;


sx = x1 + x2 + x3; sy = y1 + y2 + y3;
px = x1 * x2 * x3; py = y1 * y2 * y3;
E2x = (x1 + x2)*x3 + x1*x2;
E2y = (y1 + y2)*y3 + y1*y2;
#
A1 = x1*y1 + x2*y2 + x3*y3;
B1 = x1*y2 + x2*y3 + x3*y1;
C1 = x1*y3 + x2*y1 + x3*y2;
E2ABC = (A1 + B1)*C1 + A1*B1;
# [intermediary]
A2a = x1^2*x2 + x2^2*x3 + x3^2*x1;
A2b = x1*x2^2 + x2*x3^2 + x3*x1^2;
B2a = y1^2*y2 + y2^2*y3 + y3^2*y1;
B2b = y1*y2^2 + y2*y3^2 + y3*y1^2;


### Derivation:
p2 = toPoly.pm("A1*B1*C1 - py*(sx^3 - 3*E2x*sx) - px*(sy^3 - 3*E2y*sy) - 18*px*py +
	- 2*A2a*B2a + (E2x*sx - 3*px)*B2a + (E2y*sy - 3*py)*A2a +
	- E2x*sx*E2y*sy + 3*E2x*py*sx + 3*E2y*px*sy");
p1 = toPoly.pm("px*py - E3");
pA2a = toPoly.pm("A2a^2 - (E2x*sx - 3*px)*A2a - (6*E2x*px*sx - E2x^3 - 9*px^2 - px*sx^3)");
pB2a = toPoly.pm("B2a^2 - (E2y*sy - 3*py)*B2a - (6*E2y*py*sy - E2y^3 - 9*py^2 - py*sy^3)");

# Note: B2a cancels out:
pR = solve.lpm(p2, pA2a, pB2a, xn=c("A2a", "B2a"))
pR = drop.pm(pR[[2]]$Rez)
pR = solve.pm(p1, pR, "py")
toCoeff(pR$Rez, "px")

