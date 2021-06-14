########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Mixed Heterogeneous Symmetric
###
### draft v.0.1b-structure


###############

###############
### History ###
###############


### draft v.0.1b - v.0.1b-structure:
# - exploration of Rotations;
# - structure of solutions; [v.0.1b-structure]
### draft v.0.1a:
# - moved section on Mixed S4 systems
#   from: Poly.System.Hetero.Symmetric.S4.R
#   to this new file;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R
# - e.g. round0(), round0.p(),
#   solve.EnAll(), solveEn();


############################
############################

############################
### Mixed Hetero-Systems ###
############################

### x^2 + y^2 = R1
### x1*y1 + a2*x1*y2 + a3*y1*x2 = R2


# x1^2 + y1^2 = R1
# x2^2 + y2^2 = R1
# x1*y1 + a2*x1*y2 + a3*x2*y1 = R2
# x2*y2 + a2*x2*y1 + a3*x1*y2 = R2

### Solution:

### Special Case:
# x1 = x2; y1 = y2;
x^2 + y^2 - R1 # = 0
(a2 + a3 + 1)*x*y - R2 # = 0
# x1 = -x2; y1 = -y2;
x^2 + y^2 - R1 # = 0
(a2 + a3 - 1)*x*y + R2 # = 0

S^2 - 2*x*y - R1 # = 0

### Case:
# x1 != +/- x2; y1 != +/- y2;
# TODO:
# - are there any such solutions ???

### Solution
solve.S4P2 = function(R, a) {
	solve.part = function(div) {
		xy = R[2] / div;
		S = roots(c(1, 0, -2*xy - R[1]));
		xy = rep(xy, each=2);
		xy.diff = sqrt(R[1] - 2*xy + 0i)
		x = (S + xy.diff) / 2;
		y = (S - xy.diff) / 2;
		sol = cbind(x, y); sol = rbind(sol, sol[,2:1])
	}
	a.s = a[1] + a[2];
	sol2 = solve.part(a.s + 1);
	sol = cbind(sol2, sol2);
	# - R[2]
	R[2] = - R[2];
	sol2 = solve.part(a.s - 1);
	sol = rbind(sol, cbind(sol2, -sol2));
	return(sol)
}

### Examples
R = c(-1, 2)
a = c(3,4)
sol = solve.S4P2(R, a);
x1 = sol[,1]; x2 = sol[,3]; y1 = sol[,2]; y2 = sol[,4];

### Test
x1^2 + y1^2 # - R[1]
x2^2 + y2^2 # - R[1]
x1*y1 + a[1]*x1*y2 + a[2]*x2*y1 # - R[2]
x2*y2 + a[1]*x2*y1 + a[2]*x1*y2 # - R[2]


######################

### Simple P3

### x^3 + y^3 = R1
### a1*x1*y2 + a2*y1*x2 = R2

# x1^3 + y1^3 = R1
# x2^3 + y2^3 = R1
# a1*x1*y2 + a2*x2*y1 = R2
# a1*x2*y1 + a2*x1*y2 = R2

### Solution:

### Diff: Eq 3 - 4
(a1 - a2)*(x1*y2 - x2*y1) # = 0
# assumption: a1 != a2
x1*y2 - x2*y1 # = 0
# x1*y2 = x2*y1 = R2 / (a1+a2);

### Mult Eq 1r * 2r
x1^3*y2^3 - y1^3*x2^3 + R1*(x2^3 + y1^3) - R1^2 # = 0
R1*(x2^3 + y1^3) - R1^2 # = 0
# x2^3 + y1^3 = R1;
# x1^3 + y2^3 = R1;
# =>
(x1+y2)^3 - 3*x1*y2*(x1+y2) - R1 # = 0
S12^3 - 3*R2/(a1+a2)*S12 - R1 # = 0
(a1+a2)*S12^3 - 3*R2*S12 - (a1+a2)*R1 # = 0

# =>
(x2+y1)^3 - 3*x2*y1*(x2+y1) - R1 # = 0
S21^3 - 3*R2/(a1+a2)*S21 - R1 # = 0
(a1+a2)*S21^3 - 3*R2*S21 - (a1+a2)*R1 # = 0

### Solver:
solve.simple.S4P3 = function(R, a, debug=TRUE) {
	as = a[1] + a[2];
	S12 = roots(c(as, 0, - 3*R[2], -as*R[1]))
	S21 = roots(c(as, 0, - 3*R[2], -as*R[1]))
	if(debug) print(S12);
	solve.S2 = function(S, xy) {
		xy.diff = sqrt(S^2 - 4*xy + 0i);
		x = (S + xy.diff) / 2;
		y = (S - x);
		sol = cbind(x, y)
		return(rbind(sol, sol[,2:1]))
	}
	xy = R[2] / as;
	sol1 = solve.S2(S12, xy);
	sol2 = solve.S2(S21, xy);
	sol = cbind(sol1, sol2)
	# (x1, y2, x2, y1)
	return(sol[,c(1,4,3,2)])
}

### Examples
R = c(-1, 2)
a = c(3, 4)
sol = solve.simple.S4P3(R, a);
x1 = sol[,1]; x2 = sol[,3]; y1 = sol[,2]; y2 = sol[,4];

### Test
x1^3 + y1^3 # - R[1]
x2^3 + y2^3 # - R[1]
a[1]*x1*y2 + a[2]*x2*y1 # - R[2]
a[1]*x2*y1 + a[2]*x1*y2 # - R[2]


######################

### P3

### x^3 + y^3 = R1
### x1*y1 + a1*x1*y2 + a2*y1*x2 = R2

# x1^3 + y1^3 = R1
# x2^3 + y2^3 = R1
# x1*y1 + a1*x1*y2 + a2*x2*y1 = R2
# x2*y2 + a1*x2*y1 + a2*x1*y2 = R2

### Solution:

### Solve Pseudo-Liniar:
# (a1^2 - a2^2) * x1*y2 = (a1 - a2)*R2 - (a1*x1*y1 - a2*x2*y2)
# (a1^2 - a2^2) * x2*y1 = (a1 - a2)*R2 + (a2*x1*y1 - a1*x2*y2)

### Mult Eq 1r * 2r
x1^3*y2^3 - y1^3*x2^3 + R1*(x2^3 + y1^3) - R1^2 # = 0
# =>
# x1^3 + y2^3 # =
2*R1 - (x2^3 + y1^3)
# =>
# x2^3 = ... - y1^3
# y1^3 = R1 - x1^3 *OR*
# y2^3 = R1 - x2^3
# => Mult =>
x2^3*y1^3 + x1^3*y1^3 - R1*(...) + R1*y1^3 + (...)*x1^3 # = 0

### Diff =>
# x1^3 - y2^3 = x2^3 - y1^3

### Mult =>
x1^3*x2^3 + y1^3*y2^3 + x1^3*y2^3 + x2^3*y1^3 - R1^2 # = 0
# also:
x1^3*x2^3 - y1^3*y2^3 + R1*(y1^3 + y2^3) - R1^2 # = 0

### TODO:
# - general case;

### Special Cases:
# x1 = x2; y1 != y2;
# valid if: (a2 - a1 + 1) == 0;
(a2 - a1 + 1)*y1 - (a2 - a1 + 1)*y2 # = 0
# y1 = y2; x1 != x2;
# valid if: (a1 - a2 + 1) == 0;
(a1 - a2 + 1)*x1 - (a1 - a2 + 1)*x2 # = 0


### Solver
solve.complete.S4P3 = function(R, a) {
	a.diff = a[1] - a[2];
	if(abs(a.diff) == 1) {
		# simple solution
		m.all = unity(3, all=TRUE)
		m = m.all[2]
		if(a.diff > 0) {
			# (a1*m^j + (a2+1))*y1 = R2 / x;
			id = 1; # TODO
			div = (a[1]*m^id + (a[2]+1));
			x3 = roots(c(1, -R[1], R[2]^3 / div^3));
			x = as.vector(sapply(rootn(x3, 3), function(x) x*m.all));
			y = R[2] / x / div;
			sol = cbind(x1=x, y1=y, x2=x, y2=y*m^id);
		} else {
			# (a2*m^j + (a1+1))*x1 = R2 / y;
			id = 1; # TODO
			div = (a[2]*m^id + (a[1]+1));
			y3 = roots(c(1, -R[1], R[2]^3 / div^3));
			y = as.vector(sapply(rootn(y3, 3), function(y) y*m.all));
			x = R[2] / y / div;
			sol = cbind(x1=x, y1=y, x2=x*m^id, y2=y);
		}
	} else if(round0(a[1] + a[2]) == 0) {
		# x1*y1 + x2*y2 = 2*R2;
	}
}

### Examples:

### Special case: x1 == x2
R = c(-1, 2)
a = c(4, 3)
sol = solve.complete.S4P3(R, a)
x1 = sol[,1]; x2 = sol[,3]; y1 = sol[,2]; y2 = sol[,4];


### Special case: y1 == y2
R = c(-1, 2)
a = c(3, 4)
sol = solve.complete.S4P3(R, a)
x1 = sol[,1]; x2 = sol[,3]; y1 = sol[,2]; y2 = sol[,4];

### Test
x1^3 + y1^3 # - R[1]
x2^3 + y2^3 # - R[1]
x1*y1 + a[1]*x1*y2 + a[2]*x2*y1 # - R[2]
x2*y2 + a[1]*x2*y1 + a[2]*x1*y2 # - R[2]


### Debug
R = c(-1, 2)
a = c(3, 5)
x1 = -1.3908831512 - 2.5503258859i;
y1 =  1.4119576563 + 2.5831136271i;
x2 =  0.4925974849 + 0.7567717580i;
y2 =  0.2526662764 + 0.6178135417i;


R = c(-1, 2)
a = c(3, -3)
x1 = -0.1559428928 - 2.0642921485i;
y1 =  0.2316237980 + 2.0785818464i;
x2 = -0.7699227697 + 0.3433777972i;
y2 =  0.6635372661 - 0.7460946871i;


######################
######################

### 3D
# - possible extension to 3D sphere;

# x1^2 + y1^2 + z1^2 = R1
# x2^2 + y2^2 + z2^2 = R1
# x1*y2 + a1*y1*z2 + a2*z1*x2 = R2
# x2*y1 + a1*y2*z1 + a2*z2*x1 = R2
# x1*z2 + a3*y1*x2 + a4*z1*y2 = R3
# x2*z1 + a3*y2*x1 + a4*z2*y1 = R3

### Case:
# x1 = x2; y1 = y2; z1 = z2;
x^2 + y^2 + z^2 - R1 # = 0
x*y + a1*y*z + a2*x*z - R2 # = 0
x*z + a3*x*y + a4*y*z - R3 # = 0
### Case:
# x1 = -x2; y1 = -y2; z1 = -z2;
x^2 + y^2 + z^2 - R1 # = 0
x*y + a1*y*z + a2*x*z + R2 # = 0
x*z + a3*x*y + a4*y*z + R3 # = 0


### Variants Eqs [5]&[6]:
### V1:
# x1*y1*z1 = R3
# x2*y2*z2 = R3; # must be the same R3;
# - but does NOT preserve the negative solution;
### V2:
# x1*y1*z1 + a1*x1*y1*z2 + a2*y1*z1*x2 + a3*z1*x1*y2 = R3
# x2*y2*z2 + a1*x2*y2*z1 + a2*y2*z2*x1 + a3*z2*x2*y1 = R3
# - does NOT preserve the negative solution;
### V3:
# x1*y1*y2*z2 + a1*y1*z1*z2*x2 + a2*z1*x1*x2*y2 = R3
# x2*y2*y1*z1 + a1*y2*z2*z1*x1 + a2*z2*x2*x1*y1 = R3


#############################
#############################

#################
### Rotations ###
#################

### Order: (2,1,1)

# x1^2*x2*x3 + x2^2*x3*x4 + x3^2*x4*x1 + x4^2*x1*x2 = R1
# x1*x2^2*x3 + x2*x3^2*x4 + x3*x4^2*x1 + x4*x1^2*x2 = R2
# x1*x2*x3^2 + x2*x3*x4^2 + x3*x4*x1^2 + x4*x1*x2^2 = R3
# E3 = R4

### Solution:

### Root Structure:
# - if (x1,x2,x3,x4) is a solution,
#   then (x2,x3,x4,x1),...,(x4,x1,x2,x3) are also solutions;
# - system is decomposable: P[???] = P[4] o P[???];

### Variants:
### Variant: E4 = R4;
# - in addition, if (x1,x2,x3,x4) is a solution,
#   then {-1, i, -i} * (x1,x2,x3,x4) are also solutions;
# - system is decomposable: P[???] = P[4] o (degenerate P[4]) o P[???];
### Variant: E2 = R4;
# - in addition, if (x1,x2,x3,x4) is a solution,
#   then -1 * (x1,x2,x3,x4) is also a solution;
# - system is decomposable: P[???] = P[4] o (degenerate P[2]) o P[???];


### Sum =>
E3*S - 4*E4 - (R1+R2+R3) # = 0

### Sum(Eq[i]*Eq[i+1])
# Eq for: (x1*x2 + x2*x3 + x3*x4 + x4*x1)
4*E4^2 + E4*(E2*S^2 - E3*S - 2*E2^2 + 4*E4) - E3*E4*S + E2*E3^2 - 2*E2^2*E4 + 4*E4^2 +
	+ 2*E4*(R1+R3) + E4*((x1*x2)^2 + (x2*x3)^2 + (x3*x4)^2 + (x4*x1)^2) - (R1*R2 + R1*R3 + R2*R3) # = 0
# ((x1*x2)^2 + (x2*x3)^2 + (x3*x4)^2 + (x4*x1)^2) = (x1*x2 + x2*x3 + x3*x4 + x4*x1)^2 - 2*R2 - 4*E4
4*E4^2 + E4*(E2*S^2 - E3*S - 2*E2^2 + 4*E4) - E3*E4*S + E2*E3^2 - 2*E2^2*E4 + 4*E4^2 +
	+ 2*E4*(R1+R3) + E4*((x1*x2 + x2*x3 + x3*x4 + x4*x1)^2 - 2*R2 - 4*E4) - (R1*R2 + R1*R3 + R2*R3) # = 0
E2*E4*S^2 - 2*E3*E4*S + E2*E3^2 - 4*E2^2*E4 + 8*E4^2 +
	+ 2*E4*(R1+R3 - R2) + E4*(x1*x2 + x2*x3 + x3*x4 + x4*x1)^2 - (R1*R2 + R1*R3 + R2*R3) # = 0

### Prod =>
# Eq for: (x1*x3 + x2*x4)
E4^2*sum(x^4) + 4*E4^2*((x1*x3)^2+(x2*x4)^2) + 6*E4^2*R2 +
	+ 2*E4*(R1^2+R3^2 - 4*E4*R2 - 4*E4*((x1*x3)^2+(x2*x4)^2)) +
	+ E4*(R2*((x1*x2)^2 + (x2*x3)^2 + (x3*x4)^2 + (x4*x1)^2) - 2*R2*E4) +
	+ (4*E3*E4^2*S + E3^4 - 4*E2*E3^2*E4 + 2*E2^2*E4^2 - 4*E4^3) - R1*R2*R3 # = 0
E4^2*(S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 - 4*E4) + 4*E4^2*((x1*x3)^2+(x2*x4)^2) + 6*E4^2*R2 +
	+ 2*E4*(R1^2+R3^2 - 4*E4*R2 - 4*E4*((x1*x3 + x2*x4)^2 - 2*E4)) +
	+ E4*(R2*((x1*x2 + x2*x3 + x3*x4 + x4*x1)^2 - 2*R2 - 4*E4) - 2*R2*E4) +
	+ (4*E3*E4^2*S + E3^4 - 4*E2*E3^2*E4 + 2*E2^2*E4^2 - 4*E4^3) - R1*R2*R3 # = 0
4*E2^2*E4^2 + 8*E3*E4^2*S + E4^2*S^4 - 4*E2*E4^2*S^2 + E3^4 +
	- 8*R2*E4^2 - 4*E2*E3^2*E4 + 2*(R1^2 + R3^2 - R2^2)*E4 +
	+ R2*E4*(x1*x2 + x2*x3 + x3*x4 + x4*x1)^2 +
	- 4*E4^2*(x1*x3 + x2*x4)^2 - R1*R2*R3 # = 0
4*E2^2*E4^2 + 8*E3*E4^2*S + E4^2*S^4 - 4*E2*E4^2*S^2 - 8*R2*E4^2 +
	- 4*E2*E3^2*E4 + 2*(R1^2 + R3^2 - R2^2)*E4 + E3^4 +
	- R2*(E2*E4*S^2 - 2*E3*E4*S + E2*E3^2 - 4*E2^2*E4 + 8*E4^2 + 2*E4*(R1+R3 - R2) - (R1*R2 + R1*R3 + R2*R3)) +
	- 4*E4^2*(x1*x3 + x2*x4)^2 - R1*R2*R3 # = 0
4*E2^2*E4^2 + 8*E3*E4^2*S + E4^2*S^4 - 4*E2*E4^2*S^2 - 16*R2*E4^2 +
	- 4*E2*E3^2*E4 + 4*R2*E2^2*E4 + 2*R2*E3*E4*S - R2*E2*E4*S^2 + 2*(R1^2 + R3^2 - R2*(R1+R3))*E4 +
	+ E3^4 - R2*E2*E3^2 +
	- 4*E4^2*(x1*x3 + x2*x4)^2 +
	+ R2*(R1*R2 + R2*R3) # = 0

# E4*(x1*x2 + x2*x3 + x3*x4 + x4*x1)^2 =
- E2*E4*S^2 + 2*E3*E4*S - E2*E3^2 + 4*E2^2*E4 - 8*E4^2 - 2*E4*(R1+R3 - R2) + (R1*R2 + R1*R3 + R2*R3)
# 4*E4^2*(x1*x3 + x2*x4)^2 =
4*E2^2*E4^2 + 8*E3*E4^2*S + E4^2*S^4 - 4*E2*E4^2*S^2 - 16*R2*E4^2 +
	- 4*E2*E3^2*E4 + 4*R2*E2^2*E4 + 2*R2*E3*E4*S - R2*E2*E4*S^2 + 2*(R1^2 + R3^2 - R2*(R1+R3))*E4 +
	+ E3^4 - R2*E2*E3^2 + R2*(R1*R2 + R2*R3) # = 0
# 2*E2*E4 = sqrt(4*E4*Eq2a) + sqrt(Eq2b);
# - unknown if signs of radicals have any effect?

