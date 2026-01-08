########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S2
### Heterogeneous Symmetric
###
### draft v.0.4d


### Heterogeneous Symmetric Polynomial Systems

### 2 Variables:
### Simple System:
# x^n + P(x, y) = R
# y^n + P(y, x) = R

### Generalised System:
# P(x, y) = 0
# P(y, x) = 0
# where P(x, y) is any polynomial of 2 variables;
# e.g. P(x, y) = x^n + P2(x, y) - R;


###############
### Systems ###

### 2 Variables:
# x^n + P(x, y) = R
# y^n + P(y, x) = R

### Simple Leading Term
# P3.1) x^3 + b*y = R; [P3 => P6]
# P3.2) (x - s)^3 + b*y = R; [P3 => P6: equivalent to non-shifted]
# P3.3) x^3 + b*x*y = R; [P3 => trivial P6]
# P3.4) (x - s)^3 + b*x*y = R; (P3 => P6)
# P3.5) x^3 + b2*x*y + b1*x = R; (P3 => P6)
# P3.6) x^3 + b2*x*y + b1*y = R; (P3 => P6)
# P3.7) x^3 + b3*x*y + b2*y^2 + b1*y = R; (P3 => P6)
# P3.8) x^3 + b*(x*y)^2 = R; [P4 => P8: based on m3]
# P3.9) x^3 + b2*(x*y)^2 + b1*(x*y) = R; [P4 => P8: based on m3]
# P3.10.) x^3 + b3*x^2*y + b2*x*y^2 + b1*x = R; (P3 => P6)
### Order 4 & 5:
# P4.1a.) x^4 + b*y = R; (P6 => P12)
# P4.1b.) x^4 + b2*x*y + b1*y = R; (P6 => P12)
# P4.1c.) x^4 + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P6 => P12; if(b3 == 1) P4 => P8)
# P4.1d.) (x - s)^4 + b*y = R; (P6 => P12; simple shift of P4.1a)
# P4.1e.) x^4 + b3*y^3 + b2*y^2 + b1*y = R; (P6 => P12)
# P4.2a.) x^4 + b*x*y = R; (trivial P8)
# P4.2b.) x^4 + b2*(x*y)^2 + b1*x*y = R; (TODO: trivial P8)
# P5.1.) x^5 + b*(x+y) = 0; (TODO: based on (P4[x^4])^2)
# P5.2.) x^5 + b2*x*y + b1*(x+y) = 0; (TODO: based on (P16)^2)
# P5.3.) x^5 + b*y = R; (P10 => P20)
# P5.4.) x^5 + b*y^4 = R; (P10 => P20)
# P5.5.) x^5 + b2*y^2 + b1*y = R; (also P10)
# P5.6.) x^5 + b4*x^4 + b3*y^3 + b2*y^2 + b1*y = R; (also P10)
# Shifted Side-Chain:
# P5.5.) x^5 + b3*x^3*y^2 + b2*x^2*y = R; (P10 => P20)
# P5.6.) x^5 + b3*x^3*y^2 + b2*x^2*y + b1*x = R; (also P10)
# P5.7.) x^5 + b3*x^3*y^2 + b2*x^2*y + b1*y = R; (also P10)


### Complex Leading Term/Terms:
# - moved to other files;

### 2 Leading Terms:
# - moved to file:
#   Poly.System.Hetero.Symmetric.S2.L2.R;
# B3.1) a1*x^3 + a2*y^3 + b*x = R; (P3 => P6)
# B3.2) a1*x^3 + a2*y^3 + b*x*y = R; (P3 => P6)
# B3.3) a1*x^3 + a2*y^3 + b2*x*y + b1*x = R; (P3 => P6)

### Mixed-Leading-Term:
# - moved to file:
#   Poly.System.Hetero.Symmetric.S2.Leading.R;
# - see complete summary in the specific file;
### Mixed: Order 2+1
# M21.1) x^2*y + b*x: NO solutions (x != y);
# M21.2) x^2*y + b*y: trivial;
# M21.3) x^2*y + b2*x^2 + b1*x: trivial (trivial P2);
# M21.4) x^2*y + b3*x*y + b2*x^2 + b1*x: trivial (trivial P2);
### Mixed: Order n+1
# M31.1) x^3*y + b*x: trivial P4;
# M31.2) x^3*y + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P3 => P6; some nice)
# M32.1) TODO: x^3*y^2 + b2*x^2 + b1*x = R;
# M41.1) x^4*y + b*x; (P5 => P10)
### Mixed: Order n+3
# M43.1) x^4*y^3 + b3*x*y + b2*x^2 + b1*x = R; (trivial P2; base P7)
# M43.2) x^4*y^3 + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P3 => P6)
# M43.3) x^4*y^3 + b5*x^2*y + b4*x*y^2 + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P3 => P6)
# ...

### Mixed-Trivial: (x*y)^n
# - moved to file:
#   Poly.System.Hetero.Symmetric.S2.Lnn.R;
# Mt2.1.) x^2*y^2 + b2*x^2 + b1*x: simple (P2 => P4);
# Mt2.2.) x^2*y^2 + b3*x^2*y + b2*x^2 + b1*x; (TODO: P2 => P4)
# Mt2.3.) TODO: x^2*y^2 + b4*x^2*y + b3*x*y^2 + b2*x^2 + b1*x
### Mixed: (x*y)^n
# Mt3.1.) TODO: (x*y)^3 + ...:
#      (x*y)^3 + b*x^3 = R;
#      (x*y)^3 + b*x^4 + b2*(x*y)^2 + b1*x*y = R;
# Mt5.1) TODO: (x*y)^5 + b*x^3 = R

### Mixed: other
# - moved to file:
#   Poly.System.Hetero.Symmetric.S2.LMixed.R;
# Mm3.1.) x^3 + a*x^3*y = R; (TODO: P3 => P6)
# Mm3.2.) x^3 + a*x^3*y + b1*x*y = R;
# Mm3.3.) x^3 + a*x*y^3 = R; (TODO: P5 => P10)
#         https://www.youtube.com/watch?v=jRwiJTaAu4c
### Mixed 2 High-Power Terms:
# MT.1) a1*x^3*y + a2*x*y^3 + b*x = R; [P3 => interesting P6]


###############

###############
### History ###
###############

### [branch v.0.3]
#
### v.0.3e & v.0.3g: [cleanup]
# - moved Old History to file:
#   Poly.System.Hetero.Symmetric.S2.History.R;
# - moved section with Multiple Simple Leading terms to file:
#   Poly.System.Hetero.Symmetric.S2.L2.R;
# - moved Multiple/Mixed Leading to file:
#   Poly.System.Hetero.Symmetric.S2.LMixed.R;
# - moved Lnn (Mixed/Trivial Leading) to file:
#   Poly.System.Hetero.Symmetric.S2.Lnn.R;
# - moved section with Mixed Leading Term to new file:
#   Poly.System.Hetero.Symmetric.S2.Leading.R;
# - more refactoring & cleanup;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")


# library(polynom)
# library(pracma)

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;


### Other

polyGen.Simple = function(n, m, div=TRUE) {
	p = toPoly.pm("(R - x^n)^n - b^n*(R - b*x^m)^m");
	p$coeff = - p$coeff;
	if(div) p = div.pm(p, toPoly.pm("x^n + b*x^m - R"))$Rez;
	p = orderVars.pm(p, c("b", "R", "x"));
	p = sort.pm(p, "x");
	return(p);
}

##########################
##########################

##########################
### Polynomial Systems ###
##########################

###############
### Order 3 ###
###############

### x^3 + b*y

# x^3 + b1*y = R
# y^3 + b1*x = R

### Extensions:
# E1: x^3 + b2*x*y + b1*y = R; [P3 => P6]
# E2: x^3 + b3*(x*y)^2 + b2*x*y + b1*y = R; [P4 => P8]
# - are discussed in a separate section;

### Solution:

### Diff =>
# x*y = S^2 - b1;

### Sum =>
S^3 - 3*E2*S  + b1*S - 2*R = 0

### Eq:
S^3 - 2*b1*S + R # = 0

### Step 2:
# Solve:
# x + y = S
# x*y = S^2 - b1


### Solver:
solve.htShift = function(b, R, shift=0) {
	s = shift;
	r.sum = roots(c(1, - 6*s, - 2*(b[1]-6*s^2), - 8*s^3 + R + 3*s*b[1]))
	xy = r.sum^2 - 3*s*r.sum + 3*s^2 - b[1];
	r.diff = sqrt(r.sum^2 - 4*xy + 0i)
	x = (r.sum + r.diff)/2
	y = (r.sum - r.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	sol # TODO: include also x = y cases
}
solve.ht.S3P3 = function(R, b) {
	# using the new model of Extensions;
	b2 = if(length(b) > 1) b[2] else 0;
	b3 = if(length(b) > 2) b[3] else 0;
	coeff = c(1, 0, -2*b[1], R[1])
	if(b2 != 0 || b3 != 0) coeff = coeff + c(0, -b3, -b2, 0)
	S = roots(coeff)
	print(S) # Debug
	R1 = R[1] - b2*S - b3*S^2;
	xy = S^2 - b[1];
	len = length(S)
	# robust
	if(R[1] == 0) {
		# TODO
		r.diff = sqrt(S^2 - 4*xy + 0i);
		x = (S + r.diff)/2
	} else {
		x = sapply(seq(len), function(id)
			roots(c(R1[id], (xy[id]^2 - b[1]^2), R1[id]*(2*xy[id] - S[id]^2 + b[1]))))
		S = matrix(S, ncol=len, nrow=2, byrow=T)
	}
	y = S - x
	sol = cbind(x=as.vector(x), y=as.vector(y))
	sol # TODO: include also x = y cases
}

### Examples:

b = 3
R = 1
#
sol = solve.htShift(b, R)
x = sol[,1]; y = sol[,2];
sol # TODO: include also x = y cases

### Test
x^3 + b[1]*y
y^3 + b[1]*x

### Classic Polynomial
err = x^6 - b[1]*x^4 - 2*R*x^3 + b[1]^2*x^2 + b[1]*R*x + R^2 - b[1]^3
round0(err)


### Extensions: A1-type

#########
### Ex 2:
R = 1
b = c(1, 1)
#
sol = solve.ht.S3P3(R, b)
x = sol[,1]; y = sol[,2];
sol

### Test
x^3 + b[1]*y + b[2]*(x + y)
y^3 + b[1]*x + b[2]*(x + y)

round0.p(poly.calc(x))
err = -3 + 3*x^2 - 2*x^3 + x^6
round0(err)


#########
### Ex 3:
R = 1
b = c(1, 0, 1)
#
sol = solve.ht.S3P3(R, b)
x = sol[,1]; y = sol[,2];
sol

### Test
x^3 + b[1]*y + b[2]*(x + y) + b[3]*(x + y)^2
y^3 + b[1]*x + b[2]*(x + y) + b[3]*(x + y)^2

round0.p(poly.calc(x))
err = -1 - 2*x + 2*x^3 - x^5 + x^6
round0(err)


#####################
### Shifted Roots ###
#####################

### (x - s)^3 + b*y

# (x - s)^3 + b1*y = R
# (y - s)^3 + b1*x = R

### Equivalent:
# - trivial shift:
#  (x - s)^3 + b1*(y - s) = R - s*b1;
#  (y - s)^3 + b1*(x - s) = R - s*b1;
#   solve: x => x - s; y => y - s;
# - after shift-back: only a shift in R;

solve.htShift = function(b, R, shift=0) {
	s = shift;
	r.sum = roots(c(1, - 6*s, - 2*(b[1]-6*s^2), - 8*s^3 + R + 3*s*b[1]))
	xy = r.sum^2 - 3*s*r.sum + 3*s^2 - b[1];
	r.diff = sqrt(r.sum^2 - 4*xy + 0i)
	x = (r.sum + r.diff)/2
	y = (r.sum - r.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	sol # TODO: include also x = y cases
}

### Solution:

### Example
b = 2
R = 1
s = 1
#
sol = solve.htShift(b, R, shift=s)
x = sol[,1]; y = sol[,2];
sol

### Test
(x-s)^3 + b[1]*y
(y-s)^3 + b[1]*x

### Classic Polynomial
poly.calc(sol[,1])
round0.p(poly.calc(sol[,1] - s))

# back-shift
x = x - s
err = x^6 - b[1]*x^4 - 2*(R - b[1]*s)*x^3 + b[1]^2*x^2 + b[1]*(R - b[1]*s)*x + (R - b[1]*s)^2 - b[1]^3
round0(err)


### Example 2:
b = 2; R = 1;
s = 3/2
#
sol = solve.htShift(b, R, shift=s)
x = sol[,1]; y = sol[,2];
x = sol[,1] - s; # with shift back!
-4 - 4*x + 4*x^2 + 4*x^3 - 2*x^4 + x^6

### Example 3:
b = 2; R = 1;
s = 5/2
#
sol = solve.htShift(b, R, shift=s)
x = sol[,1]; y = sol[,2];
x = sol[,1] - s; # with shift back!
8 - 8*x + 4*x^2 + 8*x^3 - 2*x^4 + x^6

###
p = sapply((-6:6) + 1/2, function(s) print(round0.p(poly.calc(solve.htShift(2, 1, shift=s)[,1] - s))))
# [some b0] - 4*s0*x + 4*x^2 + 4*s0*x^3 - 2*x^4 + x^6 # s0 = s - 1/2

### Classic Polynomial

# - some experimentation;

### some P12s
shiftSqrt.p = function(b, R, shift) {
	s = sqrt(shift + 0i)
	r = c(solve.htShift(b, R, shift=s)[,1] - s, solve.htShift(b, R, shift=-s)[,1] + s)
	list(r=r, p=round0.p(poly.calc(r)))
}

b = 2; R = 1;
p = sapply(-6:6, function(s) print(shiftSqrt.p(b, R, shift=s/2)$p))
#
409 + 20*x - 100*x^2 - 4*x^3 - 12*x^4 - 24*x^5 - 2*x^6 + 12*x^7 + 12*x^8 - 4*x^9 - 4*x^10 + x^12 
329 + 12*x - 92*x^2 + 4*x^3 - 4*x^4 - 24*x^5 - 6*x^6 + 12*x^7 + 12*x^8 - 4*x^9 - 4*x^10 + x^12 
257 + 4*x - 84*x^2 + 12*x^3 + 4*x^4 - 24*x^5 - 10*x^6 + 12*x^7 + 12*x^8 - 4*x^9 - 4*x^10 + x^12 
193 - 4*x - 76*x^2 + 20*x^3 + 12*x^4 - 24*x^5 - 14*x^6 + 12*x^7 + 12*x^8 - 4*x^9 - 4*x^10 + x^12 
137 - 12*x - 68*x^2 + 28*x^3 + 20*x^4 - 24*x^5 - 18*x^6 + 12*x^7 + 12*x^8 - 4*x^9 - 4*x^10 + x^12 
89 - 20*x - 60*x^2 + 36*x^3 + 28*x^4 - 24*x^5 - 22*x^6 + 12*x^7 + 12*x^8 - 4*x^9 - 4*x^10 + x^12 
49 - 28*x - 52*x^2 + 44*x^3 + 36*x^4 - 24*x^5 - 26*x^6 + 12*x^7 + 12*x^8 - 4*x^9 - 4*x^10 + x^12 
17 - 36*x - 44*x^2 + 52*x^3 + 44*x^4 - 24*x^5 - 30*x^6 + 12*x^7 + 12*x^8 - 4*x^9 - 4*x^10 + x^12 
-7 - 44*x - 36*x^2 + 60*x^3 + 52*x^4 - 24*x^5 - 34*x^6 + 12*x^7 + 12*x^8 - 4*x^9 - 4*x^10 + x^12 
-23 - 52*x - 28*x^2 + 68*x^3 + 60*x^4 - 24*x^5 - 38*x^6 + 12*x^7 + 12*x^8 - 4*x^9 - 4*x^10 + x^12 
-31 - 60*x - 20*x^2 + 76*x^3 + 68*x^4 - 24*x^5 - 42*x^6 + 12*x^7 + 12*x^8 - 4*x^9 - 4*x^10 + x^12 
-31 - 68*x - 12*x^2 + 84*x^3 + 76*x^4 - 24*x^5 - 46*x^6 + 12*x^7 + 12*x^8 - 4*x^9 - 4*x^10 + x^12 
-23 - 76*x - 4*x^2 + 92*x^3 + 84*x^4 - 24*x^5 - 50*x^6 + 12*x^7 + 12*x^8 - 4*x^9 - 4*x^10 + x^12


###################

###################
### x^3 + b*x*y ###

# x^3 + b1*x*y = R
# y^3 + b1*x*y = R

# relatively trivial
# (x^3 - b[1]/2*x^2 - R)^2 + 3/4 * b[1]^2*x^4 = 0

### Solution:

### Diff =>
# x*y = (x + y)^2

### Sum =>
Z^3 - b1*Z^2 + R # = 0

solve.htxy = function(b, R) {
	x.sum = roots(c(1, - b[1], 0, R))
	xy = x.sum^2
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	sol
}

### Example:

b = 3
R = 1
#
sol = solve.htxy(b, R)
x = sol[,1]; y = sol[,2];
sol

### Test
x^3 + b[1]*x*y
y^3 + b[1]*x*y

### Classical Polynomial
round0.p(poly.calc(sol[,1]))

err = x^6 - b[1]*x^5 + b[1]^2*x^4 - 2*R*x^3 + b[1]*R*x^2 + R^2
round0(err)


### Example 2:
b = 5
#
sol = solve.htxy(b, 1) # R = 1;
x = sol[,1]; y = sol[,2];
sol
err = 1 + b[1]*x^2 - 2*x^3 + b[1]^2*x^4 - b[1]*x^5 + x^6
round0(err)


### Classical Polynomial
x^6 - b[1]*x^5 + b[1]^2*x^4 - 2*R*x^3 + b[1]*R*x^2 + R^2
(x^3 - b[1]/2*x^2 - R)^2 + 3/4 * b[1]^2*x^4 # relatively trivial


#################

#################
### Shifted Roots

### (x - s)^3 + b*x*y

# (x-s)^3 + b1*x*y = R
# (y-s)^3 + b1*x*y = R

### Solution:

### Diff =>
# x*y = Z^2 - 3*s*Z + 3*s^2

### Sum =>
Z^3 - (6*s+b1)*Z^2 + (12*s^2+3*b1*s)*Z - 8*s^3 - 3*b1*s^2 + R


solve.htxyShift = function(b, R, s) {
	r.sum = roots(c(1, - (6*s+b[1]), (12*s^2+3*b[1]*s), - 8*s^3 - 3*b[1]*s^2 + R))
	xy = r.sum^2 - 3*s*r.sum + 3*s^2
	r.diff = sqrt(r.sum^2 - 4*xy + 0i)
	x = (r.sum + r.diff)/2
	y = (r.sum - r.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	sol
}
poly.htxy = function(b, R, s) {
	sol = solve.htxyShift(b, R, s=s)
	p1 = round0.p(poly.calc(sol[,1]))
	p = round0.p(poly.calc(sol[,1] - s))
	p.coeff = c((b[1]*s^2-R)^2, b[1]*s*(b[1]*s^2 - R), b[1]*R, (2*b[1]*s^2 + b[1]^2*s - 2*R), b[1]*(s+b[1]), - b[1], 1)
	return(list(r=sol, coeff=p.coeff, p1=p1, p=p))
}

### Example
b = 2
R = 1
s = 3
#
sol = solve.htxyShift(b, R, s=s)
x = sol[,1]; y = sol[,2];
sol

### Test
(x-s)^3 + b[1]*x*y
(y-s)^3 + b[1]*x*y

### Classical Polynomial

# back-shifted:
x = x - s
err = (b[1]*s^2-R)^2 + b[1]*s*(b[1]*s^2 - R)*x + b[1]*R*x^2 + (2*b[1]*s^2 + b[1]^2*s - 2*R)*x^3 + b[1]*(s+b[1])*x^4 - b[1]*x^5 + x^6
round0(err)

p.coeff = c((b[1]*s^2-R)^2, b[1]*s*(b[1]*s^2 - R), b[1]*R, (2*b[1]*s^2 + b[1]^2*s - 2*R), b[1]*(s+b[1]), - b[1], 1)
p.coeff

round0.p(poly.calc(sol[,1]))
# back-shift
round0.p(poly.calc(sol[,1] - s))


### Example 2:
b = -1; s = 1; R = -2;
sol = poly.htxy(b, R, s=s)
x = sol$r[,1]; y = sol$r[,2];
sol
x = sol$r[,1] - s; # shift back;
err = 1 - x + 2*x^2 + 3*x^3 + x^5 + x^6
round0(err)


# back-shifted:
(x^3 + b[1]*x^2*m^1 - b[1]*s*x*m^2 + b[1]*s^2 - R)*(x^3 + b[1]*x^2*m^2 - b[1]*s*x*m + b[1]*s^2 - R)


###
b = 2; R = 1
p = sapply(-6:6, function(s) print(poly.htxy(b, R, s)$p))
# back-shifted:
# (b*s^2-R)^2 + b*s*(b*s^2 - R)*x + b*R*x^2 + (2*b*s^2 + b^2*s - 2*R)*x^3 - b*(s+b)*x^4 - b*x^5 + x^6
5041 - 852*x + 2*x^2 + 118*x^3 - 8*x^4 - 2*x^5 + x^6
2401 - 490*x + 2*x^2 + 78*x^3 - 6*x^4 - 2*x^5 + x^6
961 - 248*x + 2*x^2 + 46*x^3 - 4*x^4 - 2*x^5 + x^6
289 - 102*x + 2*x^2 + 22*x^3 - 2*x^4 - 2*x^5 + x^6
49 - 28*x + 2*x^2 + 6*x^3 - 0 - 2*x^5 + x^6
1 - 2*x + 2*x^2 - 2*x^3 + 2*x^4 - 2*x^5 + x^6
1 - 0 + 2*x^2 - 2*x^3 + 4*x^4 - 2*x^5 + x^6
1 + 2*x + 2*x^2 + 6*x^3 + 6*x^4 - 2*x^5 + x^6 
49 + 28*x + 2*x^2 + 22*x^3 + 8*x^4 - 2*x^5 + x^6 
289 + 102*x + 2*x^2 + 46*x^3 + 10*x^4 - 2*x^5 + x^6 
961 + 248*x + 2*x^2 + 78*x^3 + 12*x^4 - 2*x^5 + x^6 
2401 + 490*x + 2*x^2 + 118*x^3 + 14*x^4 - 2*x^5 + x^6 
5041 + 852*x + 2*x^2 + 166*x^3 + 16*x^4 - 2*x^5 + x^6


######################
######################

######################
### xy & x/y-Terms ###
######################

### x-Term
### x^3 + b2*x*y + b1*x
# x-Terms vs y-Terms: are equivalent;

# x^3 + b2*x*y + b1*x = R
# y^3 + b2*x*y + b1*y = R

### Solution:

### Diff =>
# x*y = S^2 + b1

### Sum =>
S^3 - b2*S^2 + b1*S - b1*b2 + R # = 0


solve.htxy = function(b, R, isX=TRUE) {
	if(length(b) < 2) b = c(0, b) # TODO: verify order
	if(isX) {
		x.sum = roots(c(1, - b[2], b[1], - b[1]*b[2] + R))
		xy = x.sum^2 + b[1]
	} else {
		x.sum = roots(c(1, - b[2], -2*b[1], b[1]*b[2] + R))
		xy = x.sum^2 - b[1]
	}
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	sol
}

### Example
b = c(2, 3)
R = 1
#
sol = solve.htxy(b, R)
x = sol[,1]; y = sol[,2]
sol

### Test
x^3 + b[2]*x*y + b[1]*x 
y^3 + b[2]*x*y + b[1]*y

### Classical Polynomial
round0.p(poly.calc(sol[,1]))

err = x^6 - b[2]*x^5 + (b[2]^2 + 2*b[1])*x^4 - (b[1]*b[2] + 2*R)*x^3 +
	+ (b[1]*b[2]^2 + R*b[2] + b[1]^2)*x^2 - 2*b[1]*R*x + R^2
round0(err)


### Example 2: 4 + 8*x - 2*x^5 + x^6
b = c(-2, 2)
R = 2
#
sol = solve.htxy(b, R)
x = sol[,1]; y = sol[,2]
sol

### Test
x^3 + b[2]*x*y + b[1]*x 
y^3 + b[2]*x*y + b[1]*y

### Classical Polynomial
round0.p(poly.calc(sol[,1]))

err = x^6 - b[2]*x^5 + (b[2]^2 + 2*b[1])*x^4 - (b[1]*b[2] + 2*R)*x^3 +
	(b[1]*b[2]^2 + R*b[2] + b[1]^2)*x^2 - 2*b[1]*R*x + R^2
round0(err)


###############

###############
### y-Term ####

### x^3 + b2*x*y + b1*y

# x^3 + b2*x*y + b1*y = R
# y^3 + b2*x*y + b1*x = R

### Solution:

### Diff =>
# x*y = (x + y)^2 - b1

### Eq:
S^3 - b2*S^2 - 2*b1*S + b1*b2 + R


### Example
b = c(-1, 1)
R = 2
#
sol = solve.htxy(b, R, isX=FALSE)
x = sol[,1]; y = sol[,2]
sol

### Test
x^3 + b[2]*x*y + b[1]*y
y^3 + b[2]*x*y + b[1]*x

### Classical Polynomial
err = (R^2 - b[1]^3) + (R*b[1] - 2*b[1]^2*b[2])*x + (R*b[2] - b[1]*b[2]^2 + b[1]^2)*x^2 +
	- 2*(R - b[1]*b[2])*x^3 - (b[1] - b[2]^2)*x^4 - b[2]*x^5 + x^6
round0(err)


round0.p(poly.calc(sol[,1]))


### Example 2
b = c(-1, -1)
R = 2
#
sol = solve.htxy(b, R, isX=FALSE)
x = sol[,1]; y = sol[,2]
sol

### Test
x^3 + b[2]*x*y + b[1]*y
y^3 + b[2]*x*y + b[1]*x

err = 5 - 2*x^3 + 2*x^4 + x^5 + x^6
round0(err)



#########################
#########################

### Mixed Side Chain:
### x^3 + b3*x*y + b2*y^2 + b1*y

# x^3 + b3*x*y + b2*y^2 + b1*y = R
# y^3 + b3*x*y + b2*x^2 + b1*x = R

### Solution:

### Diff =>
# x*y = S^2 - b2*S - b1

### Sum =>
S^3 - (b[2] + b[3])*S^2 - (2*b[1] - b[2]*b[3] + b[2]^2)*S - b[1]*b[2] + b[1]*b[3] + R


### Solver:
solve.htxy = function(b, R) {
	x.sum = roots(c(1, - (b[2] + b[3]), - (2*b[1] - b[2]*b[3] + b[2]^2), - b[1]*b[2] + b[1]*b[3] + R))
	xy =  x.sum^2 - b[2]*x.sum - b[1]
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	#
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}

### Example
b = c(3,-1,1)
R = 1
#
sol = solve.htxy(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^3 + b[3]*x*y + b[2]*y^2 + b[1]*y
y^3 + b[3]*x*y + b[2]*x^2 + b[1]*x

### Classical Polynomial
round0.p(poly.calc(sol$sol[,1]))


### Example 2:
b = c(-1,1,-1)
p = sapply(-6:6, function(r) print(solve.htxy(b, r)$p))
# (R + 1)^2 + 3*x^2 - 2*(R + 2)*x^3 + 3*x^4 + x^6

### Examples:
b = c(-1,-2,2)
p = sapply(-6:6, function(r) print(solve.htxy(b, r)$p))
# (R + 1)^2 - 3*(R + 4)*x + 33*x^2 - 2*(R - 16)*x^3 + 9*x^4 + x^6

### 10 + 6*x + 6*x^2 + x^6
b = c(2, 1, -1)
R = 2*b[3]^3
#
sol = solve.htxy(b, R)
x = sol$sol[,1]
10 + 6*x + 6*x^2 + x^6


### Classical polynomial:
coeff = c(1, -(b[2] + b[3]), -b[1] + b[2]^2 + b[3]^2, -(3*b[2]^3 + b[2]*b[3]^2 - 2*b[1]*b[2] - 2*b[1]*b[3] + 2*R),
	(R[1]*b[2] + R[1]*b[3] - b[1]*b[2]*b[3] - b[1]*b[3]^2 + b[1]^2 - b[2]^3*b[3] + b[2]^4),
	(b[1] - b[2]*b[3])*R - 2*b[1]^2*b[3] + b[1]*b[2]^3 - b[1]*b[2]^2*b[3],
	- R*(3*b[1]*b[2] + b[2]^3 - R) - b[1]^3)
coeff



###################
###################

###################
### (xy)^2 Term ###
###################

### x^3 + b*(x*y)^2

# x^3 + b1*(x*y)^2 = R
# y^3 + b1*(x*y)^2 = R

### Solution:

m3 = unity(3, all=F)

### Diff =>
# y = x*m, where m^3 = 1;
# separate equations for: m & m^2
# b1*x^4*m^2 + x^3 - R = 0
# b1*x^4*m + x^3 - R = 0

solve.xysq = function(b, R, isInverse=FALSE) {
	# isInverse: b = 1/n and parameter b = n;
	coeff = if(isInverse) c(m3^2, b[1], 0,0, - R*b[1]) else c(b[1]*m3^2, 1,0,0, - R);
	x = roots(coeff)
	y = x*m3
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	if(isInverse) coeff = c(R^2*b[1]^2, 0,0, - 2*R*b[1]^2, b[1]*R, 0, b[1]^2, - b[1], 1)
	else coeff = c(R^2, 0,0, - 2*R, b[1]*R, 0, 1, - b[1], b[1]^2) / b[1]^2
	return(list(sol=sol, coeff=coeff))
}

### Example
b = 1/2
R = 1
#
sol = solve.xysq(b, R, isInverse=F)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^3 + b[1]*(x*y)^2
y^3 + b[1]*(x*y)^2

### Classic Polynomial
err = b[1]^2*x^8 - b[1]*x^7 + x^6 + b[1]*R*x^4 - 2*R*x^3 + R^2
round0(err)

round0.p(poly.calc(x))
err = round0(4 - 8*x^3 + 2*x^4 + 4*x^6 - 2*x^7 + x^8)
err


##################
### Extensions ###

### x^3 + b2*(x*y)^2 + b1*(x*y)

# x^3 + b2*(x*y)^2 + b1*(x*y) = R
# y^3 + b2*(x*y)^2 + b1*(x*y) = R

### Solution:

m3 = unity(3, all=F)

### Diff =>
# y = x*m, where m^3 = 1;
# separate equations for: m & m^2
# b2*m^2*x^4 + x^3 + b1*m^1*x^2 - R = 0
# b2*m^1*x^4 + x^3 + b1*m^2*x^2 - R = 0

solve.xysq = function(b, R, isInverse=FALSE) {
	# TODO: verify isInverse!
	coeff = if(isInverse) c(m3^2, b[2], b[1]*b[2]*m3,0, - R*b[2]) else c(b[2]*m3^2, 1,b[1]*m3,0, - R);
	x = roots(coeff)
	y = x*m3
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	# TODO: verify Coeff!
	# b[2]^2*x^8 - b[2]*x^7 - (b[1]*b[2] - 1)*x^6 - b[1]*x^5 + (b[1]^2 + b[2]*R)*x^4 - 2*R*x^3 + b[1]*R*x^2 + R^2
	if(isInverse) coeff = c(R^2*b[2]^2, 0, b[1]*b[2]^2*R, - 2*b[2]^2*R, b[1]^2*b[2]^2 + b[2]*R, -b[1]*b[2]^2, -(b[1]*b[2]-b[2]^2), - b[2], 1)
	else coeff = c(R^2, 0, b[1]*R, - 2*R, b[1]^2 + b[2]*R, -b[1], -(b[1]*b[2]-1), - b[2], b[2]^2) / b[2]^2
	return(list(sol=sol, coeff=coeff))
}

### Example
b = c(1, 1/2)
R = 1
#
sol = solve.xysq(b, R, isInverse=F)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^3 + b[2]*(x*y)^2 + b[1]*x*y
y^3 + b[2]*(x*y)^2 + b[1]*x*y

### Classical Polynomial
round0.p(poly.calc(sol$sol[,1]))

b[2]^2*x^8 - b[2]*x^7 - (b[1]*b[2] - 1)*x^6 - b[1]*x^5 + (b[1]^2 + b[2]*R)*x^4 - 2*R*x^3 + b[1]*R*x^2 + R^2


#############################
#############################

#######################
### Base: x^n       ###
### + x^j*y^k Terms ###
#######################

### x^3 + b3*x^2*y + b2*x*y^2 + b1*x

# x^3 + b3*x^2*y + b2*x*y^2 + b1*x = R
# y^3 + b3*y^2*x + b2*y*x^2 + b1*y = R

### Solution:

### Diff =>
# (b2 - b3 + 1)*x*y = S^2 + b1 

### Eq:
(b2 - 1)*S^3 + b1*S*(b2 - 1) - R*(b2 - b3 + 1)

### TODO: verify
# Case: b2 - b3 + 1 == 0;
# Case: b2 == 1;

solve.htxy2 = function(b, R) {
	if(round0(b[3] - b[2] - 1, tol=1E-10) == 0) {
		x.sum = sqrt(-b[1] + 0i) * c(1, -1)
		# TODO: all cases;
		if(round0(b[2]+b[3] - 3) == 0) {
			xy = 0; # b2 = 1; b3 = 2;
			print("Only trivial solution!")
		} else {
			xy = (x.sum^3 + b[1]*x.sum - 2*R) / x.sum / (b[2] + b[3] - 3)
		}
	} else {
		if(b[2] == 1) print("Only trivial solution!")
		x.sum = roots(c(1, 0, b[1], - R*(b[2] - b[3] + 1)/(b[2] - 1)))
		xy = (x.sum^2 + b[1]) / (b[2] - b[3] + 1)
	}
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	sol
}

### Example:
b = c(3, 2, 1)
R = 1
#
sol = solve.htxy2(b, R)
x = sol[,1]; y = sol[,2]
sol

### Test
x^3 + b[3]*x^2*y + b[2]*x*y^2 + b[1]*x
y^3 + b[3]*y^2*x + b[2]*y*x^2 + b[1]*y

### Classic Polynomial
round0.p(poly.calc(sol[,1]))

b1 = b[1]; b2 = b[2]; b3 = b[3]; R = R[1];
#
err = R^2 +
	2*R*b[1] * (b[2] - 1)*x +
	b[1]^2*(1 - b[2])^2 *x^2 +
	R * (- 2 + b2 - 3*b2*b3 - b2*b3^2 + 2*b2^2 + 2*b2^2*b3 - b2^3 + b3 + b3^2)*x^3 +
	b[1]*(1 - b[2])^2 * (b2 - b3 + 2)*x^4 +
	(1 - b[2])^2 * (b2 - b3 + 1)*x^6
round0(err)


#############################
#############################

###############
### Order 4 ###
###############

### x^4 + b*y

### Structural Extensions:
### Extension 1: x^4 + b2*x*y + b1*y = R
### Extension 2: x^4 + b3*(x*y)^2 + b2*x*y + b1*y = R

# x^4 + b1*y = R
# y^4 + b1*x = R

### Solution:

### Diff =>
# x*y = (S^2 - b1/S)/2

### Sum =>
S^6 - 4*b1*S^3 + 4*R*S^2 - b1^2

### Simple (x*y) extensions:
### E1: x^4 + b2*x*y + b1*y = R
S^6 - 2*b2*S^4 - 4*b1*S^3 + 4*R*S^2 + 2*b1*b2*S - b1^2 # = 0
### E2: x^4 + b3*(x*y)^2 + b2*x*y + b1*y = R
(b3 - 1)*S^6 + 2*b2*S^4 - (2*b1*b3 - 4*b1)*S^3 - 4*R*S^2 - 2*b1*b2*S + b1^2 + b1^2*b3 # = 0

### Solver:
solve.ht4 = function(R, b, b.ext=0) {
	if(length(b) == 1) coeff = c(1, 0, 0, -4*b[1], 4*R, 0, - b[1]^2)
	else if(length(b) == 2) coeff = c(1, 0, -2*b[2], -4*b[1], 4*R, 2*b[1]*b[2], - b[1]^2)
	else if(length(b) == 3) coeff = - c(b[3]-1, 0, 2*b[2], -2*b[1]*b[3]+4*b[1], -4*R,
		-2*b[1]*b[2], b[1]^2 + b[1]^2*b[3]);
	# A1-type Extensions:
	b.e1 = if(length(b.ext) >= 1) b.ext[1] else 0;
	b.e2 = if(length(b.ext) >= 2) b.ext[2] else 0;
	if(b.e1 != 0 || b.e2 != 0) coeff = coeff + c(0,0, -4*b.e2, -4*b.e1, 0,0,0)
	
	r.sum = roots(coeff);
	r.sum = r.sum[r.sum != 0] # if b[3] == -1
	xy = (r.sum^2 - b[1]/r.sum)/2
	r.diff = sqrt(r.sum^2 - 4*xy + 0i)
	x = (r.sum + r.diff)/2
	y = (r.sum - r.diff)/2
	sol = cbind(x, y) # TODO: add also x = y cases;
	sol = round0(rbind(sol, sol[,2:1]))
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}

### Example 1:
R = 1
b = 2
#
sol = solve.ht4(R, b=b)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^4 + b[1]*y
y^4 + b[1]*x

### Classic Polynomial:
# - only simple case!
# - see below for extensions;
err = x^12 - b[1]*x^9 - 3*R*x^8 + b[1]^2*x^6 + 2*R*b[1]*x^5 + 3*R^2*x^4 +
	- b[1]^3*x^3 - R*b[1]^2*x^2 - R^2*b[1]*x + b[1]^4 - R^3
round0(err)

# P6 => P12;
# when b3 == 1: P4 => P8;

### Example 2: Extended version
R = 1
b = c(2, 1)
#
sol = solve.ht4(R, b=b)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^4 + b[2]*x*y + b[1]*y
y^4 + b[2]*x*y + b[1]*x


### Example 3: Extended version
R = 1
b = c(2, 1, -1)
#
sol = solve.ht4(R, b=b)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^4 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
y^4 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Classic Polynomial:
### Simple variant:
(x^4 + b[1]*x - R) *
	(x^12 - b[1]*x^9 - 3*R*x^8 + b[1]^2*x^6 + 2*R*b[1]*x^5 + 3*R^2*x^4 +
	- b[1]^3*x^3 - R*b[1]^2*x^2 - R^2*b[1]*x + b[1]^4 - R^3)

### Extension: + b2*x*y
b1 = b[1]; b2 = b[2]; R = R[1];
x^12 - b2*x^10 - b1*x^9 + (- 3*R + b2^2)*x^8 + 2*b1*b2*x^7 + (2*R*b2 + b1^2 - b2^3)*x^6 +
	+ (2*R*b1 - 3*b1*b2^2)*x^5 + (- R*b2^2 + 3*R^2 - 3*b1^2*b2)*x^4 + b1*(-2*R*b2 + b2^3 - b1^2)*x^3 +
	+ (- R*b1^2 - R^2*b2 + 3*b1^2*b2^2)*x^2 - b1*(R^2 - 3*b1^2*b2)*x - R^3 + b1^4


### Extensions: A1-type

### Ext A1: Ex 1
R = 1
b = 2
b.ext = 1
#
sol = solve.ht4(R, b=b, b.ext=b.ext)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^4 + b[1]*y + b.ext[1]*(x + y)
y^4 + b[1]*x + b.ext[1]*(x + y)


### Ext A1: Ex 2
R = 1
b = 2
b.ext = c(-1, -1)
#
sol = solve.ht4(R, b=b, b.ext=b.ext)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^4 + b[1]*y + b.ext[1]*(x + y) + b.ext[2]*(x + y)^2
y^4 + b[1]*x + b.ext[1]*(x + y) + b.ext[2]*(x + y)^2


### Ext A1: Ex 3
R = 1
b = c(1, 2)
b.ext = c(0, -1)
#
sol = solve.ht4(R, b=b, b.ext=b.ext)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^4 + b[1]*y + b[2]*x*y + b.ext[1]*(x + y) + b.ext[2]*(x + y)^2
y^4 + b[1]*x + b[2]*x*y + b.ext[1]*(x + y) + b.ext[2]*(x + y)^2


#############
### Shift ###

### (x - s)^4 + b*y

# (x-s)^4 + b1*y = R
# (y-s)^4 + b1*x = R

# simple approach:
# x - s => x; y - s => y;
# S => S - 2*s;
# R => R - b1*s;

### Solution:

# Trivial solution: x = y;

# Diff =>
# x*y = (S*(S^2 - 4*s*S + 6*s^2) - 4*s^3 - b1) / (2*S - 4*s)

### Sum =>
S^6 - 12*s*S^5 + 60*s^2*S^4 - 4*(b1 + 40*s^3)*S^3 + 4*(R + 5*b1*s + 60*s^4)*S^2 +
	- 16*(R*s + 2*b1*s^2 + 12*s^5)*S + 16*(R*s^2 + b1*s^3 + 4*s^6) - b1^2

### Solution:
solve.ht_sh.S2P4 = function(R, b, s=0) {
	coeff = c(1, -12*s, 60*s^2, -4*(b[1] + 40*s^3), 4*(R + 5*b[1]*s + 60*s^4),
		-16*(R*s + 2*b[1]*s^2 + 12*s^5), 16*(R*s^2 + b[1]*s^3 + 4*s^6) - b[1]^2)
	x.sum = roots(coeff)
	xy = (x.sum*(x.sum^2 - 4*s*x.sum + 6*s^2) - 4*s^3 - b[1]) / (2*x.sum - 4*s)
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x=x, y=y)
	return(rbind(sol, sol[,2:1]))
}

### Example
R = 1
b = 3
s = 1
#
sol = solve.ht_sh.S2P4(R, b, s)
x = sol[,1]; y = sol[,2];
sol

### Test
(x - s)^4 + b[1]*y
(y - s)^4 + b[1]*x

### Classical Polynomial: P12
round0.p(poly.calc(sol[,1]))
# shifted back
round0.p(poly.calc(sol[,1] - s))

x = x - s
err = (x^12 - b[1]*x^9 - 3*(R - b[1]*s)*x^8 + b[1]^2*x^6 - 2*b[1]*(b[1]*s - R)*x^5 +
	+ 3*(b[1]*s - R)^2*x^4 - b[1]^3*x^3 + b[1]^2*(b[1]*s - R)*x^2 +
	- b[1]*(b[1]*s - R)^2*x + (b[1]*s - R)^3 + b[1]^4)
round0(err)


################

################
### Variants ###
################

### Extended Variant:

# x^4 + b3*y^3 + b2*y^2 + b1*y = R
# y^4 + b3*x^3 + b2*x^2 + b1*x = R


### Solution:

### Diff =>
S^3 - b3*S^2 - b2*S - b1 - x*y*(2*S-b3) # = 0

### Sum =>
S^4 - 4*x*y*S^2 + 2*(x*y)^2 + b3*S^3 - 3*b3*x*y*S + b2*S^2 - 2*b2*x*y + b1*S - 2*R # = 0

### Eq S:
S^6 - b3*S^5 - 2*(b3^2 + b2)*S^4 + (b3^3 - 4*b1 - 4*b3*b2)*S^3 +
	+ (4*R - b1*b3 + 2*b3^2*b2 - 3*b2^2)*S^2 +
	+ (b3^2*b1 + b3*b2^2 - 4*b3*R - 4*b1*b2)*S + b3^2*R - b1^2 + b1*b2*b3 # = 0

### Auxiliary Eq:
# x*y*(2*S-b3) = S^3 - b3*S^2 - b2*S - b1


### Solver:
solve.S2P4Full = function(R, b, debug=TRUE) {
	b1 = b[1]; b2 = b[2]; b3 = b[3];
	coeff = c(1, - b3, - 2*b3^2 - 2*b2, b3^3 - 4*b1 - 4*b3*b2,
			4*R - b3*b1 + 2*b3^2*b2 - 3*b2^2, b3^2*b1 + b3*b2^2 - 4*b3*R - 4*b1*b2,
			b3^2*R - b1^2 + b1*b2*b3);
	S = roots(coeff);
	if(debug) print(S);
	xy = (S^3 - b3*S^2 - b2*S - b1) / (2*S - b3);
	xy.d = sqrt(S^2 - 4*xy + 0i)
	x = (S + xy.d) / 2;
	y = (S - xy.d) / 2;
	sol = cbind(x=x, y=y);
	sol = rbind(sol, sol[,2:1]);
	return(sol);
}

### Examples:

R = -1;
b = c(3,3,-1)
sol = solve.S2P4Full(R, b);
x = sol[,1]; y = sol[,2];

### Test:
round0(x^4 + b[3]*y^3 + b[2]*y^2 + b[1]*y - R) # == 0
round0(y^4 + b[3]*x^3 + b[2]*x^2 + b[1]*x - R) # == 0


### Classic Polynomial
# - see derivation;

###################
###################

### Simple Variants
### x^4 + b*x*y

# x^4 + b1*x*y = R
# y^4 + b1*x*y = R

### Solution:

# Trivial solutions: x = y *OR* x = -y;

### Diff =>
# => x =  y *OR*
#    x = -y *OR* 2*x*y = (x+y)^2;
# Case: x != +/- y =>
# 2*x*y = S^2

### Sum =>
S^4 - 2*b1*S^2 + 4*R


### Example
b = 3
R = 1
#
x.sum = roots(c(1, 0, -2*b[1], 0, 4*R))
xy = x.sum^2/2
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol

### Test
x^4 + b[1]*x*y
y^4 + b[1]*x*y

### Classical Polynomial
round0.p(poly.calc(sol[,1]))

# very trivial P8
R^2 + (b[1]^2 - 2*R)*x^4 + x^8


##############################

### Variant:
### x^4 + b2*(x*y)^2 + b1*x*y

# x^4 + b2*(x*y)^2 + b1*x*y = R
# y^4 + b2*(x*y)^2 + b1*x*y = R

### Solution:

# Trivial solutions: x = y *OR* x = -y;

### Diff =>
# same Diff:
# 2*x*y = Z^2

### Sum =>
(b2 - 1)*Z^4 + 2*b1*Z^2 - 4*R


### Example
b = c(3,2)
R = 1
#
x.sum = roots(c(b[2]-1, 0, 2*b[1], 0, -4*R))
xy = x.sum^2/2
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol

### Test
x^4 + b[2]*(x*y)^2 + b[1]*x*y
y^4 + b[2]*(x*y)^2 + b[1]*x*y

### Classical Polynomial
round0.p(poly.calc(sol[,1]))

# trivial P8;


###################################
###################################
###################################

###############
### Order 5 ###
###############

##############
### Simple ###

### x^5 + b*y

# x^5 + b*y = R
# y^5 + b*x = R

### Solution

### Diff =>
S^4 - 3*x*y*S^2 + (x*y)^2 - b # = 0

### Sum =>
S*(S^4 - 5*x*y*S^2 + 5*(x*y)^2 + b) - 2*R # = 0

### Auxiliary Eq:
# 5*x*y*S^3 = 2*S^5 - 3*b*S + R

### Eq S:
S^10 - 8*b*S^6 + 11*R*S^5 - 9*b^2*S^2 + 6*b*R*S - R^2 # = 0


### Solver:
solve.S2Ht.P5Basic = function(R, b, doPoly=TRUE) {
	x.sum = roots(c(1, 0,0,0, -8*b[1], 11*R, 0,0, - 9*b[1]^2, 6*b[1]*R, - R^2))
	xy = (2*x.sum^5 - 3*b[1]*x.sum + R) / x.sum^3 / 5
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(round0(x), round0(y))
	sol = rbind(sol, sol[,2:1]);
	p = NULL;
	if(doPoly) p = round0.p(poly.calc(sol[,1]));
	return(list(sol=sol, p=p))
}

### Examples:
R = 1
b = 3
#
sol = solve.S2Ht.P5Basic(R, b)
x = sol$sol[,1]; y = sol$sol[,2]
sol$p

### Test
x^5 + b[1]*y
y^5 + b[1]*x


### Ex 2:
R = 1
b = -1
#
sol = solve.S2Ht.P5Basic(R, b)
x = sol$sol[,1]; y = sol$sol[,2]
sol$p

### Test
x^5 + b[1]*y
y^5 + b[1]*x


### Ex 3:
R = -2
b = -3
#
sol = solve.S2Ht.P5Basic(R, b)
x = sol$sol[,1]; y = sol$sol[,2]
sol$p

### Test
x^5 + b[1]*y
y^5 + b[1]*x


### Classic Polynomial
# P[5] * P[20]
x^20 - b[1]*x^16 - 4*R*x^15 + b[1]^2*x^12 + 3*R*b[1]*x^11 + 6*R^2*x^10 - b[1]^3*x^8 +
 - 2*R*b[1]^2*x^7 - 3*R^2*b[1]*x^6 - 4*R^3*x^5 + b[1]^4*x^4 + R*b[1]^3*x^3 +
 + R^2*b[1]^2*x^2 + R^3*b[1]*x + R^4 - b[1]^5


#########
### Ex 3:
# R^4 - b[1]^5 == 0
R = -1
b = 1
#
sol = solve.S2Ht.P5Basic(R, b)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^5 + b[1]*y
y^5 + b[1]*x
# x*(x+1) * P18
err = -1 + 2*x - 3*x^2 + 4*x^3 - 3*x^5 + 5*x^6 - 6*x^7 + 6*x^8 - 3*x^10 + 4*x^11 - 4*x^12 +
	+ 4*x^13 - x^15 + x^16 - x^17 + x^18
round0(err)


### Derived System:
# - based on the conjugated roots;
R = -1
b = 1
#
sol = solve.S2Ht.P5Basic(R, b)$sol;
sol = sol[isConj.f(sol[,1], sol[,2]),];
# 4 of the roots:
x = Re(sol[,1]); y = Im(sol[,1]);

x^5 - 10*x^3*y^2 + 5*x*y^4 + b*x - R # = 0
y^4 - 10*x^2*y^2 + 5*x^4 - b # = 0

### Sum(Eq 1 + y*Eq 2) => Original system
(x + y*1i)^5 + b*(x-y*1i) - R
(x - y*1i)^5 + b*(x+y*1i) - R

# Classic Polynomial =>
12*x^5 - 20*x^3*y^2 - 3*b*x + R/2
(88*x^5 + 3*b*x - R/2)^2 - 400*x^6*(20*x^4 + b)
2^10*x^10 - 2^9*b*x^6 + 2^5*11*R*x^5 - 4*9*b^2*x^2 + 4*3*b*R*x - R^2
x = 2*x # re-scaled
x^10 - 8*b*x^6 + 11*R*x^5 - 9*b^2*x^2 + 6*b*R*x - R^2
# original Eq for S:
# S^10 - 8*b*S^6 + 11*R*S^5 - 9*b^2*S^2 + 6*b*R*S - R^2


#####################

#####################
### x^5 + b*(x+y) = R

# x^5 + b1*(x+y) = R
# y^5 + b1*(y+x) = R

### Solution:
# [trivial]

m5 = unity(5, all=F)
m4 = unity(4, all=T)

### Diff =>
x^5 - y^5 # = 0
# y = x * m5^j

### =>
x^5 + b1*(m5+1)*x # = R
# variant for R = 0:
x^4 + b1*(m5+1) # = 0

### Solver:
m5 = unity(5, all=F)
solve.S2P5Trivial = function(R, b) {
	# does NOT include case x == y;
	x = sapply(1:4, function(id) roots(c(1, 0,0,0, b[1]*(m5^id + 1), -R[1])));
	x = as.vector(x);
	y = x * rep(m5^(1:4), each=5);
	sol = cbind(x, y);
	return(sol);
}

### Example
R = 0
b = 1
#
sol = solve.S2P5Trivial(R, b);
x = sol[,1]; y = sol[,2];
sol

### Test
round0(x^5 + b[1]*(x+y))
round0(y^5 + b[1]*(y+x))

### Classic Polynomial
round0.p(poly.calc(sol[,1]))

err = 1 + 4*x^4 + 12*x^8 + 22*x^12 + 30*x^16 + 28*x^20 + 17*x^24 + 6*x^28 + x^32
round0(err)
# based on: (1 + 2*x + 4*x^2 + 3*x^3 + x^4)^2


### Derived System:
x^5 - 10*x^3*y^2 + 5*x*y^4 + 2*b*x - R # = 0
y^4 - 10*x^2*y^2 + 5*x^4 # = 0

R = -1
b = 1
#
sol = solve.S2P5Trivial(R, b);
isConj = isConj.f(sol[,1], sol[,2]);
# 4 of the solutions of the Derived system
sol = sol[isConj, ];
x = Re(sol[,1]); y = Im(sol[,1]);

### Test:
x^5 - 10*x^3*y^2 + 5*x*y^4 + 2*b*x # - R
round0(y^4 - 10*x^2*y^2 + 5*x^4) # = 0


### Classic Poly (Derived system):
24*x^5 - 40*x^3*y^2 - 2*b*x + R # = 0
# 2 variants:
16*(11 + 5*sqrt(5))*x^5 + 2*b*x - R # = 0
16*(11 - 5*sqrt(5))*x^5 + 2*b*x - R # = 0
# a somehow "cyclic" dependency;


### Example 2:
# [one of the variants]
R = 16*(11 + 5*sqrt(5));
b = - 8*(11 + 5*sqrt(5));
sol = solve.S2P5Trivial(R, b);
isConj = isConj.f(sol[,1], sol[,2]);
# actually only 2 of the solutions:
sol = sol[isConj, ];
x = Re(sol[,1]); y = Im(sol[,1]);

err = x^5 - x - 1 # = 0
round0(err)


###############################
### x^5 + b2*x*y + b1*(x+y) = 0

# x^5 + b2*x*y + b1*(x+y) = 0
# y^5 + b2*x*y + b1*(x+y) = 0

### Solution

m5 = unity(5, all=F)
m4 = unity(4, all=T)

### Diff =>
# x^5 - y^5 = 0
# y = x * m5^j

### =>
# x^5 + b2*m5*x^2 + b1*(m5+1)*x = 0
# x^4 + b2*m5*x + b1*(m5+1) = 0

solve.ht5Zero = function(b) {
	x = as.vector(sapply(1:4, function(id) roots(c(1,0,0, b[2]*m5^id, b[1]*(m5^id + 1)))))
	y = x * m5^rep(1:4, each=4)
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	sol
}

### Example
b = c(1, 1)
#
sol = solve.ht5Zero(b)
x = sol[,1]; y = sol[,2]
sol

### Test
round0(x^5 + b[2]*x*y + b[1]*(x+y))
round0(y^5 + b[2]*x*y + b[1]*(y+x))

### Classic Polynomial
round0.p(poly.calc(sol[,1]))
# (P16)^2

err = (1 + 2*x + 4*x^2 + 3*x^3 + 3*x^4 - 2*x^5 - x^6 - x^7 + 4*x^8 - x^9 + x^10 + 3*x^12 -  x^13 + x^16)^2
round0(err)

### TODO

### Example 2:
b = c(-1, 2)
#
sol = solve.ht5Zero(b)
x = sol[,1]; y = sol[,2]
sol

### Test
round0(x^5 + b[2]*x*y + b[1]*(x+y))
round0(y^5 + b[2]*x*y + b[1]*(y+x))

### Classic Polynomial
round0.p(poly.calc(sol[,1]))

err = (1 - 4*x + 16*x^2 - 24*x^3 + 14*x^4 - 4*x^5 + 4*x^6 - 8*x^7 + 4*x^8 + 2*x^9 + 4*x^10 - 3*x^12 - 2*x^13 + x^16)^2
round0(err)


####################
####################

### x^5 + b*y^4
# - efficient reduction;

# x^5 + b*y^4 = R
# y^5 + b*x^4 = R

### Solution

### Diff =>
S^4 - 3*x*y*S^2 + (x*y)^2 - b*S*(S^2 - 2*x*y) # = 0
S^4 - b*S^3 + (x*y)^2 - 3*x*y*S^2 + 2*b*x*y*S # = 0

### Diff 2:
# y*Eq 1 - x*Eq 2 =>
x*y*(x^4 - y^4) - b*(x^5 - y^5) + R*(x - y) # = 0
(x*y - b^2)*(x^4 - y^4) + R*(x - y) # = 0
(x*y - b^2)*S*(S^2 - 2*x*y) + R # = 0

### Sum =>
S*(S^4 - 5*x*y*S^2 + 5*(x*y)^2) + b*(S^4 - 4*x*y*S^2 + 2*(x*y)^2) - 2*R # = 0

### Eq S:
S^10 - b*S^9 - 3*b^2*S^8 + 2*b^3*S^7 + b^4*S^6 + 11*R*S^5 - 18*b*R*S^4 + 4*b^2*R*S^3 +
	+ 4*b^3*R*S^2 - R^2 # = 0


### Solver:

solve.S2Ht.P5Y4 = function(R, b, debug=TRUE, all=TRUE) {
	coeff = c(1, - b, - 3*b^2, 2*b^3, b^4, 11*R, - 18*b*R,
		4*b^2*R, 4*b^3*R, 0, - R^2);
	S = roots(coeff);
	if(debug) print(S);
	xy = (2*S^5 - 2*b*S^4 - b^2*S^3 + R) / (5*S^3 - 4*b*S^2 - 2*b^2*S);
	d = sqrt(S^2 - 4*xy + 0i);
	x = (S + d)/2;
	y = S - x;
	sol = cbind(x, y);
	if(all) sol = rbind(sol, sol[, c(2,1)]);
	return(sol);
}
test.S2Ht.P5Y4 = function(sol, b, R=NULL) {
	x = sol[,1]; y = sol[,2];
	err1 = x^5 + b*y^4;
	err2 = y^5 + b*x^4;
	err = rbind(err1, err2);
	err = round0(err);
	return(err);
}

### Examples:

R = 2;
b = -3;
sol = solve.S2Ht.P5Y4(R, b)

test.S2Ht.P5Y4(sol, b=b)


### Ex 2:
R = -3;
b = -1;
sol = solve.S2Ht.P5Y4(R, b)

test.S2Ht.P5Y4(sol, b=b)


### Classic Polynomial:
x = sol[,1];
x^20 - b*x^19 + b^2*x^18 - b^3*x^17 + b^4*x^16 - (b^5 + 4*R)*x^15 + b*(b^5 + 3*R)*x^14 +
	- b^7*x^13 - 2*b^2*R*x^13 + b^8*x^12 + b^3*R*x^12 - b^5*R*x^10 + 6*R^2*x^10 +
	+ 2*b^6*R*x^9 - 3*b*R^2*x^9 - 3*b^7*R*x^8 + b^2*R^2*x^8 - b^5*R^2*x^5 - 4*R^3*x^5 +
	+ 3*b^6*R^2*x^4 + b*R^3*x^4 - b^5*R^3 + R^4 # = 0


########################
########################

### Extended Side-Chains
### + b2*y^2 + b1*y

# x^5 + b2*y^2 + b1*y = R
# y^5 + b2*x^2 + b1*x = R

### Solution

### Diff =>
S^4 - 3*x*y*S^2 + (x*y)^2 - b2*S - b1 # = 0

### Diff(x*...) =>
# - with efficient reduction;
x*y*S^3 - 2*(x*y)^2*S + b2*x*y - b2*S^2 - b1*S + R # = 0
# Reduction =>
2*S^5 - 5*x*y*S^3 + b2*x*y - 3*b2*S^2 - 3*b1*S + R # = 0

### Eq S:
S^10 - 4*b2*S^7 - 8*b1*S^6 + 11*R*S^5 - 11*b2^2*S^4 - 19*b1*b2*S^3 +
	- 3*(3*b1^2 - b2*R)*S^2 + (6*b1*R + b2^3)*S - R^2 + b1*b2^2 # = 0

### Solver:

solve.S2Ht.P5Y21 = function(R, b, debug=TRUE, all=TRUE) {
	b1 = b[1]; b2 = b[2];
	coeff = c(1, 0, 0, - 4*b2, - 8*b1, 11*R, - 11*b2^2, - 19*b1*b2,
		- 9*b1^2 + 3*b2*R, 6*b1*R + b2^3, b1*b2^2 - R^2);
	S = roots(coeff);
	if(debug) print(S);
	xy = (2*S^5 - 3*b2*S^2 - 3*b1*S + R) / (5*S^3 - b2);
	d = sqrt(S^2 - 4*xy + 0i);
	x = (S + d)/2;
	y = S - x;
	sol = cbind(x, y);
	if(all) sol = rbind(sol, sol[, c(2,1)]);
	return(sol);
}
test.S2Ht.P5Y21 = function(sol, b, R=NULL) {
	x = sol[,1]; y = sol[,2];
	err1 = x^5 + b[2]*y^2 + b[1]*y;
	err2 = y^5 + b[2]*x^2 + b[1]*x;
	err = rbind(err1, err2);
	err = round0(err);
	return(err);
}

### Examples:

R = 2;
b = c(-3, -1);
sol = solve.S2Ht.P5Y21(R, b)

test.S2Ht.P5Y21(sol, b=b)


### Ex 2:
R = -3;
b = c(1, -4);
sol = solve.S2Ht.P5Y21(R, b)

test.S2Ht.P5Y21(sol, b=b)


########################
########################

### Extra-Extended Side-Chains

# x^5 + b4*x^4 + b3*y^3 + b2*y^2 + b1*y = R
# y^5 + b4*y^4 + b3*x^3 + b2*x^2 + b1*x = R

### Solution

### Diff =>
S^4 + (x*y)^2 - 3*x*y*S^2- 2*b4*x*y*S + b3*x*y + b4*S^3 - b3*S^2 - b2*S - b1 # = 0

### Diff(x*...) - (x-y) * Sum(...) =>
# - with efficient reduction & simple Reduction =>
(- 5*S^3 - 6*b4*S^2 - 2*b4^2*S + 4*b3*S + b3*b4 + b2)*x*y +
	+ 2*S^5 + 3*b4*S^4 + (b4^2 - 3*b3)*S^3 - (b3*b4 + 3*b2)*S^2 - (b2*b4 + 3*b1)*S - b1*b4 + R # = 0

### Eq S:
# - see coefficients in the solver
#   or in file Poly.System.Hetero.Symmetric.Derivation.R;


### Solver:

solve.S2Ht.P5Ch4 = function(R, b, debug=TRUE, all=TRUE) {
	if(length(b) < 4) {
		warning("Missing b coefficients! Set to 0.");
		b = c(b, rep(0, 4 - length(b)));
	}
	b1 = b[1]; b2 = b[2]; b3 = b[3]; b4 = b[4];
	coeff = c(1, 4*b4, - 2*b3 + 6*b4^2, - 4*b2 - 4*b3*b4 + 4*b4^3,
		- 8*b1 - 8*b2*b4 - 6*b3^2 - 4*b3*b4^2 + b4^4,
		11*R - 17*b1*b4 - 14*b2*b3 - 7*b2*b4^2 - 6*b3^2*b4 - b3*b4^3,
		22*b4*R - 7*b1*b3 - 14*b1*b4^2 - 11*b2^2 - 12*b2*b3*b4 - 2*b2*b4^3 + 4*b3^3 - 2*b3^2*b4^2,
		- (11*b3 - 16*b4^2)*R - 19*b1*b2 + 9*b2*b3^2 - (2*b1*b3 + 10*b2^2 - b3^3)*b4 +
			- 4*b2*b3*b4^2 - 4*b1*b4^3,
		3*b2*R - 15*b3*b4*R + 4*b4^3*R - 9*b1^2 - 15*b1*b2*b4 + 4*b1*b3^2 - b1*b3*b4^2 +
			+ 6*b2^2*b3 - 3*b2^2*b4^2 + 2*b2*b3^2*b4,
		(6*b1 + 4*b3^2 - 4*b3*b4^2)*R - (6*b1^2 - b1*b3^2 - b2^2*b3)*b4 + 5*b1*b2*b3 - 4*b1*b2*b4^2 + b2^3,
		- R^2 + 2*b1*b4*R + b2*b3*R + b3^2*b4*R - b1^2*b4^2 + b1*b2^2 + b1*b2*b3*b4);
	S = roots(coeff);
	if(debug) print(S);
	xy = (2*S^5 + 3*b4*S^4 - (3*b3 - b4^2)*S^3 - (3*b2 + b3*b4)*S^2 - (3*b1 + b2*b4)*S + R - b1*b4);
	xy = xy / (5*S^3 + 6*b4*S^2 - (4*b3 - 2*b4^2)*S - b2 - b3*b4);
	d = sqrt(S^2 - 4*xy + 0i);
	x = (S + d)/2;
	y = S - x;
	sol = cbind(x, y);
	if(all) sol = rbind(sol, sol[, c(2,1)]);
	return(sol);
}
test.S2Ht.P5Ch4 = function(sol, b, R=NULL) {
	x = sol[,1]; y = sol[,2];
	if(length(b) < 4) {
		warning("Missing b coefficients! Set to 0.");
		b = c(b, rep(0, 4 - length(b)));
	}
	err1 = x^5 + b[4]*x^4 + b[3]*y^3 + b[2]*y^2 + b[1]*y;
	err2 = y^5 + b[4]*y^4 + b[3]*x^3 + b[2]*x^2 + b[1]*x;
	err = rbind(err1, err2);
	err = round0(err);
	return(err);
}

### Examples:

R = 2;
b = c(-3, -1, 2, -1);
sol = solve.S2Ht.P5Ch4(R, b)

test.S2Ht.P5Ch4(sol, b=b)


### Ex 2:
R = -3;
b = c(1, -4, 5, -2);
sol = solve.S2Ht.P5Ch4(R, b)

test.S2Ht.P5Ch4(sol, b=b)


#######################
#######################

### Shifted Side-Chains
### x^5 + b3*x^3*y^2 + b2*x^2*y

# x^5 + b3*x^3*y^2 + b2*x^2*y = R
# y^5 + b3*y^3*x^2 + b2*y^2*x = R

### Solution

### Diff =>
S^4 - 3*x*y*S^2 + (b3 + 1)*(x*y)^2 + b2*x*y # = 0

### Diff 2:
# y*Eq 1 - x*Eq 2 =>
x*y*(x^4 - y^4) + R*(x - y) # = 0
x*y*S*(S^2 - 2*x*y) + R # = 0


### Eq S:
(b3 - 1)*S^10 + 2*b2*S^8 + (7*b3 - 11)*R*S^5 - b2*(b3 - 11)*R*S^3 +
	- 2*b2^2*R*S + (b3 + 1)^2*R^2 # = 0


### Solver:

solve.S2Ht.P5ChShift = function(R, b, debug=TRUE, all=TRUE) {
	if(length(b) < 2) stop("Wrong parameter b!");
	coeff = coeff.S2Ht.P5ChShift(R, b=b);
	S = roots(coeff);
	if(debug) print(S);
	if(length(b) >= 3) {
		b1 = b[1]; b2 = b[2]; b3 = b[3];
	} else {
		b1 = 0; b2 = b[1]; b3 = b[2];
	}
	xy = (2*S^5 + 2*b1*S + (b3 + 1)*R) / (5*S^3 - b3*S^3 - 2*b2*S);
	d = sqrt(S^2 - 4*xy + 0i);
	x = (S + d)/2;
	y = S - x;
	sol = cbind(x, y);
	if(all) sol = rbind(sol, sol[, c(2,1)]);
	return(sol);
}
coeff.S2Ht.P5ChShift = function(R, b) {
	if(length(b) == 2) {
		b2 = b[1]; b3 = b[2];
		coeff = c(b3 - 1, 0, 2*b2, 0, 0, (7*b3 - 11)*R, 0,
			- b2*(b3 - 11)*R, 0, - 2*b2^2*R, (b3 + 1)^2*R^2);
		return(coeff);
	}
	b1 = b[1]; b2 = b[2]; b3 = b[3];
	coeff = c(b3 - 1, 0, 2*b2, 0, b1*(b3 + 3), (7*b3 - 11)*R, 2*b1*b2,
		- b2*(b3 - 11)*R, 4*b1^2, (4*b1*(b3 + 1) - 2*b2^2)*R, (b3 + 1)^2*R^2);
	return(coeff);
}
test.S2Ht.P5ChShift = function(sol, b, R=NULL) {
	x = sol[,1]; y = sol[,2];
	if(length(b) >= 3) {
		b1 = b[1]; b2 = b[2]; b3 = b[3];
	} else {
		b1 = 0; b2 = b[1]; b3 = b[2];
	}
	err1 = x^5 + b3*x^3*y^2 + b2*x^2*y + b1*x;
	err2 = y^5 + b3*y^3*x^2 + b2*y^2*x + b1*y;
	err = rbind(err1, err2);
	isNum = ! is.nan(err)
	err[isNum] = round0(err[isNum]);
	return(err);
}

### Examples:

R = 2;
b = c(-2, -3);
sol = solve.S2Ht.P5ChShift(R, b)

test.S2Ht.P5ChShift(sol, b=b)


### Ex 2:
R = 3;
b = c(-1, 3);
sol = solve.S2Ht.P5ChShift(R, b)

test.S2Ht.P5ChShift(sol, b=b)


### Ex 3: Special Case
# TODO
R = 3;
b = c(-2, -1);
sol = solve.S2Ht.P5ChShift(R, b)

test.S2Ht.P5ChShift(sol, b=b)


### Ex 4:
R = -5;
b = c(-1, 2, -2);
sol = solve.S2Ht.P5ChShift(R, b)

test.S2Ht.P5ChShift(sol, b=b)


### Classic Polynomial:
b2 = b[1]; b3 = b[2]; x = sol[,1];
(b3^3 - b3^2 - b3 + 1)*x^20 - b2*(b3 - 1)^2*x^18 + b2^2*(b3 + 1)*x^16 +
	- (3*b3^3 - 4*b3^2 - 3*b3 + 4)*R*x^15 - b2^3*x^14 +
	- b2*(b3^3 - 4*b3^2 + 4*b3 - 3)*R*x^13 + b2^2*(b3^2 - 5*b3 - 2)*R*x^11 +
	+ (6 - 3*b3 - 5*b3^2 + 2*b3^3 + b3^4)*R^2*x^10 - b2^3*(b3 - 1)*R*x^9 +
	- b2*(3 - 2*b3 + 3*b3^2)*R^2*x^8 + b2^4*R*x^7 + b2^2*(4*b3 + 1)*R^2*x^6 +
	- 4*R^3*x^5 + b3*R^3*x^5 + 2*b3^2*R^3*x^5 +
	+ b2*R^3*x^3 + R^4 # = 0


##############
### Extension: + b1*x

# x^5 + b3*x^3*y^2 + b2*x^2*y + b1*x = R

### Diff =>
S^4 - 3*x*y*S^2 + (b3 + 1)*(x*y)^2 + b2*x*y + b1 # = 0

### Eq S:
(b3 - 1)*S^10 + 2*b2*S^8 + b1*(b3 + 3)*S^6 + (7*b3 - 11)*R*S^5 + 2*b1*b2*S^4 +
	- b2*(b3 - 11)*R*S^3 + 4*b1^2*S^2 + (4*b1*(b3 + 1) - 2*b2^2)*R*S +
	+ (b3 + 1)^2*R^2 # = 0

##############

##############
### Extension: + b1*y

# x^5 + b3*x^3*y^2 + b2*x^2*y + b1*y = R

### System Transform:
# b1 => - b1; (Step 1)
# R  => R - b1*S; (Step 2)

### Eq S:
(b3 - 1)*S^10 + 2*b2*S^8 - b1*(b3 + 3)*S^6 + (7*b3 - 11)*(R - b1*S)*S^5 - 2*b1*b2*S^4 +
	- b2*(b3 - 11)*(R - b1*S)*S^3 + 4*b1^2*S^2 - (4*b1*(b3 + 1) + 2*b2^2)*(R - b1*S)*S +
	+ (b3 + 1)^2*(R - b1*S)^2 # = 0
# =>
(b3 - 1)*S^10 + 2*b2*S^8 - 8*b1*(b3 - 1)*S^6 + (7*b3 - 11)*R*S^5 + b1*b2*(b3 - 13)*S^4 +
	- b2*(b3 - 11)*R*S^3 + (b1^2*(b3 + 3)^2 + 2*b1*b2^2)*S^2 +
	- 2*(b1*(b3+1)*(b3 + 3) + b2^2)*R*S + (b3 + 1)^2*R^2 # = 0

### Solver:

solve.S2Ht.P5ChShiftY = function(R, b, debug=TRUE, all=TRUE) {
	if(length(b) < 3) stop("Wrong parameter b!");
	coeff = coeff.S2Ht.P5ChShiftY(R, b=b);
	S = roots(coeff);
	if(debug) print(S);
	b1 = b[1]; b2 = b[2]; b3 = b[3];
	xy = (2*S^5 - b1*(b3 + 3)*S + (b3 + 1)*R) / (5*S^3 - b3*S^3 - 2*b2*S);
	d = sqrt(S^2 - 4*xy + 0i);
	x = (S + d)/2;
	y = S - x;
	sol = cbind(x, y);
	if(all) sol = rbind(sol, sol[, c(2,1)]);
	return(sol);
}
coeff.S2Ht.P5ChShiftY = function(R, b) {
	b1 = b[1]; b2 = b[2]; b3 = b[3];
	coeff = c(b3 - 1, 0, 2*b2, 0, -8*b1*(b3 - 1), (7*b3 - 11)*R, b1*b2*(b3 - 13),
		- b2*(b3 - 11)*R, (b1^2*(b3 + 3)^2 + 2*b1*b2^2),
		- 2*(b1*(b3+1)*(b3 + 3) + b2^2)*R, (b3 + 1)^2*R^2);
	return(coeff);
}
test.S2Ht.P5ChShiftY = function(sol, b, R=NULL) {
	x = sol[,1]; y = sol[,2];
	b1 = b[1]; b2 = b[2]; b3 = b[3];
	err1 = x^5 + b3*x^3*y^2 + b2*x^2*y + b1*y;
	err2 = y^5 + b3*y^3*x^2 + b2*y^2*x + b1*x;
	err = rbind(err1, err2);
	isNum = ! is.nan(err)
	err[isNum] = round0(err[isNum]);
	return(err);
}

### Examples:

R = 2;
b = c(-1, -3, 2);
sol = solve.S2Ht.P5ChShiftY(R, b)

test.S2Ht.P5ChShiftY(sol, b=b)


### Ex 2:
R = 3;
b = c(-5, 3, -2);
sol = solve.S2Ht.P5ChShiftY(R, b)

test.S2Ht.P5ChShiftY(sol, b=b)


###################################
###################################
###################################

