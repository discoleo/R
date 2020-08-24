
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogenous Symmetric
###
### draft v.0.1i
### & branch v.0.2b


### Heterogenous Symmetric Polynomial Systems

### 2 Variables:
# x^n + P(x, y) = R
# y^n + P(y, x) = R


###############
### Systems ###

### 2 Variables:
# x^n + P(x, y) = R
# y^n + P(y, x) = R

# 1.) x^3 + b*y = R; [P3 => P6]
# 2.) (x - s)^3 + b*y = R; [P3 => P6: equivalent to non-shifted]
# 3.) x^3 + b*x*y = R; [P3 => trivial P6]
# 4.) (x - s)^3 + b*x*y = R; [P3 => P6]
# 5.) x^3 + b2*x*y + b1*x = R; (TODO: P6)
# 6.) x^3 + b2*x*y + b1*y = R; (TODO: P6)
# 7.) x^3 + b3*x*y + b2*y^2 + b1*y = R; (TODO: P6)
# 8.) TODO: Shift for [5-7];
# 9.) x^3 + b*(x*y)^2 = R; [P3 => P6: based on m3]
# 10.) x^3 + b2*(x*y)^2 + b1*(x*y) = R; [P4 => P8: based on m3]
# 11.) TODO: Shift for [9-10];
# 12.) x^3 + b3*x^2*y + b2*x*y^2 + b1*x = R;
### 2 High-Power Terms:
# B1.) a1*x^3 + a2*y^3 + b*x = R;
# B2.) a1*x^3 + a2*y^3 + b*x*y = R;
# B2.) TODO: a1*x^3 + a2*y^3 + b2*x*y + b1*x = R;
# B3.) TODO: Shift for [B1-3];
### Mixt-High-Power:
# M1.) x^2*y + b*x: NO solutions (x != y);
# M2.) x^2*y + b*y: trivial;
# M3.) x^2*y + b2*x^2 + b1*x: trivial (trivial P2);
### Mixt (x*y)^2
# M4.) x^2*y^2 + b2*x^2 + b1*x: simple (P2 => P4);
# M5.) x^2*y^2 + b3*x^2*y + b2*x^2 + b1*x; (TODO: P2 => P4)
# M6.) TODO: x^2*y^2 + b4*x^2*y + b3*x*y^2 + b2*x^2 + b1*x
# M7.) TODO: x^3*y^2 + b2*x^2 + b1*x = R;
### Mixt (x*y)^n
# M7.) TODO: (x*y)^3 + ...:
#      (x*y)^3 + x^3 = R;
# M8.) TODO: (x*y)^5 + b*x^3 = R
### Mixt Order 3+1:
# 19.) x^3*y + b*x: trivial;
# 20.) x^4*y + b*x; (P5 => P10)
### Order 4 & 5:
# O4.1.) x^4 + b*y = R; (P6 => P12)
# O4.2.) (x - s)^4 + b*y = R; (P6 => P12)
# O4.3.) x^4 + b*x*y = R; (trivial P8)
# O4.4.) x^4 + b2*(x*y)^2 + b1*x*y = R; (TODO: trivial P8)
# O4.5.) TODO: x^4 + b3*(x*y)^3 + b2*(x*y)^2 + b1*x*y = R;
# O5.1.) x^5 + b*(x+y) = 0; (TODO: based on (P4[x^4])^2)
# O5.2.) x^5 + b2*x*y + b1*(x+y) = 0; (TODO: based on (P16)^2)



###############
### History ###

### [branch v.0.2]
### draft v.0.2b:
# - some exploration of systems with x*y*z terms;
### branch v.0.2a:
# - more work on systems with 3 variables:
#   "proper" implementation of: x[i]^2 + b*x[k];
# - TODO: robust removal of set of wrong solutions;
#   [or avoid getting superfluous solutions ???]
### branch v.0.2a-pre-a:
# - initial work on systems with 3 variables;
# - the simple cases are less rewarding;

### [branch v.0.1]
### draft v.0.1i:
# - added x^3 + b3*x*y + b2*y^2 + b1*y = R;
# - initial exploration of: x^5 + b*(x+y) = 0 variants;
# - more classical polynomials computed;
### draft v.0.1h - v.0.1h-x:
# - initial work on:
#   x^3 + b3*x^2*y + b2*x*y^2 + b1*x = R;
# - improved formatting + some fixes + completion of some cases;
### draft v.0.1f - v.0.1g:
# - added various Mixt-High-Power variants:
#   x^4*y + b*x = R;
#   x^2*y^2 + ... variants; [v.0.1g]
### draft v.0.1e:
# - added shifted version: (x-s)^4 + b*y = R;
### draft v.0.1d:
# - added variant with 2 high-power terms:
#  -- variant 1: a1*x^3 + a2*y^3 + b*x;
#  -- variant 2: a1*x^3 + a2*y^3 + b*x*y;
### draft v.0.1c:
# - added x^3 + b1*x*y + b2*x = R;
# - added x^3 + b1*x*y + b2*y = R;
# - added x^3 + b1*(x*y)^2 = R;
# - TODO:
#  -- shifted versions;
#  -- parametric polynomials;
### draft v.0.1b - v.0.1b-x:
# - added a basic xy-type: x^3 + x*y = R;
# - added also the shift (v.0.1b-sh);
# - TODO: parametric classic polynomial [DONE: in v.0.1b-x];
### draft v.0.1a-shift:
# - derivation of the classical polynomial for shifted root;
# - more interesting polynomials are generated,
#   when the shifted root is shifted back;
#   [in general not identical to non-shifted root polynomials]
### draft v.0.1a:
# - initial version: basic heterogenous systems;



library(polynom)
library(pracma)


####################

### helper Functions

### see in the other R files:
# round0(), round0.p;

#####################


###############
### Order 3 ###

### x^3 + b*y

# x^3 + b1*y = R
# y^3 + b1*x = R

### Solution:

# Diff =>
# x^3 - y^3 - b1*(x-y) = 0
# (x - y)*(x^2 + y^2 + x*y - b1) = 0
# => x = y *OR* x^2 + y^2 + x*y - b1 = 0;
# =>
# (x+y)^2 - x*y - b1 = 0
# x*y = (x+y)^2 - b1;

# Sum =>
# (x+y)^3 - 3*x*y*(x+y) + b1*(x+y) = 2*R
# S^3 - 3*(S^2 - b1)*S + b1*S - 2*R = 0
# S^3 - 2*b1*S + R = 0

### Step 2:
# Solve:
# x + y = S
# x*y = S^2 - b1

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

### Example
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

### Derivation:
# b1*y = R - x^3
# =>
# (R - x^3)^3 / b1^3 + b1*x - R = 0
# (R - x^3)^3 + b1^4*x - R*b1^3
# (x^3 - R)^3 - b1^4*x + R*b1^3
# (x^3 + b1*x - R - b1*x)^3 - b1^4*x + R*b1^3
# (x^3 + b1*x - R)^3 - 3*b1*x*(x^3 + b1*x - R)^2 + 3*b1^2*x^2*(x^3 + b1*x - R) - b1^3*x^3 - b1^4*x + R*b1^3
# let: p = (x^3 + b1*x - R)
# p^3 - 3*b1*x*p^2 + 3*b1^2*x^2*p - b1^3*p
# p*(p^2 - 3*b1*x*p + 3*b1^2*x^2 - b1^3)
# (x^3 + b1*x - R)*(x^6 - b1*x^4 - 2*R*x^3 + b1^2*x^2 + b1*R*x + R^2 - b1^3)
# =>
(x^3 - b[1]/2 * x - R)^2 + 3/4 * b[1]^2*x^2 - b[1]^3


#####################
#####################
### Shifted Roots ###

### (x - s)^3 + b*y

# (x - s)^3 + b1*y = R
# (y - s)^3 + b1*x = R

# trivial shift (only 1 liniar non-shifted term):
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
### Diff =>
# (x - s)^3 - (y - s)^3 - b1*(x-y) = 0
# (x - y)*(x^2 + x*y + y^2 - 3*s*(x+y) + 3*s^2 - b1) = 0
# => x = y *OR* x^2 + x*y + y^2 - 3*s*(x+y) + 3*s^2 - b1 = 0;
# =>
# (x+y)^2 - 3*s*(x+y) - x*y + 3*s^2 - b1 = 0
# x*y = (x+y)^2 - 3*s*(x+y) + 3*s^2 - b1;
# x*y = Z^2 - 3*s*Z + 3*s^2 - b1;

### Sum =>
# (x - s)^3 + (y - s)^3 + b1*(x+y) = 2*R
# (x+y)^3 - 3*x*y*(x+y) - 3*s*(x^2+y^2) + (3*s^2+b1)*(x+y) - 2*s^3 - 2*R = 0
# Z^3 - 3*s*Z^2 - 3*x*y*Z + 6*s*x*y + (3*s^2+b1)*Z - 2*s^3 - 2*R
# Z^3 - 3*s*Z^2 - 3*(Z^2 - 3*s*Z + 3*s^2 - b1)*Z + 6*s*(Z^2 - 3*s*Z + 3*s^2 - b1) + (3*s^2+b1)*Z - 2*s^3 - 2*R
# -2*Z^3 + 12*s*Z^2 + 4*(b1-6*s^2)*Z + 16*s^3 - 2*R - 6*s*b1
# Z^3 - 6*s*Z^2 - 2*(b1-6*s^2)*Z - 8*s^3 + R + 3*s*b1

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
### Derivation:
# b1*y = R - (x-s)^3
# =>
# (R - (x-s)^3 - s*b1)^3 / b1^3 + b1*x - R = 0
# (R - (x-s)^3 - s*b1)^3 + b1^4*x - R*b1^3
# ((x-s)^3 - R + s*b1)^3 - b1^4*x + R*b1^3
# ((x-s)^3 - R + b1*x - b1*x + s*b1)^3 - b1^4*x + R*b1^3 # p = ((x-s)^3 - R + b1*x)
# ((x-s)^3 - R + b1*x)*(p^2 - 3*(b1*x - s*b1)*p + 3*(b1*x - s*b1)^2) - (b1*x - s*b1)^3 - b1^4*x + R*b1^3
# p*(p^2 - 3*(b1*x - s*b1)*p + 3*(b1*x - s*b1)^2) - b1^3 * ((x-s)^3 + b1*x - R)
# p*(p^2 - 3*(b1*x - s*b1)*p + 3*(b1*x - s*b1)^2 - b1^3)
# p == 0 *OR* 2nd (...) == 0, where p = ((x-s)^3 - R + b1*x)
# p^2 - 3*b1*(x - s)*p + 3*b1^2*(x - s)^2 - b1^3 == 0
# (x-s)^6 + b1^2*x^2 + R^2 + 2*b1*x*(x-s)^3 - 2*R*(x-s)^3 - 2*b1*R*x - 3*b1*(x - s)*p + 3*b1^2*(x - s)^2 - b1^3
# (x-s)^6 + 2*b1*x*(x-s)^3 - 2*R*(x-s)^3 - 3*b1*(x - s)*p + b1^2*x^2 + 3*b1^2*(x - s)^2 - 2*b1*R*x + R^2 - b1^3
# (x-s)^6 + 2*b1*x*(x-s)^3 - 2*R*(x-s)^3 - 3*b1*(x-s)^4 - 3*b1^2*(x-s)*x + 4*b1^2*x^2 + 3*b1*R*(x-s) - 6*b1^2*s*x + 3*b1^2*s^2 - 2*b1*R*x + R^2 - b1^3
# (x-s)^6 - 3*b1*(x-s)^4 + 2*b1*x*(x-s)^3 - 2*R*(x-s)^3 + b1^2*x^2 - 3*b1^2*s*x + b1*R*x - 3*b1*R*s + 3*b1^2*s^2 + R^2 - b1^3
# (x-s)^6 - 3*b1*(x-s)^4 + 2*b1*x*(x-s)^3 - 2*R*x^3 + 6*R*s*x^2 + b1^2*x^2 - 6*R*s^2*x - 3*b1^2*s*x + b1*R*x - 3*b1*R*s + 3*b1^2*s^2 + 2*R*s^3 + R^2 - b1^3
# (x-s)^6 - (b1*x - 3*b1*s)*(x-s)^3 - 2*R*x^3 + 6*R*s*x^2 + b1^2*x^2 - 6*R*s^2*x - 3*b1^2*s*x + b1*R*x - 3*b1*R*s + 3*b1^2*s^2 + 2*R*s^3 + R^2 - b1^3
# (x-s)^6 - (b1*x - 3*b1*s)*(x-s)^3 - 2*R*x^3 + (6*R*s + b1^2)*x^2 - (6*R*s^2 + 3*b1^2*s - b1*R)*x - 3*b1*R*s + 3*b1^2*s^2 + 2*R*s^3 + R^2 - b1^3
#
# (x-s)^6 - b1*x^4 + (6*b1*s - 2*R)*x^3 - (12*b1*s^2 - 6*R*s - b1^2)*x^2 + (10*b1*s^3 - 6*R*s^2 - 3*b1^2*s + b1*R)*x - 3*b1*R*s + 3*b1^2*s^2 + 2*R*s^3 - 3*b1*s^4 + R^2 - b1^3


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


########################

###################
### x^3 + b*x*y ###

# x^3 + b1*x*y = R
# y^3 + b1*x*y = R

# relatively trivial
# (x^3 - b[1]/2*x^2 - R)^2 + 3/4 * b[1]^2*x^4 = 0

### Solution:

### Diff =>
# x^3 - y^3 = 0
# (x-y)*(x^2 + x*y + y^2) = 0
# Case 2:
# x^2 + x*y + y^2 = 0
# x*y = (x + y)^2

### Sum =>
# x^3 + y^3 + 2*b1*x*y - 2*R = 0
# (x+y)^3 - 3*x*y*(x+y) + 2*b1*x*y - 2*R
# Z^3 - 3*Z^3 + 2*b1*Z^2 - 2*R
# Z^3 - b1*Z^2 + R

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


### Derivation
# b1*y = (R - x^3) / x
(b[1]*y)^3 + b[1]^4*x*y - R*b[1]^3
(R - x^3)^3 / x^3 + b[1]^3*x*(R - x^3) / x - R*b[1]^3
(R - x^3)^3 + b[1]^3*x^3*(R - x^3) - R*b[1]^3*x^3
(x^3 - R)^3 + b[1]^3*x^6
(x^3 + b[1]*x^2 - R)*((x^3 - R)^2 - b[1]*x^2*(x^3 - R) + b[1]^2*x^4)
# =>
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
# (x-s)^3 - (y-s)^3 = 0
# (x-y)*((x-s)^2 + (x-s)*(y-s) + (y-s)^2) = 0
# alternatively: x - s = (y-s)*m, where m^3 = 1;
# Case: x != y
# (x-s)^2 + (x-s)*(y-s) + (y-s)^2 = 0
# (x+y - 2*s)^2 = (x-s)*(y-s)
# x*y = (x+y - 2*s)^2 + s*(x+y) - s^2
# x*y = (Z - 2*s)^2 + s*Z - s^2
# x*y = Z^2 - 3*s*Z + 3*s^2

### Sum =>
# (x-s)^3 + (y-s)^3 + 2*b1*x*y - 2*R = 0
# (x+y)^3 - 3*x*y*(x+y) - 3*s*(x^2+y^2) + 3*s^2*(x+y) + 2*b1*x*y - 2*R - 2*s^3
# Z^3 - 3*s*(Z^2 - 2*x*y) + 3*s^2*Z - 3*x*y*Z + 2*b1*x*y - 2*R - 2*s^3
# Z^3 - 3*s*Z^2 + 3*s^2*Z - 3*x*y*Z + 6*s*x*y + 2*b1*x*y - 2*R - 2*s^3
# Z^3 - 3*s*Z^2 + 3*s^2*Z -3*Z*(Z^2 - 3*s*Z + 3*s^2) + (6*s+2*b1)*(Z^2 - 3*s*Z + 3*s^2) - 2*R - 2*s^3
# Z^3 - (6*s+b1)*Z^2 + (12*s^2+3*b1*s)*Z - 8*s^3 - 3*b1*s^2 + R


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


### Derivation:
# y - s = (x-s)*m, where m^3 = 1;
# =>
# (x-s)^3 + b1*x*(s + (x-s)*m) = R
# (x-s)^3 + b1*x*(x-s)*m + b1*s*x - R = 0
round0((x-s)^3 + b[1]*x*(x-s)*m^1 + b[1]*s*x - R)
# ((x-s)^3 + b[1]*x*(x-s)*m^1 + b[1]*s*x - R)*((x-s)^3 + b[1]*x*(x-s)*m^2 + b[1]*s*x - R)
# back-shifted: (x^3 + b[1]*x^2*m^1 - b[1]*s*x*m^2 + b[1]*s^2 - R)*(x^3 + b[1]*x^2*m^2 - b[1]*s*x*m + b[1]*s^2 - R)

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

# x-Terms vs y-Terms: are equivalent;


### x-Term
### x^3 + b2*x*y + b1*x

# x^3 + b2*x*y + b1*x = R
# y^3 + b2*x*y + b1*y = R

### Solution:

### Diff =>
# x^3 - y^3 + b1*(x - y) = 0
# (x-y)*(x^2 + x*y + y^2 + b1) = 0
# Case: x != y
# x^2 + x*y + y^2 + b1 = 0
# x*y = (x + y)^2 + b1
# x*y = Z^2 + b1

### Sum =>
# x^3 + y^3 + 2*b2*x*y + b1*(x+y) - 2*R = 0
# (x+y)^3 - 3*x*y*(x+y) + 2*b2*x*y + b1*(x+y) - 2*R
# Z^3 - 3*x*y*Z + 2*b2*x*y + b1*Z - 2*R
# Z^3 - 3*Z^3 - 3*b1*Z + 2*b2*Z^2 + 2*b2*b2 + b1*Z - 2*R
# Z^3 - b2*Z^2 + b1*Z - b1*b2 + R

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


# TODO
x^6 - b[1]*x^5 + (b[1]^2 + 2*b[2])*x^4 + ...


#######################

### y-Term
### x^3 + b2*x*y + b1*y

# x^3 + b2*x*y + b1*y = R
# y^3 + b2*x*y + b1*x = R

### Solution:

### Diff =>
# x^3 - y^3 - b1*(x - y) = 0
# (x-y)*(x^2 + x*y + y^2 - b1) = 0
# Case: x != y
# x^2 + x*y + y^2 - b1 = 0
# x*y = (x + y)^2 - b1

### Sum =>
# x^3 + y^3 + 2*b2*x*y + b1*(x+y) - 2*R = 0
# (x+y)^3 - 3*x*y*(x+y) + 2*b2*x*y + b1*(x+y) - 2*R
# Z^3 - 3*x*y*Z + 2*b2*x*y + b1*Z - 2*R
# Z^3 - 3*Z^3 + 3*b1*Z + 2*b2*Z^2 - 2*b2*b2 + b1*Z - 2*R
# Z^3 - b2*Z^2 - 2*b1*Z + b1*b2 + R

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
# TODO

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

################################
### x^3 + b3*x*y + b2*y^2 + b1*y

# x^3 + b3*x*y + b2*y^2 + b1*y = R
# y^3 + b3*x*y + b2*x^2 + b1*x = R

### Solution:

### Diff =>
# x^3 - y^3 - b2*(x^2 - y^2) - b1*(x - y) = 0
# (x-y)*(x^2 + y^2 + x*y - b2*(x+y) - b1) = 0
# Case: x != y
# x^2 + y^2 + x*y - b2*(x+y) - b1 = 0
# x*y = (x + y)^2 - b2*(x+y) - b1
# x*y = Z^2 - b2*Z - b1

### Sum =>
# S^3 - 3*x*y*S + 2*b3*x*y + b2*S^2 - 2*b2*x*y + b1*S - 2*R = 0
# S^3 - 3*S^3 + 3*b2*S^2 + 3*b1*S + 2*b3*S^2 - 2*b2*b3*S - 2*b1*b3 +
#  + b2*S^2 - 2*b2*S^2 + 2*b2^2*S + 2*b1*b2 + b1*S - 2*R = 0
# -2*S^3 + 2*b2*S^2 + 2*b3*S^2 + 4*b1*S - 2*b2*b3*S + 2*b2^2*S + 2*b1*b2 - 2*b1*b3 - 2*R
# S^3 - (b2 + b3)*S^2 - (2*b1 - b2*b3 + b2^2)*S - b1*b2 + b1*b3 + R

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
round0.p(poly.calc(sol[,1]))

### TODO: classical polynomial;


### Example 2:
b = c(-1,1,-1)
p = sapply(-6:6, function(r) print(solve.htxy(b, r)$p))
# (R + 1)^2 + 3*x^2 - 2*(R + 2)*x^3 + 3*x^4 + x^6

###
b = c(-1,-2,2)
p = sapply(-6:6, function(r) print(solve.htxy(b, r)$p))
# (R + 1)^2 - 3*(R + 4)*x + 33*x^2 - 2*(R - 16)*x^3 + 9*x^4 + x^6


################

### TODO:
# - shift;



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
# x^3 - y^3 = 0
# y = x*m, where m^3 = 1;
# Case: x != y =>
# separate equations for: m & m^2
# b1*x^4*m^2 + x^3 - R = 0
# b1*x^4*m + x^3 - R = 0

solve.xysq = function(b, R, isInverse=FALSE) {
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
x = sol[,1]
err = b[1]^2*x^8 - b[1]*x^7 + x^6 + b[1]*R*x^4 - 2*R*x^3 + R^2
round0(err)

round0.p(poly.calc(x))
err = round0(4 - 8*x^3 + 2*x^4 + 4*x^6 - 2*x^7 + x^8)
err


####################
####################

### x^3 + b2*(x*y)^2 + b1*(x*y)

# x^3 + b2*(x*y)^2 + b1*(x*y) = R
# y^3 + b2*(x*y)^2 + b1*(x*y) = R

### Solution:

m3 = unity(3, all=F)

### Diff =>
# x^3 - y^3 = 0
# y = x*m, where m^3 = 1;
# Case: x != y =>
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

### x^3 + b3*x^2*y + b2*x*y^2 + b1*x

# x^3 + b3*x^2*y + b2*x*y^2 + b1*x = R
# y^3 + b3*y^2*x + b2*y*x^2 + b1*y = R

### Solution:

### Diff =>
# x^3 - y^3 + (b3 - b2)*x*y*(x-y) + b1*(x-y) = 0
# (x-y)*(x^2 + y^2 + x*y + (b3 - b2)*x*y + b1) = 0
# Case: x != y =>
# x^2 + y^2 + x*y + (b3 - b2)*x*y + b1 = 0
# (x+y)^2 + (b3 - b2 - 1)*x*y + b1
# Z^2 + (b3 - b2 - 1)*x*y + b1
# (b3 - b2 - 1)*x*y = - Z^2 - b1
# (b2 - b3 + 1)*x*y = Z^2 + b1 

### Sum =>
# x^3 + y^3 + (b3 + b2)*x*y*(x+y) + b1*(x+y) = 2*R
# Z^3 - 3*x*y*Z + (b2 + b3)*x*y*Z + b1*Z - 2*R = 0
# Z^3 + b1*Z - 2*R + (b2 + b3 - 3)*x*y*Z = 0
# Z^3 + b1*Z - 2*R + (b2 + b3 - 3)*Z*(Z^2 + b1)/(b2 - b3 + 1)
# (b2 - b3 + 1)*(Z^3 + b1*Z - 2*R) + (b2 + b3 - 3)*Z*(Z^2 + b1)
# 2*(b2 - 1)*Z^3 + 2*b1*Z*(b2 - 1) - 2*R*(b2 - b3 + 1)
# (b2 - 1)*Z^3 + b1*Z*(b2 - 1) - R*(b2 - b3 + 1)

### TODO:
# Case: b2 - b3 + 1 == 0;
# Case: b2 == 1;

### Example:
b = c(3, 2, 1)
R = 1
#
x.sum = roots(c(1, 0, b[1], - R*(b[2] - b[3] + 1)/(b[2] - 1)))
xy = (x.sum^2 + b[1]) / (b[2] - b[3] + 1)
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol

### Test
x^3 + b[3]*x^2*y + b[2]*x*y^2 + b[1]*x
y^3 + b[3]*y^2*x + b[2]*y*x^2 + b[1]*y

### Classic Polynomial
round0.p(poly.calc(sol[,1]))

### TODO


#########################
#########################

#########################
### High-Power Terms: > 1

### a1*x^3 + a2*y^3 + b1*x

# a1*x^3 + a2*y^3 + b1*x = R
# a2*x^3 + a1*y^3 + b1*y = R

### Solution

### Diff =>
# (a1 - a2)*(x^3 - y^3) = -b1*(x - y)
# (x - y)*(x^2 + y^2 + x*y + b1/(a1 - a2)) = 0
# Case x != y:
# S^2 - x*y + b1/(a1 - a2) = 0
# x*y = S^2 + b1/(a1 - a2)

### Sum =>
# (a1+a2)*(x^3 + y^3) + b1*(x+y) = 2*R
# (a1+a2)*(S^3 - 3*x*y*S) + b1*S - 2*R = 0
# (a1+a2)*(S^3 - 3*S*(S^2 + b1/(a1 - a2))) + b1*S - 2*R
# (a1+a2)*(-2*S^3 - 3*S*b1/(a1 - a2)) + b1*S - 2*R
# 2*(a1+a2)*S^3 + 3*b1*(a1+a2)/(a1 - a2)*S - b1*S + 2*R


solve.htm1 = function(b, a, R) {
	x.sum = roots(c(2*(a[1]+a[2]), 0, 3*b[1]*(a[1]+a[2])/(a[1] - a[2]) - b[1], 2*R))
	xy = x.sum^2 + b[1]/(a[1] - a[2])
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}

### Example: has Fractions
b = 1
a = c(1/2, 1/3)
R = 1
#
sol = solve.htm1(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
a[1]*x^3 + a[2]*y^3 + b[1]*x
a[2]*x^3 + a[1]*y^3 + b[1]*y

### Classic Polynomial

### TODO


### Example 2:
b = 3
a = c(1/2, -1/4)
R = 1
#
sol = solve.htm1(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
a[1]*x^3 + a[2]*y^3 + b[1]*x
a[2]*x^3 + a[1]*y^3 + b[1]*y

err = 80 - 48*x + 48*x^2 - 8*x^3 + 12*x^4 + x^6
round0(err)


################
################

### High-Power Terms: > 1

### a1*x^3 + a2*y^3 + b1*x*y

# a1*x^3 + a2*y^3 + b1*x*y = R
# a2*x^3 + a1*y^3 + b1*x*y = R

### Solution

### Diff =>
# (a1 - a2)*(x^3 - y^3) = 0
# (x - y)*(x^2 + y^2 + x*y) = 0
# Case x != y:
# S^2 - x*y = 0
# x*y = S^2

### Sum =>
# (a1+a2)*(x^3 + y^3) + 2*b1*x*y = 2*R
# (a1+a2)*(S^3 - 3*x*y*S) + 2*b1*x*y - 2*R = 0
# (a1+a2)*(S^3 - 3*S^3) + 2*b1*S^2 - 2*R = 0
# -2*(a1+a2)*S^3 + 2*b1*S^2 - 2*R
# (a1+a2)*S^3 - b1*S^2 + R

solve.htm = function(b, a, R) {
	x.sum = roots(c((a[1]+a[2]), - b[1], 0, R))
	xy = x.sum^2
	x.diff= sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}


### Example 1:
b = 3
a = c(1, 2)
R = 1
#
sol = solve.htm(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol


### Test
a[1]*x^3 + a[2]*y^3 + b[1]*x*y
a[2]*x^3 + a[1]*y^3 + b[1]*x*y

### Classic Polynomial
err = 1/9 + 1/3*x^2 - 2/3*x^3 + x^4 - x^5 + x^6
round0(err)

### TODO


### Example 2:
b = 3
a = c(1/2, -1)
R = 1
#
sol = solve.htm(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Classic Polynomial
err = 4 + 12*x^2 + 4*x^3 + 36*x^4 + 6*x^5 + x^6
round0(err)


###################################
###################################

###############
### Order 4 ###
###############

### x^4 + b*y

# x^4 + b1*y = R
# y^4 + b1*x = R

# Diff =>
# x^4 - y^4 - b1*(x-y) = 0
# (x - y)*((x+y)*(x^2+y^2) - b1) = 0
# => x = y *OR* (x+y)*(x^2+y^2) - b1 = 0;
# =>
# (x+y)^2 - 2*x*y = b1 / s
# x*y = (s^2 - b1/s)/2

### Sum =>
# x^4 + y^4 + b1*s - 2*R = 0
# s^4 - 4*x*y*s^2 + 2*(x*y)^2 + b1*s - 2*R
# s^4 - 2*(s^4 - b1*s) + (s^2 - b1/s)^2/2 + b1*s - 2*R
# -s^4 + 2*b1*s + (s^2 - b1/s)^2/2 + b1*s - 2*R
# s^4 - 3*b1*s - (s^2 - b1/s)^2/2 + 2*R
# 2*s^6 - 6*b1*s^3 - (s^3 - b1)^2 + 4*R*s^2
# s^6 - 4*b1*s^3 + 4*R*s^2 - b1^2


### Example
b = 2
R = 1
#
r.sum = roots(c(1,0,0, -4*b[1], 4*R, 0, - b[1]^2))
xy = (r.sum^2 - b[1]/r.sum)/2
r.diff = sqrt(r.sum^2 - 4*xy + 0i)
x = (r.sum + r.diff)/2
y = (r.sum - r.diff)/2
sol = cbind(x, y) # TODO: add also x = y cases;
sol

### Test
x^4 + b[1]*y
y^4 + b[1]*x

### Classical
err = x^12 - b[1]*x^9 - 3*R*x^8 + b[1]^2*x^6 + 2*R*b[1]*x^5 + 3*R^2*x^4 - b[1]^3*x^3 - R*b[1]^2*x^2 - R^2*b[1]*x + b[1]^4 - R^3
round0(err)

### Derivation:
# b1*y = R - x^4
# (R - x^4)^4/b1^4 + b1*x - R = 0
(R - x^4)^4 + b[1]^5*x - R*b[1]^4

x^16 - 4*R*x^12 + 6*R^2*x^8 - 4*R^3*x^4 + b[1]^5*x - R*b[1]^4 + R^4
(x^4 + b[1]*x - R)*(x^12 - b[1]*x^9 - 3*R*x^8 + b[1]^2*x^6 + 2*R*b[1]*x^5 + 3*R^2*x^4 - b[1]^3*x^3 - R*b[1]^2*x^2 - R^2*b[1]*x + b[1]^4 - R^3)



##############

#############
### Shift ###

### (x - s)^4 + b*y

# (x-s)^4 + b1*y = R
# (y-s)^4 + b1*x = R

### Solution:

# Trivial solution: x = y;

# Diff =>
# (x-s)^4 - (y-s)^4 - b1*(x-y) = 0
# (x - y)*((x-s)^3 + (y-s)^3 + (x-s)^2*(y-s) + (x-s)*(y-s)^2 - b1) = 0
# (x - y)*((x+y-2*s)*((x-s)^2 + (y-s)^2 - (x-s)*(y-s)) + (x+y-2*s)*(x-s)*(y-s) - b1) = 0
# (x - y)*((x+y-2*s)*(x^2 + y^2 - 2*s*(x+y) + 2*s^2 - (x-s)*(y-s) + (x-s)*(y-s)) - b1) = 0
# (x - y)*((x+y-2*s)*(x^2 + y^2 - 2*s*(x+y) + 2*s^2) - b1) = 0
# (x - y)*((x+y-2*s)*((x+y)^2 - 2*s*(x+y) - 2*x*y + 2*s^2) - b1) = 0
# (x - y)*((x+y)*((x+y)^2 - 2*s*(x+y) - 2*x*y + 2*s^2) - 2*s*(x+y)*((x+y) - 2*s) + 4*s*x*y - 4*s^3 - b1) = 0
# (x - y)*((x+y)*((x+y)^2 - 4*s*(x+y) - 2*x*y + 2*s^2 + 4*s^2) + 4*s*x*y - 4*s^3 - b1) = 0
# => x = y *OR* (x+y)*((x+y)^2 - 4*s*(x+y) - 2*x*y + 6*s^2) + 4*s*x*y - 4*s^3 - b1 = 0;
# =>
# (x+y)*((x+y)^2 - 4*s*(x+y) - 2*x*y + 6*s^2) + 4*s*x*y - 4*s^3 - b1 = 0
# Z*(Z^2 - 4*s*Z - 2*x*y + 6*s^2) + 4*s*x*y = 4*s^3 + b1
# Z*(Z^2 - 4*s*Z + 6*s^2) = (2*Z - 4*s)*x*y + 4*s^3 + b1
# x*y = (Z*(Z^2 - 4*s*Z + 6*s^2) - 4*s^3 - b1) / (2*Z - 4*s)

### Sum =>
# (x-s)^4 + (y-s)^4 + b1*(x+y) = 2*R
# x^4 + y^4 - 4*s*(x^3 + y^3) + 6*s^2*(x^2 + y^2) - 4*s^3*(x+y) + b1*(x+y) + 2*s^4 - 2*R = 0;
# x^4 + y^4 - 4*s*Z*(Z^2 - 3*x*y) + 6*s^2*(Z^2 - 2*x*y) - (4*s^3 - b1)*Z + 2*s^4 - 2*R = 0;
# x^4 + y^4 - 4*s*Z^3 + 6*s^2*Z^2 - (4*s^3 - b1)*Z + 12*(s*Z - s^2)*x*y + 2*s^4 - 2*R = 0;
# Z^4 - 4*x*y*Z^2 + 2*(x*y)^2 - 4*s*Z^3 + 6*s^2*Z^2 - (4*s^3 - b1)*Z + 12*(s*Z - s^2)*x*y + 2*s^4 - 2*R = 0;
# Z^4 - 4*s*Z^3 + 6*s^2*Z^2 - (4*s^3 - b1)*Z + 2*(x*y)^2 - 4*(Z^2 - 3*s*Z + 3*s^2)*x*y + 2*s^4 - 2*R = 0;
# =>
# Z^6 - 12*s*Z^5 + 60*s^2*Z^4 - 4*(b1 + 40*s^3)*Z^3 + 4*(R + 5*b1*s + 60*s^4)*Z^2 - 16*(R*s + 2*b1*s^2 + 12*s^5)*Z + 16*(R*s^2 + b1*s^3 + 4*s^6) - b1^2


### Example
b = 3
R = 1
s = 1
#
coeff = c(1, -12*s, 60*s^2, -4*(b[1] + 40*s^3), 4*(R + 5*b[1]*s + 60*s^4), -16*(R*s + 2*b[1]*s^2 + 12*s^5), 16*(R*s^2 + b[1]*s^3 + 4*s^6) - b[1]^2)
x.sum = roots(coeff)
xy = (x.sum*(x.sum^2 - 4*s*x.sum + 6*s^2) - 4*s^3 - b[1]) / (2*x.sum - 4*s)
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
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
+ 3*(b[1]*s - R)^2*x^4 - b[1]^3*x^3 + b[1]^2*(b[1]*s - R)*x^2 - b[1]*(b[1]*s - R)^2*x + (b[1]*s - R)^3 + b[1]^4)
round0(err)

### Derivation:
# b1*y = (-(x-s)^4 + R)
# (b1*y - b1*s)^4 + b1^5*x = b1^4*R
# =>
(-(x-s)^4 + R - b[1]*s)^4 + b[1]^5*x - b[1]^4*R
((x-s)^4 + b[1]*s - R)^4 + b[1]^5*x - b[1]^4*R
# x = x - s # shift-back;
(x^4 + b[1]*s - R)^4 + b[1]^5*(x + s) - b[1]^4*R
x^16 + 4*(b[1]*s - R)*x^12 + 6*(b[1]*s - R)^2*x^8 + 4*(b[1]*s - R)^3*x^4 + (b[1]^5)*x + (b[1]*s - R)^4 + b[1]^5*s - R*b[1]^4
(x^4 + b[1]*x + b[1]*s - R)*(x^12 - b[1]*x^9 - 3*(R - b[1]*s)*x^8 + b[1]^2*x^6 - 2*b[1]*(b[1]*s - R)*x^5 +
+ 3*(b[1]*s - R)^2*x^4 - b[1]^3*x^3 + b[1]^2*(b[1]*s - R)*x^2 - b[1]*(b[1]*s - R)^2*x + (b[1]*s - R)^3 + b[1]^4)



###################################

### x^4 + b*x*y

# x^4 + b1*x*y = R
# y^4 + b1*x*y = R

### Solution:

# Trivial solutions: x = y *OR* x = -y;

# Diff =>
# x^4 - y^4 = 0
# (x - y)*(x^3 + y^3 + x*y*(x+y)) = 0
# (x - y)*((x+y)^3 - 2*x*y*(x+y)) = 0
# (x - y)*(x + y)((x+y)^2 - 2*x*y) = 0
# => x = y *OR* x = -y *OR* 2*x*y = (x+y)^2;
# Case: x != +/- y =>
# 2*x*y = Z^2

### Sum =>
# x^4 + y^4 + 2*b1*x*y - 2*R = 0
# (x+y)^4 - 4*x*y*(x+y)^2 + 2*(x*y)^2 + 2*b1*x*y - 2*R = 0
# Z^4 - 2*Z^2*Z^2 + 2/4*Z^4 + b1*Z^2 - 2*R
# -1/2*Z^4 + b1*Z^2 - 2*R
# Z^4 - 2*b1*Z^2 + 4*R


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


#############################
### x^4 + b2*(x*y)^2 + b1*x*y

# x^4 + b2*(x*y)^2 + b1*x*y = R
# y^4 + b2*(x*y)^2 + b1*x*y = R

### Solution:

# Trivial solutions: x = y *OR* x = -y;

# Diff =>
# x^4 - y^4 = 0
# (x - y)*(x^3 + y^3 + x*y*(x+y)) = 0
# (x - y)*((x+y)^3 - 2*x*y*(x+y)) = 0
# (x - y)*(x + y)((x+y)^2 - 2*x*y) = 0
# => x = y *OR* x = -y *OR* 2*x*y = (x+y)^2;
# Case: x != +/- y =>
# 2*x*y = Z^2

### Sum =>
# x^4 + y^4 + 2*b2*(x*y)^2 + 2*b1*x*y - 2*R = 0
# (x+y)^4 - 4*x*y*(x+y)^2 + 2*(x*y)^2 + 2*b2*(x*y)^2 + 2*b1*x*y - 2*R = 0
# Z^4 - 2*Z^4 + 2/4*Z^4 + 2/4*b2*Z^4 + b1*Z^2 - 2*R
# (b2 - 1)/2*Z^4 + b1*Z^2 - 2*R
# (b2 - 1)*Z^4 + 2*b1*Z^2 - 4*R

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

##########################
### Mixed Highest Term ###
##########################

#########################
### x^j*y^k + P(x, y) ###


###################
### x^2*y + b*x ###

# x^2*y + b1*x = R
# y^2*x + b1*y = R

### Solution:

### Diff =>
# x*y*(x - y) + b1*(x-y) = 0
# (x - y)*(x*y + b1) = 0
# Case: x != y
# x*y = -b1

### Sum =>
# x*y*(x+y) + b1*(x+y) = 2*R
# (x+y)*(x*y + b1) - 2*R = 0
# but: x*y + b1 = 0
# => NO solution (x != y);



###############
### x^2*y + b*y

# x^2*y + b1*y = R
# y^2*x + b1*x = R

# very simple system

### Solution:

### Diff =>
# x*y*(x - y) - b1*(x-y) = 0
# (x - y)*(x*y - b1) = 0
# Case: x != y
# x*y = b1

### Sum =>
# x*y*(x+y) + b1*(x+y) = 2*R
# (x+y)*(x*y + b1) - 2*R = 0
# 2*b1*Z - 2*R = 0
# b1*Z = R
# Z = R / b1


### Example
b = 3
R = 1
#
x.sum = R/b[1]
xy = b[1]
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol


### Test
x^2*y + b[1]*y
y^2*x + b[1]*x


#########################
### x^2*y + b2*x^2 + b1*x

# x^2*y + b2*x^2 + b1*x = R
# y^2*x + b2*y^2 + b1*y = R

### Solution:

### Diff =>
# x*y*(x - y) + b2*(x^2 - y^2) + b1*(x-y) = 0
# (x - y)*(x*y + b2*(x+y) + b1) = 0
# Case: x != y
# x*y = -b2*(x+y) - b1
# x*y = -b2*Z - b1

### Sum =>
# x*y*(x+y) + b2*(x^2 + y^2) + b1*(x+y) = 2*R
# (x+y)*(x*y + b1) + b2*(x^2 + y^2) - 2*R = 0
# Z*(-b2*Z - b1 + b1) + b2*(Z^2 - 2*x*y) - 2*R
# -b2*Z^2 + b2*Z^2 - 2*b2*(-b2*Z - b1) - 2*R
# 2*b2*(b2*Z + b1) - 2*R
# b2^2*Z + b1*b2 - R
# Z = (R - b1*b2) / b2^2

### Example:
b = c(1, 3)
R = -1
#
x.sum = (R - b[1]*b[2]) / b[2]^2
xy = -b[2]*x.sum - b[1]
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol

### Test
x^2*y + b[2]*x^2 + b[1]*x
y^2*x + b[2]*y^2 + b[1]*y

### Classical Polynomial:
# trivial P2;
b[2]*x^2 - (R/b[2] - b[1])*x - R


###########################

###########################
### x^2*y^2 + b2*x^2 + b1*x

# x^2*y^2 + b2*x^2 + b1*x = R
# y^2*x^2 + b2*y^2 + b1*y = R

### Solution:

### Diff =>
# b2*(x^2 - y^2) + b1*(x-y) = 0
# (x - y)*(b2*(x+y) + b1) = 0
# Case: x != y
# x + y = -b1/b2
# Z = -b1/b2

### Sum =>
# 2*(x*y)^2 + b2*(x^2 + y^2) + b1*(x+y) = 2*R
# 2*(x*y)^2 + b2*(Z^2 - 2*x*y) + b1*Z - 2*R = 0
# 2*(x*y)^2 - 2*b2*x*y + b2*Z^2 + b1*Z - 2*R = 0
# (x*y)^2 - b2*x*y - R = 0


### Example:
b = c(3, 1)
R = 1
#
x.sum = - b[1]/b[2]
xy = roots(c(1, - b[2], - R))
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol


### Test
x^2*y^2 + b[2]*x^2 + b[1]*x
y^2*x^2 + b[2]*y^2 + b[1]*y

### Classical Polynomial
round0.p(poly.calc(sol[,1]))

- R + b[1]*x + (b[2] + b[1]^2/b[2]^2)*x^2 + 2*b[1]/b[2]*x^3 + x^4


######################################


######################################
### x^2*y^2 + b3*x^2*y + b2*x^2 + b1*x

# x^2*y^2 + b3*x^2*y + b2*x^2 + b1*x = R
# y^2*x^2 + b3*y^2*x + b2*y^2 + b1*y = R

### Solution:

### Diff =>
# b3*x*y*(x-y) + b2*(x^2 - y^2) + b1*(x-y) = 0
# (x - y)*(b3*x*y + b2*(x+y) + b1) = 0
# Case: x != y
# b3*x*y = - (b2*Z + b1)
# Z = -(b3*x*y + b1) / b2

### Sum =>
# 2*(x*y)^2 + b3*x*y*(x+y) + b2*(x^2 + y^2) + b1*(x+y) = 2*R
# 2*(x*y)^2 - (b2*Z + b1)*Z + b2*(Z^2 - 2*x*y) + b1*Z - 2*R = 0
# 2*(x*y)^2 - 2*b2*x*y - 2*R
# (x*y)^2 - b2*x*y - R

### Example:
b = c(3, 1, 1)
R = 1
#
xy = roots(c(1, -b[2], -R))
x.sum = -(b[3] * xy + b[1]) / b[2]
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol

### Test
x^2*y^2 + b[3]*x^2*y + b[2]*x^2 + b[1]*x
y^2*x^2 + b[3]*y^2*x + b[2]*y^2 + b[1]*y

### Classical Polynomial
round0.p(poly.calc(sol[,1]))

### TODO: classical Polynomial: P4;


#################################


#################################
### x^2*y^2 + b4*x^2*y + b3*x*y^2
### + b2*x^2 + b1*x

# x^2*y^2 + b4*x^2*y + b3*x*y^2 + b2*x^2 + b1*x = R
# y^2*x^2 + b4*y^2*x + b3*y*x^2 + b2*y^2 + b1*y = R

### Solution:

### Diff =>
# (b4 - b3)*x*y*(x-y) + b2*(x^2 - y^2) + b1*(x-y) = 0
# (x - y)*((b4 - b3)*x*y + b2*(x+y) + b1) = 0
# Case: x != y

# Sub-Case 1: b4 = b3
# x + y = - b1/b2;
# Z = - b1/b2;
#
# Sub-Case 2: b4 != b3
# (b4 - b3)*x*y = - (b2*Z + b1)
# Z = -((b4 - b3)*x*y + b1) / b2

### Sum =>
# 2*(x*y)^2 + (b3+b4)*x*y*(x+y) + b2*(x^2 + y^2) + b1*(x+y) = 2*R
#  2*(x*y)^2 + (b3+b4)*x*y*Z + b2*(Z^2 - 2*x*y) + b1*Z - 2*R = 0
# Sub-Case 1:
#  2*(x*y)^2 - b1/b2*(b3+b4)*x*y + b2*(b1^2/b2^2 - 2*x*y) - b1*b1/b2 - 2*R = 0
#  2*(x*y)^2 - b1/b2*(b3+b4)*x*y - 2*b2*x*y - 2*R = 0
#  2*b2*(x*y)^2 - (b1*(b3+b4) + 2*b2^2)*x*y - 2*b2*R = 0
# Sub-Case 2:
# ...

### TODO: Sub-Case 2;

solve.htMixt22 = function(b, R) {
	if(b[3] == b[4]) {
		x.sum = - b[1]/b[2];
		xy = roots(c(2*b[2], - (b[1]*(b[3]+b[4]) + 2*b[2]^2), - 2*b[2]*R))
	}
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	sol
}

### Example SC-1:
b = c(3,2,1,1)
R = 1
#
sol = solve.htMixt22(b, R)
x = sol[,1]; y = sol[,2];
sol

### Test
x^2*y^2 + b[4]*x^2*y + b[3]*x*y^2 + b[2]*x^2 + b[1]*x
y^2*x^2 + b[4]*y^2*x + b[3]*y*x^2 + b[2]*y^2 + b[1]*y

### Classical Polynomial
round0.p(poly.calc(sol[,1]))


### Example SC-2:
### TODO


##########################

######################
### (x*y)^3 Series ###
######################

### (x*y)^3 + b*x^3

# (x*y)^3 + b1*x^3 = R
# (x*y)^3 + b1*y^3 = R


m3 = unty(3, all=F)

### Solution:

### Diff =>
# x^3 - y^3 = 0
# y = x * m3

### =>
# x^6 + b1*x^3 = R

### Example
b = 1
R = 1
#
x = roots(c(1, b[1], - R))
x = as.vector(sapply(x, function(r) roots(c(1,0,0, -r))))
y = x*m3
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol

### Test
(x*y)^3 + b[1]*x^3
(x*y)^3 + b[1]*y^3

### Classic Polynomial:
round0.p(poly.calc(sol[,1]))
# P12 based on P4;


##########################

###############
### x^3*y + b*x

# x^3*y + b1*x = R
# y^3*x + b1*y = R

### Solution:

### Diff =>
# x*y*(x^2 - y^2) + b1*(x-y) = 0
# (x - y)*(x*y*(x+y) + b1) = 0
# Case: x != y
# x*y = -b1 / (x+y)
# x*y = -b1/Z

### Sum =>
# x*y*(x^2 + y^2) + b1*(x+y) = 2*R
# -b1/Z * (Z^2 - 2*x*y) + b1*Z - 2*R = 0
# -b1*(Z^2 - 2*x*y) + b1*Z^2 - 2*R*Z = 0
# 2*b1*x*y - 2*R*Z = 0
# R*Z - b1*x*y = 0
# R*Z + b1^2/Z = 0
# R*Z^2 + b1^2 = 0
# Z^2 = -b1^2 / R;

### Example:
b = 3
R = 1
#
x.sum = sqrt(-b[1]^2 / R + 0i)
x.sum = c(x.sum, -x.sum)
xy = -b[1]/x.sum
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol


### Test
x^3*y + b[1]*x
y^3*x + b[1]*y


#####################

###############
### Order 5 ###
###############

###############
### x^4*y + b*x

# x^4*y + b1*x = R
# y^4*x + b1*y = R

### Solution:

### Diff =>
# x*y*(x^3 - y^3) + b1*(x-y) = 0
# (x - y)*(x*y*(x^2 + y^2 + x*y) + b1) = 0
# (x - y)*(x*y*(x+y)^2 - (x*y)^2 + b1) = 0
# Case: x != y
# x*y*(x+y)^2 - (x*y)^2 + b1 = 0
# (x*y)^2 - x*y*Z^2 = b1

### Sum =>
# x*y*(x^3 + y^3) + b1*(x+y) = 2*R
# (x+y)*x*y*(x^2 + y^2 - x*y) + b1*(x+y) - 2*R = 0
# Z*x*y*(Z^2 - 3*x*y) + b1*Z - 2*R
# Z^3*x*y - 3*Z*(x*y)^2 + b1*Z - 2*R
# Z^3*x*y - 3*Z*(x*y*Z^2 + b1) + b1*Z - 2*R
# -2*Z^3*x*y - 2*b1*Z - 2*R
# Z^3*x*y + b1*Z + R
# x*y = -(b1*Z + R) / Z^3
# =>
# (b1*Z + R)^2 / Z^6 + (b1*Z + R) / Z^3 * Z^2 - b1 = 0
# (b1*Z + R)^2 + (b1*Z + R)*Z^5 - b1*Z^6 = 0
# (b1*Z + R)^2 + R*Z^5
# R*Z^5 + b1^2*Z^2 + 2*b1*R*Z + R^2

### Example:
b = 3
R = 1
#
x.sum = roots(c(R, 0, 0, b[1]^2, 2*b[1]*R, R^2))
xy = -(b[1]*x.sum + R) / x.sum^3
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol

### Test
x^4*y + b[1]*x
y^4*x + b[1]*y

### Classical Polynomial: P10
round0.p(poly.calc(sol[,1]))

err = R*x^10 + b[1]^2*x^7 - 2*R*b[1]*x^6 + R^2*x^5 - b[1]^3*x^3 + 3*R*b[1]^2*x^2 - 3*R^2*b[1]*x + R^3
round0(err)


### Derivation:
# y = (R - b[1]*x)/x^4
(R - b[1]*x)^4/x^16 * x + b[1]*(R - b[1]*x)/x^4 - R # = 0
(R - b[1]*x)^4/x^15 + b[1]*(R - b[1]*x)/x^4 - R
(R - b[1]*x)^4 + b[1]*(R - b[1]*x)*x^11 - R*x^15
R*x^15 - b[1]*(R - b[1]*x)*x^11 - (R - b[1]*x)^4
R*x^15 + b[1]^2*x^12 - b[1]*R*x^11 - (R - b[1]*x)^4

R*x^15 + b[1]^2*x^12 - R*b[1]*x^11 - b[1]^4*x^4 + 4*R*b[1]^3*x^3 - 6*R^2*b[1]^2*x^2 + 4*R^3*b[1]*x^1 - R^4
(x^5 + b[1]*x - R)*(R*x^10 + b[1]^2*x^7 - 2*R*b[1]*x^6 + R^2*x^5 - b[1]^3*x^3 + 3*R*b[1]^2*x^2 - 3*R^2*b[1]*x + R^3)


###################################
###################################

###############
### Order 5 ###
###############

#####################
### x^5 + b*(x+y) = 0

# x^5 + b1*(x+y) = 0
# y^5 + b1*(y+x) = 0

### Solution

m5 = unity(5, all=F)
m4 = unity(4, all=T)

### Diff =>
# x^5 - y^5 = 0
# y = x * m5^j

### =>
# x^5 + b1*(m5+1)*x = 0
# x^4 + b1*(m5+1) = 0

### Example
b = 1
x = as.vector(sapply(1:4, function(j) ( - b[1]*(m5^j + 1))^(1/4) * m4))
y = - (x^5 / b[1] + x)
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol

### Test
round0(x^5 + b[1]*(x+y))
round0(y^5 + b[1]*(y+x))

### Classic Polynomial
round0.p(poly.calc(sol[,1]))

err = 1 + 4*x^4 + 12*x^8 + 22*x^12 + 30*x^16 + 28*x^20 + 17*x^24 + 6*x^28 + x^32
round0(err)
# based on: (1 + 2*x + 4*x^2 + 3*x^3 + x^4)^2



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



###################################


###################################
###################################
###################################


### Heterogenous Symmetric
### Polynomial Systems: 3 Variables

### 3 Variables:
# x^n + P(x, y, z) = R
# y^n + P(y, z, x) = R
# z^n + P(z, x, y) = R



###############

###############
### Order 2 ###

### x[i]^2 + b*x[i+1]

# x^2 + b1*y = R
# y^2 + b1*z = R
# z^2 + b1*x = R

# Trivial solution: x = y = z;

### TODO:
# - find correct solution;
# - current solution involves P6,
#   which is also the classic polynomial;

### Method 1:
# Diff =>
# x^2 - y^2 = b1*(z-y)
# y^2 - z^2 = b1*(x-z)
# z^2 - x^2 = b1*(y-x)
# Prod =>
# (x+y)*(x+z)*(y+z) = (-1)*b1^3;
# S^3 - (x^3 + y^3 + z^3) + 3*b1^3 = 0;

# (x^2 + x*y + x*z + y*z)*(y+z) = - b1^3
# S[x^2y] + 2*x*y*z + b1^3 = 0;
# E2*S - 3*E3 + 2*E3 + b1^3 = 0;
# E2*S - E3 + b1^3 = 0

# Sum(x[i]*...) =>
# x^3 + y^3 + z^3 + b1*E2 = R*S
# =>
# S^3 - (R*S - b1*E2) + 3*b1^3 = 0
# S^3 - R*S + b1*E2 + 3*b1^3 = 0

# Sum =>
# x^2 + y^2 + z^2 + b1*(x+y+z) = 3*R
# S^2 - 2*E2 + b1*S - 3*R = 0;
# 2*E2 = S^2 + b1*S - 3*R;
# =>
# S^3 - R*S + b1*(S^2 + b1*S - 3*R)/2 + 3*b1^3 = 0
# 2*S^3 - 2*R*S + b1*(S^2 + b1*S - 3*R) + 6*b1^3 = 0
# 2*S^3 + b1*S^2 + (b1^2 - 2*R)*S - 3*b1*R + 6*b1^3 = 0



### Alternative
# b1*y = R - x^2
# b1^3*z = b1^2*R - (R - x^2)^2
# b1^3*z = b1^2*R - R^2 - x^4 + 2*R*x^2
# =>
# x^8 - 4*R*x^6 + (6*R^2 - 2*b1^2*R)*x^4 + 4*R^2*(b1^2 - R)*x^2 + b[1]^7*x + (b1^2*R - R^2)^2 - b[1]^6*R = 0
# (x^2 + b1*x - R) * P6;

### Example
b = 2
R = 1
#
x.sum = roots(c(2, b[1], (b[1]^2 - 2*R), - 3*b[1]*R + 6*b[1]^3))
E2 = (x.sum^2 + b[1]*x.sum - 3*R)/2
E3 = E2*x.sum + b1^3
x = sapply(1:length(x.sum), function(id) roots(c(1, -x.sum[id], E2[id], -E3[id])))
x = cbind(as.vector(x[,-1])) # TODO: remove robustly set of wrong solutions
y = (R - x^2)/b[1]
z = (R - y^2)/b[1]
sol = cbind(x,y,z)
sol

### Test
x^2 + b[1]*y 
y^2 + b[1]*z
z^2 + b[1]*x

round0.p(poly.calc(sol[,1]))


# alternative / classic
x = roots(c(1,0, - 4*R,0, (6*R^2 - 2*b[1]^2*R), 0, 4*R^2*(b[1]^2 - R), b[1]^7, (b[1]^2*R - R^2)^2 - b[1]^6*R))
y = (R - x^2)/b[1]
z = (R - y^2)/b[1]
sol = cbind(x, y, z)
sol

### Test
x^2 + b[1]*y 
y^2 + b[1]*z
z^2 + b[1]*x


######################
######################

### x[i]^2 + b*x[-i]

# x^2 + b1*y*z = R
# y^2 + b1*x*z = R
# z^2 + b1*x*y = R

### Special Case: b1 = 2
# Z^2 = 3*R

# x^2 - y^2 = b1*z*(x-y)
# x^2 - z^2 = b1*y*(x-z)
# y^2 - z^2 = b1*x*(y-z)
# => if x != y != z
#   x + y - b1*z = 0
#   x - b1*y + z = 0
# - b1*x + y + z = 0
# => x = y = z = 0;

# Case 1: x = y
# x^2 + b1*x*z = R
# z^2 + b1*x^2 = R
# =>
# b1*z = R/x - x
# b1^2*z^2 + b1^3*x^2 = b1^2*R
# x^2 - 2*R + R^2/x^2 + b1^3*x^2 - b1^2*R = 0
# (b1^3+1)*x^4 - R*(b1^2 + 2)*x^2 + R^2 = 0

b = 1
R = 2
#
x = roots(c((b[1]^3+1), 0, - R*(b[1]^2 + 2), 0, R^2))
y = x
z = (R - x^2)/y/b[1]
sol = round0(cbind(x, y, z))
sol

### Test
x^2 + b[1]*y*z
y^2 + b[1]*x*z
z^2 + b[1]*x*y


########################
########################

### Variant:
### x[i]^2 + b*(x[j] + x[k])

# x^2 + b1*(y+z) = R
# y^2 + b1*(x+z) = R
# z^2 + b1*(x+y) = R
# [a trivial system]

# Trivial solution: x = y = z;

### Solution

# Diff =>
# x^2 - y^2 = b1*(x-y)
# x^2 - z^2 = b1*(x-z)
# y^2 - z^2 = b1*(y-z)

# if x != y != z
# x + y = b1
# x + z = b1
# y + z = b1
# => x = y = z, which violates assumption;

# =>
# x = y
# x^2 + b1*(x+z) = R
# z^2 + 2*b1*x = R

# Case x != z
# x + z = b1
# x^2 + b1^2 = R

### Example:
b = 3
R = 1
#
x = sqrt(R - b[1]^2 + 0i)
x = c(x, -x)
y = x
z = b[1] - x
sol = cbind(x, y, z)
sol

### Test
x^2 + b[1]*(y+z)
y^2 + b[1]*(x+z)
z^2 + b[1]*(x+y)

####################

### Shifted Variant:
### (x[i] - s)^2 + b*(x[j] + x[k])

# (x-s)^2 + b1*(y+z) = R
# (y-s)^2 + b1*(x+z) = R
# (z-s)^2 + b1*(x+y) = R


### Solution

### Diff =>
# (x-s)^2 - (y-s)^2 = b1*(x - y)
# (x - y)*(x + y - 2*s) = b1*(x - y)
# (x - y)*(x + y - 2*s - b1) = 0
# (x - z)*(x + z - 2*s - b1) = 0
# (y - z)*(y + z - 2*s - b1) = 0

### Case x != y != z
# is NOT solvable (except in very special conditions);

### Case x = y && x != z
# x + z - 2*s - b1 = 0
# x + z = b1 + 2*s;
# =>
# (x-s)^2 + b1*(x+z) - R = 0
# x^2 - 2*s*x + s^2 + b1*(b1 + 2*s) - R
# x^2 - 2*s*x + b1^2 + 2*b1*s + s^2 - R

### Example
b = 3
s = -1
R = 1
#
x = roots(c(1, - 2*s, b[1]^2 + 2*b[1]*s + s^2 - R))
y = x
z = b[1] + 2*s - x
sol = cbind(x, y, z)
sol


### Test
(x-s)^2 + b[1]*(y+z)
(y-s)^2 + b[1]*(x+z)
(z-s)^2 + b[1]*(x+y)


####################

### "Asymmetric" Variant:
### x[i]^2 + b1*x[j] + b2*x[k]

# x^2 + b1*y + b2*z = R
# y^2 + b1*z + b2*x = R
# z^2 + b1*x + b2*y = R

### Solution

# TODO

### Sum =>
# S^2 - 2*E2 + (b1+b2)*S = 3*R

### Diff =>
# x^2 - y^2 = b2*(x - z) - b1*(y - z)
# x^2 - z^2 = b2*(y - z) - b1*(y - x)
# y^2 - z^2 = b2*(y - x) - b1*(z - x)


#########################

#########################
#########################

### High-Power Terms: > 1

### x[i]^2 + x[j]^2 + b*(x[i] + x[j])

# x^2 + y^2 + b1*(x+y) = R
# y^2 + z^2 + b1*(y+z) = R
# x^2 + z^2 + b1*(x+z) = R

# trivial solution: x = y = z;

### Solution

### Diff =>
# x^2 - z^2 = -b1*(x - z)
# (x-z)*(x + z + b1) = 0

# Case: x != y != z
# - has NO solutions;

# Case x = y && x != z
# TODO: trivial;




########################
########################

### High-Power Terms: > 1

### x[i]^2 + x[j]^2 + b*x[k]

# x^2 + y^2 + b1*z = R
# y^2 + z^2 + b1*x = R
# x^2 + z^2 + b1*y = R

# trivial solution: x = y = z;
# trivial system;

### Solution

### Diff =>
# x^2 - z^2 = b1*(x - z)
# (x-z)*(x + z - b1) = 0

# Case: x != y != z
# - has NO solutions;

# Case x = y && x != z
# z = -x + b1;
# =>
# 2*x^2 + b1*z - R = 0
# 2*x^2 - b1*x + b1^2 - R

### Example

b = 3
R = 1
#
x = roots(c(2, - b[1], b[1]^2 - R))
y = x
z = -x + b[1]
sol = cbind(x, y, z)
sol

### Test
x^2 + y^2 + b[1]*z
y^2 + z^2 + b[1]*x
x^2 + z^2 + b[1]*y


########################
########################

### Problems:
# - Difference works well for systems with 2 variables;
# - but it does NOT work well in systems with 3 variables;

###########
### Order 3
### x[i]^3 + b*(x[j] + x[k])

# x^3 + b1*(y+z) = R
# y^3 + b1*(x+z) = R
# z^3 + b1*(x+y) = R

# Trivial solution: x = y = z;
# Trivial system;

### Solution

# Diff =>
# x^3 - y^3 = b1*(x-y)
# x^3 - z^3 = b1*(x-z)
# y^3 - z^3 = b1*(y-z)

# (x-y)*(x^2 + y^2 + x*y - b1) = 0

# Case: x != y != z
# x^2 + y^2 + x*y - b1 = 0
# x^2 + z^2 + x*z - b1 = 0
# y^2 + z^2 + y*z - b1 = 0
# E2 = 2*S^2/3 - b1
# Diff =>
# y^2 - z^2 + x*(y - z) = 0
# Case x != y != z:
# x + y + z = 0
# => E2 = -b1

### Sum =>
# x^3 + y^3 + z^3 + 2*b1*S = 3*R, where S = 0
# x^3 + y^3 + z^3 = 3*R

# x^3 + y^3 + z^3 = S^3 - 3*E2*S + 3*E3
# =>
# E3 = R


### Example 1:
b = 3
R = 1
#
x = roots(c(1,0, -b[1], -R))
y = as.vector(sapply(x, function(x) roots(c(1, x, x^2 - b[1]))))
x = rep(x, each=2)
z = -x-y
sol = cbind(x,y,z)
sol

### Test
x^3 + b[1]*(y+z)
y^3 + b[1]*(x+z)
z^3 + b[1]*(x+y)

# trivial:
(x^3 - b*x - R)^2


########################

###############
### Order 2 ###

### x[i]^2*x[j] + b*Sum

# x^2*y + b1*(x+y+z) = R
# y^2*z + b1*(x+y+z) = R
# z^2*x + b1*(x+y+z) = R


# Diff =>
# x^2 = y*z
# y^2 = x*z
# z^2 = x*y

# y = x^2/z
# => x^4/z^2 = x*z
# => x^3 = z^3
# => y^3 = z^3

# Case x != y != z
# y = x*m
# z = x*m^2
# => 1 + m + m^2 = 0!
# x^3*m = R

m3 = unity(3, all=FALSE)

### Example 1:

b = 3
R = 1
#
x = (R/m3)^(1/3) * c(1, m3, m3^2)
y = x * m3
z = x * m3^2
sol = cbind(x,y,z)
sol

### Test
x^2*y + b[1]*(x+y+z)
y^2*z + b[1]*(x+y+z)
z^2*x + b[1]*(x+y+z)

### Classical Polynomial
x^9 - R^3

#######################

### Shifted
### x[i]^2 * (x[j] - shift) + b*Sum

# x^2*(y - s) + b1*(x+y+z) = R
# y^2*(z - s) + b1*(x+y+z) = R
# z^2*(x - s) + b1*(x+y+z) = R

### Solution

# Diff =>
# x^2*(y - s) = y^2*(z - s)
# y^2*(z - s) = z^2*(x - s)
# z^2*(x - s) = x^2*(y - s)
# =>
# x^4*(y-s)^2 = y^2*z^2*(x-s)*(z-s)

### TODO



#####################

### Shifted
### x[i]^2 * (x[i] - shift) + b*Sum

# x^2*(x - s) + b1*(x+y+z) = R
# y^2*(y - s) + b1*(x+y+z) = R
# z^2*(z - s) + b1*(x+y+z) = R

### Solution

# - trivial solution: x = y = z;

### Diff =>
# x^2*(x-s) = y^2*(y-s)
# y^2*(y-s) = z^2*(z-s)
# z^2*(z-s) = x^2*(x-s)
# =>
# x^3 - y^3 - s*(x^2 - y^2) = 0
# (x-y)*(x^2 + y^2 + x*y - s*(x+y)) = 0

# Case: x != y != z
# x^2 + y^2 + x*y - s*(x+y) = 0
# x^2 + z^2 + x*z - s*(x+z) = 0
# y^2 + z^2 + y*z - s*(y+z) = 0
# =>
# y^2 + x*y - s*(x+y) = z^2 + x*z - s*(x+z)
# y^2 - z^2 + x*(y-z) - s*(y-z) = 0
# (y-z)*(x + y + z - s) = 0
# x + y + z = s

### Sum =>
# x^3 + y^3 + z^3 - s*(x^2 + y^2 + z^2) + 3*b1*S = 3*R
# x^3 + y^3 + z^3 = S^3 - 3*E2*S + 3*E3 =>
# S^3 - 3*E2*S + 3*E3 - s*(S^2 - 2*E2) + 3*b1*S - 3*R = 0
# S = s =>
# - s*E2 + 3*b1*s - 3*R + 3*E3 = 0
# E3 =  s*E2/3 - b1*s + R

# E2 = 2/3*S^2 - 2/3*s*S
# E2 = 0
# E3 = s*E2/3 - b*s + R
# E3 = - b1*s + R

### Example:

b = 3
R = 1
s = 2
#
x = roots(c(1, -s, 0, b[1]*s - R))
yz.sum = s - x
yz = yz.sum^2 - s*yz.sum
yz.diff = sqrt(yz.sum^2 - 4*yz + 0i)
y = (yz.sum + yz.diff)/2
z = (yz.sum - yz.diff)/2
sol = cbind(x,y,z)
sol

### Test
x^2*(x - s) + b[1]*(x+y+z)
y^2*(y - s) + b[1]*(x+y+z)
z^2*(z - s) + b[1]*(x+y+z)

### Classical Polynomial
# x => trivial polynomial: P3;
# y, z => (P3)^2

round0.p(poly.calc(sol[,2:3]))
round0.p(poly.calc(sol[,2:3] - s))

# (x^3 - 2*x^2 + 5)^2
err = 25 - 20*x^2 + 10*x^3 + 4*x^4 - 4*x^5 + x^6
round0(err)
x = x - s # shift back;
err = 25 + 40*x + 56*x^2 + 42*x^3 + 24*x^4 + 8*x^5 + x^6
round0(err)



########################

###############
### Order 3 ###

### x^3 + b3*x*y*z + b[j]*Sum[j]

# x^3 + b3*x*y*z + b2*(x^2+y^2+z^2) + b1*(x+y+z) = R
# y^3 + b3*x*y*z + b2*(x^2+y^2+z^2) + b1*(x+y+z) = R
# z^3 + b3*x*y*z + b2*(x^2+y^2+z^2) + b1*(x+y+z) = R

m3 = unity(3, all=F)

### Solution:

# Trivial system;

### Diff =>
# x^3 - y^3 = 0
# x^3 - z^3 = 0
# Case: x != y != z =>
# y = x*m3
# z = x*m3^2

### =>
# x^3 + b3*x^3 = R
# x^3 = R / (b3 + 1)

### TODO


########################

### x^3 + b2*x^2*y*z + b1*x

# x^3 + b2*x^2*y*z + b1*x = R
# y^3 + b2*x*y^2*z + b1*y = R
# z^3 + b2*x*y*z^2 + b1*z = R

### Solution:

### Diff =>
# x^3 - y^3 + b2*x*y*z*(x-y) + b1*(x-y) = 0
# (x-y)*(x^2 + y^2 + x*y + b2*x*y*z + b1) = 0
# (y-z)*(y^2 + z^2 + y*z + b2*x*y*z + b1) = 0
# (z-x)*(x^2 + z^2 + x*z + b2*x*y*z + b1) = 0
# Case: x != y != z =>
# x^2 + y^2 + x*y + b2*x*y*z + b1 = 0
# y^2 + z^2 + y*z + b2*x*y*z + b1 = 0
# x^2 + z^2 + x*z + b2*x*y*z + b1 = 0
# =>
# x + y + z = 0

### Sum =>
# x^3 + y^3 + z^3 + b2*x*y*z*(x+y+z) + b1*(x+y+z) = 3*R
# x^3 + y^3 + z^3 = 3*R
# (x+y+z)^3 - E2*(x+y+z) - 3*x*y*z = 3*R
# x*y*z = -R

### Sum(x^2) =>
# 2*(x^2 + y^2 + z^2) + E2 + 3*b2*E3 + 3*b1 = 0
# 2*(x+y+z)^2 - 4*E2 + E2 - 3*b2*R + 3*b1 = 0
# E2 + b2*R - b1 = 0
# E2 = -b2*R + b1

### Example
b = c(1, 3)
R = 1
#
x = roots(c(1,0, -b[2]*R + b[1], R))
yz = -R / x;
yz.sum = -x;
yz.diff = sqrt(yz.sum^2 - 4*yz + 0i)
y = (yz.sum + yz.diff)/2
z = (yz.sum - yz.diff)/2
sol = cbind(x, y, z)
sol = rbind(sol, sol[,c(2,3,1)], sol[,c(3,1,2)]) # add all variants

### Test
x^3 + b[2]*x^2*y*z + b[1]*x
y^3 + b[2]*x*y^2*z + b[1]*y
z^3 + b[2]*x*y*z^2 + b[1]*z

### Classic Polynomial

# Trivial: P9 = (P3)^3

### TODO


