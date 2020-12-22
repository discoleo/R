
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogenous Symmetric
###
### draft v.0.3a-pre


### Heterogenous Symmetric Polynomial Systems

### 2 Variables:
### Simple System:
# x^n + P(x, y) = R
# y^n + P(y, x) = R

### Generalised System:
# P(x, y) = 0
# P(y, x) = 0
# where P(x, y) is any polynomial of 2 variables;


###############
### Systems ###

### 2 Variables:
# x^n + P(x, y) = R
# y^n + P(y, x) = R

# 1.) x^3 + b*y = R; [P3 => P6]
# 2.) (x - s)^3 + b*y = R; [P3 => P6: equivalent to non-shifted]
# 3.) x^3 + b*x*y = R; [P3 => trivial P6]
# 4.) (x - s)^3 + b*x*y = R; [P3 => P6]
# 5.) x^3 + b2*x*y + b1*x = R; (P3 => P6)
# 6.) x^3 + b2*x*y + b1*y = R; (TODO: P6)
# 7.) x^3 + b3*x*y + b2*y^2 + b1*y = R; (TODO: P6)
# 8.) TODO: Shift for [5-7];
# 9.) x^3 + b*(x*y)^2 = R; [P4 => P8: based on m3]
# 10.) x^3 + b2*(x*y)^2 + b1*(x*y) = R; [P4 => P8: based on m3]
# 11.) TODO: Shift for [9-10];
# 12.) x^3 + b3*x^2*y + b2*x*y^2 + b1*x = R;
### 2 High-Power Terms:
# B1.) a1*x^3 + a2*y^3 + b*x = R; (TODO: P3 => P6)
# B2.) a1*x^3 + a2*y^3 + b*x*y = R;
# B2.) a1*x^3 + a2*y^3 + b2*x*y + b1*x = R; (TODO: P3 => P6)
# B3.) TODO: Shift for [B1-3];
### Mixt-High-Power: Order 2+1
# M1.) x^2*y + b*x: NO solutions (x != y);
# M2.) x^2*y + b*y: trivial;
# M3.) x^2*y + b2*x^2 + b1*x: trivial (trivial P2);
### Mixt: (x*y)^2
# M4.) x^2*y^2 + b2*x^2 + b1*x: simple (P2 => P4);
# M5.) x^2*y^2 + b3*x^2*y + b2*x^2 + b1*x; (TODO: P2 => P4)
# M6.) TODO: x^2*y^2 + b4*x^2*y + b3*x*y^2 + b2*x^2 + b1*x
# M7.) TODO: x^3*y^2 + b2*x^2 + b1*x = R;
### Mixt: other
# Mm3.1.) x^3 + a*x^3*y = R; (TODO: P3 => P6)
# Mm3.2.) x^3 + a*x^3*y + b1*x*y = R;
# Mm3.3.) x^3 + a*x*y^3 = R; (TODO: P5 => P10)
#         https://www.youtube.com/watch?v=jRwiJTaAu4c
### Mixt: (x*y)^n
# M7.) TODO: (x*y)^3 + ...:
#      (x*y)^3 + x^3 = R;
# M8.) TODO: (x*y)^5 + b*x^3 = R
### Mixt: Order n+1:
# M31.1) x^3*y + b*x: trivial P4;
# M31.2) x^3*y + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P3 => P6; some nice)
# M41.1) x^4*y + b*x; (P5 => P10)
# M43.1) x^4*y^3 + b3*x*y + b2*x^2 + b1*x = R; (trivial P2; base P7)
# M43.2) x^4*y^3 + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P3 => P6)
# M43.3) x^4*y^3 + b5*x^2*y + b4*x*y^2 + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P3 => P6)
### Mixt 2 High-Power Terms:
# MT.1) a1*x^3*y + a2*x*y^3 + b*x + R; [TODO: interesting P6]
### Order 4 & 5:
# O4.1a.) x^4 + b*y = R; (P6 => P12)
# O4.1b.) x^4 + b2*x*y + b1*y = R; (TODO: P6 => P12)
# O4.1c.) x^4 + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P6 => P12; if(b3 == 1) P4 => P8)
# O4.1d.) (x - s)^4 + b*y = R; (P6 => P12; simple shift of O4.1a)
# O4.2a.) x^4 + b*x*y = R; (trivial P8)
# O4.2b.) x^4 + b2*(x*y)^2 + b1*x*y = R; (TODO: trivial P8)
# O4.2c.) TODO: x^4 + b3*(x*y)^3 + b2*(x*y)^2 + b1*x*y = R;
# O5.1.) x^5 + b*(x+y) = 0; (TODO: based on (P4[x^4])^2)
# O5.2.) x^5 + b2*x*y + b1*(x+y) = 0; (TODO: based on (P16)^2)
# O5.3.) x^5 + b*y = R; (P10 => P20)



###############
### History ###

### [branch v.0.3]
#
### v.0.3a-pre:
# - more work on classic polynomials;
# - TODO: cleanup;

### [branch v.0.2]
### P3 Systems:
# - moved to new file:
#   Poly.System.Hetero.Symmetric.S3.R;
### draft v.0.2.d:
# - more work on x[i]^3 + b*x[i+1] = R;
# - added: x^2*y*z + b*y = R;
# - various formatting improvements;
### draft v.0.2.c:
# - solved: x[i]^2 + s*x[i] + b*x[i+1] = R;
# - TODO: correct various bugs [DONE];
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
### draft v.0.1n:
# - solved:
#   -- x^3 + a*x^3*y = R; (TODO: P3 => P6)
#   -- x^3 + a*x*y^3 = R; (TODO: P5 => P10)
### draft v.0.1m - v.0.1.m-ext:
# - various extensions:
#   -- x^4*y^3 Series;
#   -- x^3*y Series [v.0.1.m-ext];
# - more work on older P6 polynomials;
### draft v.0.1l:
# - solved/extension:
#   a1*x^3 + a2*y^3 + b2*x*y + b1*x = R; (TODO: P3 => P6)
# - improved formatting;
### draft v.0.1k & v.0.1k-x:
# - solved: x^5 + b*y = R;
# - worked out various older issues;
# - P6 for x^3 + b2*x*y + b1*x [DONE in v.0.1k-x];
# - some extensions:
#   O4.1b.) x^4 + b2*x*y + b1*y = R; (TODO: P6 => P12)
#   O4.1c.) x^4 + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P6 => P12)
### draft v.0.1j:
# - added a1*x^3*y + a2*x*y^3 + b*x = R;
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
#  -- parametric polynomials [DONE for 1st variant];
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

### Extensions:
# E1: x^3 + b2*x*y + b1*y = R; [P3 => P6]
# E2: x^3 + b3*(x*y)^2 + b2*x*y + b1*y = R; [P4 => P8]
# - are discussed in a separate section;

### Solution:

### Diff =>
# x^3 - y^3 - b1*(x-y) = 0
# (x - y)*(x^2 + y^2 + x*y - b1) = 0
# => x = y *OR* x^2 + y^2 + x*y - b1 = 0;
# =>
# (x+y)^2 - x*y - b1 = 0
# x*y = (x+y)^2 - b1;
# x*y = S^2 - b1;

### Sum =>
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
### Shifted Roots ###
#####################

### (x - s)^3 + b*y

# (x - s)^3 + b1*y = R
# (y - s)^3 + b1*x = R

### Equivalent:
# - trivial shift (only 1 liniar non-shifted term):
#  (x - s)^3 + b1*(y - s) = R - s*b1;
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


###################

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
######################

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

err = x^6 - b[2]*x^5 + (b[2]^2 + 2*b[1])*x^4 - (b[1]*b[2] + 2*R)*x^3 + (b[1]*b[2]^2 + R*b[2] + b[1]^2)*x^2 - 2*b[1]*R*x + R^2
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

err = x^6 - b[2]*x^5 + (b[2]^2 + 2*b[1])*x^4 - (b[1]*b[2] + 2*R)*x^3 + (b[1]*b[2]^2 + R*b[2] + b[1]^2)*x^2 - 2*b[1]*R*x + R^2
round0(err)

### Derivation:
# Z^3 - b2*Z^2 + b1*Z - b1*b2 + R
# x[1]+x[2]+x[3]+x[4]+x[5]+x[6] = b2;
# (x[1]+x[4])*(x[2]+x[5]+x[3]+x[6]) + (x[2]+x[5])*(x[3]+x[6]) = b1
# x[1]*(x[2]+x[3] + 0 +x[5]+x[6]) + x[2]*(x[3]+x[4] + 0 + x[6]) +
#   + x[3]*(x[4]+x[5] + 0) + x[4]*(x[5]+x[6]) + x[5]*x[6] = b1;
# (x[1]+x[4])*(x[2]+x[5])*(x[3]+x[6]) = b[1]*b[2] - R
# E2 = b[1] + x[1]*x[4] + x[2]*x[5] + x[3]*x[6]
#    = b[1] + (x[1]+x[4])^2 + (x[2]+x[5])^2 + (x[3]+x[6])^2 + 3*b[1]
#    = 4*b[1]+ b[2]^2 - 2*b[1]
#    = 2*b1] + b[2]^2;
# E3 = x[1]*x[2]*(x[3]+_x[4]+_x[5]+x[6]) + x[1]*x[3]*(_x[4]+x[5]+_x[6]) + x[1]*x[4]*(_x[5]+_x[6]) + x[1]*x[5]*x[6] +
#    + x[2]*x[3]*(x[4]+_x[5]+_x[6]) + x[2]*x[4]*(_x[5]+x[6]) + _x[2]*x[5]*x[6] +
#    + x[3]*x[4]*(x[5]+_x[6]) + _x[3]*x[5]*x[6] + x[4]*x[5]*x[6]
#    = x[1]*x[4]*(x[2]+x[3]+x[5]+x[6]) + x[2]*x[5]*(x[1]+x[3]+x[4]+x[6]) + x[3]*x[6]*(x[1]+x[2]+x[4]+x[5])
#      + (x[1]+x[4])*(x[2]+x[5])*(x[3]+x[6])
#    = (S1^2 + b[1])*(b[2] - S1) + ...
#    = -S1^3 + b[2]*S1^2 - b[1]*S1 + b[1]*b[2] + ... # Note: Sum = + 3*b[1]*b[2]
#    = -(b[2]^3 - 3*b[1]*b[2] + 3*(b[1]*b[2] - R)) + b[2]*(b[2]^2 - 2*b[1]) - b[1]*b[2] + 3*b[1]*b[2] + b[1]*b[2] - R
#    = b[1]*b[2] + 2*R
# E6 = (S1^2 + b[1])*(S2^2 + b[1])*(S3^2 + b[1])
#    = (S1^2*S2^2 + b[1]*(S1^2 + S2^2) + b[1]^2)*(S3^2 + b[1])
#    = (S1*S2*S3)^2 + b[1]^2*(S1^2 + S2^2 + S3^2) + b[1]*(S1^2*S2^2 + S1^2*S3^2 + S2^2*S3^2) + b[1]^3
#    = (b[1]*b[2] - R)^2 + b[1]^2*(b[2]^2 - 2*b[1]) + b[1]*(b[1]^2 - 2*b[2]*(b[1]*b[2] - R)) + b[1]^3
#    = R^2
# E5 = S1*(S2^2 + b[1])*(S3^2 + b[1]) + S2*(S1^2 + b[1])*(S3^2 + b[1]) + S3*(S1^2 + b[1])*(S2^2 + b[1])
#    = S1*S2^2*S3^2 + b[1]*S1*(S2^2 + S3^2) + b[1]^2*S1 + ...
#    = S1*S2*S3*E2 + b[1]^2*(S1+S2+S3) + b[1]*(S1*S2^2 + S1*S3^2 + S2*S1^2 + S2*S3^2 + S3*S1^2 + S3*S2^2)
#    = b[1]*(b[1]*b[2] - R) + b[1]^2*b[2] + b[1]*(b[2]*b[1] - 3*(b[1]*b[2] - R))
#    = 2*b[1]*R
# E4 = x[1]*x[4]*(x[2]+x[5])*(x[3]+x[6]) + x[2]*x[5]*(x[1]+x[4])*(x[3]+x[6]) + x[3]*x[6]*(x[1]+x[4])*(x[2]+x[5]) +
#      + x[1]*x[2]*x[4]*x[5] + x[1]*x[3]*x[4]*x[6] + x[2]*x[3]*x[5]*x[6]
#    = S2*S3*(S1^2 + b[1]) + S1*S3*(S2^2 + b[1]) + S1*S2*(S3^2 + b[1]) +
#      + (S1^2 + b[1])*(S2^2 + b[1]) + (S1^2 + b[1])*(S3^2 + b[1]) + (S2^2 + b[1])*(S3^2 + b[1])
#    = E3*(S1 + S2 + S3) + b[1]*E2 + (S1^2*S2^2 + ... + 2*b[1]*(S1^2 + S2^2 + S3^2) + 3*b[1]^2)
#    = b[2]*(b[1]*b[2] - R) + b[1]^2 + (E2^2 - 2*E3*(S1 + S2 + S3)) + 2*b[1]*(b[2]^2 - 2*b[1]) + 3*b[1]^2
#    = 3*b[1]*b[2]^2 - R*b[2] + b[1]^2 - 2*(b[1]*b[2] - R)*b[2]
#    = b[1]*b[2]^2 + R*b[2] + b[1]^2
x^6 - b[2]*x^5 + (b[2]^2 + 2*b[1])*x^4 - (b[1]*b[2] + 2*R)*x^3 + (b[1]*b[2]^2 + R*b[2] + b[1]^2)*x^2 - 2*b[1]*R*x + R^2


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
# x*y = S^2 - b2*S - b1

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


### TODO: classical polynomial;
c(1, -(b[2] + b[3]), -b[1] + b[2]^2 + b[3]^2, -(3*b[2]^3 + b[2]*b[3]^2 - 2*b[1]*b[2] - 2*b[1]*b[3] + 2*R),
	NA,
	(b[1] - b[2]*b[3])*R - 2*b[1]^2*b[3] + b[1]*b[2]^3 - b[1]*b[2]^2*b[3], - R*(3*b[1]*b[2] + b[2]^3 - R) - b[1]^3)

### Derivation
# S^3 - (b[2] + b[3])*S^2 - (2*b[1] - b[2]*b[3] + b[2]^2)*S - b[1]*b[2] + b[1]*b[3] + R
# x*y = S^2 - b2*S - b1

# E1 = b[2] + b[3]
# E2 = - (2*b[1] - b[2]*b[3] + b[2]^2) + (S1^2 + S2^2 + S3^2) - b[2]*E1 - 3*b[1]
#    = - (2*b[1] - b[2]*b[3] + b[2]^2) + E1^2 - 2*E2_P3 - b[2]*E1 - 3*b[1]
#    = - (2*b[1] - b[2]*b[3] + b[2]^2) + E1^2 + 2*(2*b[1] - b[2]*b[3] + b[2]^2) - b[2]*E1 - 3*b[1]
#    = -b[1] - b[2]*b[3] + b[2]^2 + E1^2 - b[2]*E1
#    = -b[1] + b[2]*b[3] + 2*b[2]^2 + b[3]^2 - b[2]*E1
#    = -b[1] + b[2]^2 + b[3]^2;
# E3 = E1*(P[1]+P[2]+P[3]) - P[1]*S[1] - P[2]*S[2] - P[3]*S[3] + E3_P3
#    = (b[2] + b[3])*(S1^2+S2^2+S3^2 - b[2]*E1 - 3*b[1]) +
#      - S[1]*(S1^2 - b2*S1 - b1) - S[2]*(S2^2 - b2*S2 - b1) - S[3]*(S3^2 - b2*S3 - b1) + b1*b2 - b1*b3 - R
#    = (b[2] + b[3])*(E1^2 - 2*E2_P3 - b[2]*E1 - 3*b[1]) +
#      - (S[1]^3+S[2]^3+S[3]^3) + b[2]*(E1^2 - 2*E2_P3) + b[1]*E1   + b[1]*b[2] - b[1]*b[3] - R
#    = (b[2] + b[3])*(2*b[2]^2 + b[3]^2 + b[1] - b[2]*b[3]) +
#      - (E1^3 - 3*E1*E2_P3 + 3*E3_P3) + b[2]*(E1^2 - 2*E2_P3) + b[1]*E1   + b[1]*b[2] - b[1]*b[3] - R
#    = (b[2] + b[3])*(-2*b[2]^2 - 5*b[1]) + # + 3*E2_P3
#      - 3*E3_P3 + b[2]*(E1^2 - 2*E2_P3) + b[1]*E1   + b[1]*b[2] - b[1]*b[3] - R
#    = (b[2] + b[3])*(-2*b[2]^2 - 4*b[1] + b[2]*E1) - 2*b[2]*E2_P3 - 2*b[1]*b[2] + 2*b[1]*b[3] + 2*R
#    = (b[2] + b[3])*(b[2]^2 - 4*b[1] + b[2]*b[3]) - 2*b[2]*E2_P3 - 2*b[1]*b[2] + 2*b[1]*b[3] + 2*R
#    = b[2]*(3*b[2]^2 - b[2]*b[3]) + b[3]*(b[2]^2 - 4*b[1] + b[2]*b[3]) - 2*b[1]*b[2] + 2*b[1]*b[3] + 2*R
#    = 3*b[2]^3 + b[2]*b[3]^2 - 2*b[1]*b[2] - 2*b[1]*b[3] + 2*R

# E6 = (S1^2 - b[2]*S1 - b[1])*(S2^2 - b[2]*S2 - b[1])*(S3^2 - b[2]*S3 - b[1])
#    = - b[2]*E3_P3*(S1*S2 + S1*S3 + S2*S3) + b[2]^2*E3_P3*(S1+S2+S3) - b[2]^3*E3_P3 + E3_P3^2 - b[1]*b[2]^2*(S1*S2 + S1*S3 + S2*S3) +
#      + b[1]*b[2]*(S1^2*S2 + S1^2*S3 + S1*S2^2 + S2^2*S3 + S1*S3^2 + S2*S3^2) - b[1]^2*b[2]*(S1+S2+S3) +
#      - b[1]*(S1^2*S2^2 + S1^2*S3^2 + S2^2*S3^2) + b[1]^2*(S1^2 + S2^2 + S3^2) - b[1]^3
#    = - b[2]*E3_P3*E2_P3 + b[2]^2*E3_P3*E1 - b[2]^3*E3_P3 + E3_P3^2 - b[1]*b[2]^2*E2_P3 +
#      + b[1]*b[2]*(E1*E2_P3 - 3*E3_P3) - b[1]^2*b[2]*E1 +
#      - b[1]*(E2_P3^2 - 2*E3_P3*E1) + b[1]^2*(E1^2 - 2*E2_P3) - b[1]^3
#    = - b[2]*E3_P3*E2_P3 + b[2]^2*E3_P3*(b[2] + b[3]) - b[2]^3*E3_P3 + E3_P3^2 - b[1]*b[2]^2*E2_P3 +
#      + b[1]*b[2]*((b[2] + b[3])*E2_P3 - 3*E3_P3) - b[1]^2*b[2]*(b[2] + b[3]) +
#      - b[1]*(E2_P3^2 - 2*E3_P3*(b[2] + b[3])) + b[1]^2*((b[2] + b[3])^2 - 2*E2_P3) - b[1]^3
#    = - b[2]*E3_P3*E2_P3 + b[2]^2*b[3]*E3_P3 + E3_P3^2 +
#      + b[1]*b[2]*b[3]*E2_P3 - b[1]*b[2]*E3_P3 +
#      - b[1]*(E2_P3^2 - 2*b[3]*E3_P3) + b[1]^2*(b[3]^2 + b[2]*b[3] - 2*E2_P3) - b[1]^3
#    = E3_P3*(b[1]*b[2] + 2*b[1]*b[3] + b[2]^3 + E3_P3) +
#      - b[1]*b[2]*b[3]*(3*b[1] - b[2]*b[3] + b[2]^2) +
#      - b[1]*(2*b[1] - b[2]*b[3] + b[2]^2)^2 + b[1]^2*b[3]^2 + 2*b[1]^2*b[2]^2 + 3*b[1]^3
#    = - R*(3*b[1]*b[2] + b[2]^3 - R) - b[1]^3

# E5 = S1*P2*P3 + S2*P1*P3 + S3*P1*P2
#    = S1*(S2^2 - b[2]*S2 - b[1])*(S3^2 - b[2]*S3 - b[1]) + S2*(S1^2 - b[2]*S1 - b[1])*(S3^2 - b[2]*S3 - b[1]) +
#      + S3*(S1^2 - b[2]*S1 - b[1])*(S2^2 - b[2]*S2 - b[1])
#    = E3_P3*(S1*S2 + S1*S3 + S2*S3) - 2*E3_P3*b[2]*(S1 + S2 + S3) + 3*b[2]^2*E3_P3 + 2*b[1]*b[2]*(S1*S2 + S1*S3 + S2*S3) +
#      - b[1]*(S1^2*S2 + S1^2*S3 + S1*S2^2 + S2^2*S3 + S1*S3^2 + S2*S3^2) + b[1]^2*(S1 + S2 + S3)
#    = E3_P3*E2_P3 - 2*b[2]*E3_P3*E1 + 3*b[2]^2*E3_P3 + 2*b[1]*b[2]*E2_P3 - b[1]*(E1*E2_P3 - 3*E3_P3) + b[1]^2*E1
#    = (-b[2]*b[3] + b[1])*(b[1]*b[2] - b[1]*b[3] - R) +
#      - b[1]^2*b[2] + 3*b[1]^2*b[3] - b[1]*b[2]^3 + 2*b[1]*b[2]^2*b[3] - b[1]*b[2]*b[3]^2
#    = -b[2]*b[3]*(b[1]*b[2] - b[1]*b[3] - R) + b[1]*(b[1]*b[2] - b[1]*b[3] - R) +
#      - b[1]^2*b[2] + 3*b[1]^2*b[3] - b[1]*b[2]^3 + 2*b[1]*b[2]^2*b[3] - b[1]*b[2]*b[3]^2
#    = (b[2]*b[3] - b[1])*R + 2*b[1]^2*b[3] - b[1]*b[2]^3 + b[1]*b[2]^2*b[3]


# E4 = ...; # TODO



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

#######################
### Base: x^n       ###
### + x^j*y^k Terms ###
#######################

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
# Case: b2 - b3 + 1 == 0; TODO: verify;
# Case: b2 == 1;

solve.htxy2 = function(b, R) {
	if(abs(b[3] - b[2] - 1) < 1E-12) {
		x.sum = sqrt(-b[1] + 0i) * c(1, -1)
		xy = (x.sum^3 + b[1]*x.sum - 2*R) / x.sum / (b[2] + b[3] - 3)
	} else {
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

### TODO


#############################
#############################

#############################
### High-Power Terms: > 1 ###
#############################

##########################
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


#########################
### High-Power Terms: > 1

### Variant:
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

### Variant:
### a1*x^3 + a2*y^3 + b2*x*y + b1*x

# a1*x^3 + a2*y^3 + b2*x*y + b1*x = R
# a2*x^3 + a1*y^3 + b2*x*y + b1*y = R

### Solution

### Diff =>
# (a1 - a2)*(x^3 - y^3) = -b1*(x - y)
# (x - y)*(x^2 + y^2 + x*y + b1/(a1 - a2)) = 0
# Case: x != y =>
# S^2 - x*y + b1/(a1 - a2) = 0
# x*y = S^2 + b1/(a1 - a2)

### Sum =>
# (a1+a2)*(x^3 + y^3) + 2*b2*x*y + b1*(x+y) = 2*R
# (a1+a2)*(S^3 - 3*x*y*S) + 2*b2*(S^2 + b1/(a1 - a2)) + b1*S - 2*R = 0
# (a1+a2)*(S^3 - 3*S*(S^2 + b1/(a1 - a2))) + 2*b2*S^2 + 2*b1*b2/(a1 - a2) + b1*S - 2*R
# (a1+a2)*(-2*S^3 - 3*S*b1/(a1 - a2)) + 2*b2*S^2 + b1*S + 2*b1*b2/(a1 - a2) - 2*R
# 2*(a1+a2)*S^3 - 2*b2*S^2 + 3*b1*(a1+a2)/(a1 - a2)*S - b1*S - 2*b1*b2/(a1 - a2) + 2*R
# 2*(a1-a2)*(a1+a2)*S^3 - 2*b2*(a1-a2)*S^2 + 3*b1*(a1+a2)*S - b1*(a1-a2)*S - 2*b1*b2 + 2*(a1-a2)*R

solve.ht2a = function(b, a, R) {
	a.s = a[1]+a[2]; a.d = a[1]-a[2];
	coeffs = c(2*a.d*a.s, -2*b[2]*a.d, 3*b[1]*a.s - b[1]*a.d, -2*b[1]*b[2] + 2*a.d*R)
	r.sum = roots(coeffs)
	xy = r.sum^2 + b[1]/(a[1] - a[2])
	r.diff = sqrt(r.sum^2 - 4*xy + 0i)
	x = (r.sum + r.diff)/2
	y = (r.sum - r.diff)/2
	sol = cbind(x, y) # TODO: add also x = y cases;
	sol = round0(rbind(sol, sol[,2:1]))
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}

### Example:
a = c(2/3, 1/3)
b = c(3, -1)
R = 1
#
sol = solve.ht2a(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
a[1]*x^3 + a[2]*y^3 + b[2]*x*y + b[1]*x
a[2]*x^3 + a[1]*y^3 + b[2]*x*y + b[1]*y


### Example 2:
a = c(2/3, -1/3)
b = c(1, -3)
R = 1
#
sol = solve.ht2a(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
a[1]*x^3 + a[2]*y^3 + b[2]*x*y + b[1]*x
a[2]*x^3 + a[1]*y^3 + b[2]*x*y + b[1]*y


###################################
###################################

###############
### Order 4 ###
###############

### x^4 + b*y
### Extension 1: x^4 + b2*xy + b1*y = R
### Extension 2: x^4 + b3*(x*y)^2 + b2*x*y + b1*y = R

# x^4 + b1*y = R
# y^4 + b1*x = R

### Solution:

### Diff =>
# x^4 - y^4 - b1*(x-y) = 0
# (x - y)*((x+y)*(x^2+y^2) - b1) = 0
# => x = y *OR* (x+y)*(x^2+y^2) - b1 = 0;
# Case: x != y =>
# (x+y)^2 - 2*x*y = b1 / S
# x*y = (S^2 - b1/S)/2

### Sum =>
# x^4 + y^4 + b1*S - 2*R = 0
# S^4 - 4*x*y*S^2 + 2*(x*y)^2 + b1*S - 2*R
# S^4 - 2*(S^4 - b1*S) + (S^2 - b1/S)^2/2 + b1*S - 2*R
# -S^4 + 2*b1*S + (S^2 - b1/S)^2/2 + b1*S - 2*R
# S^4 - 3*b1*S - (S^2 - b1/S)^2/2 + 2*R
# 2*S^6 - 6*b1*S^3 - (S^3 - b1)^2 + 4*R*S^2
# S^6 - 4*b1*S^3 + 4*R*S^2 - b1^2

### Simple extensions:
### E1: x^4 + b2*x*y + b1*y = R
# S^6 - 2*b2*S^4 - 4*b1*S^3 + 4*R*S^2 + 2*b1*b2*S - b1^2 = 0
### E2: x^4 + b3*(x*y)^2 + b2*x*y + b1*y = R
# S^6 - b3*(S^3 - b1)^2 - 2*b2*S^4 - 4*b1*S^3 + 4*R*S^2 + 2*b1*b2*S - b1^2 = 0
# (1 - b3)*S^6 - 2*b2*S^4 - (4*b1 - 2*b1*b3)*S^3 + 4*R*S^2 + 2*b1*b2*S - b1^2 - b1^2*b3 = 0
# (b3 - 1)*S^6 + 2*b2*S^4 - (2*b1*b3 - 4*b1)*S^3 - 4*R*S^2 - 2*b1*b2*S + b1^2 + b1^2*b3 = 0

solve.ht4 = function(b, R) {
	r.sum = if(length(b) == 1) { roots(c(1,0,0, -4*b[1], 4*R, 0, - b[1]^2));
		} else if(length(b) == 2) { roots(c(1,0,-2*b[2], -4*b[1], 4*R, 2*b[1]*b[2], - b[1]^2));
		} else roots(c(b[3]-1, 0, 2*b[2], -2*b[1]*b[3]+4*b[1], -4*R, -2*b[1]*b[2], b[1]^2 + b[1]^2*b[3]));
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
b = 2
R = 1
#
sol = solve.ht4(b, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^4 + b[1]*y
y^4 + b[1]*x

### Classical Polynomial
err = x^12 - b[1]*x^9 - 3*R*x^8 + b[1]^2*x^6 + 2*R*b[1]*x^5 + 3*R^2*x^4 - b[1]^3*x^3 - R*b[1]^2*x^2 - R^2*b[1]*x + b[1]^4 - R^3
round0(err)
# - only for simple case!

# P6 => P12;
# when b3 == 1: P4 => P8;

### Example 2: Extended version
b = c(2, 1)
R = 1
#
sol = solve.ht4(b, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^4 + b[2]*x*y + b[1]*y
y^4 + b[2]*x*y + b[1]*x


### Example 3: Extended version
b = c(2, 1, -1)
R = 1
#
sol = solve.ht4(b, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^4 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
y^4 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Classical Polynomial
### Derivation:
# b1*y = R - x^4
# (R - x^4)^4/b1^4 + b1*x - R = 0
(R - x^4)^4 + b[1]^5*x - R*b[1]^4

x^16 - 4*R*x^12 + 6*R^2*x^8 - 4*R^3*x^4 + b[1]^5*x - R*b[1]^4 + R^4
(x^4 + b[1]*x - R)*(x^12 - b[1]*x^9 - 3*R*x^8 + b[1]^2*x^6 + 2*R*b[1]*x^5 + 3*R^2*x^4 - b[1]^3*x^3 - R*b[1]^2*x^2 - R^2*b[1]*x + b[1]^4 - R^3)

# TODO: for extended version;


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

###############
### x^4 + b*x*y

# x^4 + b1*x*y = R
# y^4 + b1*x*y = R

### Solution:

# Trivial solutions: x = y *OR* x = -y;

### Diff =>
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

### Diff =>
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
#########################

###################
### x^2*y + b*x ###

# x^2*y + b1*x = R
# y^2*x + b1*y = R

### Solution: *NO*

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

### NO Solution!


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

### Extension:
# x^2*y + b3*x*y + b2*x^2 + b1*x = R

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

### Extensions:
### E1: x^2*y + b3*x*y + b2*x^2 + b1*x = R
### Sum =>
# b2^2*Z - b3*(b2*Z + b1) + b1*b2 - R = 0
# (b2^2 - b2*b3)*Z + b1*b2 - b1*b3 - R = 0

solve.ht21 = function(b, R) {
	if(length(b) == 2) {
		r.sum = (R - b[1]*b[2]) / b[2]^2
	} else {
		r.sum = (R - b[1]*b[2] + b[1]*b[3]) / (b[2]^2 - b[2]*b[3])
	}
	xy = -b[2]*r.sum - b[1]
	r.diff = sqrt(r.sum^2 - 4*xy + 0i)
	x = (r.sum + r.diff)/2
	y = (r.sum - r.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	sol
}

### Example:
b = c(1, 3)
R = -1
#
sol = solve.ht21(b, R)
x = sol[,1]; y = sol[,2];
sol

### Test
x^2*y + b[2]*x^2 + b[1]*x
y^2*x + b[2]*y^2 + b[1]*y

### Classical Polynomial:
# trivial P2;
b[2]*x^2 - (R/b[2] - b[1])*x - R


### Example 2: Extension
b = c(1, 3, 1)
R = -1
#
sol = solve.ht21(b, R)
x = sol[,1]; y = sol[,2];
sol

### Test
x^2*y + b[3]*x*y + b[2]*x^2 + b[1]*x
y^2*x + b[3]*x*y + b[2]*y^2 + b[1]*y


############################

##################
### x^3*y Term ###
##################

### x^3*y + b3*(x*y)^2 + b2*x*y + b1*y

# x^3*y + b3*(x*y)^2 + b2*x*y + b1*y = R
# y^3*x + b3*(x*y)^2 + b2*x*y + b1*x = R

### Solution:

### Diff =>
# x*y*(x - y)*(x + y) - b1*(x-y) = 0
# (x - y)*(x*y*(x+y) - b1) = 0
# Case: x != y
# x*y = b1 / (x+y)
# x*y = b1 / S;

### Sum =>
# x*y*(x^2 + y^2) + 2*b3*(x*y)^2 + 2*b2*x*y + b1*(x+y) = 2*R
# x*y*(S^2 - 2*x*y) + 2*b3*(x*y)^2 + 2*b2*x*y + b1*S - 2*R = 0
# x*y*S^2*(S^2 - 2*x*y) + 2*b3*(x*y)^2*S^2 + 2*b2*x*y*S^2 + b1*S^3 - 2*R*S^2 = 0
# b1*(S^3 - 2*b1) + 2*b1^2*b3 + 2*b1*b2*S + b1*S^3 - 2*R*S^2 = 0
# b1*S^3 - R*S^2 + b1*b2*S + b1^2*b3 - b1^2 = 0

solve.htx3y = function(b, R) {
	if(R == 0 && b[1] == 0) return(list(sol=NA, p=NA))
	r.sum = roots(c(b[1], - R, b[1]*b[2], b[1]^2*b[3] - b[1]^2))
	r.sum = round0(r.sum)
	r.sum = r.sum[ r.sum != 0 ] # avoid division by 0
	xy = b[1] / r.sum
	r.diff = sqrt(r.sum^2 - 4*xy + 0i)
	x = (r.sum + r.diff)/2
	y = (r.sum - r.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}

### Example:
b = c(1/2, -1, 2)
R = 1
#
sol = solve.htx3y(b, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^3*y + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
y^3*x + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x
#
err = -0.25 + 3*x - 4*x^2 - 2*x^5 + x^6
round0(err)


### Example 2:
b = c(1,-1,-2)
R = 1
#
sol = solve.htx3y(b, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^3*y + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
y^3*x + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Example 3:
b = c(1, -3, 2)
R = 1
#
sol = solve.htx3y(b, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
err = -1 + 7*x - 12*x^2 + x^3 - x^5 + x^6
round0(err)


### Examples:
R = 1
p = sapply(-6:6, function(r) print(solve.htx3y(c(1, r, 2), R)$p))
#
b = c(1, -4, 2)
sol = solve.htx3y(b, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol
-1 + 9*x - 19*x^2 - x^5 + x^6

###
p = sapply(-6:6, function(r) print(solve.htx3y(c(r, -3, 2), r)$p))
#
r = -3; sol = solve.htx3y(c(r, -3, 2), r); x = sol$sol[,1]
-9 - 21*x - 15*x^3 - x^5 + x^6


###################
### Simple variant:
### x^3*y + b*x

# x^3*y + b1*x = R
# y^3*x + b1*y = R

# trivial P4;

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


############################

############################
### x^4*y^3 + b*x Series ###
############################

### Simple variant:
### x^4*y^3 + b2*(x+y)^2 + b1*x

# x^4*y^3 + b2*(x+y)^2 + b1*x = R
# y^3*x^4 + b2*(x+y)^2 + b1*y = R

# trivial polynomial;

### Solution:

### Diff =>
# (x*y)^3*(x - y) + b1*(x-y) = 0
# (x - y)*((x*y)^3 + b1) = 0
# Case: x != y
# x*y = (-b1)^(1/3) * m3.all

### Sum =>
# (x*y)^3*(x+y) + 2*b2*(x+y)^2 + b1*(x+y) = 2*R
# -b1*S + 2*b2*S^2 + b1*S - 2*R = 0
# b2*S^2 - R = 0
# S^2 = R / b2

m3.all = unity(3, all=T)

### Example
b = c(3, 1)
R = 1
#
xy = (-b[1] + 0i)^(1/3) * m3.all
x.sum = sqrt(R / b[2] + 0i)
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol

### Test
x^4*y^3 + b[2]*(x+y)^2 + b[1]*x
y^4*x^3 + b[2]*(x+y)^2 + b[1]*y

### Classical Polynomial:
round0.p(poly.calc(sol[,1]))
# trivial
err = -3 - x^3 + 3*x^4 - 3*x^5 + x^6
round0(err)

# trivial:
-b[1] + x^3*(x^3 - sqrt(R/b[2]))^3


##################
### Extensions ###

### Extension 1:

### x^4*y^3 + b3*x*y + b2*x^2 + b1*x

# x^4*y^3 + b3*x*y + b2*x^2 + b1*x = R
# y^3*x^4 + b3*x*y + b2*y^2 + b1*y = R

### Solution:

# "Trivial" solution:: x = y
# x^7 + (b2+b3)*x^2 + b1*x - R = 0
# &
# Trivial remaining P2;

### Diff =>
# (x*y)^3*(x - y) + b2*(x-y)*(x+y) + b1*(x-y) = 0
# (x - y)*((x*y)^3 + b2*(x+y) + b1) = 0
# Case: x != y
# (x*y)^3 = - (b2*S + b1);

### Sum =>
# (x*y)^3*(x+y) + 2*b3*x*y + b2*(x^2 + y^2) + b1*(x+y) = 2*R
# - (b2*S + b1)*S + 2*b3*x*y + b2*(S^2 - 2*x*y) + b1*S - 2*R = 0
# 2*b3*x*y - 2*b2*x*y - 2*R = 0
# (b3-b2)*x*y = R
# x*y = R / (b3-b2)
# =>
# (x*y)^3 = R^3 / (b3-b2)^3
# b2*S + b1 + R^3 / (b3-b2)^3 = 0

### Example
b = c(3,2,1)
R = 1
#
xy = R / (b[3] - b[2])
x.sum = - (b[1] + R^3 / (b[3] - b[2])^3) / b[2]
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
x = sol[,1]; y = sol[,2]
sol

### Test
x^4*y^3 + b[3]*x*y + b[2]*x^2 + b[1]*x
x^3*y^4 + b[3]*x*y + b[2]*y^2 + b[1]*y

# trivial P2 polynomial;

# all roots
x2 = roots(c(1, 0,0,0,0, (b[2]+b[3]), b[1], - R))
x = c(x, x2); y = c(y, x2)


###############
### Extension 2

### x^4*y^3 + b3*(x*y)^2 + b2*x*y + b1*y

# x^4*y^3 + b3*(x*y)^2 + b2*x*y + b1*y = R
# y^3*x^4 + b3*(x*y)^2 + b2*x*y + b1*x = R

### Solution:

# "Trivial" solution: x = y
# x^7 + b3*x^4 + b2*x^2 + b1*x - R = 0

### Diff =>
# (x*y)^3*(x - y) - b1*(x-y) = 0
# (x - y)*((x*y)^3 - b1) = 0
# Case: x != y
# (x*y)^3 = b1;

### Sum =>
# (x*y)^3*(x+y) + 2*b3*(x*y)^2 + 2*b2*x*y + b1*(x+y) = 2*R
# ((x*y)^3 + b1)*S = 2*R - 2*b3*(x*y)^2 - 2*b2*x*y
# 2*b1*S = 2*R - 2*b3*(x*y)^2 - 2*b2*x*y
# b1*S = R - b3*(x*y)^2 - b2*x*y

solve.htMixt = function(b, R) {
	if(length(b) == 3) {
		xy = roots(c(1,0,0, -b[1]))
		r.sum = (R - b[3]*xy^2 - b[2]*xy) / b[1]
	} else {
		if(length(b) == 4) b = c(b, 0);
		xy = roots(c(1,0, b[5]-b[4], -b[1]))
		div = (b[4]*xy + b[1])
		r.sum = (R - b[3]*xy^2 - b[2]*xy) / div
	}
	r.diff = sqrt(r.sum^2 - 4*xy + 0i)
	x = (r.sum + r.diff)/2
	y = (r.sum - r.diff)/2
	x = round0(x)
	sol = cbind(x, y)
	sol = sol[x != 0, ] # x == 0 => y -> Inf!
	sol = rbind(sol, sol[,2:1])
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}


### Example:
# b1 == 1 is trivial;
b = c(1, 2, 3)
R = 1
#
sol =solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^4*y^3 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
x^3*y^4 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Example 2:
b = c(1/2, 2, -3)
R = 1
#
sol =solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^4*y^3 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
x^3*y^4 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Example 3:
b = c(2,2,-2)
R = 2
#
sol =solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^4*y^3 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
x^3*y^4 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Classic Polynomial
round0.p(poly.calc(sol[,1]))

# TODO


###############
### Extension 3

### x^4*y^3 + b5*x^2*y + b4*x*y^2 + b3*(x*y)^2 + b2*x*y + b1*y

# x^4*y^3 + b5*x^2*y + b4*x*y^2 + b3*(x*y)^2 + b2*x*y + b1*y = R
# y^3*x^4 + b5*x^2*y + b4*x*y^2 + b3*(x*y)^2 + b2*x*y + b1*x = R

### Solution:

# "Trivial" solution: x = y
# x^7 + b3*x^4 + (b4+b5)*x^3 + b2*x^2 + b1*x - R = 0

### Diff =>
# (x*y)^3*(x - y) + (b5-b4)*x*y*(x-y) - b1*(x-y) = 0
# (x - y)*((x*y)^3 + (b5-b4)*x*y - b1) = 0
# Case: x != y
# (x*y)^3 + (b5-b4)*x*y - b1 = 0;

### Sum =>
# (x*y)^3*(x+y) + (b4+b5)*x*y*(x+y) + 2*b3*(x*y)^2 + 2*b2*x*y + b1*(x+y) = 2*R
# ((x*y)^3 + (b4+b5)*x*y + b1)*S = 2*R - 2*b3*(x*y)^2 - 2*b2*x*y
# (b4*x*y + b1)*S = R - b3*(x*y)^2 - b2*x*y


### Example:
b = c(1, 2, 3, -1, -2)
R = 1
#
sol =solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^4*y^3 + b[5]*x^2*y + b[4]*x*y^2 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
x^3*y^4 + b[5]*x*y^2 + b[4]*x^2*y + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Example 2:
b = c(1, 0, -1, -1, 1)
R = 1
#
sol =solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^4*y^3 + b[5]*x^2*y + b[4]*x*y^2 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
x^3*y^4 + b[5]*x*y^2 + b[4]*x^2*y + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x
#
err = 1 - 4*x - 2*x^2 + 2*x^3 - 2*x^5 + x^6
round0(err)


### Example 3:
b = c(1,1,1,0,-1)
R = 1
#
sol =solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
err = 1 + 4*x + 6*x^2 - 4*x^4 - x^5 + x^6
round0(err)


### Classic Polynomial
# TODO


####################
####################

####################
### x^2*y^2 Term ###
####################

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
###  + b2*x^2 + b1*x

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
##########################

######################
### (x*y)^3 Series ###
######################

### (x*y)^3 + b*x^3

# (x*y)^3 + b1*x^3 = R
# (x*y)^3 + b1*y^3 = R


m3 = unity(3, all=F)

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

##########################
### 2 Highest Terms:   ###
###   Simple + Mixt    ###
##########################

###############
### x^3 + x^3*y

# x^3 + a*x^3*y = R
# y^3 + a*y^3*x = R

### Extension:
# x^3 + a*x^3*y + b1*x*y = R

### Solution:

### Diff =>
# (x-y)*(x^2 + y^2 + x*y + a*x*y*(x+y)) = 0
# Case: x != 0 =>
# x^2 + y^2 + x*y + a*x*y*(x+y) = 0
# S^2 - x*y + a*x*y*S = 0
# x*y = - S^2 / (a*S - 1)

### Sum =>
# x^3 + y^3 + a*x*y*(x^2 + y^2) - 2*R = 0
# S^3 - 3*x*y*S + a*x*y*(S^2 - 2*x*y) - 2*R
# S^3 + 3*S^3 / (a*S - 1) - a*S^2*(S^2 - 2*x*y) / (a*S - 1) - 2*R
# a*S^4 - S^3 + 3*S^3 - a*S^2*(S^2 - 2*x*y) - 2*R*(a*S - 1)
# 2*S^3 + 2*a*S^2*x*y - 2*R*a*S + 2*R
# 2*S^3 - 2*a*S^4 / (a*S - 1) - 2*R*a*S + 2*R
# 2*a*S^4 - 2*S^3 - 2*a*S^4 - 2*R*a^2*S^2 + 2*R*a*S + 2*a*R*S - 2*R
# - 2*S^3 - 2*R*a^2*S^2 + 4*a*R*S - 2*R
# S^3 + R*a^2*S^2 - 2*a*R*S + R

### Extension: + b1*x*y
# (b1*a+1)*S^3 + (a^2*R-b1)*S^2 - 2*a*R*S + R

solve.x3y = function(a, R, b=0, type=1) {
	# type = 1: simple;
	if(type == 2) {
		coeff = c(b[1]*a[1]+1, a[1]^2*R - b[1], - 2*a[1]*R, R)
	} else if(type == 1) {
		coeff = c(1, R*a[1]^2, - 2*a[1]*R, R)
	} else if(type == 5) {
		coeff = c(a[1]^2, - a[1], -1, - a[1]^2*R, - 2*a[1]*R, - R)
	}
	x.s = roots(coeff)
	xy = if(type < 5) x.s^2 / (1 - a[1]*x.s) else x.s^2 / (a[1]*x.s + 1);
	x.d = sqrt(x.s^2 - 4*xy + 0i)
	x = (x.s + x.d)/2
	y = (x.s - x.d)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}

### Example
a = 3
R = 1
#
sol = solve.x3y(a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^3 + a[1]*x^3*y
y^3 + a[1]*y^3*x

### Classic Polynomial
round0.p(poly.calc(sol$sol[,1]))
err = 1 + 3*x - 2*x^3 - 3*x^4 + 9*x^5 + x^6
round0(err)

###
R = 1
l = sapply(-6:6, function(a) print(solve.x3y(a, R)$p))

### Extension

### Example 2: Extension
a = -3
b = c(1/2)
R = 2
#
sol = solve.x3y(a, R=R, b=b, type=2)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^3 + a[1]*x^3*y + b[1]*x*y
y^3 + a[1]*y^3*x + b[1]*x*y



################
### Alternative:
### x^3 + x*y^3

# x^3 + a*x*y^3 = R
# y^3 + a*x^3*y = R

# Cheol-Hyun Cho:
# Fukaya category for Landau-Ginzburg orbifolds and Berglund-Hübsch
# homological mirror symmetry for curve singularities
# [from 53:30] https://www.youtube.com/watch?v=jRwiJTaAu4c

### TODO: add x*y & (x*y)^2 terms;

### Diff =>
# (x-y)*(x^2 + y^2 + x*y - a*x*y*(x+y)) = 0
# Case: x != 0 =>
# x^2 + y^2 + x*y - a*x*y*(x+y) = 0
# S^2 - x*y - a*x*y*S = 0
# x*y = S^2 / (a*S + 1)

### Sum =>
# x^3 + y^3 + a*x*y*(x^2 + y^2) - 2*R = 0
# S^3 - 3*x*y*S + a*x*y*(S^2 - 2*x*y) - 2*R
# S^3 - 3*S^3 / (a*S + 1) + a*S^2*(S^2 - 2*x*y) / (a*S + 1) - 2*R
# a*S^4 + S^3 - 3*S^3 + a*S^2*(S^2 - 2*x*y) - 2*R*(a*S + 1)
# a*S^4 - S^3 - a*S^2*x*y - a*R*S - R
# a*S^4 - S^3 - a*S^4/(a*S + 1) - a*R*S - R
# a^2*S^5 - a*S^4 - S^3 - a^2*R*S^2 - 2*a*R*S - R

### Example
a = 1
R = 1
#
sol = solve.x3y(a, R, type=5)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^3 + a[1]*x*y^3
y^3 + a[1]*y*x^3

### Classic Polynomial
round0.p(poly.calc(sol$sol[,1]))
err = 1 - 3*x + 3*x^2 - 3*x^3 + 4*x^4 - 3*x^5 + 3*x^6 - 2*x^7 + x^8 - x^9 + x^10
round0(err)


####################
####################


####################
### Mixt Order 5 ###
####################

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


############################

############################
### 2 Mixed Highest Term ###
############################

###############################
### a1*x^j*y^k + a2*x^k*y^j ###
###############################

### a1*x^3*y + a2*x*y^3 + b*x

# a1*x^3*y + a2*x*y^3 + b1*x = R
# a1*x*y^3 + a2*x^3*y + b1*y = R

### Solution

### Diff =>
# (a1-a2)*x*y*(x^2 - y^2) + b1*(x-y) = 0
# (x-y)*((a1-a2)*x*y*Z + b1) = 0
# Case: x != y
# (a1-a2)*x*y*Z + b1 = 0
# x*y*Z = -b1/(a1-a2);

### Sum =>
# (a1+a2)*x*y*(x^2 + y^2) + b1*(x+y) = 2*R
# (a1+a2)*x*y*(Z^2 - 2*x*y) + b1*Z - 2*R = 0
# (a1+a2)*x*y*Z*(Z^3 - 2*x*y*Z) + b1*Z^3 - 2*R*Z^2 = 0
# =>
# -(a1+a2)*b1/(a1-a2)*(Z^3 + 2*b1/(a1-a2)) + b1*Z^3 - 2*R*Z^2
# -(a1+a2)*b1*((a1-a2)*Z^3 + 2*b1) + b1*(a1-a2)^2*Z^3 - 2*(a1-a2)^2*R*Z^2
# -b1*(a1+a2)*(a1-a2)*Z^3 + b1*(a1-a2)^2*Z^3 - 2*(a1-a2)^2*R*Z^2 - 2*b1^2*(a1+a2)
# 2*b1*a2*(a1-a2)*Z^3 + 2*(a1-a2)^2*R*Z^2 + 2*b1^2*(a1+a2)
# b1*a2*(a1-a2)*Z^3 + (a1-a2)^2*R*Z^2 + b1^2*(a1+a2)

solve.ht2a = function(b, a, R) {
	x.sum = roots(c(b[1]*a[2]*(a[1]-a[2]), (a[1]-a[2])^2*R, 0, b[1]^2*(a[1]+a[2])))
	x.sum = x.sum[x.sum != 0]
	if(length(x.sum) == 0) { return(list(sol=NA, p=NA)); }
	xy = -b[1]/(a[1]-a[2]) / x.sum
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}

### Example 1:
b = 3
a = c(1,2)
R = 1
#
sol = solve.ht2a(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
a[1]*x^3*y + a[2]*x*y^3 + b[1]*x
a[1]*x*y^3 + a[2]*x^3*y + b[1]*y

### Classic Polynomial
round0.p(poly.calc(sol$sol[,1]))


### Example 2:
b = 3
a = c(1,2)
R = 6
#
sol = solve.ht2a(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
a[1]*x^3*y + a[2]*x*y^3 + b[1]*x
a[1]*x*y^3 + a[2]*x^3*y + b[1]*y

### Classic Polynomial
round0.p(poly.calc(sol$sol[,1]))


###
a = c(1,2)
p = sapply(-6:6, function(r) print(solve.ht2a(6*r, a, 6)$p))
p = sapply(-6:6, function(r) print(solve.ht2a(6*r, c(1,1/2), 9)$p)) # interesting
p = sapply(-6:6, function(r) print(solve.ht2a(2*r, c(-1/3,2/3), 12)$p)) # interesting
p = sapply(-6:6, function(r) print(solve.ht2a(r, c(-1/3,2/3), 12)$p)) # interesting
# type: x - x^3 + x^5


a1 = a[1]; a2 = a[2]; b1 = b[1]; R = R[1];
#
(a1*x^4 + a2*x^4 + b1*x - R) *
	((- a2^3*b1^3) +
	(- 3*R^2*a1*a2^3 + 3*R^2*a1^2*a2^2 - R^2*a1^3*a2 + R^2*a2^4)*x +
	(- 3*R*a1^2*a2^2*b1 + 2*R*a1^3*a2*b1 + R*a2^4*b1)*x^2 +
	(a1*a2^3*b1^2 + 2*a1^2*a2^2*b1^2 - a1^3*a2*b1^2 - 2*a2^4*b1^2)*x^3 +
	(- 2*R*a1*a2^4 + 2*R*a1^3*a2^2 - R*a1^4*a2 + R*a2^5)*x^5 +
	(a1*a2^4*b1 + a1^2*a2^3*b1 - a1^3*a2^2*b1 - a2^5*b1)*x^6 )


### full: un-decomposed
(R*a2^3*b1^3) +
(3*R^3*a1*a2^3 - 3*R^3*a1^2*a2^2 + R^3*a1^3*a2 - R^3*a2^4 - a2^3*b1^4)*x +
(- 3*R^2*a1*a2^3*b1 + 6*R^2*a1^2*a2^2*b1 - 3*R^2*a1^3*a2*b1)*x^2 +
(- R*a1*a2^3*b1^2 - 5*R*a1^2*a2^2*b1^2 + 3*R*a1^3*a2*b1^2 + 3*R*a2^4*b1^2)*x^3 +
(2*a1^2*a2^2*b1^3 - a1^3*a2*b1^3 - 3*a2^4*b1^3)*x^4 +
(- 2*R*a1*a2^4*b1 - 4*R*a1^2*a2^3*b1 + 2*R*a1^3*a2^2*b1 + R*a1^4*a2*b1 + 3*R*a2^5*b1)*x^6 +
(4*a1^2*a2^3*b1^2 - a1^4*a2*b1^2 - 3*a2^5*b1^2)*x^7 +
(- R*a1*a2^5 - 2*R*a1^2*a2^4 + 2*R*a1^3*a2^3 + R*a1^4*a2^2 - R*a1^5*a2 + R*a2^6)*x^9 +
(2*a1^2*a2^4*b1 - a1^4*a2^2*b1 - a2^6*b1)*x^10



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


###############
###############

#############
### x^5 + b*y

# x^5 + b*y = R
# y^5 + b*x = R

### Solution

### Diff =>
# x^5 - y^5 - b*(x-y) = 0
# (x-y)*(x^4 + y^4 + x*y*(x^2 + y^2) + (x*y)^2 - b) = 0
# Case: x != y =>
# x^4 + y^4 + x*y*(x^2 + y^2) + (x*y)^2 - b = 0
# S^4 - 4*x*y*(x+y)^2 + 2*(x*y)^2 + x*y*(S^2 - 2*x*y) + (x*y)^2 - b = 0
# S^4 - 3*x*y*S^2 + (x*y)^2 - b = 0

### Sum =>
# x^5 + y^5 + b*(x+y) = 2*R
# (x+y)*(x^4 + y^4 - x*y*(x^2 + y^2) + (x*y)^2 + b) - 2*R = 0
# S*(S^4 - 4*x*y*S^2 + 2*(x*y)^2 - x*y*(S^2 - 2*x*y) + (x*y)^2 + b) - 2*R = 0
# S*(S^4 - 5*x*y*S^2 + 5*(x*y)^2 + b) - 2*R = 0

### Sum - 5*S*Diff =>
# S*(S^4 - 5*x*y*S^2 + 5*(x*y)^2 + b - 5*S^4 + 15*x*y*S^2 - 5*(x*y)^2 + 5*b) - 2*R = 0
# S*(-4*S^4 + 10*x*y*S^2 + 6*b) - 2*R = 0
# S*(4*S^4 - 10*x*y*S^2 - 6*b) + 2*R = 0
# S*(2*S^4 - 5*x*y*S^2 - 3*b) + R = 0
# 5*x*y*S^3 = 2*S^5 - 3*b*S + R

### => Diff =>
# S^4 - 3*x*y*S^2 + (x*y)^2 - b = 0 # * S^6
# S^10 - 3*x*y*S^8 + (x*y)^2*S^6 - b*S^6 = 0
# 25*S^10 - 75*x*y*S^8 + 25*(x*y)^2*S^6 - 25*b*S^6 = 0
# 25*S^10 - 15*S^5*(2*S^5 - 3*b*S + R) + (2*S^5 - 3*b*S + R)^2 - 25*b*S^6 = 0
# 5*S^10 - (2*S^5 - 3*b*S + R)^2 - 20*b*S^6 + 15*R*S^5 = 0
# 5*S^10 - (4*S^10 + 9*b^2*S^2 + R^2 - 12*b*S^6 + 4*R*S^5 - 6*b*R*S) - 20*b*S^6 + 15*R*S^5 = 0
# S^10 - 8*b*S^6 + 11*R*S^5 - 9*b^2*S^2 + 6*b*R*S - R^2 = 0

solve.ht5Basic = function(b, R) {
	x.sum = roots(c(1, 0,0,0, -8*b[1], 11*R, 0,0, - 9*b[1]^2, 6*b[1]*R, - R^2))
	xy = (2*x.sum^5 - 3*b[1]*x.sum + R) / x.sum^3 / 5
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(round0(x), round0(y))
	sol = rbind(sol, sol[,2:1])
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}

### Example:
b = 3
R = 1
#
sol = solve.ht5Basic(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^5 + b[1]*y
y^5 + b[1]*x


### Example 2:
b = -1
R = 1
#
sol = solve.ht5Basic(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^5 + b[1]*y
y^5 + b[1]*x

### Classic Polynomial
x^20 - b[1]*x^16 - 4*R*x^15 + b[1]^2*x^12 + 3*R*b[1]*x^11 + 6*R^2*x^10 - b[1]^3*x^8 +
 - 2*R*b[1]^2*x^7 - 3*R^2*b[1]*x^6 - 4*R^3*x^5 + b[1]^4*x^4 + R*b[1]^3*x^3 +
 + R^2*b[1]^2*x^2 + R^3*b[1]*x + R^4 - b[1]^5

### Example 3: R^4 - b[1]^5 == 0
b = 1
R = -1
#
sol = solve.ht5Basic(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^5 + b[1]*y
y^5 + b[1]*x
# x*(x+1) * P18
err = -1 + 2*x - 3*x^2 + 4*x^3 - 3*x^5 + 5*x^6 - 6*x^7 + 6*x^8 - 3*x^10 + 4*x^11 - 4*x^12 + 4*x^13 - x^15 + x^16 - x^17 + x^18
round0(err)
 
 
 
### Derivation:
# y = (-x^5 + R)/b[1]
(x^5 - R)^5 - b[1]^6*x + R*b[1]^5
#
x^25 - 5*R*x^20 + 10*R^2*x^15 - 10*R^3*x^10 + 5*R^4*x^5 - b[1]^6*x + (R*b[1]^5 - R^5)
(x^5 + b[1]*x - R)*(x^20 - b[1]*x^16 - 4*R*x^15 + b[1]^2*x^12 + 3*R*b[1]*x^11 + 6*R^2*x^10 - b[1]^3*x^8 +
 - 2*R*b[1]^2*x^7 - 3*R^2*b[1]*x^6 - 4*R^3*x^5 + b[1]^4*x^4 + R*b[1]^3*x^3 +
 + R^2*b[1]^2*x^2 + R^3*b[1]*x + R^4 - b[1]^5)



###################################


