
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogenous Symmetric
###
### draft v.0.1d
### & branch v.0.2a


### Heterogenous Symmetric Polynomial Systems

### 2 Variables:
# x^n + P(x, y) = R
# y^n + P(y, x) = R


###############
### History ###

### branch v.0.2a:
# - more work on systems with 3 variables:
#   "proper" implementation of: x[i]*2 + b*x[k];
# - TODO: robust removal of set of wrong solutions;
#   [or avoid getting superfluous solutions ???]
### branch v.0.2a-pre-a:
# - initial work on systems with 3 variables;
# - the simple cases are less rewarding;
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
#   when shifted root is shifted back;
#   [in general not identical to non-shifted root polynomials]
### draft v.0.1a:
# - initial version: basic heterogenous systems;



library(polynom)
library(pracma)


###############

###############
### Order 3 ###

### x^3 + b*y

# x^3 + b1*y = R
# y^3 + b1*x = R

# Diff =>
# x^3 - y^3 - b1*(x-y) = 0
# (x - y)*(x^2 + x*y + y^2 - b1) = 0
# => x = y *OR* x^2 + x*y + y^2 - b1 = 0;
# =>
# (x+y)^2 - x*y - b1 = 0
# x*y = (x+y)^2 - b1;

# Sum =>
# (x+y)^3 - 3*x*y*(x+y) + b1*(x+y) = 2*R
# s^3 - 3*(s^2 - b1)*s + b1*s - 2*R = 0
# s^3 - 2*b1*s + R = 0

### Step 2:
# Solve:
# x + y = s
# x*y = s^2 - b1

### Example
b = 3
R = 1
#
r.sum = roots(c(1,0, -2*b[1], R))
xy = r.sum^2 - b[1]
diff = sqrt(r.sum^2 - 4*xy + 0i)
x = (r.sum + diff)/2
y = (r.sum - diff)/2
sol = cbind(x, y)
sol # TODO: include also x = y cases

### Test
x^3 + b[1]*y
y^3 + b[1]*x
# Classic
err = x^6 - b[1]*x^4 - 2*R*x^3 + b[1]^2*x^2 + b[1]*R*x + R^2 - b[1]^3
round0(err)

### Classic
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

###################

### Shifted Roots

### (x - s)^3 + b*y

# (x - s)^3 + b1*y = R
# (y - s)^3 + b1*x = R

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
# Diff =>
# (x - s)^3 - (y - s)^3 - b1*(x-y) = 0
# (x - y)*(x^2 + x*y + y^2 - 3*s*(x+y) + 3*s^2 - b1) = 0
# => x = y *OR* x^2 + x*y + y^2 - 3*s*(x+y) + 3*s^2 - b1 = 0;
# =>
# (x+y)^2 - 3*s*(x+y) - x*y + 3*s^2 - b1 = 0
# x*y = (x+y)^2 - 3*s*(x+y) + 3*s^2 - b1;
# x*y = Z^2 - 3*s*Z + 3*s^2 - b1;

# Sum =>
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
r.sum = roots(c(1, - 6*s, - 2*(b[1]-6*s^2), - 8*s^3 + R + 3*s*b[1]))
xy = r.sum^2 - 3*s*r.sum + 3*s^2 - b[1];
r.diff = sqrt(r.sum^2 - 4*xy + 0i)
x = (r.sum + r.diff)/2
y = (r.sum - r.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol # TODO: include also x = y cases

sol = solve.htShift(b, R, shift=s)
x = sol[,1]; y = sol[,2];
sol

### Test
(x-s)^3 + b[1]*y
(y-s)^3 + b[1]*x

### TODO:
# - classic + polynomial P6;
poly.calc(sol[,1])
round0.p(poly.calc(sol[,1] - s))

###
b = 2; R = 1;
s = 3/2
#
sol = solve.htShift(b, R, shift=s)
x = sol[,1]; y = sol[,2];
x = sol[,1] - s; # with shift back!
-4 - 4*x + 4*x^2 + 4*x^3 - 2*x^4 + x^6

###
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


### Classic
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
########################

###############
### Order 4 ###

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
# TODO:
# b1*y = R - x^4
# (R - x^4)^4/b1^4 + b1*x - R = 0
# (R - x^4)^4 + b1^5*x - R*b1^4


########################
########################

###############
### Order 3 ###

### x^3 + b*x*y

# x^3 + b1*x*y = R
# y^3 + b1*x*y = R

# Diff =>
# x^3 - y^3 = 0
# (x-y)*(x^2 + x*y + y^2) = 0
# Case 2:
# x^2 + x*y + y^2 = 0
# x*y = (x + y)^2
# Sum =>
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

b = 5
#
sol = solve.htxy(b, 1) # R = 1;
x = sol[,1]; y = sol[,2];
sol
err = 1 + b[1]*x^2 - 2*x^3 + b[1]^2*x^4 - b[1]*x^5 + x^6
round0(err)


###################

### Shifted Roots

### (x - s)^3 + b*x*y

# (x-s)^3 + b1*x*y = R
# (y-s)^3 + b1*x*y = R

# Diff =>
# (x-s)^3 - (y-s)^3 = 0
# (x-y)*((x-s)^2 + (x-s)*(y-s) + (y-s)^2) = 0
# alternatively: x - s = (y-s)*m, where m^3 = 1;
# Case 2:
# (x-s)^2 + (x-s)*(y-s) + (y-s)^2 = 0
# (x+y - 2*s)^2 = (x-s)*(y-s)
# x*y = (x+y - 2*s)^2 + s*(x+y) - s^2
# x*y = (Z - 2*s)^2 + s*Z - s^2
# x*y = Z^2 - 3*s*Z + 3*s^2 ###
# Sum =>
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

# Test
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


################
################

##################
### xy & x/y-Terms

### x-Term
### x^3 + b1*x*y + b2*x

# x^3 + b1*x*y + b2*x = R
# y^3 + b1*x*y + b2*y = R

### Solution:
# Diff =>
# x^3 - y^3 + b2*(x - y) = 0
# (x-y)*(x^2 + x*y + y^2 + b2) = 0
# Case 2:
# x^2 + x*y + y^2 + b2 = 0
# x*y = (x + y)^2 + b2
# Sum =>
# x^3 + y^3 + 2*b1*x*y + b2*(x+y) - 2*R = 0
# (x+y)^3 - 3*x*y*(x+y) + 2*b1*x*y + b2*(x+y) - 2*R
# Z^3 - 3*x*y*Z + 2*b1*x*y + b2*Z - 2*R
# Z^3 - 3*Z^3 - 3*b2*Z + 2*b1*Z^2 + 2*b1*b2 + b2*Z - 2*R
# Z^3 - b1*Z^2 + b2*Z - b1*b2 + R

solve.htxy = function(b, R, isX=TRUE) {
	if(length(b) < 2) b = c(b, 0)
	if(isX) {
		x.sum = roots(c(1, - b[1], b[2], - b[1]*b[2] + R))
		xy = x.sum^2 + b[2]
	} else {
		x.sum = roots(c(1, - b[1], -2*b[2], b[1]*b[2] + R))
		xy = x.sum^2 - b[2]
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
x^3 + b[1]*x*y + b[2]*x 
y^3 + b[1]*x*y + b[2]*y

### Classical Polynomial
# TODO

round0.p(poly.calc(sol[,1]))

################

### y-Term
### x^3 + b1*x*y + b2*y

# x^3 + b1*x*y + b2*y = R
# y^3 + b1*x*y + b2*x = R

### Solution:
# Diff =>
# x^3 - y^3 - b2*(x - y) = 0
# (x-y)*(x^2 + x*y + y^2 - b2) = 0
# Case 2:
# x^2 + x*y + y^2 - b2 = 0
# x*y = (x + y)^2 - b2
# Sum =>
# x^3 + y^3 + 2*b1*x*y + b2*(x+y) - 2*R = 0
# (x+y)^3 - 3*x*y*(x+y) + 2*b1*x*y + b2*(x+y) - 2*R
# Z^3 - 3*x*y*Z + 2*b1*x*y + b2*Z - 2*R
# Z^3 - 3*Z^3 + 3*b2*Z + 2*b1*Z^2 - 2*b1*b2 + b2*Z - 2*R
# Z^3 - b1*Z^2 - 2*b2*Z + b1*b2 + R

### Example
b = c(-1, 1)
R = 2
#
sol = solve.htxy(b, R, isX=FALSE)
x = sol[,1]; y = sol[,2]
sol

### Test
x^3 + b[1]*x*y + b[2]*y
y^3 + b[1]*x*y + b[2]*x

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
x^3 + b[1]*x*y + b[2]*y
y^3 + b[1]*x*y + b[2]*x

err = 5 - 2*x^3 + 2*x^4 + x^5 + x^6
round0(err)


################

### TODO:
# - shift;



################
################

### (xy)^2 Term
### x^3 + b1*(x*y)^2

# x^3 + b1*(x*y)^2 = R
# y^3 + b1*(x*y)^2 = R

m3 = unity(3, all=F)

# Diff =>
# x^3 - y^3 = 0
# y = x*m, where m^3 = 1;
# Case 2:
# separate equations for: m & m^2
# b1*x^4*m^2 + x^3 - R = 0
# b1*x^4*m + x^3 - R = 0

solve.xysq = function(b, R, isInverse=FALSE) {
	coeff = if(isInverse) c(m3^2, b[1], 0,0, - R*b[1]) else c(b[1]*m3^2, 1,0,0, - R);
	x = roots(coeff)
	y = x*m
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


################
################

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

###########
### Order 2

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
