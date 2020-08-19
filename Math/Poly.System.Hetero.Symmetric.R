
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogenous Symmetric
###
### draft v.0.1b-sh


###############
### History ###

### draft v.0.1b - v.0.1b-sh:
# - added a basic xy-type: x^3 + x*y = R;
# - added also the shift;
# - TODO: parametric classic polynomial;
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
	return(list(r=sol, p1=p1, p=p))
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
round0.p(poly.calc(sol[,1]))
# back-shift
round0.p(poly.calc(sol[,1] - s))

### TODO:
# - classical polynomial: parametric;

###
b = 2; R = 1
p = sapply(-6:6, function(s) print(poly.htxy(b, R, s)$p))
# [] - []*x + b*x^2 + b*[]*x^3 - b*(s+b)*x^4 - b*x^5 + x^6
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


