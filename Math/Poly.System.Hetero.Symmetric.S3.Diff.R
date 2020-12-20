
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Hetero-Symmetric Differences
###
### draft v.0.1b


### Hetero-Symmetric Differences
### Polynomial Systems: 3 Variables

### Example:
x^n - y^n + b*x*y = R
y^n - z^n + b*y*z = R
z^n - x^n + b*z*x = R

#####################

###############
### History ###

### draft v.0.1b:
# - solved: x^3 - y^3 + b1*x*y = R;
### draft v.0.1a-ext2:
# - extension A1 power 2: + b[3]*(x+y+z)^2;
### draft v.0.1a:
# - moved [Difference] section to this new file;
### old file: Poly.System.Hetero.Symmetric.S3.R
### draft v.0.2b - v.0.2b-ext:
# - first concepts / solved [v.0.2b-sol]:
#   x^2 - y^2 + b*x*y = R;
# - solved extension A1 (Pow 1): + b[2]*(x+y+z); [v.0.2b-ext]


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R

##########################


##########################
### Difference Types   ###

# quasi-[Negative Correlations]


###############
### Order 2 ###
###############

# x^2 - y^2 + b*x*y = R
# y^2 - z^2 + b*y*z = R
# z^2 - x^2 + b*x*z = R

### Solution:

### Case 1:
det = sqrt(R[1]/b[1])
x = y = z = c(det, -det)

### Case 2:
# (x, y, z) NOT equal;

### Sum =>
b*E2 - 3*R # = 0
E2 = 3*R[1] / b[1]

### Sum(z^2*...) =>
# b*x*y*z*S = R*(x^2+y^2+z^2)
b*E3*S - R*(S^2 - 2*E2) # = 0
b*E3*S = R*(S^2 - 2*E2)
b*E3*S = R*(S^2 - 6*R / b)

### Sum(x*y*...) =>
(x^3*y - x*y^3 - x^3*z + x*z^3 + y^3*z - y*z^3) + b*(x^2*y^2 + x^2*z^2 + y^2*z^2) - R*E2 # = 0
E2*(S^2 - 2*E2) - E3*S - 2*(x*y^3 + x^3*z + y*z^3) + b*(E2^2 - 2*E3*S) - R*E2 # = 0
2*(x*y^3 + x^3*z + y*z^3) - E2*S^2 + 2*b*E3*S + E3*S + (2-b)*E2^2 + R*E2 # = 0
2*(x*y^3 + x^3*z + y*z^3) - E2*S^2 + 2*(R*S^2 - 2*R*E2) + E3*S + (2-b)*E2^2 + R*E2 # = 0
2*b*(x*y^3 + x^3*z + y*z^3) - b*E2*S^2 + 2*b*R*S^2 - 4*R*b*E2 + b*E3*S + (2-b)*b*E2^2 + R*b*E2 # = 0
b^2*(x*y^3 + x^3*z + y*z^3) + b^2*R*S^2 - b*R*S^2 + 6*R^2 - 9*b*R^2 # = 0
### Eq 1: b^2*(x*y^3 + x^3*z + y*z^3) = - (b^2*R*S^2 - b*R*S^2 + 6*R^2 - 9*b*R^2)
### (x^3*y+y^3*z+z^3*x)* =>
2*b^2*(x^4*y^4+x^4*z^4+y^4*z^4 + E3*(x^5+y^5+z^5) + E3^2*E2) +
	+ (2*b^2*R*S^2 - 2*b*R*S^2 + 12*R^2 - 18*b*R^2)*(x^3*y+y^3*z+z^3*x) # = 0
2*b^2*(4*E2*E3^2 + E2^4 - 4*E2^2*E3*S + 2*E3^2*S^2 + E3*(S^5 - 5*E2*S^3 + 5*E3*S^2 + 5*E2^2*S - 5*E3*E2) + E3^2*E2) +
	+ (2*b^2*R*S^2 - 2*b*R*S^2 + 12*R^2 - 18*b*R^2)*(x^3*y+y^3*z+z^3*x) # = 0
2*b^2*(E3*S^5 - 5*E2*E3*S^3 + 7*E3^2*S^2 + E2^2*E3*S + E2^4) +
	+ (2*b^2*R*S^2 - 2*b*R*S^2 + 12*R^2 - 18*b*R^2)*(x^3*y+y^3*z+z^3*x) # = 0
# =>
2*b*(R*(S^2 - 2*E2)*S^4 - 5*E2*R*(S^2 - 2*E2)*S^2 + 7*R*(S^2 - 2*E2)*E3*S + E2^2*R*(S^2 - 2*E2) + b*E2^4) +
	+ (2*b^2*R*S^2 - 2*b*R*S^2 + 12*R^2 - 18*b*R^2)*(x^3*y+y^3*z+z^3*x) # = 0
2*(b*R*S^6 - 14*R^2*S^4 + 5*R^2*E2*S^2 + 31*R^2*E2^2) +
	+ (2*b^2*R*S^2 - 2*b*R*S^2 + 12*R^2 - 18*b*R^2)*(x^3*y+y^3*z+z^3*x) # = 0
(b^3*R*S^6 - 14*b^2*R^2*S^4 + 15*R^3*b*S^2 + 279*R^4) +
	+ b^2*(b^2*R*S^2 - b*R*S^2 + 6*R^2 - 9*b*R^2)*(x^3*y+y^3*z+z^3*x) # = 0
### (Eq 1) * (...) + Eq 2:
b^2*(b^2*R*S^2 - b*R*S^2 + 6*R^2 - 9*b*R^2)*(x*y^3 + x^3*z + y*z^3 + x^3*y + y^3*z + z^3*x) +
	+ (b^3*R*S^6 - 14*b^2*R^2*S^4 + 15*R^3*b*S^2 + 279*R^4) +
	+ (b^2*R*S^2 - b*R*S^2 + 6*R^2 - 9*b*R^2)*(b^2*R*S^2 - b*R*S^2 + 6*R^2 - 9*b*R^2)
b^2*(b^2*R*S^2 - b*R*S^2 + 6*R^2 - 9*b*R^2)*(E2*(S^2 - 2*E2) - E3*S) +
	+ (b^3*R*S^6 - 14*b^2*R^2*S^4 + 15*R^3*b*S^2 + 279*R^4) +
	+ (b^2*R*S^2 - b*R*S^2 + 6*R^2 - 9*b*R^2)*(b^2*R*S^2 - b*R*S^2 + 6*R^2 - 9*b*R^2)
(b^2*R*S^2 - b*R*S^2 + 6*R^2 - 9*b*R^2)*(2*b*R*S^2 - 12*R^2) +
	+ (b^3*R*S^6 - 14*b^2*R^2*S^4 + 15*R^3*b*S^2 + 279*R^4) +
	+ (b^2*R*S^2 - b*R*S^2 + 6*R^2 - 9*b*R^2)*(b^2*R*S^2 - b*R*S^2 + 6*R^2 - 9*b*R^2)
b^3*S^6 + R*(b^4 - 15*b^2)*S^4 + R^2*(27*b - 18*b^3)*S^2 + 81*R^3*(b^2 + 3)
(b*S^2 - 9*R)*(b^2*S^4 + R*(b^3 - 6*b)*S^2 - 9*R^2*(b^2 + 3))
(b*S^2 - 9*R)^2 * (b*S^2 + R*(b^2 + 3))

### Extension A1: pow 1:
(b[1]*S^2 + 9*b[2]*S - 9*R)^2 * (b[1]*S^2 - b[2]*(b[1]^2 + 3)*S + R*(b[1]^2 + 3))
### Extension A1: pow 2:
((b[1] + 9*b[3])*S^2 + 9*b[2]*S - 9*R)^2 * ((b[1] - b[3]*(b[1]^2 + 3))*S^2 - b[2]*(b[1]^2 + 3)*S + R*(b[1]^2 + 3))


### Solution

solve.Ht3Diff = function(R, b) {
	if(length(b) == 1) {
		coeff = c(b[1], 0, R[1]*(b[1]^2 + 3))
	} else if(length(b) == 2) {
		coeff = c(b[1], - b[2]*(b[1]^2 + 3), R[1]*(b[1]^2 + 3))
	} else if(length(b) == 3) {
		coeff = c((b[1] - b[3]*(b[1]^2 + 3)), - b[2]*(b[1]^2 + 3), R[1]*(b[1]^2 + 3))
	}
	det = sqrt(coeff[2]^2 - 4*coeff[1]*coeff[3] + 0i)
	S = c( - coeff[2] + det, - coeff[2] - det) / 2 / coeff[1]
	print(S)
	#
	b2 = if(length(b) > 1) b[2] else 0; # Ext A1: pow 1;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A3: pow 2;
	R1 = R[1] - b2*S - b3*S^2
	E2 = 3*R1 / b[1]
	E3 = R1*(S^2 - 2*E2) / b[1] / S
	### x
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], - E3[id])))
	len = length(S)
	S  = rep(S,  each=3)
	E3 = rep(E3, each=3)
	R1 = rep(R1, each=3)
	isZero = round0(x) == 0
	x = x[ ! isZero] # TODO: changes length of x => S;
	yz = E3/x
	yz.s = S - x
	### robust: includes Ext A1: powers 1 & 2;
	y2 = (yz.s^2 + R1 - (2 + b[1])*yz) / 2
	y = (R1 - x^2 + y2) / b[1] / x
	z = yz.s - y;
	sol = cbind(x=x, y=y, z=z)
	### x = 0
	if(any(isZero)) {
		print("TODO: x == 0")
	}
	### x == y == z
	if(length(b) < 2) {
		x = y = z = c(1,-1) * sqrt(R[1] / b[1] + 0i)
	} else if(length(b) == 2) {
		x = y = z = roots(c(b[1], 3*b[2], -R[1]))
	} else if(length(b) == 3) {
		x = y = z = roots(c(b[1] + 9*b[3], 3*b[2], -R[1]))
	}
	sol = rbind(sol, cbind(x,y,z))
	return(sol)
}

### Examples:

R = 1
b = 1
#
sol = solve.Ht3Diff(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 - y^2 + b[1]*x*y # - R
y^2 - z^2 + b[1]*y*z # - R
z^2 - x^2 + b[1]*x*z # - R

round0.p(poly.calc(x[1:6]))

#########
### Ex 2:
R = 1
b = 5
#
sol = solve.Ht3Diff(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 - y^2 + b[1]*x*y # - R
y^2 - z^2 + b[1]*y*z # - R
z^2 - x^2 + b[1]*x*z # - R

round0.p(poly.calc(x[1:6])) * 5

###############
### Extensions:

R = 1
b = c(1, -2)
#
sol = solve.Ht3Diff(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 - y^2 + b[1]*x*y + b[2]*(x+y+z) # - R
y^2 - z^2 + b[1]*y*z + b[2]*(x+y+z) # - R
z^2 - x^2 + b[1]*x*z + b[2]*(x+y+z) # - R

round0.p(poly.calc(x[1:6]))

err = 25 + 60*x - 131*x^2 - 284*x^3 - 38*x^4 + 8*x^5 + x^6
round0(err)

##########
### Ext 2:

R = 1
b = c(1, -2, 1)
#
sol = solve.Ht3Diff(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 - y^2 + b[1]*x*y + b[3]*(x+y+z)^2 + b[2]*(x+y+z) # - R
y^2 - z^2 + b[1]*y*z + b[3]*(x+y+z)^2 + b[2]*(x+y+z) # - R
z^2 - x^2 + b[1]*x*z + b[3]*(x+y+z)^2 + b[2]*(x+y+z) # - R

round0.p(poly.calc(x[1:6]))

err = 25 + 60*x - 131*x^2 - 284*x^3 - 38*x^4 + 8*x^5 + x^6
round0(err)


### Debug
x =  1.6510934088i
y = -1.2738905550i
z = -2.3772028539i
#
S = x+y+z
E3 = x*y*z
E2 = x*y + x*z + y*z

# [old][redundant]
### Sum((x^2+y^2)*...) =>
# b*(x^3*y + x*y^3 + x^3*z + x*z^3 + y^3*z + y*z^3) = 2*R*(x^2+y^2+z^2)
b*(E2*(S^2 - 2*E2) - E3*S) - 2*R*(S^2 - 2*E2) # = 0
b*E2*S^2 - 2*b*E2^2 - b*E3*S - 2*R*S^2 + 4*R*E2 # = 0
b*E2*S^2 - 2*b*E2^2 - R*(S^2 - 6*R / b) - 2*R*S^2 + 4*R*E2 # = 0
b^2*E2*S^2 - 2*b^2*E2^2 - R*(b*S^2 - 6*R) - 2*b*R*S^2 + 4*R*b*E2 # = 0
3*b*R*S^2 - 3*b*R*S^2 # = 0


###########################

###############
### Order 3 ###
###############

# x^3 - y^3 + b*x*y = R
# y^3 - z^3 + b*y*z = R
# z^3 - x^3 + b*x*z = R

### Solution:

### Case 1:
det = sqrt(R[1]/b[1])
x = y = z = c(det, -det)

### Case 2:
# (x, y, z) NOT equal;

### Sum =>
b*E2 - 3*R # = 0
E2 = 3*R[1] / b[1]

### Sum(z^3*...) =>
# b*x*y*z*(x^2 + y^2 + z^2) = R*(x^3 + y^3 + z^3)
b*E3*(S^2 - 2*E2) - R*(S^3 - 3*E2*S + 3*E3) # = 0
(b*S^2 - 2*b*E2 - 3*R)*E3 - R*(S^3 - 3*E2*S) # = 0
(b*S^2 - 9*R)*E3 - R*(S^3 - 3*E2*S) # = 0
b*(b*S^2 - 9*R)*E3 - R*(b*S^3 - 9*R*S) # = 0

### Eq 3:
# x^3 - y^3 = R - b[1]*x*y
### sq =>
x^6 + y^6 - 2*x^3*y^3 - (R^2 + b[1]^2*x^2*y^2 - 2*b[1]*R*x*y) # = 0
### Sum() =>
2*(x^6 + y^6 + z^6) - 2*(x^3*y^3 + x^3*z^3 + y^3*z^3) - b[1]^2*(x^2*y^2 + x^2*z^2 + y^2*z^2) + 2*b[1]*R*E2 - 3*R^2
2*(- 2*E2^3 + 3*E3^2 - 12*E2*E3*S + 9*E2^2*S^2 + 6*E3*S^3 - 6*E2*S^4 + S^6) +
	- 2*(E2^3 - 3*E3*(E2*S - E3)) - b[1]^2*(E2^2 - 2*E3*S) + 2*b[1]*R*E2 - 3*R^2
- 6*E2^3 - b[1]^2*E2^2 - 18*E2*E3*S + 18*E2^2*S^2 + 12*E3*S^3 - 12*E2*S^4 + 2*S^6 + 2*b[1]^2*E3*S + 2*b[1]*R*E2 - 3*R^2
- 6*b[1]^3*E2^3 - b[1]^5*E2^2 - 18*b[1]^3*E2*E3*S + 18*b[1]^3*E2^2*S^2 +
	+ 12*b[1]^3*E3*S^3 - 12*b[1]^3*E2*S^4 + 2*b[1]^3*S^6 + 2*b[1]^5*E3*S + 2*b[1]^4*R*E2 - 3*b[1]^3*R^2
(12*b[1]^3*S^3 + 2*b[1]^5*S - 54*b[1]^2*R*S)*E3 +
	+ 2*b[1]^3*S^6 - 36*b[1]^2*R*S^4 + 2*81*b[1]*R^2*S^2 - 6*27*R^3 - 6*b[1]^3*R^2
R*(12*b[1]^3*S^3 + 2*b[1]^5*S - 54*b[1]^2*R*S)*(b[1]*S^3 - 9*R*S) +
	+ b[1]*(b[1]*S^2 - 9*R)*(2*b[1]^3*S^6 - 36*b[1]^2*R*S^4 + 2*81*b[1]*R^2*S^2 - 6*27*R^3 - 6*b[1]^3*R^2)

### =>
(27*R^3*b1^4 + 729*R^4*b1) +
	+ (- 12*R^2*b1^5 - 567*R^3*b1^2)*S^2 +
	+ (R*b1^6 + 162*R^2*b1^3)*S^4 +
	+ (- 21*R*b1^4)*S^6 + b1^5*S^8
(b[1]*S^2 - 9*R) * (b[1]^4*S^6 - 12*R*b[1]^3*S^4 + (R*b[1]^5 + 54*R^2*b[1]^2)*S^2 - (3*R^2*b[1]^4 + 81*R^3*b[1]))
(b[1]*S^2 - 9*R) * (b[1]*S^2 - 3*R) * (b[1]^3*S^4 - 9*R*b[1]^2*S^2 + (R*b[1]^4 + 27*R^2*b[1]))


coeff = c(b[1]^3, 0, - 9*R[1]*b[1]^2, 0, (R[1]*b[1]^4 + 27*R[1]^2*b[1]))


### Solution
solve.Ht3DiffP3 = function(R, b) {
	if(b[1] == 0) stop("Cannot solve when b[1] == 0!")
	if(length(b) == 1) {
		coeff = c(b[1]^3, 0, - 9*R[1]*b[1]^2, 0, (R[1]*b[1]^4 + 27*R[1]^2*b[1]))
	} else if(length(b) == 2) {
		# TODO
	} else if(length(b) == 3) {
		# TODO
	}
	S = roots(coeff)
	print(S)
	#
	b2 = if(length(b) > 1) b[2] else 0; # Ext A1: pow 1;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A3: pow 2;
	R1 = R[1] - b2*S - b3*S^2
	E2 = 3*R1 / b[1]
	E3 = R1*(b[1]*S^3 - 9*R1*S) / b[1] / (b[1]*S^2 - 9*R1)
	### x
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], - E3[id])))
	len = length(S)
	S  = rep(S,  each=3)
	E3 = rep(E3, each=3)
	R1 = rep(R1, each=3)
	isZero = round0(x) == 0
	# x = x[ ! isZero] # TODO: changes length of x => S;
	yz = E3/x
	yz.s = S - x
	### robust: includes Ext A1: powers 1 & 2;
	y3 = (yz.s^3 + R1 - (3*yz.s + b[1])*yz) / 2
	y = (R1 - x^3 + y3) / b[1] / x
	z = yz.s - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	### x = 0
	if(any(isZero)) {
		print("TODO: x == 0")
	}
	### x == y == z
	if(length(b) < 2) {
		x = y = z = c(1,-1) * sqrt(R[1] / b[1] + 0i)
	} else if(length(b) == 2) {
		x = y = z = roots(c(b[1], 3*b[2], -R[1]))
	} else if(length(b) == 3) {
		x = y = z = roots(c(b[1] + 9*b[3], 3*b[2], -R[1]))
	}
	sol = rbind(sol, cbind(x,y,z))
	return(sol)
}


### Examples:
R = 2
b = 3
#
sol = solve.Ht3DiffP3(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^3 - y^3 + b[1]*x*y # - R[1]
y^3 - z^3 + b[1]*y*z # - R[1]
z^3 - x^3 + b[1]*x*z # - R[1]


### Debug
x = 1.2583507026 + 1.1639791181i
y = 0.2282819654 - 1.0195743279i
z = 0.4163440587 + 0.6438339512i
#
S = x+y+z
E3 = x*y*z
E2 = x*y + x*z + y*z


