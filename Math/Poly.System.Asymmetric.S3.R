########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Asymmetric S3: Simple / Basic
###
### draft v.0.1b-fix


##########################
### Asymmetric Systems ###

### Simple Asymmetric Systems
### Basic Types

# - some simple Models;


####################

###############
### History ###
###############

### draft v.0.1b - v.0.1b-fix:
# - entanglement with roots of unity:
#   Order 1: x + m*y + m^2*z = 0;
#   Order 2: x^2 + m*y^2 + m^2*z^2 = 0;
# - improved versions; [v.0.1b-fix]
### draft v.0.1a:
# - initial draft;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R
# - e.g. round0(), round0.p(),
#   solve.EnAll(), solveEn();


########################

###############
### Order 2 ###
###############

# x^2 + b1*x + b2*S = R1
# y^2 + b1*y + b2*S = R2
# z^2 + b1*z + b2*S = R3

### Solution

### Sum =>
x^2 + y^2 + z^2 + (b1 + 3*b2)*S - (R1 + R2 + R3) # = 0
S^2 + (b1 + 3*b2)*S - 2*E2 - (R1 + R2 + R3) # = 0
# 2*E2 = S^2 + (b1 + 3*b2)*S - (R1 + R2 + R3)

### Prod(R[i] - b2*S)
(x^2 + b1*x)*(y^2 + b1*y)*(z^2 + b1*z) + (b2*S - R1)*(b2*S - R2)*(b2*S - R3) # = 0
E3^2 + b1*E3*E2 + b1^2*E3*S + b1^3*E3 + (b2*S - R1)*(b2*S - R2)*(b2*S - R3)
E3^2 + b1*E3*E2 + b1^2*E3*S + b1^3*E3 +
	+ (b2^3*S^3 - b2^2*(R1+R2+R3)*S^2 + b2*(R1*R2+R1*R3+R2*R3)*S - R1*R2*R3)

### Prod =>
(x^2 + b1*x + b2*S)*(y^2 + b1*y + b2*S)*(z^2 + b1*z + b2*S) - R1*R2*R3 # = 0
b1*E2*E3 + E3*b1^2*S + E3*b1^3 + E3^2 + b1*b2*(E2*S - 3*E3)*S + b1^2*b2*E2*S +
	+ b2*(E2^2 - 2*E3*S)*S + b2^2*(S^2 - 2*E2)*S^2 + b2^3*S^3 + b1*b2^2*S^3 - R1*R2*R3
b1*E2*E3 - 2*b2*E3*S^2 + b1^2*E3*S - 3*b1*b2*E3*S + b1^3*E3 + E3^2 +
	+ b1*b2*E2*S^2 - 2*b2^2*E2*S^2 + b1^2*b2*E2*S + b2*E2^2*S +
	+ b2^2*S^4 + b2^3*S^3 + b1*b2^2*S^3 - R1*R2*R3

### Rel:
SubstE3 = - 4*R1*R2*S*b2 - 4*R1*R3*S*b2 - 4*R1*S*b1^2*b2 - 8*R1*S^2*b1*b2 + 4*R1*S^2*b2^2 - 4*R1*S^3*b2 + 2*R1^2*S*b2 - 4*R2*R3*S*b2 - 4*R2*S*b1^2*b2 - 8*R2*S^2*b1*b2 + 4*R2*S^2*b2^2 - 4*R2*S^3*b2 + 2*R2^2*S*b2 - 4*R3*S*b1^2*b2 - 8*R3*S^2*b1*b2 + 4*R3*S^2*b2^2 - 4*R3*S^3*b2 + 2*R3^2*S*b2 + 12*S^2*b1^2*b2^2 + 4*S^2*b1^3*b2 + 24*S^3*b1*b2^2 + 10*S^3*b1^2*b2 - 6*S^3*b2^3 + 8*S^4*b1*b2 + 12*S^4*b2^2 + 2*S^5*b2;
DivE3 = - 24*S*b1*b2 - 16*S^2*b2;
E3 = - SubstE3 / DivE3;


### Eq: actual P[8]
# - seems NO win over direct solution;
# - embedded in S are a set of 8 real solutions and a larger set of false solutions:
#   the consequences of this are not yet fully understood;
(- 84*R1*R2*R3*b1^2 + 4*R1*R2*R3^2 + 8*R1*R2*b1^4 + 4*R1*R2^2*R3 + 10*R1*R2^2*b1^2 - 4*R1*R2^3 + 8*R1*R3*b1^4 + 10*R1*R3^2*b1^2 - 4*R1*R3^3 - 24*R1*b1^6 + 4*R1^2*R2*R3 + 10*R1^2*R2*b1^2 + 6*R1^2*R2^2 + 10*R1^2*R3*b1^2 + 6*R1^2*R3^2 + 28*R1^2*b1^4 - 4*R1^3*R2 - 4*R1^3*R3 - 10*R1^3*b1^2 + R1^4 + 8*R2*R3*b1^4 + 10*R2*R3^2*b1^2 - 4*R2*R3^3 - 24*R2*b1^6 + 10*R2^2*R3*b1^2 + 6*R2^2*R3^2 + 28*R2^2*b1^4 - 4*R2^3*R3 - 10*R2^3*b1^2 + R2^4 - 24*R3*b1^6 + 28*R3^2*b1^4 - 10*R3^3*b1^2 + R3^4)*S^2 +
(- 120*R1*R2*R3*b1 - 24*R1*R2*R3*b2 + 44*R1*R2*b1^2*b2 + 36*R1*R2*b1^3 + 12*R1*R2^2*b1 - 4*R1*R2^2*b2 + 44*R1*R3*b1^2*b2 + 36*R1*R3*b1^3 + 12*R1*R3^2*b1 - 4*R1*R3^2*b2 - 72*R1*b1^4*b2 - 120*R1*b1^5 + 12*R1^2*R2*b1 - 4*R1^2*R2*b2 + 12*R1^2*R3*b1 - 4*R1^2*R3*b2 + 10*R1^2*b1^2*b2 + 78*R1^2*b1^3 - 12*R1^3*b1 + 4*R1^3*b2 + 44*R2*R3*b1^2*b2 + 36*R2*R3*b1^3 + 12*R2*R3^2*b1 - 4*R2*R3^2*b2 - 72*R2*b1^4*b2 - 120*R2*b1^5 + 12*R2^2*R3*b1 - 4*R2^2*R3*b2 + 10*R2^2*b1^2*b2 + 78*R2^2*b1^3 - 12*R2^3*b1 + 4*R2^3*b2 - 72*R3*b1^4*b2 - 120*R3*b1^5 + 10*R3^2*b1^2*b2 + 78*R3^2*b1^3 - 12*R3^3*b1 + 4*R3^3*b2 + 72*b1^6*b2 + 24*b1^7)*S^3 +
(- 40*R1*R2*R3 + 72*R1*R2*b1*b2 + 48*R1*R2*b1^2 + 20*R1*R2*b2^2 + 4*R1*R2^2 + 72*R1*R3*b1*b2 + 48*R1*R3*b1^2 + 20*R1*R3*b2^2 + 4*R1*R3^2 - 54*R1*b1^2*b2^2 - 228*R1*b1^3*b2 - 238*R1*b1^4 + 4*R1^2*R2 + 4*R1^2*R3 + 12*R1^2*b1*b2 + 80*R1^2*b1^2 - 2*R1^2*b2^2 - 4*R1^3 + 72*R2*R3*b1*b2 + 48*R2*R3*b1^2 + 20*R2*R3*b2^2 + 4*R2*R3^2 - 54*R2*b1^2*b2^2 - 228*R2*b1^3*b2 - 238*R2*b1^4 + 4*R2^2*R3 + 12*R2^2*b1*b2 + 80*R2^2*b1^2 - 2*R2^2*b2^2 - 4*R2^3 - 54*R3*b1^2*b2^2 - 228*R3*b1^3*b2 - 238*R3*b1^4 + 12*R3^2*b1*b2 + 80*R3^2*b1^2 - 2*R3^2*b2^2 - 4*R3^3 + 108*b1^4*b2^2 + 360*b1^5*b2 + 116*b1^6)*S^4 +
(24*R1*R2*b1 + 24*R1*R2*b2 + 24*R1*R3*b1 + 24*R1*R3*b2 - 84*R1*b1*b2^2 - 256*R1*b1^2*b2 - 240*R1*b1^3 - 12*R1*b2^3 + 36*R1^2*b1 + 4*R1^2*b2 + 24*R2*R3*b1 + 24*R2*R3*b2 - 84*R2*b1*b2^2 - 256*R2*b1^2*b2 - 240*R2*b1^3 - 12*R2*b2^3 + 36*R2^2*b1 + 4*R2^2*b2 - 84*R3*b1*b2^2 - 256*R3*b1^2*b2 - 240*R3*b1^3 - 12*R3*b2^3 + 36*R3^2*b1 + 4*R3^2*b2 + 54*b1^2*b2^3 + 342*b1^3*b2^2 + 714*b1^4*b2 + 234*b1^5)*S^5 +
(4*R1*R2 + 4*R1*R3 - 120*R1*b1*b2 - 130*R1*b1^2 - 28*R1*b2^2 + 6*R1^2 + 4*R2*R3 - 120*R2*b1*b2 - 130*R2*b1^2 - 28*R2*b2^2 + 6*R2^2 - 120*R3*b1*b2 - 130*R3*b1^2 - 28*R3*b2^2 + 6*R3^2 + 84*b1*b2^3 + 384*b1^2*b2^2 + 720*b1^3*b2 + 255*b1^4 + 9*b2^4)*S^6 +
(- 36*R1*b1 - 20*R1*b2 - 36*R2*b1 - 20*R2*b2 - 36*R3*b1 - 20*R3*b2 + 180*b1*b2^2 + 390*b1^2*b2 + 162*b1^3 + 28*b2^3)*S^7 +
(- 4*R1 - 4*R2 - 4*R3 + 108*b1*b2 + 60*b1^2 + 30*b2^2)*S^8 +
(12*b1 + 12*b2)*S^9 +
(1)*S^10


### Example:

R = c(1,2,3)
b = c(2,3)
R1 = R[1]; R2 = R[2]; R3 = R[3];
b1 = b[1]; b2 = b[2];
#
coeff = c(
	(- 84*R1*R2*R3*b1^2 + 4*R1*R2*R3^2 + 8*R1*R2*b1^4 + 4*R1*R2^2*R3 + 10*R1*R2^2*b1^2 - 4*R1*R2^3 +
		8*R1*R3*b1^4 + 10*R1*R3^2*b1^2 - 4*R1*R3^3 - 24*R1*b1^6 + 4*R1^2*R2*R3 + 10*R1^2*R2*b1^2 +
		6*R1^2*R2^2 + 10*R1^2*R3*b1^2 + 6*R1^2*R3^2 + 28*R1^2*b1^4 - 4*R1^3*R2 - 4*R1^3*R3 - 10*R1^3*b1^2 +
		R1^4 + 8*R2*R3*b1^4 + 10*R2*R3^2*b1^2 - 4*R2*R3^3 - 24*R2*b1^6 + 10*R2^2*R3*b1^2 + 6*R2^2*R3^2 +
		28*R2^2*b1^4 - 4*R2^3*R3 - 10*R2^3*b1^2 + R2^4 - 24*R3*b1^6 + 28*R3^2*b1^4 - 10*R3^3*b1^2 + R3^4),
	(- 120*R1*R2*R3*b1 - 24*R1*R2*R3*b2 + 44*R1*R2*b1^2*b2 + 36*R1*R2*b1^3 + 12*R1*R2^2*b1 - 4*R1*R2^2*b2 +
		44*R1*R3*b1^2*b2 + 36*R1*R3*b1^3 + 12*R1*R3^2*b1 - 4*R1*R3^2*b2 - 72*R1*b1^4*b2 - 120*R1*b1^5 +
		12*R1^2*R2*b1 - 4*R1^2*R2*b2 + 12*R1^2*R3*b1 - 4*R1^2*R3*b2 + 10*R1^2*b1^2*b2 + 78*R1^2*b1^3 +
		- 12*R1^3*b1 + 4*R1^3*b2 + 44*R2*R3*b1^2*b2 + 36*R2*R3*b1^3 + 12*R2*R3^2*b1 - 4*R2*R3^2*b2 +
		- 72*R2*b1^4*b2 - 120*R2*b1^5 + 12*R2^2*R3*b1 - 4*R2^2*R3*b2 + 10*R2^2*b1^2*b2 + 78*R2^2*b1^3 +
		- 12*R2^3*b1 + 4*R2^3*b2 - 72*R3*b1^4*b2 - 120*R3*b1^5 + 10*R3^2*b1^2*b2 + 78*R3^2*b1^3 +
		- 12*R3^3*b1 + 4*R3^3*b2 + 72*b1^6*b2 + 24*b1^7), # *S +
	(- 40*R1*R2*R3 + 72*R1*R2*b1*b2 + 48*R1*R2*b1^2 + 20*R1*R2*b2^2 + 4*R1*R2^2 + 72*R1*R3*b1*b2 +
		48*R1*R3*b1^2 + 20*R1*R3*b2^2 + 4*R1*R3^2 - 54*R1*b1^2*b2^2 - 228*R1*b1^3*b2 - 238*R1*b1^4 +
		4*R1^2*R2 + 4*R1^2*R3 + 12*R1^2*b1*b2 + 80*R1^2*b1^2 - 2*R1^2*b2^2 - 4*R1^3 + 72*R2*R3*b1*b2 +
		48*R2*R3*b1^2 + 20*R2*R3*b2^2 + 4*R2*R3^2 - 54*R2*b1^2*b2^2 - 228*R2*b1^3*b2 - 238*R2*b1^4 +
		4*R2^2*R3 + 12*R2^2*b1*b2 + 80*R2^2*b1^2 - 2*R2^2*b2^2 - 4*R2^3 - 54*R3*b1^2*b2^2 - 228*R3*b1^3*b2 +
		- 238*R3*b1^4 + 12*R3^2*b1*b2 + 80*R3^2*b1^2 - 2*R3^2*b2^2 - 4*R3^3 + 108*b1^4*b2^2 +
		360*b1^5*b2 + 116*b1^6), # *S^2 +
	(24*R1*R2*b1 + 24*R1*R2*b2 + 24*R1*R3*b1 + 24*R1*R3*b2 - 84*R1*b1*b2^2 - 256*R1*b1^2*b2 - 240*R1*b1^3 +
		- 12*R1*b2^3 + 36*R1^2*b1 + 4*R1^2*b2 + 24*R2*R3*b1 + 24*R2*R3*b2 - 84*R2*b1*b2^2 - 256*R2*b1^2*b2 +
		- 240*R2*b1^3 - 12*R2*b2^3 + 36*R2^2*b1 + 4*R2^2*b2 - 84*R3*b1*b2^2 - 256*R3*b1^2*b2 - 240*R3*b1^3 +
		- 12*R3*b2^3 + 36*R3^2*b1 + 4*R3^2*b2 + 54*b1^2*b2^3 + 342*b1^3*b2^2 + 714*b1^4*b2 + 234*b1^5), # *S^3 +
	(6*R1^2 + 6*R2^2 + 6*R3^2 + 4*R1*R2 + 4*R1*R3 + 4*R2*R3 - 120*b1*b2*(R1 + R2 + R3) +
		 - 130*b1^2*(R1 + R2 + R3) - 28*b2^2*(R1 + R2 + R3) +
		84*b1*b2^3 + 384*b1^2*b2^2 + 720*b1^3*b2 + 255*b1^4 + 9*b2^4), # *S^4
	(-(36*b1 + 20*b2)*(R1 + R2 + R3) + 180*b1*b2^2 + 390*b1^2*b2 +
		162*b1^3 + 28*b2^3), # *S^5
	(- 4*(R1 + R2 + R3) + 108*b1*b2 + 60*b1^2 + 30*b2^2), # *S^6
	(12*b1 + 12*b2), 1 # S^8
)
S = roots(rev(coeff))
print(S)
E2 = (S^2 + (b1 + 3*b2)*S - (R1 + R2 + R3)) / 2;
SubstE3 = - 4*R1*R2*S*b2 - 4*R1*R3*S*b2 - 4*R1*S*b1^2*b2 - 8*R1*S^2*b1*b2 + 4*R1*S^2*b2^2 - 4*R1*S^3*b2 + 2*R1^2*S*b2 - 4*R2*R3*S*b2 - 4*R2*S*b1^2*b2 - 8*R2*S^2*b1*b2 + 4*R2*S^2*b2^2 - 4*R2*S^3*b2 + 2*R2^2*S*b2 - 4*R3*S*b1^2*b2 - 8*R3*S^2*b1*b2 + 4*R3*S^2*b2^2 - 4*R3*S^3*b2 + 2*R3^2*S*b2 + 12*S^2*b1^2*b2^2 + 4*S^2*b1^3*b2 + 24*S^3*b1*b2^2 + 10*S^3*b1^2*b2 - 6*S^3*b2^3 + 8*S^4*b1*b2 + 12*S^4*b2^2 + 2*S^5*b2;
DivE3 = - 24*S*b1*b2 - 16*S^2*b2;
E3 = - SubstE3 / DivE3;
### TODO: robust & all roots;
### x
# - unstable & fails!
# x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
# len = length(S)
# S = matrix(S, ncol=len, nrow=3, byrow=T)
# - fails to compute many of the roots;
x = sapply(seq_along(S), function(id) roots(c(1, b1, b2*S[id] - R1)))
### y, z
y = sapply(seq_along(S), function(id) roots(c(1, b1, b2*S[id] - R2)))
z = sapply(seq_along(S), function(id) roots(c(1, b1, b2*S[id] - R3)))
len = length(S)
S = matrix(S, ncol=len, nrow=2, byrow=T)
# x = matrix(x, ncol=len, nrow=2, byrow=T)
### Alternative:
x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
sol = solve.EnAll(x, n=3);
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 + b1*x + b2*(x+y+z) # - R1
y^2 + b1*y + b2*(x+y+z) # - R2
z^2 + b1*z + b2*(x+y+z) # - R3

isTrue1 = round0(x^2 + b1*x + b2*(x+y+z) - R1) == 0
isTrue2 = round0(y^2 + b1*y + b2*(x+y+z) - R2) == 0
isTrue3 = round0(z^2 + b1*z + b2*(x+y+z) - R3) == 0
isTrue = apply(cbind(isTrue1, isTrue2, isTrue3), 1, all)
(x+y+z)[isTrue] # every S-root;
sol[isTrue,]

poly.calc(sol[isTrue,1])
poly.calc(sol[isTrue,2])
poly.calc(sol[isTrue,3])


### Debug
b1 = -2; b2 = 3;
x = sqrt(2); y = sqrt(3); z = sqrt(5);
S = x + y + z; E2 = x*y + x*z + y*z; E3 = x*y*z;
R1 = x^2 + b1*x + b2*S;
R2 = y^2 + b1*y + b2*S;
R3 = z^2 + b1*z + b2*S;


### Classic Polynomial:
(- 4*R1*R2*b1*b2^3 - 4*R1*R3*b1*b2^3 + 4*R1*b1^2*b2^4 + 2*R1*b1^3*b2^3 - 2*R1^2*R2*b2^2 - 2*R1^2*R3*b2^2 + 8*R1^2*b1*b2^3 + 5*R1^2*b1^2*b2^2 + 4*R1^3*b1*b2 + 4*R1^3*b2^2 + R1^4 - 2*R2*R3*b2^4 - 2*R2*b1^2*b2^4 + R2^2*b2^4 - 2*R3*b1^2*b2^4 + R3^2*b2^4) +
(4*R1*R2*b1*b2^2 + 4*R1*R2*b2^3 + 4*R1*R3*b1*b2^2 + 4*R1*R3*b2^3 - 8*R1*b1*b2^4 - 26*R1*b1^2*b2^3 - 10*R1*b1^3*b2^2 - 24*R1^2*b1*b2^2 - 12*R1^2*b1^2*b2 - 8*R1^2*b2^3 - 4*R1^3*b1 - 4*R1^3*b2 + 4*R2*b1*b2^4 + 4*R2*b1^2*b2^3 + 4*R3*b1*b2^4 + 4*R3*b1^2*b2^3 - 6*b1^3*b2^4 - 2*b1^4*b2^3)*x^1 +
(4*R1*R2*b2^2 + 4*R1*R3*b2^2 + 12*R1*b1*b2^3 + 26*R1*b1^2*b2^2 + 12*R1*b1^3*b2 + 4*R1*b2^4 + 6*R1^2*b1^2 - 6*R1^2*b2^2 - 4*R1^3 - 2*R2*b1^2*b2^2 - 2*R2*b2^4 - 2*R3*b1^2*b2^2 - 2*R3*b2^4 + 9*b1^2*b2^4 + 16*b1^3*b2^3 + 5*b1^4*b2^2)*x^2 +
(36*R1*b1*b2^2 + 12*R1*b1^2*b2 - 4*R1*b1^3 + 12*R1*b2^3 + 12*R1^2*b1 + 12*R1^2*b2 - 4*R2*b1*b2^2 - 4*R2*b2^3 - 4*R3*b1*b2^2 - 4*R3*b2^3 + 6*b1^2*b2^3 - 6*b1^3*b2^2 - 4*b1^4*b2)*x^3 +
(- 12*R1*b1*b2 - 12*R1*b1^2 + 6*R1^2 - 2*R2*b2^2 - 2*R3*b2^2 - 16*b1*b2^3 - 25*b1^2*b2^2 - 8*b1^3*b2 + b1^4 - 3*b2^4)*x^4 +
(- 12*R1*b1 - 12*R1*b2 - 12*b1*b2^2 + 4*b1^3 - 4*b2^3)*x^5 +
(- 4*R1 + (b1 + b2)*(6*b1 + 2*b2))*x^6 +
4*(b1 + b2)*x^7 +
(1)*x^8


##########################
##########################

###################
### Curiosities ###
###################

### Omega-Entanglements

# - entanglement with roots of unity:
#   m^3 = 1;
# - Eqs 2 & 3: can be anything symmetric;

###############
### Order 1 ###
###############

# x + m*y + m^2*z = 0
# x*y + x*z + y*z = R2
# x*y*z = R3

### Solution:

### Eq 1 =>
(x + m*y + m^2*z)*(m^2*x + m*y + z) # = 0
m^2*(x^2 + y^2 + z^2) + (m+1)*(x*y + x*z + y*z) # = 0
m^2*(S^2 - 2*E2) + (m+1)*E2 # = 0
m^2*S^2 + (m+1 - 2*m^2)*E2 # = 0
m^2*S^2 - 3*m^2*E2 # = 0
S^2 - 3*E2 # = 0
### Eq:
# S^2 = 3*E2

### Solver:
solve.omega.P1 = function(R, debug=TRUE) {
	S = sqrt(3*R[1])
	S = c(S, -S);
	if(debug) print(S);
	x = sapply(S, function(S) roots(c(1, -S, R[1], -R[2])))
	# TODO: still NOT robust!
	m = unity(3, all=F);
	yz.s = S - x;
	yz.ms = - x;
	y = (m^2*yz.s - yz.ms) / (m^2 - m);
	z = (yz.s - y)
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	return(sol);
}
test.omega.P1 = function(sol, R=0, pow=1) {
	m = unity(3, all=F);
	x = sol[,1]; y = sol[,2]; z = sol[,3];
	err1 = x^pow + m*y^pow + m^2*z^pow # = 0
	err2 = x*y + x*z + y*z # - R2
	err3 = x*y*z # - R3
	err = cbind(err1, err2, err3)
	err = round0(err)
	err
}

### Examples

R = c(1, -1)
sol = solve.omega.P1(R)

### Test
test.omega.P1(sol, R);

x = sol[,1]
round0.p(poly.calc(sol[ ,1]))


###############

###############
### Order 2 ###
###############

# x^2 + m*y^2 + m^2*z^2 = 0
# x*y + x*z + y*z = R2
# x*y*z = R3

### Solution:

### Eq 1 =>
(x^2 + m*y^2 + m^2*z^2)*(m^2*x^2 + m*y^2 + z^2) # = 0
m^2*(x^4 + y^4 + z^4) + (m+1)*((x*y)^2 + (x*z)^2 + (y*z)^2) # = 0
m^2*(S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2) + (m+1)*(E2^2 - 2*E3*S) # = 0
m^2*S^4 - 4*E2*m^2*S^2 + 4*m^2*E3*S + 2*m^2*E2^2 + (m+1)*E2^2 - 2*(m+1)*E3*S # = 0
m^2*S^4 - 4*E2*m^2*S^2 - 2*(m+1 - 2*m^2)*E3*S + (2*m^2 + m + 1)*E2^2 # = 0
m^2*S^4 - 4*E2*m^2*S^2 + 6*m^2*E3*S + m^2*E2^2 # = 0
### Eq:
S^4 - 4*E2*S^2 + 6*E3*S + E2^2 # = 0

### Solver:
solve.omega.P2 = function(R) {
	coeff = c(1, 0, - 4*R[1], 6*R[2], R[1]^2)
	S = roots(coeff)
	x = sapply(S, function(S) roots(c(1, -S, R[1], -R[2])))
	R2 = R[1]; # !!
	m = unity(3, all=F);
	# TODO: robust
	# sol = solve.EnAll(x, n=3, max.perm=1)
	yz.s = S - x;
	yz = R2 - x*yz.s;
	yz.ms = - x^2;
	y2 = ((yz.s^2 - 2*yz)*m^2 - yz.ms) / (m^2 - m);
	y = sqrt(y2 + 0i);
	y = c(y, -y); yz.s = c(yz.s, yz.s); x = c(x, x)
	z = (yz.s - y);
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	isOK = round0(sol[,2]*sol[,3] - as.vector(c(yz, yz))) == 0;
	sol = sol[isOK,] # still needed;
	return(sol);
}

### Examples

R = c(1, -1)
sol = solve.omega.P2(R)

### Test
test.omega.P1(sol, R, pow=2);

x = sol[,1]
round0.p(poly.calc(sol))


##################

#############
### Variants:

# x + m*y + m^2*z = 0
# x^2 + y^2 + z^2 = R2
# x*y*z = R3

### Solution:

### Eq 1 =>
# S^2 = 3*E2

### Eq 2 =>
S^2 - 2*E2 - R2 # = 0
E2 - R2 # = 0
# E2 = R2

# - the same as the simple Order 1 version;


