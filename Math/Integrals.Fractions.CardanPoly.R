
### Leonard Mada
###
### Integrals: Polynomial Fractions
### Cardan-Type Polynomials
###
### draft 0.1


############

### History

# - this is the nicer decomposition, using conjugate roots;
# - a previous solution used all the individual order 1 polynomials:
#  -- the old approach yields a compact solution as well;
#  -- but the current approach seems better;


############


### Examples
# 1/(x^5 - 5*c*x^3 + 5*c^2*x - 2*d)

n = 5 # b0 is currently limited to n = 5!
# Parameters: free to change
c = 1
d = 3
# Roots of unity
m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
n.half = (n-1)/2
m.m = matrix(m^(1:n.half), ncol=1)
m.m = cbind(m.m, 1/m.m)
m.all = m^(0:(n-1))
m.sum = apply(m.m, 1, sum)
# Roots
det = sqrt(d^2 - c^n)
p = (d + det)^(1/n)
q = (d - det)^(1/n)
r = p*m.all + q/m.all
r
# Coefficents of Partial Fractions
b0 = 1 / (1/(r[1]*r[1]) + 2/(r[2]*r[5]) + 2/(r[3]*r[4]) ) / r[1] / (2*d)
b = -2 * b0 * r[1] # ALL b are the same;
a = -b0 / m.sum[n.half:1]
c(a, b0, b) # displays only 1 coeff. b
# Test
x = 3 # any value - for testing the fraction;
#
1/(x^5 - 5*c*x^3 + 5*c^2*x - 2*d) # ==
b0/(x-p-q) + sum( (a*x + b) / ((x - p*m.m[,1] - q*m.m[,2]) * (x - p*m.m[,2] - q*m.m[,1])) )



#################
### integrals ###

### TODO: trivial



##################
##################

##################
### Derivation ###

# some of the steps used in the derivation:


R = p + q
R11 = (m*p + m*q + m^4*p + m^4*q) # = R*(m + m^4)
R12 = (m^2*p + m^2*q + m^3*p + m^3*q) # = R*(m^2 + m^3)
R21 = (m*p^2 + m*q^2 + m^4*p^2 + m^4*q^2) # = (R^2 - 2*c)*(m + m^4)
R22 = (m^2*p^2 + m^2*q^2 + m^3*p^2 + m^3*q^2) # = (R^2 - 2*c)*(m^2 + m^3)

# =>

a1 + a2 + b0 # = 0
a1*R11 + a2*R12 + b0*R + b1 + b2 # = 0
- a1*(R21 + c - c*m^2 - c*m^3) - a2*(R22 + c - c*m - c*m^4) + b0*(R^2 - 5*c) + b1*R11 + b2*R12 # = 0
a1*(c*R12 - R^3 + 3*c*R) + a2*(c*R11 - R^3 + 3*c*R) + b0*(R^3 - 5*c*R) - b1*(R21 + c - c*m^2 - c*m^3) - b2*(R22 + c - c*m - c*m^4) # = 0
b0 * (R^4 - 5*c*R^2 + 5*c^2) + b1 * (c*R12 - R^3 + 3*c*R) + b2 * (c*R11 - R^3 + 3*c*R) # = 1


### Solution to Linear System:

c.m = matrix(
c(1,1,1,0,0,
R11, R12, R, 1, 1,
- (R21 + c - c*m^2 - c*m^3), - (R22 + c - c*m - c*m^4), (R^2 - 5*c), R11, R12,
(c*R12 - R^3 + 3*c*R), (c*R11 - R^3 + 3*c*R), (R^3 - 5*c*R), - (R21 + c - c*m^2 - c*m^3), - (R22 + c - c*m - c*m^4),
0, 0, (R^4 - 5*c*R^2 + 5*c^2), (c*R12 - R^3 + 3*c*R), (c*R11 - R^3 + 3*c*R)
), ncol=5, byrow=T)

c.m
c.sol = solve(c.m, c(rep(0,4), 1))
c.sol

###
a = c.sol[1:2]
b0 = c.sol[3]
b = c.sol[4:5]

# the parametric solution can be derived as well;
# [see section: Examples]
