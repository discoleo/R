########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Polynomial Generators: Tests
###
### draft v.0.1a


### Tests for Polynomial Generators


####################
####################

### helper functions

source("Polynomials.Helper.Generators.R")


### fast load:
# source("Polynomials.Helper.Generators.Test.R")


######################
######################

### Tests

### Class 1 Poly:
b = c(0,1,-1,0,2)
K = 2;
p = toPoly.Class1S.pm(b, 5)
r = sum((K^(1/5))^seq(4) * b[-1])

print.p(p, "x")
eval.pm(p, c(K, r));


#################

### Class 2 Poly:

### P[4]
n = 5
m = unity(n, all=FALSE)
p = toPoly.Class2.pm(n-1)

### Ex 1:
s = c(1,2,3)
p2 = replace.pm(p, s, paste0("s", 0:2))
print.p(p2, "x")

x = sum(s * m^(0:2))
round0(eval.pm(p2, x))
x^4 + x^3 + 16*x^2 - 4*x + 41

### Ex 2:
s = c(1,-3,3)
p2 = replace.pm(p, s, paste0("s", 0:2))
print.p(p2, "x")

x = sum(s * m^(0:2))
round0(eval.pm(p2, x))


### P[6]
n = 7
m = unity(n, all=FALSE)
pow = c(1,2,5)
p = toPoly.Class2.pm(n-1, s.id=pow)

### Ex 1:
s = c(1,-2,3)
p2 = replace.pm(p, s, paste0("s", pow))
print.p(p2, "x")

x = sum(s * m^pow)
round0(eval.pm(p2, x))


##################

### Class 3

n = 3
m = 2*cos(2*pi/(2*n+1) * seq(n));
p = toPoly.Class3.pm(n)
p

### Ex 1:
s = c(1,-2,3)
p2 = replace.pm(p, s, paste0("s", 0:2))
print.p(p2, "x")

x = sum(s[1], s[-1] * m[-3])
round0(eval.pm(p2, x))

##################

### S2 Order 4:

p1 = toPoly.pm("x^4 + c2*x^2*y^2 + c1*x*y + b2*y^2+b1*y+b0")
p2 = permute.pm(p1)
#
pR = solve.pm(p1, p2, xn="y")
str(pR)
# pR$Rez contains 376 monomials
max(pR$Rez$x)
# Order 22: slightly above the correct Order 16!
top.pm(pR$Rez, "x")
# ugly factorization: 3 monomials with leading power of x;
# c2^3*(c2^2 - 1)^2 * x^22;

