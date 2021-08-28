########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: Helper Functions
### Tests


### fast load:
# source("Polynomials.Helper.R")


### this file:
# source("Polynomials.Helper.Tests.R")


#######################
#######################

#############
### Tests ###
#############

### Multi-variable Multiplication

# (x^3 + b1*x - R)^3
pTest = data.frame(
	x = c(3,1,0),
	b1 = c(0,1,0),
	R = c(0,0,1),
	coeff = c(1,1,-1)
)
p = toPoly.pm("x^3 + b1*x - R")
diff.pm(p, pTest)


### Test
mult.pm(p)

p.v = pow.pm(p, 3)
p.v

print.p(p.v[,c(2,3,4,1)])

### eval
R = 2; b1 = 3; x = -5;
#
eval.pm(p.v, c(R, x, b1))
(x^3 + b1*x - R)^3


###################
###################

### Advanced Parser

p1 = toPoly.pm("(x+1)^3")
p1

toPoly.pm(paste("(x+", seq(1,6), ")", collapse="*"))

p2 = toPoly.pm("(x+a+b)^3")
p2

toPoly.pm("p1(x = x-1) * p2(x = x-a)")
# x^3 * (x+b)^3


### Power n
n = 3
p = toPoly.pm("x^n + b*x - R")
p


###
n = 3
m = 2 # "m" clashes with an internal variable;
p = toPoly.pm("x^(n+m) + b*x - R")
p
