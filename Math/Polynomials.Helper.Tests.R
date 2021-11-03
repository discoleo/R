########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: Helper Functions
### Tests


### fast load:
source("Polynomials.Helper.R")


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


### p^2
mult.pm(p)

### p^3
p.v = pow.pm(p, 3)
p.v

print.pm(p.v[,c(2,3,4,1)])

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

# (x+1)*(x+2)*...*(x+6)
pR = toPoly.pm(paste("(x+", seq(1,6), ")", collapse="*"))
pR
eval.pm(pR, -1)
eval.pm(pR, -6)

p2 = toPoly.pm("(x+a+b)^3")
p2 = sort.pm(p2, "x", xn2= c("a", "b"))
p2

### x^3
toPoly.pm("p1(x = x-1)")

### x^3 * (x+b)^3
toPoly.pm("p1(x = x-1) * p2(x = x-a)")

###
f = function(p1) toPoly.pm("p1(x = x-1)")
f(p1)
# (x-1)^3
f(toPoly.pm("x^3"))


### Power n
n = 3
p = toPoly.pm("x^n + b*x - R")
p


###
# [resolved] clashes with internal variable;
n = 3
m = 2
p = toPoly.pm("x^(n+m) + b*x - R")
p

f = function(m) toPoly.pm("x^(n+m) + b*x - R")
m = -1
f(0)
f(2)
f(3)

########################
########################

### D
n = 3
p1 = toPoly.pm("x^n + b1*x + b0")
p2 = toPoly.pm("x^n + c2*x^2 + c1*x")
p = dp.exp.pm(list(Poly=p1, Exp=p2))
p$Poly = sort.pm(p$Poly, "x", xn2=c("c2","b1","c1"))
print.pm(p$Poly, leading="x")

# 3*x^5 + 2*c2*x^4 + c1*x^3 + 3*b1*x^3 + 3*x^2 + 3*b0*x^2 + 2*c2*b1*x^2 +
#	+ 2*c2*b0*x + c1*b1*x + c1*b0 + b1
