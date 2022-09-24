########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
###  == Derivation ==
###  Type: L2 V1
###
### draft v.0.1a


####################
####################

### Helper Functions

source("Polynomials.Helper.R")


####################
####################

### History

# draft v.0.1a:
# - moved specific code from file:
#   Poly.System.Hetero.Symmetric.S4.Derivation.R; [stage: draft v.0.3b]


################################
################################

################################
### Mixed Leading Terms: L2  ###
### Simple Chain: Type V1    ###
################################

###############
### Order 2 ###
###############

### x1^2*x2^2 + b*x3 = R

# x1^2*x2^2 + b*x3 = R
# x2^2*x3^2 + b*x4 = R
# x3^2*x4^2 + b*x1 = R
# x4^2*x1^2 + b*x2 = R

# "Llamas are much bigger than frogs!"

### Solution:

### Case: all x[i] different;

### Eq 1:
# x1^2*x2^2 = R - b*x3 => Prod =>
### Eq 1:
E4^4 - b^4*E4 + b^3*R*E3 - b^2*R^2*E2 + b*R^3*S - R^4 # = 0

### Eq 2:
### x1^2*x2^2 = R - b*x3 => Prod(Eq1, Eq 3) =>
E4^2 - R^2 + b*R*(x1+x3) - b^2*x1*x3 # = 0
E4^2 - R^2 + b*R*(x2+x4) - b^2*x2*x4 # = 0
### Sum =>
2*E4^2 - b^2*(x1*x3 + x2*x4) + b*R*S - 2*R^2 # = 0
# b^2*(x1*x3 + x2*x4) = 2*E4^2 + b*R*S - 2*R^2

### Sum =>
((x1*x2)^2+(x2*x3)^2+(x3*x4)^2+(x4*x1)^2) + b*S - 4*R # = 0
E2^2 - 2*E3*S + 2*E4 + b*S - 4*R - (2*E4^2 + b*R*S - 2*R^2)^2 / b^4 + 2*E4 # = 0
b^4*E2^2 - 2*b^4*E3*S + 4*b^4*E4 + b^5*S - 4*b^4*R - (2*E4^2 + b*R*S - 2*R^2)^2 # = 0
b^4*E2^2 - 2*b^4*E3*S + 4*b^4*E4 + b^5*S - 4*b^4*R +
	- (4*E4^4 + b^2*R^2*S^2 + 4*R^4 + 4*b*R*E4^2*S - 8*R^2*E4^2 - 4*b*R^3*S) # = 0
### Eq 2:
4*E4^4 + 4*b*R*E4^2*S - 8*R^2*E4^2 - 4*b^4*E4 - b^4*E2^2 + 2*b^4*E3*S +
	+ b^2*R^2*S^2 - b^5*S - 4*b*R^3*S + 4*b^4*R + 4*R^4 # = 0
### Eq 2 (simplified)
# Diff: Eq 2 - 4*Eq 1 =>
4*b*R*E4^2*S - 8*R^2*E4^2 - 4*b^3*R*E3 - b^4*E2^2 + 4*b^2*R^2*E2 + 2*b^4*E3*S +
	- 4*b*R^3*S - b^5*S + b^2*R^2*S^2 - 4*b*R^3*S + 4*b^4*R + 8*R^4 # = 0

### Eq 3:
### Sum(x3^2*...) =>
((x1*x2*x3)^2+(x2*x3*x4)^2+(x3*x4*x1)^2+(x4*x1*x2)^2) + b*(S^3 - 3*E2*S + 3*E3) - R*(S^2 - 2*E2) # = 0
(E3^2 - 2*E2*E4) + b*(S^3 - 3*E2*S + 3*E3) - R*(S^2 - 2*E2) # = 0
### Eq 3:
E3^2 - 2*E2*E4 + 3*b*E3 - 3*b*E2*S + 2*R*E2 + b*S^3 - R*S^2 # = 0

### Eq 4:
### Sum(x3*x4*...) =>
E4*(E2 - (x1*x3+x2*x4)) + b*(x1^2*x2+x2^2*x3+x3^2*x4+x4^2*x1) - R*(E2 - (x1*x3+x2*x4)) # = 0
b^2*E2*E4 - E4*b^2*(x1*x3+x2*x4) + R*b^2*(x1*x3+x2*x4) - b^2*R*E2 + b^3*(x1^2*x2+x2^2*x3+x3^2*x4+x4^2*x1) # = 0
b^2*E2*E4 - E4*(2*E4^2 + b*R*S - 2*R^2) + R*(2*E4^2 + b*R*S - 2*R^2) - b^2*R*E2 +
	+ b^3*(x1^2*x2+x2^2*x3+x3^2*x4+x4^2*x1) # = 0
2*E4^3 - 2*R*E4^2 - b^2*E2*E4 + b*R*E4*S - 2*R^2*E4 + b^2*R*E2 - b*R^2*S + 2*R^3 +
	- b^3*(x1^2*x2+x2^2*x3+x3^2*x4+x4^2*x1) # = 0
# b^3*(x1^2*x2+x2^2*x3+x3^2*x4+x4^2*x1) =
2*E4^3 - 2*R*E4^2 - b^2*E2*E4 + b*R*E4*S - 2*R^2*E4 + b^2*R*E2 - b*R^2*S + 2*R^3
# =>
b^3*((E2 - (x1*x3+x2*x4))*S - 2*E3 - (x1*x2^2+x2*x3^2+x3*x4^2+x4*x1^2)) +
	- (2*E4^3 - 2*R*E4^2 - b^2*E2*E4 + b*R*E4*S - 2*R^2*E4 + b^2*R*E2 - b*R^2*S + 2*R^3) # = 0
b^3*(x1*x3+x2*x4)*S + b^3*(x1*x2^2+x2*x3^2+x3*x4^2+x4*x1^2) +
	+ (2*E4^3 - 2*R*E4^2 - b^2*E2*E4 + b*R*E4*S - 2*R^2*E4 + 2*b^3*E3 - b^3*E2*S + b^2*R*E2 - b*R^2*S + 2*R^3) # = 0
b^3*(x1*x2^2+x2*x3^2+x3*x4^2+x4*x1^2) +
	+ (2*E4^3 + 2*b*E4^2*S - 2*R*E4^2 - b^2*E2*E4 + b*R*E4*S - 2*R^2*E4 + 2*b^3*E3 - b^3*E2*S +
		+ b^2*R*E2 + b^2*R*S^2 - 3*b*R^2*S + 2*R^3) # = 0
### Sum(x4^2*...) =>
(E3^2 - 2*E2*E4) + b*(x1*x2^2+x2*x3^2+x3*x4^2+x4*x1^2) - R*(S^2 - 2*E2) # = 0
### =>
b^2*E3^2 + 2*b^2*R*E2 - b^2*R*S^2 - 2*b^2*E2*E4 +
	- (2*E4^3 + 2*b*E4^2*S - 2*R*E4^2 - b^2*E2*E4 + b*R*E4*S - 2*R^2*E4 + 2*b^3*E3 - b^3*E2*S +
		+ b^2*R*E2 + b^2*R*S^2 - 3*b*R^2*S + 2*R^3) # = 0
### Eq 4:
2*E4^3 + 2*b*E4^2*S - 2*R*E4^2 + b^2*E2*E4 + b*R*E4*S - 2*R^2*E4 - b^2*E3^2 + 2*b^3*E3 - b^3*E2*S +
	- b^2*R*E2 + 2*b^2*R*S^2 - 3*b*R^2*S + 2*R^3 # = 0


### Eqs:
### Eq 1:
E4^4 - b^4*E4 + b^3*R*E3 - b^2*R^2*E2 + b*R^3*S - R^4 # = 0
### Eq 2:
4*b*R*E4^2*S - 8*R^2*E4^2 - 4*b^3*R*E3 - b^4*E2^2 + 4*b^2*R^2*E2 + 2*b^4*E3*S +
	- 4*b*R^3*S - b^5*S + b^2*R^2*S^2 - 4*b*R^3*S + 4*b^4*R + 8*R^4 # = 0
### Eq 3:
E3^2 + 3*b*E3 - 2*E2*E4 - 3*b*E2*S + 2*R*E2 + b*S^3 - R*S^2 # = 0
### Eq 4:
2*E4^3 + 2*b*E4^2*S - 2*R*E4^2 + b^2*E2*E4 + b*R*E4*S - 2*R^2*E4 - b^2*E3^2 + 2*b^3*E3 - b^3*E2*S +
	- b^2*R*E2 + 2*b^2*R*S^2 - 3*b*R^2*S + 2*R^3 # = 0

### TODO: solve;


### Test:
x1^2*x2^2 + b*x3 # - R
x2^2*x3^2 + b*x4 # - R
x3^2*x4^2 + b*x1 # - R
x4^2*x1^2 + b*x2 # - R


R = -2
b = 3
# Special Sub-Case type:
x1 = 1.3715291620 + 0.2900893695i;
x2 = 0.4884236632 + 1.7041597366i;
x3 = 1.3715291620 - 0.2900893695i;
x4 = 0.4884236632 - 1.7041597366i;
x = c(x1, x2, x3, x4);
E = debug.E(x)
S = E$S; E2 = E$E2; E3 = E$E3; E4 = E$E4;


### Classic Polynomial:
# x1^2*x2^2 = R - b*x3
(R - b*x3)*x3^2 + b*x1^2*x4 - x1^2*R
b*x3^3 - R*x3^2 - b*x1^2*x4 + x1^2*R
# b*x3^3*x4^2 + b^2*x1*x3 - b*R*x3
R*x3^2*x4^2 + b^2*x1*x3 - b*R*x3 + b*x1^2*x4^3 - R*x1^2*x4^2
# =>
(b^2*x1 - b*R)*x3 + b*x1^2*x4^3 - R*x1^2*x4^2 - b*R*x1 + R^2
# =>
(b*x1^2*x4^3 - R*x1^2*x4^2 - b*R*x1 + R^2)^2*x4^2 + (b*x1 - R)*(b^2*x1 - b*R)^2 # = 0
b^2*x1^4*x4^8 - 2*b*R*x1^4*x4^7 + R^2*x1^4*x4^6 + 2*b*R^2*x1^2*x4^5 - 2*b^2*R*x1^3*x4^5 - 2*R^3*x1^2*x4^4 +
	+ 2*b*R^2*x1^3*x4^4 + R^4*x4^2 + b^2*R^2*x1^2*x4^2 - 2*b*R^3*x1*x4^2 - b^2*R^3 + b^5*x1^3 +
	+ 3*b^3*R^2*x1 - 3*b^4*R*x1^2

### Eq 1 =>
(b^2*x1 - b*R)*x1^2*x2^2 - b*(b*x1^2*x4^3 - R*x1^2*x4^2 - b*R*x1 + R^2) - (b^2*x1 - b*R)*R # = 0
(b*x1^3 - R*x1^2)*x2^2 - (b*x1^2*x4^3 - R*x1^2*x4^2) # = 0
### Eq 4 =>
b^2*x2^2 - (x4^2*x1^2 - R)^2 # = 0
b^2*(b*x1^2*x4^3 - R*x1^2*x4^2) - (b*x1^3 - R*x1^2)*(x4^2*x1^2 - R)^2 # = 0
R*x1^6*x4^4 - b*x1^7*x4^4 + b^3*x1^2*x4^3 - b^2*R*x1^2*x4^2 - 2*R^2*x1^4*x4^2 + 2*b*R*x1^5*x4^2 +
	+ R^3*x1^2 - b*R^2*x1^3
# TODO

x3^2*(b*x3^3 - R*x3^2 + R*x1^2)^2 + b^3*x1^5 - b^2*R*x1^4 # = 0

pX.gen = function(x1, x4, b=0, R=0, coeff=1) data.frame(x1=x1, x4=x4, b=b, R=R, coeff=coeff);
p14 = data.frame(
	x1 = c(2, 2, 1, 0),
	x4 = c(3, 2, 0, 0),
	b  = c(1, 0, 1, 0),
	R  = c(0, 1, 1, 2),
	coeff = c(1,-1,-1,1)
)
p1a = data.frame(x1=c(1,0), b=c(1,0), R=c(0,1), coeff=c(1,-1));
p1b = pow.pm(p1a, 3); p1b$b = p1b$b + 2;
#
p14sq = pow.pm(p14, 2);
p14sq$x4 = p14sq$x4 + 2;
p14r = sum.pm(p14sq, p1b);
p14r = sort.pm(p14r, xn="x4");
p14r = p14r[,c("b", "R", "x1", "x4", "coeff")]
print.p(p14r, leading="x4")

#
p4a = data.frame(
	x1 = 2, x4 = c(3, 2), b = c(3, 2), R = c(0, 1),
	coeff = c(1,-1)
)
p4b = data.frame(x1=c(3, 2), b=c(1, 0), R=c(0,1), coeff=c(1,-1))
p4c = data.frame(x1=c(2, 0), x4=c(2, 0), R=c(0, 1), coeff=c(1,-1))
#
p4 = mult.pm(pow.pm(p4c, 2), p4b);
p4 = diff.pm(p4a, p4);
p4 = sort.pm(p4, xn="x4");
p4 = p4[,c("b", "R", "x1", "x4", "coeff")]
print.p(p4, leading="x4")

#
p4r = diff.pm(mult.pm(p14r, diff.pm(pX.gen(2,0,R=1), pX.gen(3,0,b=1,))),
	mult.pm(p4, pX.gen(x1=0, x4=4, b=2)))
# p4r$coeff = - p4r$coeff;
p4r$x1 = p4r$x1 - min(p4r$x1);
p4r = sort.pm(p4r, c(4,2), xn="x4");
p4r = p4r[,c("b", "R", "x1", "x4", "coeff")]
print.p(p4r, leading="x4")

2*b^2*R*x1^5*x4^7 - 2*b*R^2*x1^4*x4^7 - b^5*x4^7 - b*R^2*x1^5*x4^6 - 2*b^3*R*x1^3*x4^6 + R^3*x1^4*x4^6 +
	+ 2*b^2*R^2*x1^2*x4^6 + b^4*R*x4^6 + 2*b^3*R*x1^4*x4^5 - 4*b^2*R^2*x1^3*x4^5 + 2*b*R^3*x1^2*x4^5 +
	- 2*b^2*R^2*x1^4*x4^4 + 4*b*R^3*x1^3*x4^4 + b^3*R^2*x1*x4^4 - 2*R^4*x1^2*x4^4 - b^2*R^3*x4^4 +
	- b^3*R^2*x1^3*x4^2 + 3*b^2*R^3*x1^2*x4^2 - 3*b*R^4*x1*x4^2 + R^5*x4^2 - b^6*x1^4 + 4*b^5*R*x1^3 +
	- 6*b^4*R^2*x1^2 + 4*b^3*R^3*x1 - b^2*R^4


ct1 = data.frame(coeff=1);
p4_1  = replace.pm(replace.pm(p4, ct1, "R", 1), ct1, "b", 1);
p4r_1 = replace.pm(replace.pm(p4r, ct1, "R", 1), ct1, "b", 1);
solve.pm(p4_1, p4r_1, "x4", stop.at=1)


library(gmp)
p4_1  = p4; p4_1$coeff = as.bigz(p4_1$coeff);
p4r_1 = p4r; p4r_1$coeff = as.bigz(p4r_1$coeff);
l = solve.pm(p4_1, p4r_1, "x4", stop.at=NULL)

l = read.csv("P.S2l22.x1.csv", colClasses=c(rep("numeric", 4), "character"))
l$coeff = as.bigz(l$coeff)
head(l)

p4x = l[l$x4 == 1, -4]
head(p4x)

p40 = l[l$x4 == 0, -4]
p40$coeff = - p40$coeff
head(p40)

pr = replace.fr.bigpm(p4_1, p40, p4x, "x4", 1)

