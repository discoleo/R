########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: P[5]
### v.0.1b


### Experimental Approaches to P[5]
# - evaluating solutions for P[5];


###############
### History ###
###############

### draft v.0.1a:
# - moved to this file from file:
#   Polynomials.Derived.R;

####################

### Helper Functions

source("Polynomials.Helper.R")


#########################
#########################

### Debug:
R = 1
x = roots(c(1,0,0,0,-1,-R));
a1 = Re(x[2]); b1 = Im(x[2]);
a2 = Re(x[4]); b2 = Im(x[4]);

### Real root
a3 = -2*(a1+a2);

### Split Complex
a = a1; b = b1; bb = b*b; ad = 2*a;
a^5 - 10*a^3*b^2 + 5*a*b^4 - a - R # = 0
5*a^4 - 10*a^2*b^2 + b^4 - 1 # = 0

# =>
24*a^5 - 40*a^3*bb - 4*a + R # = 0
5*a^4 - 10*a^2*bb + bb^2 - 1 # = 0
# =>
1024*a^10 + 192*a^6 + 352*R*a^5 - 16*a^2 + 8*R*a - R^2 # = 0
# =>
ad^10 + 3*ad^6 + 11*R*ad^5 - 4*ad^2 + 4*R*ad - R^2 # = 0

# 2 of the roots are known (based on the original P[5]):
# c(2*a1, 2*a2, ...);
xx = roots(c(1, 0,0,0, 3, 11*R, 0,0,-4,4*R,-R^2))
bp = (3*xx^5 - 8*xx + 4*R) / (20*xx^3); # bp = b^2;
tmp = poly.calc0(c(xx/2 + sqrt(bp)*1i, xx/2 - sqrt(bp)*1i), digits=4)
# as.pm("(x^5 - x - R[])^4")
# TODO: explore ways to use this property;

# Test:
print.pm(div.pm(tmp, as.pm("x^5 - x - R[]"), by="x"))
-(x+1)^3 + 3*x^5*(x+1)^2 - 3*x^10*(x+1) + x^15
(x^5 - x - 1)^3

# Relations to original P[5]:
id = c(2,1,2,5,4,2,4,4,1,2)
5*xx^3*(2*x[id] - xx)^2 + 3*xx^5 - 8*xx + 4*R # = 0
2*xx^5 - 5*xx^4*x[id] + 5*xx^3*x[id]^2 - 2*xx + R # = 0


# Derivation:
p1 = as.pm("24*a^5 - 40*a^3*bb - 4*a + R")
p2 = as.pm("5*a^4 - 10*a^2*bb + bb^2 - 1")
pR = solve.pm(p1, p2, by="bb")


################

### E2
(a1^2 + b1^2) + (a2^2 + b2^2) - a3^2 + 4*a1*a2 # = 0
3*(a1^2 + a2^2) + 4*a1*a2 - b1^2 - b2^2 # = 0

### E3
2*a1*(a2^2+b2^2) + 2*a2*(a1^2+b1^2) - 2*(a1+a2)*(a1^2+b1^2+a2^2+b2^2 + 4*a1*a2) # = 0
a1*b2^2 + a2*b1^2 - (a1+a2)*(a1^2+a2^2 + b1^2+b2^2) - 3*a1^2*a2 - 3*a1*a2^2 # = 0
a1^3 + a2^3 + a1*b1^2 + a2*b2^2 + 4*a1^2*a2 + 4*a1*a2^2 # = 0

### E4
4*a2*(a1+a2)*(a1^2+b1^2) + 4*a1*(a1+a2)*(a2^2+b2^2) - (a1^2+b1^2)*(a2^2+b2^2) - 1 # = 0
2*R*a2*(a1^2+b1^2) + 2*R*a1*(a2^2+b2^2) + (a1^2+b1^2)^2*(a2^2+b2^2)^2 + (a1^2+b1^2)*(a2^2+b2^2) # = 0

### E5
2*(a1+a2)*(a1^2+b1^2)*(a2^2+b2^2) + R # = 0

### TODO


#########################
#########################

### Experimental

### x^5 - x = R

### roots: r1, Conj(r1) =>
(a1+b1*1i)^5 - (a1+b1*1i) - R # = 0
(a1-b1*1i)^5 - (a1-b1*1i) - R # = 0

### Diff =>
5*a1^4 - 10*a1^2*b1^2 + b1^4 - 1 # = 0
5*a2^4 - 10*a2^2*b2^2 + b2^4 - 1 # = 0

### Sum =>
a1^5 - 10*a1^3*b1^2 + 5*a1*b1^4 - a1 - R # = 0
a2^5 - 10*a2^3*b2^2 + 5*a2*b2^4 - a2 - R # = 0

### Diff Eq(r1) - Eq(r2):
(a1+b1*1i)^4 + (a2+b2*1i)^4 + (a1+b1*1i)*(a2+b2*1i)*(a1^2+a2^2-b1^2-b2^2 + 2*(a1*b1+a2*b2)*1i) +
	+ (a1+b1*1i)^2*(a2+b2*1i)^2 - 1 # = 0
a1^4 + a2^4 + b1^4 + b2^4 - 6*(a1^2*b1^2 + a2^2*b2^2) + 4*(a1^3*b1 + a2^3*b2 - a1*b1^3 - a2*b2^3)*1i +
	+ (a1*a2 - b1*b2 + (a1*b2+a2*b1)*1i)*(a1^2+a2^2-b1^2-b2^2 + 2*(a1*b1+a2*b2)*1i) +
	+ (a1*a2 - b1*b2 + (a1*b2+a2*b1)*1i)^2 - 1 # = 0
# Re =>
a1^4 + a2^4 + b1^4 + b2^4 - 6*(a1^2*b1^2 + a2^2*b2^2) +
	+ (a1*a2 - b1*b2)*(a1^2+a2^2 - (b1^2+b2^2)) - 2*(a1*b2 + a2*b1)*(a1*b1 + a2*b2) +
	+ (a1*a2 - b1*b2)^2 - (a1*b2 + a2*b1)^2 - 1 # = 0
# TODO + Im()

# [alternative] =>
a1^5 + 5*a1^4*b1*1i - 10*a1^3*b1^2 - 10*a1^2*b1^3*1i + 5*a1*b1^4 + b1^5*1i - (a1+b1*1i) +
	- a2^5 - 5*a2^4*b2*1i + 10*a2^3*b2^2 + 10*a2^2*b2^3*1i - 5*a2*b2^4 - b2^5*1i + (a2+b2*1i) # = 0
a1^5 - a2^5 - 10*a1^3*b1^2 + 10*a2^3*b2^2 + 5*a1*b1^4 - 5*a2*b2^4 - (a1-a2) # = 0


### Debug:
R = 1
x = roots(c(1,0,0,0,-1,-R));
a1 = Re(x[2]); b1 = Im(x[2]);
a2 = Re(x[4]); b2 = Im(x[4]);


####################
####################

### MPFR

### Example: Polynomial Roots with mpfr

# library(Rmpfr) # loaded automatically
source("Polynomials.Helper.mpfr.R")
source("Polynomials.Helper.BigNumbers.R")

### Example 1:
p = as.pm("x^5 - x - 1")
# Complex root:
pz = replace.pm(p, as.pm("a + b*i"), "x")
pz = replace.pm(pz, -1, "i", pow = 2)
p1 = drop.pm(pz[pz$i == 0, ])
p2 = pz[pz$i != 0, ]; p2$i = 0;
p2 = drop.pm(p2)

r = roots.pm(p)

p1a = dp.pm(p1, by="a")
p1b = dp.pm(p1, by="b")
p2a = dp.pm(p2, by="a")
p2b = dp.pm(p2, by="b")
p1a$coeff = mpfr(p1a$coeff, 240);
p1b$coeff = mpfr(p1b$coeff, 240);
p2a$coeff = mpfr(p2a$coeff, 240);
p2b$coeff = mpfr(p2b$coeff, 240);
#
p1$coeff = mpfr(p1$coeff, 240);
p2$coeff = mpfr(p2$coeff, 240);

a = mpfr(Re(r[1]), 240);
b = mpfr(Im(r[1]), 240);
#
for(step_i in seq(3)) {
v1a = eval.pm(p1a, list(a=a, b=b))
v1b = eval.pm(p1b, list(a=a, b=b))
v2a = eval.pm(p2a, list(a=a, b=b))
v2b = eval.pm(p2b, list(a=a, b=b))
v1  = eval.pm(p1, list(a=a, b=b))
v2  = eval.pm(p2, list(a=a, b=b))
#
div = - (v1a*v2b - v2a*v1b);
dx  = (v2b*v1 - v1b*v2) / div;
dy  = (v2a*v1 - v1a*v2) / - div;
# Test:
print(v1a*dx + v1b*dy)
print(v2a*dx + v2b*dy)
# Update:
a = a + dx; b = b + dy;
}


### Other:
z1 = as.pm("a - 76488*sc"); z1$coeff = mpfr(z1$coeff, 240);
z2 = as.pm("b - 35247*sc"); z2$coeff = mpfr(z2$coeff, 240);
sc = mpfr("1E-5", 240);
z1 = replace.pm(z1, sc, "sc");
z2 = replace.pm(z2, sc, "sc");
p1$coeff = mpfr(p1$coeff, 240);
p2$coeff = mpfr(p2$coeff, 240);
#
p1 = replace.pm(p1, z1, "a", tol=0);
p1 = replace.pm(p1, z2, "b", tol=0);
p2 = replace.pm(p2, z1, "a");
p2 = replace.pm(p2, z2, "b");

# TODO: solve system & iterate for higher precision;

