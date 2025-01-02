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

