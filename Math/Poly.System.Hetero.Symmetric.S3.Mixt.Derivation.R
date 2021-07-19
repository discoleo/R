########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Heterogeneous Symmetric S3
### Mixed Type: Derivation
###
### draft v.0.3a


#####################

###############
### History ###
###############


### draft v.0.3a:
# - solved: Mixed Ht S3P31 + Symmetric P7;
#   x^7 + y^7 + z^7 = R2;
### draft v.0.2b:
# - solved: Mixed Ht S3P21 + Symmetric P3;
#   x^3 + y^3 + z^3 = R2;
### draft v.0.2a:
# - moved Derivation to this new file;
# - from file:
#   Poly.System.Hetero.Symmetric.S3.Mixt.R;
### [old]
### draft v.0.2k


######################
######################

##########################
### Mixed Systems      ###
### Type 1: 1 Rotation ###
##########################

###############
### Order 2 ###
###############

### n = 2
x*y^2 + y*z^2 + z*x^2 - R1 # = 0
x*y + y*z + z*x - R2 # = 0
x*y*z - R3 # = 0

### Solution:

### Eq 1: * x
x^2*y^2 + x*y*z^2 + z*x^3 - R1*x # = 0
R3^2/z^2 + R3*z + z*x^3 - R1*x # = 0 # *z^2
R3^2 + R3*z^3 + x^3*z^3 - R1*x*z^2 # = 0
# similar:
R3^2 + R3*x^3 + x^3*y^3 - R1*y*x^2 # = 0
R3^2 + R3*y^3 + y^3*z^3 - R1*z*y^2 # = 0

### Sum =>
R3*(x^3+y^3+z^3) + (x^3*y^3+x^3*z^3+y^3*z^3) - R1*(x^2*y + y^2*z + x*z^2) + 3*R3^2 # = 0
R3*(x^3+y^3+z^3) + (x^3*y^3+x^3*z^3+y^3*z^3) +
	- R1*(x^2*y + x^2*z + x*y^2 + y^2*z + x*z^2 + y*z^2) + R1^2 + 3*R3^2 # = 0
R3*(x^3+y^3+z^3) + (x^3*y^3+x^3*z^3+y^3*z^3) +
	- R1*((x^2+y^2+z^2)*(x+y+z) - (x^3+y^3+z^3)) + R1^2 + 3*R3^2 # = 0
(R1+R3)*(x^3+y^3+z^3) + (x^3*y^3+x^3*z^3+y^3*z^3) +
	- R1*((S^2 - 2*R2)*S) + R1^2 + 3*R3^2 # = 0
(R1+R3)*(S^3 - 3*R2*S + 3*R3) + (R2^3 - 3*R3*(R2*S - R3)) +
	- R1*((S^2 - 2*R2)*S) + R1^2 + 3*R3^2 # = 0
R3*S^3 - (R1+6*R3)*R2*S + R1^2 + R2^3 + 9*R3^2 + 3*R1*R3 # = 0

### Eq:
E3*S^3 - (R1+6*E3)*E2*S + R1^2 + E2^3 + 9*E3^2 + 3*R1*E3 # = 0

###############

###############
### Extensions:

### Extension A1:
### x*y^2 + y*z^2 + z*x^2 + b1*(x+y+z) = R1
### Extension A2:
### x*y + x*z + y*z + b2*(x+y+z) = R2;
### Extension A3:
### x*y*z + b3*(x+y+z) = R3;

### * x*z^2 =>
x^2*y^2*z^2 + x*y*z^4 + z^3*x^3 + b1*(x+y+z)*x*z^2 - R1*x*z^2 # = 0
R3^2 + R3*z^3 + z^3*x^3 + b1*S*x*z^2 - R1*x*z^2 # = 0
### Sum() =>
3*R3^2 + R3*(x^3+y^3+z^3) + (x^3*z^3+x^3*y^3+y^3*z^3) +
	+ (b1*S - R1)*(x*z^2 + x^2*y + y^2*z) # = 0
### Sum - R1*initial Eq
R3*(x^3+y^3+z^3) + (x^3*z^3+x^3*y^3+y^3*z^3) +
	(b1*S - R1)*(x*z^2 + x^2*y + y^2*z + x*y^2 + y*z^2 + z*x^2) + # = (E2*S - 3*E3)
	+ b1*S*(b1*S - R1) - R1*(b1*S - R1) + 3*R3^2 # = 0
R3*(S^3 - 3*R2*S + 3*R3) + (R2^3 - 3*R3*(R2*S - R3)) +
	(b1*S - R1)*((x^2+y^2+z^2)*(x+y+z) - (x^3+y^3+z^3)) +
	+ b1^2*S^2 - 2*R1*b1*S + R1^2 + 3*R3^2 # = 0
R3*S^3 + (b1*S - R1)*((x^2+y^2+z^2)*(x+y+z) - (x^3+y^3+z^3)) +
	+ b1^2*S^2 - 6*R2*R3*S - 2*R1*b1*S + R1^2 + R2^3 + 9*R3^2 # = 0
R3*S^3 + (b1*S - R1)*(R2*S - 3*R3) +
	+ b1^2*S^2 - 6*R2*R3*S - 2*R1*b1*S + R1^2 + R2^3 + 9*R3^2 # = 0
R3*S^3 + (b1^2 + b1*R2)*S^2 - (R1*R2 + 3*b1*R3 + 6*R2*R3 + 2*b1*R1)*S +
	+ R1^2 + R2^3 + 9*R3^2 + 3*R1*R3 # = 0

### Extension A3:
# - includes A1, but A2 was missed;
#   [see function solve.Ht3() for complete variant]
# E3 = R3 - b3*S
E3*S^3 + (b1^2 + b1*R2)*S^2 - (R1*R2 + 3*b1*E3 + 6*R2*E3 + 2*b1*R1)*S +
	+ R1^2 + R2^3 + 9*E3^2 + 3*R1*E3 # = 0
(R3 - b3*S)*S^3 + (b1^2 + b1*R2)*S^2 - (R1*R2 + 3*b1*(R3 - b3*S) + 6*R2*(R3 - b3*S) + 2*b1*R1)*S +
	+ R1^2 + R2^3 + 9*(R3 - b3*S)^2 + 3*R1*(R3 - b3*S) # = 0
-b3*S^4 + R3*S^3 + (b1^2 + 3*b1*b3 + 9*b3^2 + b1*R2 + 6*b3*R2)*S^2 +
	- (R1*R2 + 3*b1*R3 + 6*R2*R3 + 2*b1*R1 + 3*b3*R1 + 18*b3*R3)*S +
	+ R1^2 + R2^3 + 9*R3^2 + 3*R1*R3 # = 0


######################
######################

###############
### Order 3 ###
###############

x*y^3 + y*z^3 + z*x^3 - R1 # = 0
x*y + y*z + z*x - R2 # = 0
x*y*z - R3 # = 0

### Solution:

### Eq 1: * x*z^3
x^2*y^3*z^3 + x*y*z^6 + z^4*x^4 - R1*x*z^3 # = 0
R3^2*y*z + R3*z^5 + x^4*z^4 - R1*x*z^3 # = 0
# similar:
R3^2*x*z + R3*x^5 + x^4*y^4 - R1*y*x^3 # = 0
R3^2*x*y + R3*y^5 + y^4*z^4 - R1*z*y^3 # = 0

### Sum =>
R3^2*(x*y+x*z+y*z) + R3*(x^5+y^5+z^5) +
	+ (x^4*y^4+x^4*z^4+y^4*z^4) - R1*(y*x^3+z*y^3+x*z^3) # = 0
R3^2*R2 + R3*(S^5 - 5*R2*S^3 + 5*R3*S^2 + 5*R2^2*S  - 5*R2*R3) +
	+ (x^4*y^4+x^4*z^4+y^4*z^4) - R1*(y*x^3+z*y^3+x*z^3) # = 0
R3*S^5 - 5*R2*R3*S^3 + 5*R3^2*S^2 + 5*R2^2*R3*S +
	+ (4*R2*R3^2 + R2^4 - 4*R2^2*R3*S + 2*R3^2*S^2) - R1*(y*x^3+z*y^3+x*z^3) - 4*R2*R3^2 # = 0
R3*S^5 - 5*R2*R3*S^3 + 7*R3^2*S^2 + R2^2*R3*S +
	- R1*(y*x^3+z*y^3+x*z^3) + R2^4 # = 0
### Sum - R1*Initial_Eq =>
R3*S^5 - 5*R2*R3*S^3 + 7*R3^2*S^2 + R2^2*R3*S +
	- R1*(y*x^3+z*x^3+x*y^3+z*y^3+x*z^3+y*z^3) + R1^2 + R2^4 # = 0
R3*S^5 - 5*R2*R3*S^3 + 7*R3^2*S^2 + R2^2*R3*S +
	- R1*(R2*(S^2 - 2*R2) - R3*S) + R1^2 + R2^4 # = 0
R3*S^5 - 5*R2*R3*S^3 + (7*R3^2 - R1*R2)*S^2 + (R2^2*R3 + R1*R3)*S +
	+ R1^2 + 2*R1*R2^2 + R2^4 # = 0


#############################
#############################

### Generalization:
### x^p*y^n + y^p*z^n + z^p*x^n

#####################
### Higher Powers ###
### p > 1         ###
#####################

### p = 2
### x^2*y^n + y^2*z^n + z^2*x^n = R1

### Order 3: n = 3, p = 2
x^2*y^3 + y^2*z^3 + z^2*x^3 - R1 # = 0
x*y + y*z + z*x - R2 # = 0
x*y*z - R3 # = 0

### Solution:

### Eq 1: * x^3*y^2
x^5*y^5 + x^3*y^5*z^3 + x^6*y^2*z^2 - R1*x^3*y^2 # = 0
x^5*y^5 + R3^3*y + R3^2*x^4 - R1*x^3*y^2 # = 0
# similar:
x^5*z^5 + R3^3*x + R3^2*z^4 - R1*x^2*z^3 # = 0
y^5*z^5 + R3^3*z + R3^2*y^4 - R1*y^3*z^2 # = 0

### Sum =>
R3^3*(x+y+z) + R3^2*(x^4+y^4+z^4) +
	+ (x^5*y^5+x^5*z^5+y^5*z^5) - R1*(y^2*x^3+z^2*y^3+x^2*z^3) # = 0
R3^3*S + R3^2*(S^4 - 4*R2*S^2 + 4*R3*S + 2*R2^2) +
	+ (5*R2^2*R3^2 + R2^5 - 5*R2^3*R3*S - 5*R3^3*S + 5*R2*R3^2*S^2) +
	- R1*(y^2*x^3+z^2*y^3+x^2*z^3) # = 0
R3^2*S^4 + R2*R3^2*S^2 - 5*R2^3*R3*S +
	- R1*(y^2*x^3+z^2*y^3+x^2*z^3) + R2^5 + 7*R2^2*R3^2 # = 0

# Sum() - R1*(initial Eq) + R1^2 =>
R3^2*S^4 + R2*R3^2*S^2 - 5*R2^3*R3*S +
	- R1*(y^2*x^3+z^2*y^3+x^2*z^3 + x^2*y^3 + y^2*z^3 + z^2*x^3) +
	+ R1^2 + R2^5 + 7*R2^2*R3^2 # = 0
R3^2*S^4 + R2*R3^2*S^2 - 5*R2^3*R3*S +
	- R1*((R2^2 - 2*R3*S)*S - R3*R2) +
	+ R1^2 + R2^5 + 7*R2^2*R3^2 # = 0
R3^2*S^4 + (2*R1*R3 + R2*R3^2)*S^2 - (R1*R2^2 + 5*R2^3*R3)*S +
	+ R1^2 + R2^5 + 7*R2^2*R3^2 + R1*R2*R3 # = 0


#############################
#############################

################
### Variants ###
################

### Mixed with Symmetric
### P2+1 & Symmetric Order 3
### x^3 + y^3 + z^3 = R2

### Ht Order 2+1:
x*y^2 + y*z^2 + z*x^2 - R1 # = 0
x^3 + y^3 + z^3 - R2 # = 0
x*y*z - R3 # = 0

### Solution:

### Eq 2:
S^3 - 3*E2*S + 3*E3 - R2 # = 0
# =>
# 3*E2*S = S^3 + 3*E3 - R2

### Eq 1:
E3*S^3 - (R1+6*E3)*E2*S + R1^2 + E2^3 + 9*E3^2 + 3*R1*E3 # = 0
# =>
27*E3*S^6 - 27*(R1+6*E3)*E2*S^4 + 27*R1^2*S^3 + 27*E2^3*S^3 + 27*9*E3^2*S^3 + 81*R1*E3*S^3 # = 0
27*E3*S^6 - 9*(R1+6*E3)*(S^6 + 3*E3*S^3 - R2*S^3) + 27*R1^2*S^3 + (S^3 + 3*E3 - R2)^3 +
	+ 27*9*E3^2*S^3 + 81*R1*E3*S^3 # = 0
-27*E3*S^6 - 9*R1*S^6 +
	+ (S^3 + 3*E3 - R2)^3 +
	+ 81*E3^2*S^3 + 54*R1*E3*S^3 + 54*R2*E3*S^3 + 27*R1^2*S^3 + 9*R1*R2*S^3 # = 0
S^9 - (9*R1 + 3*R2 + 18*E3)*S^6 + (3*R2^2 + 36*R2*E3 + 108*E3^2 + 9*R1*R2 + 54*R1*E3 + 27*R1^2)*S^3 +
	- R2^3 + 27*E3^3 + 9*R2^2*E3 - 27*R2*E3^2


#############################
#############################

################
### Variants ###
################

### Mixed with Symmetric
### P3+1 & Symmetric Order 7
### x^7 + y^7 + z^7 = R2

### Ht Order 3+1:
x*y^3 + y*z^3 + z*x^3 - R1 # = 0
x^7 + y^7 + z^7 - R2 # = 0
x*y*z - R3 # = 0


### Solution:

### Eq 1:
E3*S^5 - 5*E2*E3*S^3 + (7*E3^2 - R1*E2)*S^2 + (E2^2*E3 + R1*E3)*S +
	+ R1^2 + 2*R1*E2^2 + E2^4 # = 0

R3*S^5 - 5*E2*R3*S^3 + (7*R3^2 - R1*E2)*S^2 + (E2^2*R3 + R1*R3)*S +
	+ R1^2 + 2*R1*E2^2 + E2^4 # = 0

### Eq 2:
S^7 - 7*E2*S^5 + 7*E3*S^4 + 14*E2^2*S^3 - 21*E3*E2*S^2 - 7*E2^3*S + 7*E3^2*S + 7*E3*E2^2 - R2 # = 0

### Eq S:
S^28 - 56*R3*S^25 - 28*R1*S^24 + 1190*R3^2*S^22 - 4*R2*S^21 + 1148*R3*R1*S^21 + 294*R1^2*S^20 +
	- 11564*R3^3*S^19 + 560*R2*R3*S^18 - 15386*R3^2*R1*S^18 - 210*R2*R1*S^17 - 9114*R3*R1^2*S^17 +
	+ 47383*R3^4*S^16 - 1274*R1^3*S^16 - 13370*R2*R3^2*S^15 + 62818*R3^3*R1*S^15 + 6*R2^2*S^14 +
	+ 476*R2*R3*R1*S^14 + 81928*R3^2*R1^2*S^14 - 58996*R3^5*S^13 + 784*R2*R1^2*S^13 + 19208*R3*R1^3*S^13 +
	+ 91728*R2*R3^3*S^12 + 93982*R3^4*R1*S^12 + 6517*R1^4*S^12 +
	+ 1792*R2^2*R3*S^11 + 41062*R2*R3^2*R1*S^11 - 173558*R3^3*R1^2*S^11 +
	+ 148862*R3^6*S^10 + 161*R2^2*R1*S^10 + 10682*R2*R3*R1^2*S^10 - 5831*R3^2*R1^3*S^10 +
	- 48118*R2*R3^4*S^9 - 461678*R3^5*R1*S^9 + 5978*R2*R1^3*S^9 - 25382*R3*R1^4*S^9 +
	+ 12194*R2^2*R3^2*S^8 - 72128*R2*R3^3*R1*S^8 + 124509*R3^4*R1^2*S^8 + 4802*R1^5*S^8 +
	- 4*R2^3*S^7 - 67228*R3^7*S^7 + 3150*R2^2*R3*R1*S^7 + 9016*R2*R3^2*R1^2*S^7 - 92610*R3^3*R1^3*S^7 +
	- 104272*R2*R3^5*S^6 + 177674*R3^6*R1*S^6 + 1323*R2^2*R1^2*S^6 - 4802*R2*R3*R1^3*S^6 + 50421*R3^2*R1^4*S^6 +
	+ 6272*R2^2*R3^3*S^5 + 26068*R2*R3^4*R1*S^5 - 57624*R3^5*R1^2*S^5 + 686*R2*R1^4*S^5 - 4802*R3*R1^5*S^5 +
	+ 105*R2^3*R3*S^4 + 117649*R3^8*S^4 - 1666*R2^2*R3^2*R1*S^4 - 1715*R2*R3^3*R1^2*S^4 +
		+ 72030*R3^4*R1^3*S^4 + 2401*R1^6*S^4 +
	- 4802*R2*R3^6*S^3 + 77*R2^3*R1*S^3 - 33614*R3^7*R1*S^3 + 833*R2^2*R3*R1^2*S^3 +
		+ 1029*R2*R3^2*R1^3*S^3 - 7203*R3^3*R1^4*S^3 +
	+ 735*R2^2*R3^4*S^2 + 10290*R2*R3^5*R1*S^2 + 36015*R3^6*R1^2*S^2 + 98*R2^2*R1^3*S^2 +
		+ 1372*R2*R3*R1^4*S^2 + 4802*R3^2*R1^5*S^2 +
	- 14*R2^3*R3^2*S - 294*R2^2*R3^3*R1*S - 2058*R2*R3^4*R1^2*S - 4802*R3^5*R1^3*S +
	+ R2^4 + 28*R2^3*R3*R1 + 294*R2^2*R3^2*R1^2 + 1372*R2*R3^3*R1^3 + 2401*R3^4*R1^4

### Derivation:
library(gmp)
#
p1 = toPoly.pm("E3*S^5 - 5*E2*E3*S^3 + 7*E3^2*S^2 - R1*E2*S^2 + E2^2*E3*S + R1*E3*S + R1^2 + 2*R1*E2^2 + E2^4")
p2 = toPoly.pm("S^7 - 7*E2*S^5 + 7*E3*S^4 + 14*E2^2*S^3 - 21*E3*E2*S^2 - 7*E2^3*S + 7*E3^2*S + 7*E3*E2^2 - R2")
p1$coeff = as.bigz(p1$coeff)
p2$coeff = as.bigz(p2$coeff)
#
pR = solve.pm(p1, p2, "E2")
div = gcd.vpm(pR$Rez, as.bigz(0));
print(div)
pR$Rez$coeff = pR$Rez$coeff / div;
denominator(pR$Rez$coeff);
pR$Rez$coeff = as.bigz(pR$Rez$coeff);
id = order( - pR$Rez$S); pR$Rez = pR$Rez[id,];
print.coeff(pR$Rez, "S")
#
xgcd = gcd.vpm(pR$x0, as.bigz(0));
xgcd = gcd.vpm(pR$div, xgcd)
if(xgcd != 1) {
	pR$x0$coeff  = as.bigz(pR$x0$coeff / xgcd);
	pR$div$coeff = as.bigz(pR$div$coeff / xgcd);
}
id = order( - pR$x0$S); pR$x0 = pR$x0[id,];
id = order( - pR$div$S); pR$div = pR$div[id,];
print.p(pR$x0, "S")
print.p(pR$div, "S")
# TODO: implement multi-variable
pR2 = factorize.p(pR$Rez, xn="S")
#
9*S^12 + 12*R3*S^9 + 12*R1*S^8 + 10*R3^2*S^6 + 8*R1*R3*S^5 + 4*R1^2*S^4 + 4*R3^3*S^3 + 4*R1*R3^2*S^2 + R3^4
#
pDiv = toPoly.pm("9*S^12 + 12*R3*S^9 + 12*R1*S^8 + 10*R3^2*S^6 + 8*R1*R3*S^5 + 4*R1^2*S^4 + 4*R3^3*S^3 + 4*R1*R3^2*S^2 + R3^4")
pDiv$coeff = as.bigz(pDiv$coeff);
names(pDiv)[match("R3", names(pDiv))] = "E3";
pR2 = div.pm(pR$Rez, pDiv, "S")
id = order( - pR2$Rez$S); pR2$Rez = pR2$Rez[id,];
print.coeff(pR2$Rez, "S")

pS = pR2$Rez;
names(pS)[match("E3", names(pS))] = "R3";
print.p(pS, "S")


solve.S3P21SymmP7 = function(R, debug=TRUE) {
	coeff = coeff.S3P21SymmP7(R);
	S = roots(coeff);
	if(debug) print(S);
	len = length(S);
	E2 = E2.S3P21SymmP7(S, R);
	E3 = R[3]; E3 = rep(E3, len);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	S = rep(S, each=3); E2 = rep(E2, each=3); E3 = rep(E3, each=3);
	yz.s = S - x; yz = E3 / x;
	# robust
	R1 = R[1]; R2 = E2; R3 = E3;
	x3 = if(R1 == 0) {
		# with chain rule!
		# x3 = (5*R[3]*S^4 - 15*R[2]*R[3]*S^2 + 14*R[3]^2*S + R[2]^2*R[3]) # * dS/dR
		dS = - R2*S^2 + R3*S + 2*R2^2;
		x3 = - dS;
		x3
	} else {
		x3 = (R3*S^5 - 5*R2*R3*S^3 + 7*R3^2*S^2 + R2^2*R3*S + R2^4) / R1;
		x3
	};
	yz.d = (x3 - R1) / (x^3 + yz*yz.s - x*(yz.s^2 - yz));
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
}
coeff.S3P21SymmP7 = function(R) {
	R1 = R[1]; R2 = R[2]; E3 = R[3];
	coeff = c(1, 0, 0, - 56*E3, - 28*R1, 0, 1190*E3^2, - 4*R2 + 1148*E3*R1, 294*R1^2,
		- 11564*E3^3, 560*R2*E3 - 15386*E3^2*R1, - 210*R2*R1 - 9114*E3*R1^2,
		47383*E3^4 - 1274*R1^3, - 13370*R2*E3^2 + 62818*E3^3*R1,
		6*R2^2 + 476*R2*E3*R1 + 81928*E3^2*R1^2, - 58996*E3^5 + 784*R2*R1^2 + 19208*E3*R1^3,
		91728*R2*E3^3 + 93982*E3^4*R1 + 6517*R1^4,
		1792*R2^2*E3 + 41062*R2*E3^2*R1 - 173558*E3^3*R1^2,
		148862*E3^6 + 161*R2^2*R1 + 10682*R2*E3*R1^2 - 5831*E3^2*R1^3,
		- 48118*R2*E3^4 - 461678*E3^5*R1 + 5978*R2*R1^3 - 25382*E3*R1^4,
		12194*R2^2*E3^2 - 72128*R2*E3^3*R1 + 124509*E3^4*R1^2 + 4802*R1^5,
		- 4*R2^3 - 67228*E3^7 + 3150*R2^2*E3*R1 + 9016*R2*E3^2*R1^2 - 92610*E3^3*R1^3,
		- 104272*R2*E3^5 + 177674*E3^6*R1 + 1323*R2^2*R1^2 - 4802*R2*E3*R1^3 + 50421*E3^2*R1^4,
		6272*R2^2*E3^3 + 26068*R2*E3^4*R1 - 57624*E3^5*R1^2 + 686*R2*R1^4 - 4802*E3*R1^5,
		105*R2^3*E3 + 117649*E3^8 - 1666*R2^2*E3^2*R1 - 1715*R2*E3^3*R1^2 + 72030*E3^4*R1^3 + 2401*R1^6,
		- 4802*R2*E3^6 + 77*R2^3*R1 - 33614*E3^7*R1 + 833*R2^2*E3*R1^2 + 1029*R2*E3^2*R1^3 - 7203*E3^3*R1^4,
		735*R2^2*E3^4 + 10290*R2*E3^5*R1 + 36015*E3^6*R1^2 + 98*R2^2*R1^3 + 1372*R2*E3*R1^4 + 4802*E3^2*R1^5,
		- 14*R2^3*E3^2 - 294*R2^2*E3^3*R1 - 2058*R2*E3^4*R1^2 - 4802*E3^5*R1^3,
		R2^4 + 28*R2^3*E3*R1 + 294*R2^2*E3^2*R1^2 + 1372*R2*E3^3*R1^3 + 2401*E3^4*R1^4);
	return(coeff)
}
E2.S3P21SymmP7 = function(S, R) {
	R1 = R[1]; R2 = R[2]; E3 = R[3];
	px0 = 5*S^17 - 57*E3*S^14 + 42*R1*S^13 - 343*E3^2*S^11 - 3*R2*S^10 - 49*E3*R1*S^10 - 175*R1^2*S^9 +
		+ 2541*E3^3*S^8 - 89*R2*E3*S^7 - 574*E3^2*R1*S^7 - 42*R2*R1*S^6 + 245*E3*R1^2*S^6 - 245*E3^4*S^5 +
		- 147*R1^3*S^5 - 637*E3^3*R1*S^4 - 2*R2^2*S^3 - 343*E3^5*S^2 - 21*R2*R1^2*S^2 - 98*E3*R1^3*S^2 +
		+ 7*R2*E3^3*S + 49*E3^4*R1*S - R2^2*E3 - 14*R2*E3^2*R1 - 49*E3^3*R1^2;
	pDiv = 22*S^15 - 448*E3*S^12 + 84*R1*S^11 + 2030*E3^2*S^9 + 26*R2*S^8 + 98*E3*R1*S^8 - 98*R1^2*S^7 +
		- 392*E3^3*S^6 + 154*R2*E3*S^5 - 1176*E3^2*R1*S^5 + 14*R2*R1*S^4 + 294*E3*R1^2*S^4 - 686*E3^4*S^3 +
		- 98*R1^3*S^3 + 28*R2*E3^2*S^2 + 196*E3^3*R1*S^2 + R2^2*S - 49*E3^2*R1^2*S;
	return(px0 / pDiv);
}

### Test
R = c(-1,2,-2)
sol = solve.S3P21SymmP7(R);
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
round0(x*y^3 + y*z^3 + z*x^3) # - R1
round0(x^7 + y^7 + z^7) # - R2
x*y*z # - R3 # = 0

round0.p(poly.calc(x), tol=0.5)


### Debug:
R = c(-1, 2, -2); R1 = R[1]; R2 = R[2]; R3 = R[3];
x = -1.2774896227 + 0.5473941212i;
y =  1.2113934262 + 0.6387638767i;
z =  1.0473828355 - 0.0844137651i;
S = (x+y+z); E2 = x*y+x*z+y*z; E3 = R[3];


#############################
#############################

###########################
### Mixed Systems       ###
### Type 2: 2 Rotations ###
###########################

### Dual Eq: p = 1
### x*y^n + y*z^n + z*x^n = R1
### x*z^n + y*x^n + z*y^n = R2

###############
### Order 2 ###
###############

x*y^2 + y*z^2 + z*x^2 - R1 # = 0
x*z^2 + y*x^2 + z*y^2 - R2 # = 0
x*y*z - R3 # = 0

### Solution:

### Eq 1 + Eq 2 =>
x*y^2 + y*z^2 + z*x^2 + x*z^2 + y*x^2 + z*y^2 - R1 - R2 # = 0
E2*S - 3*E3 - R1 - R2 # = 0
# E2 = (R1 + R2 + 3*E3) / S

### see Section "Simple System":
R3*S^3 - (R1+6*R3)*E2*S + R1^2 + E2^3 + 9*R3^2 + 3*R1*R3 # = 0
R3*S^3 - (R1+6*R3)*(R1 + R2 + 3*R3) + R1^2 + (3*R3 + R1 + R2)^3 / S^3 + 9*R3^2 + 3*R1*R3 # = 0
R3*S^6 - ((R1+6*R3)*(R1 + R2 + 3*R3) - R1^2 - 9*R3^2 - 3*R1*R3)*S^3 + (R1 + R2 + 3*R3)^3 # = 0
R3*S^6 - (R1*R2 + 6*R1*R3 + 6*R2*R3 + 9*R3^2)*S^3 + (R1 + R2 + 3*R3)^3 # = 0

### Extension M3:
R3*S^5 - (R1*R2 + 6*R2*R3/S + 9*R3^2/S^2 + 6*R1*R3/S)*S^3 + (R1 + R2 + 3*R3/S)^3 # = 0
R3*S^8 - (R1*R2 + 6*R2*R3/S + 9*R3^2/S^2 + 6*R1*R3/S)*S^6 + ((R1 + R2)*S + 3*R3)^3 # = 0
R3*S^8 - R1*R2*S^6 - 6*R3*(R1 + R2)*S^5 - 9*R3^2*S^4 + (R1 + R2)^3*S^3 +
	+ 9*R3*(R1 + R2)^2*S^2 + 27*R3^2*(R1 + R2)*S + 27*R3^3 # = 0


##################
##################
##################

##################
### Resonances ###
##################

######################
### Roots of Unity ###
######################

### Order: k + n
# - Base-Eq with terms: x^k*y^n, y^k*z^n, z^k*x^n
# - roots of unity: m^p = 1;
# - reusing: "x" = power of m in x;

k*x + n*y # = 0 (mod p);
k*y + n*z # = 0 (mod p);
k*z + n*x # = 0 (mod p);

### Sum =>
(k+n)*(x+y+z) # = 0 (mod p);

### Reduction =>
2*(k*y + n*z) + (k*z + n*x) # = 0 (mod p);
2*k*y + (2*n+k)*z) + n*x # = 0 (mod p);
### Diff((2*n+k)*(x+y+z) - ...) =>
(2*n-k)*y + (n+k)*x # = 0 (mod p);
### Reduction with Eq 1:
# n*Eq... - (2*n-k)*Eq_1  =>
(n*(n+k) - k*(2*n-k))*x # = 0 (mod p);
(n^2 + k^2 - n*k)*x # = 0 (mod p);

### Trivial Powers
# p = (k+n) # or one of its Divisors;

### Non-Trivial Powers
# p = (n^2 + k^2 - n*k) # or one of its Divisors;

### Combinations:
# p = Divisors of (k+n)*(n^2 + k^2 - n*k)

