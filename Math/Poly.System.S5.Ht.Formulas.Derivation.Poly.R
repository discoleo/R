########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S5: Hetero-Symmetric
### Derivation / Intermediary Polynomials
###
### draft v.0.1b


### Derivation:
# - Robust formulas for x2 & x3;
# - Intermediary Polynomials;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


###############

### Derivation:

p1 = toPoly.pm("x2 + x3 + x4 + x5 - s");
p2 = toPoly.pm("x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1 - E11a");
p3 = toPoly.pm("x1*x3 + x2*x4 + x3*x5 + x4*x1 + x5*x2 - E11b");
p4 = toPoly.pm("x2*x3*x4 + x2*x3*x5 + x2*x4*x5 + x3*x4*x5 - e3")

pR = solve.lpm(p1,p2,p3,p4, xn=c("x5", "x4"))

px1 = pR[[2]]$Rez
px1 = sort.pm(px1, c("x2", "x3"))
print.pm(px1, lead=NA)

px2 = pR[[3]]
px2 = sort.pm(px2, c("x2", "x3"))
print.pm(px2, lead=NA)

###
x2^4 - s*x2^3 + e2*x2^2 - e3*x2 + e4 # = 0
x3^4 - s*x3^3 + e2*x3^2 - e3*x3 + e4 # = 0

### px1:
x2^4 + x3^4 + 3*x3*x2^3 + 2*x3^3*x2 + 4*x3^2*x2^2 - 2*s*x2^3 + x2^3*x1 - s*x3^3 - 2*x3^3*x1 +
	- (4*s + x1)*x3*x2^2 - (3*s + 2*x1)*x3^2*x2 +
	+ (s^2 - 2*s*x1 + x1^2 + 2*E11b)*x2^2 +
	+ (x1^2 + E11a + 2*E11b)*x3^2 + (s^2 + 3*E11b)*x2*x3 +
	+ (s^2*x1 - s*x1^2 - 2*s*E11b + x1*E11b)*x2 +
	+ (s^2*x1 + 2*s*x1^2 - 2*E11a*x1 - s*E11b - 3*x1*E11b)*x3 +
	- s*x1^3 + E11a*x1^2 + x1^2*E11b - s*x1*E11b + E11b^2 # = 0

b00 = - s*x1^3 + E11a*x1^2 + x1^2*E11b - s*x1*E11b + E11b^2;
b10 = s^2*x1 - s*x1^2 - 2*s*E11b + x1*E11b;
b01 = s^2*x1 + 2*s*x1^2 - 2*E11a*x1 - s*E11b - 3*x1*E11b;
b20 = s^2 - 2*s*x1 + x1^2 + 2*E11b;
b02 = x1^2 + E11a + 2*E11b;
b11 = s^2 + 3*E11b;
#
b20r = b20 - e2; b02r = b02 - e2;
b10r = b10 + e3; b01r = b01 + e3;
b00r = b00 - 2*e4; b11r = b11;


x2^4 + x3^4 + 3*x3*x2^3 + 2*x3^3*x2 + 4*x3^2*x2^2 - 2*s*x2^3 + x2^3*x1 - s*x3^3 - 2*x3^3*x1 +
	- (4*s + x1)*x3*x2^2 - (3*s + 2*x1)*x3^2*x2 +
	+ b20*x2^2 + b02*x3^2 + b11*x2*x3 +
	+ b10*x2 + b01*x3 + b00 # = 0
# Subst =>
3*x3*x2^3 + 2*x3^3*x2 + 4*x3^2*x2^2 - s*x2^3 + x2^3*x1 - 2*x3^3*x1 +
	- (4*s + x1)*x3*x2^2 - (3*s + 2*x1)*x3^2*x2 +
	+ (b22 - e2)*x2^2 + (b32 - e2)*x3^2 + bm2*x2*x3 +
	+ (b21 + e3)*x2 + (b31 + e3)*x3 + b00 - 2*e4 # = 0
# =>
3*x3*x2^3 + 2*x3^3*x2 + 4*x3^2*x2^2 - s*x2^3 + x2^3*x1 - 2*x3^3*x1 +
	- (4*s + x1)*x3*x2^2 - (3*s + 2*x1)*x3^2*x2 +
	+ b20r*x2^2 + b02r*x3^2 + b11r*x2*x3 +
	+ b10r*x2 + b01r*x3 + b00r # = 0

### px2:
x2^5 + 4*x3*x2^4 + 2*x3^4*x2 + 6*x3^2*x2^3 + 5*x3^3*x2^2 - 2*s*x2^4 + x2^4*x1 +
	- 6*s*x3*x2^3 + 2*x3*x2^3*x1 - 3*s*x3^3*x2 - 2*x3^3*x2*x1 - (6*s + x1)*x3^2*x2^2 +
	+ (s^2 - 2*s*x1 + 2*E11b)*x2^3 + (E11b - s*x1)*x3^3 +
	+ (2*s^2 - 4*s*x1 + 5*E11b)*x3*x2^2 +
	+ (s^2 - s*x1 + 4*E11b)*x3^2*x2 +
	+ (s^2*x1 - 2*s*E11b + E11b*x1)*x2^2 +
	+ (s^2*x1 + s*x1^2 - s*E11b - E11b*x1 + e3)*x3^2 +
	+ (2*s^2*x1 - 3*s*E11b)*x2*x3 +
	+ (E11b^2 - s*E11b*x1)*x2 +
	- (s*E11b*x1 - E11b^2 + 2*e3*x1)*x3 + e3*x1^2 # = 0

c00 = e3*x1^2;
c10 = E11b^2 - s*E11b*x1;
c01 = - (s*E11b*x1 - E11b^2 + 2*e3*x1);
c20 = s^2*x1 - 2*s*E11b + E11b*x1;
c02 = s^2*x1 + s*x1^2 - s*E11b - E11b*x1 + e3;
c11 = 2*s^2*x1 - 3*s*E11b;
c21 = 2*s^2 - 4*s*x1 + 5*E11b;
c12 = s^2 - s*x1 + 4*E11b;
c30 = s^2 - 2*s*x1 + 2*E11b;
c03 = E11b - s*x1;
#
c00r = c00 - e4*(x1 - s);
c10r = c10 - 3*e4 - e3*s + e3*x1;
c01r = c01 - 4*e4;
c11r = c11 + 6*e3; c02r = c02;
c20r = c20 + e3 + e2*s - e2*x1;
c21r = c21 - 4*e2; c12r = c12 - 2*e2;
c30r = c30 - s^2 - e2 + x1*s; c03r = c03;

x2^5 + 4*x3*x2^4 + 2*x3^4*x2 + 6*x3^2*x2^3 + 5*x3^3*x2^2 - 2*s*x2^4 + x2^4*x1 +
	- 6*s*x3*x2^3 + 2*x3*x2^3*x1 - 3*s*x3^3*x2 - 2*x3^3*x2*x1 - (6*s + x1)*x3^2*x2^2 +
	+ c30*x2^3 + c03*x3^3 + c21*x3*x2^2 + c12*x3^2*x2 +
	+ c20*x2^2 + c02*x3^2 + c11*x2*x3 +
	+ c10*x2 + c01*x3 + c00 # = 0
# Subst =>
6*x3^2*x2^3 + 5*x3^3*x2^2 +
	- 2*s*x3*x2^3 + 2*x3*x2^3*x1 - s*x3^3*x2 - 2*x3^3*x2*x1 - (6*s + x1)*x3^2*x2^2 +
	+ (c30 - s^2 - e2 + x1*s)*x2^3 + c03*x3^3 + (c21 - 4*e2)*x3*x2^2 + (c12 - 2*e2)*x3^2*x2 +
	+ (c20 + e3 + e2*s - e2*x1)*x2^2 + c02*x3^2 + (c11 + 6*e3)*x2*x3 +
	+ (c10 - 3*e4 - e3*s + e3*x1)*x2 + (c01 - 4*e4)*x3 + c00 - e4*(x1 - s) # = 0
# =>
6*x3^2*x2^3 + 5*x3^3*x2^2 +
	- 2*s*x3*x2^3 + 2*x3*x2^3*x1 - s*x3^3*x2 - 2*x3^3*x2*x1 - (6*s + x1)*x3^2*x2^2 +
	+ c30r*x2^3 + c03r*x3^3 + c21r*x3*x2^2 + c12r*x3^2*x2 +
	+ c20r*x2^2 + c02r*x3^2 + c11r*x2*x3 +
	+ c10r*x2 + c01r*x3 + c00r # = 0

###
p1 = toPoly.pm("3*x3*x2^3 + 2*x3^3*x2 + 4*x3^2*x2^2 - s*x2^3 + x2^3*x1 - 2*x3^3*x1 +
	- (4*s + x1)*x3*x2^2 - (3*s + 2*x1)*x3^2*x2 +
	+ b20r*x2^2 + b02r*x3^2 + b11r*x2*x3 +
	+ b10r*x2 + b01r*x3 + b00r");
p2 = toPoly.pm("6*x3^2*x2^3 + 5*x3^3*x2^2 +
	- 2*s*x3*x2^3 + 2*x3*x2^3*x1 - s*x3^3*x2 - 2*x3^3*x2*x1 - (6*s + x1)*x3^2*x2^2 +
	+ c30r*x2^3 + c03r*x3^3 + c21r*x3*x2^2 + c12r*x3^2*x2 +
	+ c20r*x2^2 + c02r*x3^2 + c11r*x2*x3 +
	+ c10r*x2 + c01r*x3 + c00r")

pR = solve.pm(p1, p2, "x3")

# x2^5
- 15*x3*x2^5 + (- 8*x3^2 + 15*x3*x1 + 19*x3*s)*x2^4 +
	+ (4*x1 + 7*s)*x2^3*x3^2 +
	+ (2*c21r - 3*c03r - 6*x1^2 - 5*x1*s - 4*s^2 - 5*b11r)*x2^3*x3 +
	# x2^2
	+ (2*c12r - 4*c03r - 2*x1^2 + 4*x1*s - 3*s^2 - 5*b02r)*x2^2*x3^2 +
	+ (2*c11r - 2*c21r*x1 + c03r*x1 + 4*c03r*s - 5*b01r + 2*x1*b11r + s*b11r)*x2^2*x3 +
	# x2^1
	+ (2*c02r - 2*c12r*x1 + 2*c03r*x1 + 3*c03r*s + 2*x1*b02r + s*b02r)*x3^2*x2 +
	+ (2*c01r - 2*c11r*x1 + 2*x1*b01r + s*b01r - c03r*b11r)*x2*x3 +
	- (2*c02r*x1 + c03r*b02r)*x3^2 - (2*c01r*x1 + c03r*b01r)*x3 +
	# x2-Only:
	5*(s - x1)*x2^5 + (2*c30r + 2*x1^2 - x1*s - s^2 - 5*b20r)*x2^4 +
	+ (2*c20r - c03r*x1 - 2*c30r*x1 + c03r*s - 5*b10r + 2*x1*b20r + s*b20r)*x2^3 +
	+ (2*c10r - 2*c20r*x1 - 5*b00r + 2*x1*b10r + s*b10r - c03r*b20r)*x2^2 +
	+ (2*c00r - 2*c10r*x1 + 2*x1*b00r + s*b00r - c03r*b10r)*x2 +
	# B0
	- 2*c00r*x1 - c03r*b00r # = 0

d00 = - 2*c00r*x1 - c03r*b00r;
d10 = 2*c00r - 2*c10r*x1 + 2*x1*b00r + s*b00r - c03r*b10r;
d20 = 2*c10r - 2*c20r*x1 - 5*b00r + 2*x1*b10r + s*b10r - c03r*b20r;
d30 = 2*c20r - c03r*x1 - 2*c30r*x1 + c03r*s - 5*b10r + 2*x1*b20r + s*b20r;
d40 = 2*c30r + 2*x1^2 - x1*s - s^2 - 5*b20r;
d50 = 5*(s - x1);
# x3:
d01 = - (2*c01r*x1 + c03r*b01r);
d02 = - (2*c02r*x1 + c03r*b02r);
d11 = 2*c01r - 2*c11r*x1 + 2*x1*b01r + s*b01r - c03r*b11r;
d12 = 2*c02r - 2*c12r*x1 + 2*c03r*x1 + 3*c03r*s + 2*x1*b02r + s*b02r;
d21 = 2*c11r - 2*c21r*x1 + c03r*x1 + 4*c03r*s - 5*b01r + 2*x1*b11r + s*b11r;
d22 = 2*c12r - 4*c03r - 2*x1^2 + 4*x1*s - 3*s^2 - 5*b02r;
d31 = 2*c21r - 3*c03r - 6*x1^2 - 5*x1*s - 4*s^2 - 5*b11r;
#
d01r = d01 - 15*x1*e4 - 4*s*e4;
d11r = d11 + 15*e4 + 15*x1*e3 + 4*s*e3;
d02r = d02 + 8*e4; d12r = d12 - 8*e3;
d21r = d21 - 15*e3 - 15*x1*e2 - 4*s*e2;
d31r = d31 + 15*e2 + 15*x1*s + 4*s^2;
d22r = d22 + 8*e2;

(4*x1 + 7*s)*x2^3*x3^2 - 8*x2^4*x3^2 +
	- 15*x3*x2^5 + (15*x1 + 19*s)*x2^4*x3 +
	+ d31*x2^3*x3 + d22*x2^2*x3^2 + d21*x2^2*x3 +
	+ d12*x3^2*x2 + d11*x2*x3 + d02*x3^2 + d01*x3 +
	# x2-Only:
	d50*x2^5 + d40*x2^4 + d30*x2^3 + d20*x2^2 + d10*x2 + d00 # = 0
# Subst =>
(4*x1 - s)*x2^3*x3^2 +
	+ (d31 + 15*e2 + 15*x1*s + 4*s^2)*x2^3*x3 + (d22 + 8*e2)*x2^2*x3^2 +
	+ (d21 - 15*e3 - 15*x1*e2 - 4*s*e2)*x2^2*x3 +
	+ (d12 - 8*e3)*x3^2*x2 + (d11 + 15*e4 + 15*x1*e3 + 4*s*e3)*x2*x3 +
	+ (d02 + 8*e4)*x3^2 + (d01 - 15*x1*e4 - 4*s*e4)*x3 +
	# x2-Only:
	d50*x2^5 + d40*x2^4 + d30*x2^3 + d20*x2^2 + d10*x2 + d00 # = 0
# =>
(4*x1 - s)*x2^3*x3^2 +
	+ d31r*x2^3*x3 + d22r*x2^2*x3^2 + d21r*x2^2*x3 +
	+ d12r*x3^2*x2 + d11r*x2*x3 + d02r*x3^2 + d01r*x3 +
	# x2-Only:
	d50*x2^5 + d40*x2^4 + d30*x2^3 + d20*x2^2 + d10*x2 + d00 # = 0

###
2*x2*x3^3 - 2*x1*x3^3 + (4*x2^2 - 2*x2*x1 - 3*x2*s + b02r)*x3^2 +
	+ (3*x2^3 - x2^2*x1 - 4*x2^2*s + b11r*x2 + b01r)*x3 +
	- (s - x1)*x2^3 + b20r*x2^2 + b10r*x2 + b00r # = 0


########
p1 = toPoly.pm("2*x2*x3^3 - 2*x1*x3^3 + (4*x2^2 - 2*x2*x1 - 3*x2*s + b02r)*x3^2 +
	+ (3*x2^3 - x2^2*x1 - 4*x2^2*s + b11r*x2 + b01r)*x3 +
	- (s - x1)*x2^3 + b20r*x2^2 + b10r*x2 + b00r");
p2 = toPoly.pm("(4*x1 - s)*x2^3*x3^2 +
	+ d31r*x2^3*x3 + d22r*x2^2*x3^2 + d21r*x2^2*x3 +
	+ d12r*x3^2*x2 + d11r*x2*x3 + d02r*x3^2 + d01r*x3 +
	# x2-Only:
	d50*x2^5 + d40*x2^4 + d30*x2^3 + d20*x2^2 + d10*x2 + d00");

pR = solve.pm(p1, p2, "x3", stop.at=1)
pR = pR[[2]];
isX3 = pR$x3 > 0;
pR = pR[isX3, ]; pR$x3 = 0; pR = drop.pm(pR);
toCoeff(pR, "x2", addC=T)

# pR = replace.pm(pR, toPoly.pm("s*x2^3 - e2*x2^2 + e3*x2 - e4"), "x2", pow=4);

# [requires the full coefficients]
(dq9*x2^9 + dq8*x2^8 + dq7*x2^7 + dq6*x2^6 + dq5*x2^5 + dq4*x2^4 + dq3*x2^3 + dq2*x2^2 + dq1*x2 + dq0)*x3 +
	v10*x2^10 + v9*x2^9 + v8*x2^8 + v7*x2^7 + v6*x2^6 + v5*x2^5 + v4*x2^4 + v3*x2^3 + v2*x2^2 + v1*x2 + v0 # = 0

###
xCoeff = coeff.S5HtMixed.x2(x1, cbind(s, e2, e3, e4), E11a, E11b);
v = xCoeff$v; dq = xCoeff$dq;
v0r = v[,1]; v1r = v[,2]; v2r = v[,3]; v3r = v[,4];
dq0r = dq[,1]; dq1r = dq[,2]; dq2r = dq[,3]; dq3r = dq[,4];


(dq3r*x2^3 + dq2r*x2^2 + dq1r*x2 + dq0r)*x3 +
	v3r*x2^3 + v2r*x2^2 + v1r*x2 + v0r # = 0

### Reduction:
px2 = toPoly.pm("s*x2^3 - e2*x2^2 + e3*x2 - e4")
p1 = toPoly.pm("dq9*x2^9 + dq8*x2^8 + dq7*x2^7 + dq6*x2^6 + dq5*x2^5 + dq4*x2^4 + dq3*x2^3 + dq2*x2^2 +
	+ dq1*x2 + dq0")
for(i in seq(5)) p1 = replace.pm(p1, px2, "x2", pow=4)


p2 = toPoly.pm("v10*x2^10 + v9*x2^9 + v8*x2^8 + v7*x2^7 + v6*x2^6 + v5*x2^5 + v4*x2^4 + v3*x2^3 + v2*x2^2 +
	+ v1*x2 + v0")
for(i in seq(5)) p2 = replace.pm(p2, px2, "x2", pow=4)

### Derivation x2:
# TODO

p1 = toPoly.pm("(dq3r*x2^3 + dq2r*x2^2 + dq1r*x2 + dq0r)*x3 +
	v3r*x2^3 + v2r*x2^2 + v1r*x2 + v0r");
p2 = toPoly.pm("2*x2*x3^3 - 2*x1*x3^3 + (4*x2^2 - 2*x2*x1 - 3*x2*s + b02r)*x3^2 +
	+ (3*x2^3 - x2^2*x1 - 4*x2^2*s + b11r*x2 + b01r)*x3 +
	- (s - x1)*x2^3 + b20r*x2^2 + b10r*x2 + b00r");

pR = solve.pm(p1, p2, xn="x3")
pR = pR$Rez;
str(pR)
# substantial:
cat(paste0("f", 0:12, " = ", rev(toCoeff(pR, "x2", print=F)), ";\n"))



p2 = toPoly.pm("f12*x2^12 + f11*x2^11 + f10*x2^10 + f9*x2^9 + f8*x2^8 + f7*x2^7 + f6*x2^6 + f5*x2^5 +
	+ f4*x2^4 + f3*x2^3 + f2*x2^2 + f1*x2 + f0")

for(i in seq(6)) p2 = replace.pm(p2, px2, "x2", pow=4)
#
cat(paste0("f", 0:3, "r = ", rev(toCoeff(p2, "x2", print=F)), ";\n"))

### Test:
err = f3r*x2^3 + f2r*x2^2 + f1r*x2 + f0r # = 0
# but massive overflow;
err / min(abs(c(f0r,f1r,f2r,f3r)))

###
p1 = toPoly.pm("f3r*x2^3 + f2r*x2^2 + f1r*x2 + f0r");
p2 = toPoly.pm("x2^4 - s*x2^3 + e2*x2^2 - e3*x2 + e4");
pR = solve.pm(p1, p2, xn="x2", stop.at=2)


