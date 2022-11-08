########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S5: Hetero-Symmetric
### Useful Formulas
###
### draft v.0.1a-E3


### Formulas:

# - Formulas & derivations;
# - Useful for S5 Ht Systems;

# - Applicable for systems described in:
#   TODO


### Sections

### Basic:
# A.) E11a
# B.) Higher Powers


####################

### Helper Functions


source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


### Debug
x = sqrt(c(2,3,5,7,11));
x[1] = - x[1]; x[5] = - x[5];
x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4]; x5 = x[5];

### Notation:
S  = sum(x);
E5 = prod(x);

s1 = x1 + x3; s2 = x2 + x4;
p1 = x1 * x3; p2 = x2 * x4;
ps = s1 * s2; sp = p1 + p2;
# E2 = x1*(S - x1) + x2*(x3 + x4 + x5) + x3*(x4 + x5) + x4*x5;
E2 = sp + ps + x5*(S - x5);
# E3 = x1*x2*(x3 + x4 + x5) + x3*x4*(x1 + x2 + x5) + (x1 + x2)*(x3 + x4)*x5;
E3 = p1*s2 + p2*s1 + x5*(sp + ps);
# E4 = x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5;
E4 = p1*p2 + x5*(p1*s2 + p2*s1);

### E2:
E11a = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1;
E11b = x1*x3 + x2*x4 + x3*x5 + x4*x1 + x5*x2;

E11a + E11b - E2 # = 0

### E3:
E111a = x1*x2*x3 + x2*x3*x4 + x3*x4*x5 + x4*x5*x1 + x5*x1*x2;
E111b = x1*x2*x4 + x2*x3*x5 + x3*x4*x1 + x4*x5*x2 + x5*x1*x3;

E111a + E111b - E3 # = 0

### Note:
# - the remaining cyclic permutations (E2) equal E11b & E11a;
#   Perm(S5, by = 3) = rev(E11b) = E11b;
#   Perm(S5, by = 4) = rev(E11a) = E11a;


### Mixed C2 & C3-Decomposition

s1 = x1 + x2; p1 = x1 * x2;
s2 = S - s1; e2 = (x3 + x4)*x5 + x3*x4; e3 = x3*x4*x5;

### Transformed P[5] System:
s1 + s2 - S # = 0
s1*s2 + p1 + e2 - E2 # = 0
s1*e2 + s2*p1 + e3 - E3 # = 0
s1*e3 + p1*e2 - E4 # = 0
p1*e3 - E5 # = 0

### Note:
# Order: 5! / (2*6) = 120/12 = 10;
# - much better than 5!, but still not enough;


#######################

A1 = E11a*E111a + E11b*E111b;
B1 = E11b*E111a + E11a*E111b;

### Sum:
A1 + B1 - E2*E3 # = 0

### Mult:
A1 * B1 - E11a*E11b*(E111a^2 + E111b^2) - E111a*E111b*(E11a^2 + E11b^2) # = 0
A1 * B1 - E11a*E11b*(E3^2 - 2*E111a*E111b) - E111a*E111b*(E2^2 - 2*E11a*E11b) # = 0
A1 * B1 - (E11a*E11b*E3^2 + E111a*E111b*E2^2) + 4*E11a*E11b*E111a*E111b # = 0

# TODO:
# - unfortunately is not yet symmetric;
#   E4321 = 60 (not 5! = 120);
#   E43111 = 5*4; OK
#   E42211 = 20 w coeff=3; 5 * choose(4,2) = 30;
#   E42211 = 10 w coeff=2; 5 * choose(4,2) = 30;
E11a*E11b*E111a*E111b


pE2a = toPoly.pm("x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1")
pE2b = toPoly.pm("x1*x3 + x2*x4 + x3*x5 + x4*x1 + x5*x2")
pE3a = toPoly.pm("x1*x2*x3 + x2*x3*x4 + x3*x4*x5 + x4*x5*x1 + x5*x1*x2")
pE3b = toPoly.pm("x1*x2*x4 + x2*x3*x5 + x3*x4*x1 + x4*x5*x2 + x5*x1*x3")
pR = prod.pm(pE2a, pE2b, pE3a, pE3b)

pc = countMonoms(pR)
pc = sort.pm(pc, "V1")
pc

# TODO
A2 = E11a*E11b*E3^2 + E111a*E111b*E2^2;
B2 = E11a*E11b*E2^2 + E111a*E111b*E3^2;

### Sum =>
A2 + B2 - (E11a*E11b + E111a*E111b)*(E2^2 + E3^2) # = 0


#######################

### Disaster:
pE1 = toPoly.pm("x1 + x2 + x3 + x4 + x5 - 1"); # S = 1; !!!
pE2a = toPoly.pm("x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1"); # E11a = 0; !!!
pE2b = toPoly.pm("x1*x3 + x2*x4 + x3*x5 + x4*x1 + x5*x2"); # E11b = 0; !!!
# E3 = 0; E4 = 0; !
pE3 = toPoly.pm("x1*x2*(x3 + x4 + x5) + x3*x4*(x1 + x2 + x5) + (x1 + x2)*x4*x5 + x3*x5*(x1 + x2)");
pE4 = toPoly.pm("x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5");

pR = solve.lpm(pE1, pE2a, pE2b, pE3, pE4, xn=c("x5", "x4", "x3", "x2"))


### Overflows!
pE1 = toPoly.pm("s1 + s2 + x5 - 1"); # S = 1; !!!
pE2 = toPoly.pm("p1 + p2 + s1*s2 + x5*(1 - x5)"); # E2 = 0; S = 1;
pE3 = toPoly.pm("p1*s2 + p2*s1 + x5*(s1*s2 + p1 + p2)");
pE4 = toPoly.pm("p1*p2 + x5*(p1*s2 + p2*s1)");
pE5 = toPoly.pm("p1*p2*x5 - 1"); # E5 = 1;

pR = solve.lpm(pE1, pE5, pE2, pE3, pE4, xn=c("x5", "s2", "s1", "p2"))
max(pR[[4]]$Rez$p1)


####################

### Standard P[5]

### quasi-redundant
# spp = sp + ps; pp4 = p1*p2;
# - Solvable, but just ordinary P[5];
pE2 = toPoly.pm("spp + x5*(S - x5) - E2");
pE3 = toPoly.pm("A + x5*spp - E3");
pE4 = toPoly.pm("pp4 + x5*A - E4");
pE5 = toPoly.pm("pp4*x5 - E5");

pE4 = toPoly.pm("x5^2*(E3 - x5*spp) + E5 - E4*x5")

pR = solve.pm(pE2, pE4, xn=c("spp"))
pR = pR$Rez;
print.pm(pR, lead="x5") # trivial: ordinary P[5];


#####################

### Test

# any permutation:
x = 2*cos(2*pi* c(1,5,3,4,2) /11)
# *OR*
m = unity(5, all=T)
k = rootn(2, 5); x = sapply(m, function(m) (k*m)^3 - (k*m)^2 + 3*k*m)
#
x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4]; x5 = x[5];
S = sum(x);
s1 = x1 + x2; p1 = x1 * x2;
s2 = S - s1; e2 = (x3 + x4)*x5 + x3*x4; e3 = x3*x4*x5;

x12 = roots(c(1, -s1, p1));
x35 = roots(c(1, -s2, e2, -e3));

# the same:
poly.calc(c(x12, x35))
poly.calc(x)


##########################
##########################

# - a great discovery in mathematics;
-27 + 109*x^3 + 189*x^5 - 114*x^6 - 654*x^8 - 110*x^9 - 567*x^10 + 570*x^11 + 355*x^12 + 1635*x^13 +  
440*x^14 + 3710*x^15 - 1140*x^16 - 1065*x^17 + 1064*x^18 - 660*x^19 - 6475*x^20 - 1984*x^21 + 1065*x^22 -  
1609*x^23 + 440*x^24 + 3332*x^25 - 570*x^26 - 355*x^27 - 654*x^28 - 110*x^29 - 189*x^30 + 114*x^31 +  
109*x^33 + 27*x^35 # = 0
# TODO: work out S, E3, E4;
27*(E11a^7 + E11b^7)*E5^2 + 81*E11a*E11b*(E11a^5 + E11b^5)*E5^2 + 27*(E11a*E11b)^2*(E11a^3 + E11b^3)*E5^2 +
	- 10*(E11a*E11b)^3*E2*E5^2 - (E11a*E11b)^6 - 54*(E11a*E11b)^5*(E11a^1 + E11b^1)*S^2 +
	# TODO: + ...
	- 5^4*E5^4*(E11b^2 + 3*E11a*E11b + E11a^2) # = 0

# some R-values are still fixed!
solve.S3HtMixed = function(R = c(0,1,0,0,1), debug=TRUE) {
	E11a = R[2]; E5 = R[5];
	# coeff = c(27, 80, 27, -10, -10, -3098, -9294, -3098);
	coeff = c(27, 81*E11a - E11a^6/E5^2, 27*E11a^2, -10*E11a^3, -10*E11a^4,
		-3125*E5^2 + 27*E11a^5, -9375*E11a*E5^2 + 81*E11a^6, -3125*E11a^2*E5^2 + 27*E11a^7);
	E11b = roots(coeff);
	if(debug) print(E11b);
	E2 = R[2] + E11b;
	x1 = sapply(E2, function(E2) roots(c(1, - R[1], E2, -R[3], R[4], -R[5])));
	x1 = as.vector(x1);
	# Robust:
	# TODO
	return(x1);
}
solve.S3HtMixed.Classic = function() {
	coeff = c(27, 0, 109, 0, 114, -189, -110, -654, -355, -570, 3332, 440, -1609, 1065, -1984, -6475, -660,
		1064, -1065, -1140, 3710, 440, 1635, 355, 570, -567, -110, -654, 0, -114, 189, 0, 109, 0, 0, -27);
	x1 = roots(coeff);
	# TODO:
	return(x1)
}

### Numerical Solution:
# library(rootSolve)
solve.S5HtMixed.Num = function(x, R=c(0,1,0,0,1)) {
	x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
	x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4]; x5 = x[5];
	s1 = x1 + x3; s2 = x2 + x4; S = s1 + s2 + x5;
	p1 = x1 * x3; p2 = x2 * x4; E5 = p1 * p2 * x5;
	ps = s1 * s2; sp = p1 + p2;
	# E2 = x1*(S - x1) + x2*(x3 + x4 + x5) + x3*(x4 + x5) + x4*x5;
	# E2 = sp + ps + x5*(S - x5);
	E3 = p1*s2 + p2*s1 + x5*(sp + ps);
	E4 = p1*p2 + x5*(p1*s2 + p2*s1);

	### E2:
	E11a = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1;
	#
	y = c(S, E11a, E3, E4, E5) - R;
	y = rbind(Re(y), Im(y));
	return(y);
}
cat.sol = function(x) {
	cat(paste0(round(x, digits=4), collapse=", ")); cat("\n");
}

### Set 1:
x0 = c(-0.70449+0.64i, 0.8913, -0.70449-0.64i, 0.2589 + 1.08i, 0.2589 - 1.08i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x)


### Set 2:
x0 = c(-1.70449-0.64i, -0.8913, -1.70449+0.64i, -0.2589 + 1.08i, -0.2589 - 1.08i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 3:
x0 = c(0.05-0.94i, -2.8913, 0.05+0.04i, -0.2589 + 1.08i, -0.2589 - 1.08i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 4:
x0 = c(-1.274+0.729i, 1.23-0.47i, 0.58-0.67i, 0.178+0.725i, -0.71 - 0.31i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 5:
x0 = c(-1.274-0.729i, 1.23+0.47i, 0.58+0.67i, 0.178-0.725i, -0.71 + 0.31i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 6:
x0 = c(-0.8108+1.5014i, 0.7763-1.6039i, 0.6605-0.1778i, -0.5008-0.4342i, -0.1252+0.7145i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 7:
x0 = c(-0.8108-1.5014i, 0.7763+1.6039i, 0.6605+0.1778i, -0.5008+0.4342i, -0.1252-0.7145i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)
x.all = matrix(x.all, nc=5, byrow=T)

round0(poly.calc(x.all)) * 27
poly.calc(apply(x.all, 1, function(x) sum(x * x[c(3,4,5,1,2)]))) * 27


#################

### Case 2:

### Set 1:
R2 = c(0,1,0,0,2)
x0 = c(-0.8392+0.7287i, 1.051, -0.8392-0.7287i, 0.3137+1.2009i, 0.3137-1.2009i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0, R=R2)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x)

### Set 2:
x0 = c(-0.42-0.6586i, 0.7387+0i, -0.42+0.6586i, 0.0507+2.1058i, 0.0507-2.1058i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0, R=R2)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x.all, x)

### Set 3:
x0 = c(0.392+0.8686i, 1.4881+0i, 0.392-0.8686i, -1.136-0.4352i, -1.136+0.4352i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0, R=R2)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x.all, x)

### Set 4:
x0 = c(-0.1185+0.8347i, -0.591-0.4927i, 0.7697-0.2215i, 0.9223-1.7751i, -0.9826+1.6545i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0, R=R2)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)

### Set 5:
x0 = c(-0.1185-0.8347i, -0.591+0.4927i, 0.7697+0.2215i, 0.9223+1.7751i, -0.9826-1.6545i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0, R=R2)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)

### Set 6:
x0 = c(0.218+0.8128i, -0.8281-0.33i, -1.526+0.7887i, 1.5164-0.5302i, 0.6196-0.7413i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0, R=R2)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)

### Set 7:
x0 = c(0.218-0.8128i, -0.8281+0.33i, -1.526-0.7887i, 1.5164+0.5302i, 0.6196+0.7413i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0, R=R2)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)
x.all = matrix(x.all, nc=5, byrow=T)

round0(poly.calc(x.all)) * 27
poly.calc(apply(x.all, 1, function(x) sum(x * x[c(3,4,5,1,2)]))) * 27


-12473 - 37419*x - 12473*x^2 - 10*x^3 - 10*x^4 + 27*x^5 + 80.75*x^6 + 27*x^7


####################


x = x.all[1,]; E3 = E4 = 0;
x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4]; x5 = x[5];
s1 = x1 + x2; p1 = x1 * x2;
s2 = x3 + x4 + x5; e2 = (x3 + x4)*x5 + x3*x4; e3 = x3*x4*x5;
S = s1 + s2; E5 = p1*e3;
E11a = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1;
E11b = x1*x3 + x2*x4 + x3*x5 + x4*x1 + x5*x2;
E2 = E11a + E11b;


###
### Transformed P[5] System:
s1 + s2 - S # = 0
#
s1*S - s1^2 + p1 + e2 - E2 # = 0
s1*e2 - s1*p1 + S*p1 + e3 - E3 # = 0
s1*e3 + p1*e2 - E4 # = 0
p1*e3 - E5 # = 0
x5^3 + s1*x5^2 + e2*x5 - e3 - S*x5^2 # = 0
#
s1*S - s1^2 + p1 + e2 - E2 # = 0
s1*p1*e2 - s1*p1^2 + S*p1^2 - p1*E3 + E5 # = 0
s1*E5 + p1^2*e2 - p1*E4 # = 0
p1*x5^3 + p1*s1*x5^2 + p1*e2*x5 - p1*S*x5^2 - E5 # = 0


p1 = toPoly.pm("s1*S - s1^2 + p1 + e2 - E2")
p2 = toPoly.pm("s1*p1*e2 - s1*p1^2 + S*p1^2 - p1*E3 + E5")
p3 = toPoly.pm("s1*E5 + p1^2*e2 - p1*E4")
p4 = toPoly.pm("p1*x5^3 + p1*s1*x5^2 + p1*e2*x5 - p1*S*x5^2 - E5")

pR1 = solve.lpm(p1, p4, p2, xn=c("p1", "e2"))
pR2 = solve.lpm(p1, p4, p3, xn=c("p1", "e2"))
pR1 = pR1[[2]]$Rez; pR1$coeff = - pR1$coeff;
pR2 = pR2[[2]]$Rez; pR2$coeff = - pR2$coeff;
table(pR2$s1)

tmp = gcd.pm(pR1, pR2, by="s1")
pR2 = diff.pm(pR2, mult.pm(pR1, toPoly.pm("x5^3")))

# Note: coeff a == 0!
x5^2*(S - x5)*(x5^5 - S*x5^4 + E2*x5^3 - E3*x5^2 + E4*x5 - E5)*s1^2 +
	- x5^2*(S - x5)^2*(x5^5 - S*x5^4 + E2*x5^3 - E3*x5^2 + E4*x5 - E5)*s1 +
	- E5^2 + 2*E4*E5*x5 - E4^2*x5^2 - E5*E2*S*x5^2 + E5*E2*x5^3 + E4*E2*S*x5^3 + E5*S^2*x5^3 - E4*E2*x5^4 +
	- 2*E5*S*x5^4 - E4*S^2*x5^4 + E5*x5^5 + 2*E4*S*x5^5 + E2^2*S*x5^5 - E4*x5^6 - 2*E2*S^2*x5^6 +
	+ 2*E2*S*x5^7 + S^3*x5^7 - 2*S^2*x5^8 + S*x5^9 - E2*S*x5^4*E3 - E2*x5^5*E3 + S^2*x5^5*E3 +
	- x5^7*E3 + x5^4*E3^2 # = 0

