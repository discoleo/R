########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S5: Hetero-Symmetric
### Useful Formulas
###
### draft v.0.1a


### Derivations
# - Basic derivations;
# - Numerical approaches: basic intuition;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


#######################
#######################

### Base-System

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

