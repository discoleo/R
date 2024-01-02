########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S5
### Heterogeneous Symmetric
###
### draft v.0.1d


####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


test.S5.Simple = function(sol, R=NULL, b=1, n=2) {
	if(is.null(dim(sol))) sol = natrix(sol, nrow=1);
	x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3];
	x4 = sol[,4]; x5 = sol[,5];
	err1 = x1^n + b*x2 # - R
	err2 = x2^n + b*x3 # - R
	err3 = x3^n + b*x4 # - R
	err4 = x4^n + b*x5 # - R
	err5 = x5^n + b*x1 # - R
	err = rbind(err1, err2, err3, err4, err5);
	if( ! is.null(R)) err = err - R;
	err = round0(err);
	return(err);
}

### Solve

# x = x1 root;
solve.classic.xi = function(x, R, n = 2) {
	xi = cbind(x, 0, 0, 0, 0);
	for(i in seq(4)) {
		xi[, i + 1] = R - xi[, i]^n;
	}
	return(xi)
}

### Generate Classic Polynomial:
# assumes b = 1;
classic.S5Simple.gen = function(n = 2) {
	p1 = data.frame(
		x = c(n,0), R = c(0,1),
		coeff = c(-1, 1)
	);
	bR.gen = function(pR=1) data.frame(R = pR, coeff = 1);
	bx.gen = function(pb, px=1) data.frame(x = px, coeff = 1);
	for(i in seq(3)) {
		p1 = pow.pm(p1, n);
		p1 = diff.pm(bR.gen(pR=1), p1);
	}
	p1 = pow.pm(p1, n)
	p1 = diff.pm(p1, bR.gen(pR=1))
	p1 = add.pm(p1, bx.gen(px=1))
	if(n %% 2 == 1) p1$coeff = - p1$coeff;
	p1 = sort.pm(p1, sort.coeff=c(4,2,3,1), xn="x")
	rownames(p1) = seq(nrow(p1))
	return(p1);
}

######################
######################

### System:
# x1^n + x2 = R
# ...
# x4^n + x5 = R
# x5^n + x1 = R

#############
### n = 2 ###
#############


### Derivation:
pR = classic.S5Simple.gen(n = 2)
tmp = div.pm(pR, as.pm("x^2 + x - R"), by = "x")
tmp = tmp$Rez;

str(tmp)
# 256 monomials
# print.pm(tmp, leading = "x")


#########
### R = 1
R = 1
tmpR = replace.pm(tmp, c(R = R), xn = "R")

# roots.pm(tmpR)

x = solve.classic.xi(roots.pm(tmpR), R = R)
test.S5.Simple(x)

S = sort(apply(x, 1, sum))
round(poly.calc0(S[seq(1,30, by = 5)]))

err = S^6 - S^5 - 8*S^4 + 7*S^3 + 44*S^2 - 77*S + 35
round0(err)


##########
### R = -2
R = -2
tmpR = replace.pm(tmp, c(R = R), xn = "R")

# roots.pm(tmpR)

x = solve.classic.xi(roots.pm(tmpR), R = R)
test.S5.Simple(x)

S = apply(x, 1, sum)
id = order(S);
x = x[id, ]; S = S[id];
round(poly.calc0(S[seq(1,30, by = 5)]))

err = S^6 - S^5 + 25*S^4 - 47*S^3 + 158*S^2 - 56*S + 320
round0(err)

###
S^6 - S^5 - (11*R - 3)*S^4 + (18*R - 11)*S^3 +
	+ (19*R^2 - 19*R + 44)*S^2 - (17*R^2 + 24*R + 36)*S +
	- (9*R^3 - 40*R^2 + 28*R - 32) # = 0


### System Transform:

E2 = (S^2 - apply(x^2, 1, sum)) / 2;
E3 = - (S^3 - 3*E2*S - apply(x^3, 1, sum)) / 3;
E4 = (S^4 - 4*S^2*E2 + 2*E2^2 + 4*S*E3 - apply(x^4, 1, sum)) / 4;
E5 = apply(x, 1, prod);
#
E2a = - (S^3 - 3*E2*S + 3*E3 - R*S);

### Eq 1:
# sum =>
S^2 + S - 2*E2 - 5*R # = 0

### Eq 2:
# sum(x2 * ...) =>
# E21a + S^2 - 2*E2 - R*S # = 0

# sum(x1^2 * ...) =>
# S^4 - 4*S^2*E2 + 2*E2^2 + 4*S*E3 - 4*E4 + E21a - R*(S^2 - 2*E2) # = 0
# =>
S^4 - 4*S^2*E2 + 2*E2^2 + 2*E2 + 4*E3*S - 4*E4 - S^2 - R*(S^2 - S - 2*E2) # = 0
# Simplification =>
S^4 + 2*S^3 - 10*R*S^2 - S^2 + 6*R*S - 2*S - 8*E3*S + 8*E4 - 15*R^2 + 10*R # = 0


### Eq 3:
S^5 - 5*E2*S^3 + S^3 + 5*E3*S^2 + 5*E2^2*S - 5*E4*S - 3*E2*S +
	- 5*E2*E3 + 5*E5 + 3*E3 + E2 +
	- 2*R*(S^3 - 3*E2*S + S + 3*E3) + R^2*S # = 0
# Simplification =>
S^5 - 3*S^3 - 4*R*S^3 + 4*S^2 + 38*R*S^2 - 10*E3*S^2 - 2*S - 22*R*S - 69*R^2*S + 10*E3*S + 20*E4*S +
	+ 10*R - 12*E3 - 26*R*E3 - 20*E5 # = 0
S^3 + 3*R*S^3 + 2*S^2 + 6*R*S^2 - E3*S^2 - 3*S - 10*R*S - 27*R^2*S - 3*E3*S + 6*E4*S +
	- 10*E5 + 8*E4 - 6*E3 - 13*R*E3 - 15*R^2 + 15*R # = 0


### Eq 4:
3*S^5 + 6*S^4 - 13*R*S^4 - 3*S^3 - 56*R*S^3 - 6*S^2 + 13*R*S^2 + 130*R^2*S^2 - 12*E4*S^2 - 24*E3*S^2 +
	+ 38*R*S - 51*R^2*S + 12*E4*S + 104*R*E3*S + 24*E5*S +
	- 44*R*E4 + 12*E3 + 12*E3^2 - 40*R^2 + 75*R^3 # = 0


### Eq 5:

# TODO:


###
# sum(x1 * ...) =>
# S^3 - 3*E2*S + 3*E3 + E2a - R*S # = 0

# sum(x5 * ...) =>
# E21c + E2b - R*S # = 0

# sum(x2^2 * ...) =>
# E22a + S^3 - 3*E2*S + 3*E3 - R*(S^2 - 2*E2) # = 0

# sum(x3^2 * ...) =>
# E22b + E21c - R*(S^2 - 2*E2) # = 0

# sum(x5^2 * ...) =>
# E22a + E21b - R*(S^2 - 2*E2) # = 0

# =>
E22a + E21b - R*(S^2 - 2*E2) +
	- (E22a + S^3 - 3*E2*S + 3*E3 - R*(S^2 - 2*E2)) # = 0

# sum(x1^3 * ...) =>
S^5 - 5*S^3*E2 + 5*S*E2^2 + 5*S^2*E3 - 5*E2*E3 - 5*S*E4 + 5*E5 +
	+ E31a - R*(S^3 - 3*E2*S + 3*E3) # = 0

# sum(x1*x2* ...) =>
E31a + E21c - R*E2a # = 0

# =>
# sum((x1^3 - x1*x2 + x5) * ...)
S^5 - 5*S^3*E2 + 5*S*E2^2 + 5*S^2*E3 - 5*E2*E3 - 5*S*E4 + 5*E5 +
	- R*(S^3 - 3*E2*S + 3*E3) + R*E2a + (E2 - E2a) - R*S # = 0
#
S^5 - 5*E2*S^3 + S^3 + 5*E3*S^2 + 5*E2^2*S - 5*E4*S - 3*E2*S +
	- 5*E2*E3 + 5*E5 + 3*E3 + E2 +
	- 2*R*(S^3 - 3*E2*S + S + 3*E3) + R^2*S # = 0

# sum(x1^4 * ...) =>
# S^6 - 6*S^4*E2 + 9*S^2*E2^2 - 2*E2^3 + 6*S^3*E3 - 12*S*E2*E3 + 3*E3^2 - 6*S^2*E4 + 6*E2*E4 + 6*S*E5 +
#   + E41a - R*(S^4 - 4*S^2*E2 + 2*E2^2 + 4*S*E3 - 4*E4) # = 0

# sum(x1^2*x2* ...) =>
# E41a + E22a - R*E21a # = 0
# E41a + (R*(S^2 - 2*E2) - S^3 + 3*E2*S - 3*E3) - R*(R*S - S^2 + 2*E2) # = 0

# =>
S^6 - 6*S^4*E2 + 9*S^2*E2^2 - 2*E2^3 + 6*S^3*E3 - 12*S*E2*E3 + 3*E3^2 - 6*S^2*E4 + 6*E2*E4 + 6*S*E5 +
	- ((R*(S^2 - 2*E2) - S^3 + 3*E2*S - 3*E3) - R*(R*S - S^2 + 2*E2)) +
	- R*(S^4 - 4*S^2*E2 + 2*E2^2 + 4*S*E3 - 4*E4) # = 0
S^6 - 6*S^4*E2 + 9*S^2*E2^2 - 2*E2^3 + 6*S^3*E3 - 12*S*E2*E3 + 3*E3^2 - 6*S^2*E4 + 6*E2*E4 + 6*S*E5 +
	+ S^3 - 3*E2*S + 3*E3 +
	- R*(S^4 - 4*S^2*E2 + 2*E2^2 + 4*S*E3 - 4*E4 + 2*S^2 - 4*E2 - R*S) # = 0
# Simplification =>
3*S^5 + 6*S^4 - 13*R*S^4 - 3*S^3 - 56*R*S^3 - 6*S^2 + 13*R*S^2 + 130*R^2*S^2 - 12*E4*S^2 - 24*E3*S^2 +
	+ 38*R*S - 51*R^2*S + 12*E4*S + 104*R*E3*S + 24*E5*S +
	- 44*R*E4 + 12*E3 + 12*E3^2 - 40*R^2 + 75*R^3 # = 0


# more manipulations:
# - does NOT seem redundant;
p3 = toPoly.pm("S^6 + 5*S^5 + 3*S^4 - 21*R*S^4 - 9*S^3 - 46*R*S^3 +
	+ 12*S^2 + 35*R*S^2 + 71*R^2*S^2 - 40*E4*S^2 +
	- 12*S - 80*R*S + 153*R^2*S - 40*E4*S + 80*E5*S + 60*R + 40*R^2 - 195*R^3 + 48*E4 + 104*R*E4")
p4 = toPoly.pm("S^8 + 4*S^7 + 2*S^6 - 20*R*S^6 - 28*R*S^5 + 9*S^4 - 32*R*S^4 + 70*R^2*S^4 - 48*E4*S^4 +
	- 4*S^3 - 108*R*S^3 + 204*R^2*S^3 - 32*E4*S^3 + 128*E5*S^3 +
	- 12*S^2 + 4*R*S^2 + 346*R^2*S^2 - 340*R^3*S^2 - 16*E4*S^2 + 160*R*E4*S^2 +
	+ 40*R*S + 60*R^2*S - 180*R^3*S + 32*E4*S + 96*R*E4*S +
	+ 100*R^2 - 300*R^3 + 225*R^4 + 160*R*E4 - 240*R^2*E4 + 64*E4^2")
tmp = toPoly.pm("(S^2 - S + 4 + R)*p3() - p4()");
tmp$coeff = tmp$coeff / 4;
print.pm(tmp, leading = "S")
print.pm(toPoly.pm("2*(R+1)*p3() - S*tmp()"), leading = "S")


### Cyclic Redundancy:
# =>
E22a + S^3 - 3*E2*S + 3*E3 - R*(S^2 - 2*E2) +
	+ E22b + E21c - R*(S^2 - 2*E2) - (E21c + E2b - R*S) +
	- (S^3 - 3*E2*S + 3*E3 + E2a - R*S) # = 0
#
E2^2 - 2*E3*S + 2*E4 - E2 - 2*R*(S^2 - S - 2*E2) # = 0
