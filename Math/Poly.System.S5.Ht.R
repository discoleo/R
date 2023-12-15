########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S5
### Heterogeneous Symmetric
###
### draft v.0.1b


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
	-(9*R^3 - 40*R^2 + 28*R - 32) # = 0


### System Transform:

E2 = (S^2 - apply(x^2, 1, sum)) / 2;
E3 = - (S^3 - 3*E2*S - apply(x^3, 1, sum)) / 3;
E4 = (S^4 - 4*S^2*E2 + 2*E2^2 + 4*S*E3 - apply(x^4, 1, sum)) / 4;

### Eq 1:
# sum =>
S^2 + S - 2*E2 - 5*R # = 0

### Eq 2:
# sum(x2 * ...) =>
# E21a + S^2 - 2*E2 - R*S # = 0

# sum(x1^2 * ...) =>
# S^4 - 4*S^2*E2 + 2*E2^2 + 4*S*E3 - 4*E4 + E21a - R*(S^2 - 2*E2) # = 0
# =>
S^4 - 4*S^2*E2 + 2*E2^2 + 2*E2 + 4*S*E3 - 4*E4 - S^2 - R*(S^2 - S - 2*E2) # = 0

### Eq 3:
# sum(x1 * ...) =>
# S^3 - 3*E2*S + 3*E3 + E2a - R*S # = 0

# TODO

