########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: Fractions
### Derivations
###
### draft v.0.1e-example


### Polynomial Fractions
# - various Decompositions;


######################

### Helper Functions

source("Polynomials.Helper.R")


#######################
#######################


b = c(1,1,2,3,4)
p = toPoly.pm("x^5 + b[5]*x^4 + b[4]*x^3 + b[3]*x^2 + b[2]*x + b[1]");
r = roots.pm(p)
#
dp1 = dp.pm(p, by="x")
dp2 = dp.pm(dp1, by="x")


###############
### Order 1 ###
###############

### sum( 1/(x - r) )
# = d(Q(x)) / Q(x);

### Ex 1:
x = 3
sum(1/(x - r))
eval.pm(dp1, x) / eval.pm(p, x)

####################

### Decomposition of:
### x * d(Q(x))/Q(x)

### sum( r[i] / (x - r[i]) )
# = x * d(Q(x)) / Q(x) - n;

x = 3
n = length(r)
#
sum( (r - x) / (x - r) ) + sum( x / (x - r) ) # =
- n + sum( x / (x - r) ) # =
# x * d(Q(x)) / Q(x) - n;

### Ex 1:
x = 3
n = length(r)
sum( r / (x - r) )
x * eval.pm(dp1, x) / eval.pm(p, x) - n;

####################

### Decomposition of:
### x^2 * d(Q(x))/Q(x)

### sum( r[i]^2 / (x - r[i]) )
# = x^2 * d(Q(x)) / Q(x) - n*x - S;

x = 3
n = length(r)
S = round(sum(r), 5);
#
sum( (r^2 - x^2) / (x - r) ) + sum( x^2 / (x - r) ) # =
- n*x - S + x^2*sum( 1 / (x - r) ) # =
# x^2 * d(Q(x)) / Q(x) - n*x - S;

####################

### Decomposition of:
### x^3 * d(Q(x))/Q(x)

### sum( r[i]^3 / (x - r[i]) )
# = x^3 * d(Q(x)) / Q(x) - n*x^2 - S*x - S^2 + 2*E2*S;

x = 3
n = length(r)
S = round(sum(r), 5); # b[n]
E2 = b[n-1]; # sum(apply(combn(r, 2), 2, prod))
#
sum( (r^3 - x^3) / (x - r) ) + sum( x^3 / (x - r) ) # =
- n*x^2 - S*x - (S^2 - 2*E2) + x^3*sum( 1 / (x - r) ) # =
# x^3 * d(Q(x)) / Q(x) - n*x^2 - S*x - S^2 + 2*E2*S;

# Note:
# - Higher powers: sum( r^k / (x-r) )
#   are useful for decomposition into fractions of Order 1;
#   ["pure" decompositions]

### Examples

px = function(x) (x^5 - x^4 - 5*x - 1);
r = roots(c(1, -1,0,0,-5,-1))
# Test value:
x = 4
# D(px)
(5*x^4 - 4*x^3 - 5)/px(x)
sum(1 / (x - r))
#
(x^4 + 20*x + 5)/px(x)
sum(r/(x - r))
#
(x^4 + 20*x^2 + 10*x + 1)/px(x)
sum(r^2/(x - r))
#
(20*x^2 - 10*x - 4)/px(x)
sum((r^2 - r)/(x - r))
#
(x^4 + 20*x^3 + 10*x^2 + 6*x + 1)/px(x)
sum(r^3/(x - r))
#
(20*x^3 + 10*x^2 - 14*x - 4)/px(x)
sum((r^3 - r)/(x - r))
#
(21*x^4 + 10*x^3 + 6*x^2 + 6*x + 1)/px(x)
sum(r^4/(x - r))
# (a*x + b0) / px(x)
- 2*(509*x + 152)/px(x)
sum((2*r^3 - r^2 - 51*r + 10)/(x - r))
# (a*x^2 + b0) / px(x)
(20*509*x^2 - 516) / px(x)
sum((509*(r^2 - r) - 5*(2*r^3 - r^2 - 51*r + 10))/(x - r))
sum((- 10*r^3 + 514*r^2 - 254*r - 50)/(x - r))


###############

###############
### Order 2 ###
###############

### sum( 1/((x - r[i])*(x - r[j])) )
# = d2(Q(x)) / (2 * Q(x))

### Ex 1:
x  = 3
xd = (x - r);
xF2 = combn(xd, 2);
#
sum( apply(1/xF2, 2, prod) )
eval.pm(dp2, x) / (2 * eval.pm(p, x))


#####################

### Decomposition of:
### x * d2(Q(x))/Q(x)

### sum( (r[i] + r[j]) / ((x - r[i])*(x - r[j])) )
# = (x*d2(Q(x)) - 4*d(Q(x))) / Q(x);

# sum( (r[i] + r[j] - 2*x + 2*x) / ((x - r[i])*(x - r[j])) ) # =
# - 4*sum(1 / (x - r[i])) + 2*x*sum( 1/((x - r[i])*(x - r[j])) ) # =
# x * d2(Q(x)) / Q(x) - 4*d(Q(x)) / Q(x) # =
# (x*d2(Q(x)) - 4*d(Q(x))) / Q(x);

### Ex 1:
x  = 3
xd = (x - r);
xF2 = combn(xd, 2);
#
sum( apply(xF2, 2, function(tx) (2*x - sum(tx)) / prod(tx)) )
(x*eval.pm(dp2, x) - 4*eval.pm(dp1, x)) / eval.pm(p, x)


#####################

### Decomposition of:
### x^2 * d2(Q(x))/Q(x)

### sum( (r[i]*r[j]) / ((x - r[i])*(x - r[j])) )
# = choose(n, 2) + (1/2*x^2*d2(Q(x)) - 4*x*d(Q(x))) / Q(x);

# sum( ((x-r[i])*(x-r[j]) - x^2 + x*(r[i]+r[j])) / ((x - r[i])*(x - r[j])) ) # =
# choose(n, 2) - sum( (x^2 - x*(r[i]+r[j])) / ((x - r[i])*(x - r[j])) ) # =
# choose(n, 2) - 1/2*x^2 * d2(Q(x)) / Q(x) + x*(x*d2(Q(x)) - 4*d(Q(x))) / Q(x)
# choose(n, 2) + 1/2*x^2 * d2(Q(x)) / Q(x) - 4*x*d(Q(x)) / Q(x)

### Ex 1:
x  = 3
n  = length(r)
xd = (x - r);
xF2 = combn(xd, 2);
#
sum( apply(xF2, 2, function(tx) prod(x - tx) / prod(tx)) )
choose(n, 2) + (1/2*x^2*eval.pm(dp2, x) - 4*x*eval.pm(dp1, x)) / eval.pm(p, x)


#####################

### Decomposition of:
### x^2 * d2(Q(x))/Q(x)
# - Alternative using (r[i]^2 + r[j]^2);

### sum( (r[i]^2 + r[j]^2) / ((x - r[i])*(x - r[j])) )
# = (x^2*d2(Q(x)) - (3*x+S)*d(Q(x))) / Q(x) - n;

# sum( (r[i]^2 + r[j]^2 - 2*x^2 + 2*x^2) / ((x - r[i])*(x - r[j])) ) # =
# - sum( (4*x + S - r[i]) / (x - r[i]) ) + 2*x^2*sum( 1/((x - r[i])*(x - r[j])) ) # =
# - (3*x+S)*sum(1/(x - r[i])) - n + x^2*d2(Q(x)) / Q(x) # =
# - (3*x+S)*d(Q(x)) / Q(x) + x^2*d2(Q(x)) / Q(x) - n # =
# (x^2*d2(Q(x)) - (3*x+S)*d(Q(x))) / Q(x) - n;

### Ex 1:
x  = 3
n  = length(r)
S  = - b[5];
xd = (x - r);
idF2 = combn(seq(n), 2);
#
sum( apply(idF2, 2, function(id) sum(r[id]^2) / prod(xd[id])) )
(x^2*eval.pm(dp2, x) - (3*x+S)*eval.pm(dp1, x)) / eval.pm(p, x) - n;


####################
####################

### Fraction Decomposition of P[5]

px = function(x) x^5 - x - 1;
dp = function(x) 5*x^4 - 1;
d2p = function(x) 20*x^3;
r = roots(c(1,0,0,0,-1,-1));
n = length(r);
idF2 = combn(seq(n), 2);

x = 3
# Helper:
xF2 = matrix(r[idF2], nrow=2); xdF2 = matrix((x - r)[idF2], nrow=2);
S   = round(sum(r), 5);
T1  = sum( 1 / (x - r) );
T1s = sum( r / (x - r) ) + n;
T1s2 = sum( r^2 / (x - r) );
T2s = sum( sapply(seq(ncol(idF2)), function(id) sum(xF2[ , id]) / prod(xdF2[ , id])) )
T2  = sum( sapply(seq(ncol(idF2)), function(id) 1 / prod(xdF2[ , id])) )


### x*dp / P[5]
x * dp(x) / px(x);
(5*x^5 - x) / px(x)
(4*x + 5) / px(x) + 5;
sum( r / (x - r) ) + n;

### x*d2p / P[5]
# 20*x^4
x * d2p(x) / px(x);
20*x^4 / px(x)
T2s + 4*T1

### Full Decomposition:

### 1 / P[5]
# [for this particular polynomial]
# = sum( (r[i] + r[j]) / ((x - r[i])*(x - r[j])) ) / 4;
1 / px(x);
1 / (x^5 - x -1)
(T2s + 4*T1)/4 - T1
T2s/4

### x / P[5]
# [for this particular polynomial]
x / px(x);
x / (x^5 - x -1)
(T1s - 5/4*T2s - 5)/4

### x^2 / P[5]
# [for this particular polynomial]
# x^2 * d(Q(x)) / Q(x)
x^2 * dp(x) / px(x)
x^2 * (5*x^4 - 1) / (x^5 - x - 1)
T1s2 + n*x + S
#
# div.pm(toPoly.pm("x^2*(5*x^4-1)"), toPoly.pm("x^5 - x - 1"), "x")$Rem
(4*x^2 + 5*x) / px(x)
T1s2; # S = 0;
# x^2 / Q(x)
print("x^2 / Q(x)")
x^2 / px(x)
(T1s2 - 5/4*(T1s - 5/4*T2s - 5))/4
(4*T1s2 - 5*T1s + 25/4*T2s + 25)/16

### x^3 / P[5]
# [for this particular polynomial]
T2 = sum( sapply(seq(ncol(idF2)), function(id) 1 / prod(xdF2[ , id])) )
x^3 / px(x);
x^3 / (x^5 - x -1)
T2 / 10


####################
####################

### [old Derivations]
# - from file:
#   Integrals.Fractions.CardanPoly.R;


# TODO: verify alternating signs;

### T1
### sum( r[i] / (x - r[i]) )
# = (E1*x^(n-1) - 2*E2*x^(n-2) + 3*E3*x^(n-3) - 4*E4*x^(n-4) + ... + (-1)^n * (n-1)*E[n-1]*x + (-1)^(n+1) * n*E[n]) / Q(x)

### T2
### sum( (r[i1] + r[i2]) / ((x - r[i1])*(x - r[i2])) )
# = (1*(n-1)*E1*x^(n-2) - 2*(n-2)*E2*x^(n-3) + 3*(n-3)*E3*x^(n-4) + ... + (-1)^(n-1) * (n-2)*2*E[n-2]*x + (-1)^n * (n-1)*1*E[n-1]) / Q(x)
# = sum( (-1)^(i+1) * (n-i)*i*E[i]*x^(n-i-1) ) / Q(x);
# where i = 1:(n-1);

### sum( (r[i1] * r[i2]) / ((x - r[i1])*(x - r[i2])) )
# = (1*E1*x^(n-2) - 3*E2*x^(n-3) + 6*E3*x^(n-4) + ... + (-1)^(n-1) * ...*E[n-2]*x + (-1)^n * ...*E[n-1]) / Q(x)
# = sum( (-1)^(i+1) * i*(i+1)/2 * E[i]*x^(n-i-1) ) / Q(x);
# where i = 1:(n-1);

### T3
### sum( (r[i1] + r[i2] + r[i3]) / ((x - r[i1])*(x - r[i2])*(x - r[i3])) )
# = ( choose(n-1, 2)*1* E1*x^(n-3) - ... - 10*(n-5)*E[n-5]*x^3 + 6*(n-4)*E[n-4]*x^2 - 3*(n-3)*E[n-3]*x + 1*(n-2)*E[n-2]) / Q(x)
# = sum( (-1)^(n-i-1) * choose(i, 2)*(n-i)*E[n-i]*x^(i-2) ) / Q(x);
# where i = 2:(n-1); ### Note: [i] starts at 2 & is in reverse order (vs T2)!
sapply(5:10, function(n) {r = roots(c(rep(c(1,-1), n-2),-1)); mult.pfr(r=r, k=3);})

### [...] many more;


###################

# [old approach]
expand.idgrid = function(n, k) {
	# expand grid to compute various root combinations
	id.l = rep(list(1:n), k)
	id.gr = expand.grid(id.l)
	for(i in 1:(k-1)) {
	for(j in (i+1):k) {
		id.gr = id.gr[id.gr[,i] < id.gr[,j],]
	}
	}
	return(id.gr)
}
mult.pfr = function(r, k=2, type=1, pow=1) {
	# type: 0 => 1/..., 1 => sum(r[i])/...,
	#       k => prod(r[i])/...;
	n = length(r); id.all = 1:n
	gr = if(k == 1) matrix(id.all, ncol=1) else expand.idgrid(n, k=k);
	p = rep(0, n)
	rows = 1:nrow(gr);
	for(id in rows) {
		is.id = id.all %in% gr[id,]
		r.inv = r[ ! is.id]
		# print(r[is.id])
		p.m = rev(Poly(r.inv))
		if(type == 0) {
		} else if(type == 1) {
			p.m = p.m * sum(r[is.id]^pow) # Sum
		} else if(type == k) {
			p.m = p.m * prod(r[is.id]^pow) # Product
		}
		p.m = c(p.m, rep(0, n - length(p.m)))
		p = p + p.m
	}
	round0.p(p)
}

######################

### Root Combinatorics

### TODO: move to separate file;

coeff = c(1,2,3,4,5,6)
coeff = c(1,0,0,0,1,1)
coeff = c(1,0,0,0,0,1,1)
coeff = c(1,0,0,0,1,1,1)
r = roots(coeff)

### Test
poly.calc(r)

p = sapply(0:4, function(pow) print(mult.pfr(r, k=1, type=1, pow=pow)) )
p = sapply(0:4, function(pow) print(mult.pfr(r, k=2, type=1, pow=pow)) ) # Coeffs = 1 * D([k = 1])
p = sapply(0:4, function(pow) print(mult.pfr(r, k=3, type=1, pow=pow)) ) # Coeffs = 1/2 * D([k = 2])

### sum (r[i]^pow / (x - r[i]))
mult.pfr(r, k=1, type=1, pow=0) # D[1]
mult.pfr(r, k=1, type=1)
mult.pfr(r, k=1, type=1, pow=2)
mult.pfr(r, k=1, type=1, pow=3)


### sum ( (r[i]^pow + r[j]^pow) / ((x - r[i])*(x-r[j]) )
mult.pfr(r, k=2, type=1)
mult.pfr(r, k=2, type=1, pow=2)
mult.pfr(r, k=2, type=1, pow=3)


### sum ( (r[i] * r[j])^pow) / ((x - r[i])*(x-r[j]) )
mult.pfr(r, k=2, type=2)
mult.pfr(r, k=2, type=2, pow=2)
mult.pfr(r, k=2, type=2, pow=3)

