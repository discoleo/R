###
### Fraction Rationalisation
###
### Leonard Mada
### 2018-2020
###
### draft v.0.2


### Rationalizing polynomial Fractions
### with Simple Radicals


###############
### Example ###

# a simple example
k = 3^(1/5) # this is here a fixed value!
1 / (k^4 + k^3 + k^2 + 3*k + 1) # ==
(3*k^4 + k^3 - 2*k^2 - 3*k - 1) / 14
# == 0.09406944
# for other examples: see below;


#####################

###############
### History ###
###############

### draft v.0.2:
# - major simplification;
### draft v.0.1:
# - minor edits & bug-fixes;
# draft pre-v.0.1.
# - initial upload to Github;


##############
### Theory ###

### Terminology:
# K, s[j] = parameters, e.g. integers;
# k = K^(1/n);
# x = sum(s[j] * k^j), with j = 0 to (n-1);
# The fraction 1/x can be rationalised.
# 1/x = sum(b[j] * k^j) / N,
# with j = 0 to (n-1);
# if K and s[j] are integers,
# => then b[j] and N are also integers;
# (N is chosen so that coefficients b[j] are integers)

# - various Examples are provided below;
# - the function inverse() computes the rational decomposition of 1/x;
# - if the initial coefficients are integers,
#   then the computed coefficients are also integers;


### Major Simplification:
# - skips solving the Vandermond matrix;

# 1.) 1/x can be rationalized;
# 2.) 1/x = prod(x[j]) / N;
# where:
#   x[j] = sum(s[i] * k[j]^i),
#   k[j] = k * m^j, j = seq(1, n-1), i = seq over s;
#   m^n  = 1, m is a root of unity of order n;
#   N = prod(x[j]), for j = seq(0, n-1);

### "Cyclic paradox"
# Note: the definition/computation should *NOT* be viewed as cyclic;
# - N is an integer, when s[i] and K are integers;
# - N can be actually computed efficiently (both parametrically and non-parametrically);
# - the actual computation of prod(r[j]) requires the symbolic multiplication of the (n-1) roots;


####################
####################

### Helper Functions

### Symbolic computations
source("Polynomials.Helper.R")


### Symbolic

# explicit roots
roots.Class1 = function(s, n, type="Complex") {
	# K = parameter
	len = length(s);
	if(len > n) stop("Too many parameters!");
	# TODO: All;
	type = pmatch(type, c("Complex", "Real", "All"));
	if(is.na(type)) stop("Unsupported type!");
	#
	isK = (s != 0);
	id = seq(0, len-1)[isK];
	if(type == 2) {
		r0 = data.frame(k = id, coeff=s[isK]);
		return(toPoly.pm(r0));
	}
	r0 = data.frame(k = id, m = id, coeff=s[isK])
	r0 = toPoly.pm(r0);
	#
	r = lapply(seq(n-1), function(id) {
		r = r0;
		r$m = (r$m * id) %% n;
		return(r);
	})
	return(r)
}


### Old (explicit) approach
inverse = function(p.coeffs, K, n=length(p.coeffs), tol=1E-10) {
	# p.coeffs = in ascending power order;
	# Base-Root & All Roots
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	pow = 0:(n-1)
	m = m^pow
	k = m * K^(1/5) # used to construct all roots;
	#
	if(length(p.coeffs) < n + 1) {
		p.coeffs = c(p.coeffs, rep(0, n - length(p.coeffs)))
	}
	x = sapply(1:n, function(id) sum(p.coeffs * k[id]^pow))
	N = prod(x) # Integer coeffs
	# multiply polynomials
	p.m = outer(p.coeffs, rep(1, n)) # ascending power
	#
	m.id = row(p.m) + col(p.m)
	p.m = p.m * k[1]^(m.id - 2)
	# print(p.m)
	id = ((m.id - 2) %% n) * n + col(p.m)
	# print(id)
	p.m = matrix(p.m[match(1:(n*n), id)], ncol=n, byrow=T)
	sol = solve(p.m, c(N, rep(0,4)))
	if( ! is.na(tol)) {
		sol = round0(sol, tol=tol)
	}
	# TODO: gcd(N, sol)
	# Test
	test = sapply(1:n, function(id) sum(sol * k[id]^pow) / N * x[id]) # == 1
	return(list(sol=sol, N=N, r.coeffs=p.coeffs, x=x, test=test))
}

################

################
### Examples ###


### Symbolic computation
n = 5;
s = c(1,1,-3,0,1); # c(s0, s1, s2, s3, ...)
r = roots.Class1(s, n=n)
# Prod:
p = prod.pm(r)
# Reduce K & m:
p = replace.pm(p, "K", "k", pow=n)
p = replace.pm(p, 1, "m", pow=n)
mZeroSum = data.frame(m=seq(0, n-2), coeff=-1);
p = replace.pm(p, mZeroSum, "m", pow=n-1);
# Note: *not* all k have canceled out,
# as the product skipped r[0];
print.pm(p, lead="k")
# Test
K = 3
rinv = eval.pm(p, c(K=K, k = K^(1/n)))
r0 = roots.Class1(s, n=n, type="Real")
r0 = eval.pm(r0, c(k = K^(1/n)));
# Test:
rinv / round(rinv * r0) # TODO: symbolic;
1 / r0;


### [old]

### n = 5
n = 5
# for *ALL* "roots"
m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
m = m^(0:(n-1))

#############
### Example 1
K = 3
k = m * K^(1/n)
# Note:
# - coefficients tend easy to explode;
# - examples have been choosen with smaller coefficients;
x = k^4 + k^3 + k^2 + 3*k + 1
### Rationalize: 1/x
1/x


### Solution:
sol = inverse(c(1, 3, 1, 1, 1), K, n)
sol

### Test
sapply(1:n, function(id) sum(sol$sol * k[id]^(0:(n-1))) / sol$N * x[id]) # == 1

# 1/x = (3*k^4 + k^3 - 2*k^2 - 3*k - 1) / 14 # N = 14 * 8
1/x # ==
(3*k^4 + k^3 - 2*k^2 - 3*k - 1) / 14

# Proof:
# (k^4 + k^3 + k^2 + 3*k + 1) * (3*k^4 + k^3 - 2*k^2 - 3*k - 1) = 14
# x * 1/x * 14 =
3*k^8 + 4*k^7 + 2*k^6 + 5*k^5 − 9*k^3 − 12*k^2 − 6*k − 1 # = 14
(k^5 - 3) * (5 + 2*k + 4*k^2 + 3*k^3) + 14
# Substituting: k^5 = 3
# = 15 - 1 = 14

### Long Procedure
# Polynomial Divisions: k^5 - K = 0
1/x # =
1 / (k^4 + k^3 + k^2 + 3*k + 1)
(k - 1) / ((k - 1) * (k^4 + k^3 + k^2 + 3*k + 1))
(k - 1) / (k^5 + 2*k^2 − 2*k − 1)
(k - 1) / (2*k^2 − 2*k − 1 + K) # K = 3
(k - 1) / (2*k^2 − 2*k + 2)
(k - 1) / (k^2 − k + 1) / 2 # simplifies the coefficients
(k - 1)*(k^3 + k^2 - 1) / ((k^3 + k^2 - 1)*(k^2 − k + 1)) / 2
(k - 1)*(k^3 + k^2 - 1) / (k^5 + k − 1) / 2
(k - 1)*(k^3 + k^2 - 1) / (k − 1 + K) / 2
(k - 1)*(k^3 + k^2 - 1) / (k + 2) / 2
(k - 1)*(k^3 + k^2 - 1)*(k^4 - 2*k^3 + 4*k^2 - 8*k + 16) / ((k^4 - 2*k^3 + 4*k^2 - 8*k + 16)*(k + 2)) / 2
(k - 1)*(k^3 + k^2 - 1)*(k^4 - 2*k^3 + 4*k^2 - 8*k + 16) / (k^5 + 32) / 2
(k - 1)*(k^3 + k^2 - 1)*(k^4 - 2*k^3 + 4*k^2 - 8*k + 16) / (3 + 32) / 2
# Rational fraction
(k - 1)*(k^3 + k^2 - 1)*(k^4 - 2*k^3 + 4*k^2 - 8*k + 16) / 70
# polynomial reduction
(k^8 − 2*k^7 + 3*k^6 − 7*k^5 + 15*k^4 + 2*k^3 − 4*k^2 − 24*k + 16) / 70
((k^3 - 2*k^2 + 3*k - 7) * (k^5 - 3) + (15*k^4 + 5*k^3 - 10*k^2 - 15*k - 5)) / 70
# k^5 - 3 = 0
(15*k^4 + 5*k^3 - 10*k^2 - 15*k - 5) / 70
(3*k^4 + k^3 - 2*k^2 - 3*k - 1) / 14

### Liniar System: manual Workout
# (k^4 + k^3 + k^2 + 3*k + 1) * (b0 + b1*k + b2*k^2 + b3*k^3 + b4*k^4) = N
# (  b0 +   b1 +   b2 + 3*b3 +   b4) = 0 # * k^4
# (  b0 +   b1 + 3*b2 +   b3 + K*b4) = 0 # * k^3
# (  b0 + 3*b1 +   b2 + K*b3 + K*b4) = 0 # * k^2
# (3*b0 +   b1 + K*b2 + K*b3 + K*b4) = 0 # * k
# (  b0 + K*b1 + K*b2 + K*b3 + 3*K*b4) = N


#############
### Example 2
K = 3
k = K^(1/5)
x = k^4 + 3*k^3 + 2*k^2 + k + 5

sol = inverse(c(5, 1, 2, 3, 1), 3, 5)
sol

# Test
sum(sol$sol * k^(0:(n-1))) * x / sol$N

(127 + 133*k + 447*k^2 - 55*k^3 - 273*k^4) / 3908


######################
######################

### OLD

### TODO:
# - cleanup;
# - use new version;

# install.packages("matlib")
library(matlib)

inverse.f <- function(b, K, N=1, showEq=FALSE) {
	# only for n == 5
	c.m <- matrix(
		c(b,
		tail(b, -1), K*head(b, 1),
		tail(b, -2), K*head(b, 2),
		tail(b, -3), K*head(b, 3),
		tail(b, -4), K*head(b, 4)), nrow=5)
	sol <- c(0,0,0,0, N)
	if(showEq) {
		showEqn(c.m, sol)
	}
	Solve(c.m, sol, fractions=TRUE)
	#
	k <- K^(1/5)
	x.s <- solve(c.m, sol)
	mult <- sum(x.s * k^(0:4)) * x
	print(mult)
	#
	return(x.s)
}

K = 2
N = 1
b <- c(1, 0, 0, -5, 4)
x.s <- inverse.f(b, K, N)
x.s

for(i in 1:15) inverse.f(c(1,-1,0,-5, i), K, N)

#
k <- K^(1/5)
x <- sum(b * k^(4:0)) # k^4 + k^3 + k^2 + 3*k + 1
# K = 3
x.s <- c(-1, -3, -2, 1, 3) # N=14, coefficients from Solve
# K = 2
x.s <- c(-9, -17, -16, 2, 36) # N=145, coefficients from Solve
#
sum(x.s * k^(0:4)) * x


###
K = 3
k = K^(1/5)
x = k^4 + 3*k^3 + 2*k^2 + k + 5

b <- c(1, 3, 2, 1, 5)
x.s <- c(127, 133, 447, -55, -273)
sum(x.s * k^(0:4)) * x # N = 3908

1/x = 1 / (k^4 + 3*k^3 + 2*k^2 + k + 5)
= (k - 3) / ((k - 3) * (k^4 + 3*k^3 + 2*k^2 + k + 5))
= (k - 3) / (k^5 − 7*k^3 − 5*k^2 + 2*k − 15)
= -(k - 3) / (7*k^3 + 5*k^2 - 2*k + 15 - K)
= -(k - 3)*(49*k^2 - 7*5*k + 39) / ((49*k^2 - 7*5*k + 39) * (7*k^3 + 5*k^2 - 2*k + 15 - K))
= -(k - 3)*(49*k^2 - 7*5*k + 39) / (343*k^5 + 853*k^2 − 498*k + 468)
= -(k - 3)*(49*k^2 - 7*5*k + 39) / (853*k^2 − 498*k + 468 + 343*K)
= -(k - 3)*(49*k^2 - 7*5*k + 39)*(k^3/853 + 498*k^2/853^2 - 1028937*k/853^3 - 1148327244/853^4) /
((k^3/853 + 498*k^2/853^2 - 1028937*k/853^3 - 1148327244/853^4) * (853*k^2 − 498*k + 1497))
= -(k - 3)*(49*k^2 - 7*5*k + 39)*(k^3/853 + 498*k^2/853^2 - 1028937*k/853^3 - 1148327244/853^4) /
(k^5 − 742024874205*k/853^4 − 1719045884268/853^4)
= (k - 3)*(49*k^2 - 7*5*k + 39)*(k^3/853 + 498*k^2/853^2 - 1028937*k/853^3 - 1148327244/853^4) /
(742024874205*k/853^4 + 1719045884268/853^4 - K)
# coefficients just exploded
= (k - 3)*(49*k^2 - 7*5*k + 39)*(k^3/853 + 498*k^2/853^2 - 1028937*k/853^3 - 1148327244/853^4) *
(853^4*k^4/742024874205 - 923308783101682270715*k^3/7341345519185947658427 + 8051333447413349530459657866125*k^2/363164065691292881683200598878369 - 70208332756971065129693346597147625896875*k/17965118010690530459652591536300161421371443 + 612222810135528879671954715221955543024087918828125/888704295464041380633895670848358331320955932710288521) /
((853^4*k^4/742024874205 - 923308783101682270715*k^3/7341345519185947658427 + 8051333447413349530459657866125*k^2/363164065691292881683200598878369 - 70208332756971065129693346597147625896875*k/17965118010690530459652591536300161421371443 + 612222810135528879671954715221955543024087918828125/888704295464041380633895670848358331320955932710288521) * (742024874205*k/853^4 + 1719045884268/853^4 - 3))

= (k - 3)*(49*k^2 - 7*5*k + 39)*(k^3/853 + 498*k^2/853^2 - 1028937*k/853^3 - 1148327244/853^4) *
(853^4*k^4/742024874205 - 923308783101682270715*k^3/7341345519185947658427 + 8051333447413349530459657866125*k^2/363164065691292881683200598878369 - 70208332756971065129693346597147625896875*k/17965118010690530459652591536300161421371443 + 612222810135528879671954715221955543024087918828125/888704295464041380633895670848358331320955932710288521) /
(k^5 + 50420161527450383717078056241969949442556083984375/296234765154680460211298556949452777106985310903429507)
# the rational fraction
= (k - 3)*(49*k^2 - 7*5*k + 39)*(k^3/853 + 498*k^2/853^2 - 1028937*k/853^3 - 1148327244/853^4) *
(853^4*k^4/742024874205 - 923308783101682270715*k^3/7341345519185947658427 + 8051333447413349530459657866125*k^2/363164065691292881683200598878369 - 70208332756971065129693346597147625896875*k/17965118010690530459652591536300161421371443 + 612222810135528879671954715221955543024087918828125/888704295464041380633895670848358331320955932710288521) /
(3 + 50420161527450383717078056241969949442556083984375/296234765154680460211298556949452777106985310903429507)
# 0.05695671 calculated in R with both formulas


