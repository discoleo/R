########################
###
### Leonard Mada
### [the one and only]
###
### Multi-Variable Polynomials
### Modular Arithmetic
### Derivation & Experiments
###
### draft v.0.2h




######################

### Helper Functions

source("Polynomials.Helper.R")
# - is automatically loaded in: Polynomials.Helper.R;
#   source("Polynomials.Helper.Factorize.R")


### Valid Solutions:

validVals.mod = function(p, pow) {
	table(sapply(seq(p-1), pow.mod, n=pow, mod=p))
}
validValsP2.mod = function(p, pow, e=2) {
	sapply(seq(0, (p-1)/pow - 1), function(pow) pow.mod(e, pow, mod=p));
}
validValsS.mod = function(p, pow, e=2, len=80 + 3*16) {
	v = sort(validValsP2.mod(p, pow=pow, e=e));
	print.asMatrix(v, len=len);
	return(invisible(v));
}
print.asMatrix = function(v, len=80) {
	v = format(v);
	nch = nchar(v[[1]]);
	nc  = len %/% (nch + 1);
	nr  = length(v) %/% nc;
	nM  = nr*nc;
	id  = seq(1, nM);
	m = matrix(v[id], nrow=nc, ncol=nr); # T(m)
	m = apply(m, 2, paste0, collapse=" ");
	cat(m, sep="\n");
	if(nM < length(v)) {
		m = v[seq(nM+1, length(v))];
		m = paste0(m, collapse=" ");
		cat(m, sep="\n");
	}
	invisible();
}

### Roots of Unity
unity.mod = function(n, mod, debug=FALSE) {
	if(n == 2) {
		# Nothing;
	} else if(n == 3) {
		mu = solve.ModP2(c(1,1,1), mod=mod);
		if(debug) print(mu);
		mu = mu$Sol;
	} else if(n == 4) {
		mu = solve.ModP2(c(1,0,1), mod=mod);
		if(debug) print(mu);
		mu = mu$Sol;
	} else if(n == 5) {
		mu = solve.ModP2(c(-1,1,1), mod=mod);
		if(debug) print(mu);
		mu = unlist(lapply(mu$Sol,
			function(S) solve.ModP2(c(1,-S,1), mod=mod)$Sol));
	} else {
		stop("Not yet implemented!");
	}
	#
	if(n %% 2 == 0) {
		mu = c(1, mod - 1, mu);
	} else mu = c(1, mu);
	mu = sort(mu);
	return(mu);
}


#######################
#######################


### Inv

###
inv.mod(5, 13)
r = 13 %% 5
13 - (13 * inv.mod(r, 5) - 1) / 5

###
inv.mod(5, 17)
r = 17 %% 5
17 - (17 * inv.mod(r, 5) - 1) / 5

###
inv.mod(5, 19)
r = 19 %% 5
19 - (19 * inv.mod(r, 5) - 1) / 5

###
x = 7
p = 19; # gcd(x, p) == 1;
inv.mod(x, p)
r = p %% x
p - (p * inv.mod(r, x) - 1) / x


###################
###################

# sc = some arbitrary "scale";
sqrt.mod.Experimental = function(p, sc = 2) {
	sc2 = sc^2;
	r =  p %% sc2;
	scl = if(r == 1 || r == p-1) 1 else inv.mod(r, sc2);
	r  = if(r == 1) 1 else if(r == p-1) -1 else 1;
	x2 = (p*scl - r) / sc2;
	if(r > 0) x2 = p - x2;
	(x2 - inv.mod(sc2, p))
	# x = sqrt(x2) (mod p)
	x = inv.mod(sc, p)
	x = c(x, p - x)
	return(list(x=x, xsq=x2, Mod=p));
}
test.sqrt = function(x, xsq, p) {
	if(is.list(x)) {
		xsq = x$xsq; p = x$Mod; x = x$x;
	}
	err = (x*x - xsq) %% p;
	cat(paste0("Error: ", paste0(err, collapse=", "), "\n"));
	print(xsq); print(x);
	invisible(list(x=x, xsq=xsq, Mod=p));
}


### SQRT

###
p = 101;
sc = 2;
r = test.sqrt(sqrt.mod.Experimental(p=p, sc))
(r$x * inv.mod(sc, p)) %% p # sqrt();
(r$xsq * inv.mod(sc*sc, p)) %% p # x^2


### Type: Classic Squares
p = 101
sort(unique( (seq(p-1)^2) %% p) )
x0 = seq(2, floor(sqrt(p)));
x0 = x0*x0; xinv = sapply(x0, inv.mod, mod=p);
x = c(1, x0, xinv);
sort(unique(x))
# "Other Squares"
len = length(x0)
for(id in seq(len)) {
	id2 = seq(len)[-id];
	# NOT unique due to reducible fractions:
	# like (36/81);
	x = c(x, (x0[id] * xinv[id2]) %% p );
}
sort(unique(x))


### 2^2 => sqrt(10) = 1/2 (mod 13)

# 1 (mod 13)
(4 * -3) %% 13
(10 - inv.mod(4, 13))
# x = sqrt(10) (mod 13)
x = inv.mod(2, 13)
x = c(x, 13 - x)
(x*x - 10) %% 13
print(x);

# * 3 => sqrt(12)
(10 * 9) %% 13 # 12 (mod 13)
x2 = (3 * x) %% 13
(x2*x2 - 12) %% 13
print(x2);

# * 5 => sqrt(3)
(10 * 25) %% 13 # 3 (mod 13)
x2 = (5 * x) %% 13
(x2*x2 - 3) %% 13
print(x2);

# * 7 => sqrt(9) # trivial
(10 * 49) %% 13 # 3 (mod 13)
x2 = (7 * x) %% 13
(x2*x2 - 9) %% 13
print(x2);


##########

### 2^2 => sqrt(8) = 1/2 (mod 31)
# 2^4 also possible (mod 31);

p = 31
# 1 (mod p)
(4 * 8) %% p
(8 - inv.mod(4, p))
# x = sqrt(8) (mod p)
x = inv.mod(2, p)
x = c(x, p - x)
(x*x - 8) %% p
print(x);


###
p = 67
#
r =  p %% 4;
r = if(r == 1) 1 else -1;
x2 = (p - r) / 4;
if(r > 0) x2 = p - x2;
(x2 - inv.mod(4, p))
# x = sqrt(x2) (mod p)
x = inv.mod(2, p)
x = c(x, p - x)
(x*x - x2) %% p
print(x2); print(x);


###
p = 37
#
scm = 3;
sc2 = scm^2;
r =  p %% sc2;
sc = if(r == 1 || r == p-1) 1 else inv.mod(r, sc2);
r = if(r == 1) 1 else if(r == p-1) -1 else 1;
x2 = (p*sc - r) / sc2;
if(r > 0) x2 = p - x2;
(x2 - inv.mod(sc2, p))
# x = sqrt(x2) (mod p)
x = inv.mod(scm, p)
x = c(x, p - x)
(x*x - x2) %% p
print(x2); print(x);


####################
####################

#############
### Cubes ###
#############


### Types of Primes:
sapply(primes(101), function(p) length(unique( (seq(p-1)^3) %% p )))
# Primes: 5 (mod 6)
sapply(c(53, 59, 71, 83, 89, 101), function(p) length(unique( (seq(p-1)^3) %% p )))

### All:
p = primes(1000)
p = p[p %% 6 == 5]
sapply(p, function(p) p - length(unique( (seq(p-1)^3) %% p )))
#
p = primes(1000)
p = p[p %% 6 != 5]
sapply(p, function(p) p - length(unique( (seq(p-1)^3) %% p )))


### Type: Classic Cubes
# - does NOT work well;
# top = 1/3: even worse;
p = 101
top = 1/2;
sort(unique( (seq(p-1)^3) %% p ))
x0 = seq(2, floor(p^top));
x0 = (x0*x0*x0) %% p; xinv = sapply(x0, inv.mod, mod=p);
x = c(1, x0, xinv);
sort(unique(x))
# "Other Cubes"
len = length(x0)
for(id in seq(len)) {
	id2 = seq(len)[-id];
	# NOT unique due to reducible fractions;
	x = c(x, (x0[id] * xinv[id2]) %% p );
}
sort(unique(x))


###############
###############

###############
### Power 5 ###
###############

pow = 5;
# pow = 7; # pow = 11; # pow = 13;
p = primes(500)
isSol = (p %% (2*pow) != 1);
sapply(p[isSol], function(p) p - length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
#
isSol = (p %% (2*pow) == 1);
sapply(p[isSol], function(p) p - length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
# low number of solutions:
isSol = (p %% (2*pow) == 1);
print(rbind(p[isSol],
	sapply(p[isSol], function(p) length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
))
# Ex All:
sort(unique( sapply(seq(498), pow.mod, pow, mod=499) ))
# Ex few:
table( sapply(seq(30), pow.mod, pow, mod=31) )


###############
### Power 9 ###
###############

pow = 9;
# simple: (p %% 3 != 1) seems sufficient;
p = primes(500)
isSol = (p %% 3 != 1)
sapply(p[isSol], function(p) p - length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
#
isSol = (p %% 3 == 1)
sapply(p[isSol], function(p) p - length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
# very low number of solutions:
isSol = (p %% pow == 1)
print(rbind(p[isSol],
	sapply(p[isSol], function(p) length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
))
# Ex All:
sort(unique( sapply(seq(460), pow.mod, pow, mod=461) ))


################
### Power 15 ###
################

pow = 15;
# simple: (p %% c(3, 5) != 1) seems sufficient;
p = primes(500)
isSol = (p %% 3 != 1) & (p %% 5 != 1)
sapply(p[isSol], function(p) p - length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
#
isSol = (p %% 3 == 1) | (p %% 5 == 1)
sapply(p[isSol], function(p) p - length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
# very low number of solutions:
isSol = (p %% pow == 1)
print(rbind(p[isSol],
	sapply(p[isSol], function(p) length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
))
# Ex All:
sort(unique( sapply(seq(442), pow.mod, pow, mod=443) ))


###############

### Even Powers

###############
### Power 6 ###
###############

pow = 6;
# pow = 10; # pow = 14;
p = primes(500)
# - NOT all solutions, but still 2 Categories:
# - many: ~ 1/2 of values;
# - few: each (x^pow) repeats (pow) times;
isSol = (p %% (pow) != 1);
sapply(p[isSol], function(p) p - length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
#
isSol = (p %% (pow) == 1);
sapply(p[isSol], function(p) p - length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
# low number of solutions:
isSol = (p %% (pow) == 1);
print(rbind(p[isSol],
	sapply(p[isSol], function(p) length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
))
# Ex many:
sort(unique( sapply(seq(346), pow.mod, pow, mod=347) ))
# Ex few:
table( sapply(seq(60), pow.mod, pow, mod=61) )


########################
########################

### Valid Values ###

###########
### Pow = 3
p = primes.mod(3, "Multiple")
# - only "Multiple" have multiple solutions;
countSol.mod(p)

###
validVals.mod(43, 3)
validValsS.mod(43, 3, e=2)

###
validVals.mod(61, 3)
validValsS.mod(61, 3, e=8)

###
validVals.mod(67, 3)
validValsS.mod(67, 3, e=5)

###
validVals.mod(73, 3)
validValsS.mod(73, 3, e=7)

###
validVals.mod(79, 3)
validValsS.mod(79, 3, e=12)

###
validVals.mod(97, 3)
validValsS.mod(97, 3, e=19)


###########
### Pow = 4
p = primes.mod(4, "Most")
countSol.mod(p)

# TODO: "Multiple"
# - Note:
#   "Most" = 4 roots;
#   "Multiple" = 2 roots;

###
validVals.mod(37, 4)
validValsS.mod(37, 4, e=7)

###
validVals.mod(41, 4)
validValsS.mod(41, 4, e=4)

###
validVals.mod(53, 4)
validValsS.mod(53, 4, e=13)

###
validVals.mod(61, 4)
validValsS.mod(61, 4, e=12)

###
validVals.mod(73, 4)
validValsS.mod(73, 4, e=18)


###########
### Pow = 5
p = primes.mod(5, "Multiple")
countSol.mod(p)

###
validVals.mod(41, 5)
validValsS.mod(41, 5, e=3)

###
validVals.mod(61, 5)
validValsS.mod(61, 5, e=21)

###
validVals.mod(71, 5)
validValsS.mod(71, 5, e=23)

###
validVals.mod(101, 5)
validValsS.mod(101, 5, e=32)

###
validVals.mod(131, 5)
validValsS.mod(131, 5, e=18)


###########
### Pow = 8
p = primes.mod(8, "Most")
countSol.mod(p)

###
validVals.mod(89, 8)
sort(validValsP2.mod(89, 8, e=2))

###
validVals.mod(97, 8)
sort(validValsP2.mod(97, 8, e=6))

###
validVals.mod(113, 8)
sort(validValsP2.mod(113, 8, e=7))

###
validVals.mod(137, 8)
sort(validValsP2.mod(137, 8, e=16))

###
validVals.mod(193, 8)
sort(validValsP2.mod(193, 8, e=16))

###
validVals.mod(233, 8)
sort(validValsP2.mod(233, 8, e=2))


#####################
#####################

#####################
### Roots Order n ###
#####################

###########
### Pow = 3

pow = 3
pp = primes.mod(pow, "Multiple")
countSol.mod(pp)

###
p = 43;
validVals.mod(p, pow)
validValsS.mod(p, pow, e=2)

### Roots of Unity
mu = unity.mod(pow, p);
print(mu)
(mu^pow) %% p

### Other Roots
r = 2
r = (r * mu) %% p;
print(r)
(r^pow) %% p

###
r = 3
r = (r * mu) %% p;
print(r)
(r^pow) %% p

###
r = 4
r = (r * mu) %% p;
print(r)
(r^pow) %% p


##################

### Pow = 4

pow = 4
pp = primes.mod(pow, "Most")
countSol.mod(pp)

###
p = 41
validVals.mod(p, pow)
validValsS.mod(p, pow, e=4)

### Roots of Unity
# (x+1)*(x-1)*(x^2 + 1)
mu = unity.mod(pow, p);
print(mu)
(mu^pow) %% p

### Other Roots
r = 2
r = (r * mu) %% p;
print(r)
(r^pow) %% p

###
r = 3
r = (r * mu) %% p;
print(r)
(r^pow) %% p


##################

### Pow = 5

pow = 5
pp = primes.mod(pow, "Multiple")
countSol.mod(pp)

### p = 41
p = 41;
validVals.mod(p, pow)
validValsS.mod(p, pow, e=3)

### Roots of Unity
# (x-1)*(x^4 + x^3 + x^2 + x + 1)
# - P[4] is a strictly symmetric polynomial:
#   (S^2 + S - 1)
mu = unity.mod(pow, p);
print(mu)
(mu^pow) %% p

### Other Roots
r = 2
r = (r * mu) %% p;
print(r)
(r^pow) %% p

###
r = 3
r = (r * mu) %% p;
print(r)
(r^pow) %% p


##########
### p = 61
p = 61;
validVals.mod(p, pow)
validValsS.mod(p, pow, e=21)

### Roots of Unity
mu = unity.mod(pow, p);
print(mu)
(mu^pow) %% p

### Other Roots
r = 2
r = (r * mu) %% p;
print(r)
(r^pow) %% p

###
r = 3
r = (r * mu) %% p;
print(r)
(r^pow) %% p

