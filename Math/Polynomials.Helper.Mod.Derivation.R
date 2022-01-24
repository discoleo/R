########################
###
### Leonard Mada
### [the one and only]
###
### Multi-Variable Polynomials
### Modular Arithmetic
### Derivation & Experiments
###
### draft v.0.1h-lowCount




######################

### Helper Functions

source("Polynomials.Helper.R")
# - is automatically loaded in: Polynomials.Helper.R;
#   source("Polynomials.Helper.Factorize.R")


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
p = p[p %% (2*pow) != 1]
sapply(p, function(p) p - length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
#
p = primes(500)
p = p[p %% (2*pow) == 1]
sapply(p, function(p) p - length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
#
sort(unique( sapply(seq(498), pow.mod, pow, mod=499) ))


###############
### Power 9 ###
###############

pow = 9;
# simple: (p %% 3 != 1) seems sufficient;
p = primes(500)
p = p[p %% (3) != 1]
sapply(p, function(p) p - length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
#
p = primes(500)
p = p[p %% (3) == 1]
sapply(p, function(p) p - length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
# very low number of solutions:
p = primes(500)
p = p[p %% 9 == 1]
print(rbind(p,
	sapply(p, function(p) length(unique( sapply(seq(p-1), pow.mod, pow, mod=p) )))
))
#
sort(unique( sapply(seq(460), pow.mod, pow, mod=461) ))

