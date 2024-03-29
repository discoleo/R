########################
###
### Leonard Mada
### [the one and only]
###
### Multi-Variable Polynomials
### Modular Arithmetic
### Derivation & Experiments
###
### draft v.0.2p-ex2




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
	} else if(n == 8) {
		mu = solve.ModP2(c(1,0,1), mod=mod);
		# (x^4 + 1);
		# TODO: more efficient way?
		mu4 = lapply(mu$Sol, function(r) solve.ModP2(c(mod - r,0,1), mod=mod)$Sol);
		mu$Sol  = c(mu$Sol, unlist(mu4));
		if(debug) print(mu);
		mu = mu$Sol;
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

# TODO

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
p = primes.mod(4, "Multiple")
countSol.mod(p)

# "Multiple" = 4 roots;
# "Mixed" = 2 roots;
# "Strict" = NO primes (except 2);

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
p = primes.mod(8, "Multiple")
countSol.mod(p)

###
validVals.mod(89, 8)
validValsS.mod(89, 8, e=2)

###
validVals.mod(97, 8)
validValsS.mod(97, 8, e=6)

###
validVals.mod(113, 8)
validValsS.mod(113, 8, e=7)

###
validVals.mod(137, 8)
validValsS.mod(137, 8, e=16)

###
validVals.mod(193, 8)
validValsS.mod(193, 8, e=16)

###
validVals.mod(233, 8)
validValsS.mod(233, 8, e=2)


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
pp = primes.mod(pow, "Multiple")
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
validValsS.mod(p, pow, e=3) # probably not the ideal;
validValsS.mod(p, pow, e=27)
validValsS.mod(p, pow, e=38)
# e = unity((p-1)/pow, mod=p)
# e = pow.mod(3^pow, (p-1)/pow - 1, mod=p)
# e = the 2nd formula works only with (3, 6)^5;

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
e = pow.mod(2^pow, (p-1)/pow - 1, mod=p)
print(e)

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


##########
### p = 71
p = 71;
validVals.mod(p, pow)
validValsS.mod(p, pow, e=23)
validValsS.mod(p, pow, e=34)
# only with base = 11;
e = pow.mod(11^pow, (p-1)/pow - 1, mod=p)
print(e)

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

##################
##################

### All Values
### Pow = 5

pow = 5
pp = primes.mod(pow, "Strict All")
countSol.mod(pp)


### p = 23
# x^5 = y %% p
# x^25 = y^5 %% p =>
# x^3 = y^5 %% p
# x^24 = y^40 %% p =>
# x^2 = y^40 %% p =>
# with correct sign:
# x^(3-2) = y^(5-40) %% p =>
# x = y^(-35) %% p
# x = y^9 %% p
# [old]
# x = y^20 %% p *OR* x = - y^20 %% p;
# Note: exponent 20 is only for p = 23!

p = 23
#
y = 9
x = pow.mod(y, 9, mod=p);
print(x)
# Test:
pow.mod(x, 5, mod=p); y;

###
y = 5
x = pow.mod(y, 9, mod=p);
print(x)
# Test:
pow.mod(x, 5, mod=p); y;

###
y = 6
x = pow.mod(y, 9, mod=p);
print(x)
# Test:
pow.mod(x, 5, mod=p); y;


##########
### p = 29
# x^5 = y
# x^30 = y^6 =>
# x^2 = y^6
# x^4 = y^12
# x = y^(-11) =>
# x = y^17

p = 29
powS = 17
#
y = 9
x = pow.mod(y, powS, mod=p);
print(x)
# Test:
pow.mod(x, 5, mod=p); y;

###
y = 10
x = pow.mod(y, powS, mod=p);
print(x)
# Test:
pow.mod(x, 5, mod=p); y;


###################
###################

### Pow = 8

pow = 8
pp = primes.mod(pow, "Multiple")
countSol.mod(pp)

###
p = 41
validVals.mod(p, pow)
validValsS.mod(p, pow, e=10)
# all bases seem to work;

### Roots of Unity
# (x+1)*(x-1)*(x^2 + 1)*(x^4 + 1)
mu = unity.mod(pow, p);
print(mu)
(mu^pow) %% p

### Other Roots
r = 2
r = (r * mu) %% p;
print(r)
pow.mod(r, pow, mod=p)

###
r = 3 # is a root of unity
r = (r * mu) %% p;
print(r)
pow.mod(r, pow, mod=p)

###
r = 5
r = (r * mu) %% p;
print(r)
pow.mod(r, pow, mod=p)


##############
##############

### Complex 1i

### p = 1 (mod p)
# by Fermat's sum of 2 squares =>
# p = a^2 + b^2
# a^2 = - b^2 (mod p)
# =>
# sqrt(-1) = c(a/b, -a/b);

### Note:
# sqrt(2i)  = c(1 + i, -1 - i);
# sqrt(-2i) = c(1 - i, -1 + i);

unityMinus = function(mod, detailed=FALSE) {
	# detailed: return base & (-1)^(1/4) where available;
	if(mod %% 4 == 3) return(NA);
	rn = (mod %% 8);
	if(rn == 5) {
		r = pow.mod(2, (mod-1)/4,  mod=mod);
		if(detailed) attr(r, "Mod") = list(Type=8, e=2);
		return(r);
	}
	### 9 Mod 16:
	rn = (mod %% 16);
	if(rn == 9) {
		r = pow.mod(2, (mod-1)/8,  mod=mod);
		# TODO: solve failures;
		if(r == 1 || r + 1 == mod) {
			r = pow.mod(3, (mod-1)/8,  mod=mod);
			# still a few failures!
			if(r == 1 || r + 1 == mod) { r = 3; }
			else {
				r2 = (r*r) %% mod; r4 = (r2*r2) %% mod;
				if(r4 + 1 == mod) {
					mneg = r; r = r2;
					if(detailed)
						attr(r, "Mod") = list(Type=16, Subtype = -1,
							e=3, hasUnits=TRUE, m=mneg, Pow=4);
				} else {
					if(detailed)
						attr(r, "Mod") = list(Type=16, Subtype = 0,
							e=3, hasUnits=FALSE);
				}
			}
		} else if(detailed)
			attr(r, "Mod") = list(Type=16, Subtype = 0,
				e=2, hasUnits=FALSE);
		return(r);
	}
	### 17 Mod 32:
	rn = (mod %% 32);
	if(rn == 17) {
		if(mod == 17) return(4);
		pow = (mod - 1) / 16;
		r = pow.mod(2, pow,  mod=mod);
		# TODO: there are some pow = 2 * pow;
		r2 = (r^2 %% mod);
		if(r + 1 == mod) {
			# pow/2, but pow is odd;
			r = pow.mod(3, pow, mod=mod);
			if(r == 1 || r + 1 == mod) {
				# - there are still a few extreme cases;
				r = 3;
			} else {
				r2 = (r*r) %% mod; r4 = (r2*r2) %% mod;
				# TODO: find these cases;
				if(r4 == 1) {
					r = r;
				} else if(r4 + 1 == mod) {
					r = r2;
				} else { r = r4; }
			}
		} else if((r2 %% pow) != 0) {
			r = (r*r) %% mod;
		}
		return(r);
	}
	# TODO: mod 64, 128, ...;
	stop("Not yet implemented!");
	# - fortunately, there are not many Fermat primes,
	#   even though these primes are easily solvable;
}


### MOD 16
pp = filter.mod(primes(15000), 9, mod=16)
print(pp)

i = sapply(pp, unityMinus);
r = sapply(seq(along=pp), function(id) (i[id]^2 + 1) %% pp[id]);
table(r)

tail(cbind(pp, i, r), n=10)

### SQRT(x)
p = pp[202]
i = unityMinus(p, detailed=T)
mn = attr(i, "Mod")$m; attr(i, "Mod") = NULL;
print(p)
c(mn^2 %% p, i, p-i);
### Case: sqrt(-1)
x = 30
r = pow.mod(x, (p+7)/16, mod=p)
r = (r*i) %% p; r = c(r, p - r);
r^2 %% p;
# Ex 2:
x = 2203
r = pow.mod(x, (p+7)/16, mod=p)
r = (r*i) %% p; r = c(r, p - r);
r^2 %% p;

### Case: (-1)^(1/4)
x = 25
r = pow.mod(x, (p+7)/16, mod=p)
r^2 %% p # NOT the root!
r = (r*mn) %% p; r = c(r, p - r);
r^2 %% p;
# Ex 2: (-1)^(3/4)
x = 125^2 %% p
r = pow.mod(x, (p+7)/16, mod=p)
r^2 %% p # NOT the root!
r = (r*mn*i) %% p; r = c(r, p - r);
c(x, r^2 %% p); # Squares
c(125, r); # Roots


### SQRT(x) (mod 281)
p = 281;
i = unityMinus(p)
### Ex 1:
x = 2
r = pow.mod(x, (p+7)/16, mod=p)
r^2 %% p # == - x;
r = (r * i) %% p; r = c(r, p-r);
c(x, r^2 %% p)
### Ex 2: (-1)^(1/4)
x = i;
# slightly different approach:
r = pow.mod(3, (p-1)/8, mod=p);
r = c(r, p-r); m4 = r[1];
c(x, r^2 %% p)
### Ex 3:
x = 5
r = pow.mod(x, (p+7)/16, mod=p)
r^2 %% p # NOT sqrt;
r = (r * m4 * i) %% p; r = c(r, p-r);
c(x, r^2 %% p)


### Failures:
# [old]
# still unsolved:
pp[c(12, 17)]
# already solved:
pp[c(2,3,5,6,13,19,20,21,23)]

# sometimes needs different base, e.g:
# 3, 12, or 20 or 28, or 17 (for 1097);

id = 4
p = pp[id]
i = unityMinus(p)
(i^2 %% p); pp[id];
sapply(seq(1, (p-1)/8), function(pow) pow.mod(12, 2*pow, mod=p))


##########
### MOD 32
pp = filter.mod(primes(1500), 17, mod=32)
print(pp)
# Simple failure:
pp[c(3,5,6,9,10,12,13,14)]
# Other failure:
pp[c(4,8)]

i = sapply(pp, unityMinus);
r = sapply(seq(along=pp), function(id) (i[id]^2 + 1) %% pp[id]);
table(r)

cbind(pp, i, r)

unlist(sapply(((pp[r == 2] - 1) / 16), function(x) c(factors(x), NA)))


###
id = 4
p = pp[id]
i = unityMinus(p)
print(i^2 %% p); print(p);


###
p = 17
i = unityMinus(p)
(i^2 %% p)

###
p = 113
i = unityMinus(p)
(i^2 %% p)

### => 2^(2*...)
p = 241
i = unityMinus(p)
(i^2 %% p)

### FAILS
p = 337
i = unityMinus(p)
(i^2 %% p)
# pow = 14 is sufficient;
sapply(seq(1, (p-1)/8), function(pow) pow.mod(6, 2*pow, mod=p))

### => 2^(2*...)
p = 401
i = unityMinus(p)
(i^2 %% p)

### => 2^(2*...)
p = 433
i = unityMinus(p)
(i^2 %% p)

###
p = 593
i = unityMinus(p)
(i^2 %% p)

### FAILS
p = 881
i = unityMinus(p)
(i^2 %% p)
# pow = 22 is sufficient;
sapply(seq(1, (p-1)/8), function(pow) pow.mod(22, 2*pow, mod=p))

### => 2^(2*...)
p = 977
i = unityMinus(p)
(i^2 %% p)

### => 2^(2*...)
p = 1009
i = unityMinus(p)
(i^2 %% p)


### [old]
p = 29
i = 12; (2 * inv.mod(5, p)) %% p;
i^2 %% p
x = c(1 + i, -1 - i); c(x^2, 2*i) %% p;
x = c(1 - i, -1 + i); c(x^2,-2*i) %% p;

y = 5
k = (p+3)/8
x = pow.mod(y, k, mod=p)
((x)^2 %% p); y;
((x * 2^(2*k-1))^2 %% p); y;
((x * i)^2 %% p); y;