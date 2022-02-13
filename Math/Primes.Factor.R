########################
###
### Leonard Mada
### [the one and only]
###
### Prime Factorizations
###
### draft v.0.1f

# - some experiments with Prime factorizations;


####################

library(gmp)


### Helper Functions

simplify.mod = function(x, f=c(2,3)) {
	x = x - 1;
	for(p in f) {
		while(x %% p == 0) x = x %/% p;
	}
	return(x);
}
seq.mod = function(N, p) {
	x = seq(N-1);
	for(px in p) {
		x[seq(px, N-1, by=px)] = NA;
	}
	x = as.bigz(x, mod=N);
	return(x);
}
mod.inv = function(x, N, type=1) {
	if(type == 1) {
		(x + 1/x)^N;
	} else {
		xn = x^N;
		xn + 1/xn;
	}
}
bezout.mod = function(x, N) {
	xn = x^N;
	xi = 1/x; xin = 1/xn;
	(x + xi)^N - (xn + xin);
}
pollard = function(x0, N, iter=14, pow=4) {
	x = x0;
	# xn = mod.inv(x0, N);
	for(n.cycle in seq(iter)) {
		y = x0;
		for(i in seq(2^n.cycle)) {
			x = x^pow; # x^4;
			# x = x*x - xn;
			# x = x^2 - 5*x + 5;
			f = gcd(x - y, N);
			if(f > 1) return(list(f=f, i=i, cycle=n.cycle));
		}
	}
	print("NO factors found!")
}
factorizePow.mod = function(base, mod, iter=1000) {
	x = base;
	g = gcd(x, mod);
	if(g > 1) return(g);
	for(i in seq(iter)) {
		x2 = (x*x);
		x2 = (x2*x2);
		if(x2 == base) {print(i); return(x);}
		x = x2;
		g = gcd(x-1, mod);
		if(g > 1) {print(paste0("Simple: i = ", i)); return(g);}
		g = gcd(x-base, mod);
		if(g > 1) {print(paste0("Enhanced: i = ", i)); return(g);}
	}
	return(NA)
}

# N = as.bigz(1002583) * as.bigz(3001073)
N = as.bigz("50003491") * as.bigz("84300971")

# the enhanced version functions better:
# - it requires an even number of primes: possibly pairs of similar size;
# the simple version did NOT really function on this example;
base = (((((as.bigz(2*3*5*7*11*13, mod=N) * 17*19*23*29)*31*37*41*43)
	*47*53*59*61)*67*71*73*79)*83*89*97)*101*103*107
	# again slower: *109 # and recovers again: *113
factorizePow.mod(base, N, iter=60000)

# Pollard is reasonably efficient with some specific starting values:
x = as.bigz(1024*4, mod=N);
# x = 5 * as.bigz(2, mod=N)^8 - 1
# x = as.bigz(2, mod=N); x = x + 1/x;
x; mod.inv(x, N);

pollard(x, N, pow=4)
pollard(x, N, pow=8)
pollard(x, N, pow=16)

### x + Inverse(x)
x = as.bigz(2, mod=N); x = x + 1/x;
x
pollard(x, N, pow=5)
pollard(x, N, pow=25)
pollard(x, N, pow=125) # BEST!


### Bezout
# [rather inefficient]
x = as.bigz(2, mod=N);
x = bezout.mod(x, N);
x
pollard(x, N, pow=3)
pollard(x, N, pow=27)

### Bezout
x = as.bigz(3, mod=N);
x = bezout.mod(x, N);
x
pollard(x, N, pow=2)
pollard(x, N, pow=16)
pollard(x, N, pow=64)
pollard(x, N, pow=128)


######################
######################

################
### Analysis ###
################


### Types of numbers (mod p1*p2)
# - depending on the number of solutions
#   of the equation: x + 1/x = b;

### Type 0:
# - NO solutions for given b:
#   x^2 - b*x + 1 = 0 (mod N);
# - Size: ~ 75% of numbers (parameter b);

### Type 1: Squares
# - very special: x^2 = 1 (mod N);
# - Size: 4;
#   x0 = (1, -1, (x0 == 1/x0), -(x0 == 1/x0));

### Type 2:
# - only 2 solutions of eq:
#   x^2 - b*x + 1 = 0 (mod N);
# - easy factorization of N:
#   gcd(x - 1/x, N) > 1;
# - Size: relatively rare;
#   initial approx = 2*(p1 + p2) - 12;
#   (but NOT correct; TODO: exact)

### Type 4:
# - exactly 4 solutions of eq:
#   x^2 - b*x + 1 = 0 (mod N);
# - TODO: factorization of N;


### Example 1:
p1 = 97; p2 = 47;
N = p1 * p2;

x = seq.mod(N, c(p1, p2));
x = x + (1/x);
# Types of numbers:
tbl = table(as.integer(x))
head(tbl, 20)
table(tbl)

### Squares
tbl[tbl == 1]
which(x == 1647)
# very special, but there are only 2 which are useful:
# x0 = (1, x0, - x0, - 1), where x0 == 1/x0;
# special properties: x0^2 = 1;

### Quadratic-like
x.all = as.integer(names(tbl[tbl == 2]))
x.all
as.vector(sapply(x.all, function(v) which(x == v)))

head(tbl[tbl == 2], 20)
gcd(diff(which(x == 45)), N)


### Gaps
xs = sort(as.integer(x))
head(xs, 100)
# gap is actually diff() - 1;
table(diff(xs))


##############
### Example 2:
p1 = 1229; p2 = 1951;
N = p1 * p2;

x = seq(N-1);
x[x %% p1 == 0] = NA
x[x %% p2 == 0] = NA
x = as.bigz(x, mod=N)
x = x + (1/x);
# Types of numbers:
tbl = table(as.integer(x))
head(tbl, 20)
table(tbl)
# Type 2: 5.3 / 1000;

gcd(diff(which(x == 1231)), N)


### Example 3:
p1 = 97; p2 = 59;
N = p1 * p2;

x = seq(N-1)
x[x %% p1 == 0] = NA; x[x %% p2 == 0] = NA;
x = as.bigz(x, mod=N);
x = x + (1/x);
# Types of numbers:
tbl = table(as.integer(x))
table(tbl)
head(tbl[tbl == 2], 20)

which(x == 61)
# 1299 4485
gcd(diff(which(x == 61)), N)


################

### experimental
N = as.bigz("50003491") * as.bigz("84300971")

x = as.bigz(seq(2, 200), mod=N)
x = x + 1/x;

### not yet useful:
fact.experimental = function(x, N) {
	for(id in seq(length(x))) {
		xx = x[id];
		g = gcd(xx, N);
		if(g > 1) {
			return(xx);
		}
	}
	for(id in seq(length(x))) {
		xx = x[id] + 1;
		g = gcd(xx, N);
		if(g > 1) {
			return(xx);
		}
	}
	for(id in seq(length(x))) {
		xx = x[id] - 1;
		g = gcd(xx, N);
		if(g > 1) {
			return(xx);
		}
	}
	print("NO factor!")
	return(0);
}

fact.experimental(x, N=N);



######################
######################

### Polynomial Decompositions

N = 35

### Decomposition into 5*7

a = as.bigz(2, N)
# Inverse
ai = 1/a # = (N+1) / 2;
ai

### Known
a^N - ai^N

### Unknown
(a^7 - ai^7) * (a^(4*7) + ai^(4*7) + a^(2*7) + ai^(2*7) + 1)

### Known
a^N + ai^N

### Unknown
(a^7 + ai^7) * (a^(4*7) + ai^(4*7) - a^(2*7) - ai^(2*7) + 1)


### How do we compute N ???
### System:
x = (a^7 + ai^7); y = (a^7 - ai^7);
x * ((x^2-2)^2 - x^2 + 1) # = a^N + ai^N
x * (x^4 - 5*x^2 + 5) # = a^N + ai^N

y * ((x^2-2)^2 + x^2 - 3) # = a^N - ai^N
y * ((y^2+2)^2 + y^2 + 1) # = a^N - ai^N
y * (y^4 + 5*y^2 + 5) # = a^N - ai^N

### Discrete Gradient
# - (efficient) method to compute discrete Gradient or df?

### TODO:
factorize.fast = function(x, N, len=1000, iter=100000, print=FALSE) {
	p2 = seq(len);
	tbl = p2^2
	tbl = sort(tbl)
	if(print) print(head(tbl))
	xi = x; tbl.max = max(tbl);
	for(i in seq(iter)) {
		xi = xi * xi;
		if(xi <= tbl.max) {
			if(any(xi %in% tbl)) return(list(n=i, xi=xi));
		}
	}
	print("Not yet factorized!")
	return(0)
}


N = as.bigz(1002583) * as.bigz(3001073)

x = 5 * as.bigz(2, mod=N)^8 - 1
x

factorize.fast(x, N)


###
factorize.fast = function(x1, x2, N, iter=100000, print=FALSE) {
	for(i in seq(iter)) {
		x1 = x1 * x1;
		x2 = x2 * x2;
		if(x1 - x2 == 0) {
			return(list(n=i, x1=x1, x2=x2));
		}
	}
	print("Not yet factorized!")
	return(0)
}

x1 = 3 * as.bigz(83, mod=N) + 5
x2 = x1 - 7
x1; x2;

factorize.fast(x1, x2, N)


