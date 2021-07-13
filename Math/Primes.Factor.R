########################
###
### Leonard Mada
### [the one and only]
###
### Prime Factorizations
###
### draft v.0.1a

# - some experiments with Prime factorizations;


####################

library(gmp)


### Helper Functions

mod.inv = function(x, N, type=1) {
	if(type == 1) {
		(x + 1/x)^N;
	} else {
		xn = x^N;
		xn + 1/xn;
	}
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

N = as.bigz(1002583) * as.bigz(3001073)
N = as.bigz("50003491") * as.bigz("84300971")

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


