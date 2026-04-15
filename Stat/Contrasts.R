
### Contrasts

# - Polynomial contrasts are generated in R
#   using function contr.poly;
# - Hardcoded Limit: R stops if n > 95;


####################

### Helper Functions

# Note:
# - Functions have been moved to the referenced source file;

source("Contrasts.Helper.R")

####################

### Helper Functions

# Function make.poly:
# - Inner function in contr.poly;


# Note:
# if(n %% 2 == 1) poly(n)[seq(2, n-1, by = 2), (n+1)/2]
# should be probably == 0;
# FAILS for n = 23;

sapply(seq(15, 41, by=2), \(n) {
	# n = 21: still (approx) 0;
	# BUT for n = 23: jumps to -0.122;
	n2 = (n+1)/2;
	make.poly(n)[n2, c(n-3, n-1)];
})

# Separation: Probably 0 vs Non-Zero
rg = sapply(seq(3, 40), \(id, tol = 1E-9) {
	x = make.poly(id, seq(id));
	x = abs(x[x != 0]);
	isE = x <= tol;
	if(sum(isE) > 0) {
		err = x[isE];
		x = x[! isE]
	} else return(c(0, min(x)));
	c(max(err), min(x));
})
cbind(c(NA,NA), c(NA,NA), rg)

# Note:
# - lowest abs values for n = 23 is probably
#   around 1.36E-6 / 2 = 7E-7 (and therefore larger than 4.7E-7);
# - for n = 24: 3.5E-7;


### Upper Limit

# BUT: Failure of P[22] already for n = 23;
sapply(seq(100), \(id, tol = 1E-10) {
	x = make.poly(id, seq(id));
	x = abs(x[x != 0]);
	# 1E-8 may still occur for n > 26;
	sum(x <= tol);
})


###############

### Components
# Note:
# - Signs NOT set;

# Component 1: 1 / sqrt(n);
make.poly(11)[1:5, 1]; 1 / sqrt(11);

#####################

# Component 2: Linear

# Odd & Even (but misses 1 coefficient)
n = 35;
n2 = (n-1)/2;
make.poly(n)[seq(1, n2), 2]; seq(n2, 1) / sqrt((n^3 - n) / 12);

# For even:
n = 36;
n2 = (n-1)/2; 
make.poly(n)[seq(1, n2+1), 2]; c(seq(n2, 1/2, by=-1)) / sqrt((n^3 - n) / 12);

# TODO: with sign;

########################

# Component 3: Quadratic
# Note: different Structure;


# Only with Sign:
n = 17; # n = 35; # n = 20;
make.poly(n)[1:((n+1)/2), 3]; poly.cp3(n);

# Note:
# - Algorithm based on Matrix factorization
#   is remarkably stable up to n = 162; [for Component 3]
n = 161
z  = make.poly(n)[1:((n+1)/2), 3]
ze = poly.cp3.mpfr(n);
summary(abs(z - ze))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 2.900e-21 7.708e-18 1.531e-17 2.270e-17 2.361e-17 4.534e-16 

n = 162
z  = make.poly(n)[1:((n+1)/2), 3]
n1 = mpfr(n, 240); ze = poly.cp3.mpfr(n1);
summary(abs(z - ze))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 5.888e-19 9.664e-18 2.523e-17 2.487e-17 3.503e-17 1.156e-16
min(abs(ze)) # is of Order 9.9E-4 (not very small);

# Matrix factorization fails for n >= 163;


#######################

### Component 4: Cubic

# Examples:
make.poly( 5)[1:2, 4]; poly.cp4( 5);
make.poly( 6)[1:3, 4]; poly.cp4( 6);
make.poly( 7)[1:3, 4]; poly.cp4( 7);
make.poly( 8)[1:4, 4]; poly.cp4( 8);
make.poly( 9)[1:4, 4]; poly.cp4( 9);


### Numerical Stability:
library(Rmpfr)

###
n = 161;
(make.poly(n)[seq(1, (n+1)/2), 4] - poly.cp4.mpfr(n)) |> abs() |> summary();

### n = 4...161
sapply(seq(4, 161), \(n, qq = c(6,7,8) / 8) {
	d = make.poly(n)[seq(1, (n+1)/2), 4] - poly.cp4.mpfr(n);
	as.numeric(quantile(abs(d), qq));
}) |> apply(1, quantile, c(0.5, 3/4, 1));

