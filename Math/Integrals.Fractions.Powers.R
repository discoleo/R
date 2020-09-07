
########################
###
### Leonard Mada
### [the one and only]
###
### Exact Integration of Polynomial Fractions
### - Higher Powers
###   Integral( 1 / Q(x)^p )dx
### - Polynomial fractions:
###   Integral( P(x) / Q(x)^p )dx
###
### draft v.0.1a


### History

# v.0.1a:
# - basic fractions of Order 2: P(x) / (x^2 + b*x + c)^2;
# - TODO: x^3 / (x^2 + b*x + c)^2; [easy]


### Related Material:

# for Roots of unity: Baseline
# - see file: Integrals.Fractions.Unity.R;
# for Roots of unity: Powers
# - see file: Integrals.Fractions.Unity.Powers.R;


###############
###############

### Terminology

### Q(x) = x^2 + b*x + c
# F(k, p) = x^k / (x^2 + b*x + c)^p; # only the Fraction!
# I(k, p) = Integral( x^k / (x^2 + b*x + c)^p ) dx;


################
### Solution ###


###########
### Step 1:
# - creating a System of equations;

### d(1 / Q(x))
# -(2*x + b) / Q(x)^2;

### F(0, 1)
# = Q(x) / Q(x)^2
# (x^2 + b*x + c) / Q(x)^2;

### d(x / Q(x))
# -(x^2 - c) / Q(x)^2;

### F(1, 1)
# = x * Q(x) / Q(x)^2
# (x^3 + b*x^2 + c*x) / Q(x)^2;


###########
### Step 2:
# - solving system;

### F(0, 2)
# (2*c - b^2/2) * F(0, 2) = F(0,1) + d(x / Q(x)) + b/2*d(1 / Q(x))
### I(0, 2)
# = (I(0,1) + x / Q(x) + b/2 * 1/Q(x)) / (2*c - b^2/2)

### F(1, 2)
# (b - 4*c/b) * F(0, 2) = F(0,1) + d(x / Q(x)) + 2*c/b * d(1 / Q(x))
### I(1, 2)
# = (I(0,1) + x / Q(x) + 2*c/b * 1/Q(x)) / (b - 4*c/b)


### F(2, 2)
# - d(x / Q(x)) + c * F(0, 2)
### I(2, 2)
# - x / Q(x) + c * I(0, 2) # I(0, 2) is computed above;


### F(3, 2)
# TODO



############

F.f = function(x, b, k=0, p=1) {
	# b[2]*x + b[1]
	if(k == 0) {
		r = 1 / (x^2 + b[1]*x + b[2])^p
	} else {
		r = x^k / (x^2 + b[1]*x + b[2])^p
	}
	r
}
F.range = function(lim, b, k=0, p=1) {
	F.f(lim[2], b,k,p) - F.f(lim[1], b,k,p)
}
I.f = function(lim, b, k=0, p=1) {
	integrate(F.f, lower=lim[1], upper=lim[2], b=b, k=k, p=p)
}

##############

### Test
b = c(2, 3)
lim = c(2, 4)
### I(0, 2)
I.f(lim, b, 0, 2)
(I.f(lim, b, 0, 1)$value + F.range(lim, b, 1, 1) + b[1]/2 * F.range(lim, b, 0, 1)) / (2*b[2] - b[1]^2/2)

### I(1, 2)
I.f(lim, b, 1, 2)
(I.f(lim, b, 0, 1)$value + F.range(lim, b, 1, 1) + 2*b[2]/b[1] * F.range(lim, b, 0, 1)) / (b[1] - 4*b[2]/b[1])

### I(2, 2)
I.f(lim, b, 2, 2)
b[2] * I.f(lim, b, 0, 2)$value - F.range(lim, b, 1, 1)

### I(3, 2)
I.f(lim, b, 3, 2)
# TODO



