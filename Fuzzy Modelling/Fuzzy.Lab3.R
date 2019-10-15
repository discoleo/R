#### Fuzzy Sets ####

### LAB 3: Mean & Median
### Ch 2.2.2 / page 26 [19]
###
### Leonard Mada

### draft 0.1


### Function definitions

tri.f.base = function(x, a) {
	### triangular function
	# non-vectorized function
	if(x < a[1]) {
		return(0)
	} else if(x > a[3]) {
		return(0)
	} else if(x < a[2]) {
		f = (x - a[1]) / (a[2] - a[1])
	} else {
		# a[2] ... a[3]
		f = (a[3] - x) / (a[3] - a[2])
	}
}
# we need vectorized functions for:
# curve(), integrate()
tri.f = function(x, a) {
	# vectorized function
	return(sapply(x, tri.f.base, a))
}

# compute area
area.tri.f = function(x, a) {
	return(x * tri.f(x, a))
}

### Test

# triangular cutoffs
a = c(175, 185, 200)

### PLOT
# Membership Function
curve(tri.f(x, a), from=50, to=250, col="blue")

### Mean
# = Area(x * f(x)) / area(f(x))

### a.) Integrate
# integrate: does NOT work with Inf
integrate(area.tri.f, lower=0, upper=Inf, subdivisions=1000, a)
# some tricks in R:
# precision in R: seems *NOT* to use subdivisions = 1000
integrate(area.tri.f, lower=0, upper=300, a)
integrate(area.tri.f, lower=0, upper=300, subdivisions=1000, a)
# the tolerance needs to be lowered, too!
integrate(area.tri.f, lower=0, upper=300, subdivisions=1000, a)
integrate(area.tri.f, lower=0, upper=300, subdivisions=1000, rel.tol = 1e-8, a)
# selecting the interval length = 256 also improves precision & confidence!
integrate(area.tri.f, lower=0, upper=250, subdivisions=1000, a)
integrate(area.tri.f, lower=0, upper=256, subdivisions=1000, a)
integrate(area.tri.f, lower=100, upper=356, subdivisions=1000, a)

### Mean
area.tot = integrate(area.tri.f, lower=0, upper=300, subdivisions=1000, rel.tol = 1e-8, a)
area.norm = integrate(tri.f, lower=0, upper=300, subdivisions=1000, rel.tol = 1e-8, a)
area.tot
area.norm
# we can extract the value from integrate
# using $value
tri.mean = area.tot$value / area.norm$value
tri.mean

##############

### TODO: ...

### Cos Member Function
a = c(2,3,4)

cos.base.f = function(x, a) {
	# non-vectorized function
	if(x < a[1]) {
		return(0)
	} else if(x > a[3]) {
		return(0)
	} else if(x < a[2]) {
		return(cos(pi / 2 * (x - a[1]) / (a[2] - a[1])))
	} else {
		return(cos(pi / 2 * (a[3] - x) / (a[3] - a[2])))
	}
}
cos.f = function(x, a) {
	# vectorized function
	return(sapply(x, cos.base.f, a))
}

area.cos.f = function(x, a) {
	return(x * cos.f(x, a))
}

integrate(area.cos.f, lower=min(a), upper=max(a), subdivisions=10000, a)

# TODO:
# compute mean


#### Generalized ####
base.f = function(x, a, f) {
	# non-vectorized function
	if(x < a[1]) {
		return(0)
	} else if(x > a[3]) {
		return(0)
	} else if(x < a[2]) {
		return(f(pi / 2 * (x - a[1]) / (a[2] - a[1])))
	} else {
		return(f(pi / 2 * (a[3] - x) / (a[3] - a[2])))
	}
}
f.f = function(x, a, f) {
	# vectorized function
	return(sapply(x, base.f, a, f))
}

area.f = function(x, a, f=cos) {
	return(x * f.f(x, a, f))
}

# Cos(x)
integrate(area.f, lower=min(a), upper=max(a), subdivisions=10000, a, cos)
curve(f.f(x, a, cos), from=0, to=5)

# Cos(x) + 1/2


#### Distances ####

a = c(2,3,4)
b = c(5,6,7)
limits = c(0, 15)