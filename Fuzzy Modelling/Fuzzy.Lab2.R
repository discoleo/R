
#### Fuzzy Sets ####

### LAB 2: Linear, Triangular Numbers
### Ch 2.1 / page 6 [10]
###
### Leonard Mada

### draft 0.1


### Membership Functions

### 1. Linear "fuzziness"
fz.f = function(x, lower, upper) {
	if(x <= lower) {
		return(0)
	} else if(x > upper) {
		return(1)
	}
	p = (x - lower) / (upper - lower)
	return(p)
}

#### Specific Examples ####

### Tall: between 175-185 cm
tall.fz = function(x) {
	# sapply used for compatibility with curve()
	r = sapply(x, fz.f, 175, 185)
	return(r)
}

### Very Tall
very.tall.fz = function(x) {
	r = tall.fz(x)
	return(r*r)
}

### "More or Less" Tall
moreless.tall.fz = function(x) {
	r = tall.fz(x)
	return(sqrt(r))
}


### PLOT

LIMIT_MAX = 250

curve(tall.fz, from=0, to=LIMIT_MAX)

curve(very.tall.fz, from=0, to=LIMIT_MAX, add=T, col="red")
curve(moreless.tall.fz, from=0, to=LIMIT_MAX, add=T, col="green")

#################

# TODO:
# Triangular Numbers
# cos


#################

#### EXPERIMENTAL ####

### TANH is frequently used in AI / Neural Networks

tanh(3)
curve(tanh, from=-5, to=5)

### Problem 1:
# however TANH runs between -1 and 1
# so it needs centering around 1/2 and rescaling by 1/2
# (1 + tanh()) / 2

### Problem 2:
# tanh is not continuous/smooth with -1 and +1 over a finite domain;
# therefore, the rescaled function is not continuous with 0 and 1;
# Solution for lower limit:
# p * (1 + tanh(p - 1/2)) / 2
# Solution for upper limit:
# TODO

fz.h.f = function(x, a, b, xscale=2) {
	if(x <= a) {
		return(0)
	} else if(x > b) {
		return(1)
	}
	p = (x - a) / (b - a)
	mid = 1/2
	p = p * (1 + tanh(xscale*(p - mid))) / 2
	# TODO
	# p is NOT continuous with 1
	return(p)
}

tall.h.fz = function(x) {
	r = sapply(x, fz.h.f, 175, 185)
	return(r)
}

curve(tall.h.fz, from=100, to=LIMIT_MAX)

