
#### Fuzzy Sets ####

### LAB 4: Fuzzy Rules, Degree of Fulfillment
### Ch 3 / page 28 [43]
###
### Leonard Mada

### draft 0.1

### (1) Degree of Fulfillment
### (2) Experimental: using R operators


#####################

### Generic Functions

### construct Triangular/Trapezoidal Fuzzy Numbers
fzbase.f = function(x, lim) {
	if(x < lim[1]) {
		return(0)
	} else if(x < lim[2]) {
		iLow = 1
		iTop = 2
		val = (x - lim[iLow]) / (lim[iTop] - lim[iLow])
		
		return(val)
	} else if(length(lim) == 4) {
		# Trapezoidal Set
		if(x < lim[3]) {
			return(1)
		}
		offset = 1
	} else {
		offset = 0
	}
	
	if(x < lim[3 + offset]) {
		iLow = 2 + offset
		iTop = 3 + offset
		val = 1 - (x - lim[iLow]) / (lim[iTop] - lim[iLow])
		
		return(val)
	} else {
		return(0)
	}
}
fz.f = function(x, lim) {
	# vectorized version of function
	return(sapply(x, fzbase.f, lim))
}

###################

#### Exercises ####

### 1.) Definition of Numbers

### Temperature Set
# Cold: Trapezoidal number
# 0: ...:-20;
# fz: -20:-10;
# 1: -10:0;
# fz: 0:15;
# 0: 15:...;
t.f = function(x) {
	return(fz.f(x, c(-20,-10,0,15)))
}
### Plot
curve(t.f, from=-30, to=40)


### Long Walk Set
# Long Walk: Triangular number
# 0: ...:200;
# fz: 200:1500;
# fz: 1500:4000;
# 0: 4000:...;
walk.f = function(x) {
	return(fz.f(x, c(200,1500,4000)))
}
### Plot
curve(walk.f, from=0, to=5000)


##############################

### 2.) Fuzzy Rules: Degree of fulfillment

# Exercise:
# T = 5 C, Dist = 500
# is solved below

#### DOF: Operators ####

###############
### 1.) NOT ###

not.dof = function(x, f) {
	return(1 - f(x))
}

# Example
not.dof(5, t.f)
# 0.333
not.dof(500, walk.f)
# 0.769


###############
### 2.) AND ###

and.dof = function(x, f1, f2) {
	dim = dim(x)
	if( ! is.null(dim) & length(dim) == 2) {
		# processing multiple values
		dof = f1(x[ , 1]) * f2(x[ , 2])
	} else {
		dof = f1(x[[1]]) * f2(x[[2]])
	}
	return(dof)
}

# Example
and.dof(c(5, 500), t.f, walk.f)
# 0.1538

# (NOT in the Lab, but may be practical)
# multiple values
vals = cbind(5:7, seq(500, 700, by=100))
and.dof(vals, t.f, walk.f)

### multiple args AND
and.multi.dof = function(x, fs) {
	len = length(fs)
	dof = sapply(1:len, function(id) fs[[id]](x[,id]))
	return(apply(dof, 1, prod))
}
vals = cbind(5:7, seq(500, 700, by=100))
and.multi.dof(vals, c(t.f, walk.f))


#######################
### 3.) AND (Logic) ###

and.logic.dof = function(x, fs) {
	len = length(fs)
	dim = dim(x)
	if(is.null(dim)) {
		dof = sapply(1:len, function(id) fs[[id]](x[id]))
		return(min(dof))
	} else if(length(dim) == 2) {
		# processing multiple values
		dof = sapply(1:len, function(id) fs[[id]](x[,id]))
	}
	return(apply(dof, 1, min))
}

# Example
and.logic.dof(c(5, 500), c(t.f, walk.f))
# 0.2307

# multiple values
vals = cbind(5:7, seq(500, 700, by=100))
and.logic.dof(vals, c(t.f, walk.f))


##############
### 4.) OR ###
# DOF = sum() - prod()

or.dof = function(x, f1, f2) {
	dim = dim(x)
	if( ! is.null(dim) & length(dim) == 2) {
		# processing multiple values
		dof = f1(x[ , 1]) + f2(x[ , 2]) - f1(x[ , 1]) * f2(x[ , 2])
	} else {
		dof = f1(x[1]) + f2(x[2]) - f1(x[1]) * f2(x[2])
	}
	return(dof)
}

# Example
or.dof(c(5, 500), t.f, walk.f)
# 0.744


######################
### 5.) OR (Logic) ###

or.logic.dof = function(x, fs) {
	len = length(fs)
	dim = dim(x)
	if(is.null(dim)) {
		dof = sapply(1:len, function(id) fs[[id]](x[id]))
		return(max(dof))
	} else if(length(dim) == 2) {
		# processing multiple values
		dof = sapply(1:len, function(id) fs[[id]](x[,id]))
	}
	return(apply(dof, 1, max))
}

# Example
or.logic.dof(c(5, 500), c(t.f, walk.f))
# 0.67

# multiple values
vals = cbind(5:7, seq(500, 700, by=100))
or.logic.dof(vals, c(t.f, walk.f))


###############
### 6.) XOR ###
# DOF = sum() - 2*prod()

xor.dof = function(x, f1, f2) {
	dim = dim(x)
	if( ! is.null(dim) & length(dim) == 2) {
		# processing multiple values
		dof = f1(x[ , 1]) + f2(x[ , 2]) - 2 * f1(x[ , 1]) * f2(x[ , 2])
	} else {
		dof = f1(x[1]) + f2(x[2]) - 2 * f1(x[1]) * f2(x[2])
	}
	return(dof)
}

# Example
xor.dof(c(5, 500), t.f, walk.f)
# 0.5897

# multiple values
vals = cbind(5:7, seq(500, 700, by=100))
xor.dof(vals, t.f, walk.f)


#######################
### 7.) XOR (Logic) ###
# DOF = max(min(1-f1, f2), (f1, 1-f2))

xor.logic.dof = function(x, fs) {
	# TODO: compute properly when >= 3 fs
	len = length(fs)
	dim = dim(x)
	if(is.null(dim)) {
		dof = sapply(1:len, function(id) fs[[id]](x[id]))
		# TODO
		tmp = c(min(1 - dof[[1]], dof[[2]]), min(dof[[1]], 1 - dof[[2]]))
		return(max(tmp))
	} else if(length(dim) == 2) {
		# processing multiple values
		dof = sapply(1:len, function(id) fs[[id]](x[,id]))
	}
	tmp = sapply(1:len, function(col) min(1 - dof[,col], dof[,-col]))
	return(apply(dof, 1, max))
}

# Example
xor.logic.dof(c(5, 500), c(t.f, walk.f))
# 0.67

# multiple values
vals = cbind(5:7, seq(500, 700, by=100))
xor.logic.dof(vals, c(t.f, walk.f))



########################

### EXPERIMENTAL
### with Operators in R
fz = function(x, f) {
	x.fz = list("x"=x, "f"=f)
	class(x.fz) = "fz"
	return(x.fz)
}
print.fz = function(fz) {
	print(fz$x)
}
'!.fz' = function(fz1) {
	# not quite correct
	# returns the confidence value (wrapped as a fuzzy number)
	return(fz(1 - fz1$f(fz1$x), fz1$f))
}
'&.fz' = function(fz1, fz2) {
	return(fz1$f(fz1$x) * fz2$f(fz2$x))
}
'|.fz' = function(fz1, fz2) {
	r1 = fz1$f(fz1$x)
	r2 = fz2$f(fz2$x)
	return(r1 + r2 - r1 * r2)
}
x.t = fz(5, t.f)
x.walk = fz(500, walk.f)
# print/compute various values
x.t
! x.t
x.t & x.walk
x.t | x.walk
# TODO: return ? fz-object ?

