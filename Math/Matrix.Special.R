

####################
###
### Special Matrices
###
### Leonard Mada
###
### draft v.0.10

####################

### Special Matrices

# Matrices with det(M) = 0
# - Test calculations with large numbers;
# - test algorithms to improve precision;
# Conditions:
# - for "Power" matrices: if dim(M) >= pow + 2;

# [currently only older stuff - but unpublished]


#############

### Generator Functions

### a[i,j] = n[i,j]^pow + shift
m.pow.gen = function(ord, pow=2, offset=1, shift=0) {
	maxE = offset - 1 + ord^2
	m = matrix((offset:maxE)^pow + shift, nrow=ord, byrow=T)
}

### a[i,j] = (i*j)^pow + shift
xymatrix.gen = function(sz, shift=0, pow=1, offset=c(0,0)) {
	m.x = (offset[1] + 1):(offset[1] + sz)
	m.y = (offset[2] + 1):(offset[2] + sz)
	mx = outer(m.x, m.y, function(x, y) (x*y)^pow + shift)
	return(mx)
}

### TODO: ...
# Exercises for a great Textbook in Mathematics!
zdet.gen = function(n, sz=3, type=c("xy")) {
	# generates n matrices of size=sz with det = 0
	# TODO
}

### Matrix Reductions

# these Algorithms are O(n!);
# Is someone willing to unfold algorithms & implement optimized versions?

# 1D version: for simple powers
preduce.m = function(m) {
	d = dim(m)
	for(j in c(2:d[2])) {
		for(i in d[1]:j) {
			m[ , i] = m[ , i] - m[ , i-1]
		}
	}
	return(m)
}

# 2D version: for the xy-Matrices
reduce.m = function(m) {
	m = reduce.1D.m(m)
	m = reduce.1D.m(t(m))
	return(m)
}
reduce.1D.m = function(m) {
	ord = dim(m)[[2]]
	if(ord > 3) {
		for(j in c(2:(ord-2))) {
			for(i in ord:j) {
				m[ , i] = m[ , i] - m[ , i-1]
			}
		}
	} else if(ord == 3) {
		for(j in c(2, 2)) {
			for(i in ord:2) {
				m[ , i] = m[ , i] - m[ , i-1]
			}
		}
	}
	
	m
}

##################

##################
### Power-Matrices

### Squares

### Squares, 3x3: NON-Zero
pow = 2
ord = 3 # ord < (pow + 2) - 4
m = m.pow.gen(ord, pow)
m
# Test
det(m)
t(m)
det(t(m))

### Squares, 4x4
pow = 2
ord = 4
m = m.pow.gen(ord, pow)
m
# Test
det(m)
t(m)
det(t(m))
#
m = preduce.m(m)
det(m)

### Cubes

### Cubes, 4x4: NON-Zero
ord = 4
m = m.pow.gen(ord, 3)
m
det(m)

### Cubes, 5x5
ord = 5
m = m.pow.gen(ord, 3)
m
# Test
det(m) # very BIG ERROR!
#
m = preduce.m(m)
det(m)

### Cubes, 6x6
ord = 6
m = m.pow.gen(ord, 3)
m
# test
det(m)
#
m = preduce.m(m)
det(m)

### Cubes, 6x6
ord = 6
m = m.pow.gen(ord, 3, shift=-2^3)
m
# test
det(m)
#
m = preduce.m(m)
det(m)


#############
### 4th Power

### Pow 4, 6x6: Zero
ord = 6
m = m.pow.gen(ord, 4)
m
# Test
det(m) # HUGE ERROR!
#
m = preduce.m(m)
det(m)


### Pow 4, 7x7: Zero
ord = 7
m = m.pow.gen(ord, 4)
m
# Test
det(m) # NON-TRIVIAL ERROR!
#
m = preduce.m(m)
det(m)


########################

###############
### xy-Matrices

### x*y + shift
### 3x3
sz = 3
m = xymatrix.gen(sz, shift=-1)
m
det(m)

### 3x3
sz = 3
m = xymatrix.gen(sz, shift=1)
m
det(m)

### 5x5
sz = 5
m = xymatrix.gen(sz, shift=1, pow=1)
m
det(m)

###########
### Pow = 2
### (x*y)^2 + shift

### 5x5
sz = 5
m = xymatrix.gen(sz, shift=1, pow=2)
m
# Test
det(m)
mx = reduce.m(m)
mx
det(mx)


###########
### Pow = 3
### (x*y)^3 + shift

### 5x5
sz = 5
m = xymatrix.gen(sz, shift=3, pow=3, offset=c(2,0))
m
# Test
det(m)
#
mx = reduce.1D.m(m)
mx
det(mx) # [usually not sufficient]
mx = reduce.m(m)
mx
det(mx)


### 5x5
sz = 5
m = xymatrix.gen(sz, shift=3, pow=3, offset=c(2,2))
m
# Test
det(m)
#
mx = reduce.1D.m(m)
mx
det(mx) # not sufficient
mx = reduce.m(m)
mx
det(mx)


######
sz = 5
m = xymatrix.gen(sz, shift=-58826, pow=3, offset=c(2,2))
m
# Test
det(m)
#
mx = reduce.1D.m(m)
mx
det(mx) # not sufficient
mx = reduce.m(m)
mx
det(mx)

