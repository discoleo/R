########################
###
### Leonard Mada
### [the one and only]


#################

### Elementary Zeta Functions:
# - extensions of the zeta function in higher dimensions;


### Sums of Fractions
# based on:
# Michael Penn: a nice double sum.
# https://www.youtube.com/watch?v=5KpGSMyUANU

### TODO:
# - explore relationships to zeta function;


### 1 / (x1*x2*x3 * (x1 + x2 + x3))
sum.fr3v3 = function(iter) {
	sum = 0
	for(i1 in 1:iter) {
		for(i2 in 1:iter) {
			sum = sum + sum(1 / ((iter:1) * (i1 + i2 + iter:1))) / (i1*i2);
		}
		if(i1 %% 200 == 1) {
			print(i1)
		}
	}
	return(sum)
}

### 1 / (x1*x2*x3 * (x1 + x2)*(x1 + x3)*(x2 + x3))
sum.fr3v2 = function(iter) {
	sum = 0
	for(i1 in 1:iter) {
		for(i2 in 1:iter) {
			sum = sum + sum(1 / (iter:1) /
				((i1 + iter:1)*(i2 + iter:1))) / (i1*i2) / (i1 + i2)
		}
		if(i1 %% 200 == 1) {
			print(i1)
		}
	}
	return(sum)
}

### 1 / (x1*x2*x3 * (x1 + x2)*(x1 + x3)*(x2 + x3) * (x1 + x2 + x3))
sum.fr3vall = function(iter) {
	sum = 0
	for(i1 in 1:iter) {
		for(i2 in 1:iter) {
			sum = sum + sum(1 / ((iter:1)*(i1 + i2 + iter:1)) /
				((i1 + iter:1)*(i2 + iter:1))) / (i1*i2) / (i1 + i2)
		}
		if(i1 %% 200 == 1) {
			print(i1)
		}
	}
	return(sum)
}

### 1 / ((x1 + x2)*(x1 + x3)*(x2 + x3))
sum.fr3v0 = function(iter) {
	sum = 0
	for(i1 in 1:iter) {
		for(i2 in 1:iter) {
			sum = sum + sum(1 /
				((i1 + iter:1)*(i2 + iter:1))) / (i1 + i2)
		}
		if(i1 %% 200 == 1) {
			print(i1)
		}
	}
	return(sum)
}

###
iter = 3000

s1 = sum.fr3v3(iter)
# s1 = 6.41364100193
print(s1, digits=12)


s2 = sum.fr3v2(iter)
# s2 = 0.392617632782
print(s2, digits=12)


s3 = sum.fr3vall(iter)
# s3 = 0.0884001690925
print(s3, digits=12)


s0 = sum.fr3v0(iter)
# s0 = 15.0284014084
print(s0, digits=12)


