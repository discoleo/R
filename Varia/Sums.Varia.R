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

##################

library(pracma)


##################

##################
### 3D Variant ###


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


##################

##################
### 2D Variant ###


### 1 / (x1*x2 * (x1 + x2))
sum.fr2v2 = function(iter) {
	sum = 0
	for(i1 in 1:iter) {
		sum = sum + sum(1 / ((iter:1) * (i1 + iter:1)) / i1 );
		if(i1 %% 400 == 1) {
			print(i1)
		}
	}
	return(sum)
}


### 1 / (x1^2*x2^2 * (x1 + x2))
sum.fr2v21 = function(iter) {
	sum = 0
	for(i1 in 1:iter) {
		sum = sum + sum(1 / (iter:1)^2 / ((i1 + iter:1)) / i1^2 );
		if(i1 %% 400 == 1) {
			print(i1)
		}
	}
	return(sum)
}


### 1 / (x1*x2 * (x1 + x2)^2)
sum.fr2v12 = function(iter) {
	sum = 0
	for(i1 in 1:iter) {
		sum = sum + sum(1 / (iter:1) / ((i1 + iter:1)^2) / i1 );
		if(i1 %% 400 == 1) {
			print(i1)
		}
	}
	return(sum)
}

### 1 / (x1^2*x2 * (x1 + x2))
sum.fr2vas21 = function(iter) {
	sum = 0
	for(i1 in 1:iter) {
		sum = sum + sum(1 / (iter:1) / ((i1 + iter:1)) / i1^2 );
		if(i1 %% 400 == 1) {
			print(i1)
		}
	}
	return(sum)
}

###########
iter = 5000

s22 = sum.fr2v2(iter)
print(s22, 12) / 2
print(zeta(3), 12)


s221 = sum.fr2v21(iter)
print(s221, 12)


s212 = sum.fr2v12(iter)
print(s212, 12)


s2as21 = sum.fr2vas21(iter)
print(s2as21, 12)


##################

##################
### 1D Variant ###

sum.1D0 = function(b, start=1, iter=30000) {
	sum(1 / ((start:iter)^2 + b^2)) # Note: b^2;
}

sum.1D = function(b, start=1, iter=30000) {
	len = length(b) - 1; pow = 0:len;
	sum(sapply(start:iter, function(x) 1 / sum(x^pow * b) ))
}

sum.1Dsq.exact = function(b) {
	pi / (2*b) * 1/tanh(pi*b) - 1 / (2*b*b)
}
sum.1D.exact = function(b0, b1=0, useDet=FALSE) {
	check = function(b0, b1) {
		if( ! (round(b1) == b1) ) return(FALSE);
		if(b1 %% 2 != 0 && b0 != 0) {
			det = b1^2 - 4*b0;
			if(det < 0) return(FALSE);
			det = sqrt(det);
			if(det != round(det)) return(FALSE);
			# if( (-b1 + det) %% 2 != 0 ) return(FALSE);
		}
		return(TRUE);
	}
	if( ! check(b0, b1) ) {
		stop("b1 MUST be an even integer!")
	}
	### b0 == 0
	if(b0 == 0) {
		s = sum(1/(1:b1)) / b1;
		return(s);
	}
	### b1 %% 2 != 0
	if(useDet || b1 %% 2 != 0) {
		det = sqrt(b1^2 - 4*b0);
		r1 = (-b1 - det)/2; r2 = (-b1 + det)/2;
		val = 1:det;
		s = sum(1 / (val - r2)) / det;
		return(s);
	}
	### b1 %% 2 == 0
	b0.tr = b0 - b1^2/4;
	b0.tr.sqrt = sqrt(b0.tr + 0i);
	s = pi / (2*b0.tr.sqrt) * 1/tanh(pi*b0.tr.sqrt) - 1 / (2*b0.tr);
	if(b1 == 0) return(s);
	s - sum(1/((1:(b1/2))^2 + b0.tr))
}

###########

iter = 80000

###
sum.1D0(1, iter=iter)
sum.1D(c(1^2,0,1), iter=iter)
sum.1Dsq.exact(1)


###
sum.1D0(sqrt(2), iter=iter)
sum.1D(c(2,0,1), iter=iter)
sum.1Dsq.exact(sqrt(2))


### b1 %% 2 == 0
b0 = 11
b1 = 6 # must be even
sum.1D(c(b0, b1, 1), iter=iter)
sum.1D.exact(b0, b1)


### b0 == 0
b0 = 0
b1 = 6 # must be integer
sum.1D(c(b0, b1, 1), iter=iter)
sum.1D.exact(b0, b1)


### decomposable Fraction
b0 = 2
b1 = 3
sum.1D(c(b0, b1, 1), iter=iter)
sum.1D.exact(b0, b1)

#
b0 = 4
b1 = 5
sum.1D(c(b0, b1, 1), iter=iter)
sum.1D.exact(b0, b1)

#
b0 = 5/4
b1 = 3
sum.1D(c(b0, b1, 1), iter=iter)
sum.1D.exact(b0, b1)

# non-standard
b0 = 0
b1 = 1/3
sum.1D(c(b0, b1, 1), iter=iter)
sum.1D.exact(b0, b1) # NOT yet computable
sum.1D.exact(-1/9, 0) + 9/2*(1 - log(3))
### TODO: full solution;
