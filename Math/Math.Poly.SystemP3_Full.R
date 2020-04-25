
### Leonard Mada
###
### P3 Polynomial Systems
### Solver: Exact solutions
###
### draft 0.1


### P3 System:
# x^3 + y^3 + z^3 = A
# x*y + x*z + y*z = B
# x*y*z = P


### Steps:

### Subsystem 1:
#
# X^3 - 3*B*X + 3*P - A = 0
# solve for X;

### Subsystem 2:
#
# x + y + z = X
# x*y + x*z + y*z = B
# x*y*z = P
# solve for x, y, z;


### helper functions

# dir.create(tempdir())

# safe sqrt
sqrt.c = function(x) {
	x = as.complex(x)
	x.sqrt = sqrt(x)
	return(x.sqrt)
}
# safe n-th root
nroot.c = function(x, n=3) {
	if(Im(x) == 0) {
		x = Re(x)
		if(x >= 0) {
			x = complex(re = x^(1/n), im=0)
		} else if(n %% 2 == 0) {
			x = complex(re = 0, im=(-x)^(1/n))
		} else {
			x = - (-x)^(1/n)
		}
	} else {
		x = x^(1/n)
	}
	return(x)
}
### complex/matrix round
round0 = function(m, tol=1E-10) {
	m[abs(Re(m)) < tol & abs(Im(m)) < tol] = 0
	
	isZero = (Re(m) != 0) & (abs(Re(m)) < tol)
	if(sum(isZero) > 0) {
		m[isZero] = complex(re=0, im=Im(m[isZero]))
	}
	return(m)
}

### Solution
# Cubic Polynomial
solveP3 <- function(c, d, n=3, all=FALSE) {
	det = (d^2 - c^n)
	det = sqrt.c(det)
	
	p = nroot.c(d + det, n)
	q = nroot.c(d - det, n)
	if(all) {
		m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
		id = 0:(n-1)
		m = m^id
		return(p*m + q/m)
	}
	return(p + q)
}
# Cubic System
solveP3S <- function(s, bi, pr, n=3) {
	if( n== 1) {
		d2 = s
	} else if(n == 2) {
		d2 = sqrt(s + 2*bi)
	} else if(n == 3) {
		# y^3 - 3*bi * y + 3*pr - s
		d = -(3*pr - s) / 2
		d2 = solveP3(bi, d)
		# Test
		y = d2
		err = y^3 - 3*bi * y + 3*pr - s
		print(paste("err Tri =", round0(err)))
	} else {
		print("Not yet implemented!")
		return(0)
	}
	
	### Solve P3 System:
	# y^3 - 3*B*y^2 + pr*d2*y - pr^2
	# y = shifted by B = bi/3
	# y^3 - 3*B^2*y + pr*d2*y - 2*B^3 + pr*d2*B - pr^2
	B = bi / 3
	d = B^3 - pr*(d2*B - pr)/2
	y = solveP3(B^2 - pr*d2/3, d, all=TRUE)
	print(paste("y =", round(y, 5)))
	y = y + B
	# Test
	err = y^3 - 3*B*y^2 + pr*d2*y - pr^2
	print(paste("err =", round0(err), collapse=", "))
	# x1, x2, x3
	x1 = pr/y
	syz = d2 - x1
	diff = sqrt.c(syz^2 - 4*y)
	x2 = (syz + diff)/2
	x3 = (syz - diff)/2
	p = matrix(c(x1, x2, x3), ncol=3)
	# Test
	m = 1
	test = c(
		apply(p, m, function(p) sum(p^n)),
		apply(p, m, function(p) sum(p[1] * p[-1], p[2]*p[3])),
		apply(p, m, function(p) prod(p)))
	return(list("x"=p, "test" = test))
}

############################

### Test

solveP3S(1,2,3)

solveP3S(1,2,-1)

# seems to compute all 9 roots;

