
### Leonard Mada
###
### Integrals: Exercises
### draft 0.4a


### Various Exercises

# Students are invited to solve these exercises! ;-)



#####################

#########
### Ex. 1
### Michael Penn: A trigonometric integral.
### https://www.youtube.com/watch?v=V5fjlHh9pR8

### 1.a.) partial generalization
### 1/(a*cos(x)^4 + a*sin(x)^4 - 1)

a = 7

integrate(function(x) 1/(a*cos(x)^4 + a*sin(x)^4 - 1), lower=0, upper=pi/4)
# Exact solution
1/(a-1) / sqrt(2 - 2/(a-1)) * pi/2


a = 1/3
integrate(function(x) 1/(cos(x)^4 + sin(x)^4 - a), lower=0, upper=pi/4)
# Exact solution
1/(1-a) / sqrt(2 - 2/(1/a-1)) * pi/2
1/(1-a) / sqrt(2 - 2*a/(1-a)) * pi/2


a = 1/3
n = 4 # = integer multiples of pi/4!
integrate(function(x) 1/(cos(x)^4 + sin(x)^4 - a), lower=0, upper = n * pi/4)
# Exact solution
n/(1-a) / sqrt(2 - 2*a/(1-a)) * pi/2


### TODO: full generalization
a = 1/3
upper = pi/5 # ONLY between 0 & up to pi/2 !!!
integrate(function(x) 1/(cos(x)^4 + sin(x)^4 - a), lower=0, upper = upper)
# Exact solution
1/(1-a) / sqrt(2 - 2*a/(1-a)) * (pi/2 + atan( (tan(upper) - 1/tan(upper)) / sqrt(2 - 2*a/(1-a)) ))
# TODO: as multiples of pi/2!
# I[0, pi/2 + alfa] = Exact[0, pi/2] + Exact[0, alfa] # 0 is implicit;
# Exact[0, alfa] == Exact[0, pi/2 + alfa];


#########
### 1.b.) partial generalization
### sin(2*x)*cos(2*x)/(a*cos(x)^4 + a*sin(x)^4 - 1)

a = 7

integrate(function(x) sin(2*x)*cos(2*x)/(a*cos(x)^4 + a*sin(x)^4 - 1), lower=0, upper=pi/4)
# Exact solution
-1/(2*a) * log((a - 2)/(2*a - 2))
1/(2*a) * log((2*a - 2)/(a - 2))


#########
### 1.c.) partial generalization
### sin(x)*cos(x)/(a*cos(x)^4 + a*sin(x)^4 - 1)

a = 7

upper = pi/4
integrate(function(x) sin(x)*cos(x)/(a*cos(x)^4 + a*sin(x)^4 - 1), lower=0, upper=upper)
# Exact solution
coeff = sqrt((1-1/a)/2 - 1/4)
exact = function(x) -1/(4*a) * 1/coeff * atan((cos(x)^2 - 1/2)/coeff)
#
exact(upper) - exact(0)


#########
### 1.d.) partial generalization
### sin(x)^2*cos(x)^2/(a*cos(x)^4 + a*sin(x)^4 - 1)

a = 7

integrate(function(x) sin(x)^2*cos(x)^2/(a*cos(x)^4 + a*sin(x)^4 - 1), lower=0, upper=pi/4)
# Exact solution
-pi/(8*a) + (1-1/a)/(a-1) / sqrt(2 - 2/(a-1)) * pi/4
-pi/(8*a) + 1/a / sqrt(2 - 2/(a-1)) * pi/4



### 1.e.) TODO: variants using results from Ex. 2;

a = 1/3
upper = pi/4
integrate(function(x) (sin(x) + cos(x))/(cos(x)^4 + sin(x)^4 - a), lower=0, upper=upper)
# Exact solution:
# -2/(t^4 - 2*t^2 + 2*a - 1)
# 2/(sqrt(8 - 8*a)) * ( 1/(t^2 - r[1]) - 1/(t^2 - r[2]))
exact = function(a, limit) {
	a = complex(re=a, im=0)
	r = 1 + c(-1, 1)*sqrt(2 - 2*a)
	inv.sq = 1/sqrt(-r)
	t = sin(limit) - cos(limit)
	2/(sqrt(8 - 8*a)) * (inv.sq[1]*atan(t * inv.sq[1]) - inv.sq[2]*atan(t * inv.sq[2]))
}
exact(a, upper) - exact(a, 0)



################
### 1.B.) Powers

a = 1/5

integrate(function(x) 1/(cos(x)^6 + sin(x)^6 - a), lower=0, upper=pi/4)
# Exact solution
1/(1 - a) / sqrt((4*a-1)/(a - 1)) * pi/2


#################
#################


#########
### Ex. 2

# Roots of Unity
m = complex(re=cos(2*pi/3), im=sin(2*pi/3))

### 2.a.) 1/(x^4 - x^2 + 1)
upper = Inf
integrate(function(x) 1/(x^4 - x^2 + 1), lower=0, upper=upper)
# Exact solution
exact = function(x) 1/(m^2 - m) * (m*atan(x*m) - m^2*atan(x*m^2))
exact(upper) - exact(0)
# computed
exact(Inf)
pi/2



### 2.b.1.) 1/(x^6 - x^3 + 1)
# Exact solution:
# TODO - not yet; *** very ugly Partial Fraction Decomposition ***
exact = function(x) {
	mr = c(m^(1/3), (m^2)^(1/3)) # root 3 of root of unity;
	m.sum = m + m^2 # -1
	D = 2*mr - m.sum^2 * mr / 2
	m.shift = mr / 2
	m.sq = sqrt(3) * m.shift
	log.term = function(x, mr) log(x + mr) + m.sum/2 * log((x + m*mr)*(x + m^2*mr))
	atan.term = function(x, D, m.sq, shift) D/m.sq * atan((x - shift)/m.sq)
	#
	1/3 * 1/(m^2 - m) * (
	m^2 * mr[1] * ( log.term(x, mr[1]) + atan.term(x, D[1], m.sq[1], m.shift[1]) ) -
	m * mr[2] * ( log.term(x, mr[2]) + atan.term(x, D[2], m.sq[2], m.shift[2]) )
	)
}
# Test
upper = 100
integrate(function(x) 1/(x^6 - x^3 + 1), lower=0, upper=upper)
exact(upper) - exact(0)
# TODO: simplify/rewrite formula, so that Inf works as well;
# some computed values
# TODO: Inf, ...


### Fraction decomposition
# Test:
x = 2 # arbitrary value for testing;
#
mr = c(m^(1/3), (m^2)^(1/3))
m.sum = m+m^2
#
1/(x^6 - x^3 + 1)
#
1/3 * 1/(m^2 - m) *
	( m^2 * mr[1] * (1/(x + mr[1]) + (m.sum*x + 2*mr[1])/((x + m*mr[1])*(x + m^2*mr[1])) ) -
	m * mr[2] * (1/(x + mr[2]) + (m.sum*x + 2*mr[2])/((x + m*mr[2])*(x + m^2*mr[2])) ) )


### 2.b.2.) 1/(x^6 + x^3 + 1)
# Fraction decomposition
x = 2
#
mr3 = m^(1/3) # TODO: rewrite as robust formula!
m.sum = m+m^2
#
1/(x^6 + x^3 + 1)
#
1/3 * 1/(m - m^2) *
	( m^2 * (1/(x/mr3 - 1) + (m.sum*x/mr3 - 2)/((x/mr3 - m)*(x/mr3 - m^2)) ) -
	m * (1/(x/mr3^2 - 1) + (m.sum*x/mr3^2 - 2)/((x/mr3^2 - m)*(x/mr3^2 - m^2)) ) )


########################
########################


### Basic Exercises

### 1/sin(x)
lower = pi/5
upper = pi/2
integrate(function(x) 1/sin(x), lower=lower, upper=upper)
exact = function(x) 1/2 * log((1 - cos(x))/(1 + cos(x)))
exact(upper) - exact(lower)


### 1/(sin(x) + cos(x))
lower = 0
upper = pi/2
integrate(function(x) 1/(sin(x) + cos(x)), lower=lower, upper=upper)
exact = function(x) sqrt(2)/4 * log((1 - cos(x+pi/4))/(1 + cos(x+pi/4)))
exact(upper) - exact(lower)


### 1/(sin(x) + a * cos(x))
exact = function(x, a) {
	b = 1/sqrt(a^2 + 1)
	alfa = acos(b)
	b/2 * log((1 - cos(x + alfa))/(1 + cos(x + alfa)))
}
#
a = 3
lower = 0
upper = pi/2
integrate(function(x) 1/(sin(x) + a * cos(x)), lower=lower, upper=upper)
exact(upper, a) - exact(lower, a)


### 1/(sin(x) + a)
exact = function(x, a) {
	a.sq = sqrt(complex(re = 1 - a^2, im=0))
	# with fractions; altenative: atan;
	r1 = 1/2 * 1/a.sq * log((a.sq - cos(x))/(a.sq + cos(x)))
	# use result from 1/(sin(x)^2 - a^2);
	r2 = 1/2 * 1/a.sq * log((a.sq * tan(x) - a)/(a.sq * tan(x) + a))
	r1 - r2
}
#
lower = pi/5
upper = pi/2
a = 1/2
integrate(function(x) 1/(sin(x) + a), lower=lower, upper=upper)
# integrate(function(x) a/(sin(x)^2 - a^2), lower=lower, upper=upper)
exact(upper, a) - exact(lower, a)


### 1/(sin(x) + a * cos(x) + b)
exact = function(x, a, b) {
	c2 = 1/sqrt(a^2 + 1)
	# TODO: correct issues when a < 0!
	alfa = acos(c2) # atan2(a*c2, c2) # sign(a) * acos(c2)
	c1 = b * c2
	a.sq = sqrt(complex(re = 1 - c1^2, im=0))
	# print(a.sq)
	#
	# c2/2 * log((1 - cos(x + alfa))/(1 + cos(x + alfa)))
	# with fractions; alternative: atan;
	r1 = c2/2 * 1/a.sq * log((a.sq - cos(x + alfa))/(a.sq + cos(x + alfa)))
	# use result from 1/(sin(x)^2 - a^2);
	r2 = c2/2 * 1/a.sq * log((a.sq * tan(x + alfa) - c1)/(a.sq * tan(x + alfa) + c1))
	r1 - r2
	#
	r = c2/2 * 1/a.sq * log(
		(a.sq - cos(x + alfa)) / (a.sq + cos(x + alfa)) /
		(a.sq * tan(x + alfa) - c1)* (a.sq * tan(x + alfa) + c1))
}
# TODO: correct formula when a < 0!
a = 2
b = 3
lower = 0
upper = pi/2
integrate(function(x) 1/(sin(x) + a * cos(x) + b), lower=lower, upper=upper)
exact(upper, a, b) - exact(lower, a, b)


### 1/(sin(x)^2 - a^2) & 1/(sin(x)^2 + a^2)
a = 1/2
lower = pi/5
upper = pi/3
exact = function(x, a) {
	a.sq = sqrt(complex(re = 1 - a^2, im=0))
	1/(2*a) * 1/a.sq * log((a.sq * tan(x) - a)/(a.sq * tan(x) + a))
}
# Test
integrate(function(x) 1/(sin(x)^2 - a^2), lower=lower, upper=upper)
integrate(function(t) 1/ (t^2*(1-a^2) - a^2), lower=tan(lower), upper=tan(upper))
exact(upper, a) - exact(lower, a)


integrate(function(x) 1/(sin(x)^2 + a^2), lower=lower, upper=upper)
integrate(function(t) 1/ (t^2*(1+a^2) + a^2), lower=tan(lower), upper=tan(upper))
exact(upper, a * 1i) - exact(lower, a * 1i)


###############

### Basic Fractions
# see https://github.com/discoleo/R/blob/master/Math/Integrals.Fractions.Unity.R
# [the complete theory - still TODO]

unity = function(n, sum=TRUE) {
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	if(sum) {
		half = (n-1) %/% 2
		m.m = m^(1:half)
		m.m = cbind(m.m, 1/m.m)
		sum = m.m[,1] + m.m[,2]
		return(list("m"=m.m, "sum"=sum))
	}
	return(m)
}
### ONLY n = ODD!
exact = function(x, n=3) {
	m.all = unity(n)
	m = m.all$m
	m.sum = m.all$sum
	b0 = 1/n
	a = b0 * m.sum
	m.shift = m.sum/2
	D = -2*b0 + a*m.shift
	m.sq = sqrt(1 - m.shift^2)
	#
	b0 * log(-tan(x) - 1 + 0i) +
	sum(a/2*log((tan(x) + m[,1])*(tan(x) + m[,2]))) +
	sum(D / m.sq * atan(-(tan(x) + m.shift)/m.sq))
}
pseudo.exact = function(n, x.pow, lower, upper) {
	# faking now: x^x.pow / (x^n + 1) dx
	# for exact solutions:
	# see https://github.com/discoleo/R/blob/master/Math/Integrals.Fractions.Unity.R
	integrate(function(x) tan(x)^x.pow/cos(x)^2 / (tan(x)^n + 1), lower=lower, upper=upper)
}
log.exact = function(lower, upper, n) {
	log((sin(upper)^n + cos(upper)^n) / (sin(lower)^n + cos(lower)^n)) / n
}

##########
### cos(x) / (sin(x)^3 + cos(x)^3)
#
n = 3
lower = pi/10
upper = pi/2
integrate(function(x) cos(x) / (sin(x)^3 + cos(x)^3), lower=lower, upper=upper)
exact(upper, n) - exact(lower, n)


integrate(function(x) sin(x) / (sin(x)^3 + cos(x)^3), lower=lower, upper=upper)
- exact(pi/2 - upper, n) + exact(pi/2 - lower, n)


####################
### Powers 5, 7, ...
### (trivial)

### n = 5
n = 5
lower = pi/10
upper = pi/2
integrate(function(x) cos(x)^3 / (sin(x)^5 + cos(x)^5), lower=lower, upper=upper)
exact(upper, n) - exact(lower, n)



### n = 7
n = 7
lower = pi/10
upper = pi/2
integrate(function(x) cos(x)^5 / (sin(x)^7 + cos(x)^7), lower=lower, upper=upper)
exact(upper, n) - exact(lower, n)


### Workout ALL subtypes
### A[i, j] = sin(x)^i * cos(x)^j / (sin(x)^n + cos(x)^n)

### n = 5
n = 5
lower = pi/10
upper = pi/2

### A[3,0] & A[0, 3]
integrate(function(x) sin(x)^3 / (sin(x)^5 + cos(x)^5), lower=lower, upper=upper)
pseudo.exact(n, 3, lower, upper) # see also different solution above
#
integrate(function(x) cos(x)^3 / (sin(x)^5 + cos(x)^5), lower=lower, upper=upper)
exact(upper, n) - exact(lower, n)

### A[1, 2] & A[2, 1]
integrate(function(x) sin(x)*cos(x)^2 / (sin(x)^5 + cos(x)^5), lower=lower, upper=upper)
pseudo.exact(n, 1, lower, upper)
#
integrate(function(x) sin(x)^2*cos(x) / (sin(x)^5 + cos(x)^5), lower=lower, upper=upper)
pseudo.exact(n, 2, lower, upper)

### A[1, 0] & A[0, 1] = (A[1, 2] + A[3, 0], A[2, 1] + A[0, 3])
integrate(function(x) sin(x) / (sin(x)^5 + cos(x)^5), lower=lower, upper=upper)
pseudo.exact(n, 1, lower, upper)$value + pseudo.exact(n, 3, lower, upper)$value
#
integrate(function(x) cos(x) / (sin(x)^5 + cos(x)^5), lower=lower, upper=upper)
pseudo.exact(n, 2, lower, upper)$value + pseudo.exact(n, 0, lower, upper)$value


### TODO:
### A[2, 2] = (A[4,0] + A[0,4])


### A[5, 0] + A[0, 5] = x
# sin(x)^5 + cos(x)^5 =
# = (sin(x) + cos(x))*(1 - 2*A[2,2] - A[1,1] + A[2,2])
# = (sin(x) + cos(x))*(1 - A[2,2] - A[1,1])
# = A[1,0] + A[0,1] - A[3,2] - A[2,3] - A[2,1] - A[1,2]
# = (A[1,0] + A[0,1] - A[2,1] - A[1,2]) - A[3,2] - A[2,3]
# Test
integrate(function(x) (cos(x)*sin(x)^4 - sin(x)*cos(x)^4) / (sin(x)^5 + cos(x)^5), lower=lower, upper=upper)
log.exact(lower, upper, n)
# A[3, 2] & A[2, 3]
# A[3, 2]
integrate(function(x) sin(x)^3*cos(x)^2 / (sin(x)^5 + cos(x)^5), lower=lower, upper=upper)
A32 = 1/2 * (log.exact(lower, upper, n) +
		pseudo.exact(n, 1, lower, upper)$value - pseudo.exact(n, 2, lower, upper)$value -
	(upper - lower - pseudo.exact(n, 3, lower, upper)$value - pseudo.exact(n, 0, lower, upper)$value) )
A32
# A[2, 3]
integrate(function(x) sin(x)^2*cos(x)^3 / (sin(x)^5 + cos(x)^5), lower=lower, upper=upper)
A23 = (A32 - log.exact(lower, upper, n) +
	pseudo.exact(n, 2, lower, upper)$value - pseudo.exact(n, 1, lower, upper)$value)
A23



###########################
### Fractions with Radicals

### cos(x) / ( sqrt(sin(x)^2 + 1) * (1 + sin(x)) )
lower = pi/10
upper = pi/2
integrate(function(x) cos(x) / ( sqrt(sin(x)^2 + 1) * (1 + sin(x)) ), lower=lower, upper=upper)
exact = function(x) sqrt(2)/4 * log((1 - cos(atan(sin(x))+pi/4))/(1 + cos(atan(sin(x))+pi/4)))
exact = function(x) sqrt(2)/4 * log(
	(sqrt((sin(x)^2 +1)*2) - 1 + sin(x)) /
	(sqrt((sin(x)^2 +1)*2) + 1 - sin(x)))
exact(upper) - exact(lower)

