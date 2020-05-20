
### Leonard Mada
###
### Integrals: Exercises
### draft 0.3


### Various Exercises

# Students are invited to solve these exercises1 ;-)



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
