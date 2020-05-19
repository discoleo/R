
### Leonard Mada
###
### Integrals: exercises
### draft 0.2-pre-0.2


### Various Exercises

#####################

#########
### Ex. 1
###  Michael Penn: A trigonometric integral.
### https://www.youtube.com/watch?v=V5fjlHh9pR8

### 1.a.) partial generalization
### 1/(a*cos(x)^4 + a*sin(x)^4 - 1)

a = 7

integrate(function(x) 1/(a*cos(x)^4 + a*sin(x)^4 - 1), lower=0, upper=pi/4)
# Exact solution
1/(a-1) / sqrt(2 - 2/(a-1)) * pi/2


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


###############


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
upper = Inf
integrate(function(x) 1/(x^6 - x^3 + 1), lower=0, upper=upper)
# Exact solution:
# TODO - not yet; *** very ugly Partial Fraction Decomposition ***
exact = function(x) {
	# b0 = 1/3
	mr3 = m^(1/3) # root 3 of root of unity;
	m.sum = m+m^2
	m.shift = # TODO
	D1 = 2*mr3 + m.sum*mr3/2 # TODO: /2 or NOT?
	D2 = 2*mr3^2 + m.sum*mr3^2/2
	m.sq = # TODO
	#
	1/3 * 1/(m^2 - m) *
	(m^2 * mr3 * ( log(x + mr3) + m.sum/2 * log((x + m*mr3)*(x + m^2*mr3)) + D1/m1.sq * atan((x-m1.shift)/m1.sq) ) -
	m * mr3^2 * ( log(x + mr3^2) + m.sum/2 * log((x + m*mr3^2)*(x + m^2*mr3^2)) + D2/m2.sq * atan((x-m2.shift)/m2.sq)) )
}
exact(upper) - exact(0)
# computed
exact(Inf)
# TODO

### Fraction decomposition
# Test:
x = 2 # arbitrary value for testing;
#
mr3 = m^(1/3)
m.sum = m+m^2
#
1/(x^6 - x^3 + 1)
#
1/3 * 1/(m^2 - m) *
	( m^2 * mr3 * (1/(x + mr3) + (m.sum*x + 2*mr3)/((x + m*mr3)*(x + m^2*mr3)) ) -
	m * mr3^2 * (1/(x + mr3^2) + (m.sum*x + 2*mr3^2)/((x + m*mr3^2)*(x + m^2*mr3^2)) ) )


### 2.b.2.) 1/(x^6 + x^3 + 1)
# Fraction decomposition
x = 2
#
mr3 = m^(1/3)
m.sum = m+m^2
#
1/(x^6 + x^3 + 1)
#
1/3 * 1/(m - m^2) *
	( m^2 * (1/(x/mr3 - 1) + (m.sum*x/mr3 - 2)/((x/mr3 - m)*(x/mr3 - m^2)) ) -
	m * (1/(x/mr3^2 - 1) + (m.sum*x/mr3^2 - 2)/((x/mr3^2 - m)*(x/mr3^2 - m^2)) ) )
