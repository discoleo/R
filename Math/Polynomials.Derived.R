
########################
###
### Leonard Mada
### [the one and only]
###
### Derived Polynomials

### Note:
# This is the 1st part towards introducing polynomials of Class 1,
# including a different approach to polynomials.


### Theory:
# - let P[n] be a polynomial of order n with integer/rational coefficients;
# - let r be the n roots of this polynomial;
# - let f be a polynomial function with integer coefficients;
# - then f(r) are the roots of a polynomial of order n with integer coefficients;


library(polynom)

# needed to get the roots of the base polynomial;


### Examples:

### Base polynomial:
p = polynomial(c(-1,-1,0,0,0,1))
p

x0 = solve(p)
x0

### Derived polynomials

# lets create various derived polynomials of order 5 with integer coefficients;


###
x = x0^3 + x0^2 + x0
x

poly.calc(x)
err = -4 - 14*x - 21*x^2 - 11*x^3 + x^5
round0(err)

###
x = x0^4 - x0^2 + x0^3 - x0
x

poly.calc(x)
err = -1 + 3*x - 16*x^2 + 18*x^3 - 4*x^4 + x^5
round0(err)

###
x = x0^3 - x0
x

poly.calc(x)
err = -1 + 5*x - 8*x^2 + 4*x^3 + x^5
round0(err)

###
x = x0^3 + x0
x

poly.calc(x)
err = -1 - 5*x - 8*x^2 - 4*x^3 + x^5
round0(err)

###
x = x0^4 + x0^3 + x0^2 + x0
x

poly.calc(x)
err = -1 - 5*x - 10*x^2 - 10*x^3 - 4*x^4 + x^5
round0(err)

###
x = x0^4 + x0
x

poly.calc(x)
err = -4 - 2*x + 7*x^2 + x^3 - 4*x^4 + x^5
round0(err)

###
# works because free term of the initial polynomial is +/- 1!
x = x0 + 1/x0
x

poly.calc(x)
err = -1 + 4*x - 4*x^2 - 5*x^3 + x^4 + x^5
round0(err)

###
# works because free term of the initial polynomial is +/- 1!
x = 1/x0 - x0
x

poly.calc(x)
err = 1 + 4*x + 4*x^2 + 5*x^3 + x^4 + x^5
round0(err)
