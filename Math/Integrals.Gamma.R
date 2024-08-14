
### Leonard Mada
###
### Gamma(1/n)
### Relations between the G(1/n) functions
###
### draft 0.4a

# including:
# - Gamma(1/ (2*n)) = f(Gamma(1/n))
# - explicit formulas for G(1/6), G(1/10), G(1/12), G(1/14);
# - R code to test these relationships;
# - Complex Gamma;

# Derivation:
# - formulas can be derived
#   using Legendre's duplication formula;
#   (or the equivalent higher order formulas)

#########################

###############
### History ###
###############

### draft v.0.3a - v.0.3b:
# - complex line integrals (experimental);


#########################

# Base formula:
# an example
k = 4 # parameter
#
n = 2*k + 1
gamma(1/(2*n)) # ==
gamma(1/n) * gamma(k/n) * 2^((n-1)/n) / sqrt(pi) * sin(k*pi/n)
#
gamma((2*n-1)/(2*n)) # ==
1 / (gamma(1/n) * gamma(k/n)) / 2^((n-1)/n) * sqrt(pi)^3 / (sin(pi/(2*n)) * sin(k*pi/n))


####################

####################
###
### Full Derivations
###

##########
### G(1/6)
gamma(1/6) # ==
gamma(1/3)^2 / sqrt(pi) * 2^(2/3) * sin(pi/3)
#
gamma(5/6)
1/gamma(1/3)^2 * sqrt(pi)^3 / 2^(2/3) / (sin(pi/3) * sin(pi/6))
# expanded sin()
1/gamma(1/3)^2 * sqrt(pi)^3 * 2^(4/3) / sqrt(3)


###########
### G(1/10)
gamma(1/10) # ==
gamma(1/5) * gamma(2/5) / sqrt(pi) * 2^(4/5) * sin(2*pi/5)
#
gamma(3/10)
gamma(1/5) / gamma(2/5) * sqrt(pi) / 2^(3/5) / sin(3*pi/10)
#
gamma(7/10)
gamma(2/5) / gamma(1/5) * sqrt(pi) * 2^(3/5)
#
gamma(9/10)
1/(gamma(1/5) * gamma(2/5)) * sqrt(pi)^3 / 2^(4/5) / (sin(pi/10) * sin(2*pi/5))


###########
### G(1/12)
gamma(1/12)
gamma(1/4) * gamma(1/3) / sqrt(pi) * 2^(5/4) * 3^(1/8) * sqrt(sin(pi/4) * sin(5*pi/12) * sin(pi/3))
#
gamma(5/12)
gamma(1/4) / gamma(1/3) * sqrt(pi) / 2^(1/4) * 3^(1/8) * sqrt(sin(pi/4) / (sin(5*pi/12) * sin(pi/3)))
#
gamma(7/12)
gamma(1/3) / gamma(1/4) * sqrt(pi) * 2^(1/4) / 3^(1/8) / sqrt(sin(5*pi/12) * sin(pi/4) / sin(pi/3))
#
gamma(11/12)
1/(gamma(1/4) * gamma(1/3)) * sqrt(pi)^3 / 2^(5/4) / 3^(1/8) / sqrt(sin(pi/4) * sin(5*pi/12) * sin(pi/3)) / sin(pi/12)

### derivation
gamma(1/12)*gamma(5/12) # ==
gamma(1/4) / gamma(3/4) * 2 * pi * 3^(1/4)
#
gamma(1/12) / gamma(5/12) # ==
gamma(1/6) / sqrt(pi) * 2^(5/6) * sin(5*pi/12)


###########
### G(1/14)
gamma(1/14)
gamma(1/7) * gamma(3/7) / sqrt(pi) * 2^(6/7) * sin(3*pi/7)
#
gamma(3/14)
gamma(2/7) * gamma(3/7) / sqrt(pi) / 2^(3/7) * sin(3*pi/7) / sin(3*pi/14)
# alternative
gamma(2/7) / gamma(4/7) * sqrt(pi) / 2^(3/7) / sin(3*pi/14)
#
gamma(5/14)
gamma(1/7) / gamma(2/7) * sqrt(pi) / 2^(5/7) / sin(5*pi/14)
#
gamma(9/14)
gamma(2/7) / gamma(1/7) * sqrt(pi) * 2^(5/7)
#
gamma(11/14)
1/(gamma(2/7) * gamma(3/7)) * sqrt(pi)^3 * 2^(3/7) / sin(3*pi/7) # * sin(3*pi/14) / sin(11*pi/14)
#
gamma(13/14)
1/(gamma(1/7) * gamma(3/7)) * sqrt(pi)^3 / 2^(6/7) / (sin(pi/14) * sin(3*pi/7))


###########
### G(1/18)
gamma(1/18)
gamma(1/9) * gamma(4/9) * 2^(8/9) / pi^(1/2) * sin(4*pi/9)
#
gamma(5/18)
gamma(5/9) / gamma(7/9) * 2^(4/9) * pi^(1/2)
#
gamma(7/18)
gamma(7/9) / gamma(8/9) * 2^(2/9) * pi^(1/2)
#
gamma(11/18)
gamma(2/9) / gamma(1/9) * 2^(7/9) * pi^(1/2)
gamma(4/9) / gamma(1/3) / 2^(2/9) / 3^(1/6) * pi^(1/2) / sin(2*pi/9)
#
gamma(13/18)
gamma(4/9) / gamma(2/9) * 2^(5/9) * pi^(1/2)
#
gamma(17/18)
1 / (gamma(1/9) * gamma(4/9)) / 2^(8/9) * pi^(3/2) / (sin(4*pi/9) * sin(pi/18))
##


#######
### ...


#######################
###
### Partial Derivations
###
### TODO

##########
### G(1/4)
gamma(1/4)
# ???


##########
### G(1/8)
gamma(1/8)
###
gamma(1/8) / gamma(3/8) # ==
gamma(1/4) * 2^(3/4) / sqrt(pi) * sin(3*pi/8)
# 1 more eq. needed!


#########
### G(1/9)
gamma(1/9)
###
gamma(1/9)*gamma(4/9) / gamma(2/9) # ==
gamma(1/3) * 2 * 3^(1/6) * sin(2*pi/9)
# 2 more eq. needed!


###########
### G(1/15)
gamma(1/15)
# derivation
gamma(1/15) / gamma(4/15) # ==
gamma(1/5) / gamma(2/5) * 2 * 3^(3/10) * sin(4*pi/15)
#
gamma(2/15) * gamma(7/15)
gamma(1/5) * gamma(2/5) * 2 * 3^(1/10) * sin(pi/5)
#
gamma(1/15) / gamma(2/15) * gamma(4/15) * gamma(7/15)
gamma(1/3)^2 * 2^2 * 5^(1/6) * sin(pi/3) * sin(2*pi/15)
# unfortuantely is the same as eq (3)
gamma(1/15) / gamma(2/15) * gamma(4/15) * gamma(7/15)
gamma(1/3)^2 / 2^2 * 5^(1/6) * sin(pi/3) / (sin(pi/15) * sin(4*pi/15) * sin(7*pi/15))
# prod(sin(c(1,2,4,7) * pi / 15)) * 2^4
#
# 1 more eq. needed!


#######################
###
### Missing Derivations
###
### TODO
### [more advanced mathematics is needed]

n = 3
n = 5
n = 7
# ...


#####################
#####################

### Beta
library(pracma)


### B(1/3, 1/3)
m = complex(re=cos(2*pi/3), im=sin(2*pi/3))

beta(1/3, 1/3) / 10;
integrate(function(x) x^(1/3)*(1-x)^(1/3), lower=0, upper=1)
integrate(function(x) (1-x^2)^(1/3), lower=0, upper=1)$value / 2^(2/3)
integrate(function(x) x^(5/3)*(1-x^2)^(-1/2), lower=0, upper=1)$value / 2^(2/3)
integrate(function(x) x^3*(1-x^3)^(-1/2), lower=0, upper=1)$value * 3/2 / 2^(2/3)


(line_integral(function(x) x^3*(1-x^3)^(-1/2), c(0, m)) +
	line_integral(function(x) x^3*(1-x^3)^(-1/2), c(0, m^2))) * -3/2 / 2^(2/3)
(line_integral(function(x) x^3*(1-x^3)^(-1/2), c(m, 1)) +
	line_integral(function(x) x^3*(1-x^3)^(-1/2), c(m^2, 1))) * 1/2 / 2^(2/3)
line_integral(function(x) x^3*(1-x^3)^(-1/2), c(m, m^2)) * 1/2 / 2^(2/3) * (m-m^2)

### TODO


### B(1/5, 1/5)
m5 = complex(re=cos(2*pi/5), im=sin(2*pi/5))
integrate(function(x) x^(1/5)*(1-x)^(1/5), lower=0, upper=1)
integrate(function(x) (1-x^2)^(1/5), lower=0, upper=1)$value / 2^(2/5)
integrate(function(x) x^5*(1-x^5)^(-1/2), lower=0, upper=1)$value * 5/2 / 2^(2/5)

(line_integral(function(x) x^5*(1-x^5)^(-1/2), c(0, m5)) +
	line_integral(function(x) x^5*(1-x^5)^(-1/2), c(0, m5^2)) +
	line_integral(function(x) x^5*(1-x^5)^(-1/2), c(0, m5^3)) +
	line_integral(function(x) x^5*(1-x^5)^(-1/2), c(0, m5^4))) * -5/2 / 2^(2/5)
(line_integral(function(x) x^5*(1-x^5)^(-1/2), c(1, m5)) +
	line_integral(function(x) x^5*(1-x^5)^(-1/2), c(1, m5^2)) +
	line_integral(function(x) x^5*(1-x^5)^(-1/2), c(1, m5^3)) +
	line_integral(function(x) x^5*(1-x^5)^(-1/2), c(1, m5^4))) * -1/2 / 2^(2/5)
line_integral(function(x) x^5*(1-x^5)^(-1/2), c(m5, m5^4)) * 1/2 / 2^(2/5) * (m5-m5^4) +
	line_integral(function(x) x^5*(1-x^5)^(-1/2), c(m5^2, m5^3)) * 1/2 / 2^(2/5) * (m5^2-m5^3)



#####################
#####################

#####################
### Complex Gamma ###

### Gamma(1i)
# 1. Maths 505: Let's calculate i factorial
#    https://www.youtube.com/watch?v=f4ORnke4rd8
# 2. BriTheMathGuy: |i Factorial| You Won't Believe The Outcome
#    https://www.youtube.com/watch?v=R7djflEwPwQ
#  - very introductory;
# 3. Wikipedia
#    https://en.wikipedia.org/wiki/Particular_values_of_the_gamma_function

abs(pracma::gammaz(1i))
sqrt(pi / sinh(pi))


### Gamma(1i/2)
abs(pracma::gammaz(1i/2))
sqrt(2*pi / sinh(pi/2))

# Derivation:
2 * sqrt(pi / sin(pi/2*1i))
2 * sqrt(-1i * pi / sinh(pi/2))
sqrt(pi / sinh(pi/2)) * sqrt(2)

#
- (pracma::gammaz(1i/2) + pracma::gammaz(-1i/2))
log(gamma(1/4)/gamma(3/4)) - log(pi)/4

# TODO: verify & develop


### Gamma(1i / k)
k = sqrt(7)
abs(pracma::gammaz(1i/k))
sqrt(k*pi / sinh(pi/k))

# TODO
