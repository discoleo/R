
### Leonard Mada
###
### Gamma(1/n)
### Relations between the G(1/n) functions
###
### draft 0.2

# including:
# - Gamma(1/ (2*n)) = f(Gamma(1/n))
# - explicit formulas for G(1/6), G(1/10), G(1/12), G(1/14);
# - R code to test these relationships;

# Derivation:
# - formulas can be derived
#   using Legendre's duplication formula;
#   (or the equivalent higher order formulas)

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
# 3 more eq. needed!


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

