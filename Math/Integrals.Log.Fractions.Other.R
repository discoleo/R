########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Logarithms
### Log-Fractions: Other
###
### draft v.0.1a


##################
### Logarithms ###
##################

# - definite integrals;
# - various types of Logarithms combined with fractions;


### Helper Constants
Catalan = 0.915965594177219015054603514;
# Note:
# Catalan = - I(log(x)/(x^2 + 1), lower=0, upper=1)


################
################

### I( log(x)*log(1-x^p) / x )
# Maths 505: A beautiful log-trig integral featuring an important constant
# https://www.youtube.com/watch?v=2bwvuWsQSEY
# [the intermediate integral]

###
integrate(\(x) log(x)*log(1-x^2)/x, 0, 1)
pracma::zeta(3) / 4

###
p = sqrt(5)
integrate(\(x) log(x)*log(1-x^p)/x, 0, 1)
pracma::zeta(3) / p^2


################
################

# Maths 505: Another INSANE integral!
# https://www.youtube.com/watch?v=KEDEzVqlAYU

integrate(function(x) atan(x)*log((1-x)/(1+x)), 0, 1)
pi^2/16 - pi/4*log(2) - Catalan;


# Gen 1: TODO
k = 2
integrate(function(x) atan(x)*log((k - x)/(k + x)), 0, 1)
(pi/4 - log(2)/2)*log((k-1)/(k+1)) +
	+ integrate(function(x) (x*atan(x) - log(x^2+1)/2)*(1/(k - x) + 1/(k + x)), 0, 1)$value
(pi/4 - log(2)/2)*log((k-1)/(k+1)) +
	+ integrate(function(x) k*atan(x)*(1/(k-x) - 1/(x+k)), 0, 1)$value +
	- integrate(function(x) log(x^2+1)/2*(1/(k - x) + 1/(k + x)), 0, 1)$value


###
integrate(function(x) atan(x) / x, 0, 1)
Catalan

# Varia:
x = exp(pracma::lambertWp(exp(-1)) / 2 + 1/2)
# Maximum of function:
log(x) / (x^2 + 1)
