#####################
##
## Leonard Mada
## [the one and only]
##
## Polynomial Radicals:
## Other Types


### Note:
# - Basic types are in file:
#   Integrals.Fractions.Unity.Radicals.R;


##########################
##########################

### I( sqrt(sqrt(1 + x^4) + 1) )
# Maths 505: A deceivingly tough integral
# https://www.youtube.com/watch?v=nGx8j7-KHxI
# Subst: x = sqrt(tan(x));

integrate(\(x) sqrt(sqrt(x^4 + 1) + 1), 0, 1)
integrate(\(x) (1 + x^4) / (1 - x^4)^2 * 2, 0, sqrt(tan(pi/8)))
integrate(\(x) 1/(1-x^2)^2 + 1/(1+x^2)^2, 0, sqrt(tan(pi/8)))
z = sqrt(tan(pi/8))
z/(1-z^4) + atan(z)/2 + log((1+z)/(1-z)) / 4
z/(1-z^4) + atan(z)/2 + atan(1i*z) / 2i


### I( sqrt(sqrt(1 + x^4) - 1) )
integrate(\(x) sqrt(sqrt(x^4 + 1) - 1), 0, 1)
integrate(\(x) sqrt((1-cos(x))/sin(x)) / cos(x)^2 / 2, 0, pi/4)
integrate(\(x) sqrt(sin(x/2)/cos(x/2)) / cos(x)^2 / 2, 0, pi/4)
integrate(\(x) 2 * x^2 * (1 + x^4)/(1-x^4)^2, 0, sqrt(tan(pi/8)))
integrate(\(x) x^2/(1-x^2)^2 + 1/(1+x^2) - 1/(1+x^2)^2, 0, sqrt(tan(pi/8)))
z = sqrt(tan(pi/8))
z^3/(1-z^4) + atan(z)/2 - atan(1i*z) / 2i


###
integrate(\(x) sqrt(sqrt(x^4 + 1) + x^2), 0, 1)
z = sqrt(tan(pi/8))
(z/(1-z^2) + atan(z)) / sqrt(2)

###
integrate(\(x) sqrt(sqrt(x^4 + 1) - x^2), 0, 1)
z = sqrt(tan(pi/8))
(z/(1+z^2) + atan(1i*z) / 1i) / sqrt(2)


##############

### Pow 2/3
integrate(\(x) (sqrt(x^6 + 1) + 1)^(2/3), 0, 1)
integrate(\(x) 1/3 * (cos(x/2)/sin(x/2))^(2/3) / cos(x)^2, 0, pi/4)
integrate(\(x) 2/3 * x^(-2/3) * (1 + 4*x^2/(1-x^2)^2) / (1+x^2), 0, tan(pi/8))
integrate(\(x) 2/3 * x^(-2/3) * (1+x^2) / (1-x^2)^2, 0, tan(pi/8))
integrate(\(x) 2 * (1+x^6) / (1-x^6)^2, 0, (tan(pi/8))^(1/3))
# TODO


###
integrate(\(x) (sqrt(x^6 + 1) - 1)^(2/3), 0, 1)
integrate(\(x) 1/3 * (sin(x/2)/cos(x/2))^(2/3) / cos(x)^2, 0, pi/4)
integrate(\(x) 2/3 * x^(2/3) * (1 + 4*x^2/(1-x^2)^2) / (1+x^2), 0, tan(pi/8))
integrate(\(x) 2/3 * x^(2/3) * (1+x^2) / (1-x^2)^2, 0, tan(pi/8))
integrate(\(x) 2 * x^4 * (1+x^6) / (1-x^6)^2, 0, (tan(pi/8))^(1/3))
# TODO

