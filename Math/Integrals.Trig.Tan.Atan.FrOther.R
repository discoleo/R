##################
##
## Leonard Mada
##
## Integrals: ATAN
## Fractions: Other Types

### I( atan(sqrt(x^2 + b^2)) / ((x^2 + 1) * sqrt(x^2 + b^2)) )
# Maths 505: Ahmed's monster integral: solution using Feynman's technique
# https://www.youtube.com/watch?v=Rq-bYzyoSwU


###
b = sqrt(5)
integrate(\(x) (2*x^2 + b^2) * atan(sqrt(x^2 + b^2)) /
	((x^2 + 1) * (x^2 + b^2-1) * sqrt(x^2 + b^2)), 0, 1)
pi/4 * (3/2*pi + atan(sqrt(b^2 - 1)) - 2*atan(sqrt((b^2+1)/(b^2-1))) +
	- 2*atan(sqrt((b^2+1)*(b^2-1)))) / sqrt(b^2 - 1)

# Note: the I() collapses for b^2 - 1 == 1 to:
integrate(\(x) atan(sqrt(x^2 + 2)) / ((x^2 + 1) * sqrt(x^2 + 2)), 0, 1)
pi^2/8 * (7/4 - 4/3)


###
b = sqrt(5)
integrate(\(x) atan(sqrt(x^2 + b^2)) / ((x^2 + 1) * sqrt(x^2 + b^2)), 0, 1)
pi/4 * (3/2*pi + atan(sqrt(b^2 - 1)) - 2*atan(sqrt((b^2+1)/(b^2-1))) +
	- 2*atan(sqrt((b^2+1)*(b^2-1)))) / sqrt(b^2 - 1) +
	- integrate(\(x) atan(sqrt(x^2 + b^2)) / ((x^2 + b^2-1) * sqrt(x^2 + b^2)), 0, 1)$value

### [redundant]
integrate(\(x) atan(sqrt(x^2 + b^2)) / ((x^2 + b^2-1) * sqrt(x^2 + b^2)), 0, 1)
pi/4*(3/2*pi + atan(sqrt(b^2-1)) - 2*atan(sqrt((b^2+1)*(b^2-1))) +
	- 2*atan(sqrt((b^2+1)/(b^2-1)))) / sqrt(b^2 - 1) +
	- integrate(\(x) atan(sqrt(x^2 + b^2)) / ((x^2 + 1) * sqrt(x^2 + b^2)), 0, 1)$value


# D(I): 1st I
# I( atan(t * sqrt(x^2 + b^2)) / ((x^2 + 1) * sqrt(x^2 + b^2)) ) on [0,1]
pi/4 / ((b^2-1)*t^2 + 1) +
	- t * atan(t / sqrt(b^2*t^2 + 1)) / sqrt(b^2*t^2 + 1) / ((b^2-1)*t^2 + 1)
# t => 1/t; dt => - dt / t^2;
# atan(1 / sqrt(t^2 + b^2)) / sqrt(t^2 + b^2) / (t^2 + (b^2-1))
# (pi/2 - atan(sqrt(t^2 + b^2))) / sqrt(t^2 + b^2) / (t^2 + (b^2-1))

# D(I): 2nd I
# I( atan(t * sqrt(x^2 + b^2)) / ((x^2 + b^2-1) * sqrt(x^2 + b^2)) ) on [0,1]
# I: 1 / ((x^2 + b^2-1) * (t^2*x^2 + t^2*b^2 + 1))
# ( 1 / (x^2 + b^2-1) - t^2 / (t^2*x^2 + t^2*b^2 + 1) ) / (t^2 + 1)
(pi/2 - atan(sqrt(b^2-1))) / ((t^2 + 1)*sqrt(b^2 - 1)) +
	- t * atan(t / sqrt(b^2*t^2 + 1)) / sqrt(b^2*t^2 + 1) / (t^2 + 1)


### Helper:
b = 1/sqrt(5)
b = sqrt(5)
integrate(\(x) 1 / ((x^2+1) * sqrt(x^2 + b^2)), 0, 1)
integrate(\(x) 1 / (x^2 + b^2 - 1), sqrt(b^2+1), Inf)
# b^2 > 1
pi/2/sqrt(b^2-1) - atan(sqrt((b^2+1)/(b^2-1)))/sqrt(b^2-1)
# b^2 < 1
(log(sqrt(b^2+1) + sqrt(1 - b^2)) - log(b) - log(2)/2) / sqrt(1 - b^2)

#
integrate(\(x) 1 / ((x^2 + b^2-1) * sqrt(x^2 + b^2)), 0, 1)
integrate(\(x) x / (((b^2-1)*x^2 + 1) * sqrt(b^2*x^2 + 1)), 1, Inf)
integrate(\(x) 1 / ((b^2-1)*x^2 + 1), sqrt(b^2+1), Inf)
# b^2 > 1
pi/2/sqrt(b^2-1) - atan(sqrt((b^2+1)*(b^2-1)))/sqrt(b^2-1)


###
b = sqrt(5)
integrate(\(x) atan(sqrt(x^2 + b^2)) / ((x^2 + 1) * sqrt(x^2 + b^2)), 0, 1)

# D(b)
# D( atan(sqrt(x^2 + b^2)) * sqrt(x^2 + b^2) / (x^2 + 1) )
# b / ((x^2 + 1) * (x^2 + b^2+1)) + b * atan(sqrt(x^2 + b^2)) / sqrt(x^2 + b^2) / (x^2 + 1)
# ???

### TODO:
integrate(\(x) atan(sqrt(x^2 + b^2)) / ((x^2 + 2) * sqrt(x^2 + b^2)), 0, 1)
# much work;

#################

### I( atan(sqrt(x^2 + b^2)) / (x^2 + b^2)^(3/2) )
lim = c(1/5, 1/3)
integrate(\(x) atan(sqrt(x^2 + 1)) / (x^2 + 1)^(3/2), lim[1], lim[2])
x = lim
diff(x * atan(sqrt(x^2 + 1)) / sqrt(x^2 + 1) + atan(x) - sqrt(2) * atan(x/sqrt(2)))


### Gen:
lim = c(1/5, 1/3)
b = sqrt(5)
integrate(\(x) atan(sqrt(x^2 + b^2)) / (x^2 + b^2)^(3/2), lim[1], lim[2])
x = lim
diff(x * atan(sqrt(x^2 + b^2)) / sqrt(x^2 + b^2) + b*atan(x/b) +
	- sqrt(b^2 + 1) * atan(x/sqrt(b^2 + 1))) / b^2

