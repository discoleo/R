

constEuler = 0.57721566490153286060651209008240243079;

### I( (1 - exp(-x) - exp(-1/x)) / x )
# Michael Penn: if you gamma be my constant, you gamma solve this integral.
# https://www.youtube.com/watch?v=heIrFaP05Lk

integrate(\(x) (1 - exp(-x) - exp(-1/x)) / x, 0, 1)
constEuler


#########################
#########################

### I( cos(x) / (exp(1/x) + 1) )
# Michael Penn: is this integration trick TOO POWERFUL?
# https://www.youtube.com/watch?v=PX2QXILRgsc

###
lim = 1
integrate(\(x) cos(x) / (exp(1/x) + 1), -lim, lim)
sin(lim)

### Various Variants:
lim = sqrt(pi)
integrate(\(x) 1 / (exp(1/x) + 1), -lim, lim)
lim

###
lim = 1
integrate(\(x) cos(x)^2 / (exp(1/x) + 1), -lim, lim)
sin(2*lim) / 4 + lim/2

###
lim = sqrt(pi)
integrate(\(x) sin(x)^2 / (exp(1/x) + 1), -lim, lim)
lim/2 - sin(2*lim) / 4


### Variants: exp() - 1

###
lim = sqrt(pi)
integrate(\(x) 1 / (exp(1/x) - 1), -lim, lim)
- lim

###
lim = sqrt(pi)
integrate(\(x) cos(x) / (exp(1/x) - 1), -lim, lim)
- sin(lim)


