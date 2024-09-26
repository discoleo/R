#################
##
## Trig: Formulas
## ATAN & ACOS


### Resources

# - Other Trig formulas are also found
#   in file: Integrals.Fractions.Unity.Definite.R;
# TODO: copy/move here?

#################

###
Im(acos(3/2 + 0i))
Re(log((sqrt(5) + 3) * 1i/2))
2 * log((1+sqrt(5))/2)


###
Im(acos(4/3 + 0i))
log((sqrt(7) + 4)/3)

###
Im(acos(5/3 + 0i))
log(3)
# varia:
(2*digamma(1/2) - digamma(1/2 - 1/3) - digamma(1/2 + 1/3)) / 3


###
Im(acos(2 + 0i))
log(sqrt(3) + 2)


###
Im(acos(3 + 0i))
2 * log(sqrt(2) + 1)

###
Im(acos(4 + 0i))
log((sqrt(15) + 4))

###
Im(acos(sqrt(2) + 0i))
log(sqrt(2) + 1)


### Gen:
n = sqrt(7)
Im(acos(n + 0i))
log((sqrt(n^2 - 1) + n))


##################
##################

###
n = 5
# n = 7 # ODD
#
id = seq(2, n-1, by = 2)
x = 2^(1/n); cs = cos(id*pi/n); sn = sin(id*pi/n);

sum(2 * sn * atan((1 - cs) / sn))
pi/2 / sin(pi/n)

