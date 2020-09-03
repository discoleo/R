
########################
###
### Leonard Mada
### [the one and only]
###
### Exact Integration of Polynomial Fractions
### - Roots of unity: Higher Powers
###   Integral( 1 / (x^n - 1)^p )dx
### - Roots of minus unity:
###   Integral( 1 / (x^n + 1)^p )dx
### - Polynomial fractions:
###   Integral( P(x) / (x^n - 1)^p )dx
###
### draft v.0.1a


### Roots of unity: Baseline
#   Integral( P(x) / (x^n - 1) )dx
#   Integral( P(x) / (x^n + 1) )dx
# - see file: Integrals.Fractions.Unity.R


### Terminology

# F(k, n, p) = x^k / (x^n - 1)^p; # only the Fraction!
# I(k, n, p) = Integral( x^k / (x^n - 1)^p ) dx;


#######################


################
### Examples ###
################

### I(k, n, 1)
# - is computed in Integrals.Fractions.Unity.R;


##################
### I(n, n, p + 1)

# (x^n - 1 + 1) / (x^n - 1)^(p+1)
# = 1/(x^n - 1)^p + 1/(x^n - 1)^(p+1)
# =>
# I(n, n, p+1) = I(0, n, p) + I(0, n, p+1)
# I(n, n, p+1) - I(0, n, p+1) = I(0, n, p)

# d/dx F(1, n, p)
# = (x^n - 1 - n*p*x^n) / (x^n - 1 )^(p+1)
# = -(n*p - 1)*F(n, n, p+1) - F(0, n, p+1)
# =>
# (n*p - 1)*I(n, n, p+1) + I(0, n, p+1) = F(1, n, p);
# Note: Integral = Fraction!

### I(n, n, p+1)
# I(n, n, p+1) = 1/(n*p) * (F(1, n, p) + I(0, n, p))

### I(0, n, p+1)
# I(0, n, p+1) = 1/(n*p) * (F(1, n, p) - (n*p - 1)*I(0, n, p))


########################
