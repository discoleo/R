

### Systems: Lambert


library(pracma)


### x^x = k
k = 2
x = exp(lambertWp(log(k)))
#
x^x


### x^n * log(x/b) = 1
n = 3
b = log(2)
#
x = b * exp(lambertWp(n/b^n) / n)
x^n * log(x/b) # = 1


###
x = exp(lambertWp(exp(-1)) / 2 + 1/2)
# Maximum of function:
log(x) / (x^2 + 1)

