

########################

### I( log(x) / (x + 1)^p ) on [0, 1]

# TODO: truncated Polylog function;

###
integrate(function(x) log(x)/(x+1)^2, 0, 1)
- log(2)

###
integrate(function(x) log(x)/(x+1)^3, 0, 1)
- log(2)/2 + (1/2 - 1)/2

###
integrate(function(x) log(x)/(x+1)^4, 0, 1)
- log(2)/3 + (1/2 + 1/4/2 - 1 - 1/2)/3

###
integrate(function(x) log(x)/(x+1)^5, 0, 1)
- log(2)/4 + (1/2 + 1/4/2 + 1/8/3 - 1 - 1/2 - 1/3)/4

###
integrate(function(x) log(x)/(x+1)^6, 0, 1)
- log(2)/5 + (1/2 + 1/4/2 + 1/8/3 + 1/16/4 - 1 - 1/2 - 1/3 - 1/4)/5

