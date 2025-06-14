########################
##
## Leonard Mada
## [the one and only]
##
## Integrals: Logarithms
## Fractions: Polylog / Li2
##
## draft v.0.1a


### Polylog / Li2
# - Various examples;


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;

### Helper Functions

source("Integrals.Polylog.Helper.R")

####################

### Log( DIV )

### I( log((1-x)/(x+1)) / (1 - 4*x^2) )
# on [0, Inf]
integrate(\(x) log((1-x)/(x+1)) / (1 - 4*x^2) + log(3)/(2-4*x), 0, 1/2)$value +
integrate(\(x) log((1-x)/(x+1)) / (1 - 4*x^2) + log(3)/(2-4*x), 1/2, 1)$value +
integrate(\(x) log((x-1)/(x+1)) / (1 - 4*x^2), 1, Inf)$value
pi^2/24 - pracma::polylog(-1/2, 2)/2 - (log(3/2)^2 - log(3)^2)/4;

### on [1, Inf]
integrate(\(x) log((x-1)/(x+1)) / (1 - 4*x^2), 1, Inf)
pi^2/24 - (2*log(3/2)^2 - log(3)^2) / 8 +
- (pracma::polylog(-1/2, 2) + pracma::polylog(1/3, 2) - pracma::polylog(-1/3, 2)) / 2;

### on [0, 1]
integrate(\(x) log((1-x)/(x+1)) / (1 - 4*x^2) + log(3)/(2-4*x), 0, 1/2)$value +
integrate(\(x) log((1-x)/(x+1)) / (1 - 4*x^2) + log(3)/(2-4*x), 1/2, 1)$value
log(3)^2 / 8 + (pracma::polylog(1/3, 2) - pracma::polylog(-1/3, 2)) / 2;
pi^2/12 - (pracma::polylog(1/3, 2) + log(3)^2 / 4) / 2;

### on [0, 1/2]
integrate(\(x) log((1-x)/(x+1)) / (1 - 4*x^2) + log(3)/(2-4*x), 0, 1/2)
integrate(\(x) log((x-1)/(x+1)) / (1 - 4*x^2), 1, Inf)$value - log(2) * log(3) / 4;
pi^2/24 - (2*log(3/2)^2 - log(3)^2 + 2*log(2) * log(3)) / 8 +
- (pracma::polylog(-1/2, 2) + pracma::polylog(1/3, 2) - pracma::polylog(-1/3, 2)) / 2;

# Helper:
integrate(\(x) log((1-x)) / (x^2 - 4), 0, 1)
pi^2/48 + pracma::polylog(1/3, 2)/4;

#
integrate(\(x) log((1+x)) / (x^2 - 4), 0, 1)
-3/4 * pracma::polylog(1/3, 2) + pi^2/16 - log(3)^2 * 3/8;

# Li2(-1/3)
pracma::polylog(-1/3, 2)
2*pracma::polylog(1/3, 2) - pi^2/6 + log(3)^2 / 2;

