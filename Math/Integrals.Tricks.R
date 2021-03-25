########################
###
### Leonard Mada
### [the one and only]
###
### Integral Tricks
###
### draft v.0.1b


### various Integral Tricks



################

### Example 1:
# - based on: "This trick is new to me!"
#   https://www.youtube.com/watch?v=JCuplIQ6JG4

a  = 3
b2 = 5
#
b  = c(1, b2, b2*a, a^3)
bf = c(1, b2/a, b2/a, 1)
# Note: polynomial fraction can be simplified;
integrate(function(x) sapply(x,
	function(x) (x+a)*log(x)/sum(b*x^(3:0))), lower=0, upper=Inf)

integrate(function(x) sapply(x,
	function(x) log(a)/a * (x+1)/sum(bf*x^(3:0))), lower=0, upper=Inf)


#########
### Ex 2:
a = 3
b1n = 4
b1 = 2
# b2 can be added as well;

integrate(function(x) (x^2 + b1n*x + a^2)*log(x)/
	(x^4 + b1*x^3 + b1*a^2*x + a^4), lower=0, upper=Inf)
integrate(function(x) log(a)/a * (x^2 + b1n/a*x + 1)/
	(x^4 + b1/a*x^3 + b1/a*x + 1), lower=0, upper=Inf)


#########
### Ex 3:
integrate(function(x) (x^3+8)*log(x)/(x^5 + 3*x^4 + 24*x + 32), lower=0, upper=Inf)
integrate(function(x) log(2)/2 * (x^3+1)/(x^5 + 3/2*x^4 + 3/2*x + 1), lower=0, upper=Inf)


