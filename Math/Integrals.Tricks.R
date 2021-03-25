########################
###
### Leonard Mada
### [the one and only]
###
### Integral Tricks
###
### draft v.0.1a


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

integrate(function(x) sapply(x,
	function(x) (x+a)*log(x)/sum(b*x^(3:0))), lower=0, upper=Inf)

integrate(function(x) sapply(x,
	function(x) (x+1)*log(a)/a/sum(bf*x^(3:0))), lower=0, upper=Inf)


#########
### Ex 2:
integrate(function(x) (x^3+8)*log(x)/(x^5 + 3*x^4 + 24*x + 32), lower=0, upper=Inf)
integrate(function(x) (x^3+1)*log(2)/2/(x^5 + 3/2*x^4 + 3/2*x + 1), lower=0, upper=Inf)

