########################
###
### Leonard Mada
### [the one and only]


#################

### N-th Power
# based on:
# Michael Penn: A nice limit
# https://www.youtube.com/watch?v=RoqErc0NKmE


limit.power = function(x, n=25) {
	(sum(x^(1/n)) / length(x))^n
}
mean.geom = function(x) {
	len = length(x)
	prod(x)^(1/len)
}

### Test

x = c(2,3,5)
limit.power(x)
mean.geom(x)

### Ex 2:
x = c(2,3,5,7)
limit.power(x)
mean.geom(x)

### Ex 3:
x = c(2,4,4,7)
limit.power(x, n=40)
mean.geom(x)
