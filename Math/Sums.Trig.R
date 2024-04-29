

### Helper Constants

Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;


###################

###################
### Finite Sums ###

### ODD
n = 9
id = seq(1, floor(n/2));
sn = sin(2*pi*id/n);
#
sum(sn)
1/tan(pi/n/2) / 2
#
sum(id*sn)
n/sin(pi/n) / 4

### EVEN
n = 8 # n = 14
id = seq(1, n-1);
sn = sin(pi*id/n);
#
sum(sn)
1/tan(pi/n/2)
#
sum(id*sn)
# ???


################

###
n = 9; m = 4;
# n,m = Integers
id = seq(n-1)
sum(cos(2*pi*id*m/n) * digamma(id/n))
n*log(2*sin(pi*m/n)) + Euler


#####################
#####################

#####################
### Infinite Sums ###

### sum( cos(2*n*x) / n )
# Maths 505:  My take on this on wonderful infinite series from @drpeyam
# https://www.youtube.com/watch?v=mqPTvELJPM0

n = 20000
id = seq(n)

###
x = sqrt(2)
sum(cos(2*x*id)/id)
- log(2*abs(sin(x)))

###
x = sqrt(3)
sum((-1)^id * cos(2*x*id)/id)
- log(2*abs(cos(x)))

### Ex: H(3*n) - H(n)
# 1 + 1/2 - 2/3 + ...
x = pi/3
sum(cos(2*x*id)/id)
- log(2*abs(sin(x)))
#
id = seq(30000)
sum(c(1,1,-2)/id)
log(3)


########################
########################

### sum( cos(n) / n )
# Michael Penn: An infinite cosine sum.
# https://www.youtube.com/watch?v=R_Uf78si8jk
# Note: a special case of the previous sum;

id = seq(60000)
sum(cos(id) / id)
- log(2 - 2*cos(1)) / 2
- log(2*sin(1/2))


######################
######################

### prod( (1 - 1/n^2)^(+/- n) )
# 1. Michael Penn: a nice product from Ramanujan -- featuring 3 important constants!
#    https://www.youtube.com/watch?v=LG78xvZyRzU
# 2. Michael Penn: amazing techniques for uncovering this series.
#    https://www.youtube.com/watch?v=eZD7uHlSWfc


### 1 / cos(pi/2 * x)
x  = 1/7
id = seq(0, 30000)
1 / cos(pi/2 * x) # ==
4/pi * sum( (-1)^id * (2*id+1) / ((2*id+1)^2 - x^2) )


###
id = seq(1200)
prod( (1 - 1/(2*id+1)^2)^((-1)^id * (2*id+1)) )
pi/8 * exp(4*Catalan/pi)

