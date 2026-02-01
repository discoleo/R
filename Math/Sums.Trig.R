

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


##########################
##########################

### Sum( Coth(n*pi) / n^3 )
# 1. Michael Penn: a nice Ramanujan sum
#    https://www.youtube.com/watch?v=c3B3ILEtrbk
# 2. Ce Xu. Some infinite series involving hyperbolic functions.
#    https://arxiv.org/pdf/1707.06673


### Pow = 3
id = 1:10000;
sum( 1 / tanh(id*pi) / id^3 )
7/180 * pi^3;


### Sum( Coth(n*pi) / n^7 )
id = 1:10000;
sum( 1 / tanh(id*pi) / id^7 )
(pracma::zeta(8) - pracma::zeta(4)^2 + 2*pracma::zeta(2)*pracma::zeta(6)) * 1/pi;


# Sum( Tanh(n*pi) / n^3 )
sum( tanh(id*pi) / id^3 )
# TODO


### Sum( Tanh((n-1/2)*pi) / n^3 )
sum( tanh((id-1/2)*pi) / (2*id-1)^3 )
pi^3 / 32;


### Sum( coth((2*n+1)*pi) / n^3 )
id = 1:10000;
sum( 1 / tanh((2*id+1)*pi) / id^3 )
# Check: ???
pracma::zeta(3); # almost?

### Sum( coth(2*n*pi) / n^3 )
id = 1:40000;
sum( 1 / tanh(2*id*pi) / id^3 )
# TODO


# Helper:
id = 1:2000; # Converges slowly;
z = sqrt(2);
1/tanh(pi*z);
1/pi * (1/z + 2 * sum(z / (z^2 + id^2)));


# Sum( 1 / (m * n * (m+n)) )
id = 1:5000;
sum(sapply(id, \(id1) { sum(1 / (id*id1*(id+id1 + 0.0))); } ))
2*pracma::zeta(3);

# Sum( 1 / (m^2 * n * (m+n)) )
ids = 1:20000;
id1 = 1:5000; id2 = id1^2;
sum(sapply(ids, \(ids) { sum(1 / (id2*ids*(id1+ids))); } ))
pi^4 / 72;

# Sum( 1 / (m^2 * n^2 * (m+n)) )
id  = 1:2000;
id2 = id^2;
sum(sapply(id, \(id1) { 1 / (id1^2*id2*(id1+id)); } ))
integrate(\(x) sapply(x, \(x) polylog2(x)^2 / x), 0, 1, rel.tol=1E-12)
2*pracma::zeta(2)*pracma::zeta(3) - 3*pracma::zeta(5);


# Sum( 1 / (m^2 * n^2 * (m+n)^2) )
id  = 1:2000;
id2 = id^2;
sum(sapply(id, \(id1) { 1 / (id1^2*id2*(id1+id)^2); } ))
pracma::zeta(6) / 3;

# Sum( 1 / (m^2 * n^2 * (m^2 + n^2)) )
id1 = (1:4000)^2;
id2 = id1;
zlp = sum(sapply(id1, \(id1) { 1 / (id1*id2*(id1+id2)); } ))
0.6543415088859439; # higher precision;
# TODO

