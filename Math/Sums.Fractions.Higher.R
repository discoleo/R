

### Constants
Catalan = 0.915965594177219015054603514;
Euler   = 0.57721566490153286060651209008240243079;


##########################

### sum( 1 / (j^2 + 1)^p )

###
j = seq(0, 10000)
sum( 1 / (j^2 + 1)^2 )
pi^2/4 / cosh(pi)^2 + pi/4 / tanh(pi) + 1/2
# TODO


###
j = seq(10000)
sum( j / (j^2 + 1)^2 )
(pracma::psi(1, 1 + 1i) - pracma::psi(1, 1 - 1i)) * 1i/4
(pracma::psi(1, 1i) - pracma::psi(1, - 1i)) * 1i/4


###
k = sqrt(5)
j = seq(10000)
sum( j / (j^2 + k^2)^2 )
(pracma::psi(1, k*1i) - pracma::psi(1, - k*1i)) * 1i / (4*k)
# TODO: simplify ?

# Derivation:
# - generate Polygamma function;
k = sqrt(3)
j = seq(10000)
sum( j / (j^2 + k^2)^2 )
sum( (j + k*1i)^2 / (j^2 + k^2)^2 - (j - k*1i)^2 / (j^2 + k^2)^2) / 4i / k
sum( 1 / (j - k*1i)^2 - 1 / (j + k*1i)^2) / 4i / k
sum( 1/j^2 - 1 / (j + k*1i)^2 - (1/j^2 - 1 / (j - k*1i)^2) ) / 4i / k
(pracma::psi(1, k*1i) - pracma::psi(1, - k*1i)) * 1i/ (4*k)


### Psi Order 1:
k = sqrt(3)
(pracma::psi(1, k*1i) + pracma::psi(1, 1 - k*1i))
pi^2 / sin(pi*k*1i)^2
- 4*pi^2 / (exp(pi*k) - exp(-pi*k))^2


### "Sub-Sums"
# Kölbig, K. S. (1996). The polygamma function psi^k(x) for x=1/4 and x=3/4.
# J. Comput. Appl. Math. 75 (1): 43–46.
# doi:10.1016/S0377-0427(96)00055-6

(pracma::psi(1, 1/4) - pracma::psi(1, 3/4))
16 * Catalan
#
pracma::psi(1, 1/4)
pi^2 + 8*Catalan
#
pracma::psi(1, 3/4)
pi^2 - 8*Catalan


### Psi Order 2:

### Psi(2, 1/3)
# see file: Integrals.AnalyticCont.R;
pracma::psi(2, 1/3)
- 26 * pracma::zeta(3) - pi^3 * cos(pi/3) / sin(pi/3)^3
#
pracma::psi(2, 2/3)
- 26 * pracma::zeta(3) + pi^3 * cos(pi/3) / sin(pi/3)^3

### Psi(2, 1/4)
# see either article by Kolbig,
# or see file: Integrals.AnalyticCont.R;
pracma::psi(2, 1/4)
- 56 * pracma::zeta(3) - pi^3 * cos(pi/4) / sin(pi/4)^3
#
pracma::psi(2, 3/4)
- 56 * pracma::zeta(3) + pi^3 * cos(pi/4) / sin(pi/4)^3

### Psi(2, 1/6)
# see file: Integrals.AnalyticCont.R;
pracma::psi(2, 1/6)
- 182 * pracma::zeta(3) - pi^3 * cos(pi/6) / sin(pi/6)^3
#
pracma::psi(2, 5/6)
- 182 * pracma::zeta(3) + pi^3 * cos(pi/6) / sin(pi/6)^3


### Pow = 3
### "Sub-Sums"
id = 10000

### 3*k + 1
sum(1/(3*seq(0, id) + 1)^3)
13/27 * pracma::zeta(3) + pi^3/54 * cos(pi/3) / sin(pi/3)^3
#
sum(1/(3*seq(0, id) + 2)^3)
13/27 * pracma::zeta(3) - pi^3/54 * cos(pi/3) / sin(pi/3)^3


### 4*k + 1
sum(1/(4*seq(0, id) + 1)^3)
28/64 * pracma::zeta(3) + pi^3/128 * cos(pi/4) / sin(pi/4)^3
#
sum(1/(4*seq(0, id) + 3)^3)
28/64 * pracma::zeta(3) - pi^3/128 * cos(pi/4) / sin(pi/4)^3


### 6*k + 1
sum(1/(6*seq(0, id) + 1)^3)
91/6^3 * pracma::zeta(3) + pi^3/(2*6^3) * cos(pi/6) / sin(pi/6)^3
#
sum(1/(6*seq(0, id) + 5)^3)
91/6^3 * pracma::zeta(3) - pi^3/(2*6^3) * cos(pi/6) / sin(pi/6)^3


### 12*k + 1
id = seq(20000, 0)
sum(1/(12*id + 1)^3)
# TODO: ?


# Derivation:
# - but indeterminate!
sum(1/(12*id + 1)^3, 1/(12*id + 7)^3)
91/6^3 * pracma::zeta(3) + pi^3/(2*6^3) * cos(pi/6) / sin(pi/6)^3
#
sum(1/(12*id + 5)^3, 1/(12*id + 11)^3)
91/6^3 * pracma::zeta(3) - pi^3/(2*6^3) * cos(pi/6) / sin(pi/6)^3
#
sum(1/(12*id + 1)^3, - 1/(12*id + 11)^3)
(pi^3 * cos(pi/12) / sin(pi/12)^3 + 12^3 - (12/11)^3)/2/12^3
#
sum(1/(12*id + 5)^3, - 1/(12*id + 7)^3)
(pi^3 * cos(5*pi/12) / sin(5*pi/12)^3 + (12/5)^3 - (12/7)^3)/2/12^3


##########################
##########################

### sum( (-1)^j / (2*j + 1)^p )
# Flammable Maths: An AMAZING Journey of Series Evaluation!
# Calculating Euler's Sum! [ Series pi^3/32 (-1)^k/(2k+1)^3 ]
# https://www.youtube.com/watch?v=LQi7BMWU_1Y
# - could benefit from streamlining the presentation
#   and correcting the mistakes (odd function);

# sum( (1 - (-1)^k) * sin(k*pi*x) ) =
# (x - x^2) * pi^3 / 4;

###
id = seq(0, 400000)
sum((-1)^id / (2*id+1)^1)
pi/4

###
id = seq(0, 40000)
sum((-1)^id / (2*id+1)^3)
pi^3/32

###
sum((-1)^id / (2*id+1)^5)
5/1536 * pi^5


###
# x = 1/4 =>
# - see also previous section for an alternative approach based on psi(2, 1/4);
sum((-1)^id / (4*id+1)^3) + sum((-1)^id / (4*id+3)^3)
pi^3 * 3*sqrt(2)/(32*4)

# =>
sum( 1/(8*id+1)^3, -1/(8*id+7)^3 )
pi^3 * (1 + 3*sqrt(2)/4)/(32*2)
# =>
sum( 1/(8*id+1)^3, 1/(8*id+3)^3)
28/64 * pracma::zeta(3) - pi^3/128 * cos(pi/4) / sin(pi/4)^3 +
	+ pi^3 * (1 + 3*sqrt(2)/4)/(32*2)

# x = 1/8
# but redundant, as == f(psi(2, k/16));
id =  seq(0, 40000)
sum((-1)^id * sin(pi/8) / (8*id+1)^3) +
	+ sum((-1)^id * sin(3*pi/8) / (8*id+3)^3) +
	+ sum((-1)^id * sin(3*pi/8) / (8*id+5)^3) +
	+ sum((-1)^id * sin(pi/8) / (8*id+7)^3)
pi^3 * 7 / (64*8)
# x = 3/8
sum((-1)^id * sin(pi*3/8) / (8*id+1)^3) +
	- sum((-1)^id * sin(pi/8) / (8*id+3)^3) +
	- sum((-1)^id * sin(pi/8) / (8*id+5)^3) +
	+ sum((-1)^id * sin(3*pi/8) / (8*id+7)^3)
pi^3 * 15 / (64*8)

# Note: the formulas above can be composed using all versions of:
k = 1
(pracma::psi(2, k/16) - pracma::psi(2, 1 - k/16))/2
- pi^3 * cos(k*pi/16) / sin(k*pi/16)^3
#
(pracma::zeta(3) - pracma::psi(2, k/16) / 2 - k^3) / 16^3
sum(1/(16*id + k)^3)

# TODO: find solution;


######################
######################

### Sum( 1 / (n^3 + 1) )
m = cos(pi/3) + 1i*sin(pi/3);

id = seq(0, 10000)
sum(1/(id^3 + 1))
(pracma::psi(-m)*m + pracma::psi(-1/m)/m + Euler) / 3
(pracma::psi(-m)*m - pracma::psi(-1/m)*m^2 + Euler) / 3


### Gen: Sum( 1 / (n^3 + k^3) )
k = 3^(1/3)
id = seq(0, 10000)
sum(1/(id^3 + k^3))
(pracma::psi(-k*m)*m - pracma::psi(-k/m)*m^2 - pracma::psi(k)) / (3*k^2)


### Sum( 1 / (n^5 + 1) )
m = cos(pi/5) + 1i*sin(pi/5);

id = seq(0, 10000)
sum(1/(id^5 + 1))
(pracma::psi(-m)*m - pracma::psi(-1/m)*m^4 +
	+ pracma::psi(-m^3)*m^3 - pracma::psi(-1/m^3)*m^2 + Euler) / 5


###########################
###########################

### Dilog/Polylog Functions

### Sum( z^n / n^2 )
# Michael Penn: The dilogarithm -- a favorite "special function"
# https://www.youtube.com/watch?v=Bf69XqbpnzY

id = seq(40000)
sum( 1/(2^id * id^2) )
pi^2/12 - log(2)^2/2

