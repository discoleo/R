#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## ODEs: Polynomial Types
##
## draft v.0.2a

### Polynomials of Class 1


####################
####################

### Helper Functions

# include: DE.ODE.Helper.R;
source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")

#########################

### n = 5

### Example:
x = sqrt(2); # some Test
k = (x^2-x+1)^(1/5);
# Root:
# Root:
y = (x^2+x)*k^4 - (x^2+2*x+1)*k^3 + x^2*k^2 + (x^2+x)*k;
ye = expression((x^2+x)*k^4 - (x^2+2*x+1)*k^3 + x^2*k^2 + (x^2+x)*k)[[1]];
ye = do.call(substitute, list(ye, list(k=expression((x^2-x+1)^(1/5))[[1]])));
dy = eval(D(ye, "x"), list(x=x));
# P5 Polynomial (in y):
(1 + 7*x + 21*x^2 + 28*x^3 + 6*x^4 - 18*x^5 - 10*x^6 - 69*x^7 - 163*x^8 - 122*x^9 - 122*x^10 - 184*x^11 - 114*x^12 +
  - 50*x^13 - 64*x^14 - 27*x^15 + x^16 - x^17 - x^18) +
+ (5*x + 25*x^2 + 55*x^3 + 85*x^4 + 105*x^5 + 125*x^6 + 155*x^7 + 150*x^8 + 90*x^9 + 85*x^10 + 75*x^11 + 20*x^12 +
  + 10*x^13 + 5*x^14) * y +
- 5*x*(x+1)*(x^2 - x + 1)*(2*x^6 + 3*x^5 + 3*x^4 + 2*x + 1) * y^2 +
+ y^5;

# D(P): The Answer to the Universe and Everything
- 5*(20*x^9 + 27*x^8 + 24*x^7 + 14*x^6 + 18*x^5 + 25*x^4 + 4*x^3 + 4*x + 1) * y^3 +
+ 5*(14*x^13 + 26*x^12 + 48*x^11 + 165*x^10 + 170*x^9 + 162*x^8 + 240*x^7 + 217*x^6 + 150*x^5 +
	+ 105*x^4 + 68*x^3 + 33*x^2 + 10*x + 1) * y^2 +
- (18*x^17 + 17*x^16 - 16*x^15 + 405*x^14 + 896*x^13 + 650*x^12 + 1368*x^11 + 2024*x^10 + 1220*x^9 +
	+ 1098*x^8 + 1304*x^7 + 483*x^6 + 60*x^5 + 90*x^4 - 24*x^3 - 84*x^2 - 42*x - 7) * y +
+ 5*(x^2-x+1)*(x^16 + 2*x^15 + 25*x^13 + 89*x^12 + 114*x^11 + 139*x^10 + 209*x^9 + 192*x^8 +
	+ 105*x^7 + 76*x^6 + 40*x^5 - 26*x^4 - 48*x^3 - 28*x^2 - 8*x - 1) * dy +
- 20*x*(x+1)*(x^2-x+1)*(x^10 + 2*x^9 + 4*x^8 + 14*x^7 + 15*x^6 + 14*x^5 +
	+ 16*x^4 + 16*x^3 + 11*x^2 + 5*x + 1) * y * dy +
+ 15*x*(x+1)*(x^2 - x + 1)*(2*x^6 + 3*x^5 + 3*x^4 + 2*x + 1) * y^2 * dy # = 0


### Solution & Plot:
y = function(x, n=5) {
	# root
	k = rootn(x^2-x+1, n)
	y = (x^2+x)*k^4 - (x^2+2*x+1)*k^3 + x^2*k^2 + (x^2+x)*k
	y = sapply(y, round0)
	return(y)
}
y.free = function(x, y) {
	42*x*y + 50*x*y^2 - 20*x*y^3 + 84*x^2*y + 165*x^2*y^2 + 24*x^3*y + 340*x^3*y^2 - 20*x^3*y^3 - 90*x^4*y +
 + 525*x^4*y^2 - 125*x^4*y^3 - 60*x^5*y + 750*x^5*y^2 - 90*x^5*y^3 - 483*x^6*y + 1085*x^6*y^2 - 70*x^6*y^3 - 1304*x^7*y +
 + 1200*x^7*y^2 - 120*x^7*y^3 - 1098*x^8*y + 810*x^8*y^2 - 135*x^8*y^3 - 1220*x^9*y + 850*x^9*y^2 - 100*x^9*y^3 +
 - 2024*x^10*y + 825*x^10*y^2 - 1368*x^11*y + 240*x^11*y^2 - 650*x^12*y + 130*x^12*y^2 - 896*x^13*y +
 + 70*x^13*y^2 - 405*x^14*y + 16*x^15*y - 17*x^16*y - 18*x^17*y + 7*y + 5*y^2 - 5*y^3
}
dy = function(x) {
	y.x = y(x)
	div = y.x^2 * (15*x + 30*x^2 + 15*x^4 + 75*x^5 + 45*x^6 + 30*x^7 + 45*x^8 + 45*x^9 + 30*x^10) +
		y.x * (- 20*x - 100*x^2 - 220*x^3 - 340*x^4 - 420*x^5 - 500*x^6 - 620*x^7 - 600*x^8 - 360*x^9 - 340*x^10 - 300*x^11 +
			- 80*x^12 - 40*x^13 - 20*x^14) +
		(- 5 - 35*x - 105*x^2 - 140*x^3 - 30*x^4 + 90*x^5 + 50*x^6 + 345*x^7 + 815*x^8 + 610*x^9 + 610*x^10 +
			+ 920*x^11 + 570*x^12 + 250*x^13 + 320*x^14 + 135*x^15 - 5*x^16 + 5*x^17 + 5*x^18)
	dp =  -y.free(x, y.x)
	
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
curve(y, from=-2, to=2)
# a nice ?global/local? minimum
line.tan(c(-(4:1)/5, (1:6)/5), dx=3, p=y, dp=dy)


############

### Example: Simple n = 5
x = sqrt(2) # some Test
k = (x - 1)^(1/5)
# root
y = k^4 - k^3 + k^2 + k
### P5 Polynomial (in y):
(- 22 + 48*x - 30*x^2 + 5*x^3 - x^4) +
(15 - 20*x + 5*x^3)*y^1 +
(- 10 + 20*x - 10*x^2)*y^2 + y^5
### D(P5)
48*y - 60*x*y + 15*x^2*y - 4*x^3*y - 20*y^2 + 15*x^2*y^2 + 20*y^3 - 20*x*y^3 +
 + (110 - 240*x + 150*x^2 - 25*x^3 + 5*x^4) * dy +
 + (- 60 + 80*x - 20*x^3) * y * dy +
 + (30 - 60*x + 30*x^2) * y^2 * dy

### Solution & Plot:
y = function(x, n=5) {
	# root
	k = rootn(x - 1, n)
	y = k^4 - k^3 + k^2 + k
	y = sapply(y, round0)
	return(y)
}
y.free = function(x, y) {
	48*y - 60*x*y + 15*x^2*y - 4*x^3*y - 20*y^2 + 15*x^2*y^2 + 20*y^3 - 20*x*y^3
}
dy = function(x) {
	y.x = y(x)
	div = y.x^2 * (30 - 60*x + 30*x^2) +
		y.x * (- 60 + 80*x - 20*x^3) +
		(110 - 240*x + 150*x^2 - 25*x^3 + 5*x^4)
	dp =  -y.free(x, y.x)
	
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
### Plot:
curve(y, from=-2, to=3)
# Nice global minimum
line.tan(c(-(4:1)/3, (1:6)/3), dx=3, p=y, dp=dy)


############

### Example: Simple n = 5
a = 4 # Parameter
x = sqrt(2) # some Test
k = rootn(x - a, 5)
# root
y = k^4 - k^3 + k^2 + k
### P5 Polynomial (in y):
# y^5 + b2(a)*y^2 + b1(a)*y + b0(a)
(a + 42*a*x - 3*a*x^2 + 4*a*x^3 - 21*a^2 + 3*a^2*x - 6*a^2*x^2 - a^3 + 4*a^3*x - a^4 - x - 21*x^2 + x^3 - x^4) +
(5*a - 30*a*x - 15*a*x^2 + 15*a^2 + 15*a^2*x - 5*a^3 - 5*x + 15*x^2 + 5*x^3)*y +
(20*a*x - 10*a^2 - 10*x^2)*y^2 + y^5;

### ODE based on D(P5)
- 6*a*x*y - 30*a*x*y^2 + 12*a*x^2*y + 42*a*y - 30*a*y^2 + 20*a*y^3 - 12*a^2*x*y + 3*a^2*y + 15*a^2*y^2 +
	+ 4*a^3*y - 42*x*y + 30*x*y^2 - 20*x*y^3 + 3*x^2*y + 15*x^2*y^2 - 4*x^3*y - y - 5*y^2 +
 + (- 5*a - 210*a*x + 15*a*x^2 - 20*a*x^3 + 105*a^2 - 15*a^2*x + 30*a^2*x^2 + 5*a^3 - 20*a^3*x +
	+ 5*a^4 + 5*x + 105*x^2 - 5*x^3 + 5*x^4) * dy +
 + (- 20*a + 120*a*x + 60*a*x^2 - 60*a^2 - 60*a^2*x + 20*a^3 + 20*x - 60*x^2 - 20*x^3) * y * dy +
 + (- 60*a*x + 30*a^2 + 30*x^2) * y^2 * dy # = 0
 
### Solution & Plot:
y = function(x, a, n=5) {
	# root
	k = rootn(x - a, n)
	y = k^4 - k^3 + k^2 + k
	y = sapply(y, round0)
	return(y)
}
y.free = function(x, y, a) {
	- 6*a*x*y - 30*a*x*y^2 + 12*a*x^2*y + 42*a*y - 30*a*y^2 + 20*a*y^3 - 12*a^2*x*y +
	+ 3*a^2*y + 15*a^2*y^2 + 4*a^3*y - 42*x*y + 30*x*y^2 - 20*x*y^3 + 3*x^2*y + 15*x^2*y^2 - 4*x^3*y - y - 5*y^2
}
dy = function(x, a) {
	y.x = y(x, a)
	div = y.x^2 * (- 60*a*x + 30*a^2 + 30*x^2) +
		y.x * (- 20*a + 120*a*x + 60*a*x^2 - 60*a^2 - 60*a^2*x + 20*a^3 + 20*x - 60*x^2 - 20*x^3) +
		(- 5*a - 210*a*x + 15*a*x^2 - 20*a*x^3 + 105*a^2 - 15*a^2*x + 30*a^2*x^2 + 5*a^3 - 20*a^3*x +
		+ 5*a^4 + 5*x + 105*x^2 - 5*x^3 + 5*x^4)
	dp =  -y.free(x, y.x, a)
	
	dp = if(div != 0) dp / div else Inf;
	return(dp)
}
### Plot:
curve(y(x, a=a), from=a-2, to=a+3)
# a nice global minimum
line.tan(c(a-(4:1)/3, a+(1:6)/3), dx=3, p=y, dp=dy, a=a)

