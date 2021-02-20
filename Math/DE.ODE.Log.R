
########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Logarithms
###
### draft v.0.1c-fix


### ODEs Derived from Logarithms

### Note:
# Lambert W examples:
# - should be in a different file;


#########################

### Helper functions

library(pracma)
# - may be needed to solve various equations;

# include: DE.ODE.Helper.R;
source("DE.ODE.Helper.R")


#########################
#########################


### Section A: Linear ODEs

### Derived from:
### y = F1(x)*log(P1(x)) + F2(x)*log(P2(x))

### Examples:
### y = (x^m + c)*log(x^n + a) + (x^m - c)*log(x^n + b)
# - simplifies massively & eliminates y from D2(y);
# - does NOT simplify when: F1 - F2 != constant;

### D(y)
# [not run]
m*x^(m-1)*log(x^n + a) + m*x^(m-1)*log(x^n + b) +
	+ n*x^(n-1)*(x^m + c) / (x^n + a) +
	+ n*x^(n-1)*(x^m - c) / (x^n + b)
m*x^(m-1)*log(x^n + a) + m*x^(m-1)*log(x^n + b) +
	+ n*x^(n-1)*(2*x^(m+n) + (a+b)*x^m - c*(a-b)) / ((x^n + a) * (x^n + b))
# let: T(x) = x^(n-1)*(2*x^(m+n) + (a+b)*x^m - c*(a-b)) / ((x^n + a) * (x^n + b))
# Note: + n * T(x);
m*x^(m-1)*log(x^n + a) + m*x^(m-1)*log(x^n + b) + n*T
### Solve linearly for log(P(x))
log(x^n + a) =
	((x^m - c)*dy - m*x^(m-1)*y - n*T*(x^m - c)) / - (2*c*m*x^(m-1));
log(x^n + b) =
	((x^m + c)*dy - m*x^(m-1)*y - n*T*(x^m + c)) /   (2*c*m*x^(m-1));
log(x^n + a) + log(x^n + b) = (dy - n*T) / (m*x^(m-1));

###
(x^n + a)*(x^n + b) * dy +
	- m*x^(m-1)*(x^n + a)*(x^n + b)*(log(x^n + a) + log(x^n + b)) +
	- n*x^(n-1)*(2*x^(m+n) + (a+b)*x^m - c*(a-b));
T(x) = x^(n-1)*(2*x^(m+n) + (a+b)*x^m - c*(a-b)) / ((x^n + a) * (x^n + b))

### D2(y)
(x^n + a)*(x^n + b) * d2y.x + n*x^(n-1)*(2*x^n + a + b) * dy.x +
	- (m*(m-1)*x^(m-2)*(x^n + a)*(x^n + b) + m*n*x^(m+n-2)*(2*x^n + a + b)) *
		(log(x^n + a) + log(x^n + b)) +
	- n*x^(m+n-2)*(2*m*x^n + m*(a + b)) +
	- n*x^(m+n-2)*(2*(m+n)*x^n + m*(a+b)) +
	- n*(n-1)*x^(n-2)*(2*x^(m+n) + (a+b)*x^m - c*(a-b))
(x^n + a)*(x^n + b) * d2y.x + n*x^(n-1)*(2*x^n + a + b) * dy.x +
	- (m*(m-1)*x^(m-2)*(x^n + a)*(x^n + b) + m*n*x^(m+n-2)*(2*x^n + a + b)) *
		(log(x^n + a) + log(x^n + b)) +
	- n*x^(m+n-2)*(2*(2*m+2*n-1)*x^n + (2*m+n-1)*(a+b)) +
	+ n*(n-1)*c*(a-b)*x^(n-2)
x*(x^n + a)*(x^n + b) * d2y + n*x^n*(2*x^n + a + b) * dy +
	- ((m-1)*(x^n + a)*(x^n + b) + n*x^n*(2*x^n + a + b)) *
		(dy - n*T) +
	- n*x^(m+n-1)*(2*(2*m+2*n-1)*x^n + (2*m+n-1)*(a+b)) +
	+ n*(n-1)*c*(a-b)*x^(n-1)
x*(x^n + a)*(x^n + b) * d2y +
	- (m-1)*x^n*(x^n + (a+b)) * dy  - (m-1)*a*b*dy + # - (m-1)*(x^n+a)*(x^n+b)*dy
	+ n*((m-1)*(x^n + a)*(x^n + b) + n*x^n*(2*x^n + a + b)) * T +
	- n*x^(m+n-1)*(2*(2*m+2*n-1)*x^n + (2*m+n-1)*(a+b)) +
	+ n*(n-1)*c*(a-b)*x^(n-1)
x*(x^n + a)*(x^n + b) * d2y +
	- (m-1)*(x^n + a)*(x^n + b)*dy +
	- n*x^(m+n-1)*(2*(m+2*n)*x^n + (m+n)*(a+b)) +
	+ n*(n-m)*c*(a-b)*x^(n-1) +
	+ n^2*x^n*(2*x^n + a + b) * T
# D(D(y) / x^(m-1)) = ... * x^(m-1) / (x*(x^n + a)*(x^n + b));


### Example 1:
m = 2
# TODO

### m = 2; a+b = 0;
# TODO


### Solution & Plot:
y = function(x, a, b=-a, ct=0, n=2, m=2) {
	xm = x^m; xn = x^n;
	x.prod = (xn + a)*(xn +b); x.div = (xn + a) / (xn + b);
	val = xm*log(x.prod) + ct*log(x.div);
	return(val)
}
dy = function(x, a, b=-a, ct=0, n=2, m=2) {
	xm = x^m; xn = x^n;
	# y.x = y(x, a=a, b=b, ct=ct, n=n, m=m);
	x.prod = (xn + a)*(xn + b);
	div = x;
	Tnx = n*xn*(2*xm*xn + (a+b)*xm - ct*(a-b)) / (x.prod);
	dp = m*xm*log(x.prod) + Tnx;
	dp = ifelse(div != 0, dp/div, 1); # TODO
	return(dp)
}
d2y = function(x, a, b=-a, ct=0, n=2, m=2) {
	xm = x^m; xn = x^n;
	x.prod = (xn + a)*(xn + b);
	### T
	div = x * x.prod
	T = xn*(2*xm*xn + (a+b)*xm - ct*(a-b)) / div;
	### d2y
	# y.x  =  y(x, a=a, b=b, ct=ct, n=n, m=m);
	dy.x = dy(x, a=a, b=b, ct=ct, n=n, m=m);
	d2p =
	(m-1) * x.prod * dy.x +
		+ n*(xm*(2*(m+2*n)*xn + (m+n)*(a+b)) - (n-m)*ct*(a-b)) * xn / x +
		- n^2*xn*(2*xn + a + b) * T;
	d2p / div;
}
### Plot:
n = 2; m = 2;
a = 1; b = -a;
ct = 2;
px = 6/7 + (1:4)*2/7
curve(y(x, a=a, b=b, ct=ct, n=n, m=m), from= 1, to = 2.5, ylim=c(-2,20))
sapply(px, line.tan, dx=3, p=y, dp=dy, a=a, b=b, ct=ct, n=n, m=m)
# x*log(x)-like:
curve(dy(x, a=a, b=b, ct=ct, n=n, m=m), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, a=a, b=b, ct=ct, n=n, m=m, col="orange")


### Ex 2:
n = 2; m = 1/2;
a = 1; b = 3/2;
ct = 2;
px = 6/7 + c(1, 4)*2/7
curve(y(x, a=a, b=b, ct=ct, n=n, m=m), from= 1, to = 2.5, ylim=c(1/2,7))
sapply(px, line.tan, dx=3, p=y, dp=dy, a=a, b=b, ct=ct, n=n, m=m)
#
curve(dy(x, a=a, b=b, ct=ct, n=n, m=m), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, a=a, b=b, ct=ct, n=n, m=m, col="orange")


############
### Example:
y = (x^2 + k*x)*log(x^n + a) + (x^2 - k*x)*log(x^n + b)

### D(y)
(2*x + k)*log(x^n + a) + (2*x - k)*log(x^n + b) + n*x^n*((x+k)/(x^n+a) + (x-k)/(x^n+b))
(2*x + k)*log(x^n + a) + (2*x - k)*log(x^n + b) +
	+ n*x^n*(2*x^(n+1) + (a+b)*x - k*(a-b))/((x^n+a)*(x^n+b))

### D2(y)
2*(log(x^n + a) + log(x^n + b)) +
	+ n*x^(n-1)*(2*x+k) / (x^n+a) + n*x^(n-1)*(2*x-k) / (x^n+b) +
	+ D( n*x^n*(2*x^(n+1) + (a+b)*x - k*(a-b))/((x^n+a)*(x^n+b)) );
2*(x*dy - y) / x^2 +
	- n*x^(n-1)*(2*x^(n+1) + (a+b)*x - k*(a-b))/((x^n+a)*(x^n+b)) + # from x*dy
	+ n*x^(n-1)*(2*x+k) / (x^n+a) + n*x^(n-1)*(2*x-k) / (x^n+b) +
	+ D( n*x^n*(2*x^(n+1) + (a+b)*x - k*(a-b))/((x^n+a)*(x^n+b)) );
# TODO: ...


############
### Example:
y = x^3*log(x^n + a) + x^2*log(x^n + b)

### D(y)
3*x^2*log(x^n + a) + 2*x*log(x^n + b) + n*x^(n+2)/(x^n+a) + n*x^(n+1)/(x^n+b)

### D2(y)
6*x*log(x^n + a) + 2*log(x^n + b) +
	+ 3*n*x^(n+1) / (x^n + a) + 2*n*x^n / (x^n + b) +
	+ D( n*x^(n+2)/(x^n+a) + n*x^(n+1)/(x^n+b) )
2*(2*dy - 3*y) / x^2 +
	+ 3*n*x^(n+1) / (x^n + a) + 2*n*x^n / (x^n + b) +
	+ D( n*x^(n+2)/(x^n+a) + n*x^(n+1)/(x^n+b) )


############
### Example:
y = sin(x)*log(p1) + cos(x)*log(p2)

### D(y)
cos(x)*log(p1) - sin(x)*log(p2) + sin(x) / p1 * dp1 + cos(x) / p2 * dp2

### D2(y)
- sin(x)*log(p1) - cos(x)*log(p2) +
	+ cos(x) / p1 * dp1 - sin(x) / p2 * dp2 +
	D( sin(x) / p1 * dp1 + cos(x) / p2 * dp2 )
- y +
	+ cos(x) / p1 * dp1 - sin(x) / p2 * dp2 +
	D( sin(x) / p1 * dp1 + cos(x) / p2 * dp2 )


#########################
#########################


### Section B: Non-Linear ODEs

### Derived from:
### y = (log(P1(x)))^2 + (log(P2(x)))^2

### Example:
y = (log(x^2 + a))^2 + (log(x^2 + b))^2

### D(y)
4*x*log(x^2 + a) / (x^2 + a) + 4*x*log(x^2 + b) / (x^2 + b)
### (x^2+a)*(x^2+b)*dy
4*x*(x^2+b)*log(x^2 + a) + 4*x*(x^2+a)*log(x^2 + b)

### (x^2+a)*(x^2+b)*D2(y) + D((x^2+a)*(x^2+b))*dy
4*(3*x^2 + b)*log(x^2 + a) + 4*(3*x^2 + a)*log(x^2 + b) +
	+ 8*x^2*(x^2+b) / (x^2 + a) + 8*x^2*(x^2+a) / (x^2 + b)
### Solve Linear system:
# T = 8*x^2*(x^2+b) / (x^2 + a) + 8*x^2*(x^2+a) / (x^2 + b);
log(x^2 + a) =
	(x*(x^2+a)*d2y - (3*x^2 + a)*dy - x*(x^2+a)*T) / (8*(a-b)*x^3)
log(x^2 + b) =
	(x*(x^2+b)*d2y - (3*x^2 + b)*dy - x*(x^2+b)*T) / -(8*(a-b)*x^3)
# TODO: check!

### ODE:
y +
	(x*(x^2+a)*d2y - (3*x^2 + a)*dy - x*(x^2+a)*T) *
	(x*(x^2+b)*d2y - (3*x^2 + b)*dy - x*(x^2+b)*T) / (64*x^6) # = 0

