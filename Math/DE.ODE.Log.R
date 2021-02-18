
########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Logarithms
###
### draft v.0.1a


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
	- n*x^(n-1)*(2*x^(m+n) + (a+b)*x^m - c*(a-b))

### D2(y)
# TODO: correct
(x^n + a)*(x^n + b) * d2y + n*x^(n-1)*(2*x^n + a + b) * dy +
	- (m*(m-1)*x^(m-2)*(x^n + a)*(x^n + b) + m*n*x^(m+n-2)*(2*x^n + a + b)) *
		(log(x^n + a) + log(x^n + b)) +
	- n*(n-1)*x^(n-2)*(2*x^(m+n) + (a+b)*x^m - c*(a-b)) +
	- n*x^(m+n-2)*(2*(m+n)*x^n + m*(a+b))
(x^n + a)*(x^n + b) * d2y + # TODO: Error: - n*T !
	- ((m-1)*x^(-1)*(x^n + a)*(x^n + b)) * dy +
	# TODO: Error: expand - n*T !
	- n * T * (m*(m-1)*x^(m-2)*(x^n + a)*(x^n + b) + m*n*x^(m+n-2)*(2*x^n + a + b)) +
	- n*(n-1)*x^(n-2)*(2*x^(m+n) + (a+b)*x^m - c*(a-b)) +
	- n*x^(m+n-2)*(2*(m+n)*x^n + m*(a+b))


### Example 1:
m = 2
#
x*(x^n + a)*(x^n + b) * d2y +
	- (x^n + a)*(x^n + b) * dy +
	- 2*c*(a-b)*n*(2*n+1)*x^(2*n+1) +
	- (a+b)*n*(n+1)*x^(n+1) + n*(n-1)*x^(n-1)
	- n * T * (2*(x^n + a)*(x^n + b) + 2*n*x^n*(2*x^n + a + b))

### m = 2; a+b = 0;
x*(x^n + a)*(x^n - a) * d2y +
	- (x^n + a)*(x^n - a) * dy +
	- 4*c*a*n*(2*n+1)*x^(2*n+1) + n*(n-1)*x^(n-1)

### Solution & Plot:
y = function(x, a, b=-a, ct=0, n=2, m=2) {
	xm = x^m; xn = x^n;
	x.prod = (xn + a)*(xn +b); x.div = (xn + a) / (xn + b);
	val = xm*log(x.prod) + ct*log(x.div);
	return(val)
}
dy = function(x, a, b=-a, ct=0, n=2, m=2) {
	xm = x^m; xn = x^n;
	# y.x = y(x, a=a, b=b, c=c, n=n, m=m);
	x.prod = (xn + a)*(xn + b);
	div = x;
	Tnx = n*xn*(2*xm*xn + (a+b)*xm - ct*(a-b)) / (x.prod);
	dp = m*xm*log(x.prod) + Tnx;
	dp = ifelse(div != 0, dp/div, 1); # TODO
	return(dp)
}
d2y = function(x, a, b=-a, ct=0, n=2, m=2) {
	# TODO: need to correct first formulas;
}
### Plot:
n = 2; m = 2;
a = 1; b = -a;
ct = 2;
px = 5/7 + (1:3)*3/7
curve(y(x, a=a, b=b, ct=ct, n=n, m=m), from= 1, to = 2.5, ylim=c(-2,20))
sapply(px, line.tan, dx=3, p=y, dp=dy, a=a, b=b, ct=ct, n=n, m=m)
# x*log(x)-like:
curve(dy(x, a=a, b=b, ct=ct, n=n, m=m), add=T, col="green")
sapply(px, line.tan, dx=3, p=dy, dp=d2y, a=a, b=b, ct=ct, n=n, m=m, col="orange")

