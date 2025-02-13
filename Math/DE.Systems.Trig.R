########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## DE Systems: Trigonometric
##
## draft v.0.1a


####################

### Helper Functions

library(deSolve)

source("DE.ODE.Helper.R")
source("DE.ODE.Helper.Math.R")


#####################
#####################

### Base System:
# sin(k*x*y1^n) = b3*y1*y2 + b2*y2 + b1*y1 + f1
# cos(k*x*y1^n) = c3*y1*y2 + c2*y2 + c1*y1 + f2

# Note:
# - y2 can be substituted from Eq 2 into Eq 1;

### D =>
# Note: k = constant;
# Eq 1:
k*(n*x*y1^(n-1)*dy1 + y1^n)*cos(k*x*y1^n) +
	- b3*(y2*dy1 + y1*dy2) - b2*dy2 - b1*dy1 +
	- db3*y1*y2 - db2*y2 - db1*y1 - df1 # = 0
# Eq 2:
k*(n*x*y1^(n-1)*dy1 + y1^n)*sin(k*x*y1^n) +
	+ c3*(y2*dy1 + y1*dy2) + c2*dy2 + c1*dy1 +
	+ dc3*y1*y2 + dc2*y2 + dc1*y1 + df2 # = 0

# =>
# Eq 1:
k*(n*x*y1^(n-1)*dy1 + y1^n)*(c3*y1*y2 + c2*y2 + c1*y1 + f2) +
	- b3*(y2*dy1 + y1*dy2) - b2*dy2 - b1*dy1 +
	- db3*y1*y2 - db2*y2 - db1*y1 - df1 # = 0
# Eq 2:
k*(n*x*y1^(n-1)*dy1 + y1^n)*(b3*y1*y2 + b2*y2 + b1*y1 + f1) +
	+ c3*(y2*dy1 + y1*dy2) + c2*dy2 + c1*dy1 +
	+ dc3*y1*y2 + dc2*y2 + dc1*y1 + df2 # = 0

# TODO: check;


### Special Cases:

### SC 1:
# b1 = b2 = 0;
# c1 = c2 = 0;
# Eq 1:
k*(n*x*y1^(n-1)*dy1 + y1^n)*(c3*y1*y2 + f2) +
	- b3*(y2*dy1 + y1*dy2) - db3*y1*y2 - df1 # = 0
# Eq 2:
k*(n*x*y1^(n-1)*dy1 + y1^n)*(b3*y1*y2 + f1) +
	+ c3*(y2*dy1 + y1*dy2) + dc3*y1*y2 + df2 # = 0

# Solve for dy1, dy2 =>
# (for n = 1)
# Eq 1:
k*x*((b3^2 + c3^2)*y1*y2 + (f1*b3 + f2*c3))*dy1 +
	+ k*(b3^2 + c3^2)*y2*y1^2 + (b3*dc3 - c3*db3)*y1*y2 +
	+ k*(f1*b3 + f2*c3)*y1 - c3*df1 + b3*df2 # = 0
# Eq 2:
k*x*((b3^2 + c3^2)*y1*y2 + (f1*b3 + f2*c3))*y1*dy2 +
	+ k*(x*(b3*db3 + c3*dc3) - b3^2 - c3^2)*y1^2*y2^2 +
	- (b3*dc3 - db3*c3)*y1*y2^2 +
	+ k*(x*(c3*df2 + f2*dc3 + db3*f1 + df1*b3) - f1*b3 - f2*c3)*y1*y2 +
	- (b3*df2 - c3*df1)*y2 + k*x*(f1*df1 + f2*df2) # = 0

# TODO: check;


### SC 2:
# b3 = b; c3 = -b; & SC1;
# Eq 1:
k*x*(2*b*y1*y2 + (f1 - f2))*dy1 +
	+ 2*k*b*y2*y1^2 + k*(f1 - f2)*y1 + (df1 + df2) # = 0
# Eq 2:
k*x*b*(2*b*y1*y2 + (f1 - f2))*y1*dy2 +
	+ 2*k*b*(x*db - b)*y1^2*y2^2 +
	+ k*((x*db - b)*(f1 - f2) + x*b*(df1 - df2))*y1*y2 +
	- b*(df1 + df2)*y2 + k*x*(f1*df1 + f2*df2) # = 0


### Examples:

dyf = function(x, y, parms, ...) {
	params = parms;
	k = params$k;
	b = params$b(x); db = params$db(x);
	y1 = y[1]; y2 = y[2];
	f1 = params$f1(x); df1 = params$df1(x);
	f2 = params$f2(x); df2 = params$df2(x);
	div1 = k*x*(2*b*y1*y2 + (f1 - f2));
	dy1  = 2*k*b*y2*y1^2 + k*(f1 - f2)*y1 + (df1 + df2);
	dy1  = - dy1 / div1;
	div2 = k*x*b*(2*b*y1*y2 + (f1 - f2))*y1;
	dy2  = 2*k*b*(x*db - b)*y1^2*y2^2 +
		+ k*((x*db - b)*(f1 - f2) + x*b*(df1 - df2))*y1*y2 +
		- b*(df1 + df2)*y2 + k*x*(f1*df1 + f2*df2);
	dy2 = - dy2 / div2;
	# if(abs(dy2) < 1E+10) cat(x, ", ", div1, ", ", y1,
	#	", ", dy1, ", ", dy2, "\n");
	return(list(c(dy1, dy2)));
}

bf  = \(x) x^2 / 8 + 1/32; dbf = \(x) x / 4;
f1f = \(x) (3*x + 1) / 24; df1f = \(x) 1/8;
f2f = \(x) (1/3 - 3*x + x^2 / 4) / 12; df2f = \(x) -3/12 + x/24;
params = list(k = 2, b = bf, db = dbf, f1=f1f, f2=f2f, df1=df1f, df2=df2f);

x = seq(0.25, 1, length.out = 100); x0 = x[1];
y0 = (asin((f1f(x0) + f2f(x0))/sqrt(2)) - pi/4 + 2*pi) / (k*x0);
y0 = c(y0, (sin(k*x0*y0) - f1f(x0)) / (bf(x0)*y0));
# y0 = c(2, 1);
sol = ode(y0, times = x, func = dyf, parms = params)

# Plot:
matplot(sol[,1], sol[,-1], type = "l")

# Test:
err = sin(params$k * x*sol[,2]) +
	- (bf(x)*sol[,2]*sol[,3] + f1f(x));
summary(err)

#
id = 27; x = x[id]; k = params$k;
y1 = sol[id,2]; y2 = sol[id,3]; b = bf(x); db = dbf(x);
f1 = f1f(x); df1 = df1f(x); f2 = f2f(x); df2 = df2f(x);
dt = sol[id+1,1] - sol[id,1];
dy = dyf(x, c(y1, y2), params)[[1]]; dy1 = dy[1]; dy2 = dy[2]; 
n = 1; b2=b1=c2=c1 = 0; db2=db1=dc2=dc1 = 0;
b3 = b; c3 = -b; db3 = db; dc3 = - db;


###
# source("Polynomials.Helper.R")

# n = 1
# p1 = as.pm("...")
p1 = replace.pm(p1, n, xn = "n")
p2 = replace.pm(p2, n, xn = "n")
tmp = solve.pm(p1, p2, by = "dy2")
tmp = tmp$Rez; tmp$coeff = - tmp$coeff;
as.coef.pm(tmp, "dy1")

