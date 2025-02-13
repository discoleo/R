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

source("DE.ODE.Helper.R")
source("DE.ODE.Helper.Math.R")


#####################
#####################

### Base System:
# sin(k*x*y1^n) = b3*y1*y2 + b2*y2 + b1*y1 + f1
# cos(k*x*y1^n) = c3*y1*y2 + c2*y2 + c1*y1 + f2

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
# Eq 1:
k*x*((b3^2 + c3^2)*y1*y2 + k*x*(f1*b3 + f2*c3))*dy1 +
	+ k*(b3^2 + c3^2)*y2*y1^2 + (b3*dc3 - c3*db3)*y1*y2 +
	+ k*(f1*b3 + f2*c3)*y1 - c3*df1 + b3*df2 # = 0
# Eq 2:
k*x*((b3^2 + c3^2)*y2*y1^2 + (f1*b3 + f2*c3)*y1)*dy2 +
	+ k*(x*(b3*db3 + c3*dc3) - b3^2 - c3^2)*y1^2*y2^2* +
	- (b3*dc3 - db3*c3)*y1*y2^2 +
	+ k*(x*(c3*df2 + f2*dc3 + db3*f1 + df1*b3) - f1*b3 - f2*c3)*y1*y2 +
	+ df1*y2*c3 - y2*b3*df2 + k*x*(f1*df1 + f2*df2) # = 0

# TODO: check;


###
# source("Polynomials.Helper.R")

# n = 1
# p1 = as.pm("...")
p1 = replace.pm(p1, n, xn = "n")
p2 = replace.pm(p2, n, xn = "n")
tmp = solve.pm(p1, p2, by = "dy2")

