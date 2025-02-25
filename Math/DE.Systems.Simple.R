


### Simple Coupled Systems

# dy1 = b1*y1 + b2*y2 + f1;
# dy2 = c1*y1 + c2*y2 + f2;

# =>
c1*dy1 - b1*dy2 - c1*(b2*y2 + f1) + b1*(c2*y2 + f2) # = 0
c1*dy1 - b1*dy2 + (b1*c2 - b2*c1)*y2 - c1*f1 + b1*f2 # = 0


### D(Eq Sys 2) =>

### Case: b, c = constants;
d2y2 - c1*dy1 - c2*dy2 - df2 # = 0
d2y2 - (b1 + c2)*dy2 + (b1*c2 - b2*c1)*y2 - c1*f1 + b1*f2 - df2 # = 0

# => Solve ODE of Order 2;

# TODO: check;

### Case: b, c = f(x)
d2y2 - c1*dy1 - dc1*y1 - c2*dy2 - dc2*y2 - df2 # = 0
d2y2 - (b1 + c2)*dy2 - dc1*y1 +
	+ (b1*c2 - b2*c1 - dc2)*y2 - c1*f1 + b1*f2 - df2 # = 0
c1*d2y2 - c1*(b1 + c2)*dy2 - dc1*dy2 +
	+ c1*(b1*c2 - b2*c1 - dc2)*y2 + dc1*c2*y2 +
	- c1^2*f1 + b1*c1*f2 - c1*df2 + dc1*f2 # = 0

# TODO: check;

