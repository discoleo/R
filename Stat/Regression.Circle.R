

### Regression on a Circle/Sphere

### Functionality:
# - Regression Circle (2D):
#   for points distributed in a circle or circle-segment;
# - Regression on a 3D Sphere;
# - Extendible (easily) to spheres in higher dimensions;


### TODO:

# Explore:
# 1. circular: Circular Statistics
#    https://cran.r-project.org/web/packages/circular/index.html
# 2. circularKDE: Recent Methods for Kernel Density Estimation of Circular Data
#    https://cran.r-project.org/web/packages/circularKDE/index.html
# 3. circda: Circular Data Analysis
#    https://cran.r-project.org/web/packages/circda/index.html


####################

### Helper Functions

# Regression on Circle (2D)
lm.circle = function(x, y) {
	xm = mean(x); ym = mean(y);
	dx = x - xm; dy = y - ym;
	n  = length(x);
	Sdx = sum((x^2 + y^2)*dx);
	Sdy = sum((x^2 + y^2)*dy);
	#
	divx = 2*(sum(x*dx)* sum(dy*y) - sum(x*dy) * sum(dx*y));
	x0 = (Sdx * sum(dy*y) - Sdy * sum(dx*y)) / divx;
	#
	y0 = - (Sdx * sum(dy*x) - Sdy * sum(dx*x)) / divx;
	# Radius:
	r = sum((x-x0)^2, (y-y0)^2);
	r = sqrt(r/n);
	#
	fit = list(C = c(x0, y0), R = r);
	return(fit);
}
# 3D: Regression on 3D Sphere
lm.sphere = function(x, y, z) {
	xm = mean(x); ym = mean(y); zm = mean(z);
	dx = x - xm; dy = y - ym; dz = z - zm;
	n  = length(x); # Number of data-points;
	coefXYZ = 2 * c(
		sum(dx*x), sum(dx*y), sum(dx*z),
		sum(dy*x), sum(dy*y), sum(dy*z),
		sum(dz*x), sum(dz*y), sum(dz*z)
	);
	coefXYZ = matrix(coefXYZ, nrow = 3, byrow=T);
	mY = c(sum(dx*(x^2 + y^2 + z^2)), sum(dy*(x^2 + y^2 + z^2)), sum(dz*(x^2 + y^2 + z^2)));
	cc = solve(coefXYZ, mY);
	# Radius:
	r = sum((x-x0)^2, (y-y0)^2, (z-z0)^2);
	r = sqrt(r/n);
	rez = list(C = cc, R = r);
	return(rez);
}

### Line passing through Center of Circle
# - Maximizes/Minimizes (?) distance to points;
# - BUT difficult to interpret / non-intuitive;
lm.C2line = function(x, y, xyC) {
	dx = x - xyC[1];
	dy = y - xyC[2];
	a = sum(dy^2 - dx^2) / sum(dx*dy) / 2;
	b = xyC[2] - a*xyC[1];
	# y = a*x + b;
	# Note: Reverse order to abline(coef = ...);
	return(list(a = a, b = b));
}

# Distance to Line
dist.line = function(x, y, a, b) {
	if(missing(b)) {
		b = a[[2]]; a = a[[1]];
	}
	sum((y - a*x - b)^2) / (a^2 + 1)
}

### Plot:
lines.circle = function(r, centre = c(0,0), ..., n = 33) {
	phi = seq(0, 2*pi, length.out = n);
	phi[n] = phi[1];
	x = r*sin(phi) + centre[1];
	y = r*cos(phi) + centre[2];
	lines(x, y, ...);
	invisible(cbind(x,y));
}


####################
####################

### Examples:

### Ex 1:
n = 30;
r = 10;
x0 = 5; y0 = 2;
# Random data:
phi = (2*pi)*runif(n);
dri =  runif(n, -1,1) / 10 * r;
x = x0 + (r+dri) * sin(phi);
y = y0 + (r+dri) * cos(phi);

# Fit
fit = lm.circle(x, y);
Cxy = fit$C;
phi.fit = atan2(x - Cxy[1], y - Cxy[2]);
# Mean Vector:
Vm = c(sum(sin(phi.fit)), sum(cos(phi.fit)));

# Plot:
plot(x, y, xlim = range(x) + c(-1,1), ylim=range(y) + c(-1,1) )
lines.circle(fit$R, fit$C, col = "#C28092D2")
points(fit$R * sin(phi.fit) + Cxy[1], fit$R * cos(phi.fit) + Cxy[2], col = "red")
# Mean Vector: often very limited cancellation
lines(rbind(Cxy, fit$R*Vm + Cxy), col = "purple", lwd=2)
sqrt(sum(Vm^2)); # Multiples of R;


#########
### Ex 2:
n = 30;
r = 10;
x0 = 5; y0 = 2;
# Random data:
phi = (2*pi)*runif(n, 0, 3/5);
dri =  runif(n, -1,1) / 10 * r;
x = x0 + (r+dri) * sin(phi);
y = y0 + (r+dri) * cos(phi);

# Fit
fit = lm.circle(x, y);
Cxy = fit$C;
phi.fit = atan2(x - Cxy[1], y - Cxy[2]);

# Plot:
plot(x, y, xlim = range(x) + c(-1,1), ylim=range(y) + c(-1,1), asp=1)
lines.circle(fit$R, fit$C, col = "#C28092D2")
points(fit$R * sin(phi.fit) + Cxy[1], fit$R * cos(phi.fit) + Cxy[2], col = "red")


#########
### Ex 3:
n = 30;
r = 10;
x0 = 5; y0 = 2;
# Random data:
phi = (2*pi)*runif(n, -1/4, 1/4);
dri =  runif(n, -1,1) / 10 * r;
x = x0 + (r+dri) * sin(phi);
y = y0 + (r+dri) * cos(phi);

# Fit
fit = lm.circle(x, y);
Cxy = fit$C;
phi.fit = atan2(x - Cxy[1], y - Cxy[2]);

# Plot:
plot(x, y, xlim = range(x) + c(-1,1), ylim=range(y) + c(-1,1), asp=1)
lines.circle(fit$R, fit$C, col = "#C28092D2")
points(fit$R * sin(phi.fit) + Cxy[1], fit$R * cos(phi.fit) + Cxy[2], col = "red")
# Line through Center
# - Maximizes (?) distance to line;
lc = lm.C2line(x, y, Cxy)
abline(lc$b, lc$a, col = "pink")
dist.line(x, y, lc)


### Higher-Order Correlations

## TODO: perform some experiments or some theoretical approach;


#######################
#######################

### Derivation:

### Circular Regression

#
xm = mean(x); ym = mean(y);
dx = x - xm; dy = y - ym;

# Eq 1:
sum((x-x0)^2 + (y-y0)^2 - r^2) # = 0
r^2 - sum((x-x0)^2 + (y-y0)^2) / n # = 0

# Eq 2 & 3:
sum(((x-x0)^2 + (y-y0)^2 - r^2) * (x-x0)) # = 0
sum(((x-x0)^2 + (y-y0)^2 - r^2) * (y-y0)) # = 0

# Substituting Eq 1 in Eqs 2&3:
sum(((x-x0)^2 + (y-y0)^2 - r^2) * x) # = 0
sum(((x-x0)^2 + (y-y0)^2 - r^2) * y) # = 0
# =>
sum((x-x0)^2*x + (y-y0)^2*x - xm*(x-x0)^2 - xm*(y-y0)^2) # = 0
sum((x-x0)^2*y + (y-y0)^2*y - ym*(x-x0)^2 - ym*(y-y0)^2) # = 0

# =>
sum(x^2*(x-xm) - 2*x*(x-xm)*x0 + (x-xm)*x0^2 + (x-xm)*y^2 - 2*(x-xm)*y*y0 + (x-xm)*y0^2) # = 0
sum(x^2*(y-ym) - 2*x*(y-ym)*x0 + (y-ym)*x0^2 + (y-ym)*y^2 - 2*(y-ym)*y*y0 + (y-ym)*y0^2) # = 0

# =>
sum(x^2*dx + dx*y^2 - 2*x*dx*x0 - 2*dx*y*y0) # = 0
sum(x^2*dy + dy*y^2 - 2*x*dy*x0 - 2*dy*y*y0) # = 0
# =>
sum(x^2*dx + dx*y^2 - 2*x*dx*x0) * sum(dy*y) - sum(x^2*dy + dy*y^2 - 2*x*dy*x0) * sum(dx*y) # = 0
#
divx = 2*(sum(x*dx) * sum(dy*y) - sum(x*dy) * sum(dx*y));
(sum(x^2*dx + dx*y^2) * sum(dy*y) - sum(x^2*dy + dy*y^2) * sum(dx*y)) / divx;
(sum((x^2 + y^2)*dx) * sum(dy*y) - sum((x^2 + y^2)*dy) * sum(dx*y)) / divx;
x0;
#
divy = 2*(sum(dx*y)* sum(dy*x) - sum(dy*y) * sum(dx*x));
(sum(x^2*dx + dx*y^2) * sum(dy*x) - sum(x^2*dy + dy*y^2) * sum(dx*x)) / divy;
y0;

# Check/Plot:
xx = (sum(x^2*dx + dx*y^2) * sum(dy*y) - sum(x^2*dy + dy*y^2) * sum(dx*y)) / divx;
yy = (sum(x^2*dx + dx*y^2) * sum(dy*x) - sum(x^2*dy + dy*y^2) * sum(dx*x)) / divy;
plot(x, y, xlim = range(x) + c(-1,1), ylim=range(y) + c(-1,1) )
points(r*sin(phi) + xx, r*cos(phi) +y0, col = "red")


### Line through Center
# yc = a*xc + b;
b = yc - xc*a;
a * sum((y-yc - a*(x - xc))^2) + (a^2 + 1)*sum((x-xc) * (y-yc - a*(x - xc)));
# Minimize:
sum((y-yc - a*(x - xc))^2) - 2*a*sum((y-yc - a*(x - xc))*(x - xc)) +
	+ 2*a*sum((x-xc) * (y-yc - a*(x - xc))) - (a^2 + 1)*sum((x-xc)^2);
sum((y-yc)^2 - (x-xc)^2) - 2*a*sum((x - xc)*(y-yc));
a = sum((y-yc)^2 - (x-xc)^2) / sum((x - xc)*(y-yc)) / 2;

# Distance:
sum((y-yc - a*(x-xc))^2) / (a^2 + 1);
2*a*sum((y-yc - a*(x-xc))^2) / (a^2 + 1)^2 + 2*sum((y-yc - a*(x-xc))*(x-xc)) / (a^2 + 1);
a*sum((y-yc - a*(x-xc))^2) + (a^2 + 1)*sum((y-yc - a*(x-xc))*(x-xc));


###################
###################

### 3D

library(rgl)

### Ex 1:
n = 20;
r = 10;
x0 = 5; y0 = 2; z0 = 1;
# Random data:
phi = (2*pi)*runif(n);
dri =  runif(n, -1,1) / 10 * r;
dz  = seq(0, r, length.out = 5) * 0.8;
# dz  = seq(-r, r, length.out = 5) * 0.8;
x = numeric(0); y = numeric(0);
z = rep(z0 + dz, each = n);
for(zi in dz) {
	dphi = pi * runif(1);
	rz = sqrt(r^2 - zi^2);
	x = c(x, x0 + (rz+dri) * sin(phi + dphi));
	y = c(y, y0 + (rz+dri) * cos(phi + dphi));
}

# Fit:
lm.fit = lm.sphere(x, y, z);
C0 = lm.fit$C; r0 = lm.fit$R;
px = x - C0[1]; py = y - C0[2]; pz = z - C0[3];
rp = sqrt(px^2 + py^2 + pz^2);
px = C0[1] + r0*px/rp;
py = C0[2] + r0*py/rp;
pz = C0[3] + r0*pz/rp;

# Plot:
points3d(x, y, z, size=5)
points3d(px, py, pz, size=5, col = "red", alpha = 0.7)


### Derivation:

### Spherical Regression

#
xm = mean(x); ym = mean(y); zm = mean(z);
dx = x - xm; dy = y - ym; dz = z - zm;

# Eq 1:
sum((x-x0)^2 + (y-y0)^2 + (z-z0)^2 - r^2) # = 0
r^2 - sum((x-x0)^2 + (y-y0)^2 + (z-z0)^2) / n # = 0

# Eq 2, 3, 4:
sum(((x-x0)^2 + (y-y0)^2 + (z-z0)^2 - r^2) * (x-x0)) # = 0
sum(((x-x0)^2 + (y-y0)^2 + (z-z0)^2 - r^2) * (y-y0)) # = 0
sum(((x-x0)^2 + (y-y0)^2 + (z-z0)^2 - r^2) * (z-z0)) # = 0

# Substituting Eq 1 in Eqs 2-4:
sum(((x-x0)^2 + (y-y0)^2 + (z-z0)^2 - r^2) * x) # = 0
sum(((x-x0)^2 + (y-y0)^2 + (z-z0)^2 - r^2) * y) # = 0
sum(((x-x0)^2 + (y-y0)^2 + (z-z0)^2 - r^2) * z) # = 0
# =>
sum((x-x0)^2*x + (y-y0)^2*x + (z-z0)^2*x - xm*(x-x0)^2 - xm*(y-y0)^2 - xm*(z-z0)^2) # = 0
sum((x-x0)^2*y + (y-y0)^2*y + (z-z0)^2*y - ym*(x-x0)^2 - ym*(y-y0)^2 - ym*(z-z0)^2) # = 0
sum((x-x0)^2*z + (y-y0)^2*z + (z-z0)^2*z - zm*(x-x0)^2 - zm*(y-y0)^2 - zm*(z-z0)^2) # = 0

=>
sum(dx*(x^2 + y^2 + z^2) - 2*dx*(x*x0 + y*y0 + z*z0)) # = 0
sum(dy*(x^2 + y^2 + z^2) - 2*dy*(x*x0 + y*y0 + z*z0)) # = 0
sum(dz*(x^2 + y^2 + z^2) - 2*dz*(x*x0 + y*y0 + z*z0)) # = 0

coefXYZ = 2 * c(
	sum(dx*x), sum(dx*y), sum(dx*z),
	sum(dy*x), sum(dy*y), sum(dy*z),
	sum(dz*x), sum(dz*y), sum(dz*z)
	);
coefXYZ = matrix(coefXYZ, nrow = 3, byrow=T);
mY = c(sum(dx*(x^2 + y^2 + z^2)), sum(dy*(x^2 + y^2 + z^2)), sum(dz*(x^2 + y^2 + z^2)));
cc = solve(coefXYZ, mY);
print(cc);
