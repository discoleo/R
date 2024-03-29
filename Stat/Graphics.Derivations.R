

source("Polynomials.Helper.R")

### Intersection of 2 Lines:
p1 = toPoly.pm("(xA1 - xA2)*t1 - (xB1 - xB2)*t2 + xA2 - xB2");
p2 = toPoly.pm("(yA1 - yA2)*t1 - (yB1 - yB2)*t2 + yA2 - yB2");
pR = solve.pm(p1, p2, "t2")
pR$Rez = sort.pm(pR$Rez, "t1")
print.pm(pR$Rez, do.sort=FALSE, leading="t1")

### Parametric Eq:
((yB2 - yB1)*(xA2 - xA1) + (yA2 - yA1)*(xB1 - xB2))*t1 +
	+ (yA2 - yB1)*xB2 - (yB2 - yB1)*xA2 + (yB2 - yA2)*xB1 # = 0

### Intersect Function

intersect.lines = function(xA, yA, xB, yB) {
	xA1 = xA[1]; xA2 = xA[2]; xB1 = xB[1]; xB2 = xB[2];
	yA1 = yA[1]; yA2 = yA[2]; yB1 = yB[1]; yB2 = yB[2];
	#
	yB12 = yB1 - yB2;
	yA12 = yA1 - yA2;
	#
	div = yB12*(xA1 - xA2) - yA12*(xB1 - xB2);
	# TODO: remaining special cases;
	if(div == 0) {
		# Parallel lines
		d = (xA1*yA2 - xA2*yA1)*(xB1 - xB2) - (xB1*yB2 - xB2*yB1)*(xA1 - xA2);
		if(d == 0) {
			# Overlapping lines
			t1 = c(xB1 - xA2, xB2 - xA2) / (xA1 - xA2);
			if(t1[1] > t1[2]) { tmp = t1[1]; t1[1] = t1[2]; t1[2] = tmp; }
			x = NA; y = NA; t1c = NA;
			if(t1[1] < 0) {
				if(t1[2] >= 0) {
					t1c = c(0, min(1, t1[2]));
				}
			} else if(t1[1] <= 1) {
				t1c = c(t1[1], min(1, t1[2]));
			}
			if( ! is.na(t1c[1])) {
				x = xA1*t1c + (1-t1c)*xA2;
				y = yA1*t1c + (1-t1c)*yA2;
			}
			return(list(x=x, y=y, t1=t1, t2=Inf, n=2));
		} else {
			return(list(x=NA, y=NA, t1=Inf, t2=Inf, n=0));
		}
	}
	# Computations:
	t1  = (yA2 - yB1)*xB2 + yB12*xA2 + (yB2 - yA2)*xB1;
	t1  = - t1 / div;
	#
	div = yB12;
	t2  = yA2 - yB2 + (yA1 - yA2)*t1;
	# TODO: div == 0
	t2  = t2 / div;
	#
	x = xA1*t1 + (1-t1)*xA2;
	y = yA1*t1 + (1-t1)*yA2;
	return(list(x=x, y=y, t1=t1, t2=t2, n=1));
}

### Test
xA = c(1,6); yA = c(1,8);
xB = c(1,7); yB = c(4,2);
plot.new(); plot.window(xlim=c(-2,10), ylim=c(-2,10))
lines(xA, yA); lines(xB, yB, col="blue");
p = intersect.lines(xA, yA, xB, yB)
points(p$x, p$y, col="green")

###
xA = c(1,3); yA = c(1,5);
xB = c(2,7); yB = c(3,13);
plot.new(); plot.window(xlim=c(-2,10), ylim=c(-2,10))
lines(xA, yA); lines(xB, yB, col="blue");
p = intersect.lines(xA, yA, xB, yB)
points(p$x, p$y, col="green")

###
xA = c(1,3); yA = c(1,5);
xB = c(4,7); yB = c(7,13);
plot.new(); plot.window(xlim=c(-2,10), ylim=c(-2,10))
lines(xA, yA); lines(xB, yB, col="blue");
p = intersect.lines(xA, yA, xB, yB)
points(p$x, p$y, col="green")

###
xA = c(1,7); yA = c(1,13);
xB = c(3,5); yB = c(5,9);
plot.new(); plot.window(xlim=c(-2,10), ylim=c(-2,10))
lines(xA, yA); lines(xB, yB, col="blue");
p = intersect.lines(xA, yA, xB, yB)
points(p$x, p$y, col="green")


############################

### Ellipse

### Eq:
# th = angle of rotation;
x^2*(cos(th)^2/a^2 + sin(th)^2/b^2) + y^2*(sin(th)^2/a^2 + cos(th)^2/b^2) +
	+ (1/a^2 - 1/b^2)*x*y*sin(2*th) - r^2 # = 0


# naive (mathematical) plot:
plotEllipse = function(a, b, theta=0, center=c(0,0), r=1, N=128, ...) {
	B0 = cos(theta)^2/a^2 + sin(theta)^2/b^2;
	B1 = (1/a^2 - 1/b^2)*sin(2*theta);
	B2 = sin(theta)^2/a^2 + cos(theta)^2/b^2;
	B1sq = B1^2;
	x0 = 2*r*sqrt(B2 / (4*B2*B0 - B1sq));
	x = seq(-x0, x0, length.out = N);
	y = lapply(x, function(x) {
		# numerical stability:
		d = B1sq*x^2 - 4*B2*(B0*x^2 - r^2);
		d = if(d <= 0) 0 else sqrt(d);
		y = c(- B1*x + d, -B1*x - d) / (2*B2);
		return(y);
	})
	y = do.call(rbind, y);
	y = c(y[,1], rev(y[,2]));
	x = c(x, rev(x));
	lines(center[1] + x, center[2] + y, ...);
}

solveEllipse = function(p1, p2, slope=0, r=1) {
	if(slope == 0) {
		div = (p2[2]^2 - p1[2]^2) / (p1[1]^2 - p2[1]^2);
		a = sqrt(p1[1]^2 + p1[2]^2 / div) / r;
		b = a * sqrt(div);
		return(list(a=a, b=b));
	}
	# TODO
}


### Test
# plot.base = open a new plot window;
theta = pi/4
plot.base(xlim=c(-10,10), ylim=c(-10,10))
shape::plotellipse(rx=3, ry=1, angle= theta*180 / pi, col="red")
plotEllipse(3, 1, theta=theta, col="green")


###
theta = pi/7
plot.base(xlim=c(-10,10), ylim=c(-10,10))
shape::plotellipse(rx=3, ry=1, angle= theta*180 / pi, col="red")
plotEllipse(3, 1, theta=theta, col="green")


###
theta = pi/7
center = c(1,2)
plot.base(xlim=c(-10,10), ylim=c(-10,10))
shape::plotellipse(rx=3, ry=1, angle= theta*180 / pi, mid=center, col="red")
plotEllipse(3, 1, theta=theta, center=center, col="green")


### Solve

### Ex 1:
p1 = c(0,3)
p2 = c(1, 3/2)
sol = solveEllipse(p1, p2)
#
plot.base(xlim=c(-10,10), ylim=c(-10,10))
shape::plotellipse(rx=sol$a, ry=sol$b, angle=0, lcol="red")
plotEllipse(sol$a, sol$b, theta=0, col="green")
px = t(cbind(p1, p2))
points(px, col="blue")


### Ex 2:
p1 = c(-2,-1)
p2 = c(1, 3/2)
sol = solveEllipse(p1, p2)
#
plot.base(xlim=c(-10,10), ylim=c(-10,10))
shape::plotellipse(rx=sol$a, ry=sol$b, angle=0, lcol="red")
plotEllipse(sol$a, sol$b, theta=0, col="green")
px = t(cbind(p1, p2))
points(px, col="blue")

