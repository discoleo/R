

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

