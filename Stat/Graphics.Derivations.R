

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

# TODO: special cases;
intersect.lines = function(xA, yA, xB, yB) {
	xA1 = xA[1]; xA2 = xA[2]; xB1 = xB[1]; xB2 = xB[2];
	yA1 = yA[1]; yA2 = yA[2]; yB1 = yB[1]; yB2 = yB[2];
	#
	div = (yB2 - yB1)*(xA2 - xA1) + (yA2 - yA1)*(xB1 - xB2);
	t1  = (yA2 - yB1)*xB2 - (yB2 - yB1)*xA2 + (yB2 - yA2)*xB1;
	t1  = - t1 / div;
	#
	div = yB[1] - yB[2];
	t2  = yA2 - yB2 + (yA1 - yA2)*t1;
	t2  = t2 / div;
	#
	x = xA1*t1 + (1-t1)*xA2;
	y = yA1*t1 + (1-t1)*yA2;
	return(list(x=x, y=y, t1=t1, t2=t2));
}
