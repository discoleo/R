########################
###
### Leonard Mada
### [the one and only]
###
### A Famous Geometry
### Optimization Problem
###
### draft v.0.3b
### [very early draft]


###############
### History ###

### draft v.0.3 - v.0.3b:
# - started work & estimations on the type 2 variant;
# - modified type 2 variant (v.0.3b);
#   [but needs correction]
### draft v.0.2:
# - improved estimation for Case 1;
### draft v.0.1:
# - initial description of cases;
# - basic plot;


#######################

### Problem:
### Maximize the radius of the inscribed Circle!

### Given:
# - a circle C1 of radius r1 = 2;
# - two squares of side length 2, respectively 1, are positioned
#   inside the circle C1 such that both touch the circle
#   (either with one or with 2 of their vertices);
# - the 2 squares slide on one of their sides;
# - a circle C2 is inscribed in C1, such that:
#   C2 is tangent to C1 and to all squares that it touches;
### Task:
# - Maximize the radius of the inscribed Circle (C2).

######################

### Graphics:

par.old = par(mfrow=c(1, 2))

### Square: Type 1
plot(0,0, xlim=c(-3,3), ylim=c(-2,4))
### Square 1:
rect(-1,0,1,2)
### Circle 1:
x = seq(-2, 2, by=0.01)
lines(x, -sqrt(4 - x^2) + sqrt(3), col="green")
lines(x, sqrt(4 - x^2) + sqrt(3), col="green")
points(0, sqrt(3), col="green")
### Square 2:
sq.x = sqrt(4 - (3-sqrt(3))^2) # x = approx. 1.5;
rect(- sq.x, 2, - sq.x + 1, 3, border="blue")
### Circle 2:
# Note: students are invited to submit their solutions ahead of my solution!
# TODO: exact computations;
r2 = 0.8462; c2.cx = r2 - sq.x + 1; c2.cy = r2;
c2.x = seq(c2.cx - r2, c2.cx + r2, by=0.01)
lines(c2.x, -sqrt(r2^2 - (c2.x-c2.cx)^2) + 2 + c2.cy, col="orange")
lines(c2.x, sqrt(r2^2 - (c2.x-c2.cx)^2) + 2 + c2.cy, col="orange")


### Square: Type 2
plot(0,0, xlim=c(-3,3), ylim=c(-2,4))
### Square 1:
lines(c(0, sqrt(2), 0, -sqrt(2), 0), c(-sqrt(2), 0, sqrt(2), 0, -sqrt(2)))
### Circle 1:
x = seq(-2, 2, by=0.01)
lines(x, -sqrt(4 - x^2) + 2 - sqrt(2), col="green")
lines(x, sqrt(4 - x^2) + 2 - sqrt(2), col="green")
points(0, 2 - sqrt(2), col="green")
### TODO: Square 2 + Circle 2;
dx = sqrt(2); dx2 = 0.195; dx2h = dx2/sqrt(2); sh = 1/sqrt(2);
lines(c(dx2h - dx, dx2h - dx - sh, dx2h - dx, dx2h - dx + sh, dx2h - dx),
	c(dx2h, dx2h + sh, dx2h + 2*sh, dx2h + sh, dx2h), col="blue")

# Q: Which circle will be larger?
# Q: Can we optimize the size even further?

par(par.old)

### Square: modified Type 2
# TODO: verify if boundaries are correct!
plot(0,0, xlim=c(-3,3), ylim=c(-2,4))
### Square 1:
a.s = 2/sqrt(13); a.c = 3/sqrt(13);
ab.s = 1/2 + 3/4*sqrt(3/13); ab.c = sqrt(1 - ab.s^2);
p4.s = 1/sqrt(2*13); p4.c = 5*p4.s;
lines(c(0, -2*ab.s, -2*sqrt(2)*p4.s, 2*ab.c, 0), c(0, 2*ab.c, 2*sqrt(2)*p4.c, 2*ab.s, 0) - 2)
### Circle 1:
x = seq(-2, 2, by=0.01)
lines(x, -sqrt(4 - x^2), col="green")
lines(x, sqrt(4 - x^2), col="green")
points(0, 0, col="green")
### TODO: Square 2 + Circle 2;


par(par.old)


#######################

################
### Analysis ###

### Why did the students fail to solve this problem?

# 1.) Why did the students fail to design any meaningful approaches to this problem?

# 2.) Why did the students NOT see the various topologies (types/configurations)?

# 3.) Why were the students not able to actually solve for the radius?

# 4.) Why were the students not able to optimize the radius?


#######################
#######################

### Experimental & Exploratory Approaches

# includes various other unrelated calculations

### unrelated half-circle
# cz.h = pi * roots(c(1, 2*(5-sqrt(3)-sq.x), (sq.x-3)*(sq.x+1) + (2-sqrt(3))^2))

### unrelated circle
dx = 2 * pi * roots(c(2, -4*sqrt(2), (3-2*sqrt(2))^2 + 1))
