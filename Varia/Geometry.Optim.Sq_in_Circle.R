########################
###
### Leonard Mada
### [the one and only]
###
### A Famous Geometry
### Optimization Problem
###
### draft v.0.1
### [very early draft]


### Problem:
### Maximize the radius of the inscribed Circle!


### Graphics:

par.old = par(mfrow=c(1, 2))

### Square: Type 1
plot(0,0, xlim=c(-3,3), ylim=c(-2,4))
### Square 1:
rect(-1,0,1,2)
x = seq(-2, 2, by=0.01)
lines(x, -sqrt(4 - x^2) + sqrt(3), col="green")
lines(x, sqrt(4 - x^2) + sqrt(3), col="green")
points(0, sqrt(3), col="green")
### Suqare 2:
sq.x = sqrt(4 - (3-sqrt(3))^2) # x = approx. 1.5;
rect(- sq.x, 2, - sq.x + 1, 3, border="blue")
### Circle 2:
# TODO: exact computations
# Note: students are invited to submit their solutions ahead of my solution!
r2 = 0.8496; c2.cx = 0.3298; c2.cy = 0.85;
c2.x = seq(c2.cx - r2, c2.cx + r2, by=0.01)
lines(c2.x, -sqrt(r2^2 - (c2.x-c2.cx)^2) + 2 + c2.cy, col="orange")
lines(c2.x, sqrt(r2^2 - (c2.x-c2.cx)^2) + 2 + c2.cy, col="orange")


### Square: Type 2
plot(0,0, xlim=c(-3,3), ylim=c(-2,4))

lines(c(0, sqrt(2), 0, -sqrt(2), 0), c(-sqrt(2), 0, sqrt(2), 0, -sqrt(2)))
lines(x, -sqrt(4 - x^2) + 2 - sqrt(2), col="green")
lines(x, sqrt(4 - x^2) + 2 - sqrt(2), col="green")
points(0, 2 - sqrt(2), col="green")

# Q: Which circle will be larger?
# Q: Can we optimize the size even further?


par(par.old)
