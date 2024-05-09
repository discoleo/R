

source("Chem3D.R")


### Examples

### Basic Math & Transformations

### Normal
p1 = c(1,5,7); p2 = c(3,5,7);
pL = matrix(c(-4,6,-5,8,-3,8), nrow=2)
test.rotate.ortho(p1, pL)
test.rotate.ortho(p2, pL, col.point = "#D68016")


### Shift Triangle
p3 = matrix(c(1,-4,6, 5,-5,8, 7,-3,8), nrow=3)
d = 3 * seq(4)
test.shift.poly(d, p3)

###
p4 = matrix(c(1,-4,-3, 5,-5,-6, 7,-3,1), nrow=3)
p4 = rbind(p4, cbind(0.8,-0.4,0.6) %*% p4)
d = 3 * seq(4)
test.shift.poly(d, p4)



### Rotate point
p = c(2,-3,17)
pL = matrix(c(-4,6,-5,8,-3,8), nrow=2)

test.rotate.point(p, pL)
test.rotate.point(c(1,1,-10), pL)


### Arbitrary Normal

###
p = matrix(c(1,5,1,3,2,-3), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)

###
pL = matrix(c(-4,6,-5,8,-3,8), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)

### Special Case:
pL = matrix(c(-4,4,-1,8,-5,4), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)

### Special Case:
pL = matrix(c(-4,1,-1,5,-5,1), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)

###
pL = matrix(c(-4,1,-5,3,-5,2), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)


#############

### Examples: Cylinder

###
p = matrix(c(1,5,1,3,2,-3), nrow=2)
m = mesh.cylinder(1, p)
points3d(m$V)
lines3d(p, col = "blue")
lines3d(rbind(p[1,], p[1,] + m$N1 * 5), col = "red")
lines3d(rbind(p[1,], p[1,] + m$N2 * 5), col = "purple")

###
wire3d(mesh3d(vertices = t(m$V), triangles = m$M))

# TODO: mesh with higher resolution near boundary;


###
p = matrix(c(1,5,1,3,2,3), nrow=2)
colnames(p) = c("x", "y", "z")
test.cylinder.line3d(1, p)


#############

### Examples: Tetrahedron

close3d()
open3d()
Th4.base(N = 32)


##################

### Torus: Cyclo-6

# Top NMR specialists have gathered in an emergency meeting
# to clarify which conformation is adopted: chair or boat?
close3d()
open3d()
test.cyclo.cylinder(19)

