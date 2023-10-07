

### More 3D Shapes


####################

### Helper Functions

library(rgl)

dist.xyz = function(xyz, xyz0 = c(0,0,0)) {
	sqrt(sum((xyz - xyz0)^2))
}
runif.sphere.fib = function(n = 20) {
	# golden angle in radians
	phi = pi * (sqrt(5) - 1);
	
	y = 1 - (2 * seq(n)) / n;  # y goes from 1 to -1
	radius = sqrt(1 - y * y);  # radius at y
	theta = phi * seq(n);  # golden angle increment
	x = cos(theta) * radius;
	z = sin(theta) * radius;
	pp = cbind(x, y, z);
    return(pp);
}

shift.3d = function(xyz, d, center = c(0,0,0)) {
	xyz = xyz - rep(center, each = nrow(xyz));
	L   = dist.xyz(xyz[1,]);
	tt  = 1 + d / L;
	xyz = xyz * tt + rep(center, each = nrow(xyz));
	return(xyz);
}


###
pp = runif.sphere.fib()

close3d()
open3d()
spheres3d(c(0,0,0), r = 1, alpha = 0.2)
spheres3d(pp, r = 1/8)
spheres3d(shift.3d(pp, d = 1/4), r = 1/8)


###
n = 32
pp1 = runif.sphere.fib(n=n)
pp2 = shift.3d(pp1, d = 2)

close3d()
open3d()
spheres3d(c(0,0,0), r = 1, alpha = 0.2)
tmp = sapply(seq(nrow(pp1)), function(id) {
	cyl = cylinder3d(rbind(pp1[id,], pp2[id,]), radius = c(1/2, 1/8), sides = 12,
		alpha=0.9)
	# shade3d(addNormals(subdivision3d(cyl)))
	shade3d(cyl)
})
spheres3d(pp2, r = 1/8)

