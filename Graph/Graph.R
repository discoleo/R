
#### Some "Graph" Utility Functions ####
#### v 0.1
####
#### Leonard Mada

### Limitations
# - does NOT generate (V, E)-structures;

### A.) Generate Random Points on a Grid
# e.g. for a traveling salesman problem

grid.gen = function(n, sizeGrid, offset=c(0,0)) {
	i = 0:(sizeGrid[1] * sizeGrid[2] - 1)
	locations.arr = sample(i, n)
	locations.df = data.frame(
		x = locations.arr %/% sizeGrid[1],
		y = locations.arr %% sizeGrid[1])
	locations.df$x = locations.df$x + offset[1]
	locations.df$y = locations.df$y + offset[2]
	return(locations.df)
}

### Examples

gr = grid.gen(10, c(10,10))
plot(gr$x, gr$y, pch=16)


gr = grid.gen(95, c(10,10))
plot(gr$x, gr$y, pch=16)



gr = grid.gen(25, c(10,10))
plot(gr$x, gr$y, pch=16)


### B.) Plot some special "Graphs"

### special cases of Cover

poly.geom = function(v, r = 1, rot=pi/2, addP=c(0,0), add=FALSE) {
	# v = number of peripheral vertices
	# rot = pi/2 or pi/v
	i = 1:v
	i = i - 1
	a = i * 2 * pi / v + rot
	x = r * cos(a)
	y = r * sin(a)
	if(add) {
		points(x, y, pch=16)
	} else {
		plot(x, y, pch=16)
	}
	lines(c(x, x[1]), c(y, y[1]))
	
	if(length(addP) == 2) {
		points(addP[1], addP[2], pch=16)
		i = i + 1
		for(p in i) {
			lines(c(addP[1],x[p]), c(addP[2], y[p]))
		}
	}
	
	return(data.frame(x=x, y=y))
}

### Examples

poly.geom(5)

### assymetric Center
poly.geom(5, addP=c(0.5, 0))

### Other sizes
poly.geom(4, rot=pi/4)

poly.geom(8, rot=pi/8)

### Overlap
poly.geom(7)
poly.geom(5, rot=-pi/2, r=0.5, add=T)

### Overlap
poly.geom(8, rot=pi/8)
poly.geom(8, rot=pi/8, r=0.5, add=T)
