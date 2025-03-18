

### Test Functions

# Functions for Various Tests


####################

### Helper Functions

# source("Chem3D.R")


####################

### Projection of Point on Plane

test.proj.plane = function(p, x, y, z,
		col = c("#9664D6", "#4464E8", "#FF6432", "#F2E248"),
		size = c(5,6), alpha = c(0.5, 0.325), scale = 1.5,
		test.math = TRUE, verbose = TRUE, tol = 1E-8) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	pp = proj.plane3d(p, x=x, y=y, z=z, verbose=verbose);
	p2 = rbind(p, pp$P);
	# Test:
	if(test.math) {
		pA = c(x[1], y[1], z[1]);
		d2 = dist.xyz(p2)^2 + dist.xyz(rbind(pp$P, pA))^2;
		d2 = d2 - dist.xyz(rbind(p, pA))^2;
		if(abs(d2) > tol) warning("Not right triangle: ", d2);
	}
	# Plot:
	polygon3d(x, y, z, col = col[1], alpha = alpha[1]);
	if(scale != 1) {
		xyz = expand.polygon3d(scale, x, y, z);
		polygon3d(xyz, col = col[4], alpha = alpha[2]);
	}
	lines3d(p2, size = size[1], col = col[2]);
	points3d(p2, size = size[2], col = col[3]);
	invisible(pp);
}


### Shifts

# p = PolygonL the first 3 points are used;
test.shift.poly = function(d, p, expand = 2, alpha = 0.75,
		alpha.expand = 0.25, shift.expand = 0.25) {
	pN  = eigen.plane(p[1:3, ]);
	len = nrow(p);
	# Expand base:
	pExp = expand.polygon3d(expand, p) - shift.expand*rep(pN$N, each = len);
	polygon3d(pExp, col = "#3296FF", alpha = alpha.expand);
	polygon3d(p, col = "blue", alpha = alpha);
	for(di in d) {
		pS = p + di*rep(pN$N, each = len);
		polygon3d(pS, col = "red", alpha = alpha);
	}
	invisible(pN);
}

test.rotate.ortho = function(p, pL, size.point = 6, col.point = "orange") {
	pR = rotate.ortho3d(p, pL);
	# Plot:
	lines3d(pL);
	# Lines: Projection & Orthogonal Rotation
	lines3d(rbind(p, pR$P), col = "blue");
	lines3d(rbind(pR$T1, pR$P), col = "red");
	points3d(p, col = col.point, size = size.point);
	points3d(pR$T1, col = col.point, size = size.point);
	invisible(pR);
}

test.rotate.point = function(p, pL, r = 1, n = 32,
		col.line = "blue", col.point = "red", size.point = 4) {
	lines3d(pL);
	N  = eigen.xy3D(p, pL);
	r  = sqrt(sum((p - N$P)^2));
	th = seq(0, 2*pi, length.out = n);
	for(phi in th) {
		pR = rotate.point3dN2(p, r=r, phi=phi, N=N);
		pR = rbind(N$P, pR);
		lines3d(pR, col = col.line);
		points3d(pR, size = size.point, col = col.point);
	}
	invisible(N);
}

test.eigen.lineAny = function(p, d = 1, both = TRUE, rev = FALSE,
		normalize = TRUE, col.point = "red", size.point = 6) {
	N  = eigen.lineAny(p, normalize=normalize);
	if(rev) N$N = - N$N;
	pL = rbind(p[1,], p[1,] + d * N$N);
	lines3d(p);
	lines3d(pL, col = "blue");
	points3d(p[1,], col = col.point, size = size.point);
	# Test:
	dd = sum((p[1,] + N$N - p[2,])^2, - N$N^2, - (p[1,] - p[2,])^2);
	if(abs(dd) > 1E-8) {
		cat("Error: Not Orthogonal! Diff = ", dd, "\n");
	}
	#
	if(both) {
		test.eigen.lineAny(p[2:1,], d=d, normalize=normalize,
			both = FALSE, rev = ! rev);
		# Note: same N is added to both ends;
		pL2 = p + rep(d*N$N, each=2);
		lines3d(pL2, col = "green");
	}
	invisible(N);
}


### Minimal Distance btw 2 Lines

# Minimal Rectangle
test.rect.simple = function(L, d = 1, y = d, t.extend = c(-1, 2),
		col = c("#000000", "#64D032", "#64D032", "#D032D0"),
		alpha = c(1, 0.5), size.points = 4) {
	if(size.points > 0) points3d(L, col = col[[1]], size = size.points)
	lines3d(L, col = col[[1]], alpha = alpha[[1]]);
	lines3d(L + c(0,0,y,y,0,0), col = col[[2]], alpha = alpha[[1]]);
	lines3d(L + c(0,0,y,y,d,d), col = col[[3]], alpha = alpha[[1]]);
	lines3d(L + c(0,0,0,0,d,d), col = col[[4]], alpha = alpha[[1]]);
	if(! is.null(t.extend)) {
		Lext = L;
		Lext[1,] = (1 - t.extend[1]) * L[1,] + t.extend[1] * L[2,];
		Lext[2,] = (1 - t.extend[2]) * L[1,] + t.extend[2] * L[2,];
		test.rect.simple(Lext, d=d, y=y, t.extend = NULL, size.points = 0,
			col=col, alpha = alpha[[2]]);
	}
}

# L  = Distance from L to L0;
# L0 = Base line;
test.lines.minDist.base = function(L, L0 = NULL, add = FALSE,
		col = c("#D03232"), verbose = TRUE, y.rect = 1, tol = 1E-8) {
	L2 = L;
	if(is.null(L0)) {
		L1 = rbind(c(0,0,0), c(2,0,0));
	} else L1 = L0;
	LEN = length(L2);
	if(LEN > 1 && length(col) == 1) col = rep(col, LEN);
	#
	if(add == FALSE) { close3d(); test.rect.simple(L1, y = y.rect); }
	if(LEN == 0) return();
	#
	d = sapply(seq_along(L2), function(id) {
		L2 = L2[[id]];
		d = dist.lines3d(L1, L2, tol=tol, verbose = verbose);
		lines3d(L2, col = col[id]);
		return(d);
	});
	#
	return(d);
}

### Special Cases:
# Special: Orthogonal Intersection
# L = Base line (can be NULL);
test.lines.minDist.Special1 = function(L = NULL, z = c(1,-1),
		col = c("#D03232"),
		add = FALSE, verbose = TRUE, tol = 1E-8) {
	pz = z[1]; pz2 = z[2];
	L2 = list(
		rbind(c(0,0,pz), c(0,0,pz2)),
		rbind(c(1,0,pz), c(1,0,pz2)),
		rbind(c(2,0,pz), c(2,0,pz2)),
		rbind(c(-1,0,pz), c(-1,0,pz2))
	);
	d = test.lines.minDist.base(L2, L, add=add, col=col,
		verbose=verbose, tol=tol);
	return(d);
}
# Special: Intersection, but NOT Orthogonal
test.lines.minDist.Special2 = function(L = NULL, z = c(1,-1),
		col = c("#D03232"),
		add = FALSE, verbose = TRUE, tol = 1E-8) {
	pz = z[1]; pz2 = z[2];
	L2 = list(
		rbind(c(1,0,pz), c(-1,0,pz2)),
		rbind(c(2,0,pz), c(1,0,pz2)),
		rbind(c(3,0,pz), c(-1,0,pz2)),
		rbind(c(3,0,pz), c(1,0,pz2))
	);
	d = test.lines.minDist.base(L2, L, add=add, col=col,
		verbose=verbose, tol=tol);
	return(d);
}
### Special: Ortho
test.lines.minDist.SpecialOrtho = function(L = NULL, z = c(1,-1),
		col = c("#D03232"),
		add = FALSE, verbose = TRUE, tol = 1E-8) {
	pz = z[1]; pz2 = z[2];
	L2 = list(
		rbind(c(0,1,pz), c(0,1,pz2)),
		rbind(c(1,1,pz), c(1,1,pz2)),
		rbind(c(2,1,pz), c(2,1,pz2)),
		rbind(c(2.5,1,pz), c(2.5,1,pz2))
	);
	d = test.lines.minDist.base(L2, L, add=add, col=col,
		verbose=verbose, tol=tol);
	return(d);
}
### General Case
test.lines.minDist.General = function(L = NULL, z = c(1,-1), y = 1,
		col = c("#D03232"),
		add = FALSE, verbose = TRUE, tol = 1E-8) {
	pz = z[1]; pz2 = z[2]; y1 = y[1];
	L2 = list(
		rbind(c(1,y1,pz), c(-3/2,y1,pz2)),
		rbind(c(3/2,y1,pz), c(-1,y1,pz2)),
		# Interesting case:
		rbind(c(2,y1,pz), c(-1/2,y1,pz2)),
		rbind(c(2,y1,pz), c(4,y1,pz2))
	);
	d = test.lines.minDist.base(L2, L, add=add, col=col,
		verbose=verbose, tol=tol, y.rect = y1);
	return(d);
}

##################
##################

### 3D Shapes

### Cylinder:
test.cylinder.line3d = function(r, x, y = NULL, z = NULL,
		col = "#8032B2", sides = 16, alpha = 0.5, lwd.line = 4) {
	# cyl = cylinder.line3d(r=r, x=x, y=y, z=z);
	if(lwd.line > 0) {
		lines3d(x, y, z, lwd = lwd.line, size = lwd.line);
	}
	if( ! is.null(y)) {
		x = cbind(x=x, y=y, z=z);
	}
	shade3d(cylinder3d(center = x, radius = r,
		col=col, sides=sides, alpha=alpha))
}


### Icosahedron

test.icosa.vset = function(m, r.ball = 0.75, lwd.line = 4, size.points = 4,
		col = "#8032B2", col.balls = "#F03244", col.fill = "#808080",
		alpha = c(0.75, 0.5)) {
	m = mesh.icosa.vset(v);
	test.icosa.mesh(m, r.ball=r.ball, lwd.line=lwd.line, size.points=size.points,
		col=col, col.balls=col.balls, col.fill=col.fill, alpha=alpha);
	invisible(m);
}
test.icosa.mesh = function(m, r.ball = 0.75, lwd.line = 4, size.points = 4,
		col = "#8032B2", col.balls = "#F03244", col.fill = "#808080",
		alpha = c(0.75, 0.5)) {
	# Sides:
	if( ! is.null(col.fill)) {
		shade3d(mesh3d(vertices = t(m$V), triangles = m$M),
			col = col.fill, alpha=alpha[2]);
	}
	# BackBone:
	wire3d(mesh3d(vertices = t(m$V), triangles = m$M), col=col, size = lwd.line);
	# Balls:
	if(r.ball == 0) {
		points3d(m$V, col = col.balls, size = size.points);
	} else {
		spheres3d(x = m$V, radius = r.ball,
			color = col.balls, alpha=alpha[1]);
	}
}

test.icosa = function(p, dH = NULL, dV = NULL, test.dist = TRUE,
		col.line = "blue", col.p = c("red", "orange"), size = 5) {
	ic = icosa(p, dH=dH, dV=dV, detailed = TRUE);
	plot.bb.icosa(ic, p=p, col.line=col.line, col.p=col.p, size=size);
	# Test:
	if(test.dist) {
		d = test.math.icosa(ic);
		cat("Sides:\n"); print(d);
	}
	invisible(ic)
}
test.math.icosa = function(v) {
	d = v[c(1:5, 1:5),] - v[c(6:10, 10, 6:9),];
	d = rbind(d, v[1:5, ] - rep(v[11,], each=5));
	d = rbind(d, v[6:10,] - rep(v[12,], each=5));
	d = apply(d, 1, function(d) sqrt(sum(d^2)));
	d = matrix(d, ncol = 5, byrow = TRUE);
	rownames(d) = c("Lat 1", "Lat 2", "Top", "Bottom");
	return(d);
}


test.cyclo.cylinder = function(n = 15, R = 3) {
	obj = cyclo.cylinder(n=n, R=R);
	lapply(obj, function(obj) lines3d(obj));
	invisible(obj);
}

