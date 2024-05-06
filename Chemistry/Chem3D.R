

library(rgl)

# rgl >= 1.3.1


###############

### Helper

dist.xyz = function(x, y, z) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	sqrt(diff(x)^2 + diff(y)^2 + diff(z)^2);
}

center.xyz = function(x, y, z) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	len = length(x);
	x = sum(x) / len;
	y = sum(y) / len;
	z = sum(z) / len;
	return(c(x,y,z));
}

expand.polygon3d = function(d, x, y, z, is.rel = TRUE) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	cx = mean(x); cy = mean(y); cz = mean(z);
	if( ! is.rel) {
		dx  = x - cx; dy = y - cy; dz = z - cz;
		div = sqrt(dx^2 + dy^2 + dz*2);
		d = d / div + 1; # 1 unit from centre by default;
	}
	x = d*x + (1-d)*cx;
	y = d*y + (1-d)*cy;
	z = d*z + (1-d)*cz;
	p = cbind(x=x, y=y, z=z);
	return(p);
}


##########################

### Orthogonal Projection/Rotation

# p = point which will be rotated by pi/2;
# (x,y,z) = coordinates of 2 points defining line;
### Out:
# T1,T2 = rotated point (2 solutions);
# P = projected point on line;

eigen.plane = function(x, y = NULL, z = NULL, normalize = TRUE) {
	if(is.null(y)) {
		p = x[1,]; x = x[-1,];
		y = x[,2]; z = x[,3]; x = x[,1];
	} else {
		p = c(x[1], y[1], z[1]);
		x = x[-1]; y = y[-1]; z = z[-1];
	}
	# pP = Projection of P on the line;
	# pP = proj.line3d(p, x, y, z);
	# Inline variant:
	dx = x[2] - x[1]; dx0 = p[1] - x[1];
	dy = y[2] - y[1]; dy0 = p[2] - y[1];
	dz = z[2] - z[1]; dz0 = p[3] - z[1];
	tt = (dx*dx0 + dy*dy0 + dz*dz0);
	tt = tt / (dx^2 + dy^2 + dz^2);
	t1 = 1 - tt;
	px = t1*x[1] + tt*x[2];
	py = t1*y[1] + tt*y[2];
	pz = t1*z[1] + tt*z[2];
	# Rotation:
	dx0 = p[1] - px; dx1 = x[1] - px;
	dy0 = p[2] - py; dy1 = y[1] - py;
	dz0 = p[3] - pz; dz1 = z[1] - pz;
	R = sqrt(dx0^2 + dy0^2 + dz0^2);
	dxy = dx0*dy1 - dx1*dy0;
	dxz = dx0*dz1 - dx1*dz0;
	dyz = dy0*dz1 - dy1*dz0;
	div = sqrt(dxy^2 + dxz^2 + dyz^2);
	# Note: both roots of "R" are valid;
	# TODO: check if dyz = 0, 
	dxT = dyz * R / div;
	dyT = - dxz * dxT / dyz;
	dzT =   dxy * dxT / dyz;
	# Note: - dT is also a valid solution;
	N = c(dxT, dyT, dzT);
	if(normalize) {
		# TODO: check;
		d = sqrt(dxT^2 + dyT^2 + dzT^2);
		N = N / d;
	}
	lst = list(N = N, P = c(px,py,pz));
	return(lst);
}
eigen.xy3D = function(p, x, y = NULL, z = NULL, normalize = TRUE) {
	if(is.null(y)) {
		xyz = x;
	} else {
		xyz = cbind(x,y,z);
	}
	xyz = rbind(p, xyz);
	N = eigen.plane(xyz, normalize=normalize);
	pP = N$P;
	Nx = pP - p;
	if(normalize) {
		Nx = Nx / sqrt(sum(Nx^2));
	}
	lst = list(Ny = N$N, Nx = Nx, P = pP);
	return(lst);
}
# Any Normal
eigen.lineAny = function(x, y = NULL, z = NULL, normalize = TRUE) {
	if(is.null(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	nx = x[2] - x[1]; ny = y[2] - y[1]; nz = z[2] - z[1];
	if(normalize) {
		d  = sqrt(nx*nx + ny*ny + nz*nz);
		nx = nx/d; ny = ny/d; nz = nz/d;
	}
	dn1 = nx - nz; dn2 = ny - nz;
	if(abs(dn1) < 1E-8) {
		# stop("Special Case: Div by 0!");
		cat("Special Case\n");
		# Normalized: if normalize == TRUE;
		fr = list(N = c(- nx, ny, nz));
		return(fr);
	}
	# Arbitrary Normal:
	if(normalize) {
		div = sqrt(2 * (dn1^2 + dn2^2 - dn1*dn2));
		fx  = dn1 / div;
	} else fx = dn1;
	# fx = dn1 / div;
	fy = - fx * dn2 / dn1;
	fz = - (fx + fy);
	fr = list(N = c(fx, fy, fz));
	return(fr);
}

### Rotate by angle phi
rotate.point3d = function(p, phi, x, y = NULL, z = NULL) {
	N = eigen.xy3D(p=p, x=x, y=y, z=z, normalize = FALSE);
	pP = N$P; # Projection of p on line;
	pR = pP + cos(phi) * r * N$Nx + sin(phi) * r * N$Ny;
	return(pR);
}
### Rotate by angle phi
# N = Normals given;
rotate.point3dN2 = function(p, r, phi, N) {
	pP = N$P; # Projection of p on line;
	pR = pP + cos(phi) * r * N$Nx + sin(phi) * r * N$Ny;
	return(pR);
}

rotate.ortho3d = function(p, x, y = NULL, z = NULL) {
	if(is.null(y)) {
		xyz = x;
	} else {
		xyz = cbind(x,y,z);
	}
	xyz = rbind(p, xyz);
	N  = eigen.plane(xyz, normalize = FALSE);
	pp = N$P; px = pp[1]; py = pp[2]; pz = pp[3];
	dT = N$N; dxT = dT[1]; dyT = dT[2]; dzT = dT[3];
	# Note: - dT is also a valid solution;
	xT1 = px + dxT; yT1 = py + dyT; zT1 = pz + dzT;
	xT2 = px - dxT; yT2 = py - dyT; zT2 = pz - dzT;
	lst = list(T1 = c(xT1,yT1,zT1), T2 = c(xT2,yT2,zT2),
		P = c(px,py,pz));
	return(lst);
}

### Tests

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
	if(both) {
		test.eigen.lineAny(p[2:1,], d=d, normalize=normalize,
			both = FALSE, rev = ! rev);
		# Note: same N is added to both ends;
		pL2 = p + rep(d*N$N, each=2);
		lines3d(pL2, col = "green");
		if(abs(dist.xyz(pL2) - dist.xyz(p)) > 1E-8) {
			cat("Error: Distances differ!");
		}
	}
	invisible(N);
}

### Examples:

###
p1 = c(1,5,7); p2 = c(3,5,7);
pL = matrix(c(-4,6,-5,8,-3,8), nrow=2)
test.rotate.ortho(p1, pL)
test.rotate.ortho(p2, pL, col.point = "#D68016")


### Shift Triangle
p3 = matrix(c(1,-4,6, 5,-5,8, 7,-3,8), nrow=3)
d = 3 * seq(4)

pN = eigen.plane(p3)

pExp = expand.polygon3d(2, p3) - 0.25*rep(pN$N, each = 3);
polygon3d(pExp, col = "#3296FF", alpha = 0.25)
polygon3d(p3, col = "blue", alpha = 0.75)
for(di in d) {
	pS = p3 + di*rep(pN$N, each = 3)
	polygon3d(pS, col = "red", alpha = 0.75)
}

### Rotate point
p = c(2,-3,17)
pL = matrix(c(-4,6,-5,8,-3,8), nrow=2)

test.rotate.point(p, pL)
test.rotate.point(c(1,1,-10), pL)


### Arbitrary Normal

###
pL = matrix(c(-4,6,-5,8,-3,8), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)

### Special Case:
pL = matrix(c(-4,4,-5,8,-5,3), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)

### Special Case:
pL = matrix(c(-4,4,-5,1,-5,3), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)

###
pL = matrix(c(-4,1,-5,3,-5,3), nrow=2)
d  = 3;
N = test.eigen.lineAny(pL, d=d, normalize = FALSE)
N = test.eigen.lineAny(pL, d=d)


################

################
### Cylinder ###

# Note:
# - rgl::cylinder3d already implements the cylinder,
#   including proper orientation based on centre-points;
# C = Center-point from where it starts;
cylinder.line3d = function(r, x, y = NULL, z = NULL) {
	if(is.null(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	N = c(x[2] - x[1], y[2] - y[1], z[2] - z[1]);
	d = sqrt(sum(N^2));
	N = N / d;
	lst = list(R = r, L = d, N = N, C = c(x[1],y[1],z[1]));
	return(lst);
}

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

p = matrix(c(1,5,1,3,2,3), nrow=2)
colnames(p) = c("x", "y", "z")
test.cylinder.line3d(1, p)


### Tetrahedron
Th4 = function(r = 1, center = c(0,0,0), phi = 0) {
	rT = r; r = r * 2 * sqrt(2) / 3;
	# Base Triangle:
	pB = if(phi == 0) {
			r2 = - r/2; r3 = - r2 * sqrt(3);
			cbind(c(r,r2,r2) + center[1], c(0,r3,-r3) + center[2], center[3]);
		} else {
			phi3 = c(0, 2*pi/3, 4*pi/3) + phi;
			x = r * cos(phi3) + center[1];
			y = r * sin(phi3) + center[2];
			cbind(x, y, center[3]);
		}
	pT = c(0, 0, r * sqrt(2)) + center; # Top
	pV = rbind(pB, pT);
	pC = c(0, 0, r/4 * sqrt(2)) + center; # Center
	return(list(V = pV, C = pC, R = rT));
}

plot.Th4 = function(v, lwd = 2, id.base = NULL,
		fill.plane = "#A0A0A064", N = 20, alpha = 0.2) {
	# Lines from Center to Vertices:
	pV = v$V;
	pC = v$C;
	for(i in seq(4)) {
		tmp = pV[i,];
		lines3d(c(tmp[1], pC[1]), c(tmp[2], pC[2]), c(tmp[3], pC[3]), lwd=lwd);
	}
	# Draw Base-plane:
	if( ! is.null(id.base)) {
		# TODO:
		id = if(id.base == 1) c(1,2,3)
			else if(id.base == 2) c(1,2,4)
			else if(id.base == 2) c(1,3,4) else c(2,3,4);
		pB = as.data.frame(pV[id, ]);
		names(pB) = c("x","y","z");
		### Transparency: set with alpha!
		polygon3d(pB, col = fill.plane, fill = TRUE, alpha = alpha);
		# fake transparency:
		trf = function(x, tt) {
			x1 = tt * x[1] + (1 - tt) * x[2];
			x2 = tt * x[1] + (1 - tt) * x[3];
			x  = c(x1, x2);
		}
		# Hashlines in Base-Plane:
		if(N > 0)
		for(tt in seq(0, N - 1)/N) {
			x = trf(pB$x, tt); y = trf(pB$y, tt);
			z = trf(pB$z, tt);
			lines3d(x, y, z, lwd = lwd, col = fill.plane);
		}
	}
}

### Basic Tetrahedron
# Note: transparency works only by specifying alpha!
Th4.base = function(lwd = 2, fill.plane = "#A0A0A064", N = 20, alpha = 0.2) {
	# Base Triangle:
	pB = list(c(1,0,0), c(-1/2, sqrt(3)/2, 0), c(-1/2, -sqrt(3)/2, 0))
	pT = c(0, 0, sqrt(2)); # Top
	pC = c(0, 0, 1/4 * sqrt(2)); # Center
	# Lines from Center to Vertices:
	for(i in seq(3)) {
		tmp = pB[[i]];
		lines3d(c(tmp[1], pC[1]), c(tmp[2], pC[2]), c(tmp[3], pC[3]), lwd=lwd);
	}
	tmp = pT;
	lines3d(c(tmp[1], pC[1]), c(tmp[2], pC[2]), c(tmp[3], pC[3]), lwd=lwd);
	# Draw Base-plane:
	if( ! is.null(fill.plane)) {
		pB = data.frame(do.call("rbind", pB));
		names(pB) = c("x","y","z");
		# Bug in old rgl:::triangulateSimple
		# polygon3d(pB, col = fill.plane);
		# polygon3d(rbind(pB, pB[1,]), ...);
		### Transparency: set with alpha!
		polygon3d(pB, col = fill.plane, fill = TRUE, alpha = alpha);
		# fake transparency:
		trf = function(x, tt) {
			x1 = tt * x[1] + (1 - tt) * x[2];
			x2 = tt * x[1] + (1 - tt) * x[3];
			x  = c(x1, x2);
		}
		# Hashlines in Base-Plane:
		if(N > 0)
		for(tt in seq(0, N - 1)/N) {
			x = trf(pB$x, tt); y = trf(pB$y, tt); z = c(0, 0);
			lines3d(x, y, z, lwd = lwd, col = fill.plane);
		}
	}
}

close3d()
open3d()
Th4.base(N = 32)


###
# xyz    = c(a,b,c)
# center = center of tetrahedron;
Th4 = function(xyz, center = c(0,0,0), r = 1) {
	# Eq: a*(x-x0) + b*(y-y0) + c*(z-z0) = 0
	# TODO
}


##################

### Torus: Cyclo-6
cycloCylinder = function(n, R = 3) {
	# Th:
	# th = pi - atan(4); cs = - 1/sqrt(17);
	# hb = dh / (1 + 1/sqrt(17));
	dz  = R * sin(pi/n) / (1 + 1/sqrt(17));
	# TODO:
	ds  = dz / 1.2;
	x1 = R * cos(2*pi*seq(n)/n); x1 = c(x1, x1[1]);
	y1 = R * sin(2*pi*seq(n)/n); y1 = c(y1, y1[1]);
	xB = 1.25 * R * cos(2*pi*seq(n)/n + pi/n);
	yB = 1.25 * R * sin(2*pi*seq(n)/n + pi/n);
	xS = 0.75 * R * cos(2*pi*seq(n)/n + pi/n);
	yS = 0.75 * R * sin(2*pi*seq(n)/n + pi/n);
	asV = function(x, xB, xS, id) {
		c(x[id], x[id], xB[id], x[id + 1], x[id + 1], xS[id], x[id]);
	}
l	tmp = lapply(seq(n), function(id) {
		lines3d(
			x = asV(x1, xB, xS, id=id),
			y = asV(y1, yB, yS, id=id),
			z = c(0, dz, dz + ds, dz, 0, - ds));
	})
	invisible()
}


# Top NMR specialists have gathered in an emergency
# to clarify which conformation is adopted: chair or boat?
close3d()
open3d()
cycloCylinder(19)


