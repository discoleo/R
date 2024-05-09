
### Chem 3D
# - Exploratory / Experimental code;
# - Intended to be used within Rpdb:
#   https://github.com/discoleo/Rpdb


library(rgl)

# rgl >= 1.3.1

# Math functions: are independent of rgl;
# Plot & Test functions: require rgl;


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

# Both Normals: but starting in arbitrary direction;
# Out: P = projection of p0 (identical to p1);
eigen.lineAnyN2 = function(x, y = NULL, z = NULL, normalize = TRUE) {
	if(is.null(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	# Normals:
	N1 = eigen.lineAny(x=x, y=y, z=z, normalize=normalize);
	N1 = N1$N;
	p1 = c(x[1], y[1], z[1]);
	p0 = p1 + N1;
	p2 = c(x[2], y = y[2], z = z[2]);
	N2 = eigen.projPoint(p0, p1, x = p2);
	N2 = N2$N;
	#
	lst = list(N1 = N1, N2 = N2, P0 = p0, P = p1);
	return(lst);
}


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
	tt = tt / (dx*dx + dy*dy + dz*dz);
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
# Projection Point: given
# - p0 projects on pP;
eigen.projPoint = function(p0, pP, x, y = NULL, z = NULL,
		normalize = TRUE, verbose = TRUE) {
	if(is.null(y)) {
		if(inherits(x, "matrix")) x = x[1,];
		y = x[2]; z = x[3]; x = x[1];
	}
	# Rotation:
	dx0 = p0[1] - pP[1]; dx1 = x[1] - pP[1];
	dy0 = p0[2] - pP[2]; dy1 = y[1] - pP[2];
	dz0 = p0[3] - pP[3]; dz1 = z[1] - pP[3];
	R = sqrt(dx0^2 + dy0^2 + dz0^2);
	dxy = dx0*dy1 - dx1*dy0;
	dxz = dx0*dz1 - dx1*dz0;
	dyz = dy0*dz1 - dy1*dz0;
	div = sqrt(dxy^2 + dxz^2 + dyz^2);
	# Note: both roots of "R" are valid;
	# TODO: check if dyz = 0; (already solved ?)
	NN = R / div;
	dxT =   dyz * NN;
	dyT = - dxz * NN;
	dzT =   dxy * NN;
	# Note: - dT is also a valid solution;
	N = c(dxT, dyT, dzT);
	if(normalize) {
		# TODO: may be already normalized; (check)
		d = sqrt(dxT^2 + dyT^2 + dzT^2);
		if(verbose) cat("Normalize by: ", d, "\n");
		N = N / d;
	}
	lst = list(N = N);
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
	if(abs(dn2) < 1E-8) {
		# stop("Special Case: Div by 0!");
		cat("Special Case\n");
		# Normalized: if normalize == TRUE;
		# TODO: check normalization;
		fr = list(N = c(0, - ny, nz));
		return(fr);
	}
	# Arbitrary Normal:
	if(normalize) {
		div = sqrt(2 * (dn1^2 + dn2^2 - dn1*dn2));
		fx  = dn2 / div;
	} else fx = dn2;
	# fx = dn1 / div;
	fy = - fx * dn1 / dn2;
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

# Vertexes & Mesh:
mesh.cylinder = function(r, x, y = NULL, z = NULL, nL = 12, nR = 16,
		type = c("Full", "Alternating")) {
	V = mesh.vertex.cylinder(r=r, x=x, y=y, z=z, nL=nL, nR=nR);
	M = mesh.cylinder.vset(V, nL=nL, nR=nR, type=type);
	invisible(M);
}

mesh.cylinder.vset = function(V, nL = 12, nR = 16,
		type = c("Full", "Alternating")) {
	# TODO: Mesh;
	len = attr(V$V, "length");
	isOdd = len[1] > len[2];
	type = match.arg(type);
	# Note: Circle has nR+1 points;
	idS2 = len[1] * (nR + 1);
	idM = matrix(0, nrow = 3, ncol = 2*nR*len[2]);
	for(id in seq(len[2])) {
		id2 = 2*(id - 1);
		idS = id2 * nR + 1;
		idE = idS + nR - 1;
		idSlice = idS:idE;
		# Note: Circle has nR+1 points;
		idV = (id - 1)*(nR + 1) + 1;
		idV = idV:(idV + nR - 1);
		idM[, idSlice] = rbind(idV, idV + 1, idV + idS2);
		idSlice = idSlice + nR;
		idM[, idSlice] = rbind(idV + 1, idV + idS2 + 1, idV + idS2);
	}
	if(type == "Full") {
		last = len[2];
		# if(isOdd) last = last - 1;
		idM2 = matrix(0, nrow = 3, ncol = 2*nR*last);
		idS2 = idS2 - nR - 1; # previous Circle;
		for(id in seq(last)) {
			id2 = 2*(id - 1);
			idS = id2 * nR + 1;
			idE = idS + nR - 1;
			idSlice = idS:idE;
			# Note: Circle has nR+1 points;
			idV = (id - 1)*(nR + 1) + nR + 2;
			idV = idV:(idV + nR - 1);
			# print("St1"); print(idSlice)
			idM2[, idSlice] = rbind(idV + 1, idV + idS2 + 1, idV + idS2);
			idSlice = idSlice + nR;
			# print("St2"); print(idSlice)
			idM2[, idSlice] = rbind(idV, idV + 1, idV + idS2);
		}
		idM = cbind(idM, idM2);
	}
	V$M = idM;
	return(V);
}

# Vertexes:
mesh.vertex.cylinder = function(r, x, y = NULL, z = NULL, nL = 12, nR = 16) {
	if( ! is.null(y)) {
		xyz = cbind(x, y, z);
	} else xyz = x;
	# Centers
	tt = seq(0, 1, length.out = nL + 1);
	tt = cbind(1 - tt, tt);
	pC = tt %*% xyz;
	# Normals:
	N = eigen.lineAnyN2(xyz, normalize = TRUE);
	N1 = N$N1;
	N2 = N$N2;
	# Start:
	phi1 = seq(0, 2*pi, length.out = nR + 1);
	cyl1 = lapply(seq(1, nL + 1, by = 2), function(id) {
		cc = pC[id, ];
		dx = r*(cos(phi1) * N1[1] + sin(phi1) * N2[1]) + cc[1];
		dy = r*(cos(phi1) * N1[2] + sin(phi1) * N2[2]) + cc[2];
		dz = r*(cos(phi1) * N1[3] + sin(phi1) * N2[3]) + cc[3];
		cbind(dx, dy, dz);
	});
	len1 = length(cyl1);
	cyl1 = do.call(rbind, cyl1);
	#
	phi2 = phi1 + pi / nR;
	cyl2 = lapply(seq(2, nL + 1, by = 2), function(id) {
		cc = pC[id, ];
		dx = r*(cos(phi2) * N1[1] + sin(phi2) * N2[1]) + cc[1];
		dy = r*(cos(phi2) * N1[2] + sin(phi2) * N2[2]) + cc[2];
		dz = r*(cos(phi2) * N1[3] + sin(phi2) * N2[3]) + cc[3];
		cbind(dx, dy, dz);
	});
	len2 = length(cyl2);
	cyl2 = do.call(rbind, cyl2);
	V = rbind(cyl1, cyl2);
	attr(V, "length") = c(len1, len2);
	return(list(V = V, N1 = N1, N2 = N2));
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


###############

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

###
# xyz    = c(a,b,c)
# center = center of tetrahedron;
Th4 = function(xyz, center = c(0,0,0), r = 1) {
	# Eq: a*(x-x0) + b*(y-y0) + c*(z-z0) = 0
	# TODO
}


##################

### Torus: Cyclo-6
# Poly-cyclic compound around a cylinder;
cyclo.cylinder = function(n, R = 3) {
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
	obj = lapply(seq(n), function(id) {
		cbind(
			x = asV(x1, xB, xS, id=id),
			y = asV(y1, yB, yS, id=id),
			z = c(0, dz, dz + ds, dz, 0, - ds, 0));
	});
	invisible(obj);
}

test.cyclo.cylinder = function(n = 15, R = 3) {
	obj = cyclo.cylinder(n=n, R=R);
	lapply(obj, function(obj) lines3d(obj));
	invisible(obj);
}

