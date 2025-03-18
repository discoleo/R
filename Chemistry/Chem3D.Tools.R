

####################

### Helper Functions


### Distances
dist.xyz = function(x, y, z) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	sqrt(diff(x)^2 + diff(y)^2 + diff(z)^2);
}
dist.l3d = function(x) {
	# TODO: process each segment of line;
	d = x[1:3] - x[4:6];
	sqrt(sum(d^2));
}
dist.p3d = function(p1, p2) {
	d = p2 - p1;
	sqrt(sum(d^2));
}


### Center
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
# x is converted to Matrix: 3 x nCol;
center.p3d = function(x) {
	dim(x) = c(3, length(x) / 3);
	apply(x, 1, mean);
}


######################

### Viewport Rotations
# (non-destructive rotations)
# - Apply rotation on the rgl viewport;

rotate.xz = function(alpha = pi/2) {
	m = matrix(c(cos(alpha),0,sin(alpha),0,
		0,1,0,0,
		-sin(alpha),0,cos(alpha),0,
		0,0,0,1), nrow = 4);
	old = par3d("userMatrix");
	par3d(userMatrix = m %*% old);
	highlevel();
	invisible(old);
}
rotate.xy = function(alpha = pi/6) {
	m = matrix(c(cos(alpha),sin(alpha),0,0,
		-sin(alpha),cos(alpha),0,0,
		0,0,1,0, 0,0,0,1), nrow = 4);
	old = par3d("userMatrix");
	par3d(userMatrix = m %*% old);
	highlevel();
	invisible(old);
}


###################

###################
### Projections ###

# - Orthogonal projections;

### Projection of Point on Line
# p = point which will be projected on (x,y,z) line;
proj.line3d = function(p, x, y, z) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	dx = x[2] - x[1]; dx0 = p[1] - x[1];
	dy = y[2] - y[1]; dy0 = p[2] - y[1]; # p[2] = yP;
	dz = z[2] - z[1]; dz0 = p[3] - z[1];
	tt = (dx*dx0 + dy*dy0 + dz*dz0);
	tt = tt / (dx^2 + dy^2 + dz^2);
	t1 = 1 - tt;
	px = t1*x[1] + tt*x[2];
	py = t1*y[1] + tt*y[2];
	pz = t1*z[1] + tt*z[2];
	# TODO: compute d or not?
	# d = sqrt(sum((p - c(px, py, pz))^2));
	lst = list(P = c(px, py, pz), t = tt);
	return(lst);
}

### Projection of Point on Plane
# p = 3D Point;
# (x, y, z) = 3 points defining the plane;
proj.plane3d = function(p, x, y, z, verbose = FALSE, tol = 1E-8) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	len = c(length(x), length(y), length(z));
	stopifnot(len == 3);
	#
	id = c(1,2);
	p1 = proj.line3d(p, x[id], y[id], z[id]);
	t1 = p1$t;
	p1 = p1$P;
	pC = c(x[3], y[3], z[3]);
	pT = proj.line3d(pC, x[id], y[id], z[id]);
	tT = pT$t; dT = 1 - tT;
	if(abs(dT) < tol) {
		if(verbose) print("Case: Div 0!");
		p2 = pC * t1 + c(x[1], y[1], z[1]) * (1-t1);
	} else if(abs(t1 - 1) < tol) {
		if(verbose) print("Case: Special!");
		p2 = pC + c(x[2], y[2], z[2]) - pT$P;
	} else {
		tt = (1 - t1) / dT;
		p2 = pC * tt + c(x[2], y[2], z[2]) * (1-tt);
	}
	pp = proj.line3d(p, rbind(p1, p2));
	return(list(P = pp$P, pP12 = p1));
}


### Decompositions into Normals

eigen.lineN2 = function(p, x, y = NULL, z = NULL, normalize = TRUE) {
	N = eigen.xy3D(p, x=x, y=y, z=z, normalize=normalize);
	N = list(N1 = N$Nx, N2 = N$Ny, P = N$P);
	return(N);
}

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


### Generate 2 Normals
# (x,y,z) = coordinates of 3 points defining plane;
# p = 1st point in (x,y,z);
#     will be rotated by pi/2 around the remaining line;
### Out:
# N = normal to plane (x,y,z);
# p[out] = projected point on line;
# Note:
# - T1,T2 = rotated point (2 valid solutions);
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
### Project point p on line (x,y,z)
# Returns:
# - Normals & Projection Point;
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
# - pP = Projection point (known);
# - p0 = point which projects on pP;
# - x, y, z = forms with pP a line on which p0 is projected;
# Note: avoids recomputing pP if it is already known;
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
	# TODO: check if dyz == 0; (already solved ?)
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
# x, y, z = given line;
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


#################
### Rotations ###

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


#####################
### Intersections ###

dist.lines3d = function(xyz, xyz0, tol = 1E-8, verbose = TRUE) {
	p1 = proj.line3d(xyz[1,], xyz0);
	d1 = sqrt(sum((xyz[1,] - p1$P)^2));
	if(d1 <= tol) {
		# Lines intersect
		if(verbose) cat("Lines intersect at P1.\n");
		return(0);
	}
	p2 = proj.line3d(xyz[2,], xyz0);
	d2 = sqrt(sum((xyz[2,] - p2$P)^2));
	if(d2 <= tol) {
		# Lines intersect
		if(verbose) cat("Lines intersect at P2.\n");
		return(0);
	}
	# Special Cases:
	d = sqrt(sum((p1$P - p2$P)^2));
	# Orthogonal lines:
	if(d <= tol) {
		# Both points project on same point;
		# Project point back:
		pO = proj.line3d(p1$P, xyz);
		dO = sqrt(sum((p1$P - pO$P)^2));
		if(dO <= tol) {
			# Lines intersect
			if(verbose) cat("Lines intersect midway.\n");
			return(0);
		}
		if(verbose) cat("Lines are orthogonal.\n");
		return(dO);
	}
	# Equidistant Projections
	if(abs(d1 - d2) <= tol) {
		pM = (p1$P + p2$P) / 2;
		pO = proj.line3d(pM, xyz);
		dO = sqrt(sum((pM - pO$P)^2));
		if(dO <= tol) {
			# Lines intersect
			if(verbose) cat("Lines intersect midway.\n");
			return(0);
		}
		# Note: lines may be parallel,
		# but is not essential for distance;
		if(verbose) cat("Special Lines but NO intersection.\n");
		return(dO);
	}
	# General Case: Minimize d
	dd = xyz[2,] - p2$P;
	dy = xyz[1,] - p1$P;
	ds = dy - dd;
	div = sum(ds*ds);
	tt  = - sum(dd * ds) / div;
	tt1 = 1 - tt;
	pL1 = tt * xyz[1,] + tt1 * xyz[2,];
	pL2 = tt * p1$P + tt1 * p2$P;
	d = sqrt(sum((pL1 - pL2)^2));
	if(d <= tol) {
		if(verbose) cat("Intersection: General case.\n");
		d = 0;
	} else if(verbose) cat("General case.\n");
	return(d);
}

# Shortest Distance between 2 lines:
# - includes points delimiting that line segment;
# - NOT to be confused with proj.line3d;
proj.lines3d = function(xyz, xyz0, tol = 1E-8) {
	p1 = proj.line3d(xyz[1,], xyz0);
	d = sqrt(sum((xyz[1,] - p$P)^2));
	if(d <= tol) {
		# Lines intersect
		lst = list(P1 = xyz[1,], p2 = p$P, d = 0,
			doIntersect = TRUE, isParallel = FALSE);
		return(lst);
	}
	stop("TODO");
	# TODO
}

