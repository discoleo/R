
# PDB Data Tools
# Exploratory Analysis


#####################

### Helper Packages & Functions

### Rpdb

# Note:
# - Function proj.line3d is available only in the version on GitHub!

# library(Rpdb)

# GitHub: 
# https://github.com/discoleo/Rpdb
#
# path = "..."
# devtools::document(path)
# devtools::load_all(path)

######################

### PDB Files
# setwd("...")

### PDB File: 4q02.pdb
# KRAS

x = read.pdb("4q02.pdb")

# GDP is a residue
unique(x$atoms$resname)

# Range:
isGDP = x$atoms$resname == "GDP"
bb = lapply(x$atoms[isGDP, c("x1", "x2", "x3")], range)
print(bb)

### Phosphate Groups
isP = isGDP & (x$atoms$elename %in% c("PA", "PB"))
bbP = x$atoms[isP, c("x1", "x2", "x3")];

### Distance from GDP-Channel
dist.atoms = function(p, data) {
	if(nrow(p) != 2) stop("p must be a line segment!");
	data = data$atoms[, c("x1", "x2", "x3")];
	pp = proj.line3d.matrix(data, p);
	dx = pp$x - data$x1;
	dy = pp$y - data$x2;
	dz = pp$z - data$x3;
	dd = sqrt(dx*dx + dy*dy + dz*dz);
	pp$d = dd;
	return(pp);
}

pp = dist.atoms(bbP, x)
# Distance from GDP-Channel
isWater = x$atoms$resname == "HOH"
isProt  = ! (isWater | isGDP)
plot(pp$d[isProt], pp$t[isProt])
points(pp$d[isWater], pp$t[isWater], col="#A032C096")
points(pp$d[isGDP], pp$t[isGDP], col="#FF0000A0")


### Plot 3D

# TODO: better visualization;

isWater = x$atoms$resname == "HOH"
isClose = pp$d <= 7 & (pp$t > -4 & pp$t < 4) & ! isGDP & ! isWater;
atoms = x$atoms[isClose, c("x1", "x2", "x3")];
# Phosphate P:
points3d(bbP, col = "blue", size = 9)
points3d(x$atoms[isGDP, c("x1", "x2", "x3")], size = 6, col = "red")
points3d(atoms)

