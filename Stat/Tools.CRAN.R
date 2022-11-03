########################
###
### Leonard Mada
### [the one and only]
###
### Tools: Packages & CRAN
###
### draft v.0.2d


# this file:
# source("Tools.CRAN.R");


###############
### History ###
###############

### draft v.0.2d:
# - [refactor] moved code to file:
#   Tools.Format.String.R;
### draft v.0.2c:
# - splitDF: split long df into multiple columns;
### draft v.0.2b - v.0.2b-fix:
# - scroll arbitrary text;
### draft v.0.2a:
# - moved examples to separate file:
#   Tools.CRAN.Examples.R;
### draft v.0.1m - v.0.1n-fix2:
# - better word wrap;
# - more formatting options: cut(sep.h="-");
# - [fixed] length => nchar; [v.0.1n-fix]
# - full handling of multi-char sep.h; [v.0.1n-fix2]
### draft v.0.1l - v.0.1l-fix:
# - more examples;
# - [fixed] crash with only 1 record;
### draft v.0.1k:
# - [minor] bug fix;
# - added more examples & comments;


#######################

### Search CRAN
library(pkgsearch)


source("Tools.Format.String.R")

### Info about Packages
# - locally installed packages;

# Basic Info:
info.pkg = function(pkg=NULL, fields=c("Repository", "Description")) {
	if(is.null(pkg)) { pkg = installed.packages(fields=fields); }
	else {
		all.pkg = installed.packages();
		pkg = all.pkg[all.pkg[,1] %in% pkg, ];
	}
	p = pkg;
	p = as.data.frame(p);
	p = p[ , c("Package", "Version", "Built", fields, "Imports")];
	return(p);
}
# Imported packages:
imports.pkg = function(pkg=NULL, sort=TRUE) {
	p = info.pkg(pkg);
	### Imported packages
	imp = lapply(p$Imports, function(s) strsplit(s, "[,][ ]*"))
	imp = unlist(imp)
	imp = imp[ ! is.na(imp)]
	# Cleanup:
	imp = sub("[ \n\r\t]*+\\([-,. >=0-9\n\t\r]++\\) *+$", "", imp, perl=TRUE)
	imp = sub("^[ \n\r\t]++", "", imp, perl=TRUE);
	# Tabulate:
	tbl = as.data.frame(table(imp), stringsAsFactors=FALSE);
	names(tbl)[1] = "Name";
	if(sort) {
		id = order(tbl$Freq, decreasing=TRUE);
		tbl = tbl[id,];
	}
	return(tbl);
}

# Package Size:
size.f.pkg = function(path=NULL, exclude=NULL) {
	if(is.null(path)) path = R.home("library");
	xd = list.dirs(path = path, full.names = FALSE, recursive = FALSE);
	if( ! is.null(exclude)) {
		isExclude = grepl(exclude, xd);
		xd = xd[ ! isExclude];
	}
	size.f = function(p) {
		p = paste0(path, "/", p);
		sum(file.info(list.files(path=p, pattern=NULL,
			full.names = TRUE, all.files = TRUE, recursive = TRUE))$size);
	}
	sapply(xd, size.f);
}

size.pkg = function(path=NULL, sort=TRUE, file="Packages.Size.csv") {
	x = size.f.pkg(path=path);
	x = as.data.frame(x);
	names(x) = "Size"
	x$Name = rownames(x);
	# Order
	if(sort) {
		id = order(x$Size, decreasing=TRUE)
		x = x[id,];
	}
	if( ! is.null(file)) {
		if( ! is.character(file)) {
			print("Error: Size NOT written to file!");
		} else write.csv(x, file=file, row.names=FALSE);
	}
	return(x);
}

match.imports = function(pkg, x=NULL, quote=FALSE) {
	if(is.null(x)) x = info.pkg();
	if(quote) {
		pkg = paste0("\\Q", pkg, "\\E");
	}
	# TODO: Use word delimiters?
	# "(?<=^|[ \n\r\t],)"
	if(length(pkg) == 1) {
		isImport = grepl(pkg, x$Imports);
		return(x[isImport, ]);
	} else {
		# TODO: concept?
		rez = lapply(pkg, function(p) x[grepl(p, x$Imports), ]);
		return(rez);
	}
}

##################
### Formatting ###

# - moved to file:
#   Tools.Format.String.R;

### Split into columns
splitDF = function(x, ncol=2) {
	if(ncol <= 0) stop("Invalid number of columns!");
	if(ncol == 1) return(x);
	#
	len = nrow(x);
	if(len < 2) return(x);
	nr = len %/% ncol;
	#
	xr  = x[seq(1, len, by=ncol), ];
	idC = seq(2, ncol);
	for(id0 in idC) {
		tmp = x[seq(id0, len, by=ncol), ];
		nr0 = nrow(tmp);
		if(nr0 < nr) {
			if(nr0 == 0) { tmp = xr[nr, , drop=FALSE]; }
			else tmp = rbind(tmp, tmp[nr0, , drop=FALSE]);
		}
		xr = cbind(xr, tmp);
	}
	return(xr);
}


### Other / UI / Display

### Packages
scroll.pkg = function(pkg, start=1, len=15, w = c(12, 80, 16), iter=2,
		sep=" ", sep.h="-", print=TRUE) {
	if(len < 1) return();
	len  = len - 1;
	id = match(c("Package", "Description"), names(pkg));
	if(any(is.na(id))) {
		if( ! inherits(pkg, "pkg_search_result"))
			stop("Package info must contain both the name & description!");
		pkg = extract.pkg(pkg, type="Basic");
		id = c(1, 2);
	}
	pkg = cbind(pkg[, id], pkg[, - id]);
	# Column Lengths
	len.col = ncol(pkg); len.other = len.col - 2;
	w = if(len.other == 0) w[1:2]
		else if(length(w) == len.col) w
		else w[c(1,2, rep(w[3], len.other))]; # TODO: if(w[4])
	# Indent
	indent = c(list(c(" ", "   "), c("   ", "")), rep(list(""), len.other));
	# Entries
	if(start > nrow(pkg)) stop("No more entries!");
	nend = min(nrow(pkg), start + len);
	if(print) cat(c("Showing packages ", start, " to ", nend, "."), sep=c(rep("", 4), "\n"))
	cat.mlines(format.lines(pkg[seq(start, nend), ], w=w, indent=indent, iter=iter),
		sep=sep, sep.h=sep.h);
}

extract.pkg = function(x, type="Basic", print=TRUE) {
	# TODO: type
	# ex: Title, Maintainer, Date/Publication, downloads
	pkg = lapply(x$package_data, function(x)
			data.frame(
				Package = x$Package, Description = x$Description,
				Version = x$Version, Repository = x$Repository));
	pkg = do.call(rbind, pkg);
	nTotal = attr(x, "metadata")$total;
	if(print) cat(c("Found ", nTotal, " packages."), sep=c(rep("", 2), "\n"))
	return(pkg);
}

find.pkg = function(s, pkg=NULL, print=TRUE, perl=TRUE) {
	if(is.null(pkg)) pkg = info.pkg();
	isF = grepl(s, pkg$Description, perl=perl);
	pkg = pkg[isF, , drop=FALSE];
	if(print) print(paste0("Found ", nrow(pkg), " packages."))
	return(pkg);
}

###################
###################

###################
### Search CRAN ###
###################

# only simple expressions are possible:
# sep.h = row / horizontal separator;
searchCran = function(s, from=1, len=60, len.print=20, extend="*",
		sep=" ", sep.h="-") {
	if( ! is.null(extend)) s = paste0(s, extend);
	x = advanced_search(s, size=len, from=from);
	if(length(x$package_data) == 0) {
		cat("No packages found!", sep="\n");
	} else {
		scroll.pkg(x, len=len.print, sep=sep, sep.h=sep.h);
	}
	invisible(x)
}

