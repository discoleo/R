########################
###
### Leonard Mada
### [the one and only]
###
### Tools: Packages & CRAN
###
### draft v.0.1d-ext



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
size.f.pkg = function(path=NULL) {
	if(is.null(path)) path = R.home("library");
	xd = list.dirs(path = path, full.names = FALSE, recursive = FALSE);
	size.f = function(p) {
		p = paste0(path, "/", p);
		sum(file.info(list.files(path=p, pattern=".",
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
	# "(<?=^|[ \n\r\t],)"
	if(length(pkg) == 1) {
		isImport = grepl(pkg, x$Imports);
		return(x[isImport, ]);
	} else {
		# TODO: concept?
		rez = lapply(pkg, function(p) x[grepl(p, x$Imports), ]);
		return(rez);
	}
}

split.line = function(s, w=80, nL=NULL, indent = c("   ", "")) {
	if(is.null(nL)) nL = 1 + ((nchar(s) - 1) %/% w);
	if(nL == 1) return(paste0(indent[1], s));
	nMax = nchar(s); n0 = 1; n = w; dn = w %/% 2;
	if(length(indent) == 1) indent = c(indent, "");
	#
	s2 = character(nL);
	for(id in seq(nL - 1)) {
		DO_NEXT = FALSE;
		indent0 = if(id == 1) indent[1] else indent[2];
		for(npos in seq(n, n - dn)) {
			if(substr(s, npos, npos) %in% c(" ", "\n", ",", "-", ")")) {
				s2[id] = paste0(indent0, substr(s, n0, npos));
				n0 = npos + 1; n = min(nMax, n0 + w);
				DO_NEXT = TRUE;
				break;
			}
		}
		if(DO_NEXT) next;
		s2[id] = paste0(indent0, substr(s, n0, n));
		n0 = min(nMax, n + 1); n = min(nMax, n0 + w);
	}
	s2[nL] = paste0(indent[2], substr(s, n0, nMax));
	return(s2);
}
format.lines = function(x, w=80, justify="left", NL.rm=TRUE) {
	if(NL.rm) {
		for(nc in seq(ncol(x))) {
			x[, nc] = gsub("[ \t]*+\n++[ \t]*+", " ", x[, nc], perl=TRUE);
		}
	}
	# Detect Long Lines
	n  = sapply(seq(ncol(x)), function(nc) nchar(x[,nc]));
	nL = as.matrix(n);
	nL = 1 + ((nL - 1) %/% w);
	maxL = apply(nL, 1, max, na.rm=TRUE);
	nL[is.na(nL)] = 1;
	csm = cumsum(c(1, maxL));
	txt = matrix("", nrow=tail(csm, 1), ncol=ncol(x));
	for(nc in seq(ncol(x))) {
		for(nr in seq(nrow(x))) {
			nL0 = nL[nr, nc];
			txt[seq(csm[nr], length.out=nL0), nc] = split.line(x[nr, nc], w=w, nL=nL0);
		}
	}
	return(apply(txt, 2, format, justify=justify));
}


###############
###############

### Package Size
# Note: takes ages!
if(FALSE) {
	# !! setwd(...); !!
	x = size.pkg();
}
if(FALSE) {
	system.time({
		x = size.pkg(file=NULL);
	})
	# elapsed time: 509 s !!!
	# 512 Packages; 1.64 GB;
}

x = read.csv("Packages.Size.csv")


### Imports
# - much faster: but NO size;
p = info.pkg();
f = imports.pkg();


#####################
### Data Analysis ###

### Size
head(x, 20)

### Description
# - packages which do NOT import any other package;
format.lines(p[is.na(p$Imports), ][1:20, -6])


# - some are NOT Bioconductor packages;
# - TODO: filter by biocViews?
p[is.na(p$Repository), 1:4]


# No imports
table(is.na(p$Imports))
# Most imported
head(f, 20)

# imported only once:
f$Name[f$Freq == 1]


match.imports("hunspell", p)
match.imports("labeling", p)
match.imports("rpart.plot", p)

# Concept?
match.imports(c("pROC", "ROCR"), p)

