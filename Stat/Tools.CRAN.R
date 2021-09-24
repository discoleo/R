########################
###
### Leonard Mada
### [the one and only]
###
### Tools: Packages & CRAN
###
### draft v.0.1a



### Info about Packages

# - locally installed packages;

# Basic Info:
info.pkg = function(pkg=NULL) {
	if(is.null(pkg)) { pkg = installed.packages(); }
	else {
		all.pkg = installed.packages();
		pkg = all.pkg[all.pkg[,1] %in% pkg, ];
	}
	p = pkg;
	p = as.data.frame(p);
	p = p[ , c("Package", "Version", "Built", "Imports")];
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

size.pkg = function(path=NULL, sort=TRUE) {
	x = size.f.pkg(path=path);
	x = as.data.frame(x);
	names(x) = "Size"
	x$Name = rownames(x);
	# Order
	if(sort) {
		id = order(x$Size, decreasing=TRUE)
		x = x[id,];
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


###############
###############

### takes ages!
x = size.pkg();

# much faster: but NO size
p = info.pkg();
f = imports.pkg();


### Analyze data

# Size
head(x, 20)

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

