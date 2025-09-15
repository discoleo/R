########################
###
### Leonard Mada
### [the one and only]
###
### Code Tools
###
### draft v.0.2k


### Tools to Process Formulas & Expressions


###############
### History ###
###############

### draft v.0.2a - v.0.2b, v.0.2k:
# - moved Code tools to this file: Tools.Code.R;
# - moved Parser Code to file: Tools.Code.Parser.R;
### draft v.0.1g - v.0.1g-improve:
# - cut.code() into code blocks;
# - small improvements; [v.0.1g-improve]
### draft v.0.1f - v.0.1f-refactor:
# - extract code tokens from R code;
# - [refactored] uniform result;


### TODO
# - Explore: remindR;

# this file:
# source("Tools.Code.R")

### Minimalistic Parser
source("Tools.Code.Parser.R")


#########################

### GREP Utilities

# dir = sub-directory to search for files;
grepdir = function(x, dir, path = NULL, perl = TRUE,
		f.pattern = "\\.[Rr]$", warn.f = TRUE) {
	if(! is.null(path)) {
		dir = paste0(path, "/", dir);
	}
	lstFiles = list.files(dir, pattern = f.pattern, full.names = TRUE);
	if(length(lstFiles) == 0) {
		cat("NO R files!\n");
		return(invisible());
	}
	# GREP:
	lst = lapply(lstFiles, function(ff) {
		tmp = readLines(ff, warn = warn.f);
		idL = which(grepl(x, tmp, perl=perl));
		if(length(idL) == 0) return();
		return(list(File = ff, Lines = idL));
	})
	hasCode = sapply(lst, function(x) ! is.null(x));
	lst = lst[hasCode];
	return(lst);
}

#########################

# List files of type ".R"
list.filesR = function(path, pattern = NULL, full.names = FALSE,
		case.sens = FALSE, perl = TRUE, file.ext = "\\.[Rr]$") {
	x   = list.files(path, full.names = full.names);
	isR = grep(file.ext, x);
	x   = x[isR];
	if(is.null(pattern)) return(x);
	if(length(x) == 0) return(character(0));
	#
	tmp = if(full.names) basename(x) else x;
	if(! case.sens) pattern = paste0("(?i)", pattern);
	isP = grepl(pattern, tmp, perl=perl);
	x = x[isP];
	return(x);
}
# List files in R-directory
# - based on package pkgGraphR;
list.filesInR = function(path, pattern = NULL, case.sens = FALSE, perl = TRUE,
		dirR = FALSE) {
	stopifnot(is.character(path), length(path) == 1);
	x = normalizePath(path);
	if(dir.exists(x)) {
		isR = grep("^R$", list.dirs(x, full.names = FALSE));
		if(any(isR)) {
			x = paste0(x, "/R");
		} else if(dirR) {
			# Force R-dir:
			stop("The path does not contain an R directory!");
		}
		fR = list.filesR(x, full.names = TRUE,
			pattern=pattern, case.sens=case.sens, perl=perl);
		if(length(fR) == 0) {
			stop("Path does not contain any .R files!");
		}
	} else if(file.exists(x)) {
        fR = grep(".[Rr]$", x, value = TRUE);
		if(length(fR) == 0){
			stop("'path' is a file, but not an .R file.");
		}
	} else {
		stop("Path must be a valid file or directory");
	}
	return(fR);
}

### List Functions
# - lists all functions defined in the source files;
# path    = path to a file or directory;
# pattern = filter list of files;
list.functions = function(path, pattern = NULL,
		exclude = "[]$:<[`-]", verbose = FALSE, ...) {
	filesR = list.filesInR(path, pattern=pattern, ...);
	return(list.functions.files(filesR, exclude=exclude, verbose=verbose));
}
list.functions.files = function(files,
		exclude = "[]$:<[`-]", verbose = FALSE, ...) {
	allFx = sapply(files, function(x) {
        pR  = parse(x, keep.source = TRUE);
        pRD = utils::getParseData(pR) |>
			dplyr::filter(token != "COMMENT", terminal != FALSE);
        return(findFunNames(pRD));
    }, simplify = FALSE, USE.NAMES = TRUE);
	hasR  = sapply(allFx, function(x) nrow(x) > 0);
	allFx = allFx[hasR];
	# allFx = do.call(rbind, allFx);
	if(! is.null(exclude)) {
		if(verbose) {
			tmp = lapply(allFx, function(x) x$name);
			tmp = unlist(tmp);
			fxExcluded = grepl(exclude, tmp);
			tmp = tmp[fxExcluded];
			cat("Excluded functions:");
			if(length(tmp) == 0) { cat(" 0;\n"); }
			else { cat("\n"); print(tmp); }
		}
		allFx = lapply(allFx, function(x) {
			isFxName = ! grepl(exclude, x$name);
			return(x[isFxName, ]);
		});
	}
	class(allFx) = c("listFx", "list");
	return(allFx)
}

# Extract calls to functions
fun.calls = function(x, inline.exclude = FALSE) {
	stopifnot(is.character(x), length(x) == 1);
	### Functions:
	funs = list.functions(x);
	# very rudimentary mechanism;
	if(inline.exclude) {
		funs = lapply(funs, function(x) {
			x[x$nrParent == 0, ];
		});
	}
	allFuns = lapply(funs, function(x) x$name);
	names(allFuns) = NULL;
	allFuns = unique(unlist(allFuns));
	#
	nms = names(funs);
	res = lapply(nms, function(z) {
		findFunCalls(z, funs, fun.names = allFuns);
	});
	res = unlist(res);
	res = list(Calls = res, Fun = allFuns);
	class(res) = c("listCalls", "list");
	return(res);
}

duplicated.listFx = function(x, fromLast = FALSE) {
	fx = lapply(x, function(x) x$name);
	fx = unlist(fx);
	isDupl = duplicated(fx, fromLast = fromLast);
	fd = fx[isDupl];
	return(fd);
}


### Function Definitions: Info
### Out:
# nrParent:
#  0 = top-tier Function;
#  n = row of data.frame (with the parent function);
# TopParent: row of data.frame where the top parent is defined;
# - see also Issues:
#   https://github.com/discoleo/R/issues/1
findFunNames = function(x) {
	len = nrow(x) - 2;
	EMPTY = function() {
		return(data.frame(line = numeric(0), name = character(0)));
	}
	if(len <= 0) return(EMPTY());
	# res = lapply(seq(len), function(i) {
	#	if(x$token[i] == "SYMBOL" && (
	#		x$token[i+1] == "LEFT_ASSIGN" ||
	#		x$token[i+1] == "EQ_ASSIGN") &&
	#		x$token[i+2] == "FUNCTION") {
	#			return(data.frame(line = x$line1[i], name = x$text[i]));
	#		}
	# });
	# res = do.call(rbind, res);
	idF = which(x$token == "FUNCTION");
	idF = idF[idF >= 3];
	if(length(idF) == 0) return(EMPTY());
	tmp = x$token[idF - 1];
	idF = idF[tmp == "LEFT_ASSIGN" | tmp == "EQ_ASSIGN"];
	if(length(idF) == 0) return(EMPTY());
	tmp = x$token[idF - 2];
	idF = idF[tmp == "SYMBOL"];
	if(length(idF) == 0) return(EMPTY());
	idF = idF - 2;
	res = data.frame(line = x$line1[idF], name = x$text[idF],
		parent = x$parent[idF+2], Assign = x$parent[idF+1]);
	# Formals
	res$idFormS = x[idF + 3, "id"];
	formE = sapply(seq(nrow(res)), function(id) {
		idF = res$parent[id];
		ids = which(x$parent == idF);
		tmp = x[ids, c("id", "token", "text")];
		tmp = tmp$id[tmp$text == ")"];
		if(length(tmp) == 0) return(NA); # some Error
		return(tail(tmp,1));
	});
	res$idFormE = formE;
	# Function Body:
	fbSE = sapply(res$idFormE, function(id) {
		if(is.na(id)) return(c(NA, NA));
		idB = which(x$id == id) + 1;
		if(x$token[idB] != "'{'") {
			return(c(x$id[idB], NA));
		}
		idP = x$parent[idB];
		ids = which(x$parent == idP);
		return(x$id[c(idB, tail(ids,1))]);
	});
	res$idBS = fbSE[1,];
	res$idBE = fbSE[2,];
	res$hasBody = ! is.na(res$idBE);
	# End of Body: special case
	ids = which(! res$hasBody & ! is.na(res$idBS));
	idBE = sapply(ids, function(id) {
		idA = res$Assign[id];
		idE = x$id[x$parent < idA];
		return(tail(idE, 1));
	});
	res$idBE[ids] = unlist(idBE); # some Error?
	# Inline Functions:
	res = parent.default(res);
	# Top Parent
	res$TopParent = which.parent.top(res$nrParent);
	return(res);
}
# x = list of data.frames with parse-info;
parent.list = function(x) {
	for(idFile in seq_along(x)) {
		x[[idFile]]$nrParent =
			which.parent(x[[idFile]]);
	}
}
# x = data.frame with parse-info;
parent.default = function(x) {
	x$nrParent = which.parent(x);
	return(x);
}
which.parent = function(x) {
	len = nrow(x);
	if(len == 0) return(numeric(0));
	lenFx  = x$idBE - x$idBS;
	parent = sapply(seq(len), function(id) {
		idP  = x$parent[id];
		idPF = which(x$idBS < idP & x$idBE > idP);
		if(length(idPF) == 0) return(0);
		if(length(idPF) == 1) return(idPF);
		lenP = lenFx[idPF];
		idP0 = which(lenP == min(lenP));
		idPF = idPF[idP0];
		return(idPF);
	});
	return(parent);
}
# Top Parent (Root)
which.parent.top = function(x) {
	topParent = x;
	notR   = topParent > 0;
	isRoot = ! notR; # has reached Root;
	idHasP = which(notR);
	idRoot = which(isRoot); # as valid ids;
	topParent[idRoot] = idRoot; # all are valid ids;
	iter = length(topParent);
	while(length(idHasP) > 0 && iter > 0) {
		iter = iter - 1;
		idP1 = topParent[idHasP];
		topParent[idHasP] = topParent[idP1];
		idT1 = which(isRoot[idP1]);
		if(length(idT1) > 0) {
			tmp_id = idHasP[idT1];
			idHasP = idHasP[- idT1];
			isRoot[tmp_id] = TRUE;
			# topParent[tmp_id] = topParent[idP1[idT1]];
		}
	}
	if(length(idHasP) > 0) {
		warning("Data contains cycles!");
		print(idHasP);
	}
	return(topParent);
}

# list.fun = list with info about function definitions;
findFunCalls = function(x, list.fun, fun.names) {
	if(! inherits(x, "parseData")) {
		pR = parse(x, keep.source = TRUE);
		pR = utils::getParseData(pR);
	} else pR = x;
	isCall = pR$token == "SYMBOL_FUNCTION_CALL";
	pRD = pR[isCall, ];
	# TODO: do.call
	isDoCall = pRD$text == "do.call";
	pRDDC = pRD[isDoCall, ];
	pRD = pRD[! isDoCall, ];
	#
	if(is.null(fun.names)) {
		# TODO: extract from list.fun;
	}
	id  = match(pRD$text, fun.names);
	isC = ! is.na(id);
	pRD = pRD[isC, ];
	if(nrow(pRD) == 0) return(character(0));
	# Parent Fx:
	# TODO:
	# Note: line1 is NOT robust!
	fD = list.fun[[x]];
	ii = findInterval(pRD$line1, fD$line);
	funParent = fD$name[ii];
	funCalls  = setNames(pRD$text, funParent);
	return(funCalls);
}

### Graph
buildPackageGraph = function(x, unique.edges = TRUE, only.connected = FALSE,
		rep.node = NULL, unique = TRUE) 
{
	stopifnot(is.logical(unique.edges), length(unique.edges) == 1);
	stopifnot(inherits(x, "listCalls"));
	res = x$Calls; allFuns = x$Fun;
	# Unique Calls
	if(unique) {
		isDupl = duplicated(cbind(res, names(res)));
		res = res[! isDupl];
	}
	# Duplicate node:
	# - e.g. an extensively connected node;
	# TODO: do NOT replicate if called in the same function;
	if(! is.null(rep.node)) {
		for(sNode in rep.node) {
			ids = which(res == sNode);
			if(length(ids) <= 1) next;
			repFx = paste0(sNode, ".", seq(length(ids)));
			res[ids] = repFx;
			allFuns  = c(allFuns, repFx);
		}
	}
	### Graph
    edges <- dplyr::filter(data.frame(from = unname(res), to = names(res)), 
        !is.na(from), !is.na(to))
    res <- list(nodes = allFuns, edges = edges)
    if (unique.edges) {
        res$edges <- dplyr::distinct(res$edges)
    }
    if (only.connected) {
        drop <- setdiff(res$nodes, unique(c(res$edges$from, res$edges$to)))
        res$nodes <- setdiff(res$nodes, drop)
    }
    return(res)
}
graph.viz = function(x, sep = "\uB7", width = 800, height = 8000) {
	# sep can be one of, e.g.:
	# c("\uA0", "\uA4", "\uA6", "\uB0", "\uB7")
	graph = x;
	graph$nodes = gsub("\\.", sep, graph$nodes);
	graph$edges$from = gsub("\\.", sep, graph$edges$from);
	graph$edges$to   = gsub("\\.", sep, graph$edges$to);
	#
	gg = sapply(1:nrow(graph$edges), function(i) {
		paste0(graph$edges$from[i], " -> ", graph$edges$to[i])
	});
	edges = paste0(gg, collapse = "\n");
	nodes = paste0(graph$nodes, collapse = ";");
	graph.init = paste0("digraph{", "\n", "graph[rankdir = LR]", "\n",
		"node[shape = box, style = rounded, fontname = Helvetica]", "\n");
	fullGraph = paste0(graph.init, nodes, edges, "\n", "}");
	if(is.null(height)) {
		tmp = DiagrammeR::grViz(fullGraph);
	} else {
		tmp = DiagrammeR::grViz(fullGraph, width=width, height=height);
	}
	return(tmp);
}


################
### Packages ###

### List All Packages
ls.pkg = function(pkg=NULL, more.fields=FALSE, fields=c("Repository", "Description", "Imports")) {
	if(more.fields) {
		all.pkg = installed.packages(fields=fields);
	} else {
		all.pkg = installed.packages();
		fields = NULL;
	}
	if(is.null(pkg)) {
		pkg = all.pkg;
	} else {
		pkg = all.pkg[all.pkg[,1] %in% pkg, , drop = FALSE];
	}
	if(nrow(pkg) == 0) {
		warning("No packages!");
		return(data.frame("Package" = character(0)));
	}
	p = pkg;
	p = as.data.frame(p);
	p = p[ , c("Package", "Version", "Built", fields)];
	rownames(p) = seq(nrow(p));
	return(p);
}

### List all Functions in a package
# Note: only functions exported in the namespace;
ls.fun = function(pkg, exclude.C = TRUE) {
	pkg = as.character(match.call()[[2]]);
	nms = ls(getNamespace(pkg));
	if(exclude.C) {
		isC = sapply(nms, is.call.C, pkg=pkg);
		nms = nms[ ! isC];
	}
	return(nms);
}

### Args
args = function(name, default = TRUE, verbose = TRUE) {
	if(default) {
		fn = match.call()[[2]];
		if(is.function.generic(fn)) {
			fn = paste0(as.character(fn), ".default");
			name = fn;
			if(verbose) cat(fn, "\n");
		}
	}
	.Internal(args(name));
}

is.function.generic = function(name) {
	# TODO: is.function.generic();
	# - this version is a little bit ugly;
	# - S4: if(isGeneric(name));
	length(do.call(.S3methods, list(name))) > 0;
}

class.fun = function(fn, pkg) {
	fn = as.symbol(fn); pkg = as.symbol(pkg);
	fn = list(substitute(pkg ::: fn));
	# deparse
	s = do.call(class, fn);
	return(s)
}

# much faster!
is.call.C = function(FUN, pkg) {
	isC = "NativeSymbolInfo" %in% class.fun(FUN, pkg);
}
is.code.RC = function(x) {
	# Note:
	# - detailed parsing of the whole code
	#   is not efficient;
	npos = parse.RC(x);
	return(any(npos$Type == 50));
}

##############

# Note: requires Tools.Code.Parser.R
# source(Tools.Code.Parser.R)


# Parse code of specific function
parse.fun = function(fn, pkg, type=99) {
	# deparse
	s = deparse.fun(fn, pkg);
	npos = parse.simple(s, all.tokens = (type > 2));
	return(npos);
}
deparse.fun = function(fn, pkg, collapse="\n", width.cutoff=160) {
	fn = as.symbol(fn); pkg = as.symbol(pkg);
	fn = list(substitute(pkg ::: fn), width.cutoff=width.cutoff);
	# deparse
	s = paste0(do.call(deparse, fn), collapse=collapse);
	return(s)
}

# Extract Calls
# TODO:
# nms = ...
# tmp = lapply(nms, function(fn) { print(fn); parse.calls(deparse.fun(fn, "stats")); } )
parse.calls = function(x) {
	# C Calls
	C.df = parse.pC(x);
	if(any(C.df$Type == 50)) {
		cat("C Call!\n");
		return();
	}
	#
	x = parse(text = x);
	parse.f2 = function(e) {
		if(is.name(e[[1]])) {
			tmp = as.character(e[[1]]);
			if(e == "<-" || e == "=") {
				if(is.call(e[[3]])) {
					if(as.character(e[[3]][[1]]) == "function") {
						# FUN Definition
						return(as.character(e[[2]]));
					} else {
						# TODO: Call to a specific function;
					}
				}
			}
		}
		return(""); # TODO
	}
	lapply(x, parse.f2);
}
parse.RC = function(x, eol="\n") {
	return(parse.simple(x, eol=eol, all.tokens = "c"));
}


### Processing

# cut code into disjoint code blocks
cut.code = function(npos, last=0) {
	# Start
	nS = npos$nS;
	nSNext = npos$nE;
	nSNext = nSNext[nSNext != max(nSNext)];
	nS = sort(unique(c(1, nS, nSNext + 1))); # 1 = 1st Token;
	# End
	nE = npos$nE;
	nENext = npos$nS;
	nENext = nENext[nENext != min(nENext)];
	nE = sort(unique(c(nE, nENext - 1)));
	# 1st Token of code
	npos1 = npos$nS[[1]] - 1;
	if(npos1 > 0) nE = c(npos1, nE);
	# last Token of code
	if(last > 0 & (maxE <- max(nE)) < last) {
		nS = c(nS, maxE + 1);
		nE = c(nE, last);
	}
	#
	tk.df = data.frame(nS = nS, nE = nE);
	# tk.df$Type = 0;
	# TODO: classify all tokens;
	tk.df = merge(tk.df, npos[, c("nS", "Type")], by="nS", all.x=TRUE);
	#
	class(tk.df) = c("code", class(tk.df));
	return(tk.df);
}
# extract strings delimited by npos[, 1:2]
extract.str = function(s, npos, strip=FALSE, trim.regex=NULL, format.sp=TRUE) {
	if(nrow(npos) == 0) return(character(0));
	strip.FUN = if(strip) {
			function(nr) {
				if(npos$nS[[nr]] + 1 < npos$nE[[nr]]) {
					nStart = npos$nS[[nr]] + 1;
					nEnd = npos$nE[[nr]] - 1; # TODO: Error with malformed string
					return(substr(s, nStart, nEnd));
				} else {
					return("");
				}
			}
		} else function(nr) substr(s, npos$nS[[nr]], npos$nE[[nr]]);
	sR = sapply(seq(nrow(npos)), strip.FUN);
	if( ! is.null(trim.regex)) sR = gsub(trim.regex, "", sR, perl=TRUE);
	if(format.sp) sR = gsub("(?<=,|<-)(?=[^ \t\n\r])", " ", sR, perl=TRUE);
	return(sR);
}
extract.str.fun = function(fn, pkg, type=1, strip=TRUE, trim.regex=NULL) {
	# TYPE: 1 = all strings; 11 = Comments; 90 = All code tokens; 99 = All tokens
	s = deparse.fun(fn, pkg);
	npos = parse.simple(s, all.tokens = (type > 2));
	# TODO: more advanced filtering;
	isType = if(type == 90) (npos$Type > 20)
		else if(type == 99) rep(TRUE, nrow(npos)) else (npos$Type == type);
	extract.str(s, npos[isType, ], strip=strip, trim.regex=trim.regex);
}
extract.str.pkg = function(pkg, type=1, exclude.z = TRUE, strip=TRUE) {
	nms = ls(getNamespace(pkg));
	l = lapply(nms, function(fn) extract.str.fun(fn, pkg, type=type, strip=strip));
	if(exclude.z) {
		hasStr = sapply(l, function(s) length(s) >= 1);
		nms = nms[hasStr];
		l = l[hasStr];
	}
	names(l) = nms;
	return(l);
}
# Read from file
read.code = function(file, all.tokens=FALSE) {
	s = readLines(file);
	s = paste0(s, collapse="\n");
	npos = parse.simple(s, all.tokens = all.tokens);
	npos = cut.code(npos, last=nchar(s));
	s = extract.str(s, npos, strip=FALSE);
	attr(s, "type") = npos$Type;
	return(s);
}

### Formatting

# preserve "relevant" NLs;
regex.trim = function() {
	paste0("(?<=^|\n)[ \n\t\r]++|[ \t]++(?=\n)|(?<=[ \t])[ \t]++",
		"|[ \t\n\r]++$|(?<=[+*-/] )[ \t\n\r]++",
		"|(?<=,)[\n\r]++|[\n\r]++(?=,| ,)");
}
format.code = function(s, check.code=TRUE) {
	type = attr(s, "type");
	if(check.code && is.null(type)) stop("Input is NOT code!");
	# ",..." => ", ..."
	isCode = is.na(type);
	s[isCode] = gsub(",(?![ \t\n\r])", ", ", s[isCode], perl=TRUE);
	s[isCode] = gsub("[ \t]+\n", "\n", s[isCode]);
	# Comments:
	isComment = ( ! is.na(type)) & (type == 2);
	s[isComment] = gsub("^#(?![# ])", "# ", s[isComment], perl=TRUE);
	return(s);
}


####################
####################

### Extract/Summary Args

summary.args = function(e) {
	rs = lapply(e, function(e) {
		if(length(e) == 0) return(data.frame(type="NULL"));
		if(length(e) == 1) {
			if(is.symbol(e)) {
				return(data.frame(type="No Default"));
			}
		}
		if(is.call(e)) {
			if(e[[1]] == "c") return(data.frame(type="Default val"));
			if(e[[1]] == "list") return(data.frame(type="Default val"));
			if(class(e) == "if") return(data.frame(type="Code"));
			return(data.frame(type="Call"));
		}
		if(is.character(e) && e == "") return(data.frame(type="Default val: Empty"));
		return(data.frame(type="Default val"));
	});
	nm = names(e);
	rs = do.call(rbind, rs);
	rs$Name = nm; rs = rs[, c(2,1)];
	rs$type[rs$Name == "..."] = "hasDot";
	rownames(rs) = NULL;
	return(rs);
}

# Arguments for all functions in a package:
summary.all.args = function(nm) {
	f = ls(getNamespace(nm))
	r = lapply(seq_along(f), function(id) {
		fn = f[id];
		# if(substr(fn,1,1) %in% c("[", "_"))
		# DONE also: "<-"
		fn = paste0("\"", fn, "\"");
		e = parse(text=paste0("formals(", nm, ":::", fn, ")"));
		e = eval(e);
		if(is.null(e)) return(data.frame(Name=NA, type=NA, FUN=f[id]));
		a = summary.args(e);
		a$FUN = f[id];
		return(a);
	})
	
	do.call(rbind, r);
}


################

### Examples:
# - moved to file: Tools.Code.Tests.R;
