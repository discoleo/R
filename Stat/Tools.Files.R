
### List Sub-Directories
list.dirs = function(path, pattern = NULL, all.cols = FALSE, verbose = FALSE) {
	x = list.files(path, pattern=pattern, no.. = TRUE,
		include.dirs = TRUE);
	if(verbose) cat("Finished Dir!\n");
	fd = file.info(paste0(path, "/", x), extra_cols = FALSE);
	fd = fd[fd$isdir, ];
	if(! all.cols) {
		id = match(c("mode", "ctime", "atime"), names(fd));
		id = id[! is.na(id)];
		fd = fd[, - id];
	}
	return(fd);
}

### Compare 2 Directories
# - file name & file size;
# - TODO: option for hash;
match.dir = function(path1, path2, sub.path = NULL, pattern = NULL, verbose = TRUE) {
	if(! is.null(sub.path)) {
		path1 = paste0(path1, "/", sub.path);
		path2 = paste0(path2, "/", sub.path);
	}
	# Note: NO easy way to exclude directories!
	x1 = list.files(path1, pattern=pattern, no.. = TRUE,
		include.dirs = FALSE);
	if(verbose) cat("Finished Dir 1\n");
	x2 = list.files(path2, pattern=pattern, no.. = TRUE,
		include.dirs = FALSE);
	if(verbose) cat("Finished Dir 2\n");
	# Filter Dirs:
	# f1 = file.size(paste0(path1, "/", x1));
	f1 = file.info(paste0(path1, "/", x1), extra_cols = FALSE);
	isFile = ! f1$isdir;
	f1 = f1$size;
	if(verbose) {
		cat("Excluded ", length(isFile) - sum(isFile), " dirs.\n");
	}
	x1 = x1[isFile];
	f1 = f1[isFile];
	f2 = file.size(paste0(path2, "/", x2));
	# Match:
	s1 = paste0(f1, "\t", x1);
	s2 = paste0(f2, "\t", x2);
	idMatch = match(s1, s2);
	#
	f1 = data.frame(Name = x1, Size = f1);
	# f2 = data.frame(Name = x2, Size = f2);
	# dd = setdiff(f1, f2);
	f1$Match = idMatch;
	return(f1);
}

# Difference between 2 dirs
diff.dir = function(path1, path2, pattern = NULL, swap = FALSE,
		copy = FALSE, verbose = TRUE) {
	if(swap) { tmp = path1; path1 = path2; path2 = tmp; }
	FILES = match.dir(path1, path2, pattern=pattern, verbose=verbose);
	FILES = FILES[is.na(FILES$Match), ];
	if(copy) {
		fN = paste0(path1, "/", FILES$Name);
		# Manual overwrite!
		r1 = file.copy(fN, path2, overwrite = FALSE, recursive = FALSE,
			copy.mode = FALSE, copy.date = TRUE);
		fFAIL = FILES[ ! r1, ];
		if(nrow(fFAIL) > 0) {
			cat("Failed to copy:\n");
			print(fFAIL);
		}
	}
	return(FILES);
}


diff.dirs = function(x, path1, path2, verbose = c("Top", "Full", "None")) {
	verbose  = match.arg(verbose);
	verb.all = verbose == "Full";
	verb.dir = verbose != "None";
	xn = nchar(x)
	isNotRoot = xn > 0;
	p1[isNotRoot] = paste0(path1, "/", x[isNotRoot]);
	p2[isNotRoot] = paste0(path2, "/", x[isNotRoot]);
	#
	lst.diff = lapply(seq_along(x), function(id) {
		dd = diff.dir(p1[id], p2[id], verbose = verb.all);
		if(verb.dir) {
			cat("Dir: ", x[id], "\n");
			if(nrow(dd) > 0) { print(dd); }
			else {
				cat("   EQUAL!\n");
			}
		}
		return(dd);
	});
	invisible(lst.diff);
}

