

### Compare 2 Directories
# - file name & file size;
# - TODO: option for hash;
match.dir = function(path1, path2, pattern = NULL) {
	x1 = list.files(path1, pattern=pattern, no.. = TRUE);
	cat("Finished Dir 1\n");
	x2 = list.files(path2, pattern=pattern, no.. = TRUE);
	cat("Finished Dir 2\n");
	f1 = file.size(paste0(path1, "/", x1));
	f2 = file.size(paste0(path2, "/", x2));
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
		copy = FALSE) {
	if(swap) { tmp = path1; path1 = path2; path2 = tmp; }
	FILES = match.dir(path1, path2, pattern=pattern);
	FILES = FILES[is.na(FILES$Match), ];
	if(copy) {
		# TODO
	}
	return(FILES);
}
