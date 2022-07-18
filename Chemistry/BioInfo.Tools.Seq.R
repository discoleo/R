
### Bioinformatics
###
### Leonard Mada


# - Tools to process data from IEDB;

##################

is.overlap = function(v1, v2, offset = 0) {
	len = min(length(v1), length(v2));
	for(npos in seq(len - offset)) {
		if(v1[[npos + offset]] != v2[[npos]]) return(FALSE);
	}
	return(TRUE)
}
overlap = function(v1, v2, strip=4) {
	for(nStart in seq(0, strip)) {
		if(is.overlap(v1, v2, offset = nStart)) {
			return(nStart);
		}
	}
	return(-1);
}
overlapAll = function(x, strip=4) {
	len = length(x);
	n0 = data.frame(id1 = numeric(0), id2 = numeric(0), offset = numeric(0));
	for(i in seq(len - 1)) {
		for(j in seq(i + 1, len)) {
			posOverlap = overlap(x[[i]], x[[j]]);
			if(posOverlap < 0) {
				posOverlap = overlap(x[[j]], x[[i]]);
				if(posOverlap < 0) next;
				posOverlap = - posOverlap;
			}
			n0 = rbind(n0, data.frame(i, j, posOverlap));
		}
	}
	names(n0) = c("id1", "id2", "offset")
	return(n0);
}

### Extract polypeptides which overlap
# - align = align the sequences;
extract = function(id, pos, data, ...) {
	UseMethod("extract");
}
extract.default = function(id, pos, data, align=TRUE, as.matrix=TRUE, debug=TRUE) {
	r  = extract.id(id, pos, data=data,
		align=align, as.matrix=as.matrix, debug=debug);
	return(r);
}
extract.row = function(id, pos, data, align=TRUE, as.matrix=TRUE, debug=TRUE) {
	id = pos$id1[[id]];
	r  = extract.id(id, pos, data=data,
		align=align, as.matrix=as.matrix, debug=debug);
	return(r);
}
extract.id = function(id, pos, data, align=TRUE, as.matrix=TRUE, debug=TRUE) {
	isID = pos$id1 == id;
	if( ! any(isID)) return("");
	#
	base = match(TRUE, isID);
	r = data[c(pos$id1[base], pos$id2[isID])];
	# Align
	if(align) {
		off = pos$offset[isID];
		off = c(0, off);
		if(debug) print(off);
		min = min(off);
		off = off - min;
		off = sapply(off, function(x) {
			if(x == 0) "" else paste0(rep("-", x), collapse="");
		})
		r = paste0(off, r);
	}
	if(as.matrix) r = matrix(r, ncol=1);
	return(r);
}

