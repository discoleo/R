
list.fasta = function(path, pattern="\\.fasta$") {
	list.files(path, pattern=pattern, no.. = TRUE,
		include.dirs = FALSE);
}

# Identifiers in the FileName
extract.fname = function(x, pattern = "\\.Chr[^.]++(?=\\.)", offset = 5) {
	if(length(offset) == 1) offset = c(offset, 0);
	sTk = extract.reg(x, pattern=pattern, offset=offset);
	return(sTk);
}
extract.reg = function(x, pattern, offset = c(0,0), perl = TRUE) {
	npos = regexpr(pattern, x, perl=perl);
	hasStr = npos > 0;
	sTk = rep("", length(x));
	if( ! any(hasStr)) return(sTk);
	posS = npos[hasStr];
	posE = attr(npos, "match.length")[hasStr];
	posE = posE + posS - 1;
	#
	sTk[hasStr] = substr(x[hasStr], posS + offset[1], posE + offset[2]);
	return(sTk);
}

# more: extract protein abbreviation;
read.fasta = function(file, path = "", more = TRUE, print.lines = NULL) {
	fn = file;
	if(nchar(path) > 0) fn = paste0(path, "/", file)
	x  = readLines(fn);
	id = which(grepl("^>", x));
	sNms = x[id];
	# Seq ID:
	sPr  = extract.reg(sNms, "^>sp[|][^| ]*+", offset = c(4,0));
	lst  = data.frame(IDP = sPr, Start = id + 1, End = c(id[-1] - 1, length(x)));
	lst$hasSeq = (lst$End - lst$Start >= 0);
	# Seq-Length
	nAA = nchar(x); # NL are already removed;
	nAA[id] = 0;
	nAA = cumsum(nAA);
	nAA = diff(c(0, nAA[lst$End]));
	lst$Len = nAA;
	# Protein abbreviation
	if(more) {
		len = nchar(lst$IDP);
		tmp = substr(sNms, len + 6, nchar(sNms));
		PrA = extract.reg(tmp, "^[^ \t]++", offset = c(0,0));
	}
	# Gene-Name
	sG = extract.reg(sNms, "(?<=[ \t|])GN=[^| \t]++", offset = c(3,0));
	idNoG = which(nchar(sG) == 0);
	isNoG = length(idNoG) > 0;
	if(more && isNoG) {
		sG[idNoG] = PrA[idNoG];
	} else if(isNoG) {
		len = nchar(lst$IDP[idNoG]);
		tmp = sNms[idNoG];
		# LEN = 5 is 1 LESS;
		tmp = substr(tmp, len + 6, nchar(tmp));
		sG[idNoG] = extract.reg(tmp, "^[^ \t]++", offset = c(0,0));
	}
	lst$Gene = sG;
	if(more) lst$PrAbbr = PrA;
	# Debug:
	if( ! is.null(print.lines)) print(x[print.lines]);
	return(lst)
}


### Search-Tools

# data = data.frame with a Gene column;
which.gene = function(pattern, data) {
	which(grepl(pattern, data$Gene));
}

