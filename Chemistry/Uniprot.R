
### FASTA Tools


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

trim = function(x) {
	sub("^[ \t\n\r]+", "",
		sub("[ \t\n\r]+$", "", x));
}

# Regex:
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
extract.uptoreg = function(x, pattern, offset = c(0,0), perl = TRUE) {
	npos = regexpr(pattern, x, perl=perl);
	hasStr = npos > 0;
	sTk = rep("", length(x));
	if( ! any(hasStr)) return(sTk);
	posE = npos[hasStr];
	#
	sTk[hasStr] = substr(x[hasStr], offset[1], posE + offset[2]);
	return(sTk);
}

# more: extract protein abbreviation;
# Output:
# - Start = first line with AA;
# - End   = last line with AA;
# Side-Effects:
# - print.lines: Which lines of AA-seq to print; can be logical;
read.fasta = function(file, path = "", more = TRUE, print.lines = NULL) {
	fn = file;
	if(nchar(path) > 0) fn = paste0(path, "/", file);
	x  = readLines(fn);
	# Names:
	id = which(grepl("^>", x));
	sNms = x[id];
	# Seq ID:
	sPr  = extract.reg(sNms, "^>(?:sp|tr)[|][^| ]*+", offset = c(4,0));
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
		len = nchar(lst$IDP) + 6;
		tmp = substr(sNms, len, nchar(sNms));
		PrA = extract.reg(tmp, "^[^ \t]++", offset = c(0,0));
		# Protein Name:
		len = len + nchar(PrA) + 1;
		tmp = substr(sNms, len, nchar(sNms));
		PrN = extract.uptoreg(tmp, " OS=", offset = c(0,-1));
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
	if(more) {
		lst$PrAbbr = PrA;
		lst$PrN = PrN;
	}
	# Debug:
	if( ! is.null(print.lines)) {
		if(is.logical(print.lines) && length(print.lines) == 1) {
			print.lines = if(print.lines) 1:10 else numeric(0);
		}
		print(x[print.lines]);
	}
	return(lst)
}

read.fasta.pr = function(file, path = "", clean.data = TRUE) {
	fn = file;
	if(nchar(path) > 0) fn = paste0(path, "/", file);
	x  = readLines(fn);
	# Names:
	id = which(grepl("^>", x));
	sNms = x[id];
	# Seq ID:
	sPr = extract.reg(sNms, "^>(?:sp|tr)[|][^| ]*+", offset = c(4,0));
	lst = data.frame(IDP = sPr, Start = id + 1, End = c(id[-1] - 1, length(x)));
	lst$hasSeq = (lst$End - lst$Start >= 0);
	lst = lst[lst$hasSeq, ];
	if(nrow(lst) == 0) return(character(0));
	aa = sapply(seq(nrow(lst)), function(nr) {
		tmp = x[seq(lst$Start[nr], lst$End[nr])];
		tmp = paste0(tmp, collapse = "");
		return(tmp);
	});
	if(clean.data) aa = gsub("[ \n\r\t]+", "", aa);
	invisible(aa);
}


### Search-Tools

# data = data.frame with a Gene column;
which.gene = function(pattern, data) {
	which(grepl(pattern, data$Gene));
}


#################

### Quasi-BLAST

blast.quasi = function(x, data, print = TRUE, clean.data = TRUE) {
	len = nchar(x);
	if(len == 0) {
		warning("Nothing to search!");
		return(c(0,0));
	}
	if(length(data) > 1) {
		tmp = lapply(data, function(data) blast.quasi(x, data=data, print = FALSE));
		if(print) {
			lapply(seq_along(tmp), function(id) {
				tmp = tmp[[id]];
				if(length(tmp) > 1) {
					txt = paste0(tmp[-1], collapse=", ");
					txt = paste0(tmp[1], " => ", txt);
				} else txt = "";
				cat(id, ": ", txt, "\n", sep = "");
			});
		}
		return(invisible(tmp));
	}
	if(clean.data) {
		data = gsub("[ \r\n\t]+", "", data);
	}
	# to int:
	rX = charToRaw(x);
	rD = charToRaw(data);
	if(length(rD) < len) {
		warning("Not yet implemented!");
		return(c(0,1));
	}
	len1 = len - 1;
	LAST = length(rD) - len1;
	m = sapply(seq(1, LAST), function(id) {
		sum(rX == rD[seq(id, id + len1)]);
	});
	MAX = max(m);
	if(MAX == 0) return(c(0, 0));
	ids = which(m == MAX);
	return(c(MAX, ids));
}

