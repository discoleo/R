########################
###
### Leonard Mada
### [the one and only]
###
### Pubmed
### Text Mining Tools
###
### draft v.0.1q


### Text Mining Tools
# - Identify position of parenthesis in Pubmed abstracts;
# - Extract content between parenthesis;



######################

### Basic Tools

extract.regex = function(x, pattern, gr=0, perl=TRUE, simplify=TRUE, verbose=TRUE) {
	if(inherits(x, "data.frame")) stop("x should be an array!");
	r = regexec(pattern, x, perl=perl);
	if(verbose) cat("Finished Regex.\nStarting extraction.\n");
	gr = gr + 1;
	s = lapply(seq(length(x)), function(id) {
		tmp = r[[id]];
		if(tmp[1] == -1) return("");
		nStart = tmp[gr];
		substr(x[id], nStart, nStart - 1 + attr(tmp, "match.length")[gr]);
	});
	if(simplify) s = unlist(s);
	return(s)
}

# strips the Section Title
stripSection = function(x, len=25) {
	patt = paste0("(?i)^[^:]{1,", len, "}+\\:\\s*+");
	npos = regexpr(patt, x, perl=TRUE);
	npos = attr(npos, "match.length");
	hasSection = (npos > 0);
	x[hasSection] = substring(x[hasSection], 1 + npos[hasSection]);
	return(x);
}


### Sentence / Parenthesis Parser

# Parser for Parenthesis:
# sub.tokens: all tokens are currently identified/extracted;
parseParenth = function(x, sub.tokens=TRUE, warn=FALSE) {
	if(length(x) > 1) {
		return(lapply(x, parseParenth, warn=warn));
	}
	len = nchar(x);
	### Type:
	# 23 = "{", 24 = "(", 25 = "[";
	### Error:
	# 1 = No open tag;
	# 23, 24, 25 = wrong match;
	tk.df = data.frame(nS=integer(0), nE=integer(0), id=integer(0), Type=integer(0),
			Nested=logical(0), hasNested=logical(0), Err=numeric(0));
	stackType = integer(0);
	stackPos  = integer(0);
	#
	npos = 1;
	while(npos <= len) {
		s = substr(x, npos, npos);
		# Brackets
		type = 0;
		if(s == "{") {
			type = 23;
		} else if(s == "(") {
			type = 24;
		} else if(s == "[") {
			type = 25;
		} else if(s == "\uFF08") {
			# Special Parenthesis
			type = 24;
		}
		if(type > 0) {
			nr = nrow(tk.df) + 1;
			isNested = if(length(stackType) > 0) TRUE else FALSE;
			if(isNested) {
				prevStack = length(stackPos);
				tk.df$hasNested[stackPos[prevStack]] = TRUE;
			}
			stackType = c(stackType, type);
			stackPos  = c(stackPos, nr);
			tk.df[nr, ] = data.frame(npos, NA, nr, type, isNested, FALSE, 0);
			npos = npos + 1;
			next;
		}
		#
		if(s == "}") {
			type = 23;
		} else if(s == ")") {
			type = 24;
		} else if(s == "]") {
			type = 25;
		} else if(s == "\uFF09") {
			# Special Parenthesis
			type = 24;
		}
		if(type > 0) {
			# check unmatched Token
			idPrev = length(stackPos);
			if(idPrev == 0) {
				nr = nrow(tk.df) + 1;
				tk.df[nr, ] = data.frame(npos, npos, nr, type, FALSE, FALSE, 1);
			} else {
				posPrev = stackPos[idPrev];
				if(stackType[idPrev] != type) {
					# Token mismatch
					tk.df[posPrev, "Err"] = type;
					# save also position of closing parenth;
					tk.df[posPrev, "nE"]  = npos;
				} else {
					tk.df[posPrev, "nE"] = npos;
				}
				if(idPrev == 1) {
					stackPos = integer(0);
					stackType = integer(0);
				} else {
					stackPos  = stackPos[ - idPrev];
					stackType = stackType[ - idPrev];
				}
			}
			npos = npos + 1;
			next;
		}
		npos = npos + 1;
	}
	# ERROR
	idPrev = length(stackPos);
	if(idPrev > 0) {
		# process all non-closed tags;
		if(idPrev == 1) {
			nr = stackPos[idPrev];
			tk.df[nr, "Err"] = stackType[idPrev];
			tk.df[nr, "nE"]  = npos;
		} else {
			for(id in seq(idPrev, 1, by=-1)) {
				nr = stackPos[id];
				tk.df[nr, "Err"] = stackType[id];
				tk.df[nr, "nE"]  = npos;
			}
		}
	}
	if(warn && any(tk.df$Err != 0)) {
		warning("Mismatched parenthesis!");
	}
	class(tk.df) = c("code", class(tk.df));
	return(tk.df);
}

### Extract Content
# - between parenthesis;
extractParenth = function(x, pos = NULL, nested = FALSE, warn = FALSE,
		exclude.Parenth = TRUE) {
	if(is.null(pos)) {
		pos = parseParenth(x);
	}
	len = length(x);
	if(len == 0) return(character(0));
	if(warn) {
		r = lapply(seq(len), function(id) {
			extractParenth1(x[[id]], pos[[id]], nested=nested, warn = id,
				exclude.Parenth=exclude.Parenth);
		})
	} else {
		r = lapply(seq(len), function(id) {
			extractParenth1(x[[id]], pos[[id]], nested=nested, warn = NULL,
				exclude.Parenth=exclude.Parenth);
		})
	}
	return(r);
}
extractParenth1 = function(x, pos = NULL, nested=FALSE, warn=NULL,
		exclude.Parenth = TRUE) {
	if(is.null(pos)) {
		pos = parseParenth(x);
	}
	# NO parenthesis:
	if(nrow(pos) == 0) return(character(0));
	if( ! nested) pos = pos[ pos$Nested == FALSE, ];
	if( ! is.null(warn) && any(pos$Err != 0)) {
		warning("Mismatched parenthesis in: ", warn);
	}
	if(exclude.Parenth) {
		FUN = function(id) {
			nnS = pos$nS[id]; nnE = pos$nE[id];
			if(nnE > nnS + 1) { nnS = nnS + 1; nnE = nnE - 1; }
			else return("");
			substr(x, nnS, nnE);
		}
	} else {
		FUN = function(id) {
			substr(x, pos$nS[id], pos$nE[id]);
		}
	}
	sapply(seq(nrow(pos)), FUN);
}

isErrParenth = function(x) {
	sapply(seq(length(x)), function(id) {
		any(x[[id]]$Err != 0);
	})
}

countErrParenth = function(x, isErr=NULL) {
	if(is.null(isErr)) {
		isErr = isErrParenth(x);
	}
	sapply(x[isErr], function(x) {
		sum(x$Err != 0);
	})
}

hasParenth = function(x) {
	isP = sapply(x, function(x) {
		return(nrow(x) > 0);
	})
	return(isP);
}

stripParenth = function(x, pos = NULL, warn = FALSE) {
	if(is.null(pos)) pos = parseParenth(x, warn=warn);
	hasP = hasParenth(pos);
	sumP = sum(hasP);
	if(sumP == 0) return(x);
	#
	tmpStr = x[hasP]; tmpPos = pos[hasP];
	strip = function(id) {
		tmpStr = tmpStr[[id]];
		pos = tmpPos[[id]];
		pos = pos[pos$Nested == FALSE & pos$Err == 0, ];
		if(nrow(pos) == 0) return(tmpStr);
		# Note: assumes ordered positions;
		sz = nchar(tmpStr);
		nr = nrow(pos);
		nS = pos$nE + 1;
		# Head
		if(pos$nS[1] == 1) {
			if(nr >= 2) {
				nE = pos$nS[seq(2, nr)] - 1;
			} else nE = sz;
		} else {
			nS = c(1, nS);
			nE = c(pos$nS - 1, sz);
		}
		# Tail
		if(pos$nE[nr] + 1 >= sz) {
			len = length(nS);
			nS = nS[ - len];
			nE = nE[ - len];
		}
		len = length(nS);
		if(len == 0) return("");
		# does NOT function:
		# txt = substr(tmpStr, nS, nE);
		txt = sapply(seq(len), function(id) substr(tmpStr, nS[id], nE[id]));
		txt = paste0(txt, collapse="");
		return(txt);
	}
	x[hasP] = sapply(seq(sumP), strip);
	return(x);
}

summary.AbstractParenth = function(x, npos=NULL) {
	# npos = position of parenthesis;
	if(is.null(npos)) {
		npos = parseParenth(x$Abstract);
	}
	count.f = function(pred) length(unique(abstracts$PMID[pred]));
	# Abstracts with Parenthesis:
	hasP = hasParenth(npos);
	nP = count.f(hasP);
	# Abstracts with Errors:
	isErr = isErrParenth(npos);
	nErr = count.f(isErr);
	# Top-Level Parenthesis:
	nonNested = sapply(npos, function(x) {
		if(nrow(x) == 0) return(0);
		isTop = ( ! x$Nested) & (x$Err == 0);
		return(length(isTop));
	})
	nNN = sum(nonNested);
	#
	str = paste(
		c("Total Abstracts", "Non-Nested", "Errors"),
		c(nP, nNN, nErr), sep=": ");
	cat(str, sep="\n");
	invisible(str);
}

### Nested
# - extract content between parenthesis with nested parenthesis;
extractNested = function(x, pos, simplify=TRUE) {
	len = length(x);
	if(len > 1) {
		sContent = lapply(seq(len), function(id) {
			extractNested(x[[id]], pos=pos[[id]]);
		})
		if(simplify) sContent = unlist(sContent);
		return(sContent);
	}
	posNested = pos[pos$hasNested, ];
	if(nrow(posNested) == 0) return(character(0));
	sContent = extractParenth(x, pos=posNested);
	return(sContent);
}


### Wikipedia Corpus

# - file = filename;
# - n = Line buffer;
# - Records are separated by empty lines;
# - each record may contain 1 or more lines of text;
# - the file is processed sequentially,
#   thus avoiding to load huge files into memory;
read.txt.wiki = function(file, pattern, n=4000, verbose=TRUE) {
	con = file(file);
	open(con, "rt");
	extract = function(str, isWord) {
		len = length(str);
		nposR = which(nchar(str) == 0);
		last  = length(nposR);
		if(last == 0) {
			return(list(str));
		} else {
			if(nposR[[1]] > 1) nposR = c(0, nposR);
			if(nposR[[last]] < len) nposR = c(nposR, len + 1);
		}
		nposW = which(isWord);
		nposT = c(nposR, nposW);
		tkT = c(rep(FALSE, length(nposR)), rep(TRUE, length(nposW)));
		id  = order(nposT);
		nposT = nposT[id];
		tkT   = tkT[id];
		len = length(tkT);
		# lim = xor(tkT, c(tkT[-1], tkT[len]));
		isUnique = ! (tkT & c(tkT[-1], tkT[len]));
		nposT = nposT[isUnique];
		tkT   = tkT[isUnique];
		id  = which(tkT);
		idS = id - 1; idE = id + 1;
		nposS = nposT[idS];
		nposE = nposT[idE];
		len = length(nposS);
		lapply(seq(len), function(id) {
			str[seq(nposS[id] + 1, nposE[id] - 1)];
		})
	}
	#
	buffer = NULL;
	sRez   = list();
	repeat({
		input = readLines(con=con, n=n);
		# TODO: process last Record;
		if(length(input) == 0) break;
		buffer = c(buffer, input);
		# Cut before Last Record:
		# Note: Last Record may be incomplete;
		len = length(buffer);
		if(len == 0) next; # should never happen;
		for(npos in seq(len, 1, by=-1)) {
			if(nchar(buffer[npos]) == 0) break;
		}
		if(npos > 1) {
			str = buffer[1:npos];
			if(npos == len) { buffer = NULL; }
			else {
				buffer = buffer[seq(npos + 1, len)];
			}
			isWord = grepl(pattern, str);
			if(any(isWord)) {
				if(verbose) cat("Found!\n");
				sRez = c(sRez, extract(str, isWord));
			}
		}
		# Debug:
		# print(tail(buffer, 1));
	})
	close(con);
	# Last Record:
	if(length(buffer) > 0) {
		str = buffer;
		isWord = grepl(pattern, str);
		if(any(isWord)) {
			if(verbose) cat("\nFound!\n\n");
			sRez = c(sRez, extract(str, isWord));
		}
	}
	return(sRez);
}


######################

strsplitTokens = function(x, sep=";\\s*+", simplify=TRUE, perl=TRUE) {
	FUN = function(x) strsplit(x, sep, perl=perl);
	cleanFUN = function(x) {
		if(length(x) == 0) return(character(0));
		x[nzchar(x)];
	}
	#
	if(is.list(x)) {
		r = lapply(x, FUN);
		r = lapply(r, cleanFUN);
		if(simplify) {
			r = lapply(r, unlist);
		}
	} else {
		r = FUN(x);
		r = lapply(r, cleanFUN);
		if(simplify) {
			r = unlist(r);
		}
	}
	return(r);
}

# TODO:
# - split ", respectively", "ie, ";

######################

# "\d[-. %\d]*+"
is.string.numExt = function(x) {
	len = length(x);
	if(len == 1) return(is.string.numExt1(x));
	#
	rez = sapply(x, is.string.numExt1, USE.NAMES = FALSE);
	return(rez);
}
is.string.numExt1 = function(x) {
	if(is.na(x)) return(FALSE);
	if( ! nzchar(x)) return(FALSE);
	#
	x = utf8ToInt(x);
	len = length(x);
	isNum = FALSE;
	for(npos in seq(len)) {
		ch = x[npos];
		# 0 <= ch <= 9
		if(ch <= 57 && ch >= 48) {
			isNum = TRUE; next;
		}
		# " .-=<>%"
		if(ch %in% c(32, 46, 45, 61, 60, 62, 37, 160,
			183, 177, 8804, 8805)) next;
		return(FALSE);
	}
	return(isNum);
}

# "p = 0.n", "p <= 0.n"
is.string.pVal = function(x) {
	len = length(x);
	if(len == 1) return(is.string.pVal1(x));
	#
	rez = sapply(x, is.string.pVal1, USE.NAMES = FALSE);
	return(rez);
}
is.string.pVal1 = function(x) {
	if(is.na(x)) return(FALSE);
	if( ! nzchar(x)) return(FALSE);
	#
	x = utf8ToInt(x);
	len = length(x);
	npos = 1;
	while(npos < len) {
		if(x[npos] == 32 || x[npos] == 160) {
			npos = npos + 1;
		} else break;
	}
	# "p"
	if((x[npos] != 112) && (x[npos] != 80)) return(FALSE);
	npos = npos + 1;
	if(npos > len) return(FALSE);
	# " = "
	# Space:
	while(npos < len) {
		# "\u820n" = ?? used ?? as space;
		if(x[npos] == 32 || x[npos] == 160 ||
				x[npos] == 8201 || x[npos] == 8202 || x[npos] == 8203) {
			npos = npos + 1;
		} else break;
	}
	if(x[npos] != 61 && x[npos] != 60 && x[npos] != 8804 && x[npos] != 62) {
		# 62: "P > ";
		return(FALSE);
	}
	npos = npos + 1;
	if(npos > len) return(FALSE);
	if(x[npos] == 60) npos = npos + 1; # "p <= "
	if(npos > len) return(FALSE);
	# Space:
	while(npos < len) {
		# "\u820n" = ?? used ?? as space;
		if(x[npos] == 32 || x[npos] == 160 ||
				x[npos] == 8201 || x[npos] == 8202 || x[npos] == 8203) {
			npos = npos + 1;
		} else break;
	}
	# "0.n"
	if(x[npos] == 48) npos = npos + 1;
	if(x[npos] == 46 || x[npos] == 183) {
		npos = npos + 1;
	} else return(FALSE);
	while(npos <= len) {
		if(x[npos] >= 48 && x[npos] <= 57) {
			npos = npos + 1;
		} else return(FALSE); # TODO: trailing space;
	}
	return(TRUE);
}


# "N = n"
is.string.nVal = function(x) {
	len = length(x);
	if(len == 1) return(is.string.nVal1(x));
	#
	rez = sapply(x, is.string.nVal1, USE.NAMES = FALSE);
	return(rez);
}
is.string.nVal1 = function(x) {
	if(is.na(x)) return(FALSE);
	if( ! nzchar(x)) return(FALSE);
	#
	x = utf8ToInt(x);
	len = length(x);
	npos = 1;
	while(npos < len) {
		if(x[npos] == 32 || x[npos] == 160) {
			npos = npos + 1;
		} else break;
	}
	# "N"
	if((x[npos] != 110) && (x[npos] != 78)) return(FALSE);
	npos = npos + 1;
	if(npos > len) return(FALSE);
	# " = "
	# Space:
	while(npos < len) {
		# "\u820n" = ?? used ?? as space;
		if(x[npos] == 32 || x[npos] == 160 ||
				x[npos] == 8201 || x[npos] == 8202 || x[npos] == 8203 ||
				x[npos] == 8239) {
			npos = npos + 1;
		} else break;
	}
	if(x[npos] != 61) return(FALSE);
	npos = npos + 1;
	# Space:
	while(npos < len) {
		# "\u820n" = ?? used ?? as space;
		if(x[npos] == 32 || x[npos] == 160 ||
				x[npos] == 8201 || x[npos] == 8202 || x[npos] == 8203 ||
				x[npos] == 8239) {
			npos = npos + 1;
		} else break;
	}
	# "n,n"
	hasNum = FALSE;
	while(npos <= len) {
		if(x[npos] >= 48 && x[npos] <= 57) {
			npos = npos + 1;
			hasNum = TRUE;
		} else if(hasNum && x[npos] == 44) {
			npos = npos + 1; # TODO: check "(?=\d)";
		} else {
			return(FALSE); # TODO: trailing space;
		}
	}
	return(TRUE);
}


# TODO:
# - OR, RR, SMD/WMD (may have units);
# - also "95 % OR";
# - "pooled OR";
# - range, SD;

# 95% CI
is.string.CI95 = function(x) {
	regCI = paste0(
		"^95\\%[ ]*+CI", "[\\:, =\u2009\uA0]++",
		"\\d++(?:[.\uB7]\\d++)?+", "\\%?+", "[- ,;\u2009\uA0]++",
		"\\d++(?:[.\uB7]\\d++)?+", "\\%?+$")
	isCI = grepl(regCI, x, perl=TRUE);
	return(isCI);
}

# "n.n, 95% CI n-n"
is.string.ValCI95 = function(x) {
	regCI = paste0(
		"^\\d++(?:[.\uB7]\\d++)?+",
		"(?:;|,[ \uA0])[ \uA0]*+", "95\\%[ ]*+CI", "[\\:, =\u2009\uA0]++",
		"\\d++(?:[.\uB7]\\d++)?+", "\\%?+", "[- ,;\u2009\uA0]++",
		"\\d++(?:[.\uB7]\\d++)?+", "\\%?+$")
	isCI = grepl(regCI, x, perl=TRUE);
	return(isCI);
}

# "OR/RR/HR n.n, 95% CI n-n"
is.string.XR_CI95 = function(x) {
	regCI = paste0(
		"^[ORH]R", "[\\: =\u2009\uA0,]++",
		"\\d++(?:[.\uB7]\\d++)?+",
		"(?:;|,[ \uA0])[ \uA0]*+", "95\\%[ \u2009\uA0]*+CI",
		"[\\:, =\u2009\uA0]++",
		"\\d++(?:[.\uB7]\\d++)?+", "\\%?+", "[- ,;\u2009\uA0]++",
		"\\d++(?:[.\uB7]\\d++)?+", "\\%?+$")
	isCI = grepl(regCI, x, perl=TRUE);
	return(isCI);
}

# "OR/RR/HR n.n, 95% CI n-n, p = n"
is.string.XR_CI95_p = function(x) {
	regCI = paste0(
		"^(?:[ORH]R|pooled OR|pOR)", "[\\: =\u2009\uA0,]++",
		"\\d++(?:[.\uB7]\\d++)?+",
		"(?:;|,[ \uA0])[ \uA0]*+", "95\\%[ \u2009\uA0]*+CI",
		"[\\:, =\u2009\uA0]++",
		"\\d++(?:[.\uB7]\\d++)?+", "\\%?+", "[- ,;\u2009\uA0]++",
		"\\d++(?:[.\uB7]\\d++)?+", "\\%?+",
		"[\\:;, =\u2009\uA0]++", "[pP][ \uA0]*+",
		"[=<>\u2264]++[ \uA0]*+", "0?+[.\uB7]\\d++",
		"$")
	isCI = grepl(regCI, x, perl=TRUE);
	return(isCI);
}

# CRD Number
is.string.CRD = function(x) {
	regCRD = paste0(
		"^(?:[rR]egistration (?:[Nn]umber|[nN]o\\.|\\#)",
			"[ :\uA0]++)?+",
		"CRD", "[ \uA0]?+", "\\d++$")
	isCRD = grepl(regCRD, x, perl=TRUE);
	return(isCRD);
}


################

# Encode numbers => "n";
encode.num = function(x, ch="n") {
	# keep years;
	reg = paste0(
		"\\d++(?=\\.\\d)|",
		"(?<=\\.)\\d++|",
		"(?<!\\d)\\d{1,3}(?!\\d)");
	gsub(reg, ch, tmp[isNum], perl=TRUE);
}

######################

### Word Frequency
# - e.g. frequency table for keywords;
tableWords = function(x, sep=", *+", case.lower=TRUE, doSort=TRUE) {
	words = strsplit(x, split=sep, perl=TRUE);
	words = unlist(words);
	if(case.lower) {
		words = tolower(words);
	}
	wfreq = table(words);
	if(doSort) {
		id = order(wfreq, decreasing=TRUE);
		wfreq = wfreq[id];
	}
	return(wfreq);
}


