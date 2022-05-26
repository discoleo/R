########################
###
### Leonard Mada
### [the one and only]
###
### Pubmed
### Text Mining Tools
###
### draft v.0.1m


### Text Mining Tools
# - Identify position of parenthesis in Pubmed abstracts;
# - Extract content between parenthesis;



######################

### Basic Tools

extract.regex = function(x, pattern, gr=0, perl=TRUE, simplify=TRUE) {
	r = regexec(pattern, x, perl=perl);
	gr = gr + 1;
	s = lapply(seq(length(x)), function(id) {
		tmp = r[[id]];
		if(tmp[1] == -1) return("");
		substr(x[id], tmp[gr], attr(tmp, "match.length")[gr]);
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
extractParenthV = function(x, pos = NULL, nested = FALSE, warn = FALSE,
		exclude.Parenth = TRUE) {
	if(is.null(pos)) {
		pos = parseParenth(x);
	}
	len = length(x);
	if(len == 0) return(character(0));
	if(warn) {
		r = lapply(seq(len), function(id) {
			extractParenth(x[[id]], pos[[id]], nested=nested, warn = id,
				exclude.Parenth=exclude.Parenth);
		})
	} else {
		r = lapply(seq(len), function(id) {
			extractParenth(x[[id]], pos[[id]], nested=nested, warn = NULL,
				exclude.Parenth=exclude.Parenth);
		})
	}
	return(r);
}
extractParenth = function(x, pos = NULL, nested=FALSE, warn=NULL,
		exclude.Parenth = TRUE) {
	if(is.null(pos)) {
		pos = parseParenth(x);
	}
	if(nrow(pos) == 0) return(character(0));
	if( ! nested) pos = pos[ pos$Nested == FALSE, ];
	if( ! is.null(warn) && any(pos$Err != 0)) {
		warning("Mismatched parenthesis in: ", warn);
	}
	if(exclude.Parenth) {
		FUN = function(id) {
			nnS = pos$nS[id]; nnE = pos$nE[id];
			if(nnE > nnS + 1) { nnS = nnS + 1; nnE = nnE - 1; }
			else return(character(0));
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
	# " = "
	# Space:
	while(npos < len) {
		# "\u820n" = ?? used ?? as space;
		if(x[npos] == 32 || x[npos] == 160 ||
				x[npos] == 8201 || x[npos] == 8202 || x[npos] == 8203) {
			npos = npos + 1;
		} else break;
	}
	if(x[npos] != 61 && x[npos] != 60 && x[npos] != 8804) return(FALSE);
	npos = npos + 1;
	if(x[npos] == 60) npos = npos + 1; # "p <= "
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


