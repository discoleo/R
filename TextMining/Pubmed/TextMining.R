########################
###
### Leonard Mada
### [the one and only]
###
### Pubmed
### Text Mining Tools
###
### draft v.0.1f


### Text Mining Tools
# - Identify position of parenthesis in Pubmed abstracts;
# - Extract content between parenthesis;



######################

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
		# TODO: all non-closed tags;
		nr = stackPos[idPrev];
		tk.df[nr, "Err"] = stackType[idPrev];
		tk.df[nr, "nE"]  = npos;
	}
	if(warn && any(tk.df$Err != 0)) {
		warning("Mismatched parenthesis!");
	}
	class(tk.df) = c("code", class(tk.df));
	return(tk.df);
}

### Extract Content
# - between parenthesis;
extractParenth = function(x, pos = NULL, nested=FALSE, warn=NULL) {
	if(is.null(pos)) {
		pos = parseParenth(x);
	}
	if(nrow(pos) == 0) return(character(0));
	if( ! nested) pos = pos[ pos$Nested == FALSE, ];
	if( ! is.null(warn) && any(pos$Err != 0)) {
		warning("Mismatched parenthesis in: ", warn);
	}
	FUN = function(id) {
		substr(x, pos$nS[id], pos$nE[id]);
	}
	sapply(seq(nrow(pos)), FUN);
}

isErrParenth = function(x) {
	sapply(seq(length(x)), function(id) {
		any(x[[id]]$Err != 0);
	})
}

hasParenth = function(x) {
	isP = sapply(x, function(x) {
		return(nrow(x) > 0);
	})
	return(isP);
}

summary.AbstractParenth = function(x, npos=NULL) {
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


