########################
###
### Leonard Mada
### [the one and only]
###
### Code Tools
###
### draft v.0.2b


### Tools to Process Formulas & Expressions


###############
### History ###
###############

### draft v.0.2a - v.0.2b:
# - moved Code tools to this file: Tools.Code.R;
### draft v.0.1g - v.0.1g-improve:
# - cut.code() into code blocks;
# - small improvements; [v.0.1g-improve]
### draft v.0.1f - v.0.1f-refactor:
# - extract code tokens from R code;
# - [refactored] uniform result;


### TODO
# - Explore: remindR;


########################

### Parser

# minimalistic parser:
parse.simple = function(x, eol="\n", all.tokens=FALSE) {
	len = nchar(x);
	# Type:
	#  1 = "string", 2 = r"raw string", 8/9 = unmatched string (raw string),
	#  11 = "#", 12 = "#" with missing EOL,
	#  23 = "{", 24 = "(", 25 = "[";
	tk.df = data.frame(nS=integer(0), nE=integer(0), id=integer(0),
			Type=integer(0), Nested=logical(0));
	#
	is.hex = function(ch) {
		# Note: only for 1 character!
		return((ch >= "0" && ch <= "9") ||
			(ch >= "A" && ch <= "F") ||
			(ch >= "a" && ch <= "f"));
	}
	npos = 1;
	while(npos <= len) {
		s = substr(x, npos, npos);
		# State: COMMENT
		if(s == "#") {
			nS = npos;
			isEOL = FALSE;
			while(npos < len) {
				npos = npos + 1;
				if(substr(x, npos, npos) == eol) { isEOL = TRUE; break; }
			}
			nr = nrow(tk.df) + 1;
			nE = if(isEOL) npos else npos + 1;
			TYPE = if(isEOL) 11 else 12;
			tk.df[nr,] = data.frame(nS, nE, nr, TYPE, FALSE);
			npos = npos + 1; next;
		}
		# State: STRING
		if(s == "\"" || s == "'") {
			nS = npos;
			if(npos > 1 && substr(x, npos-1, npos-1) == "r") {
				# Raw String: r"(...)";
				nS = nS - 1;
				nr = nrow(tk.df) + 1;
				if(npos + 3 > len) {
					nE = len + 1;
					TYPE = 9; # type = 8 vs 9!
					tk.df[nr,] = data.frame(nS, nE, nr, TYPE, FALSE);
					warning("Unmatched String!");
					break;
				}
				npos = npos + 1;
				ch1 = substr(x, npos, npos);
				isMatched = FALSE;
				while(npos < len - 2) {
					npos = npos + 1;
					if(substr(x, npos, npos) == ch1 && substr(x, npos+1, npos+1) == s) {
						isMatched = TRUE; break;
					}
				}
				if( ! isMatched) warning("Unmatched String!");
				TYPE = if(isMatched) 2 else 9;
				nE = if(isMatched) npos else npos + 1;
				tk.df[nr,] = data.frame(nS, nE, nr, TYPE, FALSE);
				npos = npos + 1; next;
			}
			# ELSE:
			while(npos < len) {
				# Standard String
				npos = npos + 1;
				se = substr(x, npos, npos);
				# Escape:
				if(se == "\\") {
					npos = npos + 1;
					# simple escape vs Unicode:
					if(substr(x, npos, npos) != "u") next;
					len.end = min(len, npos + 4); # max length of HEX = 4;
					npos = npos + 1;
					isAllHex = TRUE;
					while(npos <= len.end) {
						se = substr(x, npos, npos);
						if( ! is.hex(se)) { isAllHex = FALSE; break; }
						npos = npos + 1;
					}
					if(isAllHex) next;
				}
				if(se == s) break;
			}
			# n.str[[2]] = c(n.str[[2]], npos);
			nr = nrow(tk.df) + 1; nE = npos;
			TYPE = 1;
			tk.df[nr,] = data.frame(nS, nE, nr, TYPE, FALSE);
			npos = npos + 1; next;
		}
		# Brackets
		if(all.tokens) {
			type = 0;
			if(s == "{") {
				type = 23;
			} else if(s == "(") {
				type = 24;
			} else if(s == "[") {
				type = 25;
				# Note: not yet "[[";
			}
			if(type > 0) {
				nr = nrow(tk.df);
				tk.df[nr + 1, ] = data.frame(npos, NA, nr + 1, type, FALSE);
				if(nr > 1) {
					for(nr.prev in seq(nr, 1, by=-1)) {
						if(tk.df$Type[nr.prev] < 20) next; # TYPE < 20
						if(is.na(tk.df$nE[nr.prev])) tk.df$Nested[nr.prev] = TRUE;
						break;
					}
				}
				npos = npos + 1; next;
			}
			#
			type = 0;
			if(s == "}") {
				type = 23;
			} else if(s == ")") {
				type = 24;
			} else if(s == "]") {
				type = 25;
				# Note: not yet "]]";
			}
			if(type > 0) {
				# TODO: check for un-matched token!
				id = tail(tk.df$id[is.na(tk.df$nE) & tk.df$Type == type], 1);
				tk.df$nE[id] = npos;
				npos = npos + 1; next;
			}
		}
		npos = npos + 1;
	}
	class(tk.df) = c("code", class(tk.df));
	return(tk.df);
}
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
parse.pC = function(x, eol="\n") {
	len = nchar(x);
	# Type:
	#  1 = "string", 2 = r"raw string", 8/9 = unmatched string (raw string),
	#  11 = "#", 12 = "#" with missing EOL,
	#  23 = "{", 24 = "(", 25 = "[";
	tk.df = data.frame(nS=integer(0), nE=integer(0), id=integer(0),
			Type=integer(0), Nested=logical(0));
	#
	is.hex = function(ch) {
		# Note: only for 1 character!
		return((ch >= "0" && ch <= "9") ||
			(ch >= "A" && ch <= "F") ||
			(ch >= "a" && ch <= "f"));
	}
	npos = 1;
	while(npos <= len) {
		s = substr(x, npos, npos);
		# State: COMMENT
		if(s == "#") {
			nS = npos;
			isEOL = FALSE;
			while(npos < len) {
				npos = npos + 1;
				if(substr(x, npos, npos) == eol) { isEOL = TRUE; break; }
			}
			nr = nrow(tk.df) + 1;
			nE = if(isEOL) npos else npos + 1;
			TYPE = if(isEOL) 11 else 12;
			tk.df[nr,] = data.frame(nS, nE, nr, TYPE, FALSE);
			npos = npos + 1; next;
		}
		# State: STRING
		if(s == "\"" || s == "'") {
			nS = npos;
			if(npos > 1 && substr(x, npos-1, npos-1) == "r") {
				# Raw String: r"(...)";
				nS = nS - 1;
				nr = nrow(tk.df) + 1;
				if(npos + 3 > len) {
					nE = len + 1;
					TYPE = 9; # type = 8 vs 9!
					tk.df[nr,] = data.frame(nS, nE, nr, TYPE, FALSE);
					warning("Unmatched String!");
					break;
				}
				npos = npos + 1;
				ch1 = substr(x, npos, npos);
				isMatched = FALSE;
				while(npos < len - 2) {
					npos = npos + 1;
					if(substr(x, npos, npos) == ch1 && substr(x, npos+1, npos+1) == s) {
						isMatched = TRUE; break;
					}
				}
				if( ! isMatched) warning("Unmatched String!");
				TYPE = if(isMatched) 2 else 9;
				nE = if(isMatched) npos else npos + 1;
				tk.df[nr,] = data.frame(nS, nE, nr, TYPE, FALSE);
				npos = npos + 1; next;
			}
			# ELSE:
			while(npos < len) {
				# Standard String
				npos = npos + 1;
				se = substr(x, npos, npos);
				# Escape:
				if(se == "\\") {
					npos = npos + 1;
					# simple escape vs Unicode:
					if(substr(x, npos, npos) != "u") next;
					len.end = min(len, npos + 4); # max length of HEX = 4;
					npos = npos + 1;
					isAllHex = TRUE;
					while(npos <= len.end) {
						se = substr(x, npos, npos);
						if( ! is.hex(se)) { isAllHex = FALSE; break; }
						npos = npos + 1;
					}
					if(isAllHex) next;
				}
				if(se == s) break;
			}
			# n.str[[2]] = c(n.str[[2]], npos);
			nr = nrow(tk.df) + 1; nE = npos;
			TYPE = 1;
			tk.df[nr,] = data.frame(nS, nE, nr, TYPE, FALSE);
			npos = npos + 1; next;
		}
		# State: Pointer
		if(s == "<") {
			sP = c("p","o","i","n","t","e","r",":"," ","0","x");
			idP = 1; nposP = npos;
			while(nposP < len && idP <= length(sP)) {
				nposP = nposP + 1;
				se = substr(x, nposP, nposP);
				if(se == sP[idP]) { idP = idP + 1; next; }
				break;
			}
			if(idP == length(sP) + 1) {
				while(nposP < len) {
					nposP = nposP + 1;
					se = substr(x, nposP, nposP);
					if(is.hex(se)) { nposP = nposP + 1; next; }
					break;
				}
				if(nposP <= len) {
					se = substr(x, nposP, nposP);
					if(se == ">") {
						nr = nrow(tk.df) + 1; nE = nposP;
						TYPE = 50;
						tk.df[nr,] = data.frame(npos, nE, nr, TYPE, FALSE);
						npos = nposP + 1; next;
					}
				}
			}
		}
		npos = npos + 1;
	}
	class(tk.df) = c("code", class(tk.df));
	return(tk.df);
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
