

##############
### Parser ###

### Minimalistic Parser:
# - Proof-of-concept parser: written entirely in R;
# - Native R-parser (base::parse, utils::getParseData) may be better;
# Possible Use-Cases:
# - Can parse code that is not fully R-compliant;
parse.simple = function(x, eol="\n", all.tokens=FALSE) {
	len = nchar(x);
	if(is.character(all.tokens)) {
		tokens = pmatch(all.tokens, c("all", "braces", "c"));
		if(is.na(tokens)) {
			stop("Wrong type of tokens!");
		} else if(tokens == 2) {
			tokens.braces = TRUE; tokens.c = FALSE;
		} else if(tokens == 3) {
			tokens.braces = FALSE; tokens.c = TRUE;
		} else {
			tokens.braces = tokens.c = all.tokens;
		}
	} else {
		tokens.braces = tokens.c = all.tokens;
	}
	# Type:
	#  1 = "string", 2 = r"raw string", 8/9 = unmatched string (raw string),
	#  11 = "#", 12 = "#" with missing EOL,
	#  23 = "{", 24 = "(", 25 = "[",
	#  50 = "<pointer: 0x...>";
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
		if(tokens.braces) {
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
		# State: Pointer
		if(tokens.c) {
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
		}
		npos = npos + 1;
	}
	class(tk.df) = c("code", class(tk.df));
	return(tk.df);
}

