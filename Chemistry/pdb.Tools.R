#################
###
### PDB Tools
###
### Leonard Mada
###
###


####################

### Helper Functions

library(Rpdb)


####################

print.str = function(x, w=80) {
	cat(strwrap(x, w), sep="\n");
}

cat.aa = function(x, w=90, sep=" ") {
	id = x$idRes;
	aa = x$aa;
	# Case: id > 3 digits;
	# LN = Length of N;
	LN = max(3, length.digits(length(id)));
	n  = w %/% (LN + nchar(sep));
	nA = (length(id) %/% n);
	if(nA > 0) for(npos in seq(nA)) {
		nStart = (npos - 1) * n;
		idAA = id[seq(nStart, nStart + n)];
		idAA = format(as.character(idAA), justify = "centre", width=LN);
		cat(idAA, sep=sep); cat("\n");
		# AA:
		strAA = aa[seq(nStart, nStart + n)];
		if(LN > 3) strAA = format(strAA, justify = "centre", width=LN);
		cat(strAA, sep=sep); cat("\n");
	}
	nStart = nA * n + 1;
	if(nStart < length(id)) {
		idAA = id[seq(nStart, length(id))];
		idAA = format(as.character(idAA), justify = "centre", width=LN);
		cat(idAA, sep=sep); cat("\n");
		# AA:
		strAA = aa[seq(nStart, length(id))];
		if(LN > 3) strAA = format(strAA, justify = "centre", width=LN);
		cat(strAA, sep=sep); cat("\n");
	}
	invisible();
}
length.digits = function(x) {
	# TODO:
	return(3);
}

extract.regex = function(x, pattern, gr=0, perl=TRUE, simplify=TRUE) {
	r = regexec(pattern, x, perl=perl);
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

title.pdb = function(x) {
	t1 = x$title[1];
	t1 = sub("TITLE +", "", t1);
	t1 = sub(" *+$", "", t1, perl=TRUE);
	t2 = x$title[-1];
	if(length(t2) > 0) {
		t2 = sub("TITLE +[0-9]+ +", "", t2);
		t2 = sub(" *+$", "", t2, perl=TRUE);
		t2 = paste0(t2, collapse=" ");
		t1 = paste0(t1, " ", t2);
	}
	class(t1) = c("str", class(t1));
	return(t1);
}

resolution = function(x, sep=" ", rm.related=TRUE) {
	res0 = x$remark;
	isRes = grepl("(?i)Resolution", res0);
	isRes = which(isRes);
	id = sub("^REMARK +", "", res0[isRes]);
	id = extract.regex(id, "^[0-9]++", perl=TRUE);
	id = unique(as.numeric(id));
	res = sapply(id, function(id) {
		patt  = paste0("^REMARK +", id, " +");
		isRes = grepl(patt, res0);
		res = res0[isRes];
		res = sub(patt, "", res);
		res = sub(" ++$", "", res, perl=TRUE);
		res = res[nchar(res) > 0];
		res = paste0(res, collapse=sep);
		return(res);
	})
	if(rm.related) {
		isRelated = grepl("^(?i)Related ", res);
		isRelated = which(isRelated);
		if(length(isRelated) > 0) {
			res = res[ - isRelated];
			id  = id[ - isRelated];
		}
	}
	attr(res, "id") = id;
	return(res)
}
# Basic Summary of Collection of PDBs
read.meta.pdb = function(path = ".", pattern="\\.ent\\.gz$", FUN=NULL, rm.res.string=TRUE) {
	files = list.files(path, pattern=pattern);
	res = lapply(files, function(name) {
		x = read.pdb(name);
		tt  = title.pdb(x);
		res = resolution(x);
		isRes = grepl("^(?i)Resolution", res);
		res = res[isRes][[1]];
		nChains = length(unique(x$atoms$chainid));
		nAA = length.chains.aa(x);
		rez = data.frame(Title=tt, Resolution=res, Chains=nChains,
			maxAA = max(nAA), nAA = paste0(nAA, collapse=", "));
		if( ! is.null(FUN)) {
			rez$FUN = FUN(x);
		}
		return(rez);
	})
	res = do.call(rbind, res);
	pdb = extract.regex(files, "(?i)^pdb([^.]++)", gr=1)
	res = cbind(pdb, res); names(res)[1] = "PDB";
	if(rm.res.string) {
		res$Resolution = sub("^(?i)Resolution[. ]++", "", res$Resolution, perl=TRUE);
	}
	return(res);
}

### Chains

chains = function(x) {
	unique(x$atoms$chainid);
}

length.chains = function(x) {
	table(x$atoms$chainid);
}
length.chains.aa = function(x, drop="HETATM") {
	chs = chains(x);
	doDrop = FALSE;
	if( ! is.null(drop) && ! is.na(drop) ) {
		if(is.logical(drop)) {
			if(drop) { drop = "HETATM"; doDrop = TRUE; }
			else drop = character(0); # BREAK!
		} else doDrop = TRUE;
	}
	nAA = sapply(chs, function(ch) {
		idCh = x$atoms$chainid == ch;
		idCh = which(idCh);
		if(doDrop) {
			tmp  = x$atoms$recname[idCh];
			idCh = idCh[ ! tmp %in% drop];
		}
		aa = unique(x$atoms[idCh, c("resid")]);
		return(length(aa));
	});
	return(nAA);
}

chain = function(ch, x, drop="HETATM") {
	idCh = x$atoms$chainid == ch;
	if( ! any(idCh)) return(data.frame(atom=character(0), id=numeric(0), idRes=numeric(0)));
	#
	atoms  = x$atoms$elename[idCh];
	idAtom = x$atoms$resname[idCh];
	idRes  = x$atoms$resid[idCh];
	rez = data.frame(atom=atoms, id=idAtom, idRes=idRes);
	if( ! is.null(drop) && ! is.na(drop) ) {
		tmp = x$atoms$recname[idCh];
		rez = rez[ ! tmp %in% drop, ];
	}
	return(rez);
}

atoms.chain = function(ch, x) {
	idCh = x$atoms$chainid == ch;
	if( ! any(idCh)) return(character(0));
	atoms = x$atoms$elename[idCh];
	names(atoms) = x$atoms$resname[idCh];
	return(atoms);
}

as.aa = function(x, ch, drop="HETATM") {
	return(unique.aa(x, ch=ch, drop=drop));
}
unique.aa = function(x, ch, drop="HETATM") {
	idCh = x$atoms$chainid == ch;
	idCh = which(idCh);
	if( ! is.null(drop) && ! is.na(drop) ) {
		if(is.logical(drop)) {
			if(drop) { drop = "HETATM"; }
			else drop = character(0); # BREAK!
		}
		tmp  = x$atoms$recname[idCh];
		idCh = idCh[ ! tmp %in% drop];
	}
	aa = unique(x$atoms[idCh, c("resid", "resname")]);
	names(aa) = c("idRes", "aa");
	return(aa);
}

### Filter Chains:
which.oligo = function(x, len=30, drop="HETATM") {
	n  = length.chains.aa(x, drop=drop);
	id = which(n <= len);
	if(length(id) == 0) return(character(0));
	ch = chains(x);
	return(ch[id]);
}
filter.oligo = function(x, len=30, drop="HETATM", maxChains=0) {
	idCh = which.oligo(x, len=len, drop=drop);
	if(length(idCh) == 0) return(data.frame(idRes = numeric(0), aa = character(0)));
	if(maxChains > 0) {
		maxChains = min(maxChains, length(idCh));
		idCh = idCh[seq(maxChains)];
	}
	aa = lapply(idCh, function(ch) as.aa(x, ch=ch, drop=drop));
	if(length(aa) == 1) aa = aa[[1]];
	return(aa);
}

#################
#################

### Simple Analysis

# Frequency of AA at every position:
# - crude summary;
summary.set.pp = function(x, sep=" ", filter=6) {
	if(is.character(sep)) {
		pp = strsplit(x, sep);
	} else {
		pp = pp = strsplit(x, sep, fixed=TRUE);
	}
	if(filter > 0) {
		id = sapply(pp, length);
		pp = pp[id >= filter];
	}
	pp = lapply(pp, function(p) {
		n = length(p);
		return(data.frame(Pos = seq(n), AA=p));
	});
	pp = do.call(rbind, pp);
	pp$Count = 1;
	pp = aggregate(Count ~ ., pp, length);
	id = order(pp$Pos, pp$AA);
	pp = pp[id,];
	return(pp);
}

library(ggplot2)

plot.aa.freq = function(x) {
	ggplot(x, mapping=aes(x=AA, y=Pos, size=Count)) +
		geom_point();
}

