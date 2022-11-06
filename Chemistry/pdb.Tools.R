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
	LN = 4; # TODO
	n  = w %/% LN;
	nA = (length(id) %/% n);
	if(nA > 0) for(npos in seq(nA)) {
		nStart = (npos - 1) * n;
		idAA = id[seq(nStart, nStart + n)];
		idAA = format(as.character(idAA), justify = "centre", width=3);
		cat(idAA, sep=sep); cat("\n");
		cat(aa[seq(nStart, nStart + n)], sep=sep); cat("\n");
	}
	nStart = nA * n + 1;
	if(nStart < length(id)) {
		idAA = id[seq(nStart, length(id))];
		idAA = format(as.character(idAA), justify = "centre", width=3);
		cat(idAA, sep=sep); cat("\n");
		cat(aa[seq(nStart, length(id))], sep=sep); cat("\n");
	}
	invisible();
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
read.meta.pdb = function(path = ".", pattern="\\.ent\\.gz$", rm.res.string=TRUE) {
	files = list.files(path, pattern=pattern);
	res = lapply(files, function(name) {
		x = read.pdb(name);
		tt  = title.pdb(x);
		res = resolution(x);
		isRes = grepl("^(?i)Resolution", res);
		res = res[isRes][[1]];
		data.frame(title=tt, resolution=res);
	})
	res = do.call(rbind, res);
	pdb = extract.regex(files, "(?i)^pdb([^.]++)", gr=1)
	res = cbind(pdb, res);
	if(rm.res.string) {
		res$resolution = sub("^(?i)Resolution[. ]++", "", res$resolution, perl=TRUE);
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

#################
#################

