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

#################
#################

