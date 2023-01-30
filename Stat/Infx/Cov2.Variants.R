

### Helper Tools

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

##############

### Mutations:
mutations = function(x, verbose=TRUE, debug=FALSE) {
	z = lapply(seq(nrow(x)), function(id) {
		mut = unlist(strsplit(x$aa_definition[id], ","));
		data.frame(V = x$lineage[id], Mutation = mut);
	})
	z = do.call(rbind, z)
	
	### Protein:
	z$P = extract.regex(z$Mutation, "^[^:]++", perl=TRUE, verbose=debug);
	
	### AA:
	z$AA = extract.regex(z$Mutation, "\\:(.++)$", perl=TRUE, gr=1, verbose=debug);
	# Test:
	if(verbose) {
		if(any(is.na(z$AA) | (nchar(z$AA) == 0)))
			warning("Extraction went wrong!");
	}
	
	z$Pos = as.numeric(extract.regex(z$AA, "[0-9]++", perl=TRUE, verbose=debug));
	z$AAi = extract.regex(z$AA, "^[A-Za-z]++", perl=TRUE, verbose=debug);
	z$AAm = extract.regex(z$AA, "[A-Za-z]++$", perl=TRUE, verbose=debug);
	
	# Polymorphisms:
	z$Polymorphism = (z$AAi == z$AAm)
	if(verbose) {
		cat("Gene Polymorphisms:\n")
		print(table(z$Polymorphism));
	}
	
	return(z);
}

###############

### Cov-2 Lineages
# download from: NCBI website
# https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/activ/?report=download_lineages

x = read.csv("Lineages.2023-01-25.csv")

# Note:
# V = only shorthand code;
z = mutations(x);

head(z)

