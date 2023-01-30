########################
###
### Leonard Mada
### [the one and only]
###
### SARS Cov-2
### Lineages: Mutations
###
### draft v.0.1b


################

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

diff.lineage = function(v1, v2, data) {
	x1 = data[data$V == v1, ];
	x2 = data[data$V == v2, ];
	# Both:
	id = which(x1$Mutation %in% x2$Mutation);
	x1 = cbind(x1, id = "1");
	x1$id[id] = "B";
	# Only 2:
	Vb = x1$Mutation[id]; # Both
	id = which(x2$Mutation %in% Vb);
	x2 = cbind(x2, id = "2");
	x2$id[id] = "B";
	x = rbind(x1, x2);
	return(x);
}

### Summaries

# TODO:
# - incorrect, as list of mutations is NOT cumulative!
countMutations = function(x) {
	x = x[x$Polymorphism == FALSE, c("V", "Mutation")];
	count = tapply(x$Mutation, x$V, function(x) {
		length(x);
	})
	nms = names(count);
	count = data.frame(V = nms, count = as.numeric(count));
	return(count);
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


### Mutations: Spike protein
# Note:
# - does NOT reflect the relative abundance of the lineages;
# - the list of mutations seems NOT cumulative;
table(z$Pos[z$P == "S" & z$Polymorphism == F])

count = countMutations(z)

head(count)
table(count$count)


### New Mutations

# TODO:
# - seems that the list of mutations is NOT cumulative!
# B.1 => B.1.1 seems cumulative
diff.lineage("B.1.1", "B.1", data=z)
# but B.1.1 => B.1.1.529 is NOT cumulative anymore;
diff.lineage("B.1.1.529", "B.1.1", data=z)
diff.lineage("B.1.1.529", "BA.2", data=z)
diff.lineage("B.1.1.529", "BA.5", data=z)

