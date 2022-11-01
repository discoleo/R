
### Peptide Tools

# Peptide => Smile
# Smile => Peptide


as.smile = function(x, n=1) {
	s = strsplit(x, split=c());
	aaCodes = getAA(n = n);
	x.smi = lapply(s, function(s) paste0(c(aaCodes[s], "O"), collapse=""));
	if(length(x.smi) == 1) x.smi = x.smi[[1]];
	return(x.smi);
}

match.smile = function(x, n=1, collapse="") {
	smi = split.pp.smile(x);
	smiCodes = unlist(getAA(n=n));
	aa  = lapply(smi, function(smi) {
		last = length(smi);
		# TODO: check with carboxy-amides;
		aaLast = sub("O$", "", smi[last]);
		smi[last] = aaLast;
		ids = match(smi, smiCodes);
		aa  = names(smiCodes)[ids];
		if( ! is.null(collapse)) aa = paste0(aa, collapse=collapse);
		return(aa);
	})
	if(length(aa) == 1) aa = aa[[1]];
	return(aa);
}
split.pp.smile = function(x) {
	strsplit(x, "(?<=C\\(=O\\))(?=N(?!\\)))", perl=TRUE);
}

aaCodes = function(n=1) {
	if( ! any(n == c(0,1,3))) stop("Unsupported codes!");
	a3 = c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met",
		"Phe", "Pro", "Pyl", "Ser", "Sec", "Thr", "Trp", "Tyr", "Val", "Asx", "Glx", "Xaa", "Xle");
	a1 = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M",
		"F", "P", "O", "S", "U", "T", "W", "Y", "V", "B", "Z", "X", "J");
	name = c("Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine", "Glutamine",
		"Glutamic acid", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine",
		"Phenylalanine", "Proline", "Pyrrolysine", "Serine", "Selenocysteine", "Threonine", "Tryptophan",
		"Tyrosine", "Valine", "Aspartic acid or Asparagine", "Glutamic acid or Glutamine",
		"Any amino acid", "Leucine or Isoleucine");
	if(n == 1) {
		aa = as.list(a1);
		names(aa) = a3;
	} else if(n == 3) {
		aa = as.list(a3);
		names(aa) = a1;
	} else {
		aa = data.frame(A3=a3, A1=a1, Name=name);
	}
	return(aa);
}

getAA = function(n=1) {
	if( ! any(n == c(1,3))) stop("Unsupported codes!");
	AA.smi = list(
		Gly = "NCC(=O)",
		Ala = "N[C@@]([H])(C)C(=O)",
		# - OH, - SH
		Ser = "N[C@@]([H])(CO)C(=O)",
		Thr = "N[C@@]([H])([C@]([H])(O)C)C(=O)",
		Cys = "N[C@@]([H])(CS)C(=O)",
		Met = "N[C@@]([H])(CCSC)C(=O)",
		# Branched
		Val = "N[C@@]([H])(C(C)C)C(=O)",
		Leu = "N[C@@]([H])(CC(C)C)C(=O)",
		Ile = "N[C@@]([H])([C@]([H])(CC)C)C(=O)",
		# Aromatic
		Phe = "N[C@@]([H])(Cc1ccccc1)C(=O)",
		Tyr = "N[C@@]([H])(Cc1ccc(O)cc1)C(=O)",
		Trp = "N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)",
		His = "N[C@@]([H])(CC1=CN=C-N1)C(=O)",
		# Acidic
		Asp = "N[C@@]([H])(CC(=O)O)C(=O)",
		Asn = "N[C@@]([H])(CC(=O)N)C(=O)",
		Glu = "N[C@@]([H])(CCC(=O)O)C(=O)",
		Gln = "N[C@@]([H])(CCC(=O)N)C(=O)",
		#
		Pro = "N1[C@@]([H])(CCC1)C(=O)",
		# Basic
		Arg = "N[C@@]([H])(CCCNC(=N)N)C(=O)",
		Lys = "N[C@@]([H])(CCCCN)C(=O)"
	);
	if(n == 1) {
		stdAA = aaCodes(n=1);
		codeAA  = stdAA[names(AA.smi)];
		names(AA.smi) = codeAA;
	}
	return(AA.smi);
}


