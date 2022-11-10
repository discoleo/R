####################
###
### Chemistry Tools
### Amino Acids
###
### Leonard Mada
###
### v.0.1a



### Amino-Acids

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

groups.aa = function(n=3) {
	if(n == 3) {
		aa = c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met",
			"Phe", "Pro", "Pyl", "Ser", "Sec", "Thr", "Trp", "Tyr", "Val", "Asx", "Glx", "Xaa", "Xle");
	} else if(n == 1) {
		aa = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M",
			"F", "P", "O", "S", "U", "T", "W", "Y", "V", "B", "Z", "X", "J");
	} else stop("Unsupported Code!");
	# Aliphatic = 1; Alcohol, Thiol = 2; Branched = 3;
	# Aromatic (includes His & Trp) = 4; Pro = 5; Acidic = 6; Basic = 7;
	type = c(1, 7, 6,6, 2, 6,6, 1, 4, 3,3, 7, 2, 4, 5, 7, 2, 2, 2, 4, 4, 3, 6, 6, NA, 3);
	type = factor(type, levels=1:7, labels=c("Aliphatic", "Alcohol, Thiol", "Branched",
		"Aromatic", "Pro", "Acidic", "Basic"));
	names(type) = aa;
	return(type);
}

as.aaType = function(x) {
	n = nchar(x[[1]]);
	type = groups.aa(n=n);
	aa   = names(type);
	isUpper = function(x) TRUE; # TODO
	if(n == 3 && isUpper(substr(x[[1]], 2, 2))) aa = toupper(aa);
	id = match(x, aa);
	return(type[id]);
}
