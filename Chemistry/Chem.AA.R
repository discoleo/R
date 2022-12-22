####################
###
### Chemistry Tools
### Amino Acids
###
### Leonard Mada
###
### v.0.1c



# TODO:
# - propensity for alpha-helix & beta-sheet;
# 1. Roseman, M. A. 1988. Hydrophilicity of polar amino acid side-chains is markedly reduced
#    by flanking peptide bonds. J. Mol. Biol. 200:513-522.
# 2. Koehl, P., and M. Levitt. 1999. Structure-based conformational preferences of amino acids.
#    Proc. Natl. Acad. Sci. U. S. A. 96:12524-12529.
# 3. Street, A. G., and S. L. Mayo. 1999. Intrinsic beta-sheet propensities result from
#    van der Waals interactions between side chains and the local backbone.
#    Proc. Natl. Acad. Sci. U. S. A. 96:9074-9076.
# 4. Pawar, A. P., K. F. Dubay, J. Zurdo, F. Chiti, M. Vendruscolo, and C. M. Dobson. 2005.
#    Prediction of “aggregation-prone” and “aggregation-susceptible” regions in proteins
#    associated with neurodegenerative diseases. J. Mol. Biol. 350:379-392.

# this file:
# source("Chem.AA.R");


### Amino-Acids

aaCodes = function(n=1) {
	if( ! any(n == c(0,1,3))) stop("Unsupported codes!");
	a3 = c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met",
		"Phe", "Pro", "Pyl", "Ser", "Sec", "Thr", "Trp", "Tyr", "Val", "Asx", "Glx", "Xaa", "Xle", "Css");
	a1 = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M",
		"F", "P", "O", "S", "U", "T", "W", "Y", "V", "B", "Z", "X", "J", "~");
	name = c("Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine", "Glutamine",
		"Glutamic acid", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine",
		"Phenylalanine", "Proline", "Pyrrolysine", "Serine", "Selenocysteine", "Threonine", "Tryptophan",
		"Tyrosine", "Valine", "Aspartic acid or Asparagine", "Glutamic acid or Glutamine",
		"Any amino acid", "Leucine or Isoleucine", "Cys-S-S-Cys");
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
			"Phe", "Pro", "Pyl", "Ser", "Sec", "Thr", "Trp", "Tyr", "Val", "Asx", "Glx",
			"Xaa", "Xle", "Css");
	} else if(n == 1) {
		aa = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M",
			"F", "P", "O", "S", "U", "T", "W", "Y", "V", "B", "Z", "X", "J", "~");
	} else stop("Unsupported Code!");
	# Aliphatic = 1; Alcohol, Thiol = 2; Branched = 3;
	# Aromatic (includes His & Trp) = 4; Pro = 5; Acidic = 6; Basic = 7;
	# Cys-S-S-Cys: behaves more like a branched/non-polar AA;
	# Sec = Se-Cysteine;
	type = c(1, 7, 6,6, 2, 6,6, 1, 4, 3,3, 7, 2, 4, 5, 7, 2, 2, 2, 4, 4, 3, 6, 6, NA, 3, 3);
	type = factor(type, levels=1:7, labels=c("Aliphatic", "Alcohol, Thiol", "Branched",
		"Aromatic", "Pro", "Acidic", "Basic"));
	names(type) = aa;
	return(type);
}

# Kyte-Doolittle Hydropathicity:
# J Mol Biol 157:105-132 1982
aaHydrophobicity = function() {
	coeff = c(0.0, # 'X'
	1.8, # 'A'
	2.5, # 'C'
	-3.5, # 'D'
	-3.5, # 'E'
	2.8, # 'F'
	-0.4, # 'G'
	-3.2, # 'H'
	4.5, # 'I'
	-3.9, # 'K'
	3.8, # 'L'
	1.9, # 'M'
	-3.5, # 'N'
	-1.6, # 'P'
	-3.5, # 'Q'
	-4.5, # 'R'
	-0.8, # 'S'
	-0.7, # 'T'
	4.2, # 'V'
	-0.9, # 'W'
	-1.3, # 'Y'
	0.0,  # '*'
	# Asn, Gln; J = Leu/Ile;
	# "~" = Cys-S-S-Cys; U = Se-Cys;
	-3.5, -3.5, 4.15, 2.5, NA
	);
	names(coeff) = c("X", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
		"S", "T", "V", "W", "Y", "*", "B", "Z", "J", "~", "U");
	return(coeff);
}
hydrophobicity = function(aa) {
	if(length(aa[1]) == 3) {
		aaCodes = aaCodes(n = 0);
		id = match(aa, aaCodes$A3);
		tmpAA = aaCodes$A1[id];
	} else tmpAA = aa;
	hydroBase = aaHydrophobicity();
	id = match(tmpAA, names(hydroBase));
	hydro = hydroBase[id];
	return(hydro);
}
hydroGroups = function() {
	aa = groups.aa(1);
	tapply(names(aa), aa, function(aa) hydrophobicity(aa));
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
