####################
###
### Chemistry Tools
### Smiles Format
###
### Leonard Mada
###
### v.0.1c


### Chemistry Tools:

# - Encode / Parse Smiles codes:
#   only for peptides;

# Peptide => Smile
# Smile => Peptide

####################

### Helper Functions

# aaCodes moved to external file!
source("Chem.AA.R")

###################

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


# 1. Minkiewicz P, Iwaniak A, Darewicz M. Annotation of Peptide Structures Using SMILES
#    and Other Chemical Codes - Practical Solutions. Molecules. 2017 Nov 27;22(12):2075.
#    https://doi.org/10.3390/molecules22122075
#    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6149970/
# 2. SMILES Tutorial
#    https://www.daylight.com/meetings/summerschool98/course/dave/smiles-intro.html

getAA = function(n=1) {
	if( ! any(n == c(1,3))) stop("Unsupported codes!");
	# [C@@H] = [C@@]([H]): should be equivalent;
	AA.smi = list(
		Gly = "NCC(=O)",
		Ala = "N[C@@H](C)C(=O)",
		# - OH, - SH
		Ser = "N[C@@H](CO)C(=O)",
		Thr = "N[C@@H]([C@H](O)C)C(=O)",
		Cys = "N[C@@H](CS)C(=O)",
		Met = "N[C@@H](CCSC)C(=O)",
		# Branched
		Val = "N[C@@H](C(C)C)C(=O)",
		Leu = "N[C@@H](CC(C)C)C(=O)",
		Ile = "N[C@@H]([C@H](CC)C)C(=O)",
		# Aromatic
		Phe = "N[C@@H](Cc1ccccc1)C(=O)",
		Tyr = "N[C@@H](Cc1ccc(O)cc1)C(=O)",
		Trp = "N[C@@H](CC(=CN2)C1=C2C=CC=C1)C(=O)",
		His = "N[C@@H](CC1=CN=C-N1)C(=O)",
		# Acidic
		Asp = "N[C@@H](CC(=O)O)C(=O)",
		Asn = "N[C@@H](CC(=O)N)C(=O)",
		Glu = "N[C@@H](CCC(=O)O)C(=O)",
		Gln = "N[C@@H](CCC(=O)N)C(=O)",
		#
		Pro = "N1[C@@H](CCC1)C(=O)",
		# Basic
		Arg = "N[C@@H](CCCNC(=N)N)C(=O)",
		Lys = "N[C@@H](CCCCN)C(=O)"
	);
	if(n == 1) {
		stdAA = aaCodes(n=1);
		codeAA  = stdAA[names(AA.smi)];
		names(AA.smi) = codeAA;
	}
	return(AA.smi);
}

