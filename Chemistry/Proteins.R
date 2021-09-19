


library(Rpdb)


### AA
AA = c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
	"PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL");

### Other

# "BMA", "MAN", "NAG", 


####################

### Helper Functions

extract.aa = function(x, info=TRUE) {
	# AA
	aa = unique(x$atoms[ , c("chainid", "resname", "resid")])

	# remove H2O
	isWater = (aa$resname == "HOH");
	if(info) print(table(isWater));
	aa = aa[ ! isWater, ];
	if(info) print(table(aa$resname));
	# true AA:
	unqAA = unique(aa$resname);
	isAA = (unqAA %in% AA);
	notAA = unqAA[ ! isAA];
	isNotAA = (aa$resname %in% notAA);
	if(info) {
		print(table(aa$resname[isNotAA]))
		print(table(aa$resname[isNotAA], aa$chainid[isNotAA]))
	}
	aa = aa[ ! isNotAA, ]
	return(aa);
}

####################

x = read.pdb(file.choose())


### Extract AA Sequence

### AA
aa = extract.aa(x);

### AAs per Chain
table(aa$resname, aa$chainid)

### Chain Length
tapply(aa$chainid, aa$chainid, length)

