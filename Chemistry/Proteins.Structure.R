


# this file:
# source("Proteins.Structure.R")


###########################

### Fi-Score
# - based on package Fiscore:
#   https://github.com/AusteKan/Fiscore/blob/main/Fiscore/R/PDB_prepare.R;
# - specific code has been extracted as individual helper-functions,
#   see functions below;
# - the initial code has been refactored substantially;

source("Proteins.Structure.FiScore.R")


########################

# Open & Clean PDB file:
read.pdb = function(file_name, lim.protein=5, ...) {
	pdb_file_temp = bio3d::read.pdb(file_name, ...);
	# clean file to remove terminal residues and ligand data
	pdb_file_temp = bio3d::clean.pdb(pdb_file_temp,
		consecutive = TRUE, force.renumber = FALSE, fix.chain = FALSE, fix.aa = TRUE,
		rm.wat = TRUE, rm.lig = TRUE, rm.h = FALSE, verbose = FALSE);
	
	# Warnings:
	pdb.seq = bio3d::pdbseq(pdb_file_temp);
	if(length(pdb.seq) == 0) { stop("The file has no amino acid residues!"); }
	if(lim.protein > 0 && length(pdb.seq) <= lim.protein) {
		stop("This is a peptide and not protein!");
	}
	return(pdb_file_temp);
}

### Structure

as.type.helix = function(type) {
	types  = unique(type);
	isType = types %in% c(1:10);
	types.id = c(1:10, types[ ! isType]);
	type.nms = c(
		'Right-handed alpha helix',
		'Right-handed omega helix',
		'Right-handed pi helix',
		'Right-handed gamma helix',
		'Right-handed 310 helix',
		'Left-handed alpha helix',
		'Left-handed omega helix',
		'Left-handed gamma helix',
		'27 ribbon/helix helix',
		'Polyproline helix',
		paste0("Other: ", types[ ! isType])
	)
	type = factor(type, levels=types.id);
	levels(type) = type.nms;
	return(type);
}
as.type.sheet = function(type) {
	types  = unique(type);
	isType = types %in% c(-1,0,1);
	types.id = c(0, 1, -1, types[ ! isType]);
	type.nms = c(
		'Parallel sheet', # TODO: separate level "Start sheet"?
		'Parallel sheet',
		'Antiparalel sheet',
		paste0("Other: ", types[ ! isType])
	)
	type = factor(type, levels=types.id);
	levels(type) = type.nms;
	return(type);
}

features.pdb = function(pdb) {
	pdb_file = pdb;
	pdb.nms = attributes(pdb_file)$names;
	feature_list = list();
	if(is.null(pdb.nms)) return(feature_list);
	
	# Reference:
	# https://www.wwpdb.org/documentation/file-format-content/format23/sect5.html
	#
	### TYPE OF HELIX
	#
	#                TYPE OF HELIX          CLASS NUMBER
	#                                       (COLUMNS 39 - 40)
	#      ---------------------------------------------------
	#      Right-handed alpha (default)       1
	#      Right-handed omega                 2
	#      Right-handed pi                    3
	#      Right-handed gamma                 4
	#      Right-handed 310                   5
	#      Left-handed alpha                  6
	#      Left-handed omega                  7
	#      Left-handed gamma                  8
	#      27 ribbon/helix                    9
	#      Polyproline                       10
	#
	#
	### TYPE OF SHEET
	#
	# The sense indicates whether strand n is parallel (sense = 1)
	# or anti-parallel (sense = -1) to strand n-1.
	# Sense is equal to zero (0) for the first strand of a sheet.
	#
	### TURNS
	#
	# Turns include those sets of residues which form beta turns, i.e.,
	# have a hydrogen bond linking (C- O)[i] to (N-H)[i + 3].
	# Turns which link residue i to i+2 (gamma-bends) may also be included.
	# Others may also be classified as turns.
	
	# Helix:
	if("helix" %in% pdb.nms) {
		# test if attribute is not NULL:
		if( ! is.null(pdb_file$helix$start)) {
			helix_df = as.data.frame(pdb_file$helix);
			type = as.vector(helix_df$type);
			helix_df$Type = as.type.helix(type);
			feature_list[["helix"]] = helix_df;
		}
	}
	# Sheets:
	if("sheet" %in% pdb.nms) {
		# prepare sheet data frame;
		# test if attribute is not NULL;
		if( ! is.null(pdb_file$sheet$start)) {
			sheet_df = as.data.frame(pdb_file$sheet);
			type = as.vector(sheet_df$sense);
			sheet_df$Type = as.type.sheet(type);
			feature_list[["sheet"]] = sheet_df;
		}
	}
	# Turns:
	if("turn" %in% pdb.nms){
		# prepare turn data frame;
		# test if attribute is not NULL;
		if( ! is.null(pdb_file$turn$start)) {
			turn_df = as.data.frame(pdb_file$turn);
			type = as.vector(turn_df$turnId);
			turn_df$Type = type;
			feature_list[["turn"]] = turn_df;
		}
	}
	
	return(feature_list);
}

torsions.pdb = function(pdb) {
	### Torsion angles
	torsion_angles = bio3d::torsion.pdb(pdb);
	# Extract torsion angle table:
	pdb_df = torsion_angles$tbl;
	# - leave rows that contain full dihedral angle visualization;
	# NOTE: terminal residues do not contain all of the angles;
	isComplete = stats::complete.cases(pdb_df[ , c("phi","psi")]);
	pdb_df = pdb_df[isComplete, , drop=FALSE];
	
	# Extract residue numbers:
	# Note: the vectorised code should be faster;
	df_resno = as.numeric(stringr::str_extract(rownames(pdb_df), "[0-9]{1,}"));
	# df_resno = as.numeric(sapply(rownames(pdb_df), function(x) {
	#	stringr::str_extract(x, "[0-9]{1,4}");
	# }));
	
	# Extract residue names:
	df_res = as.vector(stringr::str_extract(rownames(pdb_df), "[A-Z]{3}"));
	# df_res = as.vector(sapply(rownames(pdb_df), function(x) {
	#	stringr::str_extract(x, "[A-Z]{3}");
	# }));
	
	# Construct the data.frame to contain residue names and numbers
	pdb_df = cbind.data.frame(pdb_df, df_resno);
	pdb_df = cbind.data.frame(pdb_df, df_res);
	attr(pdb_df, "complete") = isComplete;
	
	return(pdb_df);
}

BFactor.pdb = function(pdb, torsions=NULL, normalize=TRUE) {
	### Torsion angles:
	if(is.null(torsions)) { pdb_df = torsions.pdb(pdb); }
	else pdb_df = torsions;
	
	### B-factor extraction:
	# Full data frame:
	# - includes dihedral angles coordinates and residue info;
	# Extracting B factor information for C alpha atom to match dihedral angles
	isComplete = attr(pdb_df, "complete");
	idCA = which(pdb$atom$elety == "CA")[isComplete];
		# & (pdb$atom$resno %in% pdb_df$"df_resno") );
	pdb_b = pdb$atom[idCA, c("resno","resid","b")];
	# Adding B factor information
	if(nrow(pdb_df) != nrow(pdb_b)) {
		# TODO: ugly ERROR;
		# - should be corrected by using isComplete;
		# names(pdb_b)[3] = "B_factor";
		# pdb_df = merge(pdb_df, pdb_b[, c("resno", "B_factor")],
		#	by.x = "df_resno", by.y = "resno");
	} else {
		pdb_df$B_factor = pdb_b$b;
	}
	if(all(pdb_df$B_factor == 0)) {
		warning("All B-factors are 0 and the analysis will be limited");
	}
	
	### B-factor normalization
	# - added as norm column;
	if(normalize) {
		pdb_df$B_normalised = MINMAX_normalisation_func(pdb_df$B_factor);
	}
	return(pdb_df);
}

### Fi-score
FiScore = function(pdb_df) {
	# - Generate Fi-score per residue and store in the data.frame;
	# - Calculate Fi-score for the whole protein and individual aa;
	
	### phi/psi normalization:
	# normalization is only scaled by SD
	psi_SD = stats::sd(pdb_df$psi);
	phi_SD = stats::sd(pdb_df$phi);
	div    = psi_SD * phi_SD;
	fi_score = pdb_df$phi * pdb_df$psi * pdb_df$B_normalised / div;
	
	# OLD code:
	# for(index in seq(nrow(pdb_df))) {
	#	fi_score = c(fi_score,
	#		pdb_df[index,'phi'] * pdb_df[index,'psi'] * pdb_df[index,'B_normalised'] / div);
	# }
	return(fi_score);
}

### Helper functions for the analysis
	
### MIN-MAX normalisation based on the input array
# input = numeric array;
# returns normalised array values;
MINMAX_normalisation_func = function(array) {
	# check for cases where all B-factor values are 0;
	rg = range(array);
	rg.diff = rg[2] - rg[1];
	if(rg.diff == 0) {
		if(rg[2] == 0) return (0);
		# TODO: rg[2] != 0;
	}
		
	return ((array - rg[1])/rg.diff);
}

