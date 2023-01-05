
### Note:
# - the code is from the R package Fiscore;
# - includes refactored code for the PDB_prepare function:
#   https://github.com/AusteKan/Fiscore/blob/main/Fiscore/R/PDB_prepare.R;
# - specific code has been moved to helper functions,
#   and is available in the file: Proteins.Structure.R;
# - dependencies on package dplyr have been removed;
# - new code should be optimized compared to initial code;
# - new code: NOT tested for bugs !!!
# - see GitHub issues:
#   1) ...
#   2) https://github.com/AusteKan/Fiscore/issues/1


#' @title PDB_prepare
#'
#' @description Function to prepare a PDB file after it was pre-processed to generate Fi-score and normalised B factor values as well as secondary structure designations
#'
#' @param file_name PDB file name to load that was split into chains, e.g. '6KZ5_A.pdb'

#' @return  returns a processed data frame with Fi-score 'Fi_score', normalised B factor values 'B_normalised' and secondary structure designations
#' @ImportFrom bio3d read.pdb
#' @ImportFrom bio3d clean.pdb
#' @ImportFrom bio3d torsion.pdb
#' @ImportFrom stringr str_extract
#' @ImportFrom stats sd
#' @ImportFrom stats complete.cases
#' @ImportFrom dplyr case_when
#' @ImportFrom dplyr mutate
#' @export
#' @examples
#' path_to_processed_PDB<- system.file("extdata", "3nf5_A.pdb", package="Fiscore")
#' # you can call PDB_prepare with the set path
#' head(PDB_prepare(path_to_processed_PDB))
PDB_prepare = function(file_name, lim.protein = 5) {
	### Fi-score scoring and PDB file processing
	
	# Prepare PDB file ------
	if(inherits(file_name, "pdb")) {
		pdb_file_temp = file_name;
	} else {
		pdb_file_temp = bio3d::read.pdb(file_name);
	}
	# clean file to remove terminal residues and ligand data
	pdb_file_temp = bio3d::clean.pdb(pdb_file_temp,
		consecutive = TRUE, force.renumber = FALSE, fix.chain = FALSE, fix.aa = TRUE,
		rm.wat = TRUE, rm.lig = TRUE, rm.h = FALSE, verbose = FALSE);
	
	# Warnings:
	pdb.seq = bio3d::pdbseq(pdb_file_temp);
	if(length(pdb.seq) == 0) { stop("The file has no amino acid residues"); }
	if(lim.protein > 0 && length(pdb.seq) <= lim.protein) {
		stop("This is a peptide and not protein");
	}
	
	### Torsion angles -------------
	pdb_df = torsions.pdb(pdb_file_temp);
	
	### B-factor extraction --------
	# Full data frame:
	# - includes dihedral angles coordinates and residue info;
	pdb_df = BFactor.pdb(pdb_file_temp, torsions = pdb_df, normalize = TRUE);
	
	
	### Fi-score --------
	# - Generate Fi-score per residue and store in the data.frame;
	# - Calculate fi_score for the whole protein and individual aa;
	pdb_df$Fi_score = FiScore(pdb_df);
	
	# Incorporate information for C[a] to indicate what secondary structure element it belongs to;
	Type_vals = rep("NA", nrow(pdb_df));
	pdb_df$Type = Type_vals;
	
	# Extract specific features
	feature_list = features.pdb(pdb_file_temp);
	
	if(length(feature_list) != 0) {
		for(feature in feature_list) {
			# feature = data frame that contains available information on helix, sheet or turn;
			for(i in 1:nrow(feature)) {
				# set feature range
				isRange = (pdb_df$df_resno >= feature$start[i]) &
					(pdb_df$df_resno <= feature$end[i]);
				# feature$Type is factor!
				pdb_df$Type[isRange] = as.character(feature$Type[i]);
				# [OLD]
				# range = feature[i,"start"]:feature[i,"end"];
				# mutate column to add what type of the secondary structure the residue belongs to
				# pdb_df = dplyr::mutate(pdb_df,
				#	Type = ifelse((pdb_df$df_resno %in% range), feature[i,"Type"], pdb_df$Type));
			}
		}
	}
	
	return(pdb_df);
}
