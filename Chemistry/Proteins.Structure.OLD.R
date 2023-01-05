

### OLD Code


### [OLD]
# - based on the initial code in:
#   https://github.com/AusteKan/Fiscore/blob/main/Fiscore/R/PDB_prepare.R;
features.pdb.old = function(pdb) {
	pdb_file = pdb;
	pdb.nms = attributes(pdb_file)$names;
	feature_list = list();
	if(is.null(pdb.nms)) return(feature_list);
	
	# Helix:
	if("helix" %in% pdb.nms) {
		# test if attribute is not NULL:
		if( ! is.null(pdb_file$helix$start)) {
			helix_df = as.data.frame(pdb_file$helix);
			type = as.vector(helix_df$type);
			helix_df$Type = dplyr::case_when(
				type == 1 ~ 'Right-handed alpha helix',
				type == 2 ~ 'Right-handed omega helix',
				type == 3 ~ 'Right-handed pi helix',
				type == 4 ~ 'Right-handed gamma helix',
				type == 5 ~ 'Right-handed 310 helix',
				type == 6 ~ 'Left-handed alpha helix',
				type == 7 ~ 'Left-handed omega helix',
				type == 8 ~ 'Left-handed gamma helix',
				type == 9 ~ '27 ribbon/helix helix',
				type == 10 ~ 'Polyproline helix',
				TRUE ~ as.character(type));
			
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
			sheet_df$Type = dplyr::case_when(
				type == 0 ~ 'Parallel sheet',
				type == 1 ~ 'Parallel sheet',
				type == -1 ~ 'Antiparalel sheet',
				TRUE ~ as.character(type));
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
