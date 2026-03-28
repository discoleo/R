
### Tools to Copy Packages between versions of R

### Note:
# - Code copied from package "installr";
# - Some Refactoring: more elegant/robust code;

warning("Not yet tested!");


### Copy Packages
# TODO: refactor;
copy.packages.between.libraries =
function (from, to, ask = FALSE, keep_old = TRUE, do_NOT_override_packages_in_new_R = TRUE) 
{
    installed_R_folders <- get.RFolders();
    installed_R_folders_TABLE <- data.frame(R_version = names(installed_R_folders), 
        Folder = installed_R_folders)
    if (ask) {
        ask_again <- T
        while (ask_again) {
            ss_R_folder_from <- ask.user.for.a.row(installed_R_folders_TABLE, 
                "From: Choose an R version/folder from which to copy the package library", 
                "FROM: Write the row number of the R version/folder FROM which to copy (or move) the package library (and then press enter):\n")
            cat("\nThank you\n")
            ss_R_folder_to <- ask.user.for.a.row(installed_R_folders_TABLE, 
                "To: Choose an R version/folder from INTO which the packages from your old R installations will be copied TO", 
                "TO: Write the row number of the R version/folder INTO which to copy (or move) the package library (and then press enter):\n")
            cat("\nThank you\n")
            from <- installed_R_folders[ss_R_folder_from];
            to   <- installed_R_folders[ss_R_folder_to];
            DECISION_text <- paste("You've chosen to move your packages from: ", 
                from, "  to: ", to)
            ask_again_12 <- ask.user.for.a.row(c("yes", "no"), 
                "Are you sure?", paste(DECISION_text, " - Is this you final decision? \n(for 'yes' press 1, and for 'no' press 2)\n"))
            ask_again <- ifelse(ask_again_12 == 1, F, T)
        }
    }
    if (missing(from)) 
        from <- installed_R_folders[2]
    if (missing(to)) 
        to <- installed_R_folders[1]
	# Path to Packages:
    from_library <- file.path(from, "library")
    to_library   <- file.path(to, "library")
    packages_in_the_from_library <- list.files(from_library)
    packages_in_the_to_library <- list.files(to_library)
    packages_to_NOT_move <- unname(installed.packages(priority = "high")[, 
        "Package"])
	# Skip packages:
    if (do_NOT_override_packages_in_new_R) 
        packages_to_NOT_move <- c(packages_to_NOT_move, packages_in_the_to_library)
    ss_packages_to_NOT_move_from <- packages_in_the_from_library %in% 
        packages_to_NOT_move
    ss_packages_to_YES_move_from <- !ss_packages_to_NOT_move_from
    pkgsToMove = packages_in_the_from_library[ss_packages_to_YES_move_from]
    paths_of_packages_to_copy_from <- file.path(from_library, 
        pkgsToMove)
    if (length(pkgsToMove) == 0) {
        cat("No packages to copy.  Goodbye :) \n")
        return(F)
    }
	# Start Copying/Moving:
	cat("-----------------------", "\n")
	cat("I am now copying ", length(pkgsToMove), " packages\n".
		"  from: ", from_library, "\n",
		"  into: ", to_library)
	cat("-----------------------", "\n")
	flush.console()
	# pkgsToMove # ???
	folders.copied = file.copy(
		from = paths_of_packages_to_copy_from, 
		to   = to_library, overwrite = !do_NOT_override_packages_in_new_R, 
		recursive = TRUE);
	cat("=====================", "\n")
    cat("Done. Finished copying all packages to the new location\n")
	flush.console();
    if (!keep_old) {
        cat("Next: we will remove the packages from the old R installation ('FROM') \n")
        deleted_packages <- unlink(paths_of_packages_to_copy_from, 
            recursive = TRUE)
        cat("Done. The old packages were deleted.\n")
    }
    return(TRUE)
}

### Helper Functions

### R-Folders
get.RFolders = function (sort_by_version = TRUE, add_version_to_name = TRUE) 
{
	path = head(strsplit(R.home(), "/|\\\\")[[1]], -1);
    pathParent  = paste(path, collapse = "/");
    itemsParent = list.files(pathParent)
    foldersR    = file.path(pathParent, itemsParent);
    versionsR   = sapply(foldersR, get.RVersion);
    if (all(is.na(versionsR))) {
        warning("Could not find any R installation on your system.",
			"(You might have installed your R version on 'c:\\R' without sub folders...");
        return(NULL);
    }
    isFolderR = ! is.na(versionsR);
    foldersR  = foldersR[isFolderR];
    versionsR = versionsR[isFolderR];
	versionsR = unlist(versionsR);
	orderR    = order.versions(versionsR);
    if (add_version_to_name)
        names(foldersR) = versionsR;
    if (sort_by_version) {
        foldersR = foldersR[orderR];
    }
    return(foldersR)
}

### R-Versions
get.RVersion = function(folder) 
{
	files = list.files(folder);
	files = gsub("patched|revised", "", files);
	nmsR  = grep("README.R-[0-9]+.[0-9]+.[0-9]+$", files);
	if (length(nmsR) == 0) 
		return(NA);
	nameREADME = files[nmsR];
	versions = substr(nameREADME, start = 10, stop = nchar(nameREADME));
	versions = sort(versions);
	versions[length(versions)];
}

order.versions = function(x, decreasing = TRUE) {
    listVersions = as.numeric.version(x);
	listVersions$decreasing = decreasing;
	id = do.call(order, listVersions);
	return(id);
}
# x = String with Versions;
as.numeric.version = function(x) 
{
	if (length(x) > 1) {
		# TODO: pad for any length;
        tmp = sapply(x, as.numeric.version1);
		tmp = t(tmp);
		tmp = list(tmp[,1], tmp[,2], tmp[,3]);
	} else {
        tmp = as.numeric.version1(x);
		tmp = as.list(tmp);
	}
	return(tmp);
}
# x = String with Version;
as.numeric.version1 = function(x) 
{
	sVersion = strsplit(x, "\\.")[[1]];
	iVersion = as.numeric(sVersion);
	len = length(iVersion);
	if(len < 3) iVersion = c(iVersion, rep(0, 3 - len));
	longVersion = iVersion[1:3];
    longVersion;
}
