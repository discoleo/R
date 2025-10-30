########################
###
### Leonard Mada
### [the one and only]
###
### Pubmed
### Tools: Search Engine
###
### draft v.0.2d


### Pubmed Tools
# - Search & Retrieve Data;


################

library(curl)
library(xml2)


# process the XML files
source("Pubmed.XML.R")

# Note:
# - Users are encouraged to read carefully
#   the Entrez documentation at:
#   https://www.ncbi.nlm.nih.gov/books/NBK25500/
#   and to provide help with the development of this Search Tool;

### Other packages:
# - Package rentrez: 'Entrez' in R
#   https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
# - Package puremoe: Pubmed Unified REtrieval for Multi-Output Exploration
#   https://cran.r-project.org/web/packages/puremoe/index.html


########################
########################

### Credentials

### eMail:
GetEMail0 = function() {
	eMail = "...";
	if(eMail == "...")
		stop("Please provide a valid eMail address,
		so that I do not get blamed when the PubMed server crashes!")
	return(eMail);
}

### App data:
GetAppData = function(email = NULL) {
	if(is.null(email)) {
		email = GetEMail();
	}
	opt = list(
		AppName = "etoolR",
		User  = "test",
		EMail = email
	);
	opt$Tool = paste0(opt$AppName, "/", opt$User);
	return(opt);
}
opt.Pubmed = GetAppData()

########################

dataBaseName = "pubmed";
baseUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";

GetCredentials = function() {
	s = paste0("email=", curl_escape(opt.Pubmed$EMail),
		"&tool=", curl_escape(opt.Pubmed$Tool));
	return(s);
}

##############
### Search ###
##############

### Pubmed Search
search.entrez = function(..., options=NULL, debug=TRUE) {
	query = list(...);
	if(length(query) == 0) {
		stop("Missing query!");
	}
	query = encodeQuery(query);
	query = paste0("term=", query);
	if(debug) print(query);
	url = paste0(baseUrl, "esearch.fcgi?db=", dataBaseName,
		"&usehistory=y", "&", query,
		GetSearchOptions(options), "&", GetCredentials());
	#
	con = curl(url)
	lns = readLines(con)
	lns = paste(lns, collapse="\n")
	close(con)
	return(lns)
}
encodeQuery = function(...) {
	query = list(...);
	if(length(query) == 1) query = query[[1]];
	fields = names(query);
	# Date:
	query = as.PubMed.date(query);
	query = lapply(query, function(s) {
		# lapply: needed to preserve inheritance for later;
		if(inherits(s, "QPubmed")) return(s);
		# s = strsplit(s, "\\s++", perl=TRUE);
		s = sapply(s, curl_escape);
		# s = paste(s, collapse="+");
		return(s);
	})
	if(is.null(fields)) {
		isSearch = rep(TRUE, length(query));
	} else {
		isSearch = (nchar(fields) == 0);
	}
	# Basic Search: Title & Abstract
	if(sum(isSearch) > 0) {
		fld_SEARCH = fieldsPubmed("SEARCH")$Field;
		query[isSearch] = sapply(query[isSearch], function(s) {
			if(inherits(s, "QPubmed")) return(s);
			paste0(s, fld_SEARCH);
		})
	}
	#
	isOther = ! isSearch;
	if(sum(isOther) > 0) {
		query[isOther] = sapply(which(isOther), function(id) {
			s = query[[id]];
			if(inherits(s, "QPubmed")) return(s);
			fld_SEARCH = fieldsPubmed(fields[id])$Field;
			paste0(s, fld_SEARCH);
		})
	}
	#
	if(length(query) > 1) {
		query = paste0(query, collapse="+AND+");
	}
	return(query);
}
as.PubMed.date = function(x, collapse = TRUE, reg.period = "[\\:]") {
	fields = names(x);
	# Date:
	idDate = which(tolower(fields) == "date");
	if(length(idDate) > 0) {
		PDAT = "[PDAT]";
		sDt  = x[idDate];
		notArray = sapply(sDt, length) < 2;
		isPeriod = grepl(reg.period, sDt);
		idPeriod = which(isPeriod & notArray);
		if(length(idPeriod) > 0) {
			# Process Date-Range:
			sDt = unlist(sDt[idPeriod]);
			sDt = strsplit(sDt, reg.period);
			sDt = lapply(sDt, function(x) {
				if(nchar(x[1]) == 4) {
					x[1] = paste0(x[1], "/01/01");
				}
				if(length(x) == 1) {
					x[2] = paste0(format(Sys.Date(), "%Y"), "/12/31");
				} else {
					nch2 = nchar(x[2]);
					if(nch2 == 4) {
						x[2] = paste0(x[2], "/12/31");
					}
				}
				sDt = curl_escape(x);
				sDt = paste0(sDt, PDAT, collapse="+%3A+");
				sDt = paste0("(", sDt, ")");
				class(sDt) = c("QPubmed", class(sDt));
				return(sDt);
			})
			# Update Period:
			idDate = idDate[idPeriod];
			for(id in seq_along(idDate)) {
				x[[idDate[id]]] = sDt[[id]];
			}
		}
		# Collapse:
		if(collapse) {
			sDt = x[idDate];
			sDt = collapse.dates(sDt);
			x[[idDate]] = NULL;
			x = c(x, sDt);
		}
	}
	return(x);
}
collapse.dates = function(x) {
	isNotQPub = sapply(x, function(x) ! inherits(x, "QPubmed"));
	idNotQPub = which(isNotQPub);
	if(length(idNotQPub) > 0) {
		tmp = x[idNotQPub];
		# Array with Multiple Dates:
		lenField = sapply(tmp, length);
		idArray  = which(lenField > 1);
		if(length(idArray) > 0) {
			tmp2 = sapply(tmp[idArray], as.date.SimplePubMed);
			tmp2 = unlist(tmp2);
			tmp2 = as.date.OrPubmed(tmp2);
			id0  = idNotQPub[idArray[1]];
			idRm = idNotQPub[idArray[-1]];
			x[[id0]] = tmp2;
			x[idRm]  = NULL;
			idNotQPub = idNotQPub[- idArray];
		}
		if(length(idNotQPub) > 0) {
			tmp = as.date.SimplePubMed(x[[idNotQPub]]);
			x[idNotQPub] = tmp;
		}
	}
	if(length(x) > 1) {
		x = as.date.OrPubmed(x);
	}
	return(x);
}
as.date.SimplePubMed = function(x, PDAT = fieldsPubmed("Date")$Field) {
	if(length(x) == 0) return(character(0));
	tmp = paste0('"', curl_escape(x), '"');
	tmp = paste0(tmp, PDAT);
	return(tmp);
}
as.date.OrPubmed = function(x) {
	sDt = paste0(x, collapse = "OR");
	sDt = paste0("(", sDt, ")");
	class(sDt) = c("QPubmed", class(sDt));
	return(sDt);
}
escapeHTML = function(s) {
	curl_escape(s);
}

### Fetch Results
search.entrez.fetch = function(key, nStart=0, max=0, type="Abstract",
		file=NULL, options=NULL, debug=TRUE) {
	if(missing(key)) stop("Key must be provided!");
	isKey = inherits(key, "Entrez");
	if(isKey) {
		sKey = c(key$Key, key$Key2);
	} else if(length(key) != 2) {
		stop("Invalid key!");
	} else sKey = key;
	#
	query = paste0("query_key=", sKey[[1]],
		"&webEnv=", sKey[[2]]);
	type = pmatch(type, c("Abstract", "ids"));
	if(type == -1) stop("Invalid type!")
	retType = if(type == 1) "abstract" else "uilist";
	#
	url = paste0(baseUrl, "efetch.fcgi?db=", dataBaseName,
		"&usehistory=y&", query, "&rettype=", retType, "&retmode=xml",
		GetSearchOptions(options), "&", GetCredentials());
	if(debug) print(url);
	### Limits
	if(nStart > 0) {
		url = paste0(url, "&retstart=", nStart);
	}
	if(max > 0) {
		url = paste0(url, "&retmax=", max);
	}
	#
	con = curl(url)
	doc = readLines(con)
	doc = paste(doc, collapse="\n")
	if( ! is.null(file)) {
		writeLines(doc, con=file);
	}
	close(con)
	return(doc)
}

queryOr = function(...) {
	query = list(...);
	nms   = names(query);
	if(is.null(nms)) {
		fld_SEARCH = fieldsPubmed("SEARCH")$Field;
		query = sapply(query, function(s) {
			paste0(escapeHTML(s), fld_SEARCH);
		})
	} else {
		stop("Not yet implemented!");
	}
	#
	if(length(query) > 1) {
		query = paste0("(", query, ")", collapse="+OR+");
		query = paste0("(", query, ")")
	}
	class(query) = c("QPubmed", class(query));
	return(query);
}

### ELink

# https://www.nlm.nih.gov/dataguide/eutilities/utilities.html#elink

# Cited by:
search.cited = function(PMID, ..., options = NULL, debug = TRUE) {
	q0 = "linkname=pubmed_pubmed_citedin";
	query = list(...);
	if(length(query) > 0) {
		query = encodeQuery(query);
		query = paste0("term=", query);
		query = paste0(q0, "&", query);
	} else query = q0;
	url = paste0(baseUrl, "elink.fcgi?dbfrom=", dataBaseName,
		"&db=", dataBaseName, "&id=", PMID,
		"&usehistory=y", "&", query,
		GetSearchOptions(options), "&", GetCredentials());
	#
	con = curl(url)
	lns = readLines(con)
	lns = paste(lns, collapse="\n")
	close(con)
	return(lns)
}


################
### Advanced ###
################

fieldsPubmed = function(opt = NULL) {
	fl = list(
		SEARCH 		= list(Name="Abstract/Title", Field="[TIAB]"),
		TIAB 		= list(Name="Abstract/Title", Field="[TIAB]"),
		TITLE 		= list(Name="Title", Field="[TITL]"),
		ABSTRACT 	= list(Name="Abstract", Field="[AB]"),
		AUTHOR 		= list(Name="Author", Field="[AUTH]"),
		JOURNAL 	= list(Name="Journal", Field="[journal]"),
		DATE 		= list(Name="Date", Field="[PDAT]"),
		TYPE		= list(Name="Article Type", Field="[PTYP]"),
		COI 		= list(Name="Conflict", Field="[COI]") # Conflict of Interest: COIS ?
	);
	if(is.null(opt)) return(fl);
	#
	opt = toupper(opt);
	idOpt = pmatch(opt, names(fl));
	if(idOpt == -1) stop("Invalid Field!");
	return(fl[[idOpt]]);
}

GetSearchOptions = function(options=NULL) {
		if(is.null(options) || nchar(options) == 0) {
			return("");
		}
		return(paste0("&", options));
}

PublicationType = function(type, caseInsensitive=TRUE) {
	pT = c("Adaptive Clinical Trial",
	"Address",
	"Autobiography",
	"Bibliography",
	"Biography",
	"Case Reports",
	"Classical Article",
	"Clinical Conference",
	"Clinical Study",
	"Clinical Trial",
	"Clinical Trial, Phase I",
	"Clinical Trial, Phase II",
	"Clinical Trial, Phase III",
	"Clinical Trial, Phase IV",
	"Clinical Trial Protocol",
	"Clinical Trial, Veterinary",
	"Collected Work",
	"Comment",
	"Comparative Study",
	"Congress",
	"Consensus Development Conference",
	"Consensus Development Conference, NIH",
	"Controlled Clinical Trial",
	"Corrected and Republished Article",
	"Dataset",
	"Dictionary",
	"Directory",
	"Duplicate Publication",
	"Editorial",
	"Electronic Supplementary Materials",
	"English Abstract",
	"Equivalence Trial",
	"Evaluation Study",
	"Expression of Concern",
	"Festschrift",
	"Government Publication",
	"Guideline",
	"Historical Article",
	"Interactive Tutorial",
	"Interview",
	"Introductory Journal Article",
	"Journal Article",
	"Lecture",
	"Legal Case",
	"Legislation",
	"Letter",
	"Meta-Analysis",
	"Multicenter Study",
	"News",
	"Newspaper Article",
	"Observational Study",
	"Observational Study, Veterinary",
	"Overall",
	"Patient Education Handout",
	"Periodical Index",
	"Personal Narrative",
	"Portrait",
	"Practice Guideline",
	"Preprint",
	"Pragmatic Clinical Trial",
	"Published Erratum",
	"Randomized Controlled Trial",
	"Randomized Controlled Trial, Veterinary",
	"Research Support, American Recovery and Reinvestment Act",
	"Research Support, N.I.H., Extramural",
	"Research Support, N.I.H., Intramural",
	"Research Support, Non-U.S. Gov't",
	"Research Support, U.S. Gov't, Non-P.H.S.",
	"Research Support, U.S. Gov't, P.H.S.",
	"Retracted Publication",
	"Retraction of Publication",
	"Review", # includes Systematic Review
	"Scientific Integrity Review",
	"Systematic Review",
	"Technical Report",
	"Twin Study",
	"Validation Study",
	"Video-Audio Media",
	"Webcast");
	#
	if(caseInsensitive) type = paste0("(?i)", type);
	sT = lapply(type, function(type) {
		isType = grepl(type, pT, perl=TRUE);
		pT[isType];
	})
	sT = unique(unlist(sT));
	if(length(sT) == 0) {
		stop("Wrong Publication Type!");
	}
	return(sT);
}

