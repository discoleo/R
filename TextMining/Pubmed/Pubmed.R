########################
###
### Leonard Mada
### [the one and only]
###
### Pubmed
### Tools: Search Engine
###
### draft v.0.2c


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
# - other packages: rentrez
#   https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html


########################
########################

### Credentials

### eMail:
GetEMail = function() {
	eMail = "...";
	if(eMail == "...")
		stop("Please provide a valid eMail address,
		so that I do not get blamed when the Pubmed server crashes!")
	return(eMail);
}

### App data:
opt.Pubmed = list(
	AppName = "etoolR",
	User  = "test",
	EMail = GetEMail()
);
opt.Pubmed$Tool = paste0(opt.Pubmed$AppName, "/", opt.Pubmed$User);

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
	query  = list(...);
	if(length(query) == 1) query = query[[1]];
	fields = names(query);
	query  = lapply(query, function(s) {
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

################
### Advanced ###
################

fieldsPubmed = function(opt = NULL) {
	fl = list(
		SEARCH 		= list(Name="Abstract/Title", Field="[TIAB]"),
		TITLE 		= list(Name="Title", Field="[TITL]"),
		ABSTRACT 	= list(Name="Abstract", Field="[AB]"),
		AUTHOR 		= list(Name="Author", Field="[AUTH]"),
		JOURNAL 	= list(Name="Journal", Field="[jour]"),
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
	"Review",
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

