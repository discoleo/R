########################
###
### Leonard Mada
### [the one and only]
###
### Pubmed
### Tools: Search Engine
###
### draft v.0.2a


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

GetSearchOptions = function(options=NULL) {
		if(is.null(options) || nchar(options) == 0) {
			return("");
		}
		return(paste0("&", options));
}

##############
### Search ###
##############

### Pubmed Search
search.entrez = function(query, options=NULL, debug=TRUE) {
	if(missing(query)) {
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
encodeQuery = function(query) {
	fields = names(query);
	query  = sapply(query, function(s) {
		s = strsplit(s, "\\s++", perl=TRUE);
		s = sapply(s, curl_escape);
		s = paste(s, collapse="+");
		return(s);
	})
	if(is.null(fields)) {
		fld_SEARCH = fieldsPubmed("SEARCH")$Field;
		query = sapply(query, function(s) {
			paste0(s, fld_SEARCH);
		})
	}
	if(length(query) > 1) {
		query = paste0(query, collapse="+AND+");
	}
	return(query);
}

### Fetch Results
search.entrez.fetch = function(key, nStart=0, type="Abstract", max=0,
		options=NULL, debug=TRUE) {
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
	close(con)
	return(doc)
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
	idOpt = pmatch(opt, names(fl));
	if(idOpt == -1) stop("Invalid Field!");
	return(fl[[idOpt]]);
}

