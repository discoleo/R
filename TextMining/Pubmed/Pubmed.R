########################
###
### Leonard Mada
### [the one and only]
###
### Pubmed
### Tools: Search Engine
###
### draft v.0.1d


### Pubmed Tools
# - Search & Retrieve Data;


################

library(curl)
library(xml2)

# Note:
# - Users are encouraged to read carefully
#   the Entrez documentation at:
#   https://www.ncbi.nlm.nih.gov/books/NBK25500/
#   and to provide help with the development of this Search Tool;


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
sAppName = "etoolR";
sUser  = "test";
sEMail = GetEMail();
sTool  = paste0(sAppName, "/", sUser);

########################

dataBaseName = "pubmed";
baseUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";

GetCredentials = function() {
	s = paste0("email=", curl_escape(sEMail),
		"&tool=", curl_escape(sTool));
	return(s);
}

GetSearchOptions = function(options=NULL) {
		if(is.null(options) || nchar(options) == 0) {
			return("");
		}
		return(paste0("&", options));
}

### Pubmed Search
search.entrez = function(query, options=NULL) {
	if(missing(query)) {
		stop("Missing query!");
	}
	query = paste0("term=", query);
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
search.entrez.fetch = function(key, nStart=0, type="Abstract", options=NULL) {
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
	print(url)
	if(nStart > 0) {
		url = paste0(url, "&retstart=", nStart);
	}
	#
	con = curl(url)
	doc = readLines(con)
	doc = paste(doc, collapse="\n")
	close(con)
	return(doc)
}


###########
### XML ###
###########

### Parse

# Parse the Search Results:
parse.entrez = function(x) {
	xml = read_xml(x);
	as.xml.numeric = function(xpath) {
		as.numeric(xml_text(xml_find_all(xml, xpath)));
	}
	count = as.xml.numeric("/eSearchResult/Count");
	len   = as.xml.numeric("/eSearchResult/RetMax");
	npos  = as.xml.numeric("/eSearchResult/RetStart");
	key   = xml_text(xml_find_first(xml, "/eSearchResult/QueryKey"));
	key2  = xml_text(xml_find_first(xml, "/eSearchResult/WebEnv"));
	r = list(Count=count, Len=len, nStart=npos, Key=key, Key2=key2);
	class(r) = c("Entrez", class(r));
	return(r)
}
# Parse fetched IDs:
parse.entrez.ids = function(x, nStart=0) {
	if(inherits(nStart, "Entrez")) {
		nStart = nStart$nStart;
	}
	xml = read_xml(x);
	ns  = xml_find_all(xml, "/IdList/Id");
	sID = xml_text(ns);
	if(nStart > 0) {
		attr(sID, "nStart") = nStart;
	}
	return(sID);
}
# Count
count.entrez.ids = function(x) {
	isXML = inherits(x, "xml_document");
	xml = if(isXML) x else read_xml(x);
	n   = xml_find_num(xml, "count(/IdList/Id)");
	return(n);
}
count.entrez.abstract = function(x) {
	isXML = inherits(x, "xml_document");
	xml = if(isXML) x else read_xml(x);
	n   = xml_find_num(xml, "count(/PubmedArticleSet/PubmedArticle)");
	return(n);
}

# Extract Titles & Year:
extractTitles = function(x) {
	isXML = inherits(x, "xml_document");
	xml = if(isXML) x else read_xml(x);
	# assumes: only 1 Title:
	pred  = "count(./MedlineCitation/Article/ArticleDate/Year)";
	as.XP = function(path)
		paste0("/PubmedArticleSet/PubmedArticle[", pred, " > 0]/MedlineCitation/", path);
	PMID   = xml_find_all(xml, as.XP("PMID"));
	sTitle = xml_find_all(xml, as.XP("Article/ArticleTitle"));
	year   = xml_find_all(xml, as.XP("Article/ArticleDate/Year[1]"));
	r = data.frame(PMID = xml_text(PMID), Year = as.numeric(xml_text(year)), Title = xml_text(sTitle));
	#
	as.XP = function(path)
		paste0("/PubmedArticleSet/PubmedArticle[", pred, " = 0]/MedlineCitation/", path);
	PMID   = xml_find_all(xml, as.XP("PMID"));
	sTitle = xml_find_all(xml, as.XP("Article/ArticleTitle"));
	r2 = data.frame(PMID = xml_text(PMID), Year = NA, Title = xml_text(sTitle));
	r  = rbind(r, r2);
	return(r);
}

# DISASTER:
extractTitles.slowBeyondAnyHope = function(x) {
	# Note: libxml2 does NOT implement XPATH 2 !!!
	# [while XPATH 1.0 has massive shortcomings!]
	isXML = inherits(x, "xml_document");
	xml = if(isXML) x else read_xml(x);
	nodes = xml_find_all(xml, "/PubmedArticleSet/PubmedArticle");
	# only first 100:
	r = lapply(nodes[1:100], function(r) {
		PMID   = xml_find_first(r, "./MedlineCitation/PMID");
		PMID   = xml_text(PMID)[1];
		sTitle = xml_find_all(r, "./MedlineCitation/Article/ArticleTitle");
		sTitle = paste(xml_text(sTitle), collapse="\n");
		year   = xml_find_first(r, "./MedlineCitation/Article/ArticleDate/Year");
		year   = as.numeric(xml_text(year)[1]);
		#
		data.frame(PMID=PMID, Title=sTitle, Year=year);
	});
	r = do.call(rbind, r);
	return(r);
}

