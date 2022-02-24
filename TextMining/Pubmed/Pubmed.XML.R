########################
###
### Leonard Mada
### [the one and only]
###
### Pubmed
### XML Tools
###
### draft v.0.1d


### XML Tools
# - Extract data from Pubmed XML files;


#################

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

extractTitles.hack = function(x, max=0, debug=TRUE) {
	# DISASTER:
	# - lapply on actual nodes: slow Beyond Any Hope !!!
	# Note: libxml2 does NOT implement XPATH 2 !!!
	# [while XPATH 1.0 has massive shortcomings!]
	isXML = inherits(x, "xml_document");
	xml = if(isXML) x else read_xml(x);
	nodes = xml_find_all(xml, "/PubmedArticleSet/PubmedArticle");
	len = length(nodes);
	# only first max-nodes:
	if(max > 0) len = min(max, len);
	r = data.frame(PMID = numeric(len), Year = numeric(len), Title = character(len));
	for(id in seq(len)) {
		nd = read_xml(as.character(nodes[[id]]));
		PMID   = xml_find_first(nd, "./MedlineCitation/PMID");
		PMID   = xml_text(PMID)[1];
		sTitle = xml_find_all(nd, "./MedlineCitation/Article/ArticleTitle");
		sTitle = paste(xml_text(sTitle), collapse="\n");
		year   = xml_find_first(nd, "./MedlineCitation/Article/ArticleDate/Year");
		year   = as.numeric(xml_text(year)[1]);
		#
		r[id, ] = data.frame(PMID=PMID, Year=year, Title=sTitle);
		if(debug && id %% 100 == 1) print(id);
	}
	# r = do.call(rbind, r);
	return(r);
}

### Authors

# Extract Authors
extractAuthors = function(x, max=3, collapse=";", filter=NULL) {
	isXML = inherits(x, "xml_document");
	xml = if(isXML) x else read_xml(x);
	#
	base = "/PubmedArticleSet/PubmedArticle";
	if( ! is.null(filter)) {
		pred = "count(./MedlineCitation/Article/AuthorList/Author)";
		if(is.logical(filter)) {
			base = if(filter) paste0(base, "[", pred, "> 0]") else base;
		} else if(is.numeric(filter)) {
			base = paste0(base, "[", pred, ">= ", filter, "]");
		} else {
			stop("Filter Option NOT supported!");
		}
	}
	nodes = xml_find_all(xml, base);
	len = length(nodes);
	nA  = lapply(seq(nodes), function(id) {
		nd = read_xml(as.character(nodes[[id]]));
		PMID  = xml_find_first(nd, "./MedlineCitation/PMID");
		PMID  = xml_text(PMID);
		count = xml_find_num(nd, "count(./MedlineCitation/Article/AuthorList/Author)");
		ndAuthors = xml_find_all(nd, "./MedlineCitation/Article/AuthorList/Author");
		if(max > 0) ndAuthors = ndAuthors[min(max, length(ndAuthors))];
		# TODO: separator;
		sAuth = paste0(xml_text(ndAuthors), collapse = collapse);
		data.frame(PMID = PMID, Count = count, Authors = sAuth);
	});
	r = do.call(rbind, nA);
	return(r)
}
# Extract Authors
countAuthors = function(x) {
	isXML = inherits(x, "xml_document");
	xml = if(isXML) x else read_xml(x);
	#
	base = "/PubmedArticleSet/PubmedArticle";
	nodes = xml_find_all(xml, base);
	len = length(nodes);
	nA  = lapply(seq(nodes), function(id) {
		nd = read_xml(as.character(nodes[[id]]));
		PMID  = xml_find_first(nd, "./MedlineCitation/PMID");
		PMID  = xml_text(PMID);
		count = xml_find_num(nd, "count(./MedlineCitation/Article/AuthorList/Author)");
		data.frame(PMID = PMID, Count = count);
	});
	r = do.call(rbind, nA);
	return(r)
}

