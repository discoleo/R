########################
###
### Leonard Mada
### [the one and only]
###
### Pubmed
### XML Tools
###
### draft v.0.1i


### XML Tools
# - Extract data from Pubmed XML files;


#################

###########
### XML ###
###########

# library(xml2)


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

# ELink Result:
parse.elink = function(x) {
	xml = read_xml(x);
	as.xml.numeric = function(xpath) {
		as.numeric(xml_text(xml_find_all(xml, xpath)));
	}
	xpath = "/eLinkResult/LinkSet";
	idSrc = as.xml.numeric(paste0(xpath, "/IdList/Id"));
	idArt = as.xml.numeric(paste0(xpath, "/LinkSetDb/Link/Id"));
	r = list(ID.Src = idSrc, Count0 = length(idSrc), # Source
		ID = idArt, Count = length(idArt));
	class(r) = c("Entrez.ELink", class(r));
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

### Exctract: Title + Abstract
# n.authors = Number of authors;
# Note:
# - Initial code: very slow (non-usable) on moderately large data-sets;
# - Alternative code: reasonably fast on 8,000 records;
# - Authors: seems slower;
extract.abs = function(x, n.authors = 1, sep = "; ", verbose = TRUE) {
	isXML = inherits(x, "xml_document");
	xml = if(isXML) x else read_xml(x);
	ppp = "/PubmedArticleSet/PubmedArticle/MedlineCitation";
	xs0 = xml_find_all(xml, ppp, flatten = FALSE);
	if(verbose) {
		cat("Found", length(xs0), "records", "\n", sep = " ");
	}
	# Note: reasonable velocity with current version of xml2;
	# Tested only on xml with 200 records!
	# New code: tested on 8000 records: ~10 seconds!
	art = lapply(xs0, function(xn) {
		xn = read_xml(as.character(xn));
		PMID  = xml_find_first(xn, "./PMID");
		PMID  = xml_text(PMID);
		xArt  = xml_find_first(xn, "./Article");
		Title = xml_find_first(xArt, "./ArticleTitle");
		Title = xml_text(Title);
		Year  = xml_find_first(xArt, "./ArticleDate/Year[1]");
		Year  = as.numeric(xml_text(Year));
		Abs0  = xml_find_all(xArt,   "./Abstract");
		Abs0  = paste0(xml_text(Abs0), collapse = "\n");
		Journ = xml_find_first(xArt, "./Journal/ISOAbbreviation");
		Journ = xml_text(Journ);
		if(n.authors == 0) {
			# NO Authors
			sAuthor = NULL;
		} else if(n.authors == 1) {
			Author = xml_find_first(xArt, "./AuthorList/Author");
			Author = paste(
				xml_text(xml_find_first(Author, "./LastName")),
				xml_text(xml_find_first(Author, "./Initials")), sep = " ");
		} else {
			# xAuthor = xml_find_first(xArt,  "./AuthorList");
			# Author  = xml_find_all(xAuthor, "./Author", flatten = FALSE);
			Author = xml_find_all(xArt, "./AuthorList/Author", flatten = FALSE);
			if(n.authors > 1) {
				n = min(n.authors, length(Author));
				Author = Author[seq(n)];
			}
			Author = sapply(Author, function(x) {
				paste(
					xml_text(xml_find_first(x, "./LastName")),
					xml_text(xml_find_first(x, "./Initials")), sep = " ");
			})
			Author = paste(Author, collapse = sep);
		}
		r = data.frame(PMID = PMID, Year = Year,
			Title = Title, Abstract = Abs0,
			Journal = Journ, Authors = Author);
	})
	art = do.call(rbind, art);
	return(art)
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

### Abstract
extractAbstract = function(x, query=NULL, sep = ": ") {
	isXML = inherits(x, "xml_document");
	xml = if(isXML) x else read_xml(x);
	#
	xpBase = "/PubmedArticleSet/PubmedArticle/MedlineCitation";
	xpText = "./Article/Abstract/AbstractText";
	xpAbst = "./Article/Abstract";
	if(is.null(query)) {
		nodes = xml_find_all(xml, xpBase);
	} else {
		nodes = lapply(query, function(query) {
			xpBase = paste0(xpBase, "[./PMID/text() = '", query, "']");
			nodes  = xml_find_all(xml, xpBase);
		})
		nodes = unlist(nodes, recursive = FALSE);
	}
	len = length(nodes);
	nA  = lapply(seq(len), function(id) {
		nd = read_xml(as.character(nodes[[id]]));
		PMID  = xml_find_first(nd, "./PMID");
		PMID  = xml_text(PMID);
		ndTxt = xml_find_all(nd, xpText);
		if(length(ndTxt) == 0) {
			sTxt = "";
			nd = xml_find_all(nd, xpAbst);
			if(length(nd) > 0) sTxt = xml_text(nd);
		} else {
			sTxt = xml_text(ndTxt);
			# Labels
			lbl = sapply(ndTxt, xml_attr, attr="Label");
			hasLbl = ! is.na(lbl);
			if(sum(hasLbl) > 0) {
				sTxt[hasLbl] = paste(lbl[hasLbl], sTxt[hasLbl], sep=sep);
			}
		}
		data.frame(PMID = PMID, Abstract = sTxt);
	});
	r = do.call(rbind, nA);
	return(r)
}


### Authors

# Extract Authors
extractAuthors = function(x, max=3, collapse=";\n", filter=NULL) {
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
		if(max > 0) ndAuthors = ndAuthors[seq(min(max, length(ndAuthors)))];
		# TODO: improve separator;
		sAuth = sapply(ndAuthors, function(ndAuthors) {
			paste0(xml_text(xml_children(ndAuthors)), collapse = ". ");
		});
		sAuth = paste0(sAuth, collapse = collapse);
		data.frame(PMID = PMID, Count = count, Authors = sAuth);
	});
	r = do.call(rbind, nA);
	return(r)
}
# Count Authors
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

### Journal
extractJournal = function(x, type="Abbreviated") {
	type = pmatch(type, c("Abbreviated", "Full"));
	if(is.na(type)) stop("Invalid type!");
	isXML = inherits(x, "xml_document");
	xml = if(isXML) x else read_xml(x);
	#
	xpBase = "/PubmedArticleSet/PubmedArticle/MedlineCitation";
	xpJrnl = if(type == 1) "./Article/Journal/ISOAbbreviation"
		else "./Article/Journal/Title";
	nodes = xml_find_all(xml, xpBase);
	len = length(nodes);
	lA  = lapply(seq(len), function(id) {
		nd = read_xml(as.character(nodes[[id]]));
		PMID  = xml_find_first(nd, "./PMID");
		PMID  = xml_text(PMID);
		sJrnl = xml_text(xml_find_first(nd, xpJrnl));
		data.frame(PMID = PMID, Journal = sJrnl);
	});
	r = do.call(rbind, lA);
	return(r)
}

### Language
extractLanguage = function(x) {
	isXML = inherits(x, "xml_document");
	xml = if(isXML) x else read_xml(x);
	#
	xpBase = "/PubmedArticleSet/PubmedArticle/MedlineCitation";
	xpLang = "./Article/Language";
	nodes = xml_find_all(xml, xpBase);
	len = length(nodes);
	lA  = lapply(seq(len), function(id) {
		nd = read_xml(as.character(nodes[[id]]));
		PMID  = xml_find_first(nd, "./PMID");
		PMID  = xml_text(PMID);
		sLang = xml_text(xml_find_first(nd, xpLang));
		data.frame(PMID = PMID, Language = sLang);
	});
	r = do.call(rbind, lA);
	return(r)
}

### Specified Node
extractGeneric = function(x, type="KeywordList/Keyword", collapse=", ") {
	isXML = inherits(x, "xml_document");
	xml = if(isXML) x else read_xml(x);
	#
	xpBase = "/PubmedArticleSet/PubmedArticle/MedlineCitation";
	xpNode = paste0("./", type);
	nodes  = xml_find_all(xml, xpBase);
	len = length(nodes);
	lA  = lapply(seq(len), function(id) {
		nd = read_xml(as.character(nodes[[id]]));
		PMID  = xml_find_first(nd, "./PMID");
		PMID  = xml_text(PMID);
		sNode = xml_text(xml_find_all(nd, xpNode));
		sNode = paste(sNode, collapse=collapse);
		data.frame(PMID = PMID, Content = sNode);
	});
	r = do.call(rbind, lA);
	return(r)
}

extractKeywords = function(x, collapse=", ") {
	keys = extractGeneric(x, type="KeywordList/Keyword", collapse=collapse);
	hasNoKey = (nchar(keys$Content) == 0);
	xp2 = "MeshHeadingList/MeshHeading/DescriptorName";
	# - does NOT extract additional descriptors;
	keys$Content[hasNoKey] =
		extractGeneric(x, type=xp2, collapse=collapse)$Content[hasNoKey];
	return(keys)
}

