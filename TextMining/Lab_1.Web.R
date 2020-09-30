
######################
###
### Data Mining: Lab 1
### Extracting articles from Pubmed / BMC
###
### Leonard Mada
###

#################
### Libraries ###

# install.packages("rvest")

library(curl)
library(rvest)

library(tm)

###########

setwd("")

#############
### Intro ###

# Articles will be extracted from:
# Pubmed & BiomedCentral (BMC)


### Other Useful Sites/Packages:

### 1.) Library: easyPubMed
# https://cran.r-project.org/web/packages/easyPubMed/vignettes/getting_started_with_easyPubMed.html

### 2.) PubmedCentral Scraper
# - Rscripts to scrape data from Pubmed;
# - but may be outdated [???];
#  -- uses rvest:html(), which is deprecated;
# https://github.com/EFavDB/PubmedCentral_Scraper

### 3.) CRAN Task View: Web Technologies
# - various packages useful to connect to web pages,
#   like curl, RCurl, etc;
# https://cran.r-project.org/web/views/WebTechnologies.html


################
################

############
### URLs ###

# we will use a specific URL, which allows direct dornload using the article ID;
URL_BMC = "https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&identifier=oai:pubmedcentral.nih.gov:"


### Article: Proteasix Ontology
id = "4893253"
URL = paste0(URL_BMC, id, "&metadataPrefix=pmc")
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4893253/


################
################

### Read page

page = read_html(URL)

################

### Processing / Extraction:

# useful selectors:
### select the body-tag
html_nodes(page, "body")
### select all sub-sections
html_nodes(page, xpath="//sec")
html_nodes(page, xpath="//sec[@id='Sec1']")

###################

### Text Extraction
# text is extracted with: html_text(_nodes_);

### Section Title
# - we will use an XPATH to select the section sub-title
sec.tit = html_text(html_nodes(page, xpath="//sec/title"))

article = data.frame(title=sec.tit, txt=NA)

### Text
# - unfortunately captures the title again;
# sec.txt = html_text(html_nodes(page, xpath="//sec"))

### Exclude <title>-tag
# sec.txt = html_text(html_nodes(page, xpath="//sec/*[not(self::title)]"))
nodes = html_nodes(page, xpath="//sec")
for(id in 1:length(nodes)) {
	str = paste(html_text(html_nodes(nodes[id], xpath="./*[not(self::title)]")), collapse="\n")
	article$txt[id] = str
	print(substr(str, 1, 500)) # print first 500 characters;
}

# article:
# - title: sub-section title;
# - txt; article text;

######################

### TODO:
# - process actual text;

