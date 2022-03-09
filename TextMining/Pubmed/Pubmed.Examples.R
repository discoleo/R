########################
###
### Leonard Mada
### [the one and only]
###
### Pubmed: Examples
###
### draft v.0.1d


##################

source("Pubmed.R")

# Pretty formatting
# - src file: in R/Stat/Tools.CRAN.R;
# - formatting functions may get moved to a different file;
source("Tools.CRAN.R")


##################

### Tests:

### Test 1:
doc = search.entrez("micrographia")
r = parse.entrez(doc)
print(r)

nStart = 40;
doc2 = search.entrez.fetch(r, nStart, type="id")
ids  = parse.entrez.ids(doc2, nStart=nStart)
count.entrez.ids(doc2)
length(ids);
print(ids)


###########
### Test 2:
doc = search.entrez("serotonin syndrome")
r = parse.entrez(doc)
print(r)

# Retrieve Abstracts
nStart = 0; nMax = 30;
doc2 = search.entrez.fetch(r, nStart, max=nMax)
writeLines(doc2, con="_Pubmed_Test.xml")

# Titles
titles = extractTitles(doc2)
scroll.txt(titles)
scroll.txt(titles, start=15, len=20)

# Authors
countAuthors(doc2)
scroll.txt(extractAuthors(doc2, filter=FALSE))

# Journals
scroll.txt(extractJournal(doc2, type="Full"))

# Abstracts
x = extractAbstract(doc2)
table(nchar(x$Abstract) == 0)
# TODO: Abstract Sections are not yet merged;
scroll.txt(x)


###########
### Test 3:
doc = search.entrez(c("serotonin syndrome", "MAO"))
r = parse.entrez(doc)
print(r)


###########
### Test 4:
doc = search.entrez(
	queryOr("Jasplakinolide", "Pipestelide", "Calyxamide",
		"Microsclerodermin", "Chondramides", "Argyrin") )
r = parse.entrez(doc)
print(r)

# Retrieve Abstracts
nStart = 0; nMax = 120;
doc2 = search.entrez.fetch(r, nStart, max=nMax)
writeLines(doc2, con="_Pubmed_Test_Cyclo.xml")

# Titles
titles = extractTitles(doc2)
scroll.txt(titles, len=20)
scroll.txt(titles, start=20, len=20)

