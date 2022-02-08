########################
###
### Leonard Mada
### [the one and only]
###
### Pubmed: Examples
###
### draft v.0.1a


##################

source("Pubmed.R")


##################

### Test:

doc = search.entrez("micrographia[Title/Abstract]")
r = parse.entrez(doc)
print(r)

nStart = 40;
doc2 = search.entrez.fetch(r, nStart, type="id")
ids  = parse.entrez.ids(doc2, nStart=nStart)
count.entrez.ids(doc2)
length(ids);
print(ids)

