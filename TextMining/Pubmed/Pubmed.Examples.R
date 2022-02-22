########################
###
### Leonard Mada
### [the one and only]
###
### Pubmed: Examples
###
### draft v.0.1b


##################

source("Pubmed.R")


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


### Test 2:
doc = search.entrez("serotonin syndrome")
r = parse.entrez(doc)
print(r)


### Test 3:
doc = search.entrez(c("serotonin syndrome", "MAO"))
r = parse.entrez(doc)
print(r)

