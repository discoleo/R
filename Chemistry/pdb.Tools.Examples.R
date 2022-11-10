

# source("pdb.Tools.R");

# also required, see in "R/Stats";
# source("Tools.Format.R");

###########################

# setwd(...)


# some pdb file;
x = read.pdb("pdb1duz.ent.gz")

tt = title.pdb(x)
# print.str(tt)
# should print automatically
tt

length.chains(x)

cat.aa(as.aa(x, "A"))


#############

### All files

meta = read.meta.pdb(FUN = function(x) paste0(filter.oligo(x, maxChains=1)$aa, collapse=" "))
names(meta)[match("FUN", names(meta))] = "PP";

# TODO:
# - save the results first;
#   (e.g. if you have a large collection of pdb files)

scroll.txt(meta, w=c(8, 56, 20, 3, 6, 15, 17))


### Plot AA Frequency
# - looks better with a few hundred PPs;
plot.aa.freq(summary.set.pp(meta$PP))

