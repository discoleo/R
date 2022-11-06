

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

meta = read.meta.pdb()

scroll.txt(meta, w=c(8, 56, 20))

