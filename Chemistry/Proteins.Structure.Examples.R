


library(bio3d)


##############


x = read.pdb(file.choose())


feat = features.pdb(x)

fi = PDB_prepare(x)

