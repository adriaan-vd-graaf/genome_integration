# Title     : TODO
# Objective : TODO
# Created by: adriaan
# Created on: 6/24/20

library('SNPknock')
args = commandArgs(trailingOnly=TRUE)

r_file = args[1]
alpha_file = args[2]
theta_file = args[3]
char_file = args[4]
geno_file = args[5]
out_file = args[6]

genotypes = read.csv(geno_file, sep='\t')
genotypes = as.matrix(genotypes)
mode(genotypes) <- "integer"

print(head(genotypes))

hmm = loadHMM(r_file, alpha_file, theta_file, char_file)
knockoff_geno = knockoffGenotypes(genotypes, hmm$r, hmm$alpha, hmm$theta)
write.table(knockoff_geno, out_file, sep="\t", row.names=FALSE, col.names=FALSE)

