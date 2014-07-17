#! /usr/bin/Rscript
setwd("~/Dropbox/Min-Suraj/Research/final implementation/code_for_calcium_stuff")
require(QUIC)

cov_md = read.csv("~/Dropbox/Min-Suraj/Research/final implementation/code_for_calcium_stuff/CXX.csv", header = F)
cov_md = as.matrix(cov_md)
N = nrow(cov_md)

rho.vec = read.csv("~/Dropbox/Min-Suraj/Research/final implementation/code_for_calcium_stuff/rho_vec.csv", header = F)
rho.vec = as.numeric(rho.vec)

inv_cov_md_est_GL = QUIC(cov_md, rho.vec)$X
inv_cov_md_est_GL=matrix(inv_cov_md_est_GL, ncol = nrow(cov_md))
write.table(inv_cov_md_est_GL, file = "TEST_inv_cov_md_est_GL_mle.csv", row.names = F,  col.names = F, sep=",")

