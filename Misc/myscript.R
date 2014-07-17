#! /usr/bin/Rscript
setwd("~/Dropbox/Min-Suraj/Research/final implementation/code_for_calcium_stuff")
require(QUIC)

cov_md = read.csv("~/Dropbox/Min-Suraj/Research/final implementation/code_for_calcium_stuff/CXX.csv", header = F)
cov_md = as.matrix(cov_md)
#rho.vec = read.csv("~/Dropbox/Min-Suraj/Research/final implementation/code_for_calcium_stuff/rho_vec.csv", header = F)
#rho.vec = as.numeric(rho.vec)

sparsity = read.csv("~/Dropbox/Min-Suraj/Research/final implementation/code_for_calcium_stuff/sparsity.csv", header = F)
sparsity = as.numeric(sparsity)

N = nrow(cov_md)/2;


#### Binary search
ub = 0.5;
lb = 0;
for (i in 1:8) {
  mid = 0.5*(ub+lb)
#  if (mid<0.04) {break}
  inv_cov_md_est_GL = QUIC(cov_md, mid)$X
  inv_cov_md_est_GL = matrix(inv_cov_md_est_GL, ncol = nrow(cov_md))
  temp = inv(inv(inv_cov_md_est_GL)[1:N,1:N])
  spa = sum(temp != 0)/(N^2)
  
  if (spa > sparsity) { 
      lb = mid }
    else if (spa<sparsity) { ub = mid}
    else {break}

}
  opt_rho = ub
write.table(opt_rho,file="optimla_rho_chosen.csv")
#fit_gl_mle = matrix(inv_cov_md_est_GL, ncol = nrow(cov_md));
write.table(inv_cov_md_est_GL, file = "TEST_inv_cov_md_est_GL_mle.csv", row.names = F,  col.names = F, sep=",")
