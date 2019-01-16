library(fda)
library(R.matlab)

setwd('/home/minjay/SuperDARN')
data = readMat("all_theta.mat")
#x = data$theta[1, ]*4
x = data$theta
# can get breaks from knt2brk function in matlab
breaks = c(0,1/2,1)*pi

bS = bsplineS(x, breaks, nderiv=1)
writeMat('deriv_all_theta.mat', bS=bS)