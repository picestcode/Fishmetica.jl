W=readdlm("/home/igor/.julia/v0.6/Fishmetica/test/W.txt")

x=readdlm("/home/igor/.julia/v0.6/Fishmetica/test/x.txt")

y=readdlm("/home/igor/.julia/v0.6/Fishmetica/test/y.txt")

yb=1974
tmax,jmax=size(x)

cst=cohorts_t1(yb,tmax,jmax,x,y,w); 

eq_csmtx(csv2cst(cst2csv2(cst)),cst)