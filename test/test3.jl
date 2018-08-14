W=readdlm("W.txt")

x=readdlm("x.txt")

y=readdlm("y.txt")

yb=1974
tmax,jmax=size(x)

cst=cohorts_t1(yb,tmax,jmax,x,y,w); 

eq_csmtx(csv2cst(cst2csv2(cst)),cst)
