#test2
yb=1974
tmax,jmax=30,16
w=zeros(tmax,jmax)
x=zeros(tmax,jmax)
y=zeros(tmax,jmax)
w=readdlm("w.txt")
y=readdlm("y.txt")
x=readdlm("x.txt")
#println(size(x))
#println(tmax," ",jmax)

#println(csmtx(yb,tmax,jmax,x))
#println(csmtx2csvs(csmtx(yb,tmax,jmax,x)))
x==csvs2csmtx(csmtx2csvs(csmtx(yb,tmax,jmax,x)))
