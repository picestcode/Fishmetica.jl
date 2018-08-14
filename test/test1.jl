#test1
F=readdlm("F.txt")

M=readdlm("M.txt")

W=readdlm("W.txt")

x=readdlm("x.txt")

y=readdlm("y.txt")

Fishmetica.testgmx(x,M,F)
Fishmetica.testgmy(x,y,M,F,W,1,1)
