#test1
F=readdlm("/home/igor/.julia/v0.6/Fishmetica/test/F.txt")

M=readdlm("/home/igor/.julia/v0.6/Fishmetica/test/M.txt")

W=readdlm("/home/igor/.julia/v0.6/Fishmetica/test/W.txt")

x=readdlm("/home/igor/.julia/v0.6/Fishmetica/test/x.txt")

y=readdlm("/home/igor/.julia/v0.6/Fishmetica/test/y.txt")

Fishmetica.testgmx(x,M,F)
Fishmetica.testgmy(x,y,M,F,W,1,1)