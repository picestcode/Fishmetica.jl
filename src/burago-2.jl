"""
returns sigma for  j
"""
function sigma(j,as,bs,jmax)
    rho=(j-as)/bs
1/(1+exp(-rho))
end


"""
returns u-function for j
"""
function ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
        rho1=(a1m-j)/b1m
        rho2=(j-a2m)/b2m
        cm+
    if j<a1m rho1^d1m
            elseif a1m<=j<=a2m 0
        elseif a2m<=j<=jmax rho2^d2m
    end
end


"""
returns wieght for t,j (according to the ovsyannikov data)
"""
function w(t,j,jmax)
   wv=[0.012 0.052 0.11 0.179 0.253 0.323 0.483 0.677 0.908 1.116 1.297 1.541 1.758 1.971 1.997 2.358 2.475];
    return wv[j]
end

"""
wtj
"""
function gwtj(w,tmax,jmax)
  wtj= fill(1.0,tmax,jmax);
    for t in [1:tmax;] 
        for j in [1:jmax;] wtj[t,j]=w(t,j,jmax) 
            end 
    end
    return wtj
end

"""
fills H with 1's
"""
function h(t,j,jmax)
    hv=fill(1.0,jmax)';
    return hv[j]
end

"""
mtj
"""
function ghtj(h,tmax,jmax)
   htj= fill(0.0,tmax,jmax);
    for t in [1:tmax;] 
        for j in [1:jmax;] htj[t,j]=h(t,j,jmax) 
            end 
    end
    return htj
end


"""
mtj
"""

function gmtj(a1m,a2m,b1m,b2m,cm,d1m,d2m,htj,tmax,jmax)
   mtj= fill(1.0,tmax,jmax);
    for t in [1:tmax;] 
        for j in [1:jmax;] mtj[t,j]= htj[t,j]*ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax) 
            end 
    end
    return mtj
end

"""
ft=0.2
"""
function gft(tmax,jmax)
    ft= fill(0.2,tmax);
    return ft
end


"""
ftj
"""
function gftj(as,bs,ft,tmax,jmax)
    ftj= fill(1.0,tmax,jmax);
    for t in [1:tmax;] 
        for j in [1:jmax;] ftj[t,j]=ft[t]*sigma(j,as,bs,jmax) 
            end 
    end
    return ftj
end

"""
x
"""

function gx(x,t0,t1,tmax,j0,j1,jmax,M,F)
    for t in [t0:t1;]
        
            for i in [0:j1-j0-1;]
            #println(t+i+1," ",j0+i+1," ",x[t+i,j0+i])
                x[t+i+1,j0+i+1]=x[t+i,j0+i]*exp(-M[t+i,j0+i]-F[t+i,j0+i])
            
            end
        
    end
    return x
end

"""
y
"""
function gy(y,x,t0,t1,tmax,j0,j1,jmax,M,F,wtj,kx,kw)
    for t in [t0:t1;]
        for i in [0:j1-j0;]           
      y[t+i,j0+i]=(kw*wtj[t+i,j0+i])*(kx*x[t+i,j0+i])*F[t+i,j0+i]/(M[t+i,j0+i]+F[t+i,j0+i])*(1- exp(-M[t+i,j0+i]-F[t+i,j0+i]))
    
        end
    end
    return y
end

"""
x,y
"""
function gm_burago(tmax,jmax,M,F,W,kx,kw,x1j)
#0
    x=fill(0.0,tmax,jmax);
   
    for j in [1:jmax;] x[1,j]=x1j[j] end
    
    for t in 1:tmax x[t,1]=x1j[1] end
   
 #1   
        for t in [1:tmax-jmax+1;]
           for i in [0:jmax-2;]
            #println(t+i+1," ",i+2,"<-", t+i," ",i+1)
            x[t+i+1,i+1+1]=x[t+i,i+1]*exp(-M[t+i,i+1]-F[t+i,i+1])
           end
        end
 
 #2   
    for t in [tmax-jmax+2:tmax-1;]
           for i in [0:tmax-t-1;]
            #println(t+i+1," ",i+2,"<-", t+i," ",i+1)
            x[t+i+1,i+1+1]=x[t+i,i+1]*exp(-M[t+i,i+1]-F[t+i,i+1])
           end
        end
   #3
    
    for j in [2:jmax;]
           for i in [0:jmax-j-1;]
           #println(i+1+1,"! ",j+i+1,"<-", i+1," ",j+i)
            x[i+1+1,j+i+1]=x[i+1,j+i]*exp(-M[i+1,j+i]-F[i+1,j+i])
           end
        end
    
    
    
 #0   
     y=fill(0.0,tmax,jmax);
  
#1
   for t in [1:tmax-jmax+1;]
        for i in [0:jmax-1;]  
            #println(t+i," ",i+1)
            Z=M[t+i,i+1]+F[t+i,i+1]
            y[t+i,i+1]=
            (kw*W[t+i,i+1])*(kx*x[t+i,i+1])*F[t+i,i+1]/Z*(1- exp(-Z))
            
        end
    end
    
   #2 
     for t in [tmax-jmax+2:tmax;]
        for i in [0:tmax-t;]  
            #println(t+i," ",i+1)
            Z=M[t+i,i+1]+F[t+i,i+1]
            y[t+i,i+1]=
            (kw*W[t+i,i+1])*(kx*x[t+i,i+1])*F[t+i,i+1]/Z*(1- exp(-Z))
             
        end
    end 
    
    
    #3
    
    for j in [2:jmax;]
           for i in [0:jmax-j;]
    
     Z=M[i+1,j+i]+F[i+1,j+i]
            y[i+1,j+i]=
            (kw*W[i+1,j+i])*(kx*x[i+1,j+i])*F[i+1,j+i]/Z*(1- exp(-Z))
            
        end
    end 
    
    #gy(y,x,1,tmax-jmax+1,tmax,1,jmax,jmax,M,F,W,kx,kw)
   
           
    return x,y
end


"""
x,y
"""

function gm(tmax,jmax,M,F,W,kx,kw,x1j,xt1)
#0
    x=fill(0.0,tmax,jmax);
   
    for j in [1:jmax;] x[1,j]=x1j[j] end
    
    for t in 1:tmax x[t,1]=xt1[t] end
   
 #1   
        for t in [1:tmax-jmax+1;]
           for i in [0:jmax-2;]
            #println(t+i+1," ",i+2,"<-", t+i," ",i+1)
            x[t+i+1,i+1+1]=x[t+i,i+1]*exp(-M[t+i,i+1]-F[t+i,i+1])
           end
        end
 
 #2   
    for t in [tmax-jmax+2:tmax-1;]
           for i in [0:tmax-t-1;]
            #println(t+i+1," ",i+2,"<-", t+i," ",i+1)
            x[t+i+1,i+1+1]=x[t+i,i+1]*exp(-M[t+i,i+1]-F[t+i,i+1])
           end
        end
   #3
    
    for j in [2:jmax;]
           for i in [0:jmax-j-1;]
           #println(i+1+1,"! ",j+i+1,"<-", i+1," ",j+i)
            x[i+1+1,j+i+1]=x[i+1,j+i]*exp(-M[i+1,j+i]-F[i+1,j+i])
           end
        end
    
    
    
 #0   
     y=fill(0.0,tmax,jmax);
  
#1
   for t in [1:tmax-jmax+1;]
        for i in [0:jmax-1;]  
            #println(t+i," ",i+1)
            Z=M[t+i,i+1]+F[t+i,i+1]
            y[t+i,i+1]=
            (kw*W[t+i,i+1])*(kx*x[t+i,i+1])*F[t+i,i+1]/Z*(1- exp(-Z))
        end
    end
    
   #2 
     for t in [tmax-jmax+2:tmax;]
        for i in [0:tmax-t;]  
            #println(t+i," ",i+1)
            Z=M[t+i,i+1]+F[t+i,i+1]
            y[t+i,i+1]=
            (kw*W[t+i,i+1])*(kx*x[t+i,i+1])*F[t+i,i+1]/Z*(1- exp(-Z))
        end
    end 
    
    
    #3
    
    for j in [2:jmax;]
           for i in [0:jmax-j;]
    
     Z=M[i+1,j+i]+F[i+1,j+i]
            y[i+1,j+i]=
            (kw*W[i+1,j+i])*(kx*x[i+1,j+i])*F[i+1,j+i]/Z*(1- exp(-Z))
        end
    end 
    
    #gy(y,x,1,tmax-jmax+1,tmax,1,jmax,jmax,M,F,W,kx,kw)
   
           
    return x,y
end



"""
adds some noise with d(m,s) (laplace, normal) to a matrix
"""
function noise2ts_burago(ts,dn,dc)
    (tmax,jmax)=size(ts.x)
      x1=fill(0.0,tmax,jmax)
    y1=fill(0.0,tmax,jmax)
    for t in 1:tmax
        for j in 1:jmax
            x1[t,j]=ts.x[t,j]*exp(rand(dn))
            y1[t,j]=ts.y[t,j]*exp(rand(dc))
        end
    end
    return cohorts_t(ts.year_b,ts.year_e,x1,y1,ts.w)
     
end


"""
constraint
"""
function buragoconstraint1(x::Vector, jmax)
    tmax=div((length(x)-9),4)
        x[4*tmax+1]-x[4*tmax+2]
end
"""
constraint
"""
function buragoconstraint2(x::Vector, jmax)
     tmax=div((length(x)-9),4)
    x[4*tmax+2]-jmax
end
"""
constraint
"""
 function buragoconstraint3(x::Vector, jmax)
     tmax=div((length(x)-9),4)
    x[4*tmax+8]-jmax   
end
"""
constraint
"""
function buragoconstraint4(x::Vector, kappa)
     tmax=div((length(x)-9),4)
    p=true
    for i in 1:3*tmax
        p=p&(1/kappa<=x[i]<=kappa)
    end 
    if p
        -1.
    else 1.
    end
end
"""
constraint
"""
function buragoconstraint5(x::Vector,i, kappa1)
    kappa1-x[i]
    end
"""
constraint
"""
function buragoconstraint6(x::Vector,i, kappa2)
    x[i]-kappa2
    end
"""
Laplace index evaluation for cohorts of cst type/Lp, p=1
    using  closures in opt: (x,grad)->buragofunc(x,cst,kx,kw,grad) 
"""
function buragofunc(x,cst,kx,kw,grad)
 
    
    (tmax,jmax)=size(cst.x)
    
    local Hl=fill(1.0,tmax,jmax)
    local ftl=fill(1,tmax)
    
   
    
    if length(grad) > 0
        grad[1] = 0
        grad[2] = 0.5/sqrt(x[2])
    end    
    
    
   
    
    
  
for t in 1:tmax
        Hl[t,1:3]=x[3*(t-1)+1:3*t]
    end
    
    local a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml,asl,bsl
    ftl=x[3*tmax+1:4*tmax]
    a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml=x[4*tmax+1:4*tmax+7]
    asl,bsl=x[4*tmax+8:4*tmax+9]
    
    
    
   
    
    #println(cst,psl,Hl,ftl)
    return Lp(cst,a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml,Hl,asl,bsl,ftl,kx,kw,1)
   
    #println("f_$count($x)=$y")

    
end



"""
Gauss index evaluation for cohorts of cst type/Lp, p=2
     using  closures in opt: (x,grad)->buragofunc(x,cst,kx,kw,grad) 
"""
#p=2
function buragofuncp2(x,cst,kx,kw,grad)
 
    
    (tmax,jmax)=size(cst.x)
    
    local Hl=fill(1.0,tmax,jmax)
    local ftl=fill(1,tmax)
    
       
    if length(grad) > 0
        grad[1] = 0
        grad[2] = 0.5/sqrt(x[2])
    end    
    
    
   
 
    
    
  
for t in 1:tmax
        Hl[t,1:3]=x[3*(t-1)+1:3*t]
    end
    
    local a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml,asl,bsl
    ftl=x[3*tmax+1:4*tmax]
    a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml=x[4*tmax+1:4*tmax+7]
    asl,bsl=x[4*tmax+8:4*tmax+9]
    
    
    
   
    
    #println(cst,psl,Hl,ftl)
    return Lp(cst,a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml,Hl,asl,bsl,ftl,kx,kw,2)
   
    #println("f_$count($x)=$y")

   
end

"""
index (Lmine2) evaluation for cohorts of cst type/Lmine2, p=2
     using  closures in opt: (x,grad)->buragofunc(x,cst,kx,kw,grad) 
"""
    
"""
function minefuncp2(x::Vector, grad::Vector)
 
   
    (tmax,jmax)=size(cst.x)
    
    local Hl=fill(1.0,tmax,jmax)
    local ftl=fill(1,tmax)
    
    
    
    if length(grad) > 0
        grad[1] = 0
        grad[2] = 0.5/sqrt(x[2])
    end    
    
    
  
 
    
    
  
for t in 1:tmax
        Hl[t,1:3]=x[3*(t-1)+1:3*t]
    end
    
    local a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml,asl,bsl
    ftl=x[3*tmax+1:4*tmax]
    a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml=x[4*tmax+1:4*tmax+7]
    asl,bsl=x[4*tmax+8:4*tmax+9]
    
    
    
   
    
    #println(cst,psl,Hl,ftl)
    return Lmine2(cst,a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml,Hl,asl,bsl,ftl,kx,kw,2)
   
    #println("f_$count($x)=$y")

    
end




"""
index (Lmine2max) evaluation for cohorts of cst type/Lmine2max, p=2
     using  closures in opt: (x,grad)->buragofunc(x,cst,kx,kw,grad) 
"""
"""
function minefuncp2max(x,cst,kx,kw,grad)
 
    global cst
    (tmax,jmax)=size(cst.x)
    
    local Hl=fill(1.0,tmax,jmax)
    local ftl=fill(1,tmax)
    
    global kx,kw
    global count
    
    if length(grad) > 0
        grad[1] = 0
        grad[2] = 0.5/sqrt(x[2])
    end    
    
    
    count::Int += 1
 
    
    
  
for t in 1:tmax
        Hl[t,1:3]=x[3*(t-1)+1:3*t]
    end
    
    local a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml,asl,bsl
    ftl=x[3*tmax+1:4*tmax]
    a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml=x[4*tmax+1:4*tmax+7]
    asl,bsl=x[4*tmax+8:4*tmax+9]
    
    
    
   
    
    #println(cst,psl,Hl,ftl)
    return Lmine2max(cst,a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml,Hl,asl,bsl,ftl,kx,kw,2)
   
    #println("f_$count($x)=$y")

    
end



"""
index (Lp1) evaluation for cohorts of cst type/Lp12, p=2
"""
function buragofunc12(x,cst,kx,kw,grad)
 
    
    (tmax,jmax)=size(cst.x)
    
    local Hl=fill(1.0,tmax,jmax)
    local ftl=fill(1,tmax)
    
    
    
    if length(grad) > 0
        grad[1] = 0
        grad[2] = 0.5/sqrt(x[2])
    end    
    
    
    
 
    
    
  
for t in 1:tmax
        Hl[t,1:3]=x[3*(t-1)+1:3*t]
    end
    
    local a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml,asl,bsl
    ftl=x[3*tmax+1:4*tmax]
    a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml=x[4*tmax+1:4*tmax+7]
    asl,bsl=x[4*tmax+8:4*tmax+9]
    
    

    
   
    
    #println(cst,psl,Hl,ftl)
    return Lp12(cst,a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml,Hl,asl,bsl,ftl,kx,kw,1)
   
    #println("f_$count($x)=$y")

    
end



"""
list of optimization algorithms
"""
ALGS=(:LN_NELDERMEAD,:LN_SBPLX,:LN_COBYLA,:LN_NEWUOA,:LN_PRAXIS,:LN_BOBYQA,:LN_AUGLAG,:LN_AUGLAG_EQ)

#


"""
optimizes func with 2 algs (a -main and a1 - local), (burago)func (Lp, Lp12, etc)  and 4 buragoconstraints 
"""
function optimizeLf(x,func,tmax,jmax,a,a1,lb,lb1,xtol,xtol1,bc4,ctol,is1,is2,is3,pop,pop1,maxe)
   
   
global count = 0# keeps track of # of function evaluations
opt = Opt(ALGS[a], 4*tmax+9)
l_opt= Opt(ALGS[a1], 4*tmax+9)
    
    
    println(ALGS[a])
    println(ALGS[a1])
    
    
    
    lower_bounds!(opt, fill(lb,4*tmax+9))
    lower_bounds!(l_opt, fill(lb1,4*tmax+9))
    xtol_rel!(opt,xtol)
    xtol_rel!(l_opt,xtol1)

min_objective!(opt, func)
    inequality_constraint!(opt, (x,g) -> buragoconstraint1(x,jmax), ctol)
    inequality_constraint!(opt, (x,g) -> buragoconstraint2(x,jmax), ctol)
    inequality_constraint!(opt, (x,g) -> buragoconstraint3(x,jmax),ctol )
    inequality_constraint!(opt, (x,g) -> buragoconstraint4(x,bc4), ctol)
#for i in 1:3*tmax
#    inequality_constraint!(opt, (x,g) -> buragoconstraint5(x,i,0.1), 1e-8)
#     inequality_constraint!(opt, (x,g) -> buragoconstraint6(x,i,5.), 1e-8)
#end

local_optimizer!(opt,l_opt)


is=fill(0.4,4*tmax+9)
    is[1:3*tmax]=fill(is1,3*tmax)#h
    is[3*tmax+1:4*tmax]=fill(is2,tmax)#f
    is[4*tmax+1:4*tmax+9]=fill(is3,9)#pm&pf


initial_step!(opt,is)
population!(opt,pop)
population!(l_opt,pop1)
NLopt.srand_time()
    maxeval!(opt, maxe)

(minf,minx,ret) = NLopt.optimize(opt,x)

end

"""
optimizes Lp with 2 algs (a -main and a1 - local), buragofunc   and 4 buragoconstraints
"""
#
function optimizeL(x,cst,kx,kw,a,a1,lb,lb1,xtol,xtol1,bc4,ctol,is1,is2,is3,pop,pop1,maxe)
   (tmax,jmax)=size(cst)
   
global count = 0# keep track of # of function evaluations
opt = Opt(ALGS[a], 4*tmax+9)
l_opt= Opt(ALGS[a1], 4*tmax+9)
    
    
   # println(ALGS[a])
    # println(ALGS[a1])
    
    
    
    lower_bounds!(opt, fill(lb,4*tmax+9))
    lower_bounds!(l_opt, fill(lb1,4*tmax+9))
    xtol_rel!(opt,xtol)
    xtol_rel!(l_opt,xtol1)

min_objective!(opt, (x,grad)->buragofunc(x,cst,kx,kw,grad))
    inequality_constraint!(opt, (x,g) -> buragoconstraint1(x,jmax), ctol)
    inequality_constraint!(opt, (x,g) -> buragoconstraint2(x,jmax), ctol)
    inequality_constraint!(opt, (x,g) -> buragoconstraint3(x,jmax),ctol )
    inequality_constraint!(opt, (x,g) -> buragoconstraint4(x,bc4), ctol)
#for i in 1:3*tmax
#    inequality_constraint!(opt, (x,g) -> buragoconstraint5(x,i,0.1), 1e-8)
#     inequality_constraint!(opt, (x,g) -> buragoconstraint6(x,i,5.), 1e-8)
#end

local_optimizer!(opt,l_opt)


is=fill(0.4,4*tmax+9)
    is[1:3*tmax]=fill(is1,3*tmax)#h
    is[3*tmax+1:4*tmax]=fill(is2,tmax)#f
    is[4*tmax+1:4*tmax+9]=fill(is3,9)#pm&pf


initial_step!(opt,is)
population!(opt,pop)
population!(l_opt,pop1)
NLopt.srand_time()
    maxeval!(opt, maxe)

(minf,minx,ret) = NLopt.optimize(opt,x)

end






"""
optimizes Lp with 1 algs (a -main), buragofunc   and 4 buragoconstraints, no lower_bounds, no x_tol, no local_optimizer?
no population
"""
function optimizeLs(x,cst,kx,kw,a,a1,lb,lb1,xtol,xtol1,bc4,ctol,is1,is2,is3,pop,pop1,maxe)
   
   global count = 0
opt = Opt(ALGS[a], 4*tmax+9)
#l_opt= Opt(ALGS[a1], 4*tmax+9)
    lower_bounds!(opt, fill(lb,4*tmax+9))
    #lower_bounds!(l_opt, fill(lb1,4*tmax+9))
    xtol_rel!(opt,xtol)
    #xtol_rel!(l_opt,xtol1)

min_objective!(opt,(x,grad)->buragofunc(x,cst,kx,kw,grad))
    inequality_constraint!(opt, (x,g) -> buragoconstraint1(x,jmax), ctol)
    inequality_constraint!(opt, (x,g) -> buragoconstraint2(x,jmax), ctol)
    inequality_constraint!(opt, (x,g) -> buragoconstraint3(x,jmax),ctol )
    inequality_constraint!(opt, (x,g) -> buragoconstraint4(x,bc4), ctol)
#for i in 1:3*tmax
#    inequality_constraint!(opt, (x,g) -> buragoconstraint5(x,i,0.1), 1e-8)
#     inequality_constraint!(opt, (x,g) -> buragoconstraint6(x,i,5.), 1e-8)
#end

#local_optimizer!(opt,l_opt)


is=fill(0.4,4*tmax+9)
    is[1:3*tmax]=fill(is1,3*tmax)#h
    is[3*tmax+1:4*tmax]=fill(is2,tmax)#f
    is[4*tmax+1:4*tmax+9]=fill(is3,9)#pm&pf


initial_step!(opt,is)
population!(opt,pop)
#population!(l_opt,pop1)
NLopt.srand_time()
    maxeval!(opt, maxe)

(minf,minx,ret) = NLopt.optimize(opt,x)

end

"""
optimizes Lp with 2 algs (a -main and a1 - local), buragofunc_v   and 4 buragoconstraints
"""
function optimizeLp_v2(x,cst,kx,kw,a,a1,lb,lb1,xtol,xtol1,bc4,ctol,is1,is2,is3,pop,pop1,maxe)
   
   
global count = 0# keep track of # of function evaluations
opt = Opt(ALGS[a], 4*tmax+9)
l_opt= Opt(ALGS[a1], 4*tmax+9)
    
    
    println(ALGS[a])
    println(ALGS[a1])
    
    
    
    lower_bounds!(opt, fill(lb,4*tmax+9))
    lower_bounds!(l_opt, fill(lb1,4*tmax+9))
    xtol_rel!(opt,xtol)
    xtol_rel!(l_opt,xtol1)

min_objective!(opt, (x,grad)->buragofunc_v(x,cst,kx,kw,grad))
    inequality_constraint!(opt, (x,g) -> buragoconstraint1(x,jmax), ctol)
    inequality_constraint!(opt, (x,g) -> buragoconstraint2(x,jmax), ctol)
    inequality_constraint!(opt, (x,g) -> buragoconstraint3(x,jmax),ctol )
    inequality_constraint!(opt, (x,g) -> buragoconstraint4(x,bc4), ctol)
#for i in 1:3*tmax
#    inequality_constraint!(opt, (x,g) -> buragoconstraint5(x,i,0.1), 1e-8)
#     inequality_constraint!(opt, (x,g) -> buragoconstraint6(x,i,5.), 1e-8)
#end

local_optimizer!(opt,l_opt)


is=fill(0.4,4*tmax+9)
    is[1:3*tmax]=fill(is1,3*tmax)#h
    is[3*tmax+1:4*tmax]=fill(is2,tmax)#f
    is[4*tmax+1:4*tmax+9]=fill(is3,9)#pm&pf


initial_step!(opt,is)
population!(opt,pop)
population!(l_opt,pop1)
NLopt.srand_time()
    maxeval!(opt, maxe)

(minf,minx,ret) = NLopt.optimize(opt,x)

end




"""
evaluates vector of unknown pars with use of given 
 ALGS=(:LN_NELDERMEAD, :LN_SBPLX, :LN_COBYLA, :LN_NEWUOA, :LN_PRAXIS, :LN_BOBYQA, :LN_AUGLAG, :LN_AUGLAG_EQ)
from NLopt and LN_AUGLAG (7),  
xx - initial approximation for pars,
k - number of iterations,
algs - list of algs
"""
function opt_v(xx,k,cst,kx,kw,algs)
@eval @everywhere  xxx=$xx
flag=true
ni=0
result=0.0
while flag
    ni+=1 
    println(ni)
    pm= pmap(i->optimizeLp_v2(xxx,cst,kx,kw,7,i,0.0,0.0,1e-20,1e-20,2.,1e-20,1.,1.,0.1,5,5,50000),algs)
   
    pm1=filter(x->typeof(x)==Tuple{Float64,Array{Float64,1},Symbol},pm)
   minfs=map(i->pm1[i][1],1:length(pm1))
   
    println(minfs)
    result=minfs
    
    if minimum(minfs)==maximum(minfs) || ni>k
         println( minimum(minfs)," ",maximum(minfs));
            flag=false 
        else
        println( minimum(minfs)," ",maximum(minfs));
        mm=minimum(minfs);
        flag=true;
            @eval @everywhere  xxx=$pm[findfirst(x->$mm==x,$minfs)][2]
    end
end
    
    return(xxx)
end


