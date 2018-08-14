module Fishmetica
 using Base.Test
 #using Distributions
#using NLopt

global ALGS


#global count


#global cst


#global jtrustsmax

global kx
global kw



global lambda




export 
gm,

gft,
sigma,
ushape,

ALGS,

set_lambda,
set_count,
set_cst,
set_jtrustsmax,

noise2ts_burago,

buragoconstraint1,
buragoconstraint2,
buragoconstraint2,
buragoconstraint3,
buragoconstraint4,
buragoconstraint5,
buragoconstraint6,

optimizeL,
optimizeLf,
optimizeL12,
optimizeLp_v2,
opt_v,
optimizeLp_v,



gm1,
ww_osp,
gg,



gmtj,
gftj,

#cohort1,
#cohort,
#cohorts,
#cohorts1,
#cohorts_m,
cohorts_t,
#cohort_mp,
#cmp2cm,
cohort_v,
cohorts_v,
cohorts_v1,
cohorts_t1,
cohorts_v2,
csmtx,
cvec,
csvecs, 

cst2csv2,
csv2cst,
csmtx2csvs,
csvs2csmtx,
eq_csmtx,



gvecs,
gvecs1,




buragofunc,
buragofuncp2,
buragofunc12,
optimizeLp_vjt_2,
buragofunc_vjt_2,
buragofunc_v,

Lpv,
Lp,
Lp1,
Lp2,
Lp12,

minefuncp2,
Lmine2,

minefuncp2max,
Lmine2max,

genZ,
genG,

sigma2f,
muf,
genmus2f,

genmus2p,


genmus2s,
genQuantile,
genQuantile1


#include("types2.jl")
include("burago-2.jl")
include("fishmetica3.jl")
include("tests.jl")
include("hmm.jl")

"""
const coefficients at main equations, default value
"""
kx=1.0


"""
const coefficients at main equations, default
"""
kw=1.0


"""
global weight for criterium, default value 
"""
lambda=0.5

"""
sets global Fishmetica.lambda to l
"""
function set_lambda(l)

global lambda
    lambda=l
    end

"""
sets global Fishmetica.count to c
"""
function set_count(c)

global count
count=c
end

"""
sets global Fishmetica.jtrustsmax to c
"""
function set_jtrustsmax(c)
"""
maximum trusted age
"""
global jtrustsmax
jtrustsmax=c
end


"""
sets global Fishmetica.cst to c
"""
function set_cst(c) 
"""
global cohorts of cst type
"""
    global cst
cst=c
end


"""
returns (minf,minx,ret)
"""
 function optimizeLp_vjt_2(x,tmax,jmax,jtrustsmax,a,a1,lb,lb1,xtol,xtol1,bc4,ctol,is1,is2,is3,pop,pop1,maxe)
   
 """
global count for optimization steps, to keep track of no. of function evaluations
"""  
global count = 0

    
    opt = Opt(ALGS[a], (jtrustsmax+1)*tmax+9)
l_opt= Opt(ALGS[a1], (jtrustsmax+1)*tmax+9)
    
    
    println(ALGS[a])
    println(ALGS[a1])
    
    
    
    lower_bounds!(opt, fill(lb,(jtrustsmax+1)*tmax+9))
    lower_bounds!(l_opt, fill(lb1,(jtrustsmax+1)*tmax+9))
    xtol_rel!(opt,xtol)
    xtol_rel!(l_opt,xtol1)

min_objective!(opt, buragofunc_vjt_2)
    inequality_constraint!(opt, (x,g) -> buragoconstraint1(x,jmax), ctol)
    inequality_constraint!(opt, (x,g) -> buragoconstraint2(x,jmax), ctol)
    inequality_constraint!(opt, (x,g) -> buragoconstraint3(x,jmax),ctol )
    inequality_constraint!(opt, (x,g) -> buragoconstraint4(x,bc4), ctol)
#for i in 1:3*tmax
#    inequality_constraint!(opt, (x,g) -> buragoconstraint5(x,i,0.1), 1e-8)
#     inequality_constraint!(opt, (x,g) -> buragoconstraint6(x,i,5.), 1e-8)
#end

local_optimizer!(opt,l_opt)


is=fill(0.4,(jtrustsmax+1)*tmax+9)
    is[1:jtrustsmax*tmax]=fill(is1,jtrustsmax*tmax)#h
    is[jtrustsmax*tmax+1:(jtrustsmax+1)*tmax]=fill(is2,tmax)#f
    is[(jtrustsmax+1)*tmax+1:(jtrustsmax+1)*tmax+9]=fill(is3,9)#pm&pf


initial_step!(opt,is)
population!(opt,pop)
population!(l_opt,pop1)
NLopt.srand_time()
    maxeval!(opt, maxe)

(minf,minx,ret) = NLopt.optimize(opt,x)

end

"""
returns Lpv
""" 
function buragofunc_vjt_2(x::Vector, grad::Vector)
 
    global cst
    (tmax,jmax)=size(cst.x)
    
    local Hl=fill(1.0,tmax,jmax)
    local ftl=fill(1.,tmax)
    
    global kx,kw
    global count
    global jtrustsmax
    
    if length(grad) > 0
        grad[1] = 0
        grad[2] = 0.5/sqrt(x[2])
    end    
    
    
    count::Int += 1
 
    
    
  
for t in 1:tmax
        Hl[t,1:jtrustsmax]=x[jtrustsmax*(t-1)+1:jtrustsmax*t]
    end
    
    local a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml,asl,bsl
    ftl=x[jtrustsmax*tmax+1:(jtrustsmax+1)*tmax]
    a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml=x[(jtrustsmax+1)*tmax+1:(jtrustsmax+1)*tmax+7]
    asl,bsl=x[(jtrustsmax+1)*tmax+8:(jtrustsmax+1)*tmax+9]
    
    
    
   
    
    #println(cst,psl,Hl,ftl)
    return Lpv(cst,a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml,Hl,asl,bsl,ftl,kx,kw,2)
   
    #println("f_$count($x)=$y")

    
end




end
