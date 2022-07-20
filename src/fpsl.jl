"""
uses risduals&kde!
returns xai(len,N) - support & wai(len,N) weights for filtering 
xai support for one-step prediction
wai weights
N - number of samples
len - cohort length
kdex, kdey - distributions of residuals (kdes)
logts1,logc1 - logs of population (state, use it just as initial conditions) & catch (observation)
alpha(len) - state (population, abundance) model coefficients
beta(len) - observation (catch) model coefficients
1 filtering with the simple Gordon, Salmond, Smith, 1993 method bootstrap !!!
"""


function bootstrapfiltering(N,len,rx,ry,logts1,logc1,alpha,beta)
    uai=zeros(len,N)
    vai=zeros(len,N)
    puai=zeros(len,N)
    pvai=zeros(len,N)
    wai=zeros(len,N)
    xai=zeros(len,N)
    
    
    
 # kde for rx and ry 
  kdex=KernelDensity.kde(vec(rx));
  kdey=KernelDensity.kde(vec(ry));  
   
 #  init 
   
     
    xai[1,:]=fill(logts1[1],N)
     wai[1,:]=fill(1/N,N)
    
    
    for a in 1:len-1
     
    uai[a,:]=sample(vec(rx),N; replace=true, ordered=true) 
      
    end
    
    
    for a in 1:len-1
        for i in 1:N
     xai[a+1,i]=xai[a,i]+alpha[a]+uai[a,i]    
       
    end
    end
    
    for a in 1:len
     for i in 1:N  
    vai[a,i]=logc1[a]-xai[a,i]-beta[a]
    pvai[a,i]=pdf(kdey, vai[a,i])
        
    end
    end
    #wai=rand(dx,len,N) #samples
    #pwai=pdf.(dx,wai) #probability

   
   
    
     for a in 2:len
  
        for i in 1:N

        if sum(pvai[a,:])==0. wai[a,i]=1/N else
           
          wai[a,i]=abs(pvai[a,i])/sum(pvai[a,:])#abs to avoid "zero" negative values
          end
                 end
        
    end
   
     
   
    
    
   
     
        
    
    return(xai,wai)
end

"""
returns xai, wai (f)
c- collection of ages subseteq 1:m=tmax+jmax-1
N - number of particles
 rlxL,rlyL - residuals 
"""


function xwais(c,N,rlxL,rlyL,csvts1, csvc1,csvZL,csvlogGL);

    
 
    lc=length(c)
    xais=Array{Float64}(undef,lc);
   wais=Array{Float64}(undef,lc); 
    wais=Array{Float64}(undef,lc); 
    #mads=Array{Float64}(undef,m);
    xais=[] 
    wais=[]
    wwais=[]
    
    #means=[]
for n in c
        #println(n)
       
        
    l=length(csvts1.vecs[n].vec); 
        
        

    logts1n=map(x->log(x),csvts1.vecs[n].vec)
    logc1n=map(x->log(x),csvc1.vecs[n].vec)
    
    alphan=-csvZL.vecs[n].vec
    betan=csvlogGL.vecs[n].vec
    
    zai=Array{Float64}(undef,l,N);
    wai=Array{Float64}(undef,l,N);
    wwai=Array{Float64}(undef,l,N);
        #println("filtering")
        xai,wai=bootstrapfiltering(N,l,rlxL,rlyL,logts1n,logc1n,alphan,betan)
      #  println("smoothing")
       # wwai=fixedintervalsmoothing(N,xai,wai,rlxL,alphan)
        push!(xais,xai)
         push!(wais,wai)
        # push!(wwais,wwai)
        #push!(mads,median(abs.(map(a->median(xai[a,:],Weights(wwai[a,:])),1:l)-csvx1L.vecs[n].vec)))
    
end
    return(xais,wais)
    #return(xais,wais,wwais)
end
        


"""
returns xai, wai (f),wwai (s)
c- collection of ages subseteq 1:m=tmax+jmax-1
N - number of particles
 rlxL,rlyL - residuals 
"""


function xwwais(c,N,rlxL,rlyL,csvts1, csvc1,csvZL,csvlogGL)
    
 
    lc=length(c)
    xais=Array{Float64}(undef,lc);
   wais=Array{Float64}(undef,lc); 
    wais=Array{Float64}(undef,lc); 
    #mads=Array{Float64}(undef,m);
    xais=[] 
    wais=[]
    wwais=[]
    
    #means=[]
for n in c
        #println(n)
       
        
    l=length(csvts1.vecs[n].vec); 
        
        

    logts1n=map(x->log(x),csvts1.vecs[n].vec)
    logc1n=map(x->log(x),csvc1.vecs[n].vec)
    
    alphan=-csvZL.vecs[n].vec
    betan=csvlogGL.vecs[n].vec
    
    zai=Array{Float64}(undef,l,N);
    wai=Array{Float64}(undef,l,N);
    wwai=Array{Float64}(undef,l,N);
        #println("filtering")
        xai,wai=bootstrapfiltering(N,l,rlxL,rlyL,logts1n,logc1n,alphan,betan)
        #println("smoothing")
       wwai=fixedintervalsmoothing(N,xai,wai,rlxL,alphan)
        push!(xais,xai)
         push!(wais,wai)
       push!(wwais,wwai)
        #push!(mads,median(abs.(map(a->median(xai[a,:],Weights(wwai[a,:])),1:l)-csvx1L.vecs[n].vec)))
    
end
    #return(xais,wais)
    return(xais,wais,wwais)
end
        

"""

  Kitagawa1987, Doucet2000, partical prediction
xai,wai - samples and wieghts for filtering
szkii=Array{Float64}(kkk+1,N) sorted samples
   spzkii=Array{Float64}(kkk+1,N) sorted weight
  pzkii=Array{Float64}(kkk+1,N)  weights
  zkii=Array{Float64}(kkk+1,N) samples
prediction for kkk>1 steps, kkk=1 for a=a^e, kkk=2 for a=a^e+1, etc.
xki sample for prediction of length kkk, 1, 2,...,kkk
wki - mass
wki[1,:] - filtering mass for a=ae
wki[k,:] - prediction mass for a=ae+k-1
alpha=alpha15[len] used in state model for all k
pdf kde(vec(rx))
"""


function kitagawaprediction(kkk,N,rx,alpha,xai,wai)
    xki=Array{Float64}(undef,kkk+1,N)
    wki=Array{Float64}(undef,kkk+1,N)
    kdex=kde(vec(rx));
    len=length(wai[:,1]);
  
    xki[1,:]=xai[len,:]
    wki[1,:]=wai[len,:]
    
    for k in 2:kkk+1
        for i in 1:N
        xki[k,i]=xki[k-1,i]+alpha+sample(vec(rx))
        wki[k,i]=pdf(kdex,xki[k,i]-xki[k-1,i]-alpha)*wki[k-1,i]
        end
    end
   s=0.
    
     for k in 2:kkk+1
        s=0.
    for j in 1:N
    s=s+sum(pdf(kdex,xki[k,j]-xki[k-1,j]-alpha)*wki[k-1,j])
        end
      wki[k,:]=   wki[k,:]/s    
    end   
    
    return(xki,wki)
end

"""


returns  wwai(len,N) - updated weights for xai (smoothing)
N - number of samples
xai - filtering samples
wai - filtering weights
rx - abundance  risiduals
alpha(len) - abundance model coefficients
no printing
"""


function fixedintervalsmoothing(N,xai,wai,rx,alpha)
    len=length(wai[:,1]);
    kdex=kde(vec(rx))
    puaij=Array{Float64}(undef,len,N,N);
    spuaj=Array{Float64}(undef,len,N);
    wwai=Array{Float64}(undef,len,N);
    
    
    #puaij 
    for a in 1:len-1
        println(a)
        for i in 1:N
                    for j in 1:N
               puaij[a,i,j]=wai[a,i]*pdf(kdex,xai[a+1,j]-xai[a,i]-alpha[a])
        end
            end
    end
    #println(all(isfinite, puaij) )  
      #println(length(filter(x->isfinite(x), puaij) ))  
        #println(length(filter(x->x>0,  puaij) ))
      #  spuaj
        
      for a in 1:len-1
        println(a)
            for j in 1:N
            spuaj[a,j]=  0.
                for k in 1:N
            spuaj[a,j]=spuaj[a,j]+puaij[a,k,j]
        end
                        end
        end
     #println(all(isfinite, spuaj) )  
       #println(length(filter(x->isfinite(x), spuaj) ))  
        #println(length(filter(x->x>0, spuaj) )) 
      #  wwai  
        
   wwai[len,:]=wai[len,:] 
    
  
    
    a=len-1
    
while a>=1
    #println(a)
        for i in 1:N
         wwai[a,i]=0.
        for j in 1:N
                 dv=if spuaj[a,j]>0  wwai[a+1,j]*puaij[a,i,j]/spuaj[a,j] else 0. 
end
             wwai[a,i]= wwai[a,i]+ dv
        end
           end
a-=1
end
   #println(all(isfinite, wwai) ) 
   
    
    
    return(wwai)
end


"""

2 Kitagawa1987, Doucet2000
uses risduals&kde!
returns liklihood p
N - number of samples
len - cohort length
rx - x residuals
xy - y residuals
logts1 - data on abundance
logc1 - data on catch
alpha(len) - state (population, abundance) model coefficients
beta(len) - jbservation (catch) model coefficients

"""

function kitagawaliklihood(N,len,rx,ry,logts1,logc1,alpha,beta)
    uai=zeros(len,N)
    vai=zeros(len,N)
    puai=zeros(len,N)
    pvai=zeros(len,N)
    wtai=zeros(len,N)
    xtai=zeros(len,N)
    pyy=zeros(len)
    
    
 # kde for rx and ry 
  kdex=kde(vec(rx));
  kdey=kde(vec(ry));  
   
 #  init 
   
     
    xtai[1,:]=fill(logts1[1],N)
     wtai[1,:]=fill(1/N,N)
    
    
    for a in 1:len
     
    uai[a,:]=sample(vec(rx),N; replace=true, ordered=true) 
      
    end
    
    
    for a in 1:len-1
        for i in 1:N
     xtai[a+1,i]=xtai[a,i]+alpha[a]+uai[a,i]    
       
    end
    end
    
    
    
    for a in 2:len
     for i in 1:N  
    
    wtai[a,i]=pdf(kdex, uai[a,i])
        
    end
    end
     #println(length(filter(x->x<1, wtai) )) 
    #wai=rand(dx,len,N) #samples
    #pwai=pdf.(dx,wai) #probability

   
   for a in 1:len
     for i in 1:N  
    
   vai[a,i]=logc1[a]-xtai[a,i]-beta[a]
            
        pvai[a,i]=pdf(kdey, vai[a,i])
    end
    end
    #println(length(filter(x->x<1, pvai) )) 
    
   
     pyy[1]=pdf(kdey, logc1[1]-logts1[1]-beta[1])
   
     for a in 2:len
  
       pyy[a]=0.
        for i in 1:N
             pyy[a]= pyy[a]+wtai[a-1,i]*pvai[a,i]
          
        end
       #println(length(filter(x->x<1, pyy) ))  
    end
    
   p=pyy[1]
     for a in 2:len
        p=p*pyy[a]
    end
    
    return(p)
end
