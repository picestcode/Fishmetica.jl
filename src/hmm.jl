"""
generates Z (full)  with M (natural mortality) and F (catch mortality)
"""
function genZ(M,F)
    (tmax,jmax)=size(M)
    Z=zeros(tmax,jmax)
    Z=M+F
    return(Z)
end

"""
generates G with M (natural mortality) and F (catch mortality), Z (full mortality)
"""
function genG(M,F)
    (tmax,jmax)=size(M)
    Z=genZ(M,F)
    G=zeros(tmax,jmax)
    for t in 1:tmax
        for j in 1:jmax
            G[t,j]=F[t,j]/Z[t,j]*(1-exp(-Z[t,j]))
        end
    end
    return(G)
end

#FILTERING

"""
filtering: returns a vector of deviations (sigma^2^f) for a cohort 
sigma^2_0 - variance  for initial age
sigma^2_e -  variance  for the evidance vector (catches)
sigma^2_x -  variance  for the state vector (population)
l - the length of the cohort
"""
function sigma2f(sigma20,sigma2x,sigma2e,l)
    s2f=zeros(l)
    s2f[1]=sigma20*sigma2e/(sigma20+sigma2e)
    if l >1
        for i in 1:l-1
        s2f[i+1]=(s2f[i]+sigma2x)*sigma2e/(s2f[i]+sigma2x+sigma2e)
        end
    end
    return(s2f)
end

"""
filtering: returns a vector of means (mu^f) for a cohort 
alpha - coefficients for the state x_{t+1}=x_t +alpha_t
beta - coefficients for the observation y_{t}=x_t +beta_t
 e - the evidance vector
mu_0 - mean for the initial age
sigma^2_0 - variance (deviation is the square root of its variance) for initial age
mu_e - mean for the evidence
sigma^2_e - variance  for the evidance  (catch)
mu_x - mean for the state
sigma^2_e - variance for the state (population)
l - the length of the cohort
"""
function muf(alpha, beta, e,mu0,sigma20,mux,sigma2x,mue,sigma2e)
    l=length(e)
    muf=zeros(l)
    muf[1]=mu0
    #muf[1]=x0+beta[1]+mue
    s2f= sigma2f(sigma20,sigma2e,sigma2x,l)
    if l>1     
        for i in 1:l-1
            muf[i+1]=(muf[i]+alpha[i]+mux)+(e[i+1]-(muf[i]+alpha[i]+mux+beta[i+1]+mue))*
        (s2f[i]+sigma2x)/(s2f[i]+sigma2x+sigma2e)
        end
    end
        return(muf)
    end


"""
filtering: returns filtered cohorts represented as vecs: csvlogts1muf - means,csvlogts1s2f - sigmas^2 (deviations),
"""
function genmus2f(ts1,c1,yb,M,F,mu0,s20,mux,s2x,muy,s2y)
    tmax,jmax=size(ts1)
    
    csvlogts1= csmtx2csvs(csmtx(yb,tmax,jmax,map(x->log(x),ts1)))
    csvlogc1= csmtx2csvs(csmtx(yb,tmax,jmax,map(x->log(x),c1)))
    csvZ= csmtx2csvs(csmtx(yb,tmax,jmax,genZ(M,F)))
    csvlogG= csmtx2csvs(csmtx(yb,tmax,jmax,map(x->log(x),genG(M,F))))
    
    csvlogts1muf= deepcopy(csvlogts1)
    csvlogts1s2f= deepcopy(csvlogts1)
   
    
    for i in 1:tmax+jmax-1
        x=csvlogts1.vecs[i].vec
        e=csvlogc1.vecs[i].vec
        alpha=-csvZ.vecs[i].vec
        beta=csvlogG.vecs[i].vec
                
        vec=muf(alpha, beta, e,x[1],s20,mux,s2x,muy,s2y)
        vecs2=sigma2f(s20,s2x,s2y,length(e))
        
        csvlogts1muf.vecs[i].vec=vec
          csvlogts1s2f.vecs[i].vec=vecs2
       
         
        
       # if i == n
        #    println(vec)
       #     println(csvlogts1f.vecs[i].vec)
      #  println(x)
            
            
       # end
     
    
    end
    return(csvlogts1muf,csvlogts1s2f)
    
end


#PREDICTION





"""
prediction: returns the vector of deviation predictions (sigma^2^p_{t+1},..., sigma^2^p_{t+k+1}) for a cohort 
sigma^2_x - variance for the state vector (population)
l - the length of the cohort
"""
function sigma2p(sigma2ft,sigma2x,k)
    sigma2pk=zeros(k+1)
    sigma2pk[1]=sigma2ft+sigma2x
    for t in 1:k
        sigma2pk[t+1]=sigma2pk[t]+sigma2x
    end
    return(sigma2pk)
end
    
"""
prediction: returns the vector of mean predictions (mu^p_{t+1}, ..., mu^p_{t+k+1}) for a cohort 
alphak=(alpha^t, ..., alpha_{t+k}) - (virtual) coefficients for x_{t+k+1}=x_{t+k} +alpha_{t+k} - 
muft - mu^f_t
mu_x - mean for the state
k >= 1 - the prediction step
"""
function mup(alphatk,mux,muft,k)
    mupk=zeros(k+1)
    mupk[1]=muft+alphatk[1]+mux
            for t in 1:k
            mupk[t+1]=mupk[t]+alphatk[t+1]+mux
        end
        return(mupk)
    end

"""
makes k predictions for all cohorts (a^e+1, .., a^e+k)  with yr>=yr2. cuts extra elements with age greater than jmax. 
returns filtered cohorts represented as vecs: csvlogts1mup - means,csvlogts1s2p - sigmas^2 (deviations),
mup and s2p (means and sigmas^2 (deviations))
"""
function genmus2p(ts1,c1,yb,M,F,mu0,s20,mux,s2x,muy,s2y,k)
    tmax,jmax=size(ts1)
    
    csvlogts1= csmtx2csvs(csmtx(yb,tmax,jmax,map(x->log(x),ts1)))
    csvlogc1= csmtx2csvs(csmtx(yb,tmax,jmax,map(x->log(x),c1)))
    csvZ= csmtx2csvs(csmtx(yb,tmax,jmax,genZ(M,F)))
    csvlogG= csmtx2csvs(csmtx(yb,tmax,jmax,map(x->log(x),genG(M,F))))
    
    csvlogts1muf= deepcopy(csvlogts1)
    csvlogts1s2f= deepcopy(csvlogts1)
    csvlogts1mufp= deepcopy(csvlogts1)
     csvlogts1mufm= deepcopy(csvlogts1)
    
    csvlogts1mup= deepcopy(csvlogts1)
    csvlogts1s2p= deepcopy(csvlogts1)
    
    

    
    
    
    
    yr1=yb
    yr0=yr1-jmax+1
    yr2=yr1+tmax-jmax
    yr3=yr1+tmax-1
    
    
    
    for i in 1:tmax+jmax-1
        x=csvlogts1.vecs[i].vec
        e=csvlogc1.vecs[i].vec
        alpha=-csvZ.vecs[i].vec
        beta=csvlogG.vecs[i].vec
        yr=csvlogts1.vecs[i].yb
                
        vec=muf(alpha, beta, e,x[1],s20,mux,s2x,muy,s2y)
        vecs2=sigma2f(s20,s2x,s2y,length(e))
        
        csvlogts1muf.vecs[i].vec=vec
        csvlogts1s2f.vecs[i].vec=vecs2
        
        
        
       # if i == n
        #    println(vec)
       #     println(csvlogts1f.vecs[i].vec)
      #  println(x)
        csvlogts1mup.vecs[i].vec=vec
        
        if k > 0
        muzs=zeros(k) 
         
        s2zs=zeros(k)
         
               
          
        
        
        
        
        
            
       # end
     if yr> yr2
        begin
            muzs[1],s2zs[1]=last(vec)+last(alpha)+mux,last(vecs2)+s2x
           
            for l in 2:k
               
                      muzs[l],s2zs[l]= 
                        if length(vec)+l<=jmax
                        muzs[l-1]+last(alpha)+mux,s2zs[l-1]+s2x
                        else
                           muzs[l],s2zs[l]= muzs[jmax-length(vec)],s2zs[jmax-length(vec)]  
                        end
            end
                       #if length(vec)+k+1 <= jmax
                        #    for l in length(vec)+k+1:jmax
                         # muzs[l],s2zs[l]=muzs[length(vec)+k],s2zs[length(vec)+k]
                       # end
                       # end
            
        
                    #println(i)
                    #println(yr)
         #println(muzs)
         #  println(s2zs)         
                    
      append!(vec,muzs)
      append!(vecs2,s2zs) 
                    if length(vec)>jmax
                        vec, vecs2= vec[1:jmax], vecs2[1:jmax]
                    end
        end  
        end
            csvlogts1mup.vecs[i].vec=vec
            csvlogts1s2p.vecs[i].vec=vecs2
           
    
    end
    end
    
    
    
    
    return(csvlogts1mup,csvlogts1s2p)
    
end
    



#SMOOTHING

"""
backward: returns a vector of deviations (sigma^2^b) for a cohort 
sigma^2_0 - variance  for initial age
sigma^2_e -  variance  for the evidance vector (catches)
sigma^2_x -  variance  for the state vector (population)
l - the length of the cohort
"""
function sigma2b(sigma20,sigma2x,sigma2e,l)
    s2b=zeros(l)
    s2f=sigma2f(sigma20,sigma2x,sigma2e,l)
    s2b[l]=sigma2x+sigma2e
    if l >1
        i=l
        while i>1
        s2b[i-1]=(sigma2x*sigma2e +s2b[i]*(sigma2x+sigma2e))/(sigma2e+s2b[i])
            i=i-1
        end
    end
    return(s2b)
end
    
"""
smoothing: makes several cohorts represented as vecs: csvlogts1mus - means,csvlogts1s2s - sigmas^2 (deviations)
"""
function genmus2s(ts1,c1,yb,M,F,mu0,s20,mux,s2x,muy,s2y)
    tmax,jmax=size(ts1)
    
    csvlogts1= csmtx2csvs(csmtx(yb,tmax,jmax,map(x->log(x),ts1)))
    csvlogc1= csmtx2csvs(csmtx(yb,tmax,jmax,map(x->log(x),c1)))
    csvZ= csmtx2csvs(csmtx(yb,tmax,jmax,genZ(M,F)))
    csvlogG= csmtx2csvs(csmtx(yb,tmax,jmax,map(x->log(x),genG(M,F))))
    
    csvlogts1muf= deepcopy(csvlogts1)
    csvlogts1s2f= deepcopy(csvlogts1)
         
    csvlogts1mub= deepcopy(csvlogts1)
    csvlogts1s2b= deepcopy(csvlogts1)
    
    csvlogts1mus= deepcopy(csvlogts1)
    csvlogts1s2s= deepcopy(csvlogts1)
    
    for i in 1:tmax+jmax-1
        x=csvlogts1.vecs[i].vec
        e=csvlogc1.vecs[i].vec
        alpha=-csvZ.vecs[i].vec
        beta=csvlogG.vecs[i].vec
        
        
        
        li=length(e)
        
        
        #forward
        
        vecmuf=muf(alpha, beta, e,x[1],s20,mux,s2x,muy,s2y)
        vecs2f=sigma2f(s20,s2x,s2y,li)
        
        csvlogts1muf.vecs[i].vec=vecmuf
          csvlogts1s2f.vecs[i].vec=vecs2f
         
        #backward
        
        vecmub=zeros(li)
       vecs2b=sigma2b(s20,s2x,s2y,li)
        
       
        
        vecmub[li]=vecmuf[li]
        
       
        
       if li >1
            vecmub[li]=e[li]-alpha[li-1]-mux-beta[li]-muy
        i1=li
             
            
        while i1>2
        vecmub[i1-1]=((vecmub[i1]-alpha[i1-2]-mux)*s2y+(e[i1-1]-alpha[i1-2]-beta[i1-1]-mux)*vecs2b[i1])/
                (s2y+vecs2b[i1])
                i1=i1-1
                
        end
    end
         
        csvlogts1mub.vecs[i].vec=vecmub
          csvlogts1s2b.vecs[i].vec=vecs2b
        
     #smoothing
        
         vecmus=zeros(li)
        vecs2s= zeros(li)
        
        vecmus[li]=vecmuf[li]
        vecs2s[li]=vecs2f[li]
        
        if li >1
            
        for i2 in 1:li-1
                 
            vecmus[i2]=(vecmuf[i2]*vecs2b[i2+1]+vecmub[i2+1]*vecs2f[i2])/(vecs2b[i2+1]+vecs2f[i2])
            vecs2s[i2]=(vecs2b[i2+1]*vecs2f[i2])/(vecs2b[i2+1]+vecs2f[i2])
                #println("i2= ",i2)
                #println("vecs2s=", vecs2s[i2],"  vecs2f=", vecs2f[i2]," vecs2b=", vecs2b[i2])
        end
            end    
         csvlogts1mus.vecs[i].vec=vecmus
          csvlogts1s2s.vecs[i].vec=vecs2s
        
        
    end
    
    
    
    
    
    
    
    return(csvlogts1mus,csvlogts1s2s)
    
end 


"""
at +p%
"""
function genQuantile(lmu,lsigma2,p)
    l=length(lmu)
    lq=zeros(l) 
    for i in 1:l
        lq[i]=if lsigma2[i]>0. 
           quantile(Distributions.Normal(lmu[i],sqrt(lsigma2[i])),p)
            else 
            lmu[i]
        end
  
        end
    return(lq)
end

"""
at -p%
"""
function genQuantile1(lmu,lsigma2,p)
    l=length(lmu)
    lq=zeros(l) 
    for i in 1:l
        lq[i]=if lsigma2[i]>0. 
           2*lmu[i]-quantile(Distributions.Normal(lmu[i],sqrt(lsigma2[i])),p)
            else 
            lmu[i]
        end
  
        end
    return(lq)
end
    