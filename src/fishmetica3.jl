"""
cohorts as 3 matrices: year_b - initial year, year_e - last year,  x - abundance, y - catch, w - weight (obsolete)
"""
mutable struct  cohorts_t  
     year_b
    year_e
    x
     y
    w
              
end

"""
detail cohort representation 
tmax - number of seasons
jmax - number of ages
yb - physical year assciated with cohort (when age=1)
jb - youngest age 
l - length
jt - trusted age
xv,yv,wv - vectors with abundunce, catch, weight
"""
mutable struct cohort_v
    tmax
    jmax
    yb
    jb
    l
    jt
     xv::Array{Float64,1}
     yv::Array{Float64,1}
     wv::Array{Float64,1}
 end

"""
detail  representation of given cohorts
tmax - number of seasons
jmax - number of ages
cs -  vector of cohort of mutable struct cohort_v
"""
mutable struct cohorts_v
       tmax
    jmax
    cs::Array{cohort_v,1}
 end


"""
detail  representation of given cohorts
tmax - number of seasons
jmax - number of ages
cs -  vector of cohort of mutable struct cohort_v
"""

mutable struct cohorts_v1
    year_b
    tmax
    jmax
    cs::Array{cohort_v,1}
 end

"""
detail cohort representation 
no- the number
yb - physical year associated with cohort (when age=1)
tmax - number of seasons
jmax - number of ages
jbeg - initial age 
len - length
jtrust - trusted age
xv,yv,wv - vectors with abundunce, catch, weight
"""
mutable struct cohort_v1
    no
    yb
    tmax
    jmax
    jbeg
    len
    jtrust
     xv::Array{Float64,1}
     yv::Array{Float64,1}
     wv::Array{Float64,1}
 end

"""
cohorts as 3 tables (=cohorts_m)
yb - physical year associated with cohort (when age=1)
tmax - number of seasons
jmax - number of ages
x,y,w- tables with abundunce, catch, weight
"""
mutable struct  cohorts_t1  
     year_b
   tmax
    jmax
    x
     y
    w
end




"""cohorts as cs: cohort_v1?
year_b - physical year associated with cohort (when age=1)
tmax - number of seasons
jmax - number of ages
"""
mutable struct cohorts_v2
    year_b
    tmax
    jmax
    cs
end

"""
cohorts related matrix
year_b - initial physical year
tmax - number of seasons
jmax - number of ages
mtx: x,y,w,m,f, etc.
"""
mutable struct csmtx
    year_b
    tmax
    jmax
    mtx
end

"""
detail mtx representation of a cohort as sequence of vecs
no- the number
yb - physical year associated with cohort (when age=1)
jbeg - youngest age 
len - length
jtrust - trusted age
vec
"""
mutable struct cvec
    no
    yb
    jbeg
    len
    jtrust
    vec::Array{Float64,1}
     
 end

"""
cohort as vecs
"""
mutable struct csvecs
    year_b
    tmax
    jmax
    vecs
end



"""
cohorts_t1->cohorts_v2

cohorts_t1:
year_b - initial physical year
tmax - number of seasons
jmax - number of ages
x,y,w -  matrices  with abundunce, catch, weight


cohorts_v2:
year_b - initial physical year
tmax - number of seasons
jmax - number of ages
cs - cohorts of mutable struct cohort_v1

cohort_v1:
no- the number
yb - physical year associated with cohort (when age=1)
tmax - number of seasons
jmax - number of ages
jbeg - initial age 
len - length
jtrust - trusted age
xv,yv,wv- vectors with abundunce, catch, weight
tmax+jmax-1 -- number of related cohorts
"""
function cst2csv2(cs::cohorts_t1)
   x=cs.x
    y=cs.y
    w=cs.w
    tmax=cs.tmax
    jmax=cs.jmax
    year1=cs.year_b
    year0=year1-jmax+1
                #println(year0)
                #println(year1)
    #csv1=Array{cohort_v1}(tmax+jmax-1)
csv1=Array{cohort_v1}(undef,tmax+jmax-1)
    
    for year in year0:year1-1
        i=year-year0+1
                #println(year)
                #println(i)
       
        csv1[i]=cohort_v1(i,year,tmax,jmax,
        jmax-i+1,i,jmax-i+1,
        map(j->x[j,jmax-i+j],1:i),
        map(j->y[j,jmax-i+j],1:i),
        map(j->w[j,jmax-i+j],1:i))
                #println(csv1[i])
    end
    
    year2=year1+tmax-jmax
                #println(year2)
    
    if tmax>=jmax
        for year in year1:year2
            i=year-year1+1
               #println(year)
               #println(i)
               #println(i-year0+year1)
            csv1[i-year0+year1]=cohort_v1(i-year0+year1,year,tmax,jmax,
            1,jmax,1,  
            map(j->x[i+j-1,j],1:jmax),
            map(j->y[i+j-1,j],1:jmax),
            map(j->w[i+j-1,j],1:jmax))
               #println( csv1[i-year0+year1])
        end
    end
     
        
    
    year3=year2+jmax-1
                #println(year3)
    
    for year in year2+1:year3
            i=year-year1+1
                #println(year)
                #println(i)
                #println(i-year1+year2)
        csv1[i-year1+year2+1]=cohort_v1(i-year1+year2+1,year,tmax,jmax,
        1,tmax-i+1,1,  
        map(j->x[i+j-1,j],1:tmax-i+1),
        map(j->y[i+j-1,j],1:tmax-i+1),
        map(j->w[i+j-1,j],1:tmax-i+1))
                #println( csv1[i-year0+year1])
    end
    
    cohorts_v2(year1,tmax,jmax,csv1)
end

"""
cohorts_v2->cohirts_t1

cohorts_t1:
year_b - initial physical year
tmax - number of seasons
jmax - number of ages
x,y,w -  matrices  with abundunce, catch, weight


cohorts_v2:
year_b - initial physical year
tmax - number of seasons
jmax - number of ages
cs - cohorts of mutable struct cohort_v1

cohort_v1:
no- the number
yb - physical year associated with cohort (when age=1)
tmax - number of seasons
jmax - number of ages
jbeg - initial age 
len - length
jtrust - trusted age
xv,yv,wv- vectors with abundunce, catch, weight
"""
function csv2cst(csv::cohorts_v2)
    
    year1=csv.year_b
    tmax=csv.tmax
    jmax=csv.jmax
    cs=csv.cs
    
    
    x=fill(1.0,tmax,jmax)
    y=fill(1.0,tmax,jmax)
    w=fill(1.0,tmax,jmax)
    
    ncs=tmax+jmax-2
    year0=year1-jmax+1
    year2=year1+tmax-jmax
    year3=year2+jmax-1
    
    
    #println(year0)
    #println(year1)
    #println(year2)
    #println(year3)
    
   #1 
    for year in year0:year1-1
        i=year-year0+1
        #println(year)
        #println(i)
      
    cv1=cs[i]
       no=cv1.no
    yb=cv1.yb
    jbeg=cv1.jbeg
    len=cv1.len
        jtrust=cv1.jtrust
        xv=cv1.xv 
    yv=cv1.yv
    wv=cv1.wv
    
    
     #println(' ',i,' ',no,' ',yb,' ',jbeg,' ',len)   
     
        for j in 1:i
            #println(j,' ',jmax-i+j)
            x[j,jmax-i+j]=xv[j]
            y[j,jmax-i+j]=yv[j]
            w[j,jmax-i+j]=wv[j]
        end         
    end
  #2
     if tmax>=jmax
        for year in year1:year2
            i=year-year0+1
            i1=year-year1+1
                 #println(year)
                 #println(i)
        
            
            cv1=cs[i]
       no=cv1.no
    yb=cv1.yb
    jbeg=cv1.jbeg
    len=cv1.len
        jtrust=cv1.jtrust
        xv=cv1.xv 
    yv=cv1.yv
    wv=cv1.wv
            
            
               #println(i,' ',no,' ',yb,' ',jbeg,' ',len)   
           
            
            for j in 1:jmax
                #println(i+j-1,' ',j)
            x[i1+j-1,j]=xv[j]
            y[i1+j-1,j]=yv[j]
            w[i1+j-1,j]=wv[j]  
           
        end
    end
        end  
    
  #3  
  for year in year2+1:year3
            i=year-year0+1
        i1=year-year1+1
      
                #println(year)
                #println(i)
       
        
        
        cv1=cs[i]
       no=cv1.no
    yb=cv1.yb
    jbeg=cv1.jbeg
    len=cv1.len
        jtrust=cv1.jtrust
        xv=cv1.xv 
    yv=cv1.yv
    wv=cv1.wv
        
             #println(i,' ',no,' ',yb,' ',jbeg,' ',len)   

            
       
        for j in 1:tmax-i1+1
             #println(i1+j-1,' ',j)
            
            
        x[i1+j-1,j]=xv[j]
            y[i1+j-1,j]=yv[j]
            w[i1+j-1,j]=wv[j]
        end
    end
     
        
    
    
    cohorts_t1(year1,tmax,jmax,x,y,w)
end

"""
cohorts related matrix (of mutable struct csmtx: x,y,w,m,f,z,g)->sequence of vectors (of mutable struct csvecs)


csrm:cohorts related matrix (of mutable struct csrm)
year_b - initial physical year
tmax - number of seasons
jmax - number of ages
x,y,w, m,f,g -  cohorts related matrices  


cvec: detail cohort representation as vec
no- the number in the sequence
yb - physical year associated with cohort (when age=1)
jbeg - initial age 
len - length
jtrust - trusted age
vec

csvecs: cohorts related mtx as a sequence of cvec
no- the number
yb - physical year associated with cohort (when age=1)
tmax - number of seasons
jmax - number of ages
jbeg - initial age 
len - length
jtrust - trusted age
xv,yv,wv- vectors with abundunce, catch, weight
tmax+jmax-1 -- number of related cohorts
"""
function csmtx2csvs(csrmtx::csmtx)
   year1=csrmtx.year_b
    tmax=csrmtx.tmax
    jmax=csrmtx.jmax
    mtx=csrmtx.mtx
    
    year0=year1-jmax+1
      year2=year1+tmax-jmax
     year3=year2+jmax-1
    #println(year0)
    #println(year1)
    #println(year2)
    #println(year3)
    vecs=Array{cvec,1}(undef,tmax+jmax-1)
    
    for year in year0:year1-1
        i=year-year0+1
        #println(year)
        #println(i)
       
        vecs[i]=cvec(i,year,
        jmax-i+1,i,argmax(map(j->mtx[j,jmax-i+j],1:i)),
        map(j->mtx[j,jmax-i+j],1:i))
        #println(vecs[i])
    end
    
  
    
    if tmax>=jmax
        for year in year1:year2
            i=year-year1+1
             #println(year)
             #println(i)
        #println(i-year0+year1)
            vecs[i-year0+year1]=cvec(i-year0+year1,year,
            1,jmax, argmax(map(j->mtx[i+j-1,j],1:jmax)),  
            map(j->mtx[i+j-1,j],1:jmax))
            #println( vecs[i-year0+year1])
        end
    end
     
        
    
   #println(3)
    
    for year in year2+1:year3
            i=year-year1+1
            #println(year)
            #println(i)
            #println(i-year1+year2)
        vecs[i-year1+year2+1]=cvec(i-year1+year2+1,year,
        1,tmax-i+1,argmax(map(j->mtx[i+j-1,j],1:tmax-i+1)),  
        map(j->mtx[i+j-1,j],1:tmax-i+1))
        #println( vecs[i-year1+year2+1])
    end
    
    csvecs(year1,tmax,jmax,vecs)
end

"""
sequence of vectors (of mutable struct csvecs)->cohorts related matrix (of mutable struct csrm) (x,y,w,m,f,z,g)


csrm:cohorts related matrix (of mutable struct csrm)
year_b - initial physical year
tmax - number of seasons
jmax - number of ages
x,y,w, m,f,g -  cohorts related matrices  


cvec: detail cohort representation as vec
no- the number in the sequence
yb - physical year associated with cohort (when age=1)
jbeg - initial age 
len - length
jtrust - trusted age
vec

csvecs: cohorts related mtx as a sequence of cvec
no- the number
yb - physical year associated with cohort (when age=1)
tmax - number of seasons
jmax - number of ages
jbeg - initial age 
len - length
jtrust - trusted age
xv,yv,wv- vectors with abundunce, catch, weight
tmax+jmax-2 -- number of related cohorts
"""
function csvs2csmtx(csvs::csvecs)
    
    year1=csvs.year_b
    tmax=csvs.tmax
    jmax=csvs.jmax
    vs=csvs.vecs
    
    
    mtx=fill(1.0,tmax,jmax)
    
    
    ncs=tmax+jmax-1
    year0=year1-jmax+1
    year2=year1+tmax-jmax
    year3=year2+jmax-1
    
    
    #println(year0)
    #println(year1)
    #println(year2)
    #println(year3)
    
   #1 
    #println(1)
    for year in year0:year1-1
        i=year-year0+1
        #println(year)
        #println(i)
      
    cv1=vs[i]
    no=cv1.no
    yb=cv1.yb
    jbeg=cv1.jbeg
    len=cv1.len
    jtrust=cv1.jtrust
    xv=cv1.vec
    
    
    
     #println(' ',i,' ',no,' ',yb,' ',jbeg,' ',len)   
     
        for j in 1:i
            #println(j,' ',jmax-i+j)
            mtx[j,jmax-i+j]=xv[j]
            
        end         
    end
  #2
    #println(2)
     if tmax>=jmax
        for year in year1:year2
            i=year-year0+1
            i1=year-year1+1
            #println(year)
            #println(i)
        
            
    cv1=vs[i]
    no=cv1.no
    yb=cv1.yb
    jbeg=cv1.jbeg
    len=cv1.len
    jtrust=cv1.jtrust
    xv=cv1.vec 
    
            
            
           #println(i,' ',no,' ',yb,' ',jbeg,' ',len)   
           
            
            for j in 1:jmax
           # println(i+j-1,' ',j)
            mtx[i1+j-1,j]=xv[j]
              
           
        end
    end
        end  
    
  #3  
  #println(3)
    for year in year2+1:year3
            i=year-year0+1
        i1=year-year1+1
      
            #println(year)
            #println(i)
       
        
        
    cv1=vs[i]
    no=cv1.no
    yb=cv1.yb
    jbeg=cv1.jbeg
    len=cv1.len
    jtrust=cv1.jtrust
     xv=cv1.vec 
   
        
            #println(i,' ',no,' ',yb,' ',jbeg,' ',len)   

            
       
        for j in 1:tmax-i1+1
        # println(i1+j-1,' ',j)
            
            
        mtx[i1+j-1,j]=xv[j]
            
        end
    end
     
        
    
    
    csmtx(year1,tmax,jmax,mtx)
end


"""
generates xv1,yv1 of csvecs mutable struct with use of xv::csvecs,yv::csvecs,wv::csvecs,zv::csvecs,gv::csvecs
first makes a deep copy, then modifies vectors 
"""
function gvecs(xv::csvecs,yv::csvecs,wv::csvecs,zv::csvecs,gv::csvecs)
    xv1=deepcopy(xv)
    yv1=deepcopy(yv)
    
    
    year1=xv.year_b
    tmax=xv.tmax
    jmax=xv.jmax
    xvecs=xv.vecs
    yvecs=yv.vecs
    wvecs=yv.vecs
    zvecs=zv.vecs
    gvecs=gv.vecs
    
    x1vecs=xv1.vecs
    y1vecs=yv1.vecs
    
    
 
    
    #println(xv1)
    
    
    
    
    # println( mvecs)
 #1   
    for i in 2:jmax-1
        xvi=xvecs[i]
        yvi=yvecs[i]
        wvi=wvecs[i]
        zvi=zvecs[i]
        gvi=gvecs[i]
        
        x1vi=x1vecs[i]
        y1vi=y1vecs[i]
      #  println(mvi)
       
        
        
         no=xvi.no
        yb=xvi.yb
        jbeg=xvi.jbeg
         len=xvi.len
        jtr=xvi.jtrust
        xi=xvi.vec
        
        yi=yvi.vec
        wi=wvi.vec
        zi=zvi.vec
        gi=gvi.vec
        
        x1i=x1vi.vec
        y1i=y1vi.vec
        
        
        #println(jbeg,' ',jtr)
        
        xijtr=xi[jtr]
              
        
        x1vecs[i].vec[jtr]=xijtr
        
        
        #println(xijtr)
        if jtr >1
        for i1 in 1:jtr-1
        #println(i1)
            x1vecs[i].vec[i1]=xijtr*exp(sum(zi[i1:jtr-1]))
        end
        end
        if jtr+1<len  
        for i1 in jtr+1:len
          #println(i1)
            x1vecs[i].vec[i1]=xijtr*exp(-sum(zi[jtr+1:i1]))
                y1vecs[i].vec[i1]= wi[i1]*x1vecs[i].vec[i1]*gi[i1]  
        end
        end
            
        #println(x1vecs[i])  
                  
    end
   #2 
   for i in jmax:tmax
        xvi=xvecs[i]
        yvi=yvecs[i]
        wvi=wvecs[i]
        zvi=zvecs[i]
        gvi=gvecs[i]
        
        x1vi=x1vecs[i]
        y1vi=y1vecs[i]
      #  println(mvi)
       
        
        
         no=xvi.no
        yb=xvi.yb
        jbeg=xvi.jbeg
         len=xvi.len
        jtr=xvi.jtrust
        xi=xvi.vec
        
        yi=yvi.vec
        wi=wvi.vec
        zi=zvi.vec
        gi=gvi.vec
        
        x1i=x1vi.vec
        y1i=y1vi.vec
        
        
        #println(jbeg,' ',jtr)
        
        xijtr=xi[jtr]
              
        
        x1vecs[i].vec[jtr]=xijtr
        
        
        #println(xijtr)
        if jtr >1
        for i1 in 1:jtr-1
         #println(i1)
            x1vecs[i].vec[i1]=xijtr*exp(sum(zi[i1:jtr-1]))
        end
        end
        if jtr+1<len  
        for i1 in jtr+1:len
          #println(i1)
            x1vecs[i].vec[i1]=xijtr*exp(-sum(zi[jtr+1:i1]))
                 y1vecs[i].vec[i1]= wi[i1]*x1vecs[i].vec[i1]*gi[i1]  
        end
        end
            
        #println(x1vecs[i])  
                  
    end 
    #3
    for i in tmax+1:tmax+jmax-2
        xvi=xvecs[i]
        yvi=yvecs[i]
        wvi=wvecs[i]
        zvi=zvecs[i]
        gvi=gvecs[i]
        
        x1vi=x1vecs[i]
        y1vi=y1vecs[i]
      #  println(mvi)
       
        
        
         no=xvi.no
        yb=xvi.yb
        jbeg=xvi.jbeg
         len=xvi.len
        jtr=xvi.jtrust
        xi=xvi.vec
        
        yi=yvi.vec
        wi=wvi.vec
        zi=zvi.vec
        gi=gvi.vec
        
        x1i=x1vi.vec
        y1i=y1vi.vec
        
        
        #println(jbeg,' ',jtr)
        
        xijtr=xi[jtr]
              
        
        x1vecs[i].vec[jtr]=xijtr
        
        
        #println(xijtr)
        if jtr >1
        for i1 in 1:jtr-1
         #println(i1)
            x1vecs[i].vec[i1]=xijtr*exp(sum(zi[i1:jtr-1]))
        end
        end
        if jtr+1<len  
        for i1 in jtr+1:len
          #println(i1)
            x1vecs[i].vec[i1]=xijtr*exp(-sum(zi[jtr+1:i1]))
                 y1vecs[i].vec[i1]= wi[i1]*x1vecs[i].vec[i1]*gi[i1]  
        end
        end
            
        #println(x1vecs[i])  
        #println(y1vecs[i]) 
                  
    end 
    
    
    return(xv1,yv1)  
    
end
    
"""
generates xv1,yv1 of csvecs mutable struct with use of xv::csvecs,yv::csvecs,wv::csvecs,zv::csvecs,gv::csvecs
first makes a deep copy, then modifies vectors 
"""
function gvecs1(xv::csvecs,yv::csvecs,wv::csvecs,zv::csvecs,gv::csvecs)
    xv1=deepcopy(xv)
    yv1=deepcopy(yv)
    
    
    year1=xv.year_b
    tmax=xv.tmax
    jmax=xv.jmax
    xvecs=xv.vecs
    yvecs=yv.vecs
    wvecs=yv.vecs
    zvecs=zv.vecs
    gvecs=gv.vecs
    
    x1vecs=xv1.vecs
    y1vecs=yv1.vecs
    
    
  
    #println(xv1)
    
    
    
    
    # println( mvecs)
 #1   
    for i in 2:jmax-1
        xvi=xvecs[i]
        yvi=yvecs[i]
        wvi=wvecs[i]
        zvi=zvecs[i]
        gvi=gvecs[i]
        
        x1vi=x1vecs[i]
        y1vi=y1vecs[i]
      #  println(mvi)
       
        
        
         no=xvi.no
        yb=xvi.yb
        jbeg=xvi.jbeg
         len=xvi.len
        jtr=xvi.jtrust
        xi=xvi.vec
        
        yi=yvi.vec
        wi=wvi.vec
        zi=zvi.vec
        gi=gvi.vec
        
        x1i=x1vi.vec
        y1i=y1vi.vec
        
        
        #println(jbeg,' ',jtr)
        
        xijtr=xi[jtr]
              
        
        x1vecs[i].vec[jtr]=xijtr
        
        
        #println(xijtr)
        if jtr >1
        for i1 in 1:jtr-1
        #println(i1)
            x1vecs[i].vec[i1]=xijtr*exp(sum(zi[i1:jtr-1]))
             
        end
        end
        if jtr+1<len  
        for i1 in jtr+1:len
          #println(i1)
            x1vecs[i].vec[i1]=xijtr*exp(-sum(zi[jtr+1:i1]))
               
        end
        end
            
        #println(x1vecs[i])  
                  
    end
    
    
    
   #2 
   for i in jmax:tmax
        xvi=xvecs[i]
        yvi=yvecs[i]
        wvi=wvecs[i]
        zvi=zvecs[i]
        gvi=gvecs[i]
        
        x1vi=x1vecs[i]
        y1vi=y1vecs[i]
      #  println(mvi)
       
        
        
         no=xvi.no
        yb=xvi.yb
        jbeg=xvi.jbeg
         len=xvi.len
        jtr=xvi.jtrust
        xi=xvi.vec
        
        yi=yvi.vec
        wi=wvi.vec
        zi=zvi.vec
        gi=gvi.vec
        
        x1i=x1vi.vec
        y1i=y1vi.vec
        
        
        #println(jbeg,' ',jtr)
        
        xijtr=xi[jtr]
              
        
        x1vecs[i].vec[jtr]=xijtr
        
        
        #println(xijtr)
        if jtr >1
        for i1 in 1:jtr-1
         #println(i1)
            x1vecs[i].vec[i1]=xijtr*exp(sum(zi[i1:jtr-1]))
                
        end
        end
        if jtr+1<len  
        for i1 in jtr+1:len
          #println(i1)
            x1vecs[i].vec[i1]=xijtr*exp(-sum(zi[jtr+1:i1]))
                 
        end
        end
            
        #println(x1vecs[i])  
                  
    end 
    #3
    for i in tmax+1:tmax+jmax-2
        xvi=xvecs[i]
        yvi=yvecs[i]
        wvi=wvecs[i]
        zvi=zvecs[i]
        gvi=gvecs[i]
        
        x1vi=x1vecs[i]
        y1vi=y1vecs[i]
      #  println(mvi)
       
        
        
         no=xvi.no
        yb=xvi.yb
        jbeg=xvi.jbeg
         len=xvi.len
        jtr=xvi.jtrust
        xi=xvi.vec
        
        yi=yvi.vec
        wi=wvi.vec
        zi=zvi.vec
        gi=gvi.vec
        
        x1i=x1vi.vec
        y1i=y1vi.vec
        
        
        #println(jbeg,' ',jtr)
        
        xijtr=xi[jtr]
              
        
        x1vecs[i].vec[jtr]=xijtr
        
        
        #println(xijtr)
        if jtr >1
        for i1 in 1:jtr-1
         #println(i1)
            x1vecs[i].vec[i1]=xijtr*exp(sum(zi[i1:jtr-1]))
                
        end
        end
        if jtr+1<len  
        for i1 in jtr+1:len
          #println(i1)
            x1vecs[i].vec[i1]=xijtr*exp(-sum(zi[jtr+1:i1]))

        end
        end
            
        #println(x1vecs[i])  
        #println(y1vecs[i]) 
                  
    end 
   
    
    
    
    
    ######
    
     for i in 1:jmax-1
        xvi=xvecs[i]
        yvi=yvecs[i]
        wvi=wvecs[i]
        zvi=zvecs[i]
        gvi=gvecs[i]
        
        x1vi=x1vecs[i]
        y1vi=y1vecs[i]
      #  println(mvi)
            
        
         no=xvi.no
        yb=xvi.yb
        jbeg=xvi.jbeg
         len=xvi.len
        jtr=xvi.jtrust
        xi=xvi.vec
        
        yi=yvi.vec
        wi=wvi.vec
        zi=zvi.vec
        gi=gvi.vec
        
        x1i=x1vi.vec
        y1i=y1vi.vec
        
       
    
        
        #println(jbeg,' ',jtr)
        
        xijtr=xi[jtr]
              
        
        x1vecs[i].vec[jtr]=xijtr
        
        
        #println(xijtr)
        
        for i1 in 1:len
        #println(i1)
           y1vecs[i].vec[i1]= wi[i1]*x1vecs[i].vec[i1]*gi[i1]
             
        end
                end
            
    
    
    
   #2 
   for i in jmax:tmax
        xvi=xvecs[i]
        yvi=yvecs[i]
        wvi=wvecs[i]
        zvi=zvecs[i]
        gvi=gvecs[i]
        
        x1vi=x1vecs[i]
        y1vi=y1vecs[i]
      #  println(mvi)
       
        
        
         no=xvi.no
        yb=xvi.yb
        jbeg=xvi.jbeg
         len=xvi.len
        jtr=xvi.jtrust
        xi=xvi.vec
        
        yi=yvi.vec
        wi=wvi.vec
        zi=zvi.vec
        gi=gvi.vec
        
        x1i=x1vi.vec
        y1i=y1vi.vec
        
        
        #println(jbeg,' ',jtr)
        
        xijtr=xi[jtr]
              
        
        x1vecs[i].vec[jtr]=xijtr
        
        
        #println(xijtr)
       for i1 in 1:len
        #println(i1)
           y1vecs[i].vec[i1]= wi[i1]*x1vecs[i].vec[i1]*gi[i1]
             
        end
                end
   
    #3
    for i in tmax+1:tmax+jmax-1
        xvi=xvecs[i]
        yvi=yvecs[i]
        wvi=wvecs[i]
        zvi=zvecs[i]
        gvi=gvecs[i]
        
        x1vi=x1vecs[i]
        y1vi=y1vecs[i]
      #  println(mvi)
       
        
        
         no=xvi.no
        yb=xvi.yb
        jbeg=xvi.jbeg
         len=xvi.len
        jtr=xvi.jtrust
        xi=xvi.vec
        
        yi=yvi.vec
        wi=wvi.vec
        zi=zvi.vec
        gi=gvi.vec
        
        x1i=x1vi.vec
        y1i=y1vi.vec
        
        
        #println(jbeg,' ',jtr)
        
        xijtr=xi[jtr]
              
        
        x1vecs[i].vec[jtr]=xijtr
        
        
       for i1 in 1:len
        #println(i1)
           y1vecs[i].vec[i1]= wi[i1]*x1vecs[i].vec[i1]*gi[i1]
             
        end
                end
    
    
    
    
    
    return(xv1,yv1)  
    
end
    


#########################################










"""
look how to evaluate closeness of data and model!!!!
as L1+L2
no prints out
"""
function Lp(cst,a1m,a2m,b1m,b2m,cm,d1m,d2m,htj,as,bs,ft,kx,kw,p)
    
    function g(hh,mm,ss,ff)
   ss*ff/(hh*mm+ss*ff)*(1-exp(-hh*mm-ss*ff))
    end
    
    
    local  x=cst.x
    
      local  y=cst.y
     local   w=cst.w
 
     (tmax,jmax)=size(x)
    
   #println(tmax," ",jmax)
   # println(a1m,a2m,b1m,b2m,cm,d1m,d2m,as,bs)
   #println(htj)
   #println(ft)
     local  L1=0.0
   
    for t=1:tmax-1
        L1t=0.0
        for j=1:jmax-1
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
            L1t=L1t+abs(log(x[t,j]/x[t+1,j+1])-htj[t,j]*mj-sj*ft[t])^p
           #println(t," ",j," ",L1t) 
        end 
        L1=L1+L1t
      
        
    end
  
   #println(L1)
      
    
     local   L2=0.0
    for t=1:tmax        
        L2t=0.0
        for j= 1:jmax
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
           
            L2t=L2t+abs(log(y[t,j])-log(kw*w[t,j]*kx*x[t,j]*g(htj[t,j],mj,sj,ft[t])))^p
            
           #println(t," ",j," ",L2t) 
        end;
        L2=L2+L2t;
               
    end
    #println(L2)
  # println(L1+L2)
    
    
    if L1+L2==Inf 10000.0 else L1+L2 end
end


"""
L as L1+L2 with dispersions for the normal distribution, p=2
no prints out
"""
function Lmine2(cst,a1m,a2m,b1m,b2m,cm,d1m,d2m,htj,as,bs,ft,kx,kw,p)
    
    function g(hh,mm,ss,ff)
   ss*ff/(hh*mm+ss*ff)*(1-exp(-hh*mm-ss*ff))
    end
    
    
    local  x=cst.x
    
      local  y=cst.y
     local   w=cst.w
 
     (tmax,jmax)=size(x)
    
   #println(tmax," ",jmax)
   # println(a1m,a2m,b1m,b2m,cm,d1m,d2m,as,bs)
   #println(htj)
   #println(ft)
     local  L1=0.0
    local  L11=0.0
    local mean1
     
    
   
    for t=1:tmax-1
        L1t=0.0
        for j=1:jmax-1
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
            L1t=L1t+
                    log(x[t,j]/x[t+1,j+1])-htj[t,j]*mj-sj*ft[t]
           #println(t," ",j," ",L1t) 
        end 
      L1=L1+L1t  
     end
        
        mean1=L1/((tmax-1)*(jmax-1))
        
     for t=1:tmax-1
        L1t=0.0
        for j=1:jmax-1
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
            L1t=L1t+
                    (log(x[t,j]/x[t+1,j+1])-htj[t,j]*mj-sj*ft[t]-mean1)^2
           #println(t," ",j," ",L1t) 
        end 
        L11=L11+L1t
      end  
    
  
   #println(L1)
      
    
     local   L2=0.0
        local  L22=0.0
        local mean2
    
    for t=1:tmax        
        L2t=0.0
        for j= 1:jmax
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
           
            L2t=L2t+
                log(y[t,j])-log(kw*w[t,j]*kx*x[t,j]*g(htj[t,j],mj,sj,ft[t]))
            
           #println(t," ",j," ",L2t) 
        end;
        L2=L2+L2t;
               
    end
    
     mean2=L2/(tmax*jmax)
    
    for t=1:tmax        
        L2t=0.0
        for j= 1:jmax
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
           
            L2t=L2t+
                    (log(y[t,j])-log(kw*w[t,j]*kx*x[t,j]*g(htj[t,j],mj,sj,ft[t]))-mean2)^2
            
           #println(t," ",j," ",L2t) 
        end;
        L22=L22+L2t;
               
    end
    
    
    
    #println(L2)
  # println(L1+L2)
    local result
    result=lambda *L11/((tmax-1)*(jmax-1)-1)+(1. - lambda)*L22/(tmax*jmax)
    
    if result==Inf 10000.0 else result end
end



"""
L as max
no prints out
"""
function Lmine2max(cst,a1m,a2m,b1m,b2m,cm,d1m,d2m,htj,as,bs,ft,kx,kw,p)
    
    function g(hh,mm,ss,ff)
   ss*ff/(hh*mm+ss*ff)*(1-exp(-hh*mm-ss*ff))
    end
    
    
    local  x=cst.x
    
      local  y=cst.y
     local   w=cst.w
 
     (tmax,jmax)=size(x)
    
   #println(tmax," ",jmax)
   # println(a1m,a2m,b1m,b2m,cm,d1m,d2m,as,bs)
   #println(htj)
   #println(ft)
     local  L1=0.0
    local  L11=0.0
    local mean1
     
    
   
    for t=1:tmax-1
        L1t=0.0
        for j=1:jmax-1
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
            L1t=L1t+
                    log(x[t,j]/x[t+1,j+1])-htj[t,j]*mj-sj*ft[t]
           #println(t," ",j," ",L1t) 
        end 
      L1=L1+L1t  
     end
        
        mean1=L1/((tmax-1)*(jmax-1))
        
     for t=1:tmax-1
        L1t=0.0
        for j=1:jmax-1
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
            L1t=L1t+
                    (log(x[t,j]/x[t+1,j+1])-htj[t,j]*mj-sj*ft[t]-mean1)^2
           #println(t," ",j," ",L1t) 
        end 
        L11=L11+L1t
      end  
    
  
   #println(L1)
      
    
     local   L2=0.0
        local  L22=0.0
        local mean2
    
    for t=1:tmax        
        L2t=0.0
        for j= 1:jmax
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
           
            L2t=L2t+
                log(y[t,j])-log(kw*w[t,j]*kx*x[t,j]*g(htj[t,j],mj,sj,ft[t]))
            
           #println(t," ",j," ",L2t) 
        end;
        L2=L2+L2t;
               
    end
    
     mean2=L2/(tmax*jmax)
    
    for t=1:tmax        
        L2t=0.0
        for j= 1:jmax
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
           
            L2t=L2t+
                    (log(y[t,j])-log(kw*w[t,j]*kx*x[t,j]*g(htj[t,j],mj,sj,ft[t]))-mean2)^2
            
           #println(t," ",j," ",L2t) 
        end;
        L22=L22+L2t;
               
    end
    
    
    
    #println(L2)
  # println(L1+L2)
    local result
    global lambda
    
    
    
    result=max(lambda*L11/((tmax-1)*(jmax-1)-1),(1. - lambda)*L22/(tmax*jmax))
    
    if result==Inf 10000.0 else result end
end





"""
look how to evaluate closeness of data and model!!!!
as L1+L2
"""
function Lp12(cst,a1m,a2m,b1m,b2m,cm,d1m,d2m,htj,as,bs,ft,kx,kw,p)
    
    function g(hh,mm,ss,ff)
   ss*ff/(hh*mm+ss*ff)*(1-exp(-hh*mm-ss*ff))
    end
    
    
    local  x=cst.x
    
      local  y=cst.y
     local   w=cst.w
 
     (tmax,jmax)=size(x)
    
   #println(tmax," ",jmax)
   # println(a1m,a2m,b1m,b2m,cm,d1m,d2m,as,bs)
   #println(htj)
   #println(ft)
     local  L1=0.0
   
    for t=1:tmax-1
        L1t=0.0
        for j=1:jmax-1
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
            L1t=L1t+abs(log(x[t,j]/x[t+1,j+1])-htj[t,j]*mj-sj*ft[t])^p
           #println(t," ",j," ",L1t) 
        end 
        L1=L1+L1t
      
        
    end
  
   println(L1)
      
    
     local   L2=0.0
    for t=1:tmax        
        L2t=0.0
        for j= 1:jmax
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
           
            L2t=L2t+abs(log(y[t,j])-log(kw*w[t,j]*kx*x[t,j]*g(htj[t,j],mj,sj,ft[t])))^p
            
           #println(t," ",j," ",L2t) 
        end;
        L2=L2+L2t;
               
    end
    println(L2)
   println(L1+L2)
    
    
    if L1+L2==Inf 10000.0 else L1+L2 end
end


function Lp1(cst,a1m,a2m,b1m,b2m,cm,d1m,d2m,htj,as,bs,ft,kx,kw,p)
    
    function g(hh,mm,ss,ff)
   ss*ff/(hh*mm+ss*ff)*(1-exp(-hh*mm-ss*ff))
    end
    
    
    local  x=cst.x
    
      local  y=cst.y
     local   w=cst.w
 
     (tmax,jmax)=size(x)
    
   #println(tmax," ",jmax)
   # println(a1m,a2m,b1m,b2m,cm,d1m,d2m,as,bs)
   #println(htj)
   #println(ft)
     local  L1=0.0
   
    for t=1:tmax-1
        L1t=0.0
        for j=1:jmax-1
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
            L1t=L1t+abs(log(x[t,j]/x[t+1,j+1])-htj[t,j]*mj-sj*ft[t])^p
           #println(t," ",j," ",L1t) 
        end 
        L1=L1+L1t
      
        
    end
  
   println(L1)
      
    
     local   L2=0.0
    for t=1:tmax        
        L2t=0.0
        for j= 1:jmax
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
           
            L2t=L2t+abs(log(y[t,j])-log(kw*w[t,j]*kx*x[t,j]*g(htj[t,j],mj,sj,ft[t])))^p
            
           #println(t," ",j," ",L2t) 
        end;
        L2=L2+L2t;
               
    end
    println(L2)
   println(L1+L2)
    
    
    if L1+L2==Inf 10000.0 else L1+L2 end
end


function Lp2(cst,a1m,a2m,b1m,b2m,cm,d1m,d2m,htj,as,bs,ft,kx,kw,p)
    
    function g(hh,mm,ss,ff)
   ss*ff/(hh*mm+ss*ff)*(1-exp(-hh*mm-ss*ff))
    end
    
    
    local  x=cst.x
    
      local  y=cst.y
     local   w=cst.w
 
     (tmax,jmax)=size(x)
    
   #println(tmax," ",jmax)
   # println(a1m,a2m,b1m,b2m,cm,d1m,d2m,as,bs)
   #println(htj)
   #println(ft)
     local  L1=0.0
   
    for t=1:tmax-1
                for j=1:jmax-1
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
            L1=L1+abs(log(x[t,j]/x[t+1,j+1])-htj[t,j]*mj-sj*ft[t])^p
           #println(t," ",j," ",L1t) 
        end 
    
      
        
    end
  
   println(L1)
      
    
     local   L2=0.0
    for t=1:tmax        
              for j= 1:jmax
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
           
            L2=L2+abs(log(y[t,j])-log(kw*w[t,j]*kx*x[t,j]*g(htj[t,j],mj,sj,ft[t])))^p
            
           #println(t," ",j," ",L2t) 
        end;
                      
    end
    println(L2)
   println(L1+L2)
    
    
    if L1+L2==Inf 10000.0 else L1+L2 end
end



"""
generates x,y with M,F,W,kx,kw,x1j,xt1
"""
function gm1(tmax,jmax,M,F,W,kx,kw,x1j,xt1)
##0
    x=fill(0.0,tmax,jmax);
   
    for j in [1:jmax;] x[1,j]=x1j[j] end
    
    for t in 1:tmax x[t,1]=xt1[t] end
   
 ##1   
        for t in [1:tmax-jmax+1;]
           for i in [0:jmax-2;]
            #println(t+i+1," ",i+2,"<-", t+i," ",i+1)
            x[t+i+1,i+1+1]=x[t+i,i+1]*exp(-M[t+i,i+1]-F[t+i,i+1])
           end
        end
 
 ##2   
    for t in [tmax-jmax+2:tmax-1;]
           for i in [0:tmax-t-1;]
            #println(t+i+1," ",i+2,"<-", t+i," ",i+1)
            x[t+i+1,i+1+1]=x[t+i,i+1]*exp(-M[t+i,i+1]-F[t+i,i+1])
           end
        end
   ##3
    
    for j in [2:jmax;]
           for i in [0:jmax-j-1;]
           #println(i+1+1,"! ",j+i+1,"<-", i+1," ",j+i)
            x[i+1+1,j+i+1]=x[i+1,j+i]*exp(-M[i+1,j+i]-F[i+1,j+i])
           end
        end
    
    
    
 ##0   
     y=fill(0.0,tmax,jmax);
  
##1
   for t in 1:tmax
        for j in 1:jmax 
            #println(t+i," ",i+1)
            Z=M[t,j]+F[t,j]
            y[t,j]=
            (kw*W[t,j])*(kx*x[t,j])*F[t,j]/Z*(1- exp(-Z))
        end
    end
           
    return x,y
end



"""
returns ovsyannikov's weight for j 
"""
function ww_osp(t,j)
   wv=[0.012 0.052 0.11 0.179 0.253 0.323 0.483 0.677 0.908 1.116 1.297 1.541 1.758 1.971 1.997 2.358 2.475];
    return wv[j]
end

"""
returns total mortality z as sum of natural m and fishing f mortality matrices
"""
function gz(m::Array{Float64,2},f::Array{Float64,2},tmax,jmax)
    z=Array{Float64,2}(tmax,jmax)
    z=m+f
    return(z)
end

"""
generates g matrix with m and f, tmax,jmax
"""
function gg(m::Array{Float64,2},f::Array{Float64,2},tmax,jmax)
    z=Array{Float64,2}(tmax,jmax)
    g=Array{Float64,2}(tmax,jmax)
    z=gz(m,f,tmax,jmax)
    for t in 1:tmax
        for j in 1:jmax
            g[t,j]=f[t,j]/z[t,j]*(1-exp(-z[t,j]))
        end
    end
    return(g)
end




"""
evaluates criteria
"""
function Lpv(cst,a1m,a2m,b1m,b2m,cm,d1m,d2m,htj,as,bs,ft,kx,kw,p)
    local  x=cst.x
    
      local  y=cst.y
     local   w=cst.w
   local year1= cst.year_b
 
     (tmax,jmax)=size(x)
    
  
     M0=gmtj(a1m,a2m,b1m,b2m,cm,d1m,d2m,htj,tmax,jmax);
    F0=gftj(as,bs,ft,tmax,jmax);
    Z0=M0+F0;
    G0=gg(M0,F0,tmax,jmax);
    
    xv=csmtx2csvs(csmtx(year1,tmax,jmax,x));
yv=csmtx2csvs(csmtx(year1,tmax,jmax,y));
    mv=csmtx2csvs(csmtx(year1,tmax,jmax,M0));
fv=csmtx2csvs(csmtx(year1,tmax,jmax,F0));
zv=csmtx2csvs(csmtx(year1,tmax,jmax,Z0));
gv=csmtx2csvs(csmtx(year1,tmax,jmax,G0));
wv=csmtx2csvs(csmtx(year1,tmax,jmax,w));
 xv0,yv0=gvecs(xv,yv,wv,zv,gv);   
    
   x11=csvs2csmtx(xv0);
y11=csvs2csmtx(yv0); 
    
    
    L=0;
for i in 1:tmax+jmax-1
    lenxv0=xv0.vecs[i].len;
   vecxv0=xv0.vecs[i].vec;
    vecyv0=yv0.vecs[i].vec;
    vecxv=xv.vecs[i].vec;
    vecyv=yv.vecs[i].vec;
           for n in 1:lenxv0
        L=L+abs(log(vecxv0[n])-log(vecxv[n]))^p+abs(log(vecyv0[n])-log(vecyv[n]))^p
    end
    
end

    
    
    if L==Inf 10000.0 else L end
end

"""
optimizes Lpv
"""
function optimizeLp_v(x,tmax,jmax,a,lb,lb1,xtol,xtol1,bc4,ctol,is1,is2,is3,pop,pop1,maxe)
   
   global count = 0
opt = Opt(ALGS[a], 4*tmax+9)
#l_opt= Opt(ALGS[a1], 4*tmax+9)
    lower_bounds!(opt, fill(lb,4*tmax+9))
    #lower_bounds!(l_opt, fill(lb1,4*tmax+9))
    xtol_rel!(opt,xtol)
    #xtol_rel!(l_opt,xtol1)

min_objective!(opt, buragofunc_v)
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
returns Lpv
"""
function buragofunc_v(x::Vector, grad::Vector)
 
    global cst
    (tmax,jmax)=size(cst.x)
    
    local Hl=fill(1.0,tmax,jmax)
    local ftl=fill(1.,tmax)
    
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
    return Lpv(cst,a1ml,a2ml,b1ml,b2ml,cml,d1ml,d2ml,Hl,asl,bsl,ftl,kx,kw,1)
   
    #println("f_$count($x)=$y")

    
end

