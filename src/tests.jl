"""
checks if direct*reverse=equal for cst2csv2 and csv2cst
"""
function runtest_cst2csv(yb,tmax,jmax,x,y,w)
cst=cohorts_t1(yb,tmax,jmax,x,y,w); 
    @test eq_cst(cst, csv2cst(cst2csv2(cst)))
end

"""
checks if generated with M, F matrix x corresponds to model relations
"""
function testgmx(x,M,F)
  tmax,jmax = size(x)
  out = true
  for t in 1:tmax-1, j in 1:jmax-1
        if t+1<=tmax&&j+1<=jmax out=out && x[t+1,j+1]==x[t,j]*exp(-M[t,j]-F[t,j])
        end
    end
  out
end





"""
checks if  2 cohorts of cvec type equal
"""
function eq_cvec(cv1,cv2)
    typeof(cv1)==  typeof(cv2)&&
    cv1.no==cv2.no &&
    cv1.yb==cv2.yb &&
    cv1.jbeg==cv2.jbeg &&
    cv1.len==cv2.len&&
    cv1.jtrust==cv2.jtrust &&
    cv1.vec==cv2.vec
       
end

"""
checks if  2 descriptions of csmtx type equal
"""
function eq_csvecs(csvs1,csvs2)
    result=typeof(csvs1)==  typeof(csvs2)&&
   csvs1.year_b==csvs2.year_b &&
   csvs1.tmax==csvs2.tmax &&
   csvs1.jmax==csvs2.jmax 
        for i in 1:length(csvs1.vecs)
           result=result && eq_cvec(csvs1.vecs[i],csvs2.vecs[i])
    end
    return(result)
end

"""
checks if csv2cst(cst2csv2(cst)=cst
"""
function runtest_csmtx2csvs(yb,ymax,jmax,mtx)
cst=cohorts_t1(yb,ymax,jmax,mtx); 
    @test eq_cst(cst, csv2cst(cst2csv2(cst)))
end

"""
checks if generated with parameters matrix y corresponds  to model relations
"""
function testgmy1(x,y,a1m,a2m,b1m,b2m,cm,d1m,d2m,htj,w,as,bs,ft,kx,kw,eps)
     function g(hh,mm,ss,ff)
   ss*ff/(hh*mm+ss*ff)*(1-exp(-hh*mm-ss*ff))
    end
  tmax,jmax = size(y)
  out = true
  for t in 1:tmax, j in 1:jmax
        if t<=tmax&&j<=jmax 
            sj= sigma(j,as,bs,jmax)
            mj=ushape(j,a1m,a2m,b1m,b2m,cm,d1m,d2m,jmax)
                    
            out=out && 
            abs(y[t,j]-kw*w[t,j]*kx*x[t,j]*g(htj[t,j],mj,sj,ft[t]))<eps
            if ! out println("t=",t,"  j=",j,
                    " ytj=",y[t,j], "     ", kw*w[t,j]*kx*x[t,j]*g(htj[t,j],mj,sj,ft[t])) end
        end
    end
  return(out)
end

"""
checks if generated with M, F, W matrix y corresponds to model relations
"""
function testgmy(x,y,M,F,W,kx,kw)
  tmax,jmax = size(y)
  out = true
  for t in 1:tmax, j in 1:jmax
        if t<=tmax&&j<=jmax out=out && 
            y[t,j]==(kw*W[t,j])*(kx*x[t,j])*F[t,j]/(M[t,j]+F[t,j])*(1- exp(-M[t,j]-F[t,j]))
        end
    end
  return(out)
end

##############################################
"""
checks if  2 cohorts of cohorts_t type equal
"""
function eq_cst(cst1,cst2)
    typeof(cst1)== typeof(cst2)&&
    cst1.year_b==cst2.year_b &&
    cst1.tmax==cst2.tmax &&
    cst1.jmax==cst2.jmax &&
    cst1.x==cst2.x&&
    cst1.y==cst2.y &&
    cst1.w==cst2.w     
end

"""
checks if  2 cohorts of cohort_v1 type equal
"""
function eq_cv1(cv1,cv2)
    typeof(cv1)==  typeof(cv2)&&
    cv1.no==cv2.no &&
    cv1.yb==cv2.yb &&
    cv1.tmax==cv2.tmax &&
    cv1.jmax==cv2.jmax &&
    cv1.jbeg==cv2.jbeg &&
    cv1.len==cv2.len&&
    cv1.jtrust==cv2.jtrust &&
    cv1.xv==cv2.xv&&
    cv1.yv==cv2.yv&&
    cv1.wv==cv2.wv     
end

"""
checks if  vectors of cohorts of cohort_v1 type equal alt
"""
function eq_cv1_alt(csv1,csv2)
        result=typeof(csv1)==typeof(csv2)
    for i in 1:length(csv1)
       result= result && eq_cv1(csv1[i],csv2[i])
    end
    result=result && length(csv1)==length(csv2)
    
end

"""
checks if  2 cohorts of cohorts_v2 type equal
"""
function eq_csv2(csv1,csv2)
    typeof(csv1)==  typeof(csv2)&&
    csv1.year_b==csv2.year_b &&
    csv1.tmax==csv2.tmax &&
    csv1.jmax==csv2.jmax &&
       eq_csv1(csv1.cs,csv2.cs)    
end

"""
checks if  2 cohorts of cohorts_t1 type equal
"""
function eq_csmtx(csm1,csm2)
    typeof(csm1)== typeof(csm2)&&
    csm1.year_b==csm2.year_b &&
    csm1.tmax==csm2.tmax &&
    csm1.jmax==csm2.jmax &&
    csm1.x==csm2.x
    csm1.y==csm2.y
    csm1.w==csm2.w
       
end
