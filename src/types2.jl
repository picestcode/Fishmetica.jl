"""
cohorts as 3 matrices: year_b - the initial year,j0 - youngest age, j1 - oldest age, jmax - maximal age,  x - abundance, y - catch, w - weight (obsolete)
"""
type cohort1
    year_b
    j0
    j1
    jmax
     x::Array{Float64,2}
     y::Array{Float64,2}
     w::Array{Float64,2}
 end

"""
cohorts as 3 matrices: year_b - the initial year,   x - abundance, y - catch, w - weight (standard)
"""

type cohort
     year_b
    
     x::Array{Float64,2}
     y::Array{Float64,2}
     w::Array{Float64,2}
 end

"""
list of cohorts (obsolete)
"""
type  cohorts  
     year_b
     tmax
     jmax
     cs::Array{cohort,1}          
end

"""
cohorts as 3 matrices: year_b - the initial year, tmax - number of seasons, jmax - maximal age,  x - abundance, y - catch, w - weight (obsolete)
"""

type  cohorts1  
     year_b
     tmax
     jmax
     cs::Array{cohort1,1}          
end

"""
cohorts as 3 matrices: year_b - the initial year, x - abundance, y - catch, w - weight (standard)
"""
type  cohorts_m  
     year_b
    x
     y
    w
              
end

"""
cohorts as 3 matrices: year_b - initial year, year_e - last year,  x - abundance, y - catch, w - weight (obsolete)
"""
type  cohorts_t  
     year_b
    year_e
    x
     y
    w
              
end

"""
cohort as model parameters: year_b - initial year, j_max - maximum age,  x1j - first row of x,  w - weight, k_x, k_w - coefficients for x and y, p_m - natural mortality parameters, h - H matrix,  p_s - - fishing mortality parameters,  ft
"""
type cohort_mp
     year_b
     j_max
     x1j::Array{Float64,2}# row
     w::Array{Float64,2}     
     k_x::Float64
     k_w::Float64    
     p_m       
     h::Array{Float64,2}   
     p_s    
    ft::Array{Float64,1}# column   
end


"""
transforms cohort from cohort_mp to  cohort_m
"""
function cmp2cm(cmp)
   
    jmax=cmp.j_max
     x1j=cmp.x1j
          kx=cmp.k_x
     kw=cmp.k_w    
     a1m1,a2m1,b1m1,b2m1,cm1,d1m1,d2m1=cmp.p_m       
     Hb=cmp.h   
     as1,bs1=cmp.p_s    
    ftb=cmp.ft 
    
    tmax=length(ftb)
    
    x,y=gm_burago(tmax,jmax,Mb,Fb,W,kx,kw,x1j)#
    
    return cohorts_m(cmp.year_b,x,y,cmp.w)
    
end

# from cohort_m to cohort_mp !