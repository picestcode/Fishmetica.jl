
Types and their transformations


<a id='Representing-cohorts-with-matrices-1'></a>

## Representing cohorts with matrices

<a id='Fishmetica.csmtx' href='#Fishmetica.csmtx'>#</a>
**`Fishmetica.csmtx`** &mdash; *Type*.



cohorts related matrix year_b - initial physical year tmax - number of seasons jmax - number of ages mtx: x,y,w,m,f, etc.


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L117-L123' class='documenter-source'>source</a><br>

<a id='Fishmetica.cohorts_t' href='#Fishmetica.cohorts_t'>#</a>
**`Fishmetica.cohorts_t`** &mdash; *Type*.



cohorts as 3 matrices: year_b - initial year, year_e - last year,  x - abundance, y - catch, w - weight (obsolete)


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L1-L3' class='documenter-source'>source</a><br>

<a id='Fishmetica.cohorts_t1' href='#Fishmetica.cohorts_t1'>#</a>
**`Fishmetica.cohorts_t1`** &mdash; *Type*.



cohorts as 3 tables (=cohorts_m) yb - physical year associated with cohort (when age=1) tmax - number of seasons jmax - number of ages x,y,w- tables with abundunce, catch, weight


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L86-L92' class='documenter-source'>source</a><br>


<a id='Representing-cohorts-with-vectors-1'></a>

## Representing cohorts with vectors

<a id='Fishmetica.cvec' href='#Fishmetica.cvec'>#</a>
**`Fishmetica.cvec`** &mdash; *Type*.



detail mtx representation of a cohort as sequence of vecs no- the number yb - physical year associated with cohort (when age=1) jbeg - youngest age  len - length jtrust - trusted age vec


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L131-L139' class='documenter-source'>source</a><br>

<a id='Fishmetica.csvecs' href='#Fishmetica.csvecs'>#</a>
**`Fishmetica.csvecs`** &mdash; *Type*.



cohort as vecs


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L150-L152' class='documenter-source'>source</a><br>

<a id='Fishmetica.cohort_v' href='#Fishmetica.cohort_v'>#</a>
**`Fishmetica.cohort_v`** &mdash; *Type*.



detail cohort representation  tmax - number of seasons jmax - number of ages yb - physical year assciated with cohort (when age=1) jb - youngest age  l - length jt - trusted age xv,yv,wv - vectors with abundunce, catch, weight


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L13-L22' class='documenter-source'>source</a><br>

<a id='Fishmetica.cohorts_v' href='#Fishmetica.cohorts_v'>#</a>
**`Fishmetica.cohorts_v`** &mdash; *Type*.



detail  representation of given cohorts tmax - number of seasons jmax - number of ages cs -  vector of cohort of type cohort_v


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L35-L40' class='documenter-source'>source</a><br>

<a id='Fishmetica.cohorts_v1' href='#Fishmetica.cohorts_v1'>#</a>
**`Fishmetica.cohorts_v1`** &mdash; *Type*.



detail  representation of given cohorts tmax - number of seasons jmax - number of ages cs -  vector of cohort of type cohort_v


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L48-L53' class='documenter-source'>source</a><br>

<a id='Fishmetica.cohort_v1' href='#Fishmetica.cohort_v1'>#</a>
**`Fishmetica.cohort_v1`** &mdash; *Type*.



detail cohort representation  no- the number yb - physical year associated with cohort (when age=1) tmax - number of seasons jmax - number of ages jbeg - initial age  len - length jtrust - trusted age xv,yv,wv - vectors with abundunce, catch, weight


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L62-L72' class='documenter-source'>source</a><br>

<a id='Fishmetica.cohorts_v2' href='#Fishmetica.cohorts_v2'>#</a>
**`Fishmetica.cohorts_v2`** &mdash; *Type*.



cohorts as cs: cohort_v1? year_b - physical year associated with cohort (when age=1) tmax - number of seasons jmax - number of ages


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L105-L110' class='documenter-source'>source</a><br>


##Generating vectors

<a id='Fishmetica.gvecs' href='#Fishmetica.gvecs'>#</a>
**`Fishmetica.gvecs`** &mdash; *Function*.



generates xv1,yv1 of csvecs type with use of xv::csvecs,yv::csvecs,wv::csvecs,zv::csvecs,gv::csvecs first makes a deep copy, then modifies vectors 


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L639-L642' class='documenter-source'>source</a><br>

<a id='Fishmetica.gvecs1' href='#Fishmetica.gvecs1'>#</a>
**`Fishmetica.gvecs1`** &mdash; *Function*.



generates xv1,yv1 of csvecs type with use of xv::csvecs,yv::csvecs,wv::csvecs,zv::csvecs,gv::csvecs first makes a deep copy, then modifies vectors 


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L844-L847' class='documenter-source'>source</a><br>


##Type transformations

<a id='Fishmetica.cst2csv2' href='#Fishmetica.cst2csv2'>#</a>
**`Fishmetica.cst2csv2`** &mdash; *Function*.



cohorts_t1->cohorts_v2

cohorts_t1: year_b - initial physical year tmax - number of seasons jmax - number of ages x,y,w -  matrices  with abundunce, catch, weight

cohorts_v2: year_b - initial physical year tmax - number of seasons jmax - number of ages cs - cohorts of type cohort_v1

cohort_v1: no- the number yb - physical year associated with cohort (when age=1) tmax - number of seasons jmax - number of ages jbeg - initial age  len - length jtrust - trusted age xv,yv,wv- vectors with abundunce, catch, weight tmax+jmax-1 – number of related cohorts


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L162-L188' class='documenter-source'>source</a><br>

<a id='Fishmetica.csv2cst' href='#Fishmetica.csv2cst'>#</a>
**`Fishmetica.csv2cst`** &mdash; *Function*.



cohorts_v2->cohirts_t1

cohorts_t1: year_b - initial physical year tmax - number of seasons jmax - number of ages x,y,w -  matrices  with abundunce, catch, weight

cohorts_v2: year_b - initial physical year tmax - number of seasons jmax - number of ages cs - cohorts of type cohort_v1

cohort_v1: no- the number yb - physical year associated with cohort (when age=1) tmax - number of seasons jmax - number of ages jbeg - initial age  len - length jtrust - trusted age xv,yv,wv- vectors with abundunce, catch, weight


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L253-L278' class='documenter-source'>source</a><br>

<a id='Fishmetica.csmtx2csvs' href='#Fishmetica.csmtx2csvs'>#</a>
**`Fishmetica.csmtx2csvs`** &mdash; *Function*.



cohorts related matrix (of type csmtx: x,y,w,m,f,z,g)->sequence of vectors (of type csvecs)

csrm:cohorts related matrix (of type csrm) year_b - initial physical year tmax - number of seasons jmax - number of ages x,y,w, m,f,g -  cohorts related matrices  

cvec: detail cohort representation as vec no- the number in the sequence yb - physical year associated with cohort (when age=1) jbeg - initial age  len - length jtrust - trusted age vec

csvecs: cohorts related mtx as a sequence of cvec no- the number yb - physical year associated with cohort (when age=1) tmax - number of seasons jmax - number of ages jbeg - initial age  len - length jtrust - trusted age xv,yv,wv- vectors with abundunce, catch, weight tmax+jmax-1 – number of related cohorts


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L401-L430' class='documenter-source'>source</a><br>

<a id='Fishmetica.csvs2csmtx' href='#Fishmetica.csvs2csmtx'>#</a>
**`Fishmetica.csvs2csmtx`** &mdash; *Function*.



sequence of vectors (of type csvecs)->cohorts related matrix (of type csrm) (x,y,w,m,f,z,g)

csrm:cohorts related matrix (of type csrm) year_b - initial physical year tmax - number of seasons jmax - number of ages x,y,w, m,f,g -  cohorts related matrices  

cvec: detail cohort representation as vec no- the number in the sequence yb - physical year associated with cohort (when age=1) jbeg - initial age  len - length jtrust - trusted age vec

csvecs: cohorts related mtx as a sequence of cvec no- the number yb - physical year associated with cohort (when age=1) tmax - number of seasons jmax - number of ages jbeg - initial age  len - length jtrust - trusted age xv,yv,wv- vectors with abundunce, catch, weight tmax+jmax-2 – number of related cohorts


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/fishmetica3.jl#L490-L519' class='documenter-source'>source</a><br>


<a id='Testing-type-transformations-1'></a>

## Testing type transformations

<a id='Fishmetica.eq_cvec' href='#Fishmetica.eq_cvec'>#</a>
**`Fishmetica.eq_cvec`** &mdash; *Function*.



checks if  2 cohorts of cvec type equal


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/tests.jl#L26-L28' class='documenter-source'>source</a><br>

<a id='Fishmetica.eq_csvecs' href='#Fishmetica.eq_csvecs'>#</a>
**`Fishmetica.eq_csvecs`** &mdash; *Function*.



checks if  2 descriptions of csmtx type equal


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/tests.jl#L40-L42' class='documenter-source'>source</a><br>

<a id='Fishmetica.eq_cst' href='#Fishmetica.eq_cst'>#</a>
**`Fishmetica.eq_cst`** &mdash; *Function*.



checks if  2 cohorts of cohorts_t type equal


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/tests.jl#L100-L102' class='documenter-source'>source</a><br>

<a id='Fishmetica.eq_cv1' href='#Fishmetica.eq_cv1'>#</a>
**`Fishmetica.eq_cv1`** &mdash; *Function*.



checks if  2 cohorts of cohort_v1 type equal


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/tests.jl#L113-L115' class='documenter-source'>source</a><br>

<a id='Fishmetica.eq_cv1_alt' href='#Fishmetica.eq_cv1_alt'>#</a>
**`Fishmetica.eq_cv1_alt`** &mdash; *Function*.



checks if  vectors of cohorts of cohort_v1 type equal alt


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/tests.jl#L130-L132' class='documenter-source'>source</a><br>

<a id='Fishmetica.eq_csv2' href='#Fishmetica.eq_csv2'>#</a>
**`Fishmetica.eq_csv2`** &mdash; *Function*.



checks if  2 cohorts of cohorts_v2 type equal


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/tests.jl#L142-L144' class='documenter-source'>source</a><br>

<a id='Fishmetica.eq_csmtx' href='#Fishmetica.eq_csmtx'>#</a>
**`Fishmetica.eq_csmtx`** &mdash; *Function*.



checks if  2 cohorts of cohorts_t1 type equal


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/tests.jl#L153-L155' class='documenter-source'>source</a><br>

<a id='Fishmetica.runtest_cst2csv' href='#Fishmetica.runtest_cst2csv'>#</a>
**`Fishmetica.runtest_cst2csv`** &mdash; *Function*.



checks if direct*reverse=equal for cst2csv2 and csv2cst


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/tests.jl#L1-L3' class='documenter-source'>source</a><br>

<a id='Fishmetica.runtest_csmtx2csvs' href='#Fishmetica.runtest_csmtx2csvs'>#</a>
**`Fishmetica.runtest_csmtx2csvs`** &mdash; *Function*.



checks if csv2cst(cst2csv2(cst)=cst


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/tests.jl#L54-L56' class='documenter-source'>source</a><br>

