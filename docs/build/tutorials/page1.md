
Fuctions to generate $\mu$'s and $\sigma^2$'s


<a id='Filtering-1'></a>

## Filtering

<a id='Fishmetica.sigma2f' href='#Fishmetica.sigma2f'>#</a>
**`Fishmetica.sigma2f`** &mdash; *Function*.



filtering: returns a vector of deviations ($sigma^2^f$) for a cohort  sigma^2_0 - variance  for initial age sigma^2_e -  variance  for the evidance vector (catches) sigma^2_x -  variance  for the state vector (population) l - the length of the cohort


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/hmm.jl#L28' class='documenter-source'>source</a><br>

<a id='Fishmetica.muf' href='#Fishmetica.muf'>#</a>
**`Fishmetica.muf`** &mdash; *Function*.



filtering: returns a vector of means (mu^f) for a cohort  alpha - coefficients for the state x_{t+1}=x_t +alpha_t beta - coefficients for the observation y_{t}=x_t +beta_t  e - the evidance vector mu_0 - mean for the initial age sigma^2_0 - variance (deviation is the square root of its variance) for initial age mu_e - mean for the evidence sigma^2_e - variance  for the evidance  (catch) mu_x - mean for the state sigma^2_e - variance for the state (population) l - the length of the cohort


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/hmm.jl#L46-L58' class='documenter-source'>source</a><br>

<a id='Fishmetica.genmus2f' href='#Fishmetica.genmus2f'>#</a>
**`Fishmetica.genmus2f`** &mdash; *Function*.



filtering: returns filtered cohorts represented as vecs: csvlogts1muf - means,csvlogts1s2f - sigmas^2 (deviations),


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/hmm.jl#L75-L77' class='documenter-source'>source</a><br>


<a id='Prediction-1'></a>

## Prediction

<a id='Fishmetica.sigma2p' href='#Fishmetica.sigma2p'>#</a>
**`Fishmetica.sigma2p`** &mdash; *Function*.



prediction: returns the vector of deviation predictions (sigma^2^p_{t+1},..., sigma^2^p_{t+k+1}) for a cohort  sigma^2_x - variance for the state vector (population) l - the length of the cohort


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/hmm.jl#L125-L129' class='documenter-source'>source</a><br>

<a id='Fishmetica.mup' href='#Fishmetica.mup'>#</a>
**`Fishmetica.mup`** &mdash; *Function*.



prediction: returns the vector of mean predictions (mu^p_{t+1}, ..., mu^p_{t+k+1}) for a cohort  alphak=(alpha^t, ..., alpha_{t+k}) - (virtual) coefficients for x_{t+k+1}=x_{t+k} +alpha_{t+k} -  muft - mu^f_t mu_x - mean for the state k >= 1 - the prediction step


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/hmm.jl#L139-L145' class='documenter-source'>source</a><br>

<a id='Fishmetica.genmus2p' href='#Fishmetica.genmus2p'>#</a>
**`Fishmetica.genmus2p`** &mdash; *Function*.



makes k predictions for all cohorts (a^e+1, .., a^e+k)  with yr>=yr2. cuts extra elements with age greater than jmax.  returns filtered cohorts represented as vecs: csvlogts1mup - means,csvlogts1s2p - sigmas^2 (deviations), mup and s2p (means and sigmas^2 (deviations))


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/hmm.jl#L155-L159' class='documenter-source'>source</a><br>


<a id='Smoothing-1'></a>

## Smoothing

<a id='Fishmetica.sigma2b' href='#Fishmetica.sigma2b'>#</a>
**`Fishmetica.sigma2b`** &mdash; *Function*.



backward: returns a vector of deviations (sigma^2^b) for a cohort  sigma^2_0 - variance  for initial age sigma^2_e -  variance  for the evidance vector (catches) sigma^2_x -  variance  for the state vector (population) l - the length of the cohort


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/hmm.jl#L275-L281' class='documenter-source'>source</a><br>

<a id='Fishmetica.genmus2s' href='#Fishmetica.genmus2s'>#</a>
**`Fishmetica.genmus2s`** &mdash; *Function*.



smoothing: makes several cohorts represented as vecs: csvlogts1mus - means,csvlogts1s2s - sigmas^2 (deviations)


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/hmm.jl#L296-L298' class='documenter-source'>source</a><br>


<a id='Auxilary-1'></a>

## Auxilary

<a id='Fishmetica.genG' href='#Fishmetica.genG'>#</a>
**`Fishmetica.genG`** &mdash; *Function*.



generates G with M (natural mortality) and F (catch mortality), Z (full mortality)


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/hmm.jl#L11-L13' class='documenter-source'>source</a><br>

<a id='Fishmetica.genZ' href='#Fishmetica.genZ'>#</a>
**`Fishmetica.genZ`** &mdash; *Function*.



generates Z (full)  with M (natural mortality) and F (catch mortality)


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/hmm.jl#L1-L3' class='documenter-source'>source</a><br>

<a id='Fishmetica.genQuantile' href='#Fishmetica.genQuantile'>#</a>
**`Fishmetica.genQuantile`** &mdash; *Function*.



at +p%


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/hmm.jl#L397-L399' class='documenter-source'>source</a><br>

<a id='Fishmetica.genQuantile1' href='#Fishmetica.genQuantile1'>#</a>
**`Fishmetica.genQuantile1`** &mdash; *Function*.



at -p%


<a target='_blank' href='https://github.com/picestcode/Fishmetica.jl/blob/028737beb977c8c74424b038d5cf0c6bb61155ec/src/hmm.jl#L414-L416' class='documenter-source'>source</a><br>

