csgelrub = function(nu){
  # CSGELRUB Gelman and Rubin convergence diagnostic.
  
  # RHAT = CSGELRUB(NU) Returns the Gelman and Rubin convergence diagnostic for
  # Markov chain Monte Carlo. The input to the function is a matrix of scalar
  # summaries NU. Each row of the matrix is a sequence of scalar summaries from the chain
  # For example, the scalar summaries might be the mean, median, etc.
  
  # W. L. and A. R. Martinez, 9/15/01
  
  k = dim(nu)[1]
  n = dim(nu)[2]
  # First find the Between Sequence variance
  museq = apply(nu,1,mean)
  mutot = mean(nu)
  # square each element and sum
  B = n/(k-1)*sum((museq-mutot)^2)
  # The logic here is that we see each point in a sequence as the mean of
  # the corresponding sequence. Hence, we take n terms here.
  # Since there are k sequences, the square term should be divided by
  # k-1 to get an unbiased estimator.
  
  # Find the within-sequence variance
  varseq = apply(nu,1,var) # variance of each row
  W = mean(varseq)
  
  varhat = (n-1)/n*W + B/n
  rhat = varhat/W
  return(rhat)
}