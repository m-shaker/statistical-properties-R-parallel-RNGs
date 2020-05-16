
## Statistical tests:

#L'Ecuyer Combined multiple-recursive generator

#Loading the necessary package
library(parallel)   
library(randtoolbox)

#Create a cluster with 4 cores
cl <- makeCluster(4)

#Setting the right parallel RNG
iseed <- 0
clusterSetRNGStream(cl = cl, iseed = iseed)

#Confirming that our generator is the L'Ecuyer Combined multiple-recursive generator
parLapply(cl, 1:4, function(n)RNGkind())

#Each stream will have 100,000 random numbers
nSims <- 100000
taskFun <- function(i){
  val <- runif(1)
  return(val)
}

#Defining matrix to hold generated sequences 
cmrg_sequences = matrix(nrow=1000, ncol=100000)

#Generate 1000 random sequences in parallel each of length 100,000 using the L'Ecuyer Combined #multiple-recursive generator
for (i in 1:1000){
  cmrg_sequences[i,] <- parSapply(cl, 1:nSims, taskFun)
}

#Stopping the cluster
stopCluster(cl)


#Now we will use the 64-Bit Linear congruential generator:
  
#64-Bit Linear congruential generator

#Loading the rTRNG library 
library(rTRNG)

#Indicating that processes should be done in parallel with 4 threads
RcppParallel::setThreadOptions(numThreads = 4L)

#Defining matrix to hold generated sequences 
lcg64_sequences = matrix(nrow=1000, ncol=100000)

#Setting up the right generator 
TRNGkind("lcg64")

#Confirming that the right generator is in use
TRNGkind()

#Generate 1000 random sequences in parallel each of length 100,000 using the 64-Bit Linear congruential #generator
for (i in 1:1000){
  lcg64_sequences[i,] <- runif_trng(100000, min = 0, max = 1, parallelGrain = 100L)
}


### Gap test R implementation using sequences generated from the L'Ecuyer Combined 
###multiple-recursive generator
r1 = rep(0,1000)
for (i in 1:1000){
  if (gap.test(cmrg_sequences[i,])$p.value < 0.05){
    r1[i] =1
  } 
}
cat("Rejection:", sum(r1)/1000)


### Frequency test R implementation using sequences generated from the L'Ecuyer 
### Combined multiple-recursive generator
r2 = rep(0,1000)
for (i in 1:1000){
  if (freq.test(cmrg_sequences[i,],1:10)$p.value < 0.05){
    r2[i] =1
  } 
}
cat("Rejection:", sum(r2)/1000)


### Serial test R implementation using sequences generated from the L'Ecuyer 
### Combined multiple-recursive generator
r3 = rep(0,1000)
for (i in 1:1000){
  if (serial.test(cmrg_sequences[i,])$p.value < 0.05){
    r3[i] =1
  } 
}
cat("Rejection:", sum(r3)/1000)


### Poker test R implementation using sequences generated from the L'Ecuyer 
### Combined multiple-recursive generator
r4 = rep(0,1000)
for (i in 1:1000){
  if (poker.test(cmrg_sequences[i,])$p.value < 0.05){
    r4[i] =1
  } 
}
cat("Rejection:", sum(r4)/1000)


### Order test R implementation using sequences generated from the L'Ecuyer 
### Combined multiple-recursive generator
r5 = rep(0,1000)
for (i in 1:1000){
  if (order.test(cmrg_sequences[i,],d=5)$p.value < 0.05){
    r5[i] =1
  } 
}
cat("Rejection:", sum(r5)/1000)


### Gap test R implementation using sequences generated from the 64-Bit Linear 
### congruential generator
l1 = rep(0,1000)
for (i in 1:1000){
  if (gap.test(lcg64_sequences[i,])$p.value < 0.05){
    l1[i] =1
  } 
}
cat("Rejection:", sum(l1)/1000)


### Frequency test R implementation using sequences generated from the 64-Bit 
### Linear congruential generator
l2 = rep(0,1000)
for (i in 1:1000){
  if (freq.test(lcg64_sequences[i,])$p.value < 0.05){
    l2[i] =1
  } 
}
cat("Rejection:", sum(l2)/1000)


### Serial test R implementation using sequences generated from the 64-Bit 
### Linear congruential generator
l3 = rep(0,1000)
for (i in 1:1000){
  if (serial.test(lcg64_sequences[i,])$p.value < 0.05){
    l3[i] =1
  } 
}
cat("Rejection:", sum(l3)/1000)


### Poker test R implementation using sequences generated from the 64-Bit Linear 
### congruential generator
l4 = rep(0,1000)
for (i in 1:1000){
  if (poker.test(lcg64_sequences[i,])$p.value < 0.05){
    l4[i] =1
  } 
}
cat("Rejection:", sum(l4)/1000)


### Order test R implementation using sequences generated from the 64-Bit Linear 
### congruential generator
l5 = rep(0,1000)
for (i in 1:1000){
  if (order.test(lcg64_sequences[i,],d=5)$p.value < 0.05){
    l5[i] =1
  } 
}
cat("Rejection:", sum(l5)/1000)
