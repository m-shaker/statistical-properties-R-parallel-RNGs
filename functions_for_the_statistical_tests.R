
## Defining functions for the statistical tests

#Gap test

gap.test <- function(u,lower =0,upper=1/2,echo=TRUE)
{
  #Compute gaps such as (g_i)=1 if (u_i) > max or (u_i)<min, 0 otherwise 
  gap <-   (!((u<=upper)&(u>=lower)))*1 
  n <- length(u)
  p <- upper - lower
  
  #Find index of zeros 
  indexzero <- (gap == 1)*1:n
  indexzero <- indexzero[indexzero != 0]
  indexzero <- c(0, indexzero, n+1)
  lindzero <- length(indexzero)
  
  #Compute sizes of zero lengths
  lengthsize <- indexzero[2:lindzero] - indexzero[2:lindzero-1] -1
  lengthsize <- lengthsize[lengthsize != 0]
  maxlen <- max(lengthsize)
  maxlen <- max(maxlen, floor((log(10^(-1))- 2*log(1-p) - log(n)) / log(p)))
  
  #Compute observed and theoretical frequencies
  obsnum <- sapply(1:maxlen, function(t) sum(lengthsize == t))
  expnum <- (1-p)^2*p^(1:maxlen)*n
  
  #Compute chisquare statistic
  residu <- (obsnum - expnum)/sqrt(expnum)
  stat <- sum(residu^2)
  pvalue <- pchisq( stat, maxlen - 1, lower.tail = FALSE)
  options(digits=2)    
  if( echo )
  {
    cat("\n\t\t\t Gap test\n")
    cat("\nchisq stat = ", stat, ", df = ",maxlen-1, ", p-value = ", pvalue, "\n", sep="")
    cat("\n\t\t (sample size : ",length(u),")\n\n", sep="")
    cat("length\tobserved freq\t\ttheoretical freq\n")
    for(i in 1:maxlen)
      cat(i,"\t\t\t", obsnum[i],"\t\t\t", expnum[i],"\n")
  }
  
  res <- list(statistic=stat, parameter=maxlen-1, 
              p.value=pvalue, observed=obsnum, 
              expected=expnum, residuals=residu) 
  return(invisible(res))
}


#Frequency test
freq.test <- function(u, seq=0:15, echo=TRUE)
{
  
  #Compute integers such as min(seq) <= integernum[i] < max(seq)
  integernum <- floor(u*length(seq) + min(seq))
  
  #Observed numbers equal to seq[i]
  obsnum <- sapply(seq, function(x) sum(integernum == x))
  
  #Expected number equal to seq[i]
  expnum <- length(u)/length(seq)    
  
  #Compute chisquare statistic
  residu <-  (obsnum - expnum)/sqrt(expnum)
  stat <- sum(residu^2)
  pvalue <- pchisq(stat, length(seq) - 1, lower.tail = FALSE)
  options(digits=2)    
  if(echo)
  {
    cat("\n\t\t\t Frequency test\n")
    cat("\nchisq stat = ", stat, ", df = ",length(seq)-1, ", p-value = ", pvalue, "\n", sep="")
    cat("\n\t\t (sample size : ",length(u),")\n\n", sep="")
    cat("\tobserved number\t",obsnum,"\n")
    cat("\texpected number\t",expnum,"\n")    
  } 
  
  res <- list(statistic = stat, parameter=length(seq)-1, 
              p.value = pvalue, observed = obsnum, 
              expected = expnum, residuals =residu ) 
  return( invisible( res ) )
}



#Serial test
serial.test <- function(u , d=8, echo=TRUE)
{
  if(length(u)/d-as.integer(length(u)/d) > 0)
    stop("the length of 'u' must be a multiple of d")
  
  #Compute pairs in {0, ..., d-1}
  pair <- matrix(floor(u*d), length(u)/2, 2)
  
  #Compute (u_i)*d+(u_i)+1 in {0, ..., d^2-1}
  poly <- pair[,1]*d+pair[,2]
  
  #Compute numbers
  obsnum <- sapply(0:(d^2 -1), function(x) sum(poly == x))
  expnum <- length(u)/(2*d^2)    
  
  #Compute chisquare statistic
  residu <- (obsnum - expnum)/sqrt(expnum) 
  stat <- sum(residu^2)
  pvalue <- pchisq(stat, d^2-1, lower.tail=FALSE)
  options(digits=2)    
  if( echo )
  {
    cat("\n\t\t\t Serial test\n")
    cat("\nchisq stat = ", stat, ", df = ",d^2-1, ", p-value = ", pvalue, "\n", sep="")
    cat("\n\t\t (sample size : ",length(u),")\n\n", sep="")
    cat("\tobserved number\t",obsnum,"\n")
    cat("\texpected number\t",expnum,"\n")    
  } 
  
  res <- list(statistic = stat, parameter=d^2-1, 
              p.value=pvalue, observed=obsnum, 
              expected=expnum, residuals=residu) 
  return(invisible(res))
}



#Poker.test 
poker.test <- function(u, nbcard=5,echo=TRUE)
{
  if(length(u)/nbcard-as.integer(length(u)/nbcard) > 0)
    stop("the length of 'u' must be a multiple of nbcard")
  
  #Number of lines  
  nbl <- length(u)/nbcard
  
  #Compute "nbcard-hands" in {0, ..., k-1}
  hands <- matrix(as.integer(floor(u*nbcard)),nbl, nbcard) 
  
  #Compute observed hands
  obshands <- .Call("doPokerTest", hands, nbl, nbcard)
  
  #Compute expected hands
  fact <- vector("numeric", nbcard+1)
  fact[1] <- 1
  for(i in 2:(nbcard+1))
    fact[i] <- fact[i-1]*(i-1)
  
  #Fact contains 0!,1!, 2! ... (2*nbcard)!
  ind <- 1:nbcard
  stirlingNum <- stirling(nbcard)[-1]
  exphands <- 1/nbcard^nbcard*fact[nbcard+1]/fact[nbcard-ind+1]*stirlingNum[ind]*nbl
  
  stat <- sum((obshands-exphands)^2/exphands)
  pvalue <- pchisq(stat, nbcard-1, lower.tail=FALSE)
  
  options(digits=2)
  if( echo )
  {
    cat("\n\t\t\t Poker test\n")
    cat("\nchisq stat = ", stat, ", df = ",nbcard-1, ", p-value = ", pvalue, "\n", sep="")
    cat("\n\t\t (sample size : ",length(u),")\n\n", sep="")
    cat("\tobserved number\t",obshands,"\n")
    cat("\texpected number\t",exphands,"\n")    
  } 
  
  res <- list(statistic=stat, parameter=nbcard-1, 
              p.value=pvalue, observed=obshands, 
              expected=exphands, residuals=(obshands-exphands)^2/exphands) 
  return(invisible(res))
}



#Order test 
order.test <- function(u, d=3,echo=TRUE)
{
  if(d > 8)
    stop("too long to compute this exponential cost problem.")
  if(d < 2)
    stop("wrong argument 'd'.")
  
  if(length(u)/d-as.integer(length(u)/d) > 0)
    stop("the length of 'u' must be a multiple of d")
  
  #Store u in a matrix
  if(is.vector(u))
    u <- matrix(u, length(u) / d, d)
  if(!is.matrix(u))
    stop("wrong argument u")
  
  #Compute observed numbers manually to save on computational time 
  factOfD <- factorial(d)
  obsnum <- vector("numeric", length=factOfD)
  if(d == 2)
  {
    obsnum[1] <- sum(u[,1] < u[,2])
    obsnum[2] <- sum(u[,2] < u[,1])
  }
  if(d == 3)
  {
    obsnum[1] <- sum(u[,1] < u[,2] & u[,2] <u[,3])
    obsnum[2] <- sum(u[,1] < u[,3] & u[,3] <u[,2])
    obsnum[3] <- sum(u[,2] < u[,1] & u[,1] <u[,3])
    obsnum[4] <- sum(u[,2] < u[,3] & u[,3] <u[,1])
    obsnum[5] <- sum(u[,3] < u[,2] & u[,2] <u[,1])
    obsnum[6] <- sum(u[,3] < u[,1] & u[,1] <u[,2])
  }
  if(d == 4)
  {
    obsnum[1] <- sum(u[,1] < u[,2] & u[,2] <u[,3] & u[,3] <u[,4])
    obsnum[2] <- sum(u[,1] < u[,3] & u[,3] <u[,2] & u[,2] <u[,4])
    obsnum[3] <- sum(u[,2] < u[,1] & u[,1] <u[,3] & u[,3] <u[,4])
    obsnum[4] <- sum(u[,2] < u[,3] & u[,3] <u[,1] & u[,1] <u[,4])
    obsnum[5] <- sum(u[,3] < u[,2] & u[,2] <u[,1] & u[,1] <u[,4])
    obsnum[6] <- sum(u[,3] < u[,1] & u[,1] <u[,2] & u[,2] <u[,4])
    obsnum[7] <- sum(u[,1] < u[,2] & u[,2] <u[,4] & u[,4] <u[,3])
    obsnum[8] <- sum(u[,1] < u[,3] & u[,3] <u[,4] & u[,4] <u[,2])
    obsnum[9] <- sum(u[,2] < u[,1] & u[,1] <u[,4] & u[,4] <u[,3])
    obsnum[10] <- sum(u[,2] < u[,3] & u[,3] <u[,4] & u[,4] <u[,1])
    obsnum[11] <- sum(u[,3] < u[,2] & u[,2] <u[,4] & u[,4] <u[,1])
    obsnum[12] <- sum(u[,3] < u[,1] & u[,1] <u[,4] & u[,4] <u[,2])
    obsnum[13] <- sum(u[,1] < u[,4] & u[,4] <u[,2] & u[,2] <u[,3])
    obsnum[14] <- sum(u[,1] < u[,4] & u[,4] <u[,3] & u[,3] <u[,2])
    obsnum[15] <- sum(u[,2] < u[,4] & u[,4] <u[,1] & u[,1] <u[,3])
    obsnum[16] <- sum(u[,2] < u[,4] & u[,4] <u[,3] & u[,3] <u[,1])
    obsnum[17] <- sum(u[,3] < u[,4] & u[,4] <u[,2] & u[,2] <u[,1])
    obsnum[18] <- sum(u[,3] < u[,4] & u[,4] <u[,1] & u[,1] <u[,2])
    obsnum[19] <- sum(u[,4] < u[,1] & u[,1] <u[,2] & u[,2] <u[,3])
    obsnum[20] <- sum(u[,4] < u[,1] & u[,1] <u[,3] & u[,3] <u[,2])
    obsnum[21] <- sum(u[,4] < u[,2] & u[,2] <u[,1] & u[,1] <u[,3])
    obsnum[22] <- sum(u[,4] < u[,2] & u[,2] <u[,3] & u[,3] <u[,1])
    obsnum[23] <- sum(u[,4] < u[,3] & u[,3] <u[,2] & u[,2] <u[,1])
    obsnum[24] <- sum(u[,4] < u[,3] & u[,3] <u[,1] & u[,1] <u[,2])
  }
  
  if(d == 5)
  {
    obsnum[1]<-sum(u[,1]<u[,2]&u[,2]<u[,3]&u[,3]<u[,4]&u[,4] <u[,5])
    obsnum[2]<-sum(u[,1]<u[,3]&u[,3]<u[,2]&u[,2]<u[,4]&u[,4] <u[,5])
    obsnum[3]<-sum(u[,2]<u[,1]&u[,1]<u[,3]&u[,3]<u[,4]&u[,4] <u[,5])
    obsnum[4]<-sum(u[,2]<u[,3]&u[,3]<u[,1]&u[,1]<u[,4]&u[,4] <u[,5])
    obsnum[5]<-sum(u[,3]<u[,2]&u[,2]<u[,1]&u[,1]<u[,4]&u[,4] <u[,5])
    obsnum[6]<-sum(u[,3]<u[,1]&u[,1]<u[,2]&u[,2]<u[,4]&u[,4] <u[,5])
    obsnum[7]<-sum(u[,1]<u[,2]&u[,2]<u[,4]&u[,4]<u[,3]&u[,3] <u[,5])
    obsnum[8]<-sum(u[,1]<u[,3]&u[,3]<u[,4]&u[,4]<u[,2]&u[,2] <u[,5])
    obsnum[9]<-sum(u[,2]<u[,1]&u[,1]<u[,4]&u[,4]<u[,3]&u[,3] <u[,5])
    obsnum[10]<-sum(u[,2]<u[,3]&u[,3]<u[,4]&u[,4]<u[,1]&u[,1] <u[,5])
    obsnum[11]<-sum(u[,3]<u[,2]&u[,2]<u[,4]&u[,4]<u[,1]&u[,1] <u[,5])
    obsnum[12]<-sum(u[,3]<u[,1]&u[,1]<u[,4]&u[,4]<u[,2]&u[,2] <u[,5])
    obsnum[13]<-sum(u[,1]<u[,4]&u[,4]<u[,2]&u[,2]<u[,3]&u[,3] <u[,5])
    obsnum[14]<-sum(u[,1]<u[,4]&u[,4]<u[,3]&u[,3]<u[,2]&u[,2] <u[,5])
    obsnum[15]<-sum(u[,2]<u[,4]&u[,4]<u[,1]&u[,1]<u[,3]&u[,3] <u[,5])
    obsnum[16]<-sum(u[,2]<u[,4]&u[,4]<u[,3]&u[,3]<u[,1]&u[,1] <u[,5])
    obsnum[17]<-sum(u[,3]<u[,4]&u[,4]<u[,2]&u[,2]<u[,1]&u[,1] <u[,5])
    obsnum[18]<-sum(u[,3]<u[,4]&u[,4]<u[,1]&u[,1]<u[,2]&u[,2] <u[,5])
    obsnum[19]<-sum(u[,4]<u[,1]&u[,1]<u[,2]&u[,2]<u[,3]&u[,3] <u[,5])
    obsnum[20]<-sum(u[,4]<u[,1]&u[,1]<u[,3]&u[,3]<u[,2]&u[,2] <u[,5])
    obsnum[21]<-sum(u[,4]<u[,2]&u[,2]<u[,1]&u[,1]<u[,3]&u[,3] <u[,5])
    obsnum[22]<-sum(u[,4]<u[,2]&u[,2]<u[,3]&u[,3]<u[,1]&u[,1] <u[,5])
    obsnum[23]<-sum(u[,4]<u[,3]&u[,3]<u[,2]&u[,2]<u[,1]&u[,1] <u[,5])
    obsnum[24]<-sum(u[,4]<u[,3]&u[,3]<u[,1]&u[,1]<u[,2]&u[,2] <u[,5])
    obsnum[25]<-sum(u[,1]<u[,2]&u[,2]<u[,3]&u[,3]<u[,5]&u[,5] <u[,4])
    obsnum[26]<-sum(u[,1]<u[,3]&u[,3]<u[,2]&u[,2]<u[,5]&u[,5] <u[,4])
    obsnum[27]<-sum(u[,2]<u[,1]&u[,1]<u[,3]&u[,3]<u[,5]&u[,5] <u[,4])
    obsnum[28]<-sum(u[,2]<u[,3]&u[,3]<u[,1]&u[,1]<u[,5]&u[,5] <u[,4])
    obsnum[29]<-sum(u[,3]<u[,2]&u[,2]<u[,1]&u[,1]<u[,5]&u[,5] <u[,4])
    obsnum[30]<-sum(u[,3]<u[,1]&u[,1]<u[,2]&u[,2]<u[,5]&u[,5] <u[,4])
    obsnum[31]<-sum(u[,1]<u[,2]&u[,2]<u[,4]&u[,4]<u[,5]&u[,5] <u[,3])
    obsnum[32]<-sum(u[,1]<u[,3]&u[,3]<u[,4]&u[,4]<u[,5]&u[,5] <u[,2])
    obsnum[33]<-sum(u[,2]<u[,1]&u[,1]<u[,4]&u[,4]<u[,5]&u[,5] <u[,3])
    obsnum[34]<-sum(u[,2]<u[,3]&u[,3]<u[,4]&u[,4]<u[,5]&u[,5] <u[,1])
    obsnum[35]<-sum(u[,3]<u[,2]&u[,2]<u[,4]&u[,4]<u[,5]&u[,5] <u[,1])
    obsnum[36]<-sum(u[,3]<u[,1]&u[,1]<u[,4]&u[,4]<u[,5]&u[,5] <u[,2])
    obsnum[37]<-sum(u[,1]<u[,4]&u[,4]<u[,2]&u[,2]<u[,5]&u[,5] <u[,3])
    obsnum[38]<-sum(u[,1]<u[,4]&u[,4]<u[,3]&u[,3]<u[,5]&u[,5] <u[,2])
    obsnum[39]<-sum(u[,2]<u[,4]&u[,4]<u[,1]&u[,1]<u[,5]&u[,5] <u[,3])
    obsnum[40]<-sum(u[,2]<u[,4]&u[,4]<u[,3]&u[,3]<u[,5]&u[,5] <u[,1])
    obsnum[41]<-sum(u[,3]<u[,4]&u[,4]<u[,2]&u[,2]<u[,5]&u[,5] <u[,1])
    obsnum[42]<-sum(u[,3]<u[,4]&u[,4]<u[,1]&u[,1]<u[,5]&u[,5] <u[,2])
    obsnum[43]<-sum(u[,4]<u[,1]&u[,1]<u[,2]&u[,2]<u[,5]&u[,5] <u[,3])
    obsnum[44]<-sum(u[,4]<u[,1]&u[,1]<u[,3]&u[,3]<u[,5]&u[,5] <u[,2])
    obsnum[45]<-sum(u[,4]<u[,2]&u[,2]<u[,1]&u[,1]<u[,5]&u[,5] <u[,3])
    obsnum[46]<-sum(u[,4]<u[,2]&u[,2]<u[,3]&u[,3]<u[,5]&u[,5] <u[,1])
    obsnum[47]<-sum(u[,4]<u[,3]&u[,3]<u[,2]&u[,2]<u[,5]&u[,5] <u[,1])
    obsnum[48]<-sum(u[,4]<u[,3]&u[,3]<u[,1]&u[,1]<u[,5]&u[,5] <u[,2])
    obsnum[49]<-sum(u[,1]<u[,2]&u[,2]<u[,5]&u[,5]<u[,3]&u[,3] <u[,4])
    obsnum[50]<-sum(u[,1]<u[,3]&u[,3]<u[,5]&u[,5]<u[,2]&u[,2] <u[,4])
    obsnum[51]<-sum(u[,2]<u[,1]&u[,1]<u[,5]&u[,5]<u[,3]&u[,3] <u[,4])
    obsnum[52]<-sum(u[,2]<u[,3]&u[,3]<u[,5]&u[,5]<u[,1]&u[,1] <u[,4])
    obsnum[53]<-sum(u[,3]<u[,2]&u[,2]<u[,5]&u[,5]<u[,1]&u[,1] <u[,4])
    obsnum[54]<-sum(u[,3]<u[,1]&u[,1]<u[,5]&u[,5]<u[,2]&u[,2] <u[,4])
    obsnum[55]<-sum(u[,1]<u[,2]&u[,2]<u[,5]&u[,5]<u[,4]&u[,4] <u[,3])
    obsnum[56]<-sum(u[,1]<u[,3]&u[,3]<u[,5]&u[,5]<u[,4]&u[,4] <u[,2])
    obsnum[57]<-sum(u[,2]<u[,1]&u[,1]<u[,5]&u[,5]<u[,4]&u[,4] <u[,3])
    obsnum[58]<-sum(u[,2]<u[,3]&u[,3]<u[,5]&u[,5]<u[,4]&u[,4] <u[,1])
    obsnum[59]<-sum(u[,3]<u[,2]&u[,2]<u[,5]&u[,5]<u[,4]&u[,4] <u[,1])
    obsnum[60]<-sum(u[,3]<u[,1]&u[,1]<u[,5]&u[,5]<u[,4]&u[,4] <u[,2])
    obsnum[61]<-sum(u[,1]<u[,4]&u[,4]<u[,5]&u[,5]<u[,2]&u[,2] <u[,3])
    obsnum[62]<-sum(u[,1]<u[,4]&u[,4]<u[,5]&u[,5]<u[,3]&u[,3] <u[,2])
    obsnum[63]<-sum(u[,2]<u[,4]&u[,4]<u[,5]&u[,5]<u[,1]&u[,1] <u[,3])
    obsnum[64]<-sum(u[,2]<u[,4]&u[,4]<u[,5]&u[,5]<u[,3]&u[,3] <u[,1])
    obsnum[65]<-sum(u[,3]<u[,4]&u[,4]<u[,5]&u[,5]<u[,2]&u[,2] <u[,1])
    obsnum[66]<-sum(u[,3]<u[,4]&u[,4]<u[,5]&u[,5]<u[,1]&u[,1] <u[,2])
    obsnum[67]<-sum(u[,4]<u[,1]&u[,1]<u[,5]&u[,5]<u[,2]&u[,2] <u[,3])
    obsnum[68]<-sum(u[,4]<u[,1]&u[,1]<u[,5]&u[,5]<u[,3]&u[,3] <u[,2])
    obsnum[69]<-sum(u[,4]<u[,2]&u[,2]<u[,5]&u[,5]<u[,1]&u[,1] <u[,3])
    obsnum[70]<-sum(u[,4]<u[,2]&u[,2]<u[,5]&u[,5]<u[,3]&u[,3] <u[,1])
    obsnum[71]<-sum(u[,4]<u[,3]&u[,3]<u[,5]&u[,5]<u[,2]&u[,2] <u[,1])
    obsnum[72]<-sum(u[,4]<u[,3]&u[,3]<u[,5]&u[,5]<u[,1]&u[,1] <u[,2])
    obsnum[73]<-sum(u[,1]<u[,5]&u[,5]<u[,2]&u[,2]<u[,3]&u[,3] <u[,4])
    obsnum[74]<-sum(u[,1]<u[,5]&u[,5]<u[,3]&u[,3]<u[,2]&u[,2] <u[,4])
    obsnum[75]<-sum(u[,2]<u[,5]&u[,5]<u[,1]&u[,1]<u[,3]&u[,3] <u[,4])
    obsnum[76]<-sum(u[,2]<u[,5]&u[,5]<u[,3]&u[,3]<u[,1]&u[,1] <u[,4])
    obsnum[77]<-sum(u[,3]<u[,5]&u[,5]<u[,2]&u[,2]<u[,1]&u[,1] <u[,4])
    obsnum[78]<-sum(u[,3]<u[,5]&u[,5]<u[,1]&u[,1]<u[,2]&u[,2] <u[,4])
    obsnum[79]<-sum(u[,1]<u[,5]&u[,5]<u[,2]&u[,2]<u[,4]&u[,4] <u[,3])
    obsnum[80]<-sum(u[,1]<u[,5]&u[,5]<u[,3]&u[,3]<u[,4]&u[,4] <u[,2])
    obsnum[81]<-sum(u[,2]<u[,5]&u[,5]<u[,1]&u[,1]<u[,4]&u[,4] <u[,3])
    obsnum[82]<-sum(u[,2]<u[,5]&u[,5]<u[,3]&u[,3]<u[,4]&u[,4] <u[,1])
    obsnum[83]<-sum(u[,3]<u[,5]&u[,5]<u[,2]&u[,2]<u[,4]&u[,4] <u[,1])
    obsnum[84]<-sum(u[,3]<u[,5]&u[,5]<u[,1]&u[,1]<u[,4]&u[,4] <u[,2])
    obsnum[85]<-sum(u[,1]<u[,5]&u[,5]<u[,4]&u[,4]<u[,2]&u[,2] <u[,3])
    obsnum[86]<-sum(u[,1]<u[,5]&u[,5]<u[,4]&u[,4]<u[,3]&u[,3] <u[,2])
    obsnum[87]<-sum(u[,2]<u[,5]&u[,5]<u[,4]&u[,4]<u[,1]&u[,1] <u[,3])
    obsnum[88]<-sum(u[,2]<u[,5]&u[,5]<u[,4]&u[,4]<u[,3]&u[,3] <u[,1])
    obsnum[89]<-sum(u[,3]<u[,5]&u[,5]<u[,4]&u[,4]<u[,2]&u[,2] <u[,1])
    obsnum[90]<-sum(u[,3]<u[,5]&u[,5]<u[,4]&u[,4]<u[,1]&u[,1] <u[,2])
    obsnum[91]<-sum(u[,4]<u[,5]&u[,5]<u[,1]&u[,1]<u[,2]&u[,2] <u[,3])
    obsnum[92]<-sum(u[,4]<u[,5]&u[,5]<u[,1]&u[,1]<u[,3]&u[,3] <u[,2])
    obsnum[93]<-sum(u[,4]<u[,5]&u[,5]<u[,2]&u[,2]<u[,1]&u[,1] <u[,3])
    obsnum[94]<-sum(u[,4]<u[,5]&u[,5]<u[,2]&u[,2]<u[,3]&u[,3] <u[,1])
    obsnum[95]<-sum(u[,4]<u[,5]&u[,5]<u[,3]&u[,3]<u[,2]&u[,2] <u[,1])
    obsnum[96]<-sum(u[,4]<u[,5]&u[,5]<u[,3]&u[,3]<u[,1]&u[,1] <u[,2])        
    obsnum[97]<-sum(u[,5]<u[,1]&u[,1]<u[,2]&u[,2]<u[,3]&u[,3] <u[,4])
    obsnum[98]<-sum(u[,5]<u[,1]&u[,1]<u[,3]&u[,3]<u[,2]&u[,2] <u[,4])
    obsnum[99]<-sum(u[,5]<u[,2]&u[,2]<u[,1]&u[,1]<u[,3]&u[,3] <u[,4])
    obsnum[100]<-sum(u[,5]<u[,2]&u[,2]<u[,3]&u[,3]<u[,1]&u[,1] <u[,4])
    obsnum[101]<-sum(u[,5]<u[,3]&u[,3]<u[,2]&u[,2]<u[,1]&u[,1] <u[,4])
    obsnum[102]<-sum(u[,5]<u[,3]&u[,3]<u[,1]&u[,1]<u[,2]&u[,2] <u[,4])
    obsnum[103]<-sum(u[,5]<u[,1]&u[,1]<u[,2]&u[,2]<u[,4]&u[,4] <u[,3])
    obsnum[104]<-sum(u[,5]<u[,1]&u[,1]<u[,3]&u[,3]<u[,4]&u[,4] <u[,2])
    obsnum[105]<-sum(u[,5]<u[,2]&u[,2]<u[,1]&u[,1]<u[,4]&u[,4] <u[,3])
    obsnum[106]<-sum(u[,5]<u[,2]&u[,2]<u[,3]&u[,3]<u[,4]&u[,4] <u[,1])
    obsnum[107]<-sum(u[,5]<u[,3]&u[,3]<u[,2]&u[,2]<u[,4]&u[,4] <u[,1])
    obsnum[108]<-sum(u[,5]<u[,3]&u[,3]<u[,1]&u[,1]<u[,4]&u[,4] <u[,2])
    obsnum[109]<-sum(u[,5]<u[,1]&u[,1]<u[,4]&u[,4]<u[,2]&u[,2] <u[,3])
    obsnum[110]<-sum(u[,5]<u[,1]&u[,1]<u[,4]&u[,4]<u[,3]&u[,3] <u[,2])
    obsnum[111]<-sum(u[,5]<u[,2]&u[,2]<u[,4]&u[,4]<u[,1]&u[,1] <u[,3])
    obsnum[112]<-sum(u[,5]<u[,2]&u[,2]<u[,4]&u[,4]<u[,3]&u[,3] <u[,1])
    obsnum[113]<-sum(u[,5]<u[,3]&u[,3]<u[,4]&u[,4]<u[,2]&u[,2] <u[,1])
    obsnum[114]<-sum(u[,5]<u[,3]&u[,3]<u[,4]&u[,4]<u[,1]&u[,1] <u[,2])
    obsnum[115]<-sum(u[,5]<u[,4]&u[,4]<u[,1]&u[,1]<u[,2]&u[,2] <u[,3])
    obsnum[116]<-sum(u[,5]<u[,4]&u[,4]<u[,1]&u[,1]<u[,3]&u[,3] <u[,2])
    obsnum[117]<-sum(u[,5]<u[,4]&u[,4]<u[,2]&u[,2]<u[,1]&u[,1] <u[,3])
    obsnum[118]<-sum(u[,5]<u[,4]&u[,4]<u[,2]&u[,2]<u[,3]&u[,3] <u[,1])
    obsnum[119]<-sum(u[,5]<u[,4]&u[,4]<u[,3]&u[,3]<u[,2]&u[,2] <u[,1])
    obsnum[120]<-sum(u[,5]<u[,4]&u[,4]<u[,3]&u[,3]<u[,1]&u[,1] <u[,2]) 
  }
  
  #If d is greather than 5, compute all the permutation recursively
  if(d > 5 && d <= 8)
  {
    mypermut <- permut(d)
    obsnum <- u[, mypermut[,1]] < u[, mypermut[,2]]
    
    #Compare columns of u 
    for(i in 3:d)
    {
      obsnum <- obsnum & u[,mypermut[,i-1]] < u[,mypermut[,i]]
    }
    if(NCOL(obsnum) == 1)
      obsnum <- sum(obsnum)
    else
      obsnum <- colSums(obsnum)
  }
  
  #Compute expected numbers
  expnum <- length(u[,1])/factOfD
  
  #Compute chisquare statistic
  residu <- (obsnum - expnum)/sqrt(expnum) 
  stat <- sum(residu^2)
  pvalue <- pchisq(stat,factOfD-1, lower.tail=FALSE)
  
  options(digits=2)    
  if( echo )
  {
    cat("\n\t\t\t Order test\n")
    cat("\nchisq stat = ", stat, ", df = ",factOfD-1, ", p-value = ", pvalue, "\n", sep="")
    cat("\n\t\t (sample size : ",length(u),")\n\n", sep="")
    if(length(obsnum) <= 1000)
      cat("\tobserved number\t",obsnum,"\n")
    else
      cat("\tobserved number\t too many to be printed\n")
    cat("\texpected number\t",expnum,"\n")    
  } 
  
  res <- list(statistic=stat, parameter=factOfD-1, 
              p.value=pvalue, observed=obsnum, 
              expected=expnum, residuals=residu) 
  return(invisible(res))
}
