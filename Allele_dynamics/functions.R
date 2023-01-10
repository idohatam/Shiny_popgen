
#Calculate change in q (new allele) under overdominance scenario.
#t selection coefficient of homozygote, s selection coefficient for heterozygote
overd <- function(q,s,t){
  p <- 1-q
  dq <- (p*q(s-(2*q*s)+t*q))/(1+(2*s*q*p)-t*(q^2))
  nq <- q +dq
  return(nq)
  
}

#Calculate the allele frequency under equilibrium in the case of overdominanace
overeq <- function(s,t){
  
  qeq <- s/(2*s-t)
  peq <- 1-qeq
}

#Allele dynamics in the case of positive recessive 
res <- function(q,t){
  p <- 1-q
  dq <- (t*p*q^2)/(1+(t*q^2))
  nq <- q+dq
  return(nq)
}

#Allele dynamics in the case of positive dominant

dom <- function(q,t){
  p <- 1-q
  
  dq <- (t*p*q^2)/(1+ q*t*(1+p))
  
  nq <- q+dq
  
  return(nq)
}

#Allele dynamics in the case of codominanace

co <- function(q,t){
  p = 1-q
  
  dq <- (t*p*q)/(1+ 2*q*t)
  
  nq <- q+dq
  
  return(nq)
  
}


#data and plots

dnp <- function(qd,td,sd=NULL, scd, g){
  a <- rep(0,g)
  generations <- seq(1, g, by = 1)
  i=1
  
  if(sc == "The new allele is recessive"){
    for (i in seq_along(a)) {
      
      if(i == 1){
        a[i] <- qd
      } else {
        a[i] <- a[i-1] + res(q = qd, t = td)
      }
      
    }
    
  }
  
  return(results)
}