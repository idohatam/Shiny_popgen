
#Calculate change in q (new allele) under overdominance scenario.
#t selection coefficient of homozygote, s selection coefficient for heterozygote
overd <- function(q,s,t){
  p <- 1-q
  dq <- (p*q*(s-(2*q*s)+t*q))/(1+(2*s*q*p)-t*(q^2))
  nq <- q +dq
  return(nq)
  
}

#Calculate the allele frequency under equilibrium in the case of overdominanace
overeq <- function(s,t){
  
  qeq <- s/(2*s-t)
  peq <- 1-qeq
  return(c(qeq,peq))
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

#allele dynamics under randome genetic drift
rgd <- function(eps, q, ng, ns) {
  # Initialize empty data frame to store results
  results <- data.frame(generation = numeric(0), A1 = numeric(0), A2 = numeric(0), simulation = numeric(0))
  
  # Loop through number of simulations
  for (i in 1:ns) {
    # Calculate proportion of A2
    p <- 1 - q
    # Initialize vectors to store proportion of each allele for each generation
    q_vec <- numeric(ng)
    p_vec <- numeric(ng)
    # Set initial proportion of alleles
    q_vec[1] <- q
    p_vec[1] <- p
    # Loop through generations
    for (j in 2:ng) {
      # Calculate change in allele proportion based on genetic drift
      q_vec[j] <- q_vec[j-1] + rnorm(1, mean = 0, sd = sqrt(q_vec[j-1]*p_vec[j-1]/eps))
      # Truncate values less than 0 or greater than 1
      q_vec[j] <- ifelse(q_vec[j] < 0, 0, ifelse(q_vec[j] > 1, 1, q_vec[j]))
      p_vec[j] <- 1 - q_vec[j]
    }
    # Add results to data frame
    results <- rbind(results, data.frame(generation = 1:ng, A1 = q_vec, A2 = p_vec, simulation = i,stringsAsFactors = F))
  }
  results$simulation <- as.factor(results$simulation)
  return(results)
}


#data and plots
dnp <- function(qd,td,sd=NULL, scd, g, ps){
  #Create a vector a that would contain the proportions of the new allele
  a <- rep(0,g)
  
  #Create the vector generations that would hold the number of generations passed
  generations <- seq(1, g, by = 1)
  
  
  i=1
  #Loop to get the change in q over time depending on the scenario
  for(i in seq_along(generations)){
    if(scd == "The new allele is recessive"){
      if(i == 1){
          a[i] <- qd
        } else {
          a[i] <- res(q = a[i-1], t = td)
        }
    } else {
      if(scd == "The new allele is dominanat"){
        if(i == 1){
          a[i] <- qd
        } else {
          a[i] <- dom(q = a[i-1], t = td)
        }
      } else {
        if(scd == "The alleles are codominant"){
          if(i == 1){
            a[i] <- qd
          } else {
            a[i] <- co(q = a[i-1], t = td)
          }
        } else {if(scd == "Overdominance of the heterozygot"|
                   scd == "Underdominance of the heterozygot"){ 
          if( i == 1){
            a[i] <- qd
          } else {
            a[i] <- overd(q = a[i-1], s = sd, t = td)}
          }
        }
      }
    }
  }
    
  
  # shrinking the dataset
  #Need to filter the heterozygot and negative values
  if(!grepl("heterozygot", scd) & td > 0){
    a <- a[1 : which(round(a,3) > 0.99)[10]]
  }else{
      if(!grepl("heterozygot", scd) & td < 0){
        a <- a[1:which(round(a,4) < 0.03)[3]]
      } else{if(grepl("Overdominance", scd)){
        
        eq <- overeq(s = sd, t = td)
        
        if(eq[1] == 1){
          a <- a[1:which((1-round(a,4)) < 0.035)[3]]
        }else{
          a <- a[1:which(round(a,2) == round(eq[1],2))[3]]
        }
      }else{
        if(grepl("Underdominance", scd) & td > 0){
          
          if(max(a) > qd){
            a <- a[1:which(round(a,2) > 0.99)[10]]  
          } else {
            a <- a[1:which(round(a,4) < 0.03)[10]]
          }
          
        } else {
          if(grepl("Underdominance", scd) & td < 0){
            a <- a[1:which(round(a,4)<0.001)[3]]
          }else{
            a <- a
          }
        }
         }
        }
    }
  A <- 1 - a
  generations <- generations[1:length(a)]
  
  #genotypes proportion
  prop_AA <- A^2
  prop_aa <- a^2
  prop_Aa <- 2*a*A
  
  #genotype size
  
  size_AA <- round(ps*prop_AA, 0)
  size_aa <- round(ps*prop_aa, 0)
  size_Aa <- round(ps*prop_Aa, 0)
  
  #make dataframes
  
  df_allele <- tibble(a,A, gen = generations) %>% 
    gather(key = "Allele", value = "Proportion", -gen)
  
  df_gtype_prop <- tibble(AA = prop_AA, aa = prop_aa, Aa = prop_Aa, gen = generations) %>%
    gather(key = "Genotype", value = "Proportion", -gen)
  
  df_gtyp_size <- tibble(AA = size_AA, aa = size_aa, Aa = size_Aa, gen = generations) %>%
    gather(key = "Genotype", value = "Size", -gen)
  
  #chi square get generation for loop
  
  gts <- df_gtyp_size %>% pivot_wider(names_from = Genotype, values_from = Size)
  prop <- as.numeric(c(prop_AA[1], prop_aa[1], prop_Aa[1]))
  
  for (gn in seq_along(gts$gen)) {
    csq2 <- chisq.test(x = round(as.numeric(gts[gn,2:4]),0), 
                       p = prop)
    if(round(csq2$statistic, 3) >= 7.373 ){
      HWE <- gn
      break
    }
    
  }
 
  
  
  
  #Plot allele
  p_alleles <- ggplot(df_allele, aes(
    x = gen, y = Proportion, color = Allele))+
    geom_line(linewidth=0.5)+
    geom_point()+
    geom_vline(xintercept = HWE, linetype = 2, color = "black")+
    theme_bw(base_size = 16)+
    scale_color_brewer(palette = "Dark2")+
    xlab("Generations")+
    ylab("Allele proportion within the population")+
    ggtitle("Change in allele proportion")
  
  #Plot gen_freq
  
  p_gen_freq <- ggplot(df_gtype_prop, aes(
    x = gen, y = Proportion, color = Genotype))+
    geom_line(linewidth=0.5)+
    geom_point()+
    geom_vline(xintercept = HWE, linetype = 2, color = "black")+
    theme_bw(base_size = 16)+
    scale_color_brewer(palette = "Dark2")+
    xlab("Generations")+
    ylab("Genotype proportion within the population")+
    ggtitle("Change in genotype proportion")
  
  #Plot size
  
  p_gen_size <- ggplot(df_gtyp_size, aes(
    x = gen, y = Size, fill = Genotype,color = Genotype))+
    geom_line(linewidth = 0.5)+
    geom_point()+
    geom_vline(xintercept = HWE, linetype = 2, color = "black")+
    theme_bw(base_size = 16)+
    scale_color_brewer(palette = "Dark2")+
    scale_fill_brewer(palette = "Dark2")+
    xlab("Generations")+
    ylab("Number of individuals")+
    ggtitle("Change in number of individuals per genotype")
  
  #return the results list
  
  results <- list(df_allele, df_gtype_prop, df_gtyp_size, 
                  p_alleles, p_gen_freq, p_gen_size)
  
  
  return(results)
}
