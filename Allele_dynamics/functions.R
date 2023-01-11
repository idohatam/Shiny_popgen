
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
dnp <- function(qd,td,sd=NULL, scd, g, ps){
  require("dplyr")
  require("tidyr")
  require("ggplot2")
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
        } else {if(scd == "Overderdominance of the heterozygot"|
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
    }else{ a <- a}
  A <- 1 - a
  generations <- generations[1:length(a)]
  
  #genotypes proportion
  prop_AA <- A^2
  prop_aa <- a^2
  prop_Aa <- 2*a*A
  
  #genotype size
  
  size_AA <- ps*prop_AA
  size_aa <- ps*prop_aa
  size_Aa <- ps*prop_Aa
  
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
