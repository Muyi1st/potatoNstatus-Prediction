library(tidyverse)
library(plotly)
library(GGally)
library(StageWise)

#load FM data
load("./data/chipsNAllSpectraits.rda")

ch <- chips %>%
  #mutate_if(is.character, as.numeric) %>%
  mutate_at(vars(23:35), as.numeric) %>%
  mutate(nlevel = factor(nlevel, levels = c("0","33","66", "100"), labels = c("5.6","73.97","147.95","229.77"))) %>% # replace nlevels with applied n
  mutate(vinedm=(subvinedry/subvinefresh), tuberdm= (subtuberdry/subtuberfresh)) %>% #DM estimation for both vine and tuber
  mutate(ncontVine=vinedm*nvine, ncontTub=tuberdm*ntuber) %>% # N concentraion estimation for both vine and tuber
  mutate(tpn=ncontTub+ncontVine) %>% # total plant N estimation
  mutate(nup=tpn/as.numeric(nlevel)) %>% # N uptake estimation
  mutate(ute=tuberdm/tpn) %>% # N utilization estimation
  mutate(nue=ute*nup) %>% # N use efficiency estimation
  dplyr::select(-vinefresh, -subvinefresh, -subvinedry, -tuberfresh, -subtuberfresh, -subtuberdry) # remove columns no longer needed

# view data distribution
par(mfrow=c(3,6))
for(i in 10:ncol(ch)){
  hist(ch[,i], xlab = '', col='steelblue',  main=names(ch[i]))
}

#funcion for data cleaning based on SD deviation
sdfx <- function(col, na.rm=T){
  #filter(is.finite(col))
  d1 <- mean(col,na.rm=T) + 2*sd(col,na.rm = T)
  d2 <- mean(col,na.rm=T) - 2*sd(col,na.rm = T)
  col <- ifelse(col>d1 | col<d2, NA, col)
  #col[col>d1 | col<d2] <- NA
  return(col)
}


#data cleaning and processing of chips NCurve
ntraits <- c(names(ch)[c(26:37)]) #
#for (i in (traits)) {
ch1 <- ch %>%
  filter(spec_date=="22_cover50" | spec_date=="23_cover50") %>%
  #fm1 <- data.frame(lapply(fm1[,c(10:34)], sdfx))
  mutate(across(everything(), ~ifelse(is.infinite(.),NA, .))) %>%
  mutate(across(all_of(ntraits),sdfx)) %>%
  mutate(redness=ifelse(redness<0 | redness>33, NA,redness),
         lightness=ifelse(lightness<0 | lightness>100, NA,lightness))

for(i in 10:ncol(ch1)){
  meanXi <- mean(ch1[,i],na.rm=TRUE)
  # naive imputation
  ch1[,i] <- ifelse(is.na(ch1[,i]),meanXi, ch1[,i])
}


# view data distribution after cleaning
par(mfrow=c(3,6))
Hmisc::hist.data.frame(ch1[,10:37])

#format data for generating blues and heritaility
ch1 %<>%
  rename(env=trial) %>%
  dplyr::select(env, id, plot,range, row, rep, nlevel, 10:37) %>%
  write_csv(., "./data/chcanopy50.csv", quote = "none")


#########obtain blues
#library(StageWise)

effects <- data.frame(name=c("row","range","rep"),
                      fixed=c(FALSE,FALSE,FALSE),
                      factor=c(TRUE,TRUE,TRUE))
effects
traits <- c(names(ch1)[c(8:35)])
mod <- vector("list", length(traits))

for (i in traits) {
  mod[[i]] <- StageWise::Stage1(filename="./data/chcanopy50.csv", traits=i, solver="spats", spline=c("row","range"),effects=effects[3,])
}

saveRDS(mod,"./data/chcanopy50allTraitsBlues.rds")


#get blues, heritability, vcov matrix and residuals
mod <- readRDS("./data/chcanopy50allTraitsBlues.rds")
ch_canopy50_blues <- data.frame(env=mod$red$blues$env, id=mod$red$blues$id,  red=mod$red$blues$BLUE,green=mod$green$blues$BLUE,
                                blue=mod$blue$blues$BLUE, nir=mod$nir$blues$BLUE, red_edge=mod$red_edge$blues$BLUE, ndvi=mod$ndvi$blues$BLUE,
                                ndre=mod$ndre$blues$BLUE, gndvi=mod$gndvi$blues$BLUE, gli=mod$gli$blues$BLUE, cire=mod$cire$blues$BLUE,
                                cig=mod$cig$blues$BLUE, ndwi=mod$ndwi$blues$BLUE, tcari.osavi=mod$tcari.osavi$blues$BLUE,nvine=mod$nvine$blues$BLUE,
                                ntuber=mod$ntuber$blues$BLUE,vinedm=mod$vinedm$blues$BLUE, tuberdm=mod$tuberdm$blues$BLUE, ncontVine=mod$ncontVine$blues$BLUE,
                                nconTuber=mod$ncontTub$blues$BLUE,tpn=mod$tpn$blues$BLUE,nup=mod$nup$blues$BLUE, ute=mod$ute$blues$BLUE, nue=mod$nue$blues$BLUE)


c50_byield <-  data.frame(env=mod$yield$blues$env, id=mod$yield$blues$id, yield=mod$yield$blues$BLUE)
c50_bred <-  data.frame(env=mod$redness$blues$env,id=mod$redness$blues$id,redness=mod$redness$blues$BLUE)
c50_blight <-  data.frame(env=mod$lightness$blues$env,id=mod$lightness$blues$id,lightness=mod$lightness$blues$BLUE)
c50_round <-  data.frame(env=mod$roundness$blues$env,id=mod$roundness$blues$id,roundness=mod$roundness$blues$BLUE)
c50_lxw <-  data.frame(env=mod$lxw$blues$env,id=mod$lxw$blues$id,lxw=mod$lxw$blues$BLUE)
data1 <- f50_round %>%
  left_join(f50_lxw, by=c("env","id")) %>%
  left_join(f50_blight, by=c("env","id")) %>%
  left_join(f50_bred, by=c("env","id")) %>%
  left_join(f50_byield, by=c("env","id"))

save(data1, file = "./data/traitBLUEschYRLRL.rdata")

ch_canopy50_H2 <- data.frame(env=mod$red$fit$env, yield=mod$yield$fit$H2, redness=mod$redness$fit$H2,lightness=mod$lightness$fit$H2,
                             roundness=mod$roundness$fit$H2,lxw=mod$lxw$fit$H2,nvine=mod$nvine$fit$H2,ntuber=mod$ntuber$fit$H2,
                             vinedm=mod$vinedm$fit$H2, tuberdm=mod$tuberdm$fit$H2, ncontVine=mod$ncontVine$fit$H2,
                             nconTuber=mod$ncontTub$fit$H2,tpn=mod$tpn$fit$H2,nup=mod$nup$fit$H2,ute=mod$ute$fit$H2,nue=mod$nue$fit$H2,
                             red=mod$red$fit$H2, green=mod$green$fit$H2,blue=mod$blue$fit$H2,
                             nir=mod$nir$fit$H2, red_edge=mod$red_edge$fit$H2, ndvi=mod$ndvi$fit$H2,ndre=mod$ndre$fit$H2,gndvi=mod$gndvi$fit$H2,
                             gli=mod$gli$fit$H2, cire=mod$cire$fit$H2, cig=mod$cig$fit$H2,ndwi=mod$ndwi$fit$H2,tcari.osavi=mod$tcari.osavi$fit$H2)




#fm_canopy50_vcov <- list(yield=mod$yield$vcov,redness=mod$redness$vcov, lightness=mod$lightness$vcov, roundness=mod$roundness$vcov
#                         ,lxw=mod$lxw$vcov)

#####saving traits data separately
save(ch_canopy50_blues, file = "./data/spectcanopy50BLUEsch.rdata")
save(fm_canopy50_vcov, file = "./data/traitVcovchRLRL.rdata")



###########heritability
h2_fc50 <- ch_canopy50_H2 %>%
  gather(key = "traits", value = "h2", -c(env)) %>%
  mutate(gs=paste0("canopy cover 50%"))




#########################
########################## Canopy cover peak flowering
ch2 <-  ch %>%
  #  mutate(yield=ifelse(spec_date=="22_cover50",yield22$yield, yield)) %>% #correct yield data
  filter(spec_date=="22_coverflower" | spec_date=="23_coverflower") %>%
  mutate(across(everything(), ~ifelse(is.infinite(.),NA, .))) %>%
  mutate(across(all_of(ntraits),sdfx)) %>%
  mutate(redness=ifelse(redness<0 | redness>33, NA,redness),
         lightness=ifelse(lightness<0 | lightness>100, NA,lightness))

for(i in 10:ncol(ch2)){
  meanXi <- mean(ch2[,i],na.rm=TRUE)
  # naive imputation
  ch2[,i] <- ifelse(is.na(ch2[,i]),meanXi, ch2[,i])
}


par(mfrow=c(3,5))
Hmisc::hist.data.frame(ch2[,10:37])


ch2 %<>%
  rename(env=trial) %>%
  dplyr::select(env, id, plot,range, row, rep, nlevel, 10:37) %>%
  write_csv(., "./data/chcanopycfl.csv", quote = "none")

#########obtain blues
effects <- data.frame(name=c("row","range","rep"),
                      fixed=c(FALSE,FALSE,FALSE),
                      factor=c(TRUE,TRUE,TRUE))
effects
traits <- c(names(ch2)[c(8:20)])
mod <- vector("list", length(traits))
#heritability.c1 <- data.frame(env=factor(), H2=integer(), trait=character())
#blues.c1 <- chips1.blues[,1:5]
#mod <- list(14)

for (i in traits) {
  mod[[i]] <- StageWise::Stage1(filename="./data/chcanopycfl.csv", traits=i, solver="spats", spline=c("row","range"),effects=effects[3,])
}

saveRDS(mod,"./data/fmcanopycoverflowerTraitsBlues.rds")


#get blues, heritability, vcov matrix and residuals
mod <- readRDS("./data/fmcanopycoverflowerTraitsBlues.rds")
ch_canopycfl_blues <- data.frame(env=mod$red$blues$env, id=mod$red$blues$id,  red=mod$red$blues$BLUE,green=mod$green$blues$BLUE,
                                 blue=mod$blue$blues$BLUE, nir=mod$nir$blues$BLUE, red_edge=mod$red_edge$blues$BLUE, ndvi=mod$ndvi$blues$BLUE,
                                 ndre=mod$ndre$blues$BLUE, gndvi=mod$gndvi$blues$BLUE, gli=mod$gli$blues$BLUE, cire=mod$cire$blues$BLUE,
                                 cig=mod$cig$blues$BLUE, ndwi=mod$ndwi$blues$BLUE,tcari.osavi=mod$tcari.osavi$blues$BLUE)


ch_canopycfl_H2 <- data.frame(env=mod$red$fit$env,red=mod$red$fit$H2, green=mod$green$fit$H2,blue=mod$blue$fit$H2,
                              nir=mod$nir$fit$H2, red_edge=mod$red_edge$fit$H2, ndvi=mod$ndvi$fit$H2,ndre=mod$ndre$fit$H2,gndvi=mod$gndvi$fit$H2,
                              gli=mod$gli$fit$H2, cire=mod$cire$fit$H2, cig=mod$cig$fit$H2,ndwi=mod$ndwi$fit$H2,tcari.osavi=mod$tcari.osavi$fit$H2)

#fm_canopycfl_vcov <- list(yield=mod$yield$vcov,redness=mod$redness$vcov, lightness=mod$lightness$vcov, roundness=mod$roundness$vcov
#                         ,lxw=mod$lxw$vcov)

#####saving traits data separately
save(ch_canopycfl_blues, file = "./data/spectcanopycflBLUEsch.rdata")
#save(fm_canopy50_vcov, file = "./data/traitVcovfmRLRL.rdata")


###########heritability
h2_fcfl <- ch_canopycfl_H2 %>%
  gather(key = "traits", value = "h2", -c(env)) %>%
  mutate(gs=paste0("canopy cover flower"))

############plot h2
h2_fcfl%>%
  mutate(traits = factor(traits, levels = unique(traits)), gs = factor(gs, levels = c("canopy cover flower"), labels = c("canopy cover flower"))) %>%
  ggplot(aes(x = traits,y = h2, fill=env)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5)+
  # geom_text(aes(x = model,y = PVE,label = paste(PVE*100, "%")),
  #           position = position_stack(vjust = 0.5))+
  labs(x = "Traits", y = "H2") + facet_grid(. ~ gs)



#########################
########################## Peak flowering
ch3 <-  ch %>%
  filter(spec_date=="22_peakFlower" | spec_date=="23_peakFlower") %>%
  mutate(across(everything(), ~ifelse(is.infinite(.),NA, .))) %>%
  mutate(across(all_of(ntraits),sdfx)) %>%
  mutate(redness=ifelse(redness<0 | redness>33, NA,redness),
         lightness=ifelse(lightness<0 | lightness>100, NA,lightness))

for(i in 10:ncol(ch3))
{
  meanXi <- mean(ch3[,i],na.rm=TRUE)
  # naive imputation
  ch3[,i] <- ifelse(is.na(ch3[,i]),meanXi, ch3[,i])
}

par(mfrow=c(3,5))
for(i in 7:19){
  hist(ch3[,i], xlab = '', col='steelblue',  main=names(ch3[i]))
}



##saving info for the chips canopy50 spectra data for downstream analysis
ch3 %<>%
  dplyr::rename(env=trial) %>%
  dplyr::select(env, id, plot,range, row, rep, nlevel, 10:37) %>%
  write_csv(., "./data/chpeakFlower.csv", quote = "none")



#########obtain blues
effects <- data.frame(name=c("row","range","rep"),
                      fixed=c(FALSE,FALSE,FALSE),
                      factor=c(TRUE,TRUE,TRUE))
effects

traits <- c(names(ch3)[8:20])
mod <- vector("list", length(traits))
#heritability.c1 <- data.frame(env=factor(), H2=integer(), trait=character())
#blues.c1 <- chips1.blues[,1:5]
#mod <- list(14)

for (i in traits) {
  mod[[i]] <- Stage1(filename="./data/chpeakFlower.csv", traits=i, solver="spats", spline=c("row","range"),effects=effects[3,])
}

#yield=mod$yield$vcov$202
saveRDS(mod3,"./data/peakFlowerspecBluesch.rds")


mod3 <- readRDS("./data/peakFlowerspecBluesch.rds")


## get blues, heritability, vcov matrix and residuals

ch_peakFlower_blues <- data.frame(env=mod$red$blues$env, id=mod$red$blues$id,  red=mod$red$blues$BLUE,green=mod$green$blues$BLUE,
                                  blue=mod$blue$blues$BLUE, nir=mod$nir$blues$BLUE, red_edge=mod$red_edge$blues$BLUE, ndvi=mod$ndvi$blues$BLUE,
                                  ndre=mod$ndre$blues$BLUE, gndvi=mod$gndvi$blues$BLUE, gli=mod$gli$blues$BLUE, cire=mod$cire$blues$BLUE,
                                  cig=mod$cig$blues$BLUE, ndwi=mod$ndwi$blues$BLUE, tcari.osavi=mod$tcari.osavi$blues$BLUE)

ch_peakFlower_H2 <- data.frame(env=mod$red$fit$env,red=mod$red$fit$H2, green=mod$green$fit$H2,blue=mod$blue$fit$H2,
                               nir=mod$nir$fit$H2, red_edge=mod$red_edge$fit$H2, ndvi=mod$ndvi$fit$H2,ndre=mod$ndre$fit$H2,gndvi=mod$gndvi$fit$H2,
                               gli=mod$gli$fit$H2, cire=mod$cire$fit$H2, cig=mod$cig$fit$H2,ndwi=mod$ndwi$fit$H2,tcari.osavi=mod$tcari.osavi$fit$H2)
# fm_canopy50_vcov <- as.matrix(env=mod3$red$vcov,  red=mod3$red$vcov, green=mod3$green$vcov, blue=mod3$blue$vcov,
#                                  nir=mod3$nir$vcov, red_edge=mod3$red_edge$vcov, ndvi=mod3$ndvi$vcov, ndre=mod3$ndre$vcov,gndvi=mod3$gndvi$vcov,
#                                  gli=mod3$gli$vcov, cire=mod3$cire$vcov, cig=mod3$cig$vcov, ndwi=mod3$ndwi$vcov,tcari.osavi=mod3$tcari.osavi$vcov)

##Save peak fl blues

save(ch_peakFlower_blues, file = "./data/ch_specPeakflBLUEs.rdata")
#save(chips_canopy50_vcov, file = "./data/traitVcovchipsPeakfl.rdata")



##plot heritability for the bands and traits

h2_fpF <- ch_peakFlower_H2 %>%
  gather(key = "traits", value = "h2", -c(env)) %>%
  mutate(gs=paste0("peak flowering"),)

##update heritability for all growth stages
#saveRDS(rbind(h2_ccflower,h2_cpF), "chipsH2traits.rds")




#########################
########################## Vine biomass
ch4 <-  ch %>%
  #  mutate(yield=ifelse(spec_date=="22_cover50",yield22$yield, yield)) %>% #correct yield data
  filter(spec_date=="22_vineBiomass" | spec_date=="23_vineBiomass") %>%
  mutate(across(everything(), ~ifelse(is.infinite(.),NA, .))) %>%
  mutate(across(all_of(ntraits),sdfx)) %>%
  mutate(redness=ifelse(redness<0 | redness>33, NA,redness),
         lightness=ifelse(lightness<0 | lightness>100, NA,lightness))

for(i in 10:ncol(ch4))
{
  meanXi <- mean(ch4[,i],na.rm=TRUE)
  # naive imputation
  ch4[,i] <- ifelse(is.na(ch4[,i]),meanXi, ch4[,i])
}
par(mfrow=c(3,5))
for(i in 8:20){
  hist(ch4[,i], xlab = '', col='steelblue',  main=names(ch4[i]))
}



## saving info for the chips canopy50 spectra data for downstream analysis
ch4 %<>%
  dplyr::rename(env=trial) %>%
  dplyr::select(env, id, plot,range, row, rep, nlevel, 10:37) %>%
  write_csv(., "./data/chvineBiomass.csv", quote = "none")


traits <- c(names(ch4)[8:20])
mod <- vector("list", length(traits))
#heritability.c1 <- data.frame(env=factor(), H2=integer(), trait=character())
#blues.c1 <- chips1.blues[,1:5]
#mod <- list(14)

for (i in traits) {
  mod[[i]] <- Stage1(filename="./data/chvineBiomass.csv", traits=i, solver="spats", spline=c("row","range"),effects=effects[3,])
}

#yield=mod$yield$vcov$202
saveRDS(mod4,"./data/vineBiomassaTraitsBluesch.rds")




mod4 <- readRDS("./data/vineBiomassTraitsBluesch.rds")


## get blues, heritability, vcov matrix and residuals

ch_vineBiomass_blues <- data.frame(env=mod$red$blues$env, id=mod$red$blues$id,  red=mod$red$blues$BLUE,green=mod$green$blues$BLUE,
                                   blue=mod$blue$blues$BLUE, nir=mod$nir$blues$BLUE, red_edge=mod$red_edge$blues$BLUE, ndvi=mod$ndvi$blues$BLUE,
                                   ndre=mod$ndre$blues$BLUE, gndvi=mod$gndvi$blues$BLUE, gli=mod$gli$blues$BLUE, cire=mod$cire$blues$BLUE,
                                   cig=mod$cig$blues$BLUE, ndwi=mod$ndwi$blues$BLUE, tcari.osavi=mod$tcari.osavi$blues$BLUE)


ch_vineBiomass_H2 <- data.frame(env=mod$red$fit$env,red=mod$red$fit$H2, green=mod$green$fit$H2,blue=mod$blue$fit$H2,
                                nir=mod$nir$fit$H2, red_edge=mod$red_edge$fit$H2, ndvi=mod$ndvi$fit$H2,ndre=mod$ndre$fit$H2,gndvi=mod$gndvi$fit$H2,
                                gli=mod$gli$fit$H2, cire=mod$cire$fit$H2, cig=mod$cig$fit$H2,ndwi=mod$ndwi$fit$H2,tcari.osavi=mod$tcari.osavi$fit$H2)
# chips_vineBiomass_vcov <- as.matrix(env=mod4$red$vcov,  red=mod4$red$vcov, green=mod4$green$vcov, blue=mod4$blue$vcov,
#                                 nir=mod4$nir$vcov, red_edge=mod4$red_edge$vcov, ndvi=mod4$ndvi$vcov, ndre=mod4$ndre$vcov,gndvi=mod4$gndvi$vcov,
#                                 gli=mod4$gli$vcov, cire=mod4$cire$vcov, cig=mod4$cig$vcov, ndwi=mod4$ndwi$vcov,tcari.osavi=mod4$tcari.osavi$vcov)

##Save Vine biomass blues

save(ch_vineBiomass_blues, file = "./data/ch_spectvinebioBLUEs.rdata")


##plot heritability for the bands and traits
h2_fvB <- fm_vineBiomass_H2 %>%
  gather(key = "traits", value = "h2", -c(env)) %>%
  mutate(gs=paste0("vine biomass"),)
#saveRDS(rbind(h2_ccflower,h2_cpF,h2_cvB), "./data/chipsH2traits.rds")




#########################
########################## Senescence 1
ch5 <-  ch %>%
  #  mutate(yield=ifelse(spec_date=="22_cover50",yield22$yield, yield)) %>% #correct yield data
  filter(spec_date=="22_senesence1" | spec_date=="23_senesence1") %>%
  mutate(across(everything(), ~ifelse(is.infinite(.),NA, .))) %>%
  mutate(across(all_of(ntraits),sdfx)) %>%
  mutate(redness=ifelse(redness<0 | redness>33, NA,redness),
         lightness=ifelse(lightness<0 | lightness>100, NA,lightness))

for(i in 10:ncol(ch5)){
  meanXi <- mean(ch5[,i],na.rm=TRUE)
  # naive imputation
  ch5[,i] <- ifelse(is.na(ch5[,i]),meanXi, ch5[,i])
}

par(mfrow=c(3,5))
for(i in 7:20){
  hist(ch5[,i], xlab = '', col='steelblue',  main=names(ch5[i]))
}


ch5 %<>%
  dplyr::rename(env=trial) %>%
  dplyr::select(env, id, plot,range, row, rep, nlevel, 10:37) %>%
  write_csv(., "./data/chsenescence1.csv", quote = "none")


traits <- c(names(ch5)[8:20])
mod <- vector("list", length(traits))
#heritability.c1 <- data.frame(env=factor(), H2=integer(), trait=character())
#blues.c1 <- chips1.blues[,1:5]
#mod <- list(14)

for (i in traits) {
  mod[[i]] <- Stage1(filename="./data/chsenescence1.csv", traits=i, solver="spats", spline=c("row","range"),effects=effects[3,])
}

#yield=mod$yield$vcov$202
saveRDS(mod5,"./data/senescence1TraitsBluesch.rds")




mod5 <- readRDS("./data/senescence1TraitsBluesch.rds")


## get blues, heritability, vcov matrix and residuals

ch_senescence1_blues <- data.frame(env=mod$red$blues$env, id=mod$red$blues$id,red=mod$red$blues$BLUE,green=mod$green$blues$BLUE,
                                   blue=mod$blue$blues$BLUE, nir=mod$nir$blues$BLUE, red_edge=mod$red_edge$blues$BLUE, ndvi=mod$ndvi$blues$BLUE,
                                   ndre=mod$ndre$blues$BLUE, gndvi=mod$gndvi$blues$BLUE, gli=mod$gli$blues$BLUE, cire=mod$cire$blues$BLUE,
                                   cig=mod$cig$blues$BLUE, ndwi=mod$ndwi$blues$BLUE, tcari.osavi=mod$tcari.osavi$blues$BLUE)


ch_senescence1_H2 <- data.frame(env=mod$red$fit$env,red=mod$red$fit$H2, green=mod$green$fit$H2,blue=mod$blue$fit$H2,
                                nir=mod$nir$fit$H2, red_edge=mod$red_edge$fit$H2, ndvi=mod$ndvi$fit$H2,ndre=mod$ndre$fit$H2,gndvi=mod$gndvi$fit$H2,
                                gli=mod$gli$fit$H2, cire=mod$cire$fit$H2, cig=mod$cig$fit$H2,ndwi=mod$ndwi$fit$H2,tcari.osavi=mod$tcari.osavi$fit$H2)
# chips_canopy50_vcov <- as.matrix(env=mod5$red$vcov,  red=mod5$red$vcov, green=mod5$green$vcov, blue=mod5$blue$vcov,
#                                 nir=mod5$nir$vcov, red_edge=mod5$red_edge$vcov, ndvi=mod5$ndvi$vcov, ndre=mod5$ndre$vcov,gndvi=mod5$gndvi$vcov,
#                                 gli=mod5$gli$vcov, cire=mod5$cire$vcov, cig=mod5$cig$vcov, ndwi=mod5$ndwi$vcov,tcari.osavi=mod5$tcari.osavi$vcov)


##Save Senesce1 blues

save(ch_senescence1_blues, file = "./data/ch_spectSenesce1BLUEs.rdata")


# plot heritability for the bands and traits

h2_fS1 <- ch_senescence1_H2 %>%
  gather(key = "traits", value = "h2", -c(env)) %>%
  mutate(gs=paste0("senescence1"),)





#####################################
#################################################################################

##plot all Broad sense heritability
#update merging
saveRDS(rbind(h2_fc50, h2_fcfl,h2_fpF,h2_fvB,h2_fS1,h2_fS2), "./data/chH2traits.rds")
saveRDS(rbind(ch_H2, h2_fS2),"../data/chH2traits.rds")
ch_H2 <- readRDS("./data/chH2traits.rds")
ch_H2$class <- paste("chips market class")




############plot h2
ch2 <- ch_H2[-(1:40),] %>%
  mutate(traits = factor(traits, levels = unique(traits)), gs = factor(gs, levels = unique(gs)),
         trial = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                        "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                        labels= c("0N_2022","33N_2022","66N_2022","100N_2022","0N_2023","33N_2023","66N_2023","100N_2023"))) %>%

  ggplot(aes(x = traits,y = h2, fill =trial)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5)+
  # geom_text(aes(x = model,y = PVE,label = paste(PVE*100, "%")),
  #           position = position_stack(vjust = 0.5))+
  labs(x = "Ntraits with bands", y = "Broad sense heritability") + facet_wrap(~gs, ncol = 3)+
  theme(axis.text = element_text(face = "bold",size = 10), axis.title = element_text(face = "bold", size=14),
        legend.title = element_text(face = "bold", size=12),legend.text = element_text(face = "bold", size=8),
        strip.text = element_text(size=16,face="bold"),axis.text.x=element_text(angle=25, vjust=0.9, hjust=0.9))
ggsave("./output/fmspectrancurveH2.png",height=10.5, width=18.5, units="in", dpi=300)


## combine all spectral blues and traits

msncdata1.5_ch <- cbind(ch_canopy50_blues[,c(1:2,16:25)] %>%
                       unite(env.id,c(env,id),sep = ":"),
                     ch_canopy50_blues[,3:15]%>%
                       rename_with( ~ paste0(.x, ".1")),
                     ch_canopycfl_blues[,3:15]%>%
                       #unite(env.id,c(env,id),sep = ":") %>%
                       rename_with( ~ paste0(.x, ".2")),
                     ch_peakFlower_blues[,3:15]%>%
                       #unite(env.id,c(env,id),sep = ":") %>%
                       rename_with( ~ paste0(.x, ".3")),
                     ch_vineBiomass_blues[,3:15]%>%
                       #unite(env.id,c(env,id),sep = ":") %>%
                       rename_with( ~ paste0(.x, ".4")),
                     ch_senescence1_blues[,3:15]%>%
                       #unite(env.id,c(env,id),sep = ":") %>%
                       rename_with( ~ paste0(.x, ".5"))
)
save(msncdata1.5_ch, file = "./data/msspecCHdatablues1.5.rda")

spectra <- msncdata1.5_ch[,]
cortab <- data.frame(spectra=NA, r2=NA, year=NA, nrate=NA, gs=NA)

# Plot distribution for total plant content
d1 <- ch_canopy50_blues %>%
  mutate(year = factor(env, levels = c("2022_MN_NCurve_Ch_0","2022_MN_NCurve_Ch_33","2022_MN_NCurve_Ch_66","2022_MN_NCRT_Chips",
                                       "2023_MN_Ncurve_Ch_0","2023_MN_Ncurve_Ch_33", "2023_MN_Ncurve_Ch_66","2023_MN_NCRT_Chips"),
                       labels= c("2022","2022","2022","2022","2023","2023","2023","2023")),
         nrate = factor(env, levels = c("2022_MN_NCurve_Ch_0","2022_MN_NCurve_Ch_33","2022_MN_NCurve_Ch_66","2022_MN_NCRT_Chips",
                                        "2023_MN_Ncurve_Ch_0","2023_MN_Ncurve_Ch_33", "2023_MN_Ncurve_Ch_66","2023_MN_NCRT_Chips"),
                        labels= c("0","33","66","100","0","33","66","100"))) %>%
  #mutate(rank_round=rank(roundness), rank_lw=rank(lxw)) %>%
  ggplot( aes(x = nrate, y = tpn, colour = year ))+
  geom_boxplot() +
  #stat_smooth(method = lm,level = 0.95, colour = "black",linetype = 2)+
  labs(x = "Nlevel %", y = "Total N content")+
  theme(axis.text = element_text(face = "bold",size = 12), axis.title = element_text(face = "bold", size=12))#+
#annotate("text", x = 60, y = 79, label = "r^2==-0.82", parse = TRUE)
ggsave("./output/chspectrancurvedist.png",height=10.5, width=18.5, units="in", dpi=300)


ch_canopy50_blues %>%
  mutate(year = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                       "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                       labels= c("2022","2022","2022","2022","2023","2023","2023","2023")),
         nrate = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                        "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                        labels= c("0","33","66","100","0","33","66","100"))) %>%
  #mutate(rank_round=rank(roundness), rank_lw=rank(lxw)) %>%
  ggplot( aes(x = nrate, y = ntuber, colour = year ))+
  geom_point() +
  #stat_smooth(method = lm,level = 0.95, colour = "black",linetype = 2)+
  labs(x = "Nlevel %", y = "N concentration (tuber)")+
  theme(axis.text = element_text(face = "bold",size = 12), axis.title = element_text(face = "bold", size=12))#+

ggsave("./output/fmspectrancurvedistr.png",height=10.5, width=18.5, units="in", dpi=300)
