library(tidyverse)
library(plotly)
library(GGally)
library(StageWise)

#load FM data
load("./data/fmNAllSpectraits.rdata")

fm <- fm %>%
  #mutate_if(is.character, as.numeric) %>%
  mutate_at(vars(23:35), as.numeric) %>%
  mutate(nlevel = factor(nlevel, levels = c("0","33","66", "100"), labels = c("5.6","62.88","125.76","196.54"))) %>% # replace nlevels with applied n
  mutate(vinedm=(subvinedry/subvinefresh), tuberdm= (subtuberdry/subtuberfresh)) %>% #DM estimation for both vine and tuber
  mutate(ncontVine=vinedm*nvine, ncontTub=tuberdm*ntuber) %>% # N concentraion estimation for both vine and tuber
  mutate(tpn=ncontTub+ncontVine) %>% # total plant N estimation
  mutate(nup=tpn/as.numeric(nlevel)) %>% # N uptake estimation
  mutate(ute=tuberdm/tpn) %>% # N utilization estimation
  mutate(nue=ute*nup) %>% # N use efficiency estimation
  dplyr::select(-vinefresh, -subvinefresh, -subvinedry, -tuberfresh, -subtuberfresh, -subtuberdry) # remove columns no longer needed

# view data distribution
par(mfrow=c(3,6))
for(i in 10:ncol(fm)){
  hist(fm[,i], xlab = '', col='steelblue',  main=names(fm[i]))
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


#data cleaning and processing of FM NCurve
ntraits <- c(names(fm)[c(26:37)]) #
#for (i in (traits)) {
fm1 <- fm %>%
  filter(spec_date=="22_cover50" | spec_date=="23_cover50") %>%
  #fm1 <- data.frame(lapply(fm1[,c(10:34)], sdfx))
  mutate(across(everything(), ~ifelse(is.infinite(.),NA, .))) %>%
  mutate(across(all_of(ntraits),sdfx)) %>%
  mutate(redness=ifelse(redness<0 | redness>33, NA,redness),
         lightness=ifelse(lightness<0 | lightness>100, NA,lightness))

for(i in 10:ncol(fm1)){
  meanXi <- mean(fm1[,i],na.rm=TRUE)
  # naive imputation
  fm1[,i] <- ifelse(is.na(fm1[,i]),meanXi, fm1[,i])
}


# view data distribution after cleaning
par(mfrow=c(3,6))
Hmisc::hist.data.frame(fm1[,10:37])

#format data for generating blues and heritaility
fm1 %<>%
  rename(env=trial) %>%
  dplyr::select(env, id, plot,range, row, rep, nlevel, 10:37) %>%
  write_csv(., "./data/fmcanopy50.csv", quote = "none")


#########obtain blues
#library(StageWise)

effects <- data.frame(name=c("row","range","rep", "nlevel"),
                      fixed=c(FALSE,FALSE,TRUE,TRUE),
                      factor=c(TRUE,TRUE,TRUE,TRUE))
effects
traits <- c(names(fm1)[c(8:35)])
mod <- vector("list", length(traits))

for (i in traits) {
  mod[[i]] <- StageWise::Stage1(filename="./data/fmcanopy50.csv", traits=i, solver="spats", spline=c("row","range"))
}

# effects <- data.frame(name=c("rep","nlevel"),
#                       fixed=c(FALSE,TRUE),
#                       factor=c(TRUE,TRUE))
# for (i in traits) {
#   mod[[i]] <- StageWise::Stage1(filename="./data/fmcanopy50.csv", traits=i, solver="asreml", effects=effects)
# }
#yield=mod$yield$vcov$202
saveRDS(mod,"./data/fmcanopy50allTraitsBlues.rds")


#get blues, heritability, vcov matrix and residuals
mod <- readRDS("./data/fmcanopy50allTraitsBlues.rds")
fm_canopy50_blues <- data.frame(env=mod$red$blues$env, id=mod$red$blues$id,  red=mod$red$blues$BLUE,green=mod$green$blues$BLUE,
                                blue=mod$blue$blues$BLUE, nir=mod$nir$blues$BLUE, red_edge=mod$red_edge$blues$BLUE, ndvi=mod$ndvi$blues$BLUE,
                                ndre=mod$ndre$blues$BLUE, gndvi=mod$gndvi$blues$BLUE, gli=mod$gli$blues$BLUE, cire=mod$cire$blues$BLUE,
                                cig=mod$cig$blues$BLUE, ndwi=mod$ndwi$blues$BLUE, tcari.osavi=mod$tcari.osavi$blues$BLUE,nvine=mod$nvine$blues$BLUE,
                                ntuber=mod$ntuber$blues$BLUE,vinedm=mod$vinedm$blues$BLUE, tuberdm=mod$tuberdm$blues$BLUE, ncontVine=mod$ncontVine$blues$BLUE,
                                nconTuber=mod$ncontTub$blues$BLUE,tpn=mod$tpn$blues$BLUE,nup=mod$nup$blues$BLUE, ute=mod$ute$blues$BLUE, nue=mod$nue$blues$BLUE)


f50_byield <-  data.frame(env=mod$yield$blues$env, id=mod$yield$blues$id, yield=mod$yield$blues$BLUE)
f50_bred <-  data.frame(env=mod$redness$blues$env,id=mod$redness$blues$id,redness=mod$redness$blues$BLUE)
f50_blight <-  data.frame(env=mod$lightness$blues$env,id=mod$lightness$blues$id,lightness=mod$lightness$blues$BLUE)
f50_round <-  data.frame(env=mod$roundness$blues$env,id=mod$roundness$blues$id,roundness=mod$roundness$blues$BLUE)
f50_lxw <-  data.frame(env=mod$lxw$blues$env,id=mod$lxw$blues$id,lxw=mod$lxw$blues$BLUE)
data1 <- f50_round %>%
  left_join(f50_lxw, by=c("env","id")) %>%
  left_join(f50_blight, by=c("env","id")) %>%
  left_join(f50_bred, by=c("env","id")) %>%
  left_join(f50_byield, by=c("env","id"))

save(data1, file = "./data/traitBLUEsfmYRLRL.rdata")

fm_canopy50_H2 <- data.frame(env=mod$red$fit$env, yield=mod$yield$fit$H2, redness=mod$redness$fit$H2,lightness=mod$lightness$fit$H2,
                             roundness=mod$roundness$fit$H2,lxw=mod$lxw$fit$H2,nvine=mod$nvine$fit$H2,ntuber=mod$ntuber$fit$H2,
                             vinedm=mod$vinedm$fit$H2, tuberdm=mod$tuberdm$fit$H2, ncontVine=mod$ncontVine$fit$H2,
                             nconTuber=mod$ncontTub$fit$H2,tpn=mod$tpn$fit$H2,nup=mod$nup$fit$H2,ute=mod$ute$fit$H2,nue=mod$nue$fit$H2,
                             red=mod$red$fit$H2, green=mod$green$fit$H2,blue=mod$blue$fit$H2,
                             nir=mod$nir$fit$H2, red_edge=mod$red_edge$fit$H2, ndvi=mod$ndvi$fit$H2,ndre=mod$ndre$fit$H2,gndvi=mod$gndvi$fit$H2,
                             gli=mod$gli$fit$H2, cire=mod$cire$fit$H2, cig=mod$cig$fit$H2,ndwi=mod$ndwi$fit$H2,tcari.osavi=mod$tcari.osavi$fit$H2)




#fm_canopy50_vcov <- list(yield=mod$yield$vcov,redness=mod$redness$vcov, lightness=mod$lightness$vcov, roundness=mod$roundness$vcov
#                         ,lxw=mod$lxw$vcov)

#####saving traits data separately
save(fm_canopy50_blues, file = "./data/spectcanopy50BLUEsfm.rdata")
save(fm_canopy50_vcov, file = "./data/traitVcovfmRLRL.rdata")



###########heritability
h2_fc50 <- fm_canopy50_H2 %>%
  gather(key = "traits", value = "h2", -c(env)) %>%
  mutate(gs=paste0("canopy cover 50%"))




#########################
########################## Canopy cover peak flowering
fm2 <-  fm %>%
  #  mutate(yield=ifelse(spec_date=="22_cover50",yield22$yield, yield)) %>% #correct yield data
  filter(spec_date=="22_coverflower" | spec_date=="23_coverflower") %>%
  mutate(across(everything(), ~ifelse(is.infinite(.),NA, .))) %>%
  mutate(across(all_of(ntraits),sdfx)) %>%
  mutate(redness=ifelse(redness<0 | redness>33, NA,redness),
         lightness=ifelse(lightness<0 | lightness>100, NA,lightness))

for(i in 10:ncol(fm2)){
  meanXi <- mean(fm2[,i],na.rm=TRUE)
  # naive imputation
  fm2[,i] <- ifelse(is.na(fm2[,i]),meanXi, fm2[,i])
}


par(mfrow=c(3,5))
Hmisc::hist.data.frame(fm2[,10:37])


fm2 %<>%
  rename(env=trial) %>%
  dplyr::select(env, id, plot,range, row, rep, nlevel, 10:37) %>%
  write_csv(fm2, "./data/fmcanopycfl.csv", quote = "none")

#########obtain blues
effects <- data.frame(name=c("row","range","rep", "nlevel"),
                      fixed=c(FALSE,FALSE,FALSE,TRUE),
                      factor=c(TRUE,TRUE,TRUE,TRUE))
effects
traits <- c(names(fm2)[c(8:20)])
mod <- vector("list", length(traits))
#heritability.c1 <- data.frame(env=factor(), H2=integer(), trait=character())
#blues.c1 <- chips1.blues[,1:5]
#mod <- list(14)

for (i in traits) {
  mod[[i]] <- StageWise::Stage1(filename="./data/fmcanopycfl.csv", traits=i, solver="spats", spline=c("row","range"),effects=effects[3,])
}

saveRDS(mod,"./data/fmcanopycoverflowerTraitsBlues.rds")


#get blues, heritability, vcov matrix and residuals
mod <- readRDS("./data/fmcanopycoverflowerTraitsBlues.rds")
fm_canopycfl_blues <- data.frame(env=mod$red$blues$env, id=mod$red$blues$id,  red=mod$red$blues$BLUE,green=mod$green$blues$BLUE,
                                 blue=mod$blue$blues$BLUE, nir=mod$nir$blues$BLUE, red_edge=mod$red_edge$blues$BLUE, ndvi=mod$ndvi$blues$BLUE,
                                 ndre=mod$ndre$blues$BLUE, gndvi=mod$gndvi$blues$BLUE, gli=mod$gli$blues$BLUE, cire=mod$cire$blues$BLUE,
                                 cig=mod$cig$blues$BLUE, ndwi=mod$ndwi$blues$BLUE,tcari.osavi=mod$tcari.osavi$blues$BLUE)


fm_canopycfl_H2 <- data.frame(env=mod$red$fit$env,red=mod$red$fit$H2, green=mod$green$fit$H2,blue=mod$blue$fit$H2,
                              nir=mod$nir$fit$H2, red_edge=mod$red_edge$fit$H2, ndvi=mod$ndvi$fit$H2,ndre=mod$ndre$fit$H2,gndvi=mod$gndvi$fit$H2,
                              gli=mod$gli$fit$H2, cire=mod$cire$fit$H2, cig=mod$cig$fit$H2,ndwi=mod$ndwi$fit$H2,tcari.osavi=mod$tcari.osavi$fit$H2)

#fm_canopycfl_vcov <- list(yield=mod$yield$vcov,redness=mod$redness$vcov, lightness=mod$lightness$vcov, roundness=mod$roundness$vcov
#                         ,lxw=mod$lxw$vcov)

#####saving traits data separately
save(fm_canopycfl_blues, file = "./data/spectcanopycflBLUEsfm.rdata")
#save(fm_canopy50_vcov, file = "./data/traitVcovfmRLRL.rdata")


###########heritability
h2_fcfl <- fm_canopycfl_H2 %>%
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
fm3 <-  fm %>%
  filter(spec_date=="22_peakFlower" | spec_date=="23_peakFlower") %>%
  mutate(across(everything(), ~ifelse(is.infinite(.),NA, .))) %>%
  mutate(across(all_of(ntraits),sdfx)) %>%
  mutate(redness=ifelse(redness<0 | redness>33, NA,redness),
         lightness=ifelse(lightness<0 | lightness>100, NA,lightness))

for(i in 10:ncol(fm3))
{
  meanXi <- mean(fm3[,i],na.rm=TRUE)
  # naive imputation
  fm3[,i] <- ifelse(is.na(fm3[,i]),meanXi, fm3[,i])
}

par(mfrow=c(3,5))
for(i in 7:19){
  hist(fm3[,i], xlab = '', col='steelblue',  main=names(fm3[i]))
}



##saving info for the chips canopy50 spectra data for downstream analysis
fm3 %<>%
  dplyr::rename(env=trial) %>%
  dplyr::select(env, id, plot,range, row, rep, nlevel, 10:34) %>%
  write_csv(., "./data/fmpeakFlower.csv", quote = "none")



#########obtain blues
effects <- data.frame(name=c("row","range","rep"),
                      fixed=c(FALSE,FALSE,FALSE),
                      factor=c(TRUE,TRUE,TRUE))
effects

traits <- c(names(fm3)[8:20])
mod <- vector("list", length(traits))
#heritability.c1 <- data.frame(env=factor(), H2=integer(), trait=character())
#blues.c1 <- chips1.blues[,1:5]
#mod <- list(14)

for (i in traits) {
  mod[[i]] <- Stage1(filename="./data/fmpeakFlower.csv", traits=i, solver="spats", spline=c("row","range"),effects=effects[3,])
}

#yield=mod$yield$vcov$202
saveRDS(mod3,"./data/peakFlowerspecBluesfm.rds")


mod3 <- readRDS("./data/peakFlowerspecBluesfm.rds")


## get blues, heritability, vcov matrix and residuals

fm_peakFlower_blues <- data.frame(env=mod$red$blues$env, id=mod$red$blues$id,  red=mod$red$blues$BLUE,green=mod$green$blues$BLUE,
                                  blue=mod$blue$blues$BLUE, nir=mod$nir$blues$BLUE, red_edge=mod$red_edge$blues$BLUE, ndvi=mod$ndvi$blues$BLUE,
                                  ndre=mod$ndre$blues$BLUE, gndvi=mod$gndvi$blues$BLUE, gli=mod$gli$blues$BLUE, cire=mod$cire$blues$BLUE,
                                  cig=mod$cig$blues$BLUE, ndwi=mod$ndwi$blues$BLUE, tcari.osavi=mod$tcari.osavi$blues$BLUE)

fm_peakFlower_H2 <- data.frame(env=mod$red$fit$env,red=mod$red$fit$H2, green=mod$green$fit$H2,blue=mod$blue$fit$H2,
                               nir=mod$nir$fit$H2, red_edge=mod$red_edge$fit$H2, ndvi=mod$ndvi$fit$H2,ndre=mod$ndre$fit$H2,gndvi=mod$gndvi$fit$H2,
                               gli=mod$gli$fit$H2, cire=mod$cire$fit$H2, cig=mod$cig$fit$H2,ndwi=mod$ndwi$fit$H2,tcari.osavi=mod$tcari.osavi$fit$H2)
# fm_canopy50_vcov <- as.matrix(env=mod3$red$vcov,  red=mod3$red$vcov, green=mod3$green$vcov, blue=mod3$blue$vcov,
#                                  nir=mod3$nir$vcov, red_edge=mod3$red_edge$vcov, ndvi=mod3$ndvi$vcov, ndre=mod3$ndre$vcov,gndvi=mod3$gndvi$vcov,
#                                  gli=mod3$gli$vcov, cire=mod3$cire$vcov, cig=mod3$cig$vcov, ndwi=mod3$ndwi$vcov,tcari.osavi=mod3$tcari.osavi$vcov)

##Save peak fl blues

save(fm_peakFlower_blues, file = "./data/fm_specPeakflBLUEs.rdata")
#save(chips_canopy50_vcov, file = "./data/traitVcovchipsPeakfl.rdata")



##plot heritability for the bands and traits

h2_fpF <- fm_peakFlower_H2 %>%
  gather(key = "traits", value = "h2", -c(env)) %>%
  mutate(gs=paste0("peak flowering"),)

##update heritability for all growth stages
#saveRDS(rbind(h2_ccflower,h2_cpF), "chipsH2traits.rds")




#########################
########################## Vine biomass
fm4 <-  fm %>%
  #  mutate(yield=ifelse(spec_date=="22_cover50",yield22$yield, yield)) %>% #correct yield data
  filter(spec_date=="22_vineBiomass" | spec_date=="23_vineBiomass") %>%
  mutate(across(everything(), ~ifelse(is.infinite(.),NA, .))) %>%
  mutate(across(all_of(ntraits),sdfx)) %>%
  mutate(redness=ifelse(redness<0 | redness>33, NA,redness),
         lightness=ifelse(lightness<0 | lightness>100, NA,lightness))

for(i in 10:ncol(fm4))
{
  meanXi <- mean(fm4[,i],na.rm=TRUE)
  # naive imputation
  fm4[,i] <- ifelse(is.na(fm4[,i]),meanXi, fm4[,i])
}
par(mfrow=c(3,5))
for(i in 8:20){
  hist(fm4[,i], xlab = '', col='steelblue',  main=names(fm4[i]))
}



## saving info for the chips canopy50 spectra data for downstream analysis
fm4 %<>%
  dplyr::rename(env=trial) %>%
  dplyr::select(env, id, plot,range, row, rep, nlevel, 10:37) %>%
  write_csv(., "./data/fmvineBiomass.csv", quote = "none")


traits <- c(names(fm4)[8:20])
mod <- vector("list", length(traits))
#heritability.c1 <- data.frame(env=factor(), H2=integer(), trait=character())
#blues.c1 <- chips1.blues[,1:5]
#mod <- list(14)

for (i in traits) {
  mod[[i]] <- Stage1(filename="./data/fmvineBiomass.csv", traits=i, solver="spats", spline=c("row","range"),effects=effects[3,])
}

#yield=mod$yield$vcov$202
saveRDS(mod4,"./data/vineBiomassaTraitsBluesfm.rds")




mod4 <- readRDS("./data/vineBiomassTraitsBluesfm.rds")


## get blues, heritability, vcov matrix and residuals

fm_vineBiomass_blues <- data.frame(env=mod$red$blues$env, id=mod$red$blues$id,  red=mod$red$blues$BLUE,green=mod$green$blues$BLUE,
                                   blue=mod$blue$blues$BLUE, nir=mod$nir$blues$BLUE, red_edge=mod$red_edge$blues$BLUE, ndvi=mod$ndvi$blues$BLUE,
                                   ndre=mod$ndre$blues$BLUE, gndvi=mod$gndvi$blues$BLUE, gli=mod$gli$blues$BLUE, cire=mod$cire$blues$BLUE,
                                   cig=mod$cig$blues$BLUE, ndwi=mod$ndwi$blues$BLUE, tcari.osavi=mod$tcari.osavi$blues$BLUE)


fm_vineBiomass_H2 <- data.frame(env=mod$red$fit$env,red=mod$red$fit$H2, green=mod$green$fit$H2,blue=mod$blue$fit$H2,
                                nir=mod$nir$fit$H2, red_edge=mod$red_edge$fit$H2, ndvi=mod$ndvi$fit$H2,ndre=mod$ndre$fit$H2,gndvi=mod$gndvi$fit$H2,
                                gli=mod$gli$fit$H2, cire=mod$cire$fit$H2, cig=mod$cig$fit$H2,ndwi=mod$ndwi$fit$H2,tcari.osavi=mod$tcari.osavi$fit$H2)
# chips_vineBiomass_vcov <- as.matrix(env=mod4$red$vcov,  red=mod4$red$vcov, green=mod4$green$vcov, blue=mod4$blue$vcov,
#                                 nir=mod4$nir$vcov, red_edge=mod4$red_edge$vcov, ndvi=mod4$ndvi$vcov, ndre=mod4$ndre$vcov,gndvi=mod4$gndvi$vcov,
#                                 gli=mod4$gli$vcov, cire=mod4$cire$vcov, cig=mod4$cig$vcov, ndwi=mod4$ndwi$vcov,tcari.osavi=mod4$tcari.osavi$vcov)

##Save Vine biomass blues

save(fm_vineBiomass_blues, file = "./data/fm_spectvinebioBLUEs.rdata")


##plot heritability for the bands and traits
h2_fvB <- fm_vineBiomass_H2 %>%
  gather(key = "traits", value = "h2", -c(env)) %>%
  mutate(gs=paste0("vine biomass"),)
#saveRDS(rbind(h2_ccflower,h2_cpF,h2_cvB), "./data/chipsH2traits.rds")




#########################
########################## Senescence 1
fm5 <-  fm %>%
  #  mutate(yield=ifelse(spec_date=="22_cover50",yield22$yield, yield)) %>% #correct yield data
  filter(spec_date=="22_senesence1" | spec_date=="23_senesence1") %>%
  mutate(across(everything(), ~ifelse(is.infinite(.),NA, .))) %>%
  mutate(across(all_of(ntraits),sdfx)) %>%
  mutate(redness=ifelse(redness<0 | redness>33, NA,redness),
         lightness=ifelse(lightness<0 | lightness>100, NA,lightness))

for(i in 10:ncol(fm5)){
  meanXi <- mean(fm5[,i],na.rm=TRUE)
  # naive imputation
  fm5[,i] <- ifelse(is.na(fm5[,i]),meanXi, fm5[,i])
}

par(mfrow=c(3,5))
for(i in 7:20){
  hist(fm5[,i], xlab = '', col='steelblue',  main=names(fm5[i]))
}


fm5 %<>%
  dplyr::rename(env=trial) %>%
  dplyr::select(env, id, plot,range, row, rep, nlevel, 10:37) %>%
  write_csv(., "./data/fmsenescence1.csv", quote = "none")

effects <- data.frame(name=c("row","range","rep", "nlevel"),
                      fixed=c(FALSE,FALSE,TRUE,TRUE),
                      factor=c(TRUE,TRUE,TRUE,TRUE))
effects

traits <- c(names(fm5)[8:20])
mod <- vector("list", length(traits))
#heritability.c1 <- data.frame(env=factor(), H2=integer(), trait=character())
#blues.c1 <- chips1.blues[,1:5]
#mod <- list(14)

for (i in traits) {
  mod[[i]] <- Stage1(filename="./data/fmsenescence1.csv", traits=i, solver="spats", spline=c("row","range"),effects=effects[3,])
}

#yield=mod$yield$vcov$202
saveRDS(mod,"./data/senescence1TraitsBluesfm.rds")




mod5 <- readRDS("./data/senescence1TraitsBluesfm.rds")


## get blues, heritability, vcov matrix and residuals

fm_senescence1_blues <- data.frame(env=mod$red$blues$env, id=mod$red$blues$id,red=mod$red$blues$BLUE,green=mod$green$blues$BLUE,
                                   blue=mod$blue$blues$BLUE, nir=mod$nir$blues$BLUE, red_edge=mod$red_edge$blues$BLUE, ndvi=mod$ndvi$blues$BLUE,
                                   ndre=mod$ndre$blues$BLUE, gndvi=mod$gndvi$blues$BLUE, gli=mod$gli$blues$BLUE, cire=mod$cire$blues$BLUE,
                                   cig=mod$cig$blues$BLUE, ndwi=mod$ndwi$blues$BLUE, tcari.osavi=mod$tcari.osavi$blues$BLUE)


fm_senescence1_H2 <- data.frame(env=mod$red$fit$env,red=mod$red$fit$H2, green=mod$green$fit$H2,blue=mod$blue$fit$H2,
                                nir=mod$nir$fit$H2, red_edge=mod$red_edge$fit$H2, ndvi=mod$ndvi$fit$H2,ndre=mod$ndre$fit$H2,gndvi=mod$gndvi$fit$H2,
                                gli=mod$gli$fit$H2, cire=mod$cire$fit$H2, cig=mod$cig$fit$H2,ndwi=mod$ndwi$fit$H2,tcari.osavi=mod$tcari.osavi$fit$H2)
# chips_canopy50_vcov <- as.matrix(env=mod5$red$vcov,  red=mod5$red$vcov, green=mod5$green$vcov, blue=mod5$blue$vcov,
#                                 nir=mod5$nir$vcov, red_edge=mod5$red_edge$vcov, ndvi=mod5$ndvi$vcov, ndre=mod5$ndre$vcov,gndvi=mod5$gndvi$vcov,
#                                 gli=mod5$gli$vcov, cire=mod5$cire$vcov, cig=mod5$cig$vcov, ndwi=mod5$ndwi$vcov,tcari.osavi=mod5$tcari.osavi$vcov)


##Save Senesce1 blues

save(fm_senescence1_blues, file = "./data/fm_spectSenesce1BLUEs.rdata")


# plot heritability for the bands and traits

h2_fS1 <- fm_senescence1_H2 %>%
  gather(key = "traits", value = "h2", -c(env)) %>%
  mutate(gs=paste0("senescence1"),)




#########################
########################## Senescence 2
fm6 <-  fm %>%
  #  mutate(yield=ifelse(spec_date=="22_cover50",yield22$yield, yield)) %>% #correct yield data
  filter(spec_date=="22_senesence2" | spec_date=="23_senesence2") %>%
  mutate(across(everything(), ~ifelse(is.infinite(.),NA, .))) %>%
  mutate(across(all_of(ntraits),sdfx)) %>%
  mutate(redness=ifelse(redness<0 | redness>33, NA,redness),
         lightness=ifelse(lightness<0 | lightness>100, NA,lightness))

for(i in 10:ncol(fm6)){
  meanXi <- mean(fm1[,i],na.rm=TRUE)
  # naive imputation
  fm6[,i] <- ifelse(is.na(fm6[,i]),meanXi, fm6[,i])
}

par(mfrow=c(3,5))
for(i in 7:19){
  hist(fm6[,i], xlab = '', col='steelblue',  main=names(fm6[i]))
}


fm6 %<>%
  dplyr::rename(env=trial) %>%
  dplyr::select(env, id, plot,range, row, rep, nlevel, 10:34) %>%
  write_csv(., "./data/fmsenesence2.csv", quote = "none")


traits <- c(names(fm6)[8:32])
mod <- vector("list", length(traits))
#heritability.c1 <- data.frame(env=factor(), H2=integer(), trait=character())
#blues.c1 <- chips1.blues[,1:5]
#mod <- list(14)

for (i in traits) {
  mod[[i]] <- StageWise::Stage1(filename="./data/fmsenesence2.csv", traits=i, solver="spats", spline=c("row","range"),effects=effects[3,])
}

#yield=mod$yield$vcov$202
saveRDS(mod,"../data/senesence2allTraitsBluesfm.rds")


mod6 <- readRDS("../data/senesence2allTraitsBluesfm.rds")


#get blues, heritability, vcov matrix and residuals

fm_senesence2_blues <- data.frame(env=mod$red$blues$env, id=mod$red$blues$id,  red=mod$red$blues$BLUE,green=mod$green$blues$BLUE,
                                  blue=mod$blue$blues$BLUE, nir=mod$nir$blues$BLUE, red_edge=mod$red_edge$blues$BLUE, ndvi=mod$ndvi$blues$BLUE,
                                  ndre=mod$ndre$blues$BLUE, gndvi=mod$gndvi$blues$BLUE, gli=mod$gli$blues$BLUE, cire=mod$cire$blues$BLUE,
                                  cig=mod$cig$blues$BLUE, ndwi=mod$ndwi$blues$BLUE, tcari.osavi=mod$tcari.osavi$blues$BLUE,nvine=mod$nvine$blues$BLUE,
                                  ntuber=mod$ntuber$blues$BLUE,vinedm=mod$vinedm$blues$BLUE, tuberdm=mod$tuberdm$blues$BLUE, ncontVine=mod$ncontVine$blues$BLUE,
                                  nconTuber=mod$ncontTub$blues$BLUE,nup=mod$nup$blues$BLUE)

fm_senesence2_H2 <- data.frame(env=mod$red$fit$env,nvine=mod$nvine$fit$H2,ntuber=mod$ntuber$fit$H2,
                               vinedm=mod$vinedm$fit$H2, tuberdm=mod$tuberdm$fit$H2, ncontVine=mod$ncontVine$fit$H2,
                               nconTuber=mod$ncontTub$fit$H2,nup=mod$nup$fit$H2,red=mod$red$fit$H2, green=mod$green$fit$H2,blue=mod$blue$fit$H2,
                               nir=mod$nir$fit$H2, red_edge=mod$red_edge$fit$H2, ndvi=mod$ndvi$fit$H2,ndre=mod$ndre$fit$H2,gndvi=mod$gndvi$fit$H2,
                               gli=mod$gli$fit$H2, cire=mod$cire$fit$H2, cig=mod$cig$fit$H2,ndwi=mod$ndwi$fit$H2,tcari.osavi=mod$tcari.osavi$fit$H2)
# chips_canopy50_vcov <- as.matrix(env=mod6$red$vcov,  red=mod6$red$vcov, green=mod6$green$vcov, blue=mod6$blue$vcov,
#                                 nir=mod6$nir$vcov, red_edge=mod6$red_edge$vcov, ndvi=mod6$ndvi$vcov, ndre=mod6$ndre$vcov,gndvi=mod6$gndvi$vcov,
#                                 gli=mod6$gli$vcov, cire=mod6$cire$vcov, cig=mod6$cig$vcov, ndwi=mod6$ndwi$vcov,tcari.osavi=mod6$tcari.osavi$vcov)


#Save Senesce2 blues

save(fm_senesence2_blues, file = "./data/fm_spectSenesce2BLUEs.rdata")


#plot heritability for the bands and traits

h2_fS2 <- fm_senesence2_H2 %>%
  gather(key = "traits", value = "h2", -c(env)) %>%
  mutate(gs=paste0("senesence2"),)



#####################################
#################################################################################

##plot all Broad sense heritability
#update merging
saveRDS(rbind(h2_fc50, h2_fcfl,h2_fpF,h2_fvB,h2_fS1,h2_fS2), "./data/fmH2traits.rds")
saveRDS(rbind(fm_H2, h2_fS2),"../data/fmH2traits.rds")
fm_H2 <- readRDS("./data/fmH2traits.rds")
fm_H2$class <- paste("fresh market class")




############plot h2
fh2 <- fm_H2[-(1:40),] %>%
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

msncdata1.5 <- cbind(fm_canopy50_blues[,c(1:2,16:25)] %>%
  unite(env.id,c(env,id),sep = ":"),
fm_canopy50_blues[,3:15]%>%
  rename_with( ~ paste0(.x, ".1")),
fm_canopycfl_blues[,3:15]%>%
  #unite(env.id,c(env,id),sep = ":") %>%
  rename_with( ~ paste0(.x, ".2")),
fm_peakFlower_blues[,3:15]%>%
  #unite(env.id,c(env,id),sep = ":") %>%
  rename_with( ~ paste0(.x, ".3")),
fm_vineBiomass_blues[,3:15]%>%
  #unite(env.id,c(env,id),sep = ":") %>%
  rename_with( ~ paste0(.x, ".4")),
fm_senescence1_blues[,3:15]%>%
  #unite(env.id,c(env,id),sep = ":") %>%
  rename_with( ~ paste0(.x, ".5"))
)
save(msncdata1.5, file = "./data/msspecFMdatablues1.5.rda")

spectra <- msncdata1.5[,]
cortab <- data.frame(spectra=NA, r2=NA, year=NA, nrate=NA, gs=NA)

# Plot distribution for total plant content
d1 <- fm_canopy50_blues %>%
  mutate(year = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                        "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                        labels= c("2022","2022","2022","2022","2023","2023","2023","2023")),
         nrate = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                        "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                        labels= c("0","33","66","100","0","33","66","100"))) %>%
  #mutate(rank_round=rank(roundness), rank_lw=rank(lxw)) %>%
  ggplot( aes(x = nrate, y = tpn, colour = year ))+
  geom_boxplot() +
  #stat_smooth(method = lm,level = 0.95, colour = "black",linetype = 2)+
  labs(x = "Nlevel %", y = "Total N content")+
  theme(axis.text = element_text(face = "bold",size = 12), axis.title = element_text(face = "bold", size=12))#+
  #annotate("text", x = 60, y = 79, label = "r^2==-0.82", parse = TRUE)
ggsave("./output/fmspectrancurvedist.png",height=10.5, width=18.5, units="in", dpi=300)


fm_canopy50_blues %>%
  mutate(year = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                       "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                       labels= c("2022","2022","2022","2022","2023","2023","2023","2023")),
         nrate = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                        "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                        labels= c("0","33","66","100","0","33","66","100"))) %>%
  #mutate(rank_round=rank(roundness), rank_lw=rank(lxw)) %>%
  ggplot( aes(x = nrate, y = ncontVine, colour = year ))+
  geom_boxplot() +
  #stat_smooth(method = lm,level = 0.95, colour = "black",linetype = 2)+
  labs(x = "Nlevel %", y = "N concentration (tuber)")+
  theme(axis.text = element_text(face = "bold",size = 12), axis.title = element_text(face = "bold", size=12))#+

ggsave("./output/fmspectrancurvedistr.png",height=10.5, width=18.5, units="in", dpi=300)

p1 <- fm_canopy50_blues %>%
  mutate(year = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                       "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                       labels= c("2022","2022","2022","2022","2023","2023","2023","2023")),
         nrate = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                        "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                        labels= c("0","33","66","100","0","33","66","100"))) %>%
              ggpairs(columns = 3:25,
             mapping = ggplot2::aes(color = year ),
             upper = list(continuous = wrap("cor", size = 3))
             )

p2 <- fm_canopycfl_blues %>%
  mutate(year = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                       "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                       labels= c("2022","2022","2022","2022","2023","2023","2023","2023")),
         nrate = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                        "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                        labels= c("0","33","66","100","0","33","66","100"))) %>%
  ggpairs(columns = 3:22,
          mapping = ggplot2::aes(color = year ),
          upper = list(continuous = wrap("cor", size = 3))
  )

p3 <- fm_peakFlower_blues %>%
  mutate(year = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                       "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                       labels= c("2022","2022","2022","2022","2023","2023","2023","2023")),
         nrate = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                        "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                        labels= c("0","33","66","100","0","33","66","100"))) %>%
  ggpairs(columns = 3:22,
          mapping = ggplot2::aes(color = year ),
          upper = list(continuous = wrap("cor", size = 3))
  )
p4 <- fm_vineBiomass_blues %>%
  mutate(year = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                       "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                       labels= c("2022","2022","2022","2022","2023","2023","2023","2023")),
         nrate = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                        "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                        labels= c("0","33","66","100","0","33","66","100"))) %>%
  ggpairs(columns = 3:22,
          mapping = ggplot2::aes(color = year ),
          upper = list(continuous = wrap("cor", size = 3))
  )
p5 <- fm_senescence1_blues %>%
  mutate(year = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                       "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                       labels= c("2022","2022","2022","2022","2023","2023","2023","2023")),
         nrate = factor(env, levels = c("2022_MN_NCurve_FM_0","2022_MN_NCurve_FM_33","2022_MN_NCurve_FM_66","2022_MN_NCRT_FM",
                                        "2023_MN_NCurve_FM_0","2023_MN_NCurve_FM_33", "2023_MN_NCurve_FM_66","2023_MN_NCRT_FM"),
                        labels= c("0","33","66","100","0","33","66","100"))) %>%
  ggpairs(columns = 3:22,
          mapping = ggplot2::aes(color = year ),
          upper = list(continuous = wrap("cor", size = 3))
  )

p <- ggpubr::ggarrange(p1, p2, p3, p4, p5, ncol=3, nrow=2)

ggsave("./output/fmspectrancurvecorall.png",height=10.5, width=18.5, units="in", dpi=300)
