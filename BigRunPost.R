require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo)

cres <- readRDS("~/selection-sims/bigsim/correctedmixl.Rds")
ures <- readRDS("~/selection-sims/bigsim/uncorrectedmixl.Rds")
apollo_modelOutput(cres)
apollo_modelOutput(ures)
db <- readRDS("~/selection-sims/bigsim/simfulldata.Rds")
head(db)

rp <- db[db$RP==1,]
sp <- db[db$SP==1,]
head(sp)

reg1 <- glm(I(choice==1) ~ w1,data=rp,family="binomial")
summary(reg1)
rp$fitted <- reg1$fitted.values %>% qlogis()
rp$fitted <- (rp$fitted - mean(rp$fitted))/sd(rp$fitted)

rp2 <- rp %>% dplyr::select(ID,fitted)
sp <- left_join(sp,rp2)
sp <- sp %>% dplyr::filter(!is.na(choice))
summary(sp$fitted)

reg2 <- glm(I(choice==3) ~ x3 + x3:I(fitted) + cost3, data=sp,family="binomial")
summary(reg2)

reg3 <- glm(I(choice2==3) ~ x3 + cost3, data=db,family="binomial")
summary(reg3)
