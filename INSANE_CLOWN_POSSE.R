require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,data.table,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo,haven,readstata13,broom,stringr,heatmaply)

source("~/selection-sims/selection_helpers.R")


mat1 <- c("a1","b1","b2","c1","c2","c3","d1","d2","d3","d4")
data1 <- c("alpha","statusquo","ben","cost")

# adds a function called "rescale" to global environment 
build_rescale_function(matrixvars=mat1,datavars=data1)

# test rescaler function
do.call(rescaler,as.list(rnorm(length(mat1)+length(data1))))

#imr <- function (x) {dnorm(-1*x)/(1-pnorm(-1*x))}



icp_rn <- read.dta13("~/selection-sims/resp-nonresp-data-for-selection-one-line_JMN.dta")
icp_rn2 <- read.dta13("~/selection-sims/ICP_survey_invitation_characteristics.dta")
icp_rn <- left_join(icp_rn,icp_rn2,by="caseid")
icp <- read.dta13("~/selection-sims/ICP_simple_2-3alt_choice_data.dta")

icp_rn <- icp_rn %>% 
  mutate(statseason = ifelse(invite_student_fall,"StudentFall",
                             ifelse(invite_student_spring,"StudentSpring",
                                    ifelse(invite_employee_fall,"EmployeeFall","EmployeeSpring"))))



icp_rn$rmnds <- icp_rn$rmnd1 + icp_rn$rmnd2 + icp_rn$rmnd3

icp_rn <- icp_rn %>% mutate(RESPONDED=ifelse(caseid %in% icp$caseid,1,0)) %>% 
  select(RESPONDED,everything())

icp_rn$invite_employee <- ifelse(icp_rn$invite_employee_fall|icp_rn$invite_employee_spring,1,0)

icp_rn <- icp_rn %>% 
  mutate(newgroup = ifelse(groupnum %in% 1:2, "A",
                            ifelse(groupnum %in% 3:8, "B",
                                   ifelse(groupnum %in% c(9,11), "C",
                                          ifelse(groupnum %in% c(10,12), "D",
                                                 ifelse(groupnum %in% c(13,15),"E","F"))))))


groupdummies <- model.matrix( ~ newgroup -1, icp_rn) %>% as.data.frame()

icp_rn <- cbind(icp_rn,groupdummies)

# promisng response models
summary(glm(RESPONDED ~ invite_employee_spring + invite_employee_fall + invite_student_fall + f1 + hav_ZIPSCORES, data=icp_rn,family="binomial"))
summary(glm(RESPONDED ~ I(invite_employee_fall + invite_student_fall), data=icp_rn,family="binomial"))

summary(glm(RESPONDED ~ factor(newgroup), data=icp_rn,family="binomial"))
summary(lm(RESPONDED ~ newgroup, data=icp_rn,family="binomial"))


# main specification from original paper
clogitmodel <- clogit(best ~ ben + cost + airtrfees+ buildfees + taxpayers + spacademic + spoffsets + statusquo + strata(choice),data=icp)
clogitmodel <- clogit(best ~ ben + cost +  statusquo + strata(choice),data=icp)
summary(clogitmodel)
# reshape choice data for apollo

icp_reshape <- icp %>% select(-stdempfees,-spcarbon) %>% 
  pivot_wider(names_from=alt,
              values_from=c(best,ben , cost, airtrfees, buildfees, taxpayers, spacademic, spoffsets, statusquo))

icp_reshape <- icp_reshape %>% mutate(av1=ifelse(is.na(statusquo_1),0,1),
                                      av2=ifelse(is.na(statusquo_2),0,1),
                                      av3=ifelse(is.na(statusquo_3),0,1))

icp_reshape <- icp_reshape %>% select(choice_id = choice,everything()) %>% 
  mutate(choice = ifelse(best_1==1 & !is.na(best_1),1,
                         ifelse(best_2==1 & !is.na(best_2),2,3))) %>% 
  select(choice_id,choice,everything())

icp_rn <- icp_rn %>% select(caseid,RESPONDED,starts_with("newgroup"),f1,hav_ZIPSCORES) %>% 
  mutate(av4=1,av5=1,choice=ifelse(RESPONDED==1,4,5)) %>% 
  mutate(choice_id = row_number() + max(icp_reshape$choice_id)) %>% 
  select(-RESPONDED) %>% 
  select(choice_id,choice,caseid,everything())

icp <- bind_rows(icp_reshape,icp_rn) %>% arrange(caseid,choice_id)

icp <- icp %>% relocate(starts_with("newgroup"),.after="statusquo_2")

icp <- icp %>% replace_na(list(av1=0,av2=0,av3=0,av4=0,av5=0))


icp <- icp %>% group_by(caseid) %>% 
  mutate(f1=max(f1,na.rm=T),
         hav_ZIPSCORES=max(hav_ZIPSCORES,na.rm=T),
         across(starts_with("newgroup"), ~ max(.x,na.rm=T)))

icp$cost_1 <- ifelse(icp$av1 == 0 & icp$av2 == 1,icp$cost_2,icp$cost_1)
icp$cost_2 <- ifelse(icp$av2 == 0 & icp$av1 == 1,icp$cost_1,icp$cost_2)
icp$ben_1 <- ifelse(icp$av1 == 0 & icp$av2 == 1,icp$ben_2,icp$ben_1)
icp$ben_2 <- ifelse(icp$av2 == 0 & icp$av1 == 1,icp$ben_1,icp$ben_2)
icp$statusquo_1 <- ifelse(icp$av1 == 0 & icp$av2 == 1,icp$statusquo_2,icp$statusquo_1)
icp$statusquo_2 <- ifelse(icp$av2 == 0 & icp$av1 == 1,icp$statusquo_1,icp$statusquo_2)
# icp$airtrfees_1 <- ifelse(icp$av1 == 0,icp$airtrfees_2,icp$airtrfees_1)
# icp$airtrfees_2 <- ifelse(icp$av2 == 0,icp$airtrfees_1,icp$airtrfees_2)

icp[is.na(icp)] <- 0

icp <- icp %>% dplyr::rename(factor1=f1)

database <- as.data.frame(icp)


### Initialise code
apollo_initialise()

### Set core controls
### NOTE: YOU MAY NEED TO LOWER NCORES BEFORE RUNNING. MOST COMPUTERS MAX OUT AT 8, IN WHICH CASE YOU SHOULD USE NO MORE THAN 7.
### ALSO: "NCORES" SHOULD REALLY BE CALLED "NTHREADS". TYPICAL TO HAVE 2 THREADS PER CORE.

apollo_control <- list(
  modelName ="Selection correction",
  modelDescr ="Mixed logit model",
  indivID   ="caseid",  
  mixing    = TRUE,
  nCores    = 15,
  panelData = TRUE
)


apollo_beta <- c(
  #response variables
  
  # a
  alpha_B = 1.14,
  alpha_C = .76,
  alpha_D = .48,
  alpha_E = .92,
  alpha_F = .87,
  
  #alpha_factor1 = -0.009,
  #alpha_hav_ZIPSCORES = 0.069,
  mu_alpha = -3.17,
  
  # choice variables
  
  mu_statusquo = -0.42,
  mu_ben = 0.0147,
  mu_cost = -0.0058,  
  # mu_airtrfees = 7,
  # mu_buildfees = 0,
  # mu_taxpayers  = 0,
  # mu_spacademic = 0,
  # mu_spoffsets = 0,
  
  beta_imr = 0,
  
  a1 = 0,
  b1 = 0, b2 = 0,
  c1 = 0, c2 = 0, c3 = 0,
  d1 = 0, d2 = 0, d3 = 0, d4 = 0,
  # e1 = 0, e2 = 0, e3 = 0, e4 = 0, e5 = 0,
  # f1 = 0, f2 = 0, f3 = 0, f4 = 0, f5 = 0, f6 = 0, 
  # g1 = 0, g2 = 0, g3 = 0, g4 = 0, g5 = 0, g6 = 0, g7 = 0,
  # h1 = 0, h2 = 0, h3 = 0, h4 = 0, h5 = 0, h6 = 0, h7 = 0, h8 = 0,
  # i1 = 0, i2 = 0, i3 = 0, i4 = 0, i5 = 0, i6 = 0, i7 = 0, i8 = 0, i9 = 0,
  

  
  # scales
  
  scale_RP =  1,
  scale_SP = 1
)

zeropars <- c(rep(0,length(apollo_beta)-2),1,1)
names(zeropars) <- names(apollo_beta)

apollo_beta <- startpars %>% startparjitter()
apollo_beta <- zeropars

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed <- c("scale_RP","scale_SP")

# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws <- list(
  interDrawsType = "halton",
  interNDraws    = 500,
  interUnifDraws = c(),
  interNormDraws = c("draws_alpha",
                     "draws_statquo",
                     "draws_ben",
                     "draws_cost"
                     # "draws_airtrfees"
                     # "draws_buildfees",
                     # "draws_taxpayers",
                     # "draws_spacademic",
                     # "draws_spoffsets"
                     ),
  intraDrawsType = "halton",
  intraNDraws    = 0,
  intraUnifDraws = c(),
  intraNormDraws = c()
)




### Create random parameters
apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  
  randcoeff = list()
  
  randcoeff[["random_alpha"]]    = mu_alpha +     a1*draws_alpha
  
  randcoeff[["random_statquo"]]  = mu_statusquo + b1*draws_alpha + b2*draws_statquo
  
  randcoeff[["random_ben"]]      = mu_ben +       c1*draws_alpha + c2*draws_statquo + c3*draws_ben
    
  randcoeff[["random_cost"]] =    mu_cost +       d1*draws_alpha + d2*draws_statquo + d3*draws_ben + d4*draws_cost 
    
  # randcoeff[["random_airtrfees"]] = mu_airtrfees+ e1*draws_alpha + e2*draws_statquo + e3*draws_ben + e4*draws_cost + e5*draws_airtrfees
  # 
  # randcoeff[["random_buildfees"]] = mu_buildfees+ f1*draws_alpha + f2*draws_statquo + f3*draws_ben + f4*draws_cost + f5*draws_airtrfees + 
  #   f6*draws_buildfees 
  # 
  # randcoeff[["random_taxpayers"]] = mu_taxpayers+ g1*draws_alpha + g2*draws_statquo + g3*draws_ben + g4*draws_cost + g5*draws_airtrfees + 
  #   g6*draws_buildfees + g7*draws_taxpayers
  #   
  # randcoeff[["random_spacademic"]]=mu_spacademic +h1*draws_alpha + h2*draws_statquo + h3*draws_ben + h4*draws_cost + h5*draws_airtrfees + 
  #   h6*draws_buildfees + h7*draws_taxpayers + h8*draws_spacademic
  #   
  # randcoeff[["random_spoffsets"]] =mu_spoffsets + i1*draws_alpha + i2*draws_statquo + i3*draws_ben + i4*draws_cost + i5*draws_airtrfees + 
  #   i6*draws_buildfees + i7*draws_taxpayers + i8*draws_spacademic + i9*draws_spoffsets

  

  
  return(randcoeff)
}





# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs <- apollo_validateInputs(silent=F)


apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Function initialisation: do not change the following three commands
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  
  # 1 - alt 1
  # 2 - alt 2
  # 3 - alt 3 (status qup)
  # 4 - respond
  # 5 - don't respond (conceptually, the status quo for the response decision)
  
  
  
  V = list()
  V[['alt1']] = ( random_statquo*statusquo_1 + random_ben*ben_1 + random_cost*cost_1 + 
    beta_imr*imr((
            mu_alpha + 
            alpha_B * newgroupB +
            alpha_C * newgroupC +
            alpha_D * newgroupD +
            alpha_E * newgroupE +
            alpha_F * newgroupF
                  )/sqrt(a1^2+(pi^2/3))))/
    rescaler(  a1,
               b1,b2,
               c1,c2,c3,
               d1,d2,d3,d4,

              
               0,
               (statusquo_1 + statusquo_2)/2,
               (ben_1 + ben_2)/2,
               (cost_1 + cost_2)/2

               )
  
  V[['alt2']] = (random_statquo*statusquo_2 + random_ben*ben_2 + random_cost*cost_2 +  
    beta_imr*imr((
                    mu_alpha + 
                    alpha_B * newgroupB +
                    alpha_C * newgroupC +
                    alpha_D * newgroupD +
                    alpha_E * newgroupE +
                    alpha_F * newgroupF)/sqrt(a1^2+(pi^2/3)) ))/
    rescaler(  
               a1,
               b1,b2,
               c1,c2,c3,
               d1,d2,d3,d4,
               
               0,
               (statusquo_1 + statusquo_2)/2,
               (ben_1 + ben_2)/2,
               (cost_1 + cost_2)/2
    )
  
  V[['alt3']] = random_statquo*statusquo_3 + random_ben*ben_3 + random_cost*cost_3 
  
  
  V[['alt4']] = random_alpha + 
    alpha_B * newgroupB +
    alpha_C * newgroupC +
    alpha_D * newgroupD +
    alpha_E * newgroupE +
    alpha_F * newgroupF
    

  V[['alt5']] = 0
  
  
  
  mnl_settings = list(
    alternatives  = c(alt1=1, alt2=2, alt3=3, alt4=4,alt5=5),
    avail         = list(alt1=av1, alt2=av2,alt3=av3,alt4=av4,alt5=av5),
    choiceVar     = choice,
    V             = lapply(V,"*",scale_RP),
    rows          = (av4==1 | av5==1),
    componentName = "MNL-RP"
  )
  
  P[['RP']] = apollo_mnl(mnl_settings, functionality)
  
  ### Compute probabilities for the SP part of the data using MNL model
  mnl_settings$V = lapply(V, "*", scale_SP)
  mnl_settings$rows = (av1==1 | av2==1 | av3 ==1)
  mnl_settings$componentName = "MNL-SP"
  
  P[['SP']] = apollo_mnl(mnl_settings, functionality)
  
  ### Combined model
  P = apollo_combineModels(P, apollo_inputs, functionality)
  
  ### Compute probabilities using MNL model
  #P[["model"]] = apollo_mnl(mnl_settings, functionality)
  
  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ### Average across inter-individual draws
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  
  # use selectionpopweight as weights
  #P = apollo_weighting(P, apollo_inputs, functionality)
  
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #


model <- apollo_estimate(apollo_beta, apollo_fixed,
                         apollo_probabilities, apollo_inputs,
                         estimate_settings=list(
                          hessianRoutine="numDeriv",
                          # hessianRoutine="none",
                           silent=F))

apollo_modelOutput(model)

# 7523.375 saddle point
# 7538.718 saddle point

#startpars <- model$estimate

mat <- un_cholesky(model$estimate[11:20])

mat <- mat

drawwtp <- mvrnorm(n=100000,mu=unname(model$estimate[6:9]),Sigma=mat)
quantile(-1*drawwtp[,3]/drawwtp[,4],c(0.025,0.5,0.975))

cov2cor(mat)

startpars <- model$estimate
startpars


