require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo,haven,readstata13)

df <- read.csv("C:/Users/Joe/Dropbox (University of Oregon)/Mitchell-Nelson-Cameron/Selection Paper/real_data/CAT_demo_nonvary_1-46_vary_47-56.csv")
df <- read.dta13("C:/Users/Joe/Dropbox (University of Oregon)/Mitchell-Nelson-Cameron/Selection Paper/real_data/CAT_stacked_selection_choice_data.dta")

vlabs <- varlabel(df)

selectlabs <-  c("dm:ZCTA pr:Renter: heat none",
                 "dm:1=Own age:18-24",
                 "dm:1=Own age:25-34",
                 "dm:1=Own age:75 and up",
                 "dm:1=Own hhld inc:125-150K",
                 "dm:ZCTA pr:Owner: heat none",
                 "dm:1=Used Windows OS",
                 "dm:1=Own hhld inc:lt 20K",
                 "dm:1=Race:Amer. Indian",
                 "dm:1=Race:Asian",
                 "dm:ZCTA pr:Age 55-59",
                 "dm:1=Race:Other",
                 "dm:District pr: vote Other party 2016",
                 "dm:ZCTA pr:Renter: heat coal",
                 "dm:1=Start: hour end. at 2:00",
                 "dm:ZCTA pr:Age 15-19",
                 "dm:1=Start: hour end. at 9:00",
                 "dm:County COVID deaths/50K Jun'21",
                 "dm:County COVID deaths/50K Mar'20",
                 "dm:1=Started survey on Tue",
                 "dm:1=Used a mobile device",
                 "dm:ZCTA pr:Hawaiian/Pac.Isl.")


df <- read_dta("C:/Users/Joe/Dropbox (University of Oregon)/Mitchell-Nelson-Cameron/Selection Paper/real_data/CAT_stacked_selection_choice_data.dta")
df$complete <- ifelse(df$completepopweight==0,0,1)
slabs <- names(df)[which(vlabs %in% selectlabs)]

clabs <- c("cost", "emissions", "cbnjobspct", "grnjobspct", "auction", "equip", "workers", "regs", "statquo")
# "write" expressions to copy/paste
# cat(paste0("a_",slabs,"= 0,\n"))
# cat(paste0("a_",slabs,"*",slabs," + \n"))

summary(vanilla <- stats::glm(best ~ cost + emissions+ cbnjobspct+ grnjobspct+ auction +equip+ workers+ regs+ statquo + 0,
                       data=df[which(df$complete==1),],family = binomial()))

vanilla$coefficients["emissions"]/vanilla$coefficients["cost"]


# reshape data for apollo

df2 <- df[,c(3,which(names(df) %in% slabs))] %>% group_by(choice) %>%  summarise(across(everything(), list(diff)))
df3 <- df[,c(3,which(names(df) %in% clabs))] %>% group_by(choice) %>%  summarise(across(everything(), list(function (x) {-1*diff(x)} )))

df4 <- left_join(df3,df2)

names(df4) <- names(df4) %>% str_remove_all("_1$")

df5 <- df %>% select(choice,alt,best) %>% filter(alt==1) %>% select(-alt) %>% unique()

df6 <- left_join(df4,df5)

df7 <- df %>% select(choice,choicenum) %>% unique()

df7 <- df7 %>% mutate(SP=ifelse(choicenum==0,0,1),RP =ifelse(choicenum==0,1,0)) %>% select(-choicenum)

df8 <- left_join(df6,df7)

df9 <- df %>% select(choice,ResponseNum,selectpopweight) %>% unique()

df10 <- left_join(df8,df9) %>% as.data.frame()
# Notes to help translate Trudy code:

# SELECTIONVARIABLES "ageb55to64_Wed ageb18to24_pzipraceasian pzipraceblack_pziprural nearfiredis1019_pintntsatel ageb65up_Wed pzipiwholes_mobile"
# mixlogit best ${SELECTIONVARIABLES}, group(choice) id(ResponseNum)

# responders make 7 choices
# non-responders make 1 choice

# recode choice column:
# 1 - don't respond
# 2 - respond
# 3 - status quo
# 4 - take the policy

# number each choice task 1 thru 6 for responders, 1 for non-responders
#df$choicenum <- c(rep(1:6,1000),rep(1,548))
# take choice "1" for each person (so everyone, responder or not, has one row)

df10 <- df10 %>% mutate(choice2 = ifelse(SP==1 & best==0,3,
                                         ifelse(SP==1 & best==1,4,
                                                ifelse(SP==0 & best==1,2,1))))

df10$av1 <- ifelse(df10$choice2 %in% 1:2,1,0)
df10$av2 <- ifelse(df10$choice2 %in% 1:2,1,0)
df10$av3 <- ifelse(df10$choice2 %in% 3:4,1,0)
df10$av4 <- ifelse(df10$choice2 %in% 3:4,1,0)

repnames <- names(df10)

database <- df10 %>% arrange(ResponseNum) %>% ungroup() %>% as.matrix()
database <- database %>% as.data.frame()

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName ="Selection correction",
  modelDescr ="Mixed logit model",
  indivID   ="ResponseNum",  
  mixing    = TRUE, 
  nCores    = 12,
  weights = "selectpopweight"
)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

## Vector of parameters, including any that are kept fixed in estimation
apollo_beta=c(
              # choice variables
  
              b_emissions = -1.5,
              b_cost = -0.0078,
              b_regs = 0.47,
              b_cbnjobspct = 4,
              b_grnjobspct = 3.4,
              b_auction = -0.1,
              b_equip = -0.25,
              b_workers = 0.52,
              mu_statquo = 0.27,
              sigma_statquo = 0.6,
              
              # selection variables
              
              a_y2021m6deaths = 0.14,
              a_y2020m3deaths = 0.65,  
              a_pzipage1519= 0.1,
              a_pzipage5559= 0.155,
              a_pzipracehwpi= -0.23,
              a_pvoteothOR16= 0.8,
              a_pzipownhtnone= 0.4,
              a_pziprenthtcoal= -3.4,
              a_pziprenthtnone= 0.173,
              a_ageb18to24= -0.7,
              a_ageb25to34= -0.4,
              a_ageb75up= -1.5,
              a_raceasian= -0.75,
              a_raceamind= -1,
              a_raceother= 0.75,
              a_inclt20= -0.3,
              a_inc125to150= 0.42,
              a_mobile= -0.3,
              a_OnWindows= 0.7,
              a_Tue= -0.4,
              a_hour2= 1,
              a_hour9= 1.25,
              mu_alpha = 2.44,
              sigma_alpha = 1,
              
              # rho and scale
              
              rho = -0.36,
              scale_RP = 1,
              scale_SP = 1)

apollo_beta <- startpars[which(names(startpars)!="rho")]

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("scale_RP","scale_SP")

# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws = list(
  interDrawsType = "halton",
  interNDraws    = 500,
  interUnifDraws = c(),
  interNormDraws = c("draws_alpha","draws_statquo"),
  intraDrawsType = "halton",
  intraNDraws    = 0,
  intraUnifDraws = c(),
  intraNormDraws = c()
)

### Create random parameters
apollo_randCoeff = function(apollo_beta, apollo_inputs){
  randcoeff = list()
  
  
  randcoeff[["alpha"]] = mu_alpha + sigma_alpha * draws_alpha
  #randcoeff[["b_statquo"]] = mu_statquo + sqrt(1-rho^2)*(sigma_statquo) * draws_statquo + rho * sigma_statquo * draws_alpha
  
  randcoeff[["b_statquo"]] = mu_statquo + sigma_statquo * draws_statquo
  
  return(randcoeff)
}

# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Function initialisation: do not change the following three commands
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  
  # 1 - don't respond
  # 2 - respond
  # 3 - status quo
  # 4 - take the policy
  
  
  V = list()
  V[['alt1']] = 0
  V[['alt2']] = alpha + 
    a_y2021m6deaths*y2021m6deaths +
    a_y2020m3deaths*y2020m3deaths +
    a_pzipage1519*pzipage1519 + 
    a_pzipage5559*pzipage5559 + 
    a_pzipracehwpi*pzipracehwpi + 
    a_pvoteothOR16*pvoteothOR16 + 
    a_pzipownhtnone*pzipownhtnone + 
    a_pziprenthtcoal*pziprenthtcoal + 
    a_pziprenthtnone*pziprenthtnone + 
    a_ageb18to24*ageb18to24 + 
    a_ageb25to34*ageb25to34 + 
    a_ageb75up*ageb75up + 
    a_raceasian*raceasian + 
    a_raceamind*raceamind + 
    a_raceother*raceother + 
    a_inclt20*inclt20 + 
    a_inc125to150*inc125to150 + 
    a_mobile*mobile + 
    a_OnWindows*OnWindows + 
    a_Tue*Tue + 
    a_hour2*hour2 + 
    a_hour9*hour9
  V[['alt3']] = 0
  V[['alt4']] = b_cost*cost + b_emissions*emissions + b_regs*regs + 
    b_cbnjobspct*cbnjobspct + b_grnjobspct*grnjobspct + b_auction*auction + b_equip*equip + b_workers*workers + 
    b_statquo*statquo
  
  
  mnl_settings = list(
    alternatives  = c(alt1=1, alt2=2, alt3=3, alt4=4),
    avail         = list(alt1=av1, alt2=av2,alt3=av3,alt4=av4),
    choiceVar     = choice2,
    V             = lapply(V,"*",scale_RP),
    rows          = (RP==1),
    componentName = "MNL-RP"
  )
  
  P[['RP']] = apollo_mnl(mnl_settings, functionality)
  
  ### Compute probabilities for the SP part of the data using MNL model
  mnl_settings$V = lapply(V, "*", scale_SP)
  mnl_settings$rows = (SP==1)
  mnl_settings$componentName = "MNL-SP"
  
  P[['SP']] = apollo_mnl(mnl_settings, functionality)
  
  ### Combined model
  P = apollo_combineModels(P, apollo_inputs, functionality)
  
  
  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ### Average across inter-individual draws
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  
  # use selectionpopweight as weights
  P = apollo_weighting(P, apollo_inputs, functionality)
  
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #

model = apollo_estimate(apollo_beta, apollo_fixed,
                        apollo_probabilities, apollo_inputs, 
                        estimate_settings=list(
                           hessianRoutine="numDeriv",
                        #  hessianRoutine="none",
                          silent=F))


# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO SCREEN)                               ----
# ----------------------------------------------------------------- #

apollo_modelOutput(model)

startpars <- model$estimate

