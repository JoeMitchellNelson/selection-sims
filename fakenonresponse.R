require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo,haven,readstata13)

df <- read.csv("C:/Users/Joe/Dropbox (University of Oregon)/Mitchell-Nelson-Cameron/Selection Paper/real_data/CAT_demo_nonvary_1-46_vary_47-56.csv")
df <- read.dta13("C:/Users/Joe/Dropbox (University of Oregon)/Mitchell-Nelson-Cameron/Selection Paper/real_data/CAT_stacked_selection_choice_data.dta")

vlabs <- varlabel(df)

selectvars <- c("selectionconst", "unempr20206", "ageb18to24", 
                "ageb25to34", "ageb75up", "pvoteothOR16", 
                "inc175to200", "inc200up", "y2020m3deaths", 
                "totfire20mi2020", "y2021m6cases", "pzipracehwpi", 
                "hour22", "totfire50mi1019", "pzipiinfor") 


df <- df %>% dplyr::filter(complete==1)
slabs <- selectvars

clabs <- c("cost", "emissions", "cbnjobspct", "grnjobspct", "auction", "equip", "workers", "regs", "statquo")
# "write" expressions to copy/paste
# cat(paste0("a_",slabs,"= 0,\n"))
# cat(paste0("a_",slabs,"*",slabs," + \n"))



# reshape data for apollo

df2 <- df[,c(3,which(names(df) %in% slabs))] %>% group_by(choice) %>%  summarise(across(everything(), list(diff)))
df3 <- df[,c(3,which(names(df) %in% clabs))] %>% group_by(choice) %>%  summarise(across(everything(), list(function (x) {-1*diff(x)} )))

df4 <- left_join(df3,df2)

names(df4) <- names(df4) %>% str_remove_all("_1$")

df5 <- df %>% select(choice,alt,best,fakenonresponse) %>% filter(alt==1) %>% select(-alt) %>% unique()

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

df10$choice2 <- ifelse(df10$fakenonresponse==1 & df10$choice2==2,1,df10$choice2)

df10$av1 <- ifelse(df10$choice2 %in% 1:2,1,0)
df10$av2 <- ifelse(df10$choice2 %in% 1:2,1,0)
df10$av3 <- ifelse(df10$choice2 %in% 3:4,1,0)
df10$av4 <- ifelse(df10$choice2 %in% 3:4,1,0)

df10 <- df10 %>% dplyr::filter(!(SP==1 & fakenonresponse==1))

repnames <- names(df10)

database <- df10 %>% arrange(ResponseNum) %>% ungroup() %>% as.matrix()
database <- database %>% as.data.frame()

database$totfire50mi1019 <- database$totfire50mi1019/1000

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
  
  mu_cost = -0.0078,
  sigma_cost = 0,
  b_emissions = -1.5,
  b_regs = 0.47,
  b_cbnjobspct = 4,
  b_grnjobspct = 3.4,
  b_auction = -0.1,
  b_equip = -0.25,
  b_workers = 0.52,
  mu_statquo = 0.27,
  sigma_statquo = 0.6,
  
  # selection variables
  
  mu_alpha = 0,
  sigma_alpha = 0.5,
  a_unempr20206 = 0,
  a_ageb18to24 = 0,
  a_ageb25to34 = 0,
  a_ageb75up = 0,
  a_pvoteothOR16 = 0,
  a_inc175to200 = 0,
  a_inc200up = 0,
  a_y2020m3deaths = 0,
  a_totfire20mi2020 = 0,
  a_y2021m6cases = 0,
  a_pzipracehwpi = 0,
  a_hour22 = 0,
  a_totfire50mi1019 = 0,
  a_pzipiinfor = 0,
  
  # rho and scale
  
  # rho = -0.36,
  sigma_alpha_statquo = 0,
  sigma_alpha_cost = 0,
  sigma_statquo_cost = 0,
  scale_RP = 1,
  scale_SP = 1)

apollo_beta <- startpars
#apollo_beta <- apollo_beta[apollo_beta!="rho"]

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
  interNormDraws = c("draws_alpha","draws_statquo","draws_cost"),
  intraDrawsType = "halton",
  intraNDraws    = 0,
  intraUnifDraws = c(),
  intraNormDraws = c()
)

### Create random parameters
apollo_randCoeff = function(apollo_beta, apollo_inputs){
  randcoeff = list()
  
  
  randcoeff[["alpha"]] = mu_alpha + sigma_alpha * draws_alpha
  randcoeff[["b_statquo"]] = mu_statquo + sigma_statquo * draws_statquo + sigma_alpha_statquo * draws_alpha
  randcoeff[["b_cost"]] = mu_cost + sigma_cost * draws_cost + sigma_alpha_cost * draws_alpha + sigma_statquo_cost * draws_statquo
  
  #randcoeff[["b_statquo"]] = mu_statquo + sqrt(1-rho^2)*(sigma_statquo) * draws_statquo + rho * sigma_statquo * draws_alpha
  
  #randcoeff[["b_statquo"]] = mu_statquo + sigma_statquo * draws_statquo
  
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
    a_unempr20206*unempr20206 + 
    a_ageb18to24*ageb18to24 + 
    a_ageb25to34*ageb25to34 + 
    a_ageb75up*ageb75up + 
    a_pvoteothOR16*pvoteothOR16 + 
    a_inc175to200*inc175to200 + 
    a_inc200up*inc200up + 
    a_y2020m3deaths*y2020m3deaths + 
    a_totfire20mi2020*totfire20mi2020 + 
    a_y2021m6cases*y2021m6cases + 
    a_pzipracehwpi*pzipracehwpi + 
    a_hour22*hour22 + 
    a_totfire50mi1019*totfire50mi1019 + 
    a_pzipiinfor*pzipiinfor
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
                        #    hessianRoutine="none",
                          silent=F))


# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO SCREEN)                               ----
# ----------------------------------------------------------------- #

apollo_modelOutput(model)

startpars <- model$estimate


summary(vanilla <- stats::glm(best ~ cost + emissions+ cbnjobspct+ grnjobspct+ auction +equip+ workers+ regs+ statquo + 0,
                              data=df[which(df$fakenonresponse==0),],family = binomial()))



summary(chocolate <- stats::glm(best ~ cost + emissions+ cbnjobspct+ grnjobspct+ auction +equip+ workers+ regs+ statquo + 0,
                                data=df,family = binomial()))

vanilla$coefficients["emissions"]/vanilla$coefficients["cost"] # with selection bias
chocolate$coefficients["emissions"]/chocolate$coefficients["cost"] # unbiased
model$estimate["b_emissions"]/model$estimate["mu_cost"] # corrected
model$estimate
