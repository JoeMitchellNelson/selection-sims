require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo)

df <- read.csv("C:/Users/Joe/Dropbox (University of Oregon)/Mitchell-Nelson-Cameron/Selection Paper/real_data/CAT_demo_nonvary_1-46_vary_47-56.csv")

summary(vanilla <- stats::glm(best ~ cost + emissions+ cbnjobspct+ grnjobspct+ auction +equip+ workers+ regs+ statquo + 0,
                              data=df[which(df$complete==1),],family = binomial()))

vanilla$coefficients["emissions"]/vanilla$coefficients["cost"]

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
df$choicenum <- c(rep(1:6,1000),rep(1,548))
# take choice "1" for each person (so everyone, responder or not, has one row)
dfrp <- df %>% filter(choicenum==1)
dfrp$choice <- dfrp$complete + 1
dfrp <- dfrp %>% mutate(av1=1,av2=1,av3=0,av4=0)
df <- df %>% mutate(av1=0,av2=0,av3=1,av4=1)
df <- df %>% filter(complete==1)

df$choice <- df$best + 3


database <- rbind(df,dfrp) %>% arrange(ResponseNum)
database$w <- database$pzippov150up_pzipinc35to50k
database$SP <- ifelse(database$choice %in% 3:4,1,0)
database$RP <- (database$SP - 1)^2

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName ="Selection correction",
  modelDescr ="Mixed logit model",
  indivID   ="ResponseNum",  
  mixing    = TRUE, 
  nCores    = 12
)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

## Vector of parameters, including any that are kept fixed in estimation
apollo_beta=c(
  b_emissions = unname(vanilla$coefficients["emissions"]),
  b_cost = unname(vanilla$coefficients["cost"]),
  b_regs = unname(vanilla$coefficients["regs"]),
  b_cbnjobspct = unname(vanilla$coefficients["cbnjobspct"]),
  b_grnjobspct = unname(vanilla$coefficients["grnjobspct"]),
  b_auction = unname(vanilla$coefficients["auction"]),
  b_equip = unname(vanilla$coefficients["equip"]),
  b_workers = unname(vanilla$coefficients["workers"]),
  mu_statquo = unname(vanilla$coefficients["statquo"]),
  sigma_statquo = 1,
  # b_ageb18to24_gendmale = -.48,
  b_ageb18to24_pzipraceasian = -0.55,
  # b_inclt20_Wed = -.77,
  # b_gendfemale_Fri = -.44,
  # b_gendfemale_raceblack = -.65,
  b_pzipiwholes_mobile = -.147,
  b_pzipraceblack_pziprural = 1.18,
  b_nearfiredis1019_pintntsatel = .000077,
  b_ageb55to64_Wed = -1.347,
  b_ageb65up_Wed = -1.025,
  mu_alpha = 2.44,
  sigma_alpha = 1,
 # sigma_emissions = 1,
  rho_statquo = -.68,
 # rho_emissions = 0,
  scale_RP = 1,
  scale_SP = 1)

 apollo_beta <- startpars

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
  randcoeff[["b_statquo"]] = mu_statquo + sqrt(1-rho_statquo^2)*(sigma_statquo) * draws_statquo + rho_statquo * sigma_statquo * draws_alpha
 # randcoeff[["b_emissions"]] = mu_emissions + sqrt(1-rho_emissions^2)*(sigma_emissions) * draws_emissions + rho_emissions * sigma_emissions * draws_alpha
  
  
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
    # b_ageb18to24_gendmale*ageb18to24_gendmale + 
    b_ageb18to24_pzipraceasian*ageb18to24_pzipraceasian + 
    #  b_inclt20_Wed*inclt20_Wed + 
    # b_gendfemale_Fri*gendfemale_Fri + 
    #  b_gendfemale_raceblack*gendfemale_raceblack + 
    b_pzipiwholes_mobile*pzipiwholes_mobile + 
    b_pzipraceblack_pziprural*pzipraceblack_pziprural + 
    b_nearfiredis1019_pintntsatel*nearfiredis1019_pintntsatel + 
    b_ageb55to64_Wed*ageb55to64_Wed + 
    b_ageb65up_Wed*ageb65up_Wed
  V[['alt3']] = 0
  V[['alt4']] = b_cost*cost + b_emissions*emissions + b_regs*regs + 
    b_cbnjobspct*cbnjobspct + b_grnjobspct*grnjobspct + b_auction*auction + b_equip*equip + b_workers*workers + 
    b_statquo*statquo
  
  
  mnl_settings = list(
    alternatives  = c(alt1=1, alt2=2, alt3=3, alt4=4),
    avail         = list(alt1=av1, alt2=av2,alt3=av3,alt4=av4),
    choiceVar     = choice,
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
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #

model = apollo_estimate(apollo_beta, apollo_fixed,
                        apollo_probabilities, apollo_inputs, 
                        estimate_settings=list(hessianRoutine="numDeriv",silent=T))


# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO SCREEN)                               ----
# ----------------------------------------------------------------- #

apollo_modelOutput(model)

startpars <- model$estimate
startpars["rho_emissions"] <- -.5
startpars["rho_statquo"] <- -.3
startpars
