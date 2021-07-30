# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #

### Clear memory
rm(list = ls())

### Load Apollo library
library(apollo)

n = 10000
J <- 2 # not currently used

# just do a mixed logit with correlated random component, 2 choices per person

database <- data.frame(ID = rep(1:n,each=2),
                       choice=NA,
                       SP = rep(0:1,n),
                       RP = rep(c(1,0),n),
                       w1 = rep(rnorm(n),each=2),
                       w2 = 0,
                       x3=rnorm(2*n),
                       cost3=rnorm(2*n),
                       x4=0,
                       cost4=0,
                       m = rep(rnorm(n),each=2),
                       epsRP = rep(rlogis(n,location=0,scale=1),each=2),
                       epsSP = rep(rlogis(n,location=0,scale=1),each=2))
database$av_2 <- 1
database$av_3 <- 1
database$av_4 <- 1

#database$w1 <- database$m + rep(rnorm(n,mean=0,sd=.25),each=2)


# parameters
beta1 <- -2 # effect for cost (choice model)
beta20 <- 4  # common component of effect for x
beta21 <- 2  # m-component of effect for x

alpha0 <- 1 # constant term, increase to boost the simulated response rate
alpha1 <- 1 # coef for w in selection equation
alpha2 <- 2 # coef for m in selection equation

database$U12 <-  alpha0 + alpha1*database$w1 + alpha2*database$m + database$epsRP
database$U34 <- beta1*(database$cost3) + (beta20 + beta21*database$m)*(database$x3) + database$epsSP


database$choice1 <- ifelse(database$RP==1 & database$U12 >= 0,1,
                           ifelse(database$RP==1 & database$U12 < 0, 2,NA))
database$choice2 <- ifelse(database$SP == 1 & database$U34 >= 0,3,
                           ifelse(database$SP==1 & database$U34 < 0,4,NA))

database$choice <- ifelse(database$RP==1,database$choice1,database$choice2)

# code to remove SP for nonresponders
database <- database %>% group_by(ID) %>% mutate(respond = sum(choice==1)) %>% ungroup()
database$choice <- ifelse(database$SP==1 & database$respond==0,NA,database$choice)
database <- database %>% dplyr::filter(!is.na(choice))
database <- database[!database$choice==1,]

# apollo doesn't like tibbles (http://www.apollochoicemodelling.com/forum/viewtopic.php?f=13&p=546)
database <- as.data.frame(database)

#######################################
#### COMPARE TO UNCORRECTED CLOGIT ####
#######################################

{
  df <- database %>% dplyr::filter(SP==1)
  df$choice <- ifelse(df$choice==3,1,0)
  df2 <- df
  df2$x3 <- 0
  df2$cost3 <- 0
  df2$choice <- (df2$choice - 1)^2
  df <- rbind(df,df2)
  df <- df[order(df$ID),]
  
  summary(uncor <- clogit(choice ~ x3 + cost3 + strata(ID),data=df))
  cat(paste0("uncorrected WTP is ",round(-1*uncor$coefficients[["x3"]]/uncor$coefficients[["cost3"]],3)))
}

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName  ="Apollo_example_4",
  modelDescr ="Two-level NL model with socio-demographics on mode choice SP data",
  indivID    ="ID"
)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

# database = read.csv("~/selection-sims/apollo_modeChoiceData.csv",header=TRUE)
# 
# ### Use only SP data
# database = subset(database,database$SP==1)
# 
# ### Create new variable with average income
# database$mean_income = mean(database$income)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta=c(b_w = 0.5,
              b_x = 0.5,
              b_cost = -0.2,
              lambda_nr  = 0.11)

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c()

### Read in starting values for at least some parameters from existing model output file
#apollo_beta=apollo_readBeta(apollo_beta,apollo_fixed,"Apollo_example_3",overwriteFixed=FALSE)

# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### List of utilities: these must use the same names as in nl_settings, order is irrelevant
  V = list()
  V[['alt2']] = -1*b_w*w1
  V[['alt3']] = b_x*x3  + b_cost*cost3
  V[['alt4']] = b_x*x4  + b_cost*cost4
  
  ### Specify nests for NL model
  nlNests      = list(root=1, respond=lambda_nr)
  
  ### Specify tree structure for NL model
  nlStructure= list()
  nlStructure[["root"]]   = c("alt2","respond")
  nlStructure[["respond"]]     = c("alt3","alt4")
  
  ### Define settings for NL model
  nl_settings <- list(
    alternatives = c(alt2=2, alt3=3, alt4=4),
    avail        = list(alt2=av_2, alt3=av_3, alt4=av_4),
    choiceVar    = choice,
    V            = V,
    nlNests      = nlNests,
    nlStructure  = nlStructure
  )
  
  ### Compute probabilities using NL model
  P[["model"]] = apollo_nl(nl_settings, functionality)
  
  ### Take product across observation for same individual
 # P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #

model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO SCREEN)                               ----
# ----------------------------------------------------------------- #

apollo_modelOutput(model)

-1*model$estimate[["b_x"]]/model$estimate[["b_cost"]]
cat(paste0("uncorrected WTP is ",round(-1*uncor$coefficients[["x3"]]/uncor$coefficients[["cost3"]],3)))

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO FILE, using model name)               ----
# ----------------------------------------------------------------- #

apollo_saveOutput(model)

# ################################################################# #
##### POST-PROCESSING                                            ####
# ################################################################# #

### Print outputs of additional diagnostics to new output file (remember to close file writing when complete)
sink(paste(model$apollo_control$modelName,"_additional_output.txt",sep=""),split=TRUE)

# ----------------------------------------------------------------- #
#---- LR TEST AGAINST MNL MODEL                                  ----
# ----------------------------------------------------------------- #

apollo_lrTest("Apollo_example_3", "Apollo_example_4")
apollo_lrTest("Apollo_example_3", model)

# ----------------------------------------------------------------- #
#---- switch off writing to file                                 ----
# ----------------------------------------------------------------- #

if(sink.number()>0) sink()