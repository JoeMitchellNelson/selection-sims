require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo)

n = 500
J <- 2 # not currently used

# just do a mixed logit with correlated random component, 2 choices per person

database <- data.frame(ID = rep(1:n,each=2),
                       choice=NA,
                       x1=runif(2*n,min=0,max=3),
                       cost1=runif(2*n,min=0,max=3),
                       w1 = rnorm(2*n),
                       x2=0,
                       cost2=0,
                       w2 = 0,
                       eps = rlogis(2*n,location=0,scale=1))

ms <- mvrnorm(n=n,mu=c(0,0),Sigma = matrix(c(1,.5,
                                             .5,1),byrow = T,nrow = 2))
database$m1 <- rep(ms[,1],each=2)
database$m2 <- rep(ms[,2],each=2)

head(database)

# parameters
beta1 <- -2 # effect for cost (choice model)
beta20 <- 4  # effect for x
betaw <- 1


database$deltaU <- beta1*(database$cost1) + 
  (beta20 + database$m1)*(database$x1) + 
  (betaw + database$m2)*(database$w1) + 
  database$eps
database$choice <- ifelse(database$deltaU > 0, 1, 2)


# apollo doesn't like tibbles (http://www.apollochoicemodelling.com/forum/viewtopic.php?f=13&p=546)
database <- as.data.frame(database)
head(database)

database <- database[-1,]

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName ="Mixed logit test",
  modelDescr ="Mixed logit model, correlated normal for x and w",
  indivID   ="ID",  
  mixing    = TRUE, 
  nCores    = 2
)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta=c(mu_x = 4,
              sigma_x = 1,
              b_cost = -2,
              mu_w = 1,
              sigma_w = 1,
              sigma_x_w = .5)

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c()

# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws = list(
  interDrawsType = "halton",
  interNDraws    = 500,
  interUnifDraws = c(),
  interNormDraws = c("draws_x","draws_w"),
  intraDrawsType = "halton",
  intraNDraws    = 0,
  intraUnifDraws = c(),
  intraNormDraws = c()
)

### Create random parameters
apollo_randCoeff = function(apollo_beta, apollo_inputs){
  randcoeff = list()
  
  randcoeff[["b_x"]] = mu_x + sigma_x * draws_x 
  randcoeff[["b_w"]] = mu_w + sigma_w * draws_w + sigma_x_w * draws_x
  
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
  V = list()
  V[['alt1']] = b_x*x1  + b_cost*cost1 + b_w*w1
  V[['alt2']] = b_x*x2  + b_cost*cost2 + b_w*w2
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives  = c(alt1=1, alt2=2),
    avail         = list(alt1=1, alt2=1),
    choiceVar     = choice,
    V             = V
  )
  
  ### Compute probabilities using MNL model
  P[['model']] = apollo_mnl(mnl_settings, functionality)
  
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
                        estimate_settings=list(hessianRoutine="maxLik"))

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO SCREEN)                               ----
# ----------------------------------------------------------------- #

apollo_modelOutput(model)
