require(pacman)
p_load(tidyverse,evd,MASS)


### Clear memory
rm(list = ls())

n = 1000
m = 100

means <- c()
medians <- c()
trues <- c()
wtpspace <- c()


for (i in 1:m) {
  
  
  
  df <- data.frame(ID = rep(1:n,each=3),                     
                   best=NA,                                
                   beta_cost = -2,
                   beta_x = rep(rnorm(n,mean=6,sd=1)),
                   x=rep(runif(n,min=0,max=3)),       
                   cost=rep(runif(n,min=0,max=3)),   
                   nu = rep(rgumbel(n*3,loc=0,scale=5)))
  
  truewtp <- mean(-df$beta_x/df$beta_cost)
  
  df$u <- df$beta_cost*df$cost + df$beta_x*df$x + df$nu
  
  df <- df %>% group_by(ID) %>% mutate(maxu = max(u)) %>% ungroup()
  df$best <- ifelse(df$maxu == df$u,1,0)
  summary(res <- clogit(best ~ x + cost + strata(ID),data=df))
  
  getwtp <- mvrnorm(n=10000,mu=res$coefficients,Sigma=res$var) %>% as.data.frame()
  getwtp$wtp <- -getwtp$x/getwtp$cost
  #cat("proportion negative",sum(getwtp$wtp<0)/10000,"\n")
  #getwtp$wtp <- ifelse(getwtp$wtp<0,0,getwtp$wtp)
  
  means <- c(means,mean(getwtp$wtp))
  medians <- c(medians,median(getwtp$wtp))
  trues <- c(trues,truewtp)
  
  
  
  ######################################################################
  ######################################################################
  ############# FIX BETA_COST = 1 AND WORK IN WTP-SPACE ################
  ######################################################################
  ######################################################################
  
  
  
  
  
  ### Initialise code
 # apollo_initialise()
  
  ### Set core controls
  apollo_control = list(
    modelName  ="Apollo_example_1",
    modelDescr ="Simple MNL model on mode choice RP data",
    indivID    ="ID",
    panelData  = F
  )
  
  # ################################################################# #
  #### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
  # ################################################################# #
  
  
  df2 <- df
  df2$choice <- rep(1:3,n)
  df2 <- df2 %>%  pivot_wider(id_cols="ID",names_from = "choice",values_from=c("x","cost","best"))
  database <- df2 %>% mutate(best = 1*best_1 + 2*best_2 + 3*best_3) %>% dplyr::select(-best_1,-best_2,-best_3)
  
  
  # ################################################################# #
  #### DEFINE MODEL PARAMETERS                                     ####
  # ################################################################# #
  
  ### Vector of parameters, including any that are kept fixed in estimation
  apollo_beta=c(b_x=3,
                b_cost=-1)
  
  ### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
  apollo_fixed = c("b_cost")
  
  # ################################################################# #
  #### GROUP AND VALIDATE INPUTS                                   ####
  # ################################################################# #
  
  apollo_inputs = apollo_validateInputs(silent=T)
  
  # ################################################################# #
  #### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
  # ################################################################# #
  
  apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
    
    ### Attach inputs and detach after function exit
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    
    ### Create list of probabilities P
    P = list()
    
    ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
    V = list()
    V[['one']]  =  b_cost * cost_1 + b_x * x_1
    V[['two']]  = b_cost * cost_2 + b_x * x_2
    V[['three']]  = b_cost * cost_3 + b_x * x_3
    
    ### Define settings for MNL model component
    mnl_settings = list(
      alternatives  = c(one=1, two=2, three=3), 
      avail         = list(one=1,two=1,three=1), 
      choiceVar     = best,
      V             = V
    )
    
    ### Compute probabilities using MNL model
    P[['model']] = apollo_mnl(mnl_settings, functionality)
    
    ### Take product across observation for same individual
    # P = apollo_panelProd(P, apollo_inputs, functionality)
    
    ### Prepare and return outputs of function
    P = apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
  
  # ################################################################# #
  #### MODEL ESTIMATION                                            ####
  # ################################################################# #
  
  model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,estimate_settings = list(silent=T))
  
  # ################################################################# #
  #### MODEL OUTPUTS                                               ####
  # ################################################################# #
  
  wtpspace <- c(wtpspace,model$estimate[1])
  
  
}

summary(means)
summary(medians)
summary(trues)
summary(wtpspace)

simres <- data.frame(
  stat=rep(c("mean_WTP","median_WTP","true_WTP"),each=m),
  value=c(means,medians,trues)
)

ggplot(simres[simres$stat %in% c("median_WTP","mean_WTP"),]) +
  geom_density(aes(x=value,group=stat,fill=stat),alpha=.5) +
  lims(x=c(-1,7))
