require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,boot,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo)

source("~/selection-sims/helper_funcs.R")

n_sims <- 200

res <- matrix(data=NA,nrow=n_sims,ncol=28) %>% as.data.frame()

names(res) <- c("seed","LL","message",
                "b_w",   "mu_x",  "sigma_beta",  "b_cost", "mu_alpha", "sigma_alpha",   "rho", "scale_RP", "scale_SP" ,
                paste0(c("b_w",   "mu_x",  "sigma_beta",  "b_cost", "mu_alpha", "sigma_alpha",   "rho", "scale_RP", "scale_SP"),"_se"),
                "maxEigen","uncorrected_wtp","uncorrected_wtp_se", "fullsample_wtp","fullsample_wtp_se", "adhoc_wtp","adhoc_wtp_se")


for (i in 175:n_sims) {
  
  ### Clear memory
  
  rm(database)
  
  
  s <- Sys.time() %>% as.numeric()
  set.seed(s)
  
  n = 2000
  J = 3

  # just do a mixed logit with correlated random component, 2 choices per person
  
  
  
  database <- data.frame(ID = rep(1:n,each=2),                       # identifies individual
                         choice=NA,                                  # to be filled in with choices
                         SP = rep(0:1,n),                            # indicator for policy choice (stated preference)
                         RP = rep(c(1,0),n),                         # indicator for participation (revealed preference)
                         w1 = rep(rnorm(n),each=2),                  # indiv-specific var that affects participation
                         w2 = 0,
                         x3=rep(runif(n,min=0,max=3),each=2),        # x attribute for first policy option
                         cost3=rep(runif(n,min=0,max=3),each=2),     # cost for first policy option
                         x4=rep(runif(n,min=0,max=3),each=2),        # x att for second policy option
                         cost4=rep(runif(n,min=0,max=3),each=2),     # cost for second opion
                         x5=0,                                       # x = 0 for status quo
                         cost5=0,                                    # cost = 0 for status quo
                         m = rep(rnorm(n),each=2),                   # unobservable, affects participation and policy preference
                         #m2 = rep(rnorm(n),each=2),                  # unobservable, affects policy preference only
                         m2 = 0,
                         nu1 = rep(rgumbel(n,loc=0,scale=1),each=2), #
                         nu2 = rep(rgumbel(n,loc=0,scale=1),each=2),   #
                         nu3 = rep(rgumbel(n,loc=0,scale=1),each=2),     # type 1 extreme value error terms, one for each alternative
                         nu4 = rep(rgumbel(n,loc=0,scale=1),each=2),   #
                         nu5 = rep(rgumbel(n,loc=0,scale=1),each=2),  #
                         noise = rep(rnorm(n),each=2)
  )
  
  # indicator variables for availability
  # (1 and 2 available together for participation decision)
  # (3, 4, 5 available together for policy decision)
  database$av_1 <- database$RP * 1
  database$av_2 <- database$RP * 1
  database$av_3 <- database$SP * 1
  database$av_4 <- database$SP * 1
  database$av_5 <- database$SP * 1
  
  
  # parameters
  beta1 <- -2 # effect for cost (choice model)
  beta20 <- 4  # common component of effect for x
  beta21 <- -2  # m-component of effect for x
  beta22 <- 2 # heterogeneity, uncorrelated with response propensity
  
  alpha0 <- 0 # constant term, increase to boost the simulated response rate
  alpha1 <- 2 # coef for w in selection equation
  alpha2 <- 2 # coef for m in selection equation
  
  # calculate indirect utilities for respond (U1) and nonrespond (U2)
  database$U1 <-  alpha0 + alpha1*database$w1 + alpha2*database$m + database$nu1
  database$U2 <-  0 + alpha1*database$w2 + 0*database$m + database$nu2
  
  # calculate indirect utilities for each of three policy alternatives (U3, U4, U5)
  database$U3 <- beta1*(database$cost3) + (beta20 + beta21*database$m + beta22*database$m2)*(database$x3) + database$nu3
  database$U4 <- beta1*(database$cost4) + (beta20 + beta21*database$m + beta22*database$m2)*(database$x4) + database$nu4
  database$U5 <- beta1*(database$cost5) + (beta20 + beta21*database$m + beta22*database$m2)*(database$x5) + database$nu5
  
  
  
  database$choice1 <- ifelse(database$U1 >= database$U2,1,2)
  database$choice2 <- ifelse(database$U3 >= pmax(database$U3,database$U4,database$U5),3,
                             ifelse(database$U4 >= pmax(database$U3,database$U4,database$U5),4,5))
  
  database$choice <- ifelse(database$RP==1,database$choice1,database$choice2)
  
  
  database <- database %>% group_by(ID) %>% mutate(respond = sum(choice==1)) %>% ungroup()
  
  ######################################################
  ### clogits for comparison ###
  ######################################################
  
  SPs <- database[database$SP==1,] %>% dplyr::select(-c(SP,RP,av_1,av_2,av_3,av_4,av_5,choice))
  RPs <- database[database$RP==1,] %>% dplyr::select(-c(SP,RP,av_1,av_2,av_3,av_4,av_5,choice))
  flatdb <- left_join(SPs,RPs) %>% as.data.frame()
  
    df <- database %>% dplyr::filter(SP==1)
    df2 <- database %>% dplyr::filter(RP==1) %>% dplyr::select(ID,w1,m,respond)
    df2$noise <- rnorm(nrow(df2))
    temp <- data.frame(ID = rep(1:n,each=J),alt=1:J + 2)
    df <- left_join(temp,df)
    df$choice <- ifelse(df$choice==df$alt,1,0)
    df$x <- ifelse(df$alt==3,df$x3,
                   ifelse(df$alt==4,df$x4,df$x5))
    df$cost <- ifelse(df$alt==3,df$cost3,
                   ifelse(df$alt==4,df$cost4,df$cost5))
    df <- df %>% dplyr::select(ID,alt,choice,x,cost,respond)
    
    df2 <- left_join(df,df2,by=c("ID","respond"))
    summary(uncor <- clogit(choice ~ x + cost + strata(ID),data=df[df$respond==1,]))
    summary(cor <- clogit(choice ~ x*m + cost*m + strata(ID),data=df2))
    
    uncor1 <- mvrnorm(n=10000,mu=uncor$coefficients,Sigma=uncor$var) %>% as.data.frame()
    uncor1$wtp <- -1*uncor1$x/uncor1$cost
    uncor_wtp <- mean(uncor1$wtp)
    uncor_se <- sqrt(1/(nrow(uncor1)-1) * sum((uncor1$wtp - uncor_wtp)^2))
    
    cor1 <- mvrnorm(n=10000,mu=cor$coefficients,Sigma=cor$var) %>% as.data.frame()
    cor1$wtp <- -1*cor1$x/cor1$cost
    cor_wtp <- mean(cor1$wtp)
    cor_se <- sqrt(1/(nrow(cor1)-1) * sum((cor1$wtp - cor_wtp)^2))

    
   # summary(adhoc2 <- clogit(choice ~ x + cost + x:fittedRP + cost:fittedRP + strata(ID),data=df2[df2$respond==1,]))
    
    # bootstrapping with 1000 replications
    results <- boot(data=flatdb, statistic=bs,
                    R=500, formula=choice ~ x + cost + x:fittedRP + cost:fittedRP + strata(ID))
    
    # view results
    bsres <- summary(results)
    
  
  
  # code to remove SP for nonresponders
  database$choice <- ifelse(database$SP==1 & database$respond==0,NA,database$choice)
  database <- database %>% dplyr::filter(!is.na(choice))
  
  
  # apollo doesn't like tibbles (http://www.apollochoicemodelling.com/forum/viewtopic.php?f=13&p=546)
  database <- as.data.frame(database)


  
  ### Initialise code
  apollo_initialise()
  
  ### Set core controls
  apollo_control = list(
    modelName ="First stab at selection correction",
    modelDescr ="Mixed logit model",
    indivID   ="ID",  
    mixing    = TRUE, 
    nCores    = 12
  )
  
  # ################################################################# #
  #### DEFINE MODEL PARAMETERS                                     ####
  # ################################################################# #
  
  ### Vector of parameters, including any that are kept fixed in estimation
  apollo_beta=c(b_w = 2,
                mu_x = 2,
                sigma_beta = 1,
                b_cost = -1,
                mu_alpha = 0,
                sigma_alpha = 1,
                rho = .7,
                scale_RP = 1,
                scale_SP = 2)
  
  ### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
  apollo_fixed = c("scale_RP","b_cost")
  
  # ################################################################# #
  #### DEFINE RANDOM COMPONENTS                                    ####
  # ################################################################# #
  
  ### Set parameters for generating draws
  apollo_draws = list(
    interDrawsType = "halton",
    interNDraws    = 1000,
    interUnifDraws = c(),
    interNormDraws = c("draws_alpha","draws_x"),
    intraDrawsType = "halton",
    intraNDraws    = 0,
    intraUnifDraws = c(),
    intraNormDraws = c()
  )
  
  ### Create random parameters
  apollo_randCoeff = function(apollo_beta, apollo_inputs){
    randcoeff = list()
    
    
    randcoeff[["alpha"]] = mu_alpha + sigma_alpha * draws_alpha
    randcoeff[["b_x"]] = mu_x + sqrt(1-rho^2)*(sigma_beta) * draws_x + rho * sigma_beta * draws_alpha
    
    
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
    V[['alt1']] = alpha + b_w*w1
    V[['alt2']] = b_w*w2
    V[['alt3']] = b_x*x3  + b_cost*cost3
    V[['alt4']] = b_x*x4  + b_cost*cost4
    V[['alt5']] = b_x*x5  + b_cost*cost5
    
    
    
    mnl_settings = list(
      alternatives  = c(alt1=1, alt2=2, alt3=3, alt4=4,alt5=5),
      avail         = list(alt1=av_1, alt2=av_2,alt3=av_3,alt4=av_4,alt5=av_5),
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
  
  res$seed[i] <- s
  res$LL[i] <- model$maximum
  res$message[i] <- ifelse(model$message=="successful convergence ",1,0)
  res[i,4:12] <- model$estimate
  res[i,13:21] <- model$robse
  res$maxEigen[i] <- ifelse(!is.null(model$eigValue),model$eigValue,NA)
  res$uncorrected_wtp[i] <- uncor_wtp
  res$uncorrected_wtp_se[i] <- uncor_se
  res$fullsample_wtp[i] <- cor_wtp
  res$fullsample_wtp_se[i] <- cor_se
  res$adhoc_wtp[i] <- bsres$original
  res$adhoc_wtp_se[i] <- bsres$bootSE
  
  
  cat(paste0("\n\n\n\t\t\t",i,"\n\n\n"))
  
  cat(paste0(mean(-1*res$mu_x/res$b_cost,na.rm=T)),"\n\n")
  cat(paste0(mean(res$rho,na.rm=T)),"\n\n")
  
  
  
  
}

ggplot(res) +
  geom_density(aes(x=fullsample_wtp),color="black",alpha=0.5) +
  geom_density(aes(x=mu_x),fill="forestgreen",color=NA,alpha=0.5) +
  geom_vline(xintercept=mean(res$mu_x,na.rm=T),color="forestgreen",linetype="dashed") +
  geom_density(aes(x=uncorrected_wtp),fill="red",color=NA,alpha=0.5) +
  geom_density(aes(x=adhoc_wtp),fill="blue",color=NA,alpha=0.5) +
  labs(x="WTP Estimate",y="Density") +
  geom_vline(xintercept=2) +
  lims(x=c(1,3))
  

resnd <- res %>% dplyr::filter(maxEigen < 0 & !is.na(mu_x_se))

summary(resnd)

sum(resnd$mu_x - 1.96*resnd$mu_x_se < 2 & resnd$mu_x + 1.96*resnd$mu_x_se > 2)/nrow(resnd)


ggplot(resnd) +
  geom_density(aes(x=fullsample_wtp),color="black",alpha=0.5) +
  geom_density(aes(x=mu_x),fill="forestgreen",color=NA,alpha=0.5) +
  geom_vline(xintercept=mean(resnd$mu_x,na.rm=T),color="forestgreen",linetype="dashed") +
  geom_density(aes(x=uncorrected_wtp),fill="red",color=NA,alpha=0.5) +
  geom_density(aes(x=adhoc_wtp),fill="blue",color=NA,alpha=0.5) +
  labs(x="WTP Estimate",y="Density") +
  geom_vline(xintercept=2) +
  lims(x=c(1,3))

#write.csv(res,"~/selection-sims/simresults-8-7-2021.csv")
