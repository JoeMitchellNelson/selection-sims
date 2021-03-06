require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo)



  
  ### Clear memory
  
  rm(database)
  
  
  s <- Sys.time() %>% as.numeric()
  set.seed(s)
  
  n = 2000

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
  database$av_1 <- database$RP * 1
  database$av_2 <- database$RP * 1
  database$av_3 <- database$SP * 1
  database$av_4 <- database$SP * 1
  
  # parameters
  beta1 <- -2 # effect for cost (choice model)
  beta20 <- 4  # common component of effect for x
  beta21 <- 2  # m-component of effect for x
  
  alpha0 <- 0 # constant term, increase to boost the simulated response rate
  alpha1 <- 2 # coef for w in selection equation
  alpha2 <- 2 # coef for m in selection equation
  
  database$U12 <-  alpha0 + alpha1*database$w1 + alpha2*database$m + database$epsRP
  database$U34 <- beta1*(database$cost3) + (beta20 + beta21*database$m)*(database$x3) + database$epsSP
  
  
  database$choice1 <- ifelse(database$RP==1 & database$U12 >= 0,1,
                             ifelse(database$RP==1 & database$U12 < 0, 2,NA))
  database$choice2 <- ifelse(database$SP == 1 & database$U34 >= 0,3,
                             ifelse(database$SP==1 & database$U34 < 0,4,NA))
  
  database$choice <- ifelse(database$RP==1,database$choice1,database$choice2)
  
  
  database <- database %>% group_by(ID) %>% mutate(respond = sum(choice==1)) %>% ungroup()
  
  ######################################################
  ### clogit if 100% response ###
  ######################################################
  {
    df <- database %>% dplyr::filter(SP==1)
    df$choice <- ifelse(df$choice==3,1,0)
    df2 <- df
    df2$x3 <- 0
    df2$cost3 <- 0
    df2$choice <- (df2$choice - 1)^2
    df <- rbind(df,df2)
    df <- df[order(df$ID),]
    summary(cor <- clogit(choice ~ x3 + cost3 + strata(ID),data=df))
  }
  
  # code to remove SP for nonresponders
  database$choice <- ifelse(database$SP==1 & database$respond==0,NA,database$choice)
  database <- database %>% dplyr::filter(!is.na(choice))
  
  
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
    cat(paste0("full sample WTP is ",round(-1*cor$coefficients[["x3"]]/cor$coefficients[["cost3"]],3)))
    cat(paste0("\nuncorrected WTP is ",round(-1*uncor$coefficients[["x3"]]/uncor$coefficients[["cost3"]],3)))
  }
  
  ### Initialise code
  apollo_initialise()
  
  ### Set core controls
  apollo_control = list(
    modelName ="First stab at selection correction",
    modelDescr ="Mixed logit model",
    indivID   ="ID",  
    mixing    = TRUE, 
    nCores    = 2
  )
  
  # ################################################################# #
  #### DEFINE MODEL PARAMETERS                                     ####
  # ################################################################# #
  
  ### Vector of parameters, including any that are kept fixed in estimation
  apollo_beta=c(b_w = 1,
                mu_x = 1,
                sigma_x = 1,
                b_cost = -1,
                mu_alpha = 0,
                sigma_alpha = 1,
                sigma_x_alpha = 1,
                scale_RP = 1,
                scale_SP = 1)
  
  ### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
  apollo_fixed = c("scale_RP","b_cost")
  
  # ################################################################# #
  #### DEFINE RANDOM COMPONENTS                                    ####
  # ################################################################# #
  
  ### Set parameters for generating draws
  apollo_draws = list(
    interDrawsType = "halton",
    interNDraws    = 501,
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
    
    randcoeff[["alpha"]] = mu_alpha + sigma_alpha * draws_alpha + sigma_x_alpha * draws_x
    randcoeff[["b_x"]] = mu_x + sigma_x * draws_x
    
    # res <- mvrnorm(n=length(unique(database$ID)),mu=c(mu_x,mu_w),Sigma=matrix(c(1,.5,.5,1),byrow=T,nrow=2))
    # 
    # 
    # randcoeff[["b_x"]] = res[1]
    # randcoeff[["b_w"]] = res[2]
    
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
    
    
    mnl_settings = list(
      alternatives  = c(alt1=1, alt2=2, alt3=3, alt4=4),
      avail         = list(alt1=av_1, alt2=av_2,alt3=av_3,alt4=av_4),
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
                          estimate_settings=list(hessianRoutine="maxLik",silent=T))
  cat(paste0("full sample WTP is ",round(-1*cor$coefficients[["x3"]]/cor$coefficients[["cost3"]],3)))
  cat(paste0("\nuncorrected WTP is ",round(-1*uncor$coefficients[["x3"]]/uncor$coefficients[["cost3"]],3)))
  cat(paste0("\ncorrected WTP is ",round(-1*model$estimate[["mu_x"]]/model$estimate[["b_cost"]],3)))
  
  uwtp <- -1*uncor$coefficients[["x3"]]/uncor$coefficients[["cost3"]]
  cwtp <- -1*cor$coefficients[["x3"]]/cor$coefficients[["cost3"]]
  
  # ################################################################# #
  #### MODEL OUTPUTS                                               ####
  # ################################################################# #
  
  # ----------------------------------------------------------------- #
  #---- FORMATTED OUTPUT (TO SCREEN)                               ----
  # ----------------------------------------------------------------- #
  
  apollo_modelOutput(model)
  
  
  res0 <- c(s,model$maximum,
            ifelse(model$message=="successful convergence ",1,0),
            model$estimate,model$robse,model$eigValue,uwtp,cwtp)
  
  
  res0
  
  # res <- as.data.frame(t(res0))
  # names(res) <- c("seed","LL","message",
  #                 names(model$estimate),
  #                 paste0(names(model$robse),"_se"),
  #                 "maxEigen","uncorrected_wtp","fullsample_wtp")
  
  res <- rbind(res,res0)
  
  cat(paste0("\n\n\n\t\t\t",i,"\n\n\n"))
  
  cat(paste0(mean(res$mu_x[res$maxEigen<0],na.rm=T)),"\n\n")
  


ggplot(res[res$maxEigen<0,]) +
  geom_density(aes(x=mu_x),fill="blue",alpha=0.5,color=NA) +
  geom_density(aes(x=uncorrected_wtp),fill="red",alpha=0.5,color=NA) +
  geom_vline(xintercept=2) +
  geom_vline(xintercept=mean(res$mu_x[res$maxEigen<0],na.rm=T),color="blue",linetype="dashed") +
  geom_vline(xintercept=mean(res$uncorrected_wtp[res$maxEigen<0],na.rm=T),color="red",linetype="dashed")

(sum(res$mu_x[res$maxEigen<0] + 1.96*res$mu_x_se[res$maxEigen<0] < 2,na.rm=T) + 
    sum(res$mu_x[res$maxEigen<0] - 1.96*res$mu_x_se[res$maxEigen<0] > 2,na.rm=T))/sum(res$message[res$maxEigen<0],na.rm=T)

(sum(res$mu_x + 1.96*res$mu_x_se < 2,na.rm=T) + 
    sum(res$mu_x - 1.96*res$mu_x_se > 2,na.rm=T))/sum(res$message,na.rm=T)

ggplot(res) +
  geom_density(aes(x=mu_x),fill="blue",alpha=0.5,color=NA) +
  geom_density(aes(x=uncorrected_wtp),fill="red",alpha=0.5,color=NA) +
  geom_vline(xintercept=2) +
  geom_vline(xintercept=mean(res$mu_x,na.rm=T),color="blue",linetype="dashed") + 
  geom_vline(xintercept=mean(res$uncorrected_wtp,na.rm=T),color="red",linetype="dashed") +
  labs(x="Estimated WTP",y="Density") +
  theme_minimal()

#write.csv(res,"~/selection-sims/sim_results.csv")
res <- read.csv("~/selection-sims/sim_results.csv")
res <- res %>% filter(message==1 & !is.na(b_w_se))

sd(sqrt(res$sigma_x_alpha^2 + res$sigma_alpha^2))

sd(abs(res$sigma_x_alpha)*abs(res$sigma_x)/(abs(res$sigma_x)*sqrt(res$sigma_alpha^2 + res$sigma_x_alpha^2)))/sqrt(nrow(res))

a <- data.frame(x1=rnorm(100000)*17.83825,x2=rnorm(100000)*1)
cor(a$x2,a$x2 + a$x1)
