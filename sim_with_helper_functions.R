require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,data.table,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo,haven,readstata13,broom,stringr)


# this code handles data generation and estimation for the sim results
# comment/uncomment the appropriate code chunks in lines ~350-370 to switch between heteroskedasticity strategies


source("~/selection-sims/selection_helpers.R")



mat1 <- c("a1","b1","b2","c1","c2","c3","d1","d2","d3","d4")


# names of correlated variables in selection equation (alpha, in this case)
# and in the choice equation: statusquo, ben, cost
data1 <- c("alpha","statusquo","att","cost")

# adds a function to the global environment called "rescaler"
# that is used to rescale the index in the choice equation
build_rescale_function(corvars=data1)


results <- readRDS("~/selection-sims/100runsof2k_2.rds")
hessians <- readRDS("~/selection-sims/100runsof2k_2_hessians.rds")

for (i in 92:150) {
  
  n <- 2000
  
  # set true/average param values
  alpha <- 0
  a_iv <- 1
  b_statquo <- -1
  b_cost <- -5
  b_att <- 5
  
  ### for the correlation matrix for data generation
  
  sigma_alpha_beta_11 <- 0.5 # alpha constant and statquo (response probability uncorrelated with marginal utility for status quo)
  sigma_alpha_beta_12 <- 0.5 # alpha constant and att (likely responders have higher m.u. for attribute)
  sigma_alpha_beta_13 <- 0.5 # alpha constant and cost (likely responders have lower m.u. for $$$)
  sigma_beta_12 <- 0 # statquo and att
  sigma_beta_13 <- 0 # statquo and cost
  sigma_beta_23 <- 0 # att and cost 
  
  # utility functions
  u_respond <- function (iv,eps_alpha,omega) {(alpha + eps_alpha) + a_iv*iv + omega}
  u_choice <- function (statquo,eps_statquo,att,eps_att,cost,eps_cost,eta) {
    (b_statquo+eps_statquo)*statquo + (b_att+eps_att)*att + (b_cost+eps_cost)*cost + eta
  }
  
  
  set.seed(i)
  
  # draw the correlated components
  Sigma <- matrix(nrow=4,ncol=4,byrow=T,
                  data=c(1,sigma_alpha_beta_11,sigma_alpha_beta_12,sigma_alpha_beta_13,
                         sigma_alpha_beta_11,1,sigma_beta_12,sigma_beta_13,
                         sigma_alpha_beta_12,sigma_beta_12,1,sigma_beta_23,
                         sigma_alpha_beta_13,sigma_beta_13,sigma_beta_23,1))
  
  # eta is the policy choice error
  # omega is the response/non-response error
  # epsilons are random components of utility parameters
  
  epsilons <- mvrnorm(n=n,mu=c(0,0,0,0),Sigma=Sigma)
  
  
  dataf <- data.frame(id=1:n,
                      eta1=rlogis(n), # independent draw
                      eta2=rlogis(n), # independent draw
                      omega=rlogis(n), # independent draw
                      eps_alpha=epsilons[,1], # alpha correlated with choice parameters
                      eps_statquo=epsilons[,2], # beta_statquo correlated with alpha
                      eps_att=epsilons[,3], # beta_att correlated with alpha
                      eps_cost=epsilons[,4], # beta_cost correlated with alpha
                      iv=rnorm(n)) # independent draw for the "instrument" that affects selection eq but not choice eq
  
  dataf <- rbind(dataf,dataf) %>% arrange(id) # two rows per person
  
  # generate alternatives (all available alternatives on a single row for apollo)
  dataf <- dataf %>% mutate(choicetype=rep(c("R","S"),n),
                            cost1=rep(runif(n,0,3),each=2),
                            cost2=rep(runif(n,0,3),each=2),
                            att1 = rep(runif(n,0,3),each=2),
                            att2 = rep(runif(n,0,3),each=2),
                            statquo=-1)
  
  # calculate utilities for r/nr and policy-choice alternatives-- u(statquo alternative)=0, u(non-response)=0
  dataf <- dataf %>% mutate(u_respond1=u_respond(iv,eps_alpha,omega),
                            u_choice1=u_choice(statquo,eps_statquo,att1,eps_att,cost1,eps_cost,eta1),
                            u_choice2=u_choice(statquo,eps_statquo,att2,eps_att,cost2,eps_cost,eta2))
  
  # choice key
  # 1 - don't respond
  # 2 - respond
  # 3 - choose statquo
  # 4 - choose alt 1
  # 5 - choose alt 2
  
  # make the choices
  dataf$choice <- ifelse(dataf$choicetype=="R" & dataf$u_respond1 < 0,1,NA)
  dataf$choice <- ifelse(dataf$choicetype=="R" & dataf$u_respond1 >= 0,2,dataf$choice)
  dataf$choice <- ifelse(dataf$choicetype=="S" & dataf$u_choice1 <= 0 & dataf$u_choice2 <= 0,3,dataf$choice)
  dataf$choice <- ifelse(dataf$choicetype=="S" & dataf$u_choice1 > 0 & dataf$u_choice1 >= dataf$u_choice2,4,dataf$choice)
  dataf$choice <- ifelse(dataf$choicetype=="S" & dataf$u_choice2 > 0 & dataf$u_choice2 > dataf$u_choice1,5,dataf$choice)
  
  
  
  # label individuals responders/nonresponders
  dataf <- dataf %>% group_by(id) %>%
    mutate(response_status = min(choice)) %>%
    ungroup %>%
    mutate(response_status = ifelse(response_status==1,"NONRESPONDER","RESPONDER"))
  
  # apollo requires "availability" dummies
  dataf$av1 <- ifelse(dataf$choice %in% 1:2,1,0)
  dataf$av2 <- ifelse(dataf$choice %in% 1:2,1,0)
  dataf$av3 <- ifelse(dataf$choice %in% 3:5,1,0)
  dataf$av4 <- ifelse(dataf$choice %in% 3:5,1,0)
  dataf$av5 <- ifelse(dataf$choice %in% 3:5,1,0)
  
  # reshape data for clogit
  
  clogitdata <- dataf %>% filter(choicetype=="S") %>% mutate(cost0=0,att0=0)
  
  clogitdata <- clogitdata %>% pivot_longer(cols=starts_with("cost"),names_to="cost_drop",values_to="cost") %>%
    pivot_longer(cols=starts_with("att"),names_to="att_drop",values_to="att") %>%
    dplyr::select(id,choice,response_status,att,cost,att_drop,cost_drop) %>%
    mutate(att_drop=str_extract(att_drop,"[[:digit:]]"),
           cost_drop=str_extract(cost_drop,"[[:digit:]]")) %>%
    filter(att_drop==cost_drop) %>%
    mutate(alt=as.numeric(att_drop)+3) %>%
    mutate(choice=ifelse(choice==alt,1,0)) %>%
    mutate(statquo=ifelse(alt %in% 4:5,-1,0))
  
  # clogit everyone
  right <- clogit(choice ~ att + cost + statquo + strata(id),data=clogitdata)
  # clogit responders only (uncorrected)
  wrong <- clogit(choice ~ att + cost + statquo + strata(id),data=clogitdata[clogitdata$response_status=="RESPONDER",])
  summary(right)
  summary(wrong)
  
  
  
  df <- dataf %>% filter(!(response_status=="NONRESPONDER" & choicetype=="S"))
  
  
  
  database <- df %>% ungroup() %>% as.data.frame()
  

  
  ########### convert iv into a dummy variable to see if it still works ############
  # database$iv <- ifelse(database$iv>0,1,0)
  ##################################################################################
  
  
  ### Initialise code
  apollo_initialise()
  
  ### Set core controls
  ### NOTE: YOU MAY NEED TO LOWER NCORES BEFORE RUNNING. MOST COMPUTERS MAX OUT AT 8, IN WHICH CASE YOU SHOULD USE NO MORE THAN 7.
  ### ALSO: "NCORES" SHOULD REALLY BE CALLED "NTHREADS". TYPICAL TO HAVE 2 THREADS PER CORE.
  
  apollo_control <- list(
    modelName ="Selection correction",
    modelDescr ="Mixed logit model",
    indivID   ="id",  
    mixing    = TRUE,
    nCores    = 15,
    panelData = TRUE
  )
  
  # ################################################################# #
  #### DEFINE MODEL PARAMETERS                                     ####
  # ################################################################# #
  
  ## Vector of parameters, including any that are kept fixed in estimation
  apollo_beta <- c(
    #response variables
    
    alpha_iv = 5,
    mu_alpha = 0,
    
    # choice variables
    
    mu_att = 26,
    mu_statquo = -4,
    mu_cost = -26,
    beta_imr = 6,
    
    
    
    a1 = 9,
    b1 = 7,
    b2 = 2,
    c1 = 1.6,
    c2 = 2,
    c3 = -2.5,
    d1 = 1,
    d2 = 2,
    d3 = 1,
    d4 = -1,
    
    # scales
    
    scale_RP =  1,
    scale_SP = 1
  )
  
  
  
  
  ### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
  apollo_fixed <<- c("scale_RP","scale_SP")
  
  # ################################################################# #
  #### DEFINE RANDOM COMPONENTS                                    ####
  # ################################################################# #
  
  ### Set parameters for generating draws
  apollo_draws <<- list(
    interDrawsType = "halton",
    interNDraws    = 500,
    interUnifDraws = c(),
    interNormDraws = c("draws_alpha","draws_statquo","draws_att","draws_cost"),
    intraDrawsType = "halton",
    intraNDraws    = 0,
    intraUnifDraws = c(),
    intraNormDraws = c()
  )
  
  
  
  
  ### Create random parameters
  apollo_randCoeff <<- function(apollo_beta, apollo_inputs){
    
    randcoeff = list()
    
    randcoeff[["random_alpha"]]   = mu_alpha +   a1*draws_alpha
    
    randcoeff[["random_statquo"]] = mu_statquo + b1*draws_alpha + b2*draws_statquo
    
    randcoeff[["random_att"]]     = mu_att +     c1*draws_alpha + c2*draws_statquo + c3*draws_att
    
    randcoeff[["random_cost"]]    = mu_cost +    d1*draws_alpha + d2*draws_statquo + d3*draws_att + d4*draws_cost
    
    
    return(randcoeff)
  }
  
  
  
  
  
  # ################################################################# #
  #### GROUP AND VALIDATE INPUTS                                   ####
  # ################################################################# #
  
  apollo_inputs <<- apollo_validateInputs(silent=F)
  
  # ################################################################# #
  #### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
  # ################################################################# #
  
  
  
  apollo_probabilities <<- function(apollo_beta, apollo_inputs, functionality="estimate"){
    
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
    V[['alt2']] = random_alpha + alpha_iv * iv
    V[['alt3']] = 0
    
    # OPTION 1: USE THIS CHUNK TO IGNORE HETEREOSKEDASTICITY
    # V[['alt4']] = (random_cost*cost1 + random_att*att1 + random_statquo*statquo +
    #                  beta_imr*imr((alpha_iv*iv + mu_alpha)/sqrt(a1^2+(pi^2/3))))
    # V[['alt5']] = (random_cost*cost2 + random_att*att2 + random_statquo*statquo +
    #                  beta_imr*imr((alpha_iv*iv + mu_alpha)/sqrt(a1^2+(pi^2/3))))
    
    # OPTION 2: USE THIS CHUNK TO USE AVERAGE OF NON-STATQUO ATTRIBUTE LEVELS
    V[['alt4']] = (random_cost*cost1 + random_att*att1 + random_statquo*statquo +
                     beta_imr*imr((alpha_iv*iv + mu_alpha)/sqrt(a1^2+(pi^2/3))))/
      rescaler(a1,b1,b2,c1,c2,c3,d1,d2,d3,d4,0,-1,(att1+att2)/2,(cost1+cost2)/2)
    
    V[['alt5']] = (random_cost*cost2 + random_att*att2 + random_statquo*statquo +
                     beta_imr*imr((alpha_iv*iv + mu_alpha)/sqrt(a1^2+(pi^2/3))))/
      rescaler(a1,b1,b2,c1,c2,c3,d1,d2,d3,d4,0,-1,(att1+att2)/2,(cost1+cost2)/2)
    
    # OPTION 3: USE THIS CHUNK FOR IJ-SPECIFIC HETEROSKEDASTICITY SCALING
    # V[['alt4']] = (random_cost*cost1 + random_att*att1 + random_statquo*statquo +
    #                  beta_imr*imr((alpha_iv*iv + mu_alpha)/sqrt(a1^2+(pi^2/3))))/
    #   stupidmatrix(a1,b1,b2,c1,c2,c3,d1,d2,d3,d4,statquo,att1,cost1)
    # V[['alt5']] = (random_cost*cost2 + random_att*att2 + random_statquo*statquo +
    #                  beta_imr*imr((alpha_iv*iv + mu_alpha)/sqrt(a1^2+(pi^2/3))))/
    #   stupidmatrix(a1,b1,b2,c1,c2,c3,d1,d2,d3,d4,statquo,att2,cost2)
    
    
    
    mnl_settings = list(
      alternatives  = c(alt1=1, alt2=2, alt3=3, alt4=4,alt5=5),
      avail         = list(alt1=av1, alt2=av2,alt3=av3,alt4=av4,alt5=av5),
      choiceVar     = choice,
      V             = lapply(V,"*",scale_RP),
      rows          = (choicetype=="R"),
      componentName = "MNL-RP"
    )
    
    P[['RP']] = apollo_mnl(mnl_settings, functionality)
    
    ### Compute probabilities for the SP part of the data using MNL model
    mnl_settings$V = lapply(V, "*", scale_SP)
    mnl_settings$rows = (choicetype=="S")
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
  
  
  #startpars <- model$estimate
  
  # ################################################################# #
  #### MODEL OUTPUTS                                               ####
  # ################################################################# #
  
  # ----------------------------------------------------------------- #
  #---- FORMATTED OUTPUT (TO SCREEN)                               ----
  # ----------------------------------------------------------------- #
  
  apollo_modelOutput(model)
  
  #save(dataf, database,clogitdata,wrong,right,model, file = "~/selection-sims/simresults20k3.RData")
  
  
  #########################################################
  ######### CODE GRAVEYARD BELOW ##########################
  #########################################################
  
  ests <- model$estimate
  
  with(as.list(ests),{
    vcov =  c(a1^2,  a1*b1,       a1*c1,             a1*d1,
              a1*b1, b1^2+b2^2,   b1*c1+b2*c2,       b1*d1+b2*d2,
              a1*c1, b1*c1+b2*c2, c1^2+c2^2+c3^2,    c1*d1+c2*d2+c3*d3,
              a1*d1, b1*d1+b2*d2, c1*d1+c2*d2+c3*d3, d1^2+d2^2+d3^2+d4^2)
    matrix(vcov,nrow=4)
  }) %>% eigen()
  
  with(as.list(ests),{
    vcov =  c(a1^2,  a1*b1,       a1*c1,             a1*d1,
              a1*b1, b1^2+b2^2,   b1*c1+b2*c2,       b1*d1+b2*d2,
              a1*c1, b1*c1+b2*c2, c1^2+c2^2+c3^2,    c1*d1+c2*d2+c3*d3,
              a1*d1, b1*d1+b2*d2, c1*d1+c2*d2+c3*d3, d1^2+d2^2+d3^2+d4^2)
    matrix(vcov,nrow=4)
  }) %>% cov2cor()
  
  
  mat <- with(as.list(ests),{
    vcov =  c(a1^2,  a1*b1,       a1*c1,             a1*d1,
              a1*b1, b1^2+b2^2,   b1*c1+b2*c2,       b1*d1+b2*d2,
              a1*c1, b1*c1+b2*c2, c1^2+c2^2+c3^2,    c1*d1+c2*d2+c3*d3,
              a1*d1, b1*d1+b2*d2, c1*d1+c2*d2+c3*d3, d1^2+d2^2+d3^2+d4^2)
    matrix(vcov,nrow=4)
  })
  
  drawwtp <- mvrnorm(n=20000,mu=ests[c("mu_alpha","mu_statquo","mu_att","mu_cost")],Sigma=mat)
  
  quantile(-1*drawwtp[,3]/drawwtp[,4],c(0.025,0.5,0.975))
  quantile(-1*drawwtp[,2]/drawwtp[,4],c(0.025,0.5,0.975))
  
  
  
  ests <- model$estimate
  
  
  
  ses <- model$robse
  names(ses) <- paste0(names(ses),"_se")
  ests <- c(ests,ses)
  ests <- ests[sort(names(ests))]
  
  a <- tidy(wrong) %>%
    mutate(term=paste0("wrong_",term)) %>%
    column_to_rownames(var="term") %>%
    t() %>%
    as.data.frame()
  
  wrong_est <- unlist(a[1,])
  wrong_ses <- unlist(a[2,])
  names(wrong_ses) <- paste0(names(wrong_ses),"_se")
  wrong_est <- c(wrong_est,wrong_ses)
  wrong_est <- wrong_est[sort(names(wrong_est))]
  
  
  ests <- c(ests,wrong_est)
  
  a <- tidy(right) %>%
    mutate(term=paste0("right_",term)) %>%
    column_to_rownames(var="term") %>%
    t() %>%
    as.data.frame()
  
  right_est <- unlist(a[1,])
  right_ses <- unlist(a[2,])
  names(right_ses) <- paste0(names(right_ses),"_se")
  right_est <- c(right_est,right_ses)
  right_est <- right_est[sort(names(right_est))]
  
  
  ests <- c(ests,right_est)
  
  ests
  
  results <- rbind(results,ests) %>% as.data.frame()
  
  hessians[[i]] <- model$hessian
  
  cat("i = ", i, ", wtp = ",-1*ests[["mu_att"]]/ests[["mu_cost"]],"\n")
  
  
}

# first is uncorrected for heterosked
# second is corrected for heterosked using average sd from both non-statquo alts
# third is trudy method as written


#saveRDS(results,"~/selection-sims/100runsof2k_2.rds")
#saveRDS(hessians,"~/selection-sims/100runsof2k_2_hessians.rds")

results$ND_hessian <- lapply(hessians,function (x) {
  ifelse(!is.na(x[1,1]),sum(eigen(x)$values>=0)==0,0)
    }
  ) %>% unlist() %>% as.numeric()

summary(-1*results$mu_att[results$ND_hessian==1]/results$mu_cost[results$ND_hessian==1])

results2 <- results %>% filter(ND_hessian==1)


wtps <- NULL

for (i in 1:nrow(results2)) {
  newrow <- with(results2[i,],
                 {
                   Sigma <- un_cholesky(c(a1,b1,b2,c1,c2,c3,d1,d2,d3,d4))
                   draws <- mvrnorm(n=100000,mu=c(mu_alpha,mu_statquo,mu_att,mu_cost),Sigma=Sigma)
                   quantile(-1*draws[,3]/draws[,4],c(0.025,.5,.975))
                 }
  )
  wtps <- rbind(wtps,newrow)
  
}

wtps <- wtps %>% as.data.frame()
summary(wtps)



ggplot(wtps) +
  geom_density(aes(x=`50%`),fill="lightblue",color="lightblue",alpha=0.5) +
  geom_density(data=results2, aes(x=-wrong_att/wrong_cost),color="pink",fill="pink",alpha=0.5) +
  geom_vline(data=results2,aes(xintercept=mean(-wrong_att/wrong_cost)),color="red") +
  geom_vline(data=results2,aes(xintercept=median(-wrong_att/wrong_cost)),color="red",linetype="dashed") +
  geom_vline(aes(xintercept=median(`50%`)),linetype="dashed",color="blue") +
  geom_vline(aes(xintercept=mean(`50%`)),color="blue") +
  geom_vline(aes(xintercept=1),color="#24a62f",size=1) +
  geom_hline(aes(yintercept=0)) +
  labs(x="",y="") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        text=element_text(size=11,  family="serif"),
        panel.grid=element_blank()) 
  
  





ggplot() + 
  geom_density(data=results2,aes(x=-mu_att/mu_cost),color="blue") + 
  geom_density(data=results2, aes(x=-wrong_att/wrong_cost),color="red") +
  geom_density(data=results2, aes(x=-right_att/right_cost),color="darkgreen") +
  geom_vline(data=results2,aes(xintercept=median(-mu_att/mu_cost)),color="blue") +
  geom_vline(data=results2,aes(xintercept=median(-wrong_att/wrong_cost)),color="red") +
  geom_vline(data=results2,aes(xintercept=median(-right_att/right_cost)),color="darkgreen")
  

ggplot() + 
  geom_density(data=results2,aes(x=-mu_statquo/mu_cost),color="blue") + 
  geom_density(data=results2, aes(x=-wrong_statquo/wrong_cost),color="red") +
  geom_density(data=results2, aes(x=-right_statquo/right_cost),color="darkgreen") +
  geom_vline(data=results2,aes(xintercept=median(-mu_statquo/mu_cost)),color="blue") +
  geom_vline(data=results2,aes(xintercept=median(-wrong_statquo/wrong_cost)),color="red") +
  geom_vline(data=results2,aes(xintercept=median(-right_statquo/right_cost)),color="darkgreen")


ggplot() + 
  geom_density(data=results,aes(x=-mu_att/mu_cost),color="blue") + 
  geom_density(data=results, aes(x=-wrong_att/wrong_cost),color="red") +
  geom_density(data=results, aes(x=-right_att/right_cost),color="darkgreen") +
  geom_vline(data=results,aes(xintercept=mean(-mu_att/mu_cost)),color="blue") +
  geom_vline(data=results,aes(xintercept=mean(-wrong_att/wrong_cost)),color="red") +
  geom_vline(data=results,aes(xintercept=mean(-right_att/right_cost)),color="darkgreen")


ggplot() + 
  geom_density(data=results,aes(x=-mu_statquo/mu_cost),color="blue") + 
  geom_density(data=results, aes(x=-wrong_statquo/wrong_cost),color="red") +
  geom_density(data=results, aes(x=-right_statquo/right_cost),color="darkgreen") +
  geom_vline(data=results,aes(xintercept=mean(-mu_statquo/mu_cost)),color="blue") +
  geom_vline(data=results,aes(xintercept=mean(-wrong_statquo/wrong_cost)),color="red") +
  geom_vline(data=results,aes(xintercept=mean(-right_statquo/right_cost)),color="darkgreen")
