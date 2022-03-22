require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,data.table,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo,haven,readstata13,broom,stringr)

n <- 1000

# set true param values

b_iv <- 1
b_statquo <- 0
b_att <- 3
b_cost <- -2

u_respond <- function (iv,omega) {b_iv*iv + omega}
u_choice <- function (statquo,att,cost,eta) {b_statquo*statquo + b_att*att + b_cost*cost + eta}
imr <- function (x) {dnorm(-1*x)/(1-pnorm(-1*x))}


cor_const_sim <- function (n,b_iv,b_statquo,b_att,b_cost) {
  

  seed <- Sys.time() %>% 
    str_extract_all("[[:digit:]]") %>% 
    unlist() %>% 
    paste0(collapse='') %>%  
    str_sub(start=-8,end=-1) %>% 
    as.numeric() %>% 
    as.integer()
  
  set.seed(seed)
  
  # draw correlated components
  Sigma <- matrix(nrow=3,ncol=3,byrow=T,
                  data=c(1,.5,0,
                         .5,1,.5,
                         0,.5,1))
  
  eta_omega_iv <- mvrnorm(n=n,mu=c(0,0,0),Sigma=Sigma)
  # reshape the random draws to logistic
  eta_omega_iv <- qlogis(pnorm(eta_omega_iv))
  
  dataf <- data.frame(id=1:n,eta=eta_omega_iv[,1],omega=eta_omega_iv[,2],iv=eta_omega_iv[,3])
  dataf <- rbind(dataf,dataf) %>% arrange(id)
  dataf <- dataf %>% mutate(choicetype=rep(c("R","S"),n),
                            att=rep(runif(n,0,3),each=2),
                            cost=rep(runif(n,0,3),each=2),
                            statquo=1)
  
  dataf <- dataf %>% mutate(u_respond=u_respond(iv,omega),
                            u_choice=u_choice(statquo,att,cost,eta))
  
  # choice key
  # 1 - don't respond
  # 2 - respond
  # 3 - choose statquo
  # 4 - choose the thing
  
  
  dataf$choice <- ifelse(dataf$choicetype=="R" & dataf$u_respond < 0,1,NA)
  dataf$choice <- ifelse(dataf$choicetype=="R" & dataf$u_respond >= 0,2,dataf$choice)
  dataf$choice <- ifelse(dataf$choicetype=="S" & dataf$u_choice < 0,3,dataf$choice)
  dataf$choice <- ifelse(dataf$choicetype=="S" & dataf$u_choice >= 0,4,dataf$choice)
  
  
  
  dataf <- dataf %>% group_by(id) %>% 
    mutate(response_status = min(choice)) %>% 
    ungroup %>% 
    mutate(response_status = ifelse(response_status==1,"NONRESPONDER","RESPONDER"))
  
  dataf$av1 <- ifelse(dataf$choice %in% 1:2,1,0)
  dataf$av2 <- ifelse(dataf$choice %in% 1:2,1,0)
  dataf$av3 <- ifelse(dataf$choice %in% 3:4,1,0)
  dataf$av4 <- ifelse(dataf$choice %in% 3:4,1,0)
  
  
  wrong <- glm(I(choice==4) ~ 0 + statquo + att + cost,data=dataf[dataf$choicetype=="S" & dataf$response_status=="RESPONDER",],family = "binomial")
  
  right <- glm(I(choice==4) ~ 0 + statquo + att + cost,data=dataf[dataf$choicetype=="S",],family = "binomial")
  
  # summary(glm(I(choice==2) ~ iv,data=df[df$choicetype=="R",],family = "binomial"))
  
  
  df <- dataf %>% filter(!(response_status=="NONRESPONDER" & choicetype=="S"))
  
  
  database <<- df %>% ungroup() %>% as.data.frame()
  
  ### Initialise code
  apollo_initialise()
  
  ### Set core controls
  apollo_control <<- list(
    modelName ="Selection correction",
    modelDescr ="Mixed logit model",
    indivID   ="id",  
    mixing    = TRUE, 
    nCores    = 12,
    panelData = TRUE
  )
  
  # ################################################################# #
  #### DEFINE MODEL PARAMETERS                                     ####
  # ################################################################# #
  
  ## Vector of parameters, including any that are kept fixed in estimation
  apollo_beta <<- c(
    # choice variables
    
    beta_cost = -43 + rnorm(1),
    beta_att = 64 + rnorm(1),
    mu_statquo = 0 + rnorm(1),
    sigma_statquo = 30 + rnorm(1),
    beta_imr = 13 + rnorm(1),
    
    
    #response variables
    
    alpha_iv = 2 + rnorm(1),
    mu_alpha = 0 + rnorm(1),
    sigma_alpha = 1 + rnorm(1),
    
    # rho and scale
    
    rho = 29 + rnorm(1),
    scale_RP = 1,
    scale_SP = 1)
  
  #apollo_beta <- startpars
  
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
    interNormDraws = c("draws_alpha","draws_statquo"),
    intraDrawsType = "halton",
    intraNDraws    = 0,
    intraUnifDraws = c(),
    intraNormDraws = c()
  )
  
  ### Create random parameters
  apollo_randCoeff <<- function(apollo_beta, apollo_inputs){
    randcoeff = list()
    
    
   # dumb <- cbind(c(draws_alpha),c(draws_statquo))
    
    
    randcoeff[["random_alpha"]] = mu_alpha + sigma_alpha * draws_alpha
    #randcoeff[["b_statquo"]] = mu_statquo + sqrt(1-rho^2)*(sigma_statquo) * draws_statquo + rho * sigma_statquo * draws_alpha
    
    randcoeff[["random_statquo"]] = mu_statquo + sigma_statquo * draws_statquo + rho * draws_alpha
    
    
    #cat(dim(randcoeff[["random_alpha"]]),class(randcoeff[["random_alpha"]]),class(randcoeff))
    
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
    #V[['alt4']] = beta_cost*cost + beta_att * att + random_statquo*statquo + beta_imr*imr0
    V[['alt4']] = beta_cost*cost + beta_att * att + random_statquo*statquo + beta_imr*imr((alpha_iv*iv + mu_alpha)/sqrt(sigma_alpha^2+(pi^2/3)))
    
    
    
    mnl_settings = list(
      alternatives  = c(alt1=1, alt2=2, alt3=3, alt4=4),
      avail         = list(alt1=av1, alt2=av2,alt3=av3,alt4=av4),
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
  
  model <<- apollo_estimate(apollo_beta, apollo_fixed,
                          apollo_probabilities, apollo_inputs, 
                          estimate_settings=list(
                          #  hessianRoutine="numDeriv",
                             hessianRoutine="none",
                            silent=F))
  
  
  # ################################################################# #
  #### MODEL OUTPUTS                                               ####
  # ################################################################# #
  
  # ----------------------------------------------------------------- #
  #---- FORMATTED OUTPUT (TO SCREEN)                               ----
  # ----------------------------------------------------------------- #
  
  #apollo_modelOutput(model)
  
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
  
  suppressWarnings(rm(database,apollo_draws,apollo_control,apollo_inputs,dataf,eta_omega_iv,model,df),classes="warning")
  
  ests
  
}


cor_const_sim2 <- function (x) {
  cor_const_sim(n=x,b_iv=1,b_statquo=0,b_att=3,b_cost=-2)
}


cor_const_sim2(500)


tic()
test <- lapply(rep(1000,120),cor_const_sim2)
toc()

res1 <- do.call(rbind, test) %>% as.data.frame()

ggplot(res1) +
  geom_density(aes(x=mu_statquo/beta_cost)) +
  geom_density(aes(x=wrong_statquo/wrong_cost),color="red")

tic()
test2 <- lapply(rep(1000,20),cor_const_sim2)
toc()

res2 <- do.call(rbind, test2) %>% as.data.frame()


tic()
test3 <- lapply(rep(1000,20),cor_const_sim2)
toc()

res3 <- do.call(rbind, test3) %>% as.data.frame()


tic()
test4 <- lapply(rep(1000,105),cor_const_sim2)
toc()

res4 <- do.call(rbind, test4) %>% as.data.frame()

tic()
test5 <- lapply(rep(1000,150),cor_const_sim2)
toc()

res5 <- do.call(rbind, test5) %>% as.data.frame()

tic()
test6 <- lapply(rep(1000,200),cor_const_sim2)
toc()

res6 <- do.call(rbind, test6) %>% as.data.frame()


results <- rbind(res1,res2,res3,res4,res5,res6)


ggplot(results) +
  geom_density(aes(x=mu_statquo/beta_cost)) +
  geom_density(aes(x=wrong_statquo/wrong_cost),color="red")

ggplot(results) +
  geom_density(aes(x=beta_imr))



mean(results$mu_statquo/results$beta_cost)

mean(results$wrong_statquo/results$wrong_cost)



