require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,data.table,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo,haven,readstata13,broom,stringr)




n <- 1000

# set true param values
alpha <- 0
a_iv <- 1
b_statquo <- 0
b_cost <- -3
b_att <- 3

### for the correlation matrix

sigma_alpha_beta_11 <- 0.5 # alpha constant and statquo
sigma_alpha_beta_12 <- 0.5 # alpha constant and att
sigma_alpha_beta_13 <- 0.5 # alpha constant and cost
sigma_alpha_iv <- 0.5 # alpha constant and iv
sigma_beta_12 <- 0 # statquo and att
sigma_beta_13 <- 0 # statquo and cost
sigma_beta_23 <- 0 # att and cost


u_respond <- function (iv,eps_alpha,omega) {(alpha + eps_alpha) + a_iv*iv + omega}
u_choice <- function (statquo,eps_statquo,att,eps_att,cost,eps_cost,eta) {
  (b_statquo+eps_statquo)*statquo + (b_att+eps_att)*att + (b_cost+eps_cost)*cost + eta
}

imr <- function (x) {dnorm(-1*x)/(1-pnorm(-1*x))}

{
# hetero_scaler <- function (sigma_statquo,sigma_att,sigma_cost,rho,rho2,rho3,statquo,att,cost) {
#   
#   sqrt(
#     
#     matrix(c(statquo,att,cost),byrow=T,ncol=3) %*%
#       
#       matrix(c(sigma_statquo^2 + rho^2, rho*rho2,             rho*rho3,
#                rho*rho2,                sigma_att^2 + rho2^2, rho2*rho3,
#                rho*rho3,                rho2*rho3,            sigma_cost^2 + rho3^2),
#              byrow = T,nrow=3,ncol=3) %*%
#       
#       matrix(c(statquo,att,cost),byrow=T,ncol=1) +
#       
#       (pi^2)/3
#     
#   ) %>% as.numeric()
# }
# 
# 
# hetero_scaler2 <- function (sigma_statquo,sigma_att,sigma_cost,rho,rho2,rho3,statquo,att,cost) {
#   sqrt(
#   statquo*((sigma_statquo^2 + rho^2)*statquo + rho*rho2*att + rho*rho3*cost) +
#     att*(rho*rho2*statquo + (sigma_att^2 + rho2^2)*att + (rho2*rho3*cost)) +
#     cost*(rho*rho3*statquo + rho2*rho3*att + (sigma_cost^2+rho3^2)*cost) +
#     (pi^2)/3
#   )   
#     
# }
# 
# hetero_scaler2(1,2,3,1,2,5,1,1,1)
# hetero_scaler(1,2,3,1,2,5,1,1,1)
}

#cor_marg_sim <- function (n,alpha,a_iv,b_statquo,b_att,b_cost) {
  
  
  seed <- Sys.time() %>%
    str_extract_all("[[:digit:]]") %>%
    unlist() %>%
    paste0(collapse='') %>%
    str_sub(start=-8,end=-1) %>%
    as.numeric() %>%
    as.integer()

  set.seed(seed)
  
  # draw correlated components
  Sigma <- matrix(nrow=4,ncol=4,byrow=T,
                  data=c(1,sigma_alpha_beta_11,sigma_alpha_beta_12,sigma_alpha_beta_13,
                         sigma_alpha_beta_11,1,sigma_beta_12,sigma_beta_13,
                         sigma_alpha_beta_12,sigma_beta_12,1,sigma_beta_23,
                         sigma_alpha_beta_13,sigma_beta_13,sigma_beta_23,1))
  
  # eta is the policy choice error
  # omega is the response/non-response error
  # iv is correlated with eta but not with omega
  # epsilon att is correlated with omega and iv but not eta
  
  epsilons <- mvrnorm(n=n,mu=c(0,0,0,0),Sigma=Sigma)
  # reshape the random draws to logistic
  #epsilons[,1:2] <- qlogis(pnorm(omega_eta_epsilons[,1:2]))
  
  dataf <- data.frame(id=1:n,
                      eta1=rlogis(n),
                      eta2=rlogis(n),
                      omega=rlogis(n),
                      eps_alpha=epsilons[,1],
                      eps_statquo=epsilons[,2],
                      eps_att=epsilons[,3],
                      eps_cost=epsilons[,4],
                      iv=rnorm(n))
  dataf <- rbind(dataf,dataf) %>% arrange(id)
  dataf <- dataf %>% mutate(choicetype=rep(c("R","S"),n),
                            cost1=rep(runif(n,0,3),each=2),
                            cost2=rep(runif(n,0,3),each=2),
                            att1 = rep(runif(n,0,3),each=2),
                            att2 = rep(runif(n,0,3),each=2),
                            statquo=-1)
  
  dataf <- dataf %>% mutate(u_respond=u_respond(iv,eps_alpha,omega),
                            u_choice1=u_choice(statquo,eps_statquo,att1,eps_att,cost1,eps_cost,eta1),
                            u_choice2=u_choice(statquo,eps_statquo,att2,eps_att,cost2,eps_cost,eta2))
  
  # choice key
  # 1 - don't respond
  # 2 - respond
  # 3 - choose statquo
  # 4 - choose alt 1
  # 5 - choose alt 2
  
  
  dataf$choice <- ifelse(dataf$choicetype=="R" & dataf$u_respond < 0,1,NA)
  dataf$choice <- ifelse(dataf$choicetype=="R" & dataf$u_respond >= 0,2,dataf$choice)
  dataf$choice <- ifelse(dataf$choicetype=="S" & dataf$u_choice1 <= 0 & dataf$u_choice2 <= 0,3,dataf$choice)
  dataf$choice <- ifelse(dataf$choicetype=="S" & dataf$u_choice1 > 0 & dataf$u_choice1 >= dataf$u_choice2,4,dataf$choice)
  dataf$choice <- ifelse(dataf$choicetype=="S" & dataf$u_choice2 > 0 & dataf$u_choice2 > dataf$u_choice1,5,dataf$choice)
  
  
  
  
  dataf <- dataf %>% group_by(id) %>% 
    mutate(response_status = min(choice)) %>% 
    ungroup %>% 
    mutate(response_status = ifelse(response_status==1,"NONRESPONDER","RESPONDER"))
  
  dataf$av1 <- ifelse(dataf$choice %in% 1:2,1,0)
  dataf$av2 <- ifelse(dataf$choice %in% 1:2,1,0)
  dataf$av3 <- ifelse(dataf$choice %in% 3:5,1,0)
  dataf$av4 <- ifelse(dataf$choice %in% 3:5,1,0)
  dataf$av5 <- ifelse(dataf$choice %in% 3:5,1,0)
  
  # reshape data for clogit
  
  clogitdata <- dataf %>% filter(choicetype=="S") %>% mutate(cost0=0,att0=0)
  
  clogitdata <- clogitdata %>% pivot_longer(cols=starts_with("cost"),names_to="cost_drop",values_to="cost") %>% 
    pivot_longer(cols=starts_with("att"),names_to="att_drop",values_to="att") %>% 
    select(id,choice,response_status,att,cost,att_drop,cost_drop) %>% 
    mutate(att_drop=str_extract(att_drop,"[[:digit:]]"),
           cost_drop=str_extract(cost_drop,"[[:digit:]]")) %>% 
    filter(att_drop==cost_drop) %>% 
    mutate(alt=as.numeric(att_drop)+3) %>% 
    mutate(choice=ifelse(choice==alt,1,0)) %>% 
    mutate(statquo=ifelse(alt %in% 4:5,-1,0))
  
  right <- clogit(choice ~ att + cost + statquo + strata(id),data=clogitdata)
  wrong <- clogit(choice ~ att + cost + statquo + strata(id),data=clogitdata[clogitdata$response_status=="RESPONDER",])
  
  summary(right)
  
summary(wrong)
  
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
    nCores    = 1,
    panelData = TRUE
  )
  
  # ################################################################# #
  #### DEFINE MODEL PARAMETERS                                     ####
  # ################################################################# #
  
  ## Vector of parameters, including any that are kept fixed in estimation
  apollo_beta <<- c(
    #response variables
    
    alpha_iv = 7 + rnorm(1),
    mu_alpha = 2 + rnorm(1),
    sigma_alpha = 3,
    
    # choice variables
    
    sigma_cost = 3,
    mu_att = 2 + rnorm(1),
    sigma_att = 3,
    mu_statquo = .4 + rnorm(1),
    sigma_statquo = 3,
    beta_imr = -1.7 + rnorm(1),
    mu_cost = -2 + rnorm(1),
    

    
    # correlations
    
    sigma_alpha_beta_11_hat = 2.2, # alpha constant and statquo
    sigma_alpha_beta_12_hat = 0.2, # alpha constant and att
    sigma_alpha_beta_13_hat = -0.1, # alpha constant and cost
    sigma_beta_12_hat = 0, # statquo and att
    sigma_beta_13_hat = 0, # statquo and cost
    sigma_beta_23_hat = 0, # att and cost
    
    # scales
    
    scale_RP =  0,
    scale_SP = 0
    )
  
  #apollo_beta <- startpars
  
  ### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
  apollo_fixed <<- c("sigma_beta_12_hat","sigma_beta_13_hat","sigma_beta_23_hat","scale_RP","scale_SP")
  
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
  
  
  cormaker <<- function (sigma_alpha,sigma_alpha_beta_11_hat,sigma_alpha_beta_12_hat,sigma_alpha_beta_13_hat,
                        sigma_statquo,sigma_beta_12_hat,sigma_beta_13_hat,
                        sigma_att,sigma_beta_23_hat,
                        sigma_cost,
                        draws_alpha,draws_statquo,draws_att,draws_cost,
                        i) {
    

    
    cormatrix <- drop(eigen(matrix(nrow=4,ncol=4,byrow=T,
                      data=c(sigma_alpha,sigma_alpha_beta_11_hat,sigma_alpha_beta_12_hat,sigma_alpha_beta_13_hat,
                             sigma_alpha_beta_11_hat,sigma_statquo,sigma_beta_12_hat,sigma_beta_13_hat,
                             sigma_alpha_beta_12_hat,sigma_beta_12_hat,sigma_att,sigma_beta_23_hat,
                             sigma_alpha_beta_13_hat,sigma_beta_13_hat,sigma_beta_23_hat,sigma_cost)),symmetric=T)$vectors %*%
           diag(sqrt(pmax(eigen(matrix(nrow=4,ncol=4,byrow=T,
                                       data=c(sigma_alpha,sigma_alpha_beta_11_hat,sigma_alpha_beta_12_hat,sigma_alpha_beta_13_hat,
                                              sigma_alpha_beta_11_hat,sigma_statquo,sigma_beta_12_hat,sigma_beta_13_hat,
                                              sigma_alpha_beta_12_hat,sigma_beta_12_hat,sigma_att,sigma_beta_23_hat,
                                              sigma_alpha_beta_13_hat,sigma_beta_13_hat,sigma_beta_23_hat,sigma_cost)))$values,0)),4))
    dumb = cbind(c(draws_alpha),c(draws_statquo),c(draws_att),c(draws_cost))
    
    dumb2 = apply(dumb,1,function (x) {cormatrix %*% x})
    
    return(matrix(data=dumb2[i,],nrow=nrow(database)))
  }
  
  
  ### Create random parameters
  apollo_randCoeff <<- function(apollo_beta, apollo_inputs){
    
    randcoeff = list()
    
    randcoeff[["random_alpha"]] = draws_alpha
    
    randcoeff[["random_statquo"]] = draws_statquo
     
    randcoeff[["random_att"]] = draws_att
    
    randcoeff[["random_cost"]] = draws_cost
    
  
    # randcoeff[["random_alpha"]] = cormaker(sigma_alpha,sigma_alpha_beta_11_hat,sigma_alpha_beta_12_hat,sigma_alpha_beta_13_hat,
    #                                        sigma_statquo,sigma_beta_12_hat,sigma_beta_13_hat,
    #                                        sigma_att,sigma_beta_23_hat,
    #                                        sigma_cost,
    #                                        draws_alpha,draws_statquo,draws_att,draws_cost,
    #                                        1) + mu_alpha
    # 
    # randcoeff[["random_statquo"]] = cormaker(sigma_alpha,sigma_alpha_beta_11_hat,sigma_alpha_beta_12_hat,sigma_alpha_beta_13_hat,
    #                                          sigma_statquo,sigma_beta_12_hat,sigma_beta_13_hat,
    #                                          sigma_att,sigma_beta_23_hat,
    #                                          sigma_cost,
    #                                          draws_alpha,draws_statquo,draws_att,draws_cost,
    #                                          2) + mu_statquo
    # 
    # randcoeff[["random_att"]] = cormaker(sigma_alpha,sigma_alpha_beta_11_hat,sigma_alpha_beta_12_hat,sigma_alpha_beta_13_hat,
    #                                      sigma_statquo,sigma_beta_12_hat,sigma_beta_13_hat,
    #                                      sigma_att,sigma_beta_23_hat,
    #                                      sigma_cost,
    #                                      draws_alpha,draws_statquo,draws_att,draws_cost,
    #                                      3) + mu_att
    # 
    # randcoeff[["random_cost"]] = cormaker(sigma_alpha,sigma_alpha_beta_11_hat,sigma_alpha_beta_12_hat,sigma_alpha_beta_13_hat,
    #                                       sigma_statquo,sigma_beta_12_hat,sigma_beta_13_hat,
    #                                       sigma_att,sigma_beta_23_hat,
    #                                       sigma_cost,
    #                                       draws_alpha,draws_statquo,draws_att,draws_cost,
    #                                       4) + mu_cost
    # 
    
    return(randcoeff)
  }
  
  # ### Create random parameters
  # apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  #   randcoeff = list()
  #   {
  #     cormaker = drop(eigen(matrix(nrow=4,ncol=4,byrow=T,
  #                                   data=c(sigma_alpha,sigma_alpha_beta_11_hat,sigma_alpha_beta_12_hat,sigma_alpha_beta_13_hat,
  #                                          sigma_alpha_beta_11_hat,sigma_statquo,sigma_beta_12_hat,sigma_beta_13_hat,
  #                                          sigma_alpha_beta_12_hat,sigma_beta_12_hat,sigma_att,sigma_beta_23_hat,
  #                                          sigma_alpha_beta_13_hat,sigma_beta_13_hat,sigma_beta_23_hat,sigma_cost)),symmetric=T)$vectors %*%
  #                        diag(sqrt(pmax(eigen(matrix(nrow=4,ncol=4,byrow=T,
  #                                                    data=c(sigma_alpha,sigma_alpha_beta_11_hat,sigma_alpha_beta_12_hat,sigma_alpha_beta_13_hat,
  #                                                           sigma_alpha_beta_11_hat,sigma_statquo,sigma_beta_12_hat,sigma_beta_13_hat,
  #                                                           sigma_alpha_beta_12_hat,sigma_beta_12_hat,sigma_att,sigma_beta_23_hat,
  #                                                           sigma_alpha_beta_13_hat,sigma_beta_13_hat,sigma_beta_23_hat,sigma_cost)))$values,0)),4))
  #     
  #     dumb = cbind(c(draws_alpha),c(draws_statquo),c(draws_att),c(draws_cost))
  #     
  #     dumb2 = apply(dumb,1,function (x) {cormaker %*% x})
  #     
  #   }
  #   randcoeff[["random_alpha"]] = matrix(data=dumb2[1,],nrow=nrow(database))
  #   
  #   randcoeff[["random_statquo"]] = matrix(data=dumb2[2,],nrow=nrow(database))
  #   
  #   randcoeff[["random_att"]] = matrix(data=dumb2[3,],nrow=nrow(database))
  #   
  #   randcoeff[["random_cost"]] = matrix(data=dumb2[4,],nrow=nrow(database))
  #   
  #   
  #   return(randcoeff)
  # }
  # 

  

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
    V[['alt2']] = cormaker(sigma_alpha,sigma_alpha_beta_11_hat,sigma_alpha_beta_12_hat,sigma_alpha_beta_13_hat,
                           sigma_statquo,sigma_beta_12_hat,sigma_beta_13_hat,
                           sigma_att,sigma_beta_23_hat,
                           sigma_cost,
                           random_alpha,random_statquo,random_att,random_cost,
                           1) + mu_alpha + alpha_iv * iv
    V[['alt3']] = 0
    #V[['alt4']] = beta_cost*cost + beta_att * att + random_statquo*statquo + beta_imr*imr0
    V[['alt4']] = (random_cost*cost1 + random_att*att1 + random_statquo*statquo + 
                     beta_imr*imr((alpha_iv*iv + mu_alpha)/sqrt(sigma_alpha^2+(pi^2/3))))
    V[['alt5']] = (random_cost*cost2 + random_att*att2 + random_statquo*statquo + 
                     beta_imr*imr((alpha_iv*iv + mu_alpha)/sqrt(sigma_alpha^2+(pi^2/3))))
    # need to divide to re-scale for heteroskedasticity
    # large costs/atts will multiply components of compound error terms

    
    
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
  
#}


cor_marg_sim2 <- function (x) {
  cor_marg_sim(n=x,alpha=alpha,a_iv=a_iv,b_statquo=b_statquo,b_att=b_att,b_cost=b_cost)
}




tic()
test <- lapply(rep(1000,20),cor_marg_sim2)
toc()

res1 <- do.call(rbind, test) %>% as.data.frame()

tic()
test <- lapply(rep(2000,50),cor_marg_sim2)
toc()

res2 <- do.call(rbind, test) %>% as.data.frame()

res <- rbind(res1,res2)


tic()
test <- lapply(rep(1000,20),cor_marg_sim2)
toc()

res3 <- do.call(rbind, test) %>% as.data.frame()


res <- rbind(res1,res2)


sd(res$right_att/res$right_cost)/sqrt(60)
