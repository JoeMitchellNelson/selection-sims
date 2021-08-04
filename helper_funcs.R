

formatclogit <- function (db,n,J) {
  df <- db %>% dplyr::filter(SP==1)
  temp <- data.frame(ID = rep(1:n,each=J),alt=1:J + 2)
  df <- left_join(temp,df)
  df$choice <- ifelse(df$choice==df$alt,1,0)
  df$x <- ifelse(df$alt==3,df$x3,
                 ifelse(df$alt==4,df$x4,df$x5))
  df$cost <- ifelse(df$alt==3,df$cost3,
                    ifelse(df$alt==4,df$cost4,df$cost5))
  df <- df %>% dplyr::select(ID,alt,choice,x,cost,respond)
  df
}



bs <- function(data, indices,formula) {
  d1 <- data[indices,] # allows boot to select sample
  d1$fittedRP <- glm(respond ~ I(m+noise),data=d1,family="binomial")$fitted.values
  d1$fittedRP <- d1$fittedRP - mean(d1$fittedRP)
  
  d2 <- data.frame(ID=rep(d1$ID,each=3),alt=rep(3:5),fittedRP=rep(d1$fittedRP,each=3))
  d2 <- left_join(d2,data,by=c("ID"))
  d2$choice <- ifelse(d2$choice2==d2$alt,1,0)
  d2$x <- ifelse(d2$alt==3,d2$x3,
                 ifelse(d2$alt==4,d2$x4,d2$x5))
  d2$cost <- ifelse(d2$alt==3,d2$cost3,
                    ifelse(d2$alt==4,d2$cost4,d2$cost5))
  
  fit <- clogit(formula, data=d2[d2$respond==1,])
  wtps <- mvrnorm(n=10000,mu=fit$coefficients,Sigma=fit$var) %>% as.data.frame()
  wtp <- mean(-1*wtps$x/wtps$cost)
  return(wtp)
}



