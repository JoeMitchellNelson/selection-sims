n=200000

df <- data.frame(id=1:n,
                 w=rnorm(n),
                 m=rnorm(n),
                 eps=rlogis(n),
                 x=rnorm(n),
                 cost=rnorm(n),
                 eps2=rlogis(n))

a0 = 1
a1 = 2
a2 = 3

b0 = 1
b1 = 4
b2 = -2
bm = 2

df$u <- a0 + a1*df$w + a2*df$m + df$eps

df$y <- ifelse(df$u>0,1,0)
head(df)

summary(a <- glm(y ~ w,data=df,family="binomial"))

df$fitted <- a$fitted.values
df$ind <- predict(a)

df$mguess <- qnorm(0.5*(df$y-df$fitted)+0.5)


df$u2 <- b0 + (b1 + bm*df$m)*df$x + b2*df$cost + df$eps2
df$y2 <- ifelse(df$u2>0,1,0)
summary(wrong <- glm(y2 ~ x + cost, data=df[which(df$y==1),],family="binomial"))

summary(better <- glm(y2 ~ x*I(mguess) + cost,data=df[which(df$y==1),],family="binomial"))
summary(everyone <- glm(y2 ~ x + cost,data=df,family="binomial"))

-1*wrong$coefficients[["x"]]/wrong$coefficients[["cost"]]
-1*better$coefficients[["x"]]/better$coefficients[["cost"]]
-1*everyone$coefficients[["x"]]/everyone$coefficients[["cost"]]





ggplot(df[sample(5000),],aes(x=mguess,y=m)) +
  geom_point(aes(x=mguess,y=m,color=factor(y)),alpha=0.4) +
  geom_smooth(method="lm",se=F)

summary(lm(m ~ mguess, data=df))

# ggplot(df) +
#   geom_jitter(aes(x=pred,y=y,color=m),width=0,height=0.01) +
#   scale_color_viridis_c()

df$indbin <- cut(df$ind,100)
df$indbin <- str_remove_all(df$indbin,"\\(|\\)|\\[|\\]") %>% 
  str_split_fixed(",",2) %>% 
  as.data.frame() %>% apply(2,as.numeric) %>% rowMeans()

df$fitbin <- cut(df$fitted,100)
df$fitbin <- str_remove_all(df$fitbin,"\\(|\\)|\\[|\\]") %>% 
  str_split_fixed(",",2) %>% 
  as.data.frame() %>% apply(2,as.numeric) %>% rowMeans()

dfg <- df %>% group_by(y,indbin) %>% summarise(mean_a2_m = mean(a2*m),n=n(),sd=sd(a2*m)) %>% filter(n>20)
dfg2 <- df %>% group_by(y,fitbin) %>% summarise(mean_a2_m = mean(a2*m),n=n(),sd=sd(a2*m)) %>% filter(n>20)

dfg$se <- dfg$sd/sqrt(dfg$n)
dfg$lb <- dfg$mean_a2_m - 1.96*dfg$se
dfg$ub <- dfg$mean_a2_m + 1.96*dfg$se

f <- function (x) {dnorm(x)/(1-pnorm(x))}

f2 <- function (x) {
  out <- rep(NA,length(x))
  for (i in 1:length(x)) {
    eps <- seq(from=-4,to=4,length.out = 100)
    out[i] <- sum(f(x[i]-eps)*dlogis(eps))*(eps[2]-eps[1])
  }
 out
}

x = 3

f(4)

ggplot(dfg) +
  geom_point(aes(x=indbin,y=mean_a2_m,color=y))
  #geom_function(fun=f2)

ggplot(dfg2) +
  geom_point(aes(x=fitbin,y=mean_a2_m,color=y)) +
  geom_point(aes(x=fitbin,y=mean_a2_m,color=y))




ggplot(dfg) +
  geom_point(aes(x=fitbin,y=mean_m,color=y)) +
  geom_smooth(aes(x=fitbin,y=mean_m,group=y,color=y),se=F)


####################################

df <- data.frame(x=rnorm(n),
                 cost=rnorm(n),
                 b_x=2*rnorm(n)+4,
                 b_cost=-2,
                 eps=rlogis(n))
df <- df %>% mutate(choice = ifelse(b_x*x + b_cost*cost + eps > 0,1,0))

summary(test <- glm(choice ~ x + cost, data=df))

-1*test$coefficients[["x"]]/test$coefficients[["cost"]]
