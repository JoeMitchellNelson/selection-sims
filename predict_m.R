n=200000

df <- data.frame(id=1:n,
                 w=rnorm(n),
                 m=rnorm(n),
                 eps=rlogis(n),
                 x=rnorm(n),
                 cost=rnorm(n),
                 eps2=rlogis(n))

a0 = 0
a1 = 2
a2 = 3

b0 = 1
b1 = 4
b2 = -2
bm = -2

df$u <- a0 + a1*df$w + a2*df$m + df$eps

df$y <- ifelse(df$u>0,1,0)
head(df)

summary(a <- glm(y ~ w,data=df,family="binomial"))

df$fitted <- a$fitted.values

df$mguess <- qnorm(0.5*(df$y-df$fitted)+0.5)


df$u2 <- b0 + (b1 + bm*df$m)*df$x + b2*df$cost + df$eps2
df$y2 <- ifelse(df$u2>0,1,0)
summary(wrong <- glm(y2 ~ x + cost, data=df[which(df$y==1),],family="binomial"))

summary(better <- glm(y2 ~ x*mguess + cost,data=df[which(df$y==1),],family="binomial"))
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

df$fitbin <- cut(df$fitted,100)
df$fitbin <- str_remove_all(df$fitbin,"\\(|\\)|\\[|\\]") %>% 
  str_split_fixed(",",2) %>% 
  as.data.frame() %>% apply(2,as.numeric) %>% rowMeans()

dfg <- df %>% group_by(y,fitbin) %>% summarise(mean_m = mean(m),n=n(),sd=sd(m)) %>% filter(n>20)
dfg$se <- dfg$sd/sqrt(dfg$n)
dfg$lb <- dfg$mean_m - 1.96*dfg$se
dfg$ub <- dfg$mean_m + 1.96*dfg$se


f <- function (x) {qnorm(x)}

ggplot(dfg) +
  geom_point(aes(x=fitbin,y=mean_m,color=y)) +
  geom_point(aes(x=fitbin,y=mean_m,color=y)) +
  geom_errorbar(aes(ymin=lb,ymax=ub,x=fitbin,color=y))

ggplot(dfg) +
  geom_point(aes(x=0.5*(y-fitbin)+0.5,y=mean_m,color=y)) +
  geom_function(fun=f)


ggplot(dfg) +
  geom_point(aes(x=fitbin,y=mean_m,color=y)) +
  geom_smooth(aes(x=fitbin,y=mean_m,group=y,color=y),se=F)
