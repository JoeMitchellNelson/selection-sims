require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo)

n=500000

f <- function (x) {
  temp <- data.frame(x1=rnorm(1000000),x2=rlogis(1000000))
  temp$x1 <- temp$x1/(sqrt(1+pi^2/3)) * (sqrt(pi^2/3))
  temp$x2 <- temp$x2/(sqrt(1+pi^2/3)) * (sqrt(pi^2/3))
  
  out <- rep(NA,length(x))
  for (i in 1:length(x)) {
    out[i] <- mean(temp$x1[which(temp$x1 + temp$x2 >= x[i])])
  }
  out
}

df <- data.frame(id=1:n,
                 w=rnorm(n),
                 m=rnorm(n),
                 eps=rlogis(n,scale=1),
                 x=rnorm(n),
                 cost=rnorm(n),
                 eps2=rlogis(n))

a0 = 1
a1 = 2
a2 = 3

b0 = 1
b1 = 4
b2 = -2
bm = 1

df$u <- a0 + a1*df$w + a2*df$m + df$eps

df$y <- ifelse(df$u>0,1,0)
head(df)

summary(a <- glm(y ~ w,data=df,family="binomial"))

df$fitted <- a$fitted.values
df$ind <- predict(a)
df$ind2 <- df$ind - mean(df$ind)

df$indbin <- cut(df$ind,200)
df$indbin <- str_remove_all(df$indbin,"\\(|\\)|\\[|\\]") %>% 
  str_split_fixed(",",2) %>% 
  as.data.frame() %>% apply(2,as.numeric) %>% rowMeans()

df$fitbin <- cut(df$fitted,200)
df$fitbin <- str_remove_all(df$fitbin,"\\(|\\)|\\[|\\]") %>% 
  str_split_fixed(",",2) %>% 
  as.data.frame() %>% apply(2,as.numeric) %>% rowMeans()

dfg <- df %>% group_by(y,indbin) %>% summarise(mean_a2_m = mean(m),n=n(),sd=sd(m)) %>% filter(n>100)
dfg2 <- df %>% group_by(y,fitbin) %>% summarise(mean_a2_m = mean(m/a2),n=n(),sd=sd(m/a2)) %>% filter(n>100)

dfg$mean_m <- f(-1*dfg$indbin)

df$mguess <- qlogis(0.5*(df$y-df$fitted)+0.5)

df <- left_join(df,dfg[,which(names(dfg) %in% c("indbin","mean_m"))])

df$u2 <- b0 + (b1 + bm*df$m)*df$x + b2*df$cost + df$eps2
df$y2 <- ifelse(df$u2>0,1,0)

summary(wrong <- glm(y2 ~ x + cost, data=df[which(df$y==1),],family="binomial"))
summary(better <- glm(y2 ~ x + x:I(mguess) + cost,data=df[which(df$y==1),],family="binomial"))
summary(everyone <- glm(y2 ~ x +  x:m + cost,data=df[sample(sum(df$y==1)),],family="binomial"))

-1*wrong$coefficients[["x"]]/wrong$coefficients[["cost"]]
-1*better$coefficients[["x"]]/better$coefficients[["cost"]]
-1*everyone$coefficients[["x"]]/everyone$coefficients[["cost"]]

cor.test(df$mean_m[df$y==1],df$m[df$y==1])
cor.test(df$mguess[df$y==1],df$m[df$y==1])
cor.test(df$mean_m[df$y==1],df$mguess[df$y==1])

b <- lm(mean_a2_m ~ f(-1*indbin) + 0,data=dfg[which(dfg$y==1),],weight=n)

c <- lm(m ~ mean_m + 0, data=df[df$y==1,])

ggplot(dfg[which(dfg$y==1),]) +
  geom_point(aes(x=-1*indbin,y=(mean_a2_m-0)/c$coefficients[1],color=n<1000)) +
  geom_function(fun=f)


mean(df$m[df$y==1])
mean(df$mguess[df$y==1])
# distribution of m's around the (scaled) guess



testbin <- sort(unique(df$indbin))[100]
testmeanm <- unique(df$mean_m[df$indbin==testbin & df$y==1])
nrow(df[df$indbin==testbin & df$y==1,])

ggplot(df[df$indbin==testbin & df$y==1,]) +
  geom_density(aes(x=m)) +
  geom_vline(xintercept=b$coefficients[1]*testmeanm)





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
