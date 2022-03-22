require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,data.table,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo,haven,readstata13,broom,stringr,ggridges)

# Hi Trudy.
# Run one of the following three lines:

load("~/selection-sims/simresults20k.Rdata") # no heterosked scaling
load("~/selection-sims/simresults20k2.Rdata") # heterosked scaling using the average attribute levels in a choice set
load("~/selection-sims/simresults20k3.Rdata") # heterosked scaling using ij-specific attribute levels

# You should get the following objects in your environment:
# clogitdata - simulated data, shaped for R's clogit function (includes choice rows for non-responders)
# database - simulated data, shaped for apollo (non-responders' choice rows are dropped)
# dataf - simulated data, shaped for apollo (including choice rows for non-responders)
# model - results of mixed logit estimation in apollo
# right - clogit results on the full set of invitees
# wrong - clogit results for responders only (no correction)

# Let's have a look

# print the results the apollo model
apollo_modelOutput(model)




# do some checks on estimated VCOV MATRIX
# (matrix order is alpha, statquo, att, cost)

ests <- model$estimate


# eigen values positive?
with(as.list(ests),{
  vcov =  c(a1^2,  a1*b1,       a1*c1,             a1*d1,
            a1*b1, b1^2+b2^2,   b1*c1+b2*c2,       b1*d1+b2*d2,
            a1*c1, b1*c1+b2*c2, c1^2+c2^2+c3^2,    c1*d1+c2*d2+c3*d3,
            a1*d1, b1*d1+b2*d2, c1*d1+c2*d2+c3*d3, d1^2+d2^2+d3^2+d4^2)
  matrix(vcov,nrow=4)
}) %>% eigen()


# convert vcov matrix to a correlation matrix
# (had expected 0.5 in first row and first col, 1s on diag, 0s elsewhere...)
with(as.list(ests),{
  vcov =  c(a1^2,  a1*b1,       a1*c1,             a1*d1,
            a1*b1, b1^2+b2^2,   b1*c1+b2*c2,       b1*d1+b2*d2,
            a1*c1, b1*c1+b2*c2, c1^2+c2^2+c3^2,    c1*d1+c2*d2+c3*d3,
            a1*d1, b1*d1+b2*d2, c1*d1+c2*d2+c3*d3, d1^2+d2^2+d3^2+d4^2)
  matrix(vcov,nrow=4)
}) %>% cov2cor()

# compare that cor matrix to the cor matrix in the simulated data
getcors <- dataf %>% filter(choicetype=="R") %>% select(starts_with("eps"))
cor(getcors)


# store the vcov matrix as "mat" and use it to draw cost, attribute, statquo parameters
mat <- with(as.list(ests),{
  vcov =  c(a1^2,  a1*b1,       a1*c1,             a1*d1,
            a1*b1, b1^2+b2^2,   b1*c1+b2*c2,       b1*d1+b2*d2,
            a1*c1, b1*c1+b2*c2, c1^2+c2^2+c3^2,    c1*d1+c2*d2+c3*d3,
            a1*d1, b1*d1+b2*d2, c1*d1+c2*d2+c3*d3, d1^2+d2^2+d3^2+d4^2)
  matrix(vcov,nrow=4)
})

drawwtp <- mvrnorm(n=20000,mu=ests[c("mu_alpha","mu_statquo","mu_att","mu_cost")],Sigma=mat)

# print estimated median and 95 CI for att WTP (should be 1)
quantile(-1*drawwtp[,3]/drawwtp[,4],c(0.025,0.5,0.975))

# print estimated median and 95 CI for statquo WTP (should be -0.2)
quantile(-1*drawwtp[,2]/drawwtp[,4],c(0.025,0.5,0.975))


# compare those to the TRUE distributions of WTPs in the simulated data (including non-responders)
singleserve <- dataf %>% filter(choicetype=="R")
responders <- singleserve %>% filter(response_status=="RESPONDER")

quantile(-1*(5+singleserve$eps_att)/(-5+singleserve$eps_cost),c(0.025,0.5,0.975)) # WTP for att
quantile(-1*(-1+singleserve$eps_statquo)/(-5+singleserve$eps_cost),c(0.025,0.5,0.975)) # WTP for statquo

# plot those distributions
# red is responders
# green is everyone
# blue is estimates
# vertical lines are medians, dashed line is for estimates

drawwtpdf <- as.data.frame(drawwtp)

estimatedwtp <- drawwtpdf %>% 
  mutate(wtp_ben = -mu_att/mu_cost) %>% 
  mutate(wtp_statquo = -mu_statquo/mu_cost) %>% 
  select(wtp_ben,wtp_statquo) %>% 
  mutate(Type="Estimate")

simulatedwtp <- singleserve %>% 
  mutate(wtp_ben=-(5+eps_att)/(-5+eps_cost)) %>% 
  mutate(wtp_statquo=-(1+eps_statquo)/(-5+eps_cost)) %>% 
  mutate(Type="Responders &\nNon-responders") %>% 
  select(wtp_ben,wtp_statquo,Type)

forplotting <- rbind(estimatedwtp,simulatedwtp)

####### WTP distributions (att) ######

ggplot() +
  geom_density(data=drawwtpdf,aes(x=-mu_att/mu_cost),fill="lightblue",color="lightblue",alpha=0.5) +
  geom_density(data=responders,aes(x=-1*(5+eps_att)/(-5+eps_cost)),fill="pink",color="pink",alpha=0.5) +
  geom_vline(xintercept=median(-1*(5+responders$eps_att)/(-5+responders$eps_cost)),color="red") +
  geom_vline(xintercept=median(-1*drawwtpdf$mu_att/drawwtpdf$mu_cost),color="blue") +
  geom_density(data=singleserve,aes(x=-1*(5+eps_att)/(-5+eps_cost)),color="#24a62f",alpha=0.5,size=1) +
  geom_vline(xintercept=median(-1*(5+singleserve$eps_att)/(-5+singleserve$eps_cost)),color="#24a62f",size=1) +
  geom_hline(aes(yintercept=0)) +
  lims(x=c(0.3,2.6)) +
  labs(x="",y="") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        text=element_text(size=11,  family="serif"),
        panel.grid=element_blank()) 


median(-drawwtpdf$mu_att/drawwtpdf$mu_cost)

####### WTP distributions (statquo) ######

ggplot() +
  geom_density(data=drawwtpdf,aes(x=-mu_statquo/mu_cost),fill="blue",alpha=0.5) +
  geom_density(data=singleserve,aes(x=-1*(-1+eps_statquo)/(-5+eps_cost)),fill="green",alpha=0.5) +
  geom_density(data=responders,aes(x=-1*(-1+eps_statquo)/(-5+eps_cost)),fill="red",alpha=0.5) +
  geom_vline(xintercept=mean(-1*(-1+singleserve$eps_statquo)/(-5+singleserve$eps_cost)),color="green") +
  geom_vline(xintercept=mean(-1*(-1+responders$eps_statquo)/(-5+responders$eps_cost)),color="red") +
  geom_vline(xintercept=mean(-1*drawwtpdf$mu_statquo/drawwtpdf$mu_cost),color="blue",linetype="dashed") +
  lims(x=c(-1.3,1.3))




##############################################
############## many small runs ###############
##############################################

res <- readRDS("~/selection-sims/100runsof2k.rds")
res <- res[2:nrow(res),]


wtpatt <- NULL
wtpstatquo <- NULL

for (i in 1:nrow(res)) {
  
  mat <- with(as.list(res[i,]),{
    vcov =  c(a1^2,  a1*b1,       a1*c1,             a1*d1,
              a1*b1, b1^2+b2^2,   b1*c1+b2*c2,       b1*d1+b2*d2,
              a1*c1, b1*c1+b2*c2, c1^2+c2^2+c3^2,    c1*d1+c2*d2+c3*d3,
              a1*d1, b1*d1+b2*d2, c1*d1+c2*d2+c3*d3, d1^2+d2^2+d3^2+d4^2)
    matrix(vcov,nrow=4)
  })
  
  drawwtp <- mvrnorm(n=50000,mu=unlist(res[i,c("mu_alpha","mu_statquo","mu_att","mu_cost")]),Sigma=mat)
  
  
  wtpatt <- rbind(wtpatt,quantile(-1*drawwtp[,3]/drawwtp[,4],c(0.025,0.5,0.975)))
  wtpstatquo <- rbind(wtpstatquo,quantile(-1*drawwtp[,2]/drawwtp[,4],c(0.025,0.5,0.975)))
  
}

wtpatt <- wtpatt %>% as.data.frame()
wtpstatquo <- wtpstatquo %>% as.data.frame()

ggplot(wtpatt) +
  geom_density(aes(x=`50%`)) +
  geom_vline(xintercept=mean(wtpatt$`50%`))

ggplot(wtpstatquo) +
  geom_density(aes(x=`50%`)) +
  geom_vline(xintercept=mean(wtpstatquo$`50%`))

ggplot(wtpatt) +
  geom_density(aes(x=`2.5%`)) +
  geom_vline(xintercept=mean(wtpatt$`2.5%`))

mean(wtpatt$`50%`)
median(wtpatt$`50%`)


mean(wtpstatquo$`50%`)
median(wtpstatquo$`50%`)

quantile(wtpatt$`50%`,c(0.025,.5,0.975))


# trudyselect <- database %>% select(id,iv) %>% unique()
# 
# totrudy <- clogitdata %>% select(-ends_with("_drop")) %>% 
#   select(id,alt,response_status,best=choice,att,cost,statquo) %>% arrange(id,alt)
# 
# totrudy <- left_join(totrudy,trudyselect)
# 
# write.csv(totrudy,"~/selection-sims/for_binarized_probit_estimation.csv")

dataf$response_status <- ifelse(dataf$response_status=="RESPONDER","Responders",dataf$response_status)
dataf$response_status <- ifelse(dataf$response_status=="NONRESPONDER","Non-responders",dataf$response_status)


p1 <- ggplot(dataf) +
  stat_density_ridges(aes(x=5+eps_att,y=response_status,fill=response_status,color=response_status),
                      alpha=1,show.legend = F,quantile_lines = TRUE, quantiles = 2) +
  labs(x="βbenefit",y="") +
  scale_fill_manual(values=c("#f0f037","#1ABDB7")) +
  scale_color_manual(values=c("#84841e","#0c5552")) +
  lims(x=c(-10,10))+ 
  theme_minimal()

p2 <- ggplot(dataf) +
  stat_density_ridges(aes(x=-5+eps_cost,y=response_status,fill=response_status,color=response_status),
                      alpha=1,show.legend = F,quantile_lines = TRUE, quantiles = 2) +
  labs(x="βcost",y="") +
  lims(x=c(-10,10))+
  scale_fill_manual(values=c("#f0f037","#1ABDB7")) +
  scale_color_manual(values=c("#84841e","#0c5552")) +
  theme_ridges(center_axis_labels = TRUE) +
  theme_minimal()

p3 <- ggplot(dataf) +
  stat_density_ridges(aes(x=-(5+eps_att)/(-5+eps_cost),y=response_status,fill=response_status,color=response_status),
                      alpha=1,show.legend = F,quantile_lines = TRUE, quantiles = 2) +
  labs(x="",y="") +
  scale_fill_manual(values=c("#f0f037","#1ABDB7")) +
  scale_color_manual(values=c("#84841e","#0c5552")) +
  theme_minimal() +
  lims(x=c(0,3))

  
p4 <- ggplot(forplotting) +
  stat_density_ridges(aes(x=wtp_ben,y=Type,fill=Type,color=Type),
                      quantile_lines = TRUE, quantiles = 2,show.legend = F) +
  scale_fill_manual(values=c("#af1abd","#85d25f")) +
  scale_color_manual(values=c("#460a4c","#355426")) +
  theme_minimal() +
  labs(x="WTP for benefit",y="") +
  lims(x=c(0,3))

p1/p2

p3/p4

ggsave("~/selection-sims/parameter_ridgeplot.png",plot=p1/p2,width=8,height=4,units="in")
ggsave("~/selection-sims/wtp_corrected_ridgeplot.png",plot=p3/p4,width=8,height=4,units="in")
