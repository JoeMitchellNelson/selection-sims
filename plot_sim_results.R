require(pacman)
p_load(pacman,tidyverse,MASS,evd,sampleSelection,foreach,boot,
       doParallel,tictoc,patchwork,matrixcalc,survival,plotly,apollo,ggridges)

#############################
# num #  beta21  #  alpha2  #
#  1  #    2     #     2    #
#  2  #    2     #    -2    #
#  3  #   -2     #     2    #
#  4  #   -2     #    -2    #
#############################

res1 <- read.csv("~/selection-sims/simresults1-pos-pos.csv")
res2 <- read.csv("~/selection-sims/simresults2-pos-neg.csv")
res3 <- read.csv("~/selection-sims/simresults3-neg-pos.csv")
res4 <- read.csv("~/selection-sims/simresults4-neg-neg.csv")

{
res1plot <- res1 %>% dplyr::select(mu_x,
                               uncorrected_wtp,
                               fullsample_wtp,
                               adhoc_wtp,
                               maxEigen) %>% 
  pivot_longer(names_to="est",cols=c(mu_x,
               uncorrected_wtp,
               fullsample_wtp,
               adhoc_wtp),values_to="WTP")

res2plot <- res2 %>% dplyr::select(mu_x,
                                   uncorrected_wtp,
                                   fullsample_wtp,
                                   adhoc_wtp,
                                   maxEigen) %>% 
  pivot_longer(names_to="est",cols=c(mu_x,
                                     uncorrected_wtp,
                                     fullsample_wtp,
                                     adhoc_wtp),values_to="WTP")

res3plot <- res3 %>% dplyr::select(mu_x,
                                   uncorrected_wtp,
                                   fullsample_wtp,
                                   adhoc_wtp,
                                   maxEigen) %>% 
  pivot_longer(names_to="est",cols=c(mu_x,
                                     uncorrected_wtp,
                                     fullsample_wtp,
                                     adhoc_wtp),values_to="WTP")

res4plot <- res4 %>% dplyr::select(mu_x,
                                   uncorrected_wtp,
                                   fullsample_wtp,
                                   adhoc_wtp,
                                   maxEigen) %>% 
  pivot_longer(names_to="est",cols=c(mu_x,
                                     uncorrected_wtp,
                                     fullsample_wtp,
                                     adhoc_wtp),values_to="WTP")
}

cols <- c("mu_x" = "#f79914", "fullsample_wtp" = "grey50", "adhoc_wtp" = "#2e14f7", "uncorrected_wtp" = "#14c6f7")
key_labs <- c("mu_x" = "Mixl corrected", "fullsample_wtp" = "True DGP", "adhoc_wtp" = "Ad hoc corrected", "uncorrected_wtp" = "Naive clogit")
var_order <- c("fullsample_wtp","mu_x","adhoc_wtp","uncorrected_wtp")

makedplot <- function (dfname) { 
  ggplot(dfname[which(dfname$maxEigen<=0),]) +
  geom_density(aes(x=WTP,group=est,fill=est),color=NA,alpha=0.5) +
  scale_fill_manual(breaks=var_order,
                    values=cols,labels=key_labs) +
  geom_vline(xintercept=2) +
  labs(x="WTP Estimate",y="Density",fill="") +
  scale_x_continuous(breaks=seq(from=1,to=2.5,by=.5),limits=c(1,2.5)) +
    scale_y_continuous(breaks=NULL,limits=c(0,10)) +
    
  theme_minimal()
}

p1 <- makedplot(res1plot) + labs(title=expression(alpha[2]*" = 2, "*beta[21]*" = 2"))
p2 <- makedplot(res2plot) + labs(title=expression(alpha[2]*" = -2, "*beta[21]*" = 2"))
p3 <- makedplot(res3plot) + labs(title=expression(alpha[2]*" = 2, "*beta[21]*" = -2"))
p4 <- makedplot(res4plot) + labs(title=expression(alpha[2]*" = -2, "*beta[21]*" = -2"))
 
combinedplot <- p1 + p2 + p3 + p4 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
combinedplot[[1]] <- combinedplot[[1]] + theme(axis.title.y = element_blank(),
                                       axis.title.x = element_blank() )
combinedplot[[2]] <- combinedplot[[2]] + theme(
                                               axis.title.x = element_blank(),
                                               axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.title.y = element_blank())
combinedplot[[3]] <- combinedplot[[3]] + theme(axis.title.y = element_blank())
combinedplot[[4]] <- combinedplot[[4]] + theme(axis.text.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.title.y = element_blank() )

combinedplot

############### ridge density ###############

makerplot <- function (dfname) {
  ggplot(dfname[which(dfname$maxEigen<=0),]) +
  geom_density_ridges(aes(x=WTP,y=est,fill=est),color="white",alpha=0.5,scale=1.5,show.legend=F) +
  scale_fill_manual(breaks=rev(var_order),
                    values=cols,labels=key_labs) +
  geom_vline(xintercept=2) +
  labs(x="WTP Estimate",y="",fill="") +
  scale_x_continuous(breaks=seq(from=1,to=2.5,by=.5),limits=c(1,2.5)) +
  scale_y_discrete(limits=var_order,expand=c(0,0),labels=key_labs) +
   theme_minimal() +
    theme(panel.grid.major.y = element_blank())
}


p1 <- makerplot(res1plot) + labs(title=expression(alpha[2]*" = 2, "*beta[21]*" = 2"))
p2 <- makerplot(res2plot) + labs(title=expression(alpha[2]*" = -2, "*beta[21]*" = 2"))
p3 <- makerplot(res3plot) + labs(title=expression(alpha[2]*" = 2, "*beta[21]*" = -2"))
p4 <- makerplot(res4plot) + labs(title=expression(alpha[2]*" = -2, "*beta[21]*" = -2"))

combinedplot <- (p1 + p2)/(p3 + p4)
combinedplot

ggsave("~/selection-sims/ridgeplot.png",plot=combinedplot,width=8,height=4,units="in")
