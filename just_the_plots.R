require(pacman)
p_load(tidyverse,plotly)

gridsearch1 <- readRDS("sim_results_no_causality.rds")
gridsearch2 <- readRDS("sim_results_with_causality.rds")

# latent problem-variable is x
# instrument is z
# "bias_remain" of 1 indicates no improvement over uncorrected estimate
# 0 indicates bias completely removed
# >1 indicates "correction" introduced additional bias

##################################################
### assumes x plays no causal role in selection ##
##################################################

fig1 <- plot_ly(gridsearch1, x = ~cov1, y = ~cov2, z = ~covz,color=~bias_remain)
fig1 <- fig1 %>% add_markers()
fig1 <- fig1 %>% layout(scene = list(xaxis = list(title = 'Cov(beta1,z)'),
                                   yaxis = list(title = 'Cov(beta2,z)'),
                                   zaxis = list(title = 'Cov(x,z)')))

fig1

###################################################################################
### assumes x and z play equally important causal roles in selection ##############
### some bias_remain values are very large, so bias measure is truncated at 2 #####
###################################################################################

fig2 <- plot_ly(gridsearch2, x = ~cov1, y = ~cov2, z = ~covz,color=~bias_remain_truncated)
fig2 <- fig2 %>% add_markers()
fig2 <- fig2 %>% layout(scene = list(xaxis = list(title = 'Cov(beta1,z)'),
                                     yaxis = list(title = 'Cov(beta2,z)'),
                                     zaxis = list(title = 'Cov(x,z)')))

fig2
