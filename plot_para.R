library(tidyverse)
plot_para <- data.frame(na.omit(opt_params))
colnames(plot_para) <- c("group","step","Q","phi","gam","mu_g","mu_f","rho_f","V","mu_rho","sig_rho")

ggplot(plot_para %>% gather(key,value,-group,-step),aes(x=step,y=value,color=as.factor(group))) +
  geom_line() + facet_wrap(.~key,scale="free")

