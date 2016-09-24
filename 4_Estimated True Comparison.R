
#------ Check model to make sure it retrieves the parameter estimates

#--- Define parameter values
# Uninfected survival probability

true <- c(phi_U, 
          phi_I, 
          p_U,
          p_I,
          beta_UI,
          beta_IU,
          gamma_U, 
          gamma_U, 
          gamma_U, 
          gamma_I,
          gamma_I,
          gamma_I)

names <- c("phi_U", 
          "phi_I", 
           "p_U",
          "p_I",
           "beta_UI",
          "beta_IU",
           "gamma_U1", 
          "gamma_U2", 
           "gamma_U3", 
           "gamma_I1",
          "gamma_I2",
           "gamma_I3")

mod.mean <- c(
              js.ms$mean$phi_U,
              js.ms$mean$phi_I,
              js.ms$mean$p_U,
              js.ms$mean$p_I,
              js.ms$mean$beta_UI,
              js.ms$mean$beta_IU,
              js.ms$mean$gamma_U[1:3],
              js.ms$mean$gamma_I[1:3])

mod.q2.5 <- c(
  js.ms$q2.5$phi_U,
  js.ms$q2.5$phi_I,
  js.ms$q2.5$p_U,
  js.ms$q2.5$p_I,
  js.ms$q2.5$beta_UI,
  js.ms$q2.5$beta_IU,
  js.ms$q2.5$gamma_U[1:3],
  js.ms$q2.5$gamma_I[1:3])

mod.q97.5 <- c(
  js.ms$q97.5$phi_U,
  js.ms$q97.5$phi_I,
  js.ms$q97.5$p_U,
  js.ms$q97.5$p_I,
  js.ms$q97.5$beta_UI,
  js.ms$q97.5$beta_IU,
  js.ms$q97.5$gamma_U[1:3],
  js.ms$q97.5$gamma_I[1:3])

dat <- data.frame(names = names, true = true, mod.mean = mod.mean, mod.q2.5 = mod.q2.5, mod.q97.5 = mod.q97.5)


library(ggplot2)

cols <- c("Truth" = "red", "Estimated" = "black")

ggplot(dat, aes(x= names, y=mod.mean, ymin=mod.q2.5, ymax=mod.q97.5))+ 
  geom_linerange(size = 1) +
  geom_point(size = 3, aes(x = names, y = mod.mean, col = "Estimated")) +
  geom_point(size = 3, aes(x = names, y = true, col = "Truth")) +
  scale_colour_manual("Values", values=cols)+
  geom_hline(yintercept = 0, lty=2) +
  coord_flip() + ylab('Parameter estimates') +
  xlab("Parameter names") +
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size = 17, color = "black"), 
        axis.title.y = element_text(size = 17, color = "black"), 
        axis.title.x =element_text(size = 17, color = "black"),
        legend.title =element_text(size = 17, color = "black"),
        legend.text =element_text(size = 17, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 





