setwd("/home/federico/Documents/uni/dottorato_pd_2010/papers_for_thesis/parfm/graphs")

install.packages("parfm")

R.Version()[["version.string"]]

library("parfm")
packageDescription("parfm", fields="Version")

data("kidney")

head(kidney)

kidney$sex <- kidney$sex - 1

############################## Model estimation #################################
mod <- parfm(Surv(time, status) ~ sex + age, cluster="id", 
             data=kidney, dist="exponential", frailty="gamma")
mod
tau(mod)

ci.parfm(mod)["sex",]

############################# Frailty prediction#################################
u <- predict(mod)

pdf("prediplot.pdf", height=6, width=8)
plot(u, sort="i", main="", cex.axis=.7)
dev.off()


####################### Comparison of different models ##########################
kidney.parfm <- select.parfm(Surv(time, status) ~ sex + age, 
    cluster="id", data=kidney,
    dist=c("exponential", "weibull", "gompertz",
           "loglogistic", "lognormal"),
    frailty=c("gamma", "ingau", "possta"))
kidney.parfm

pdf("plot.pdf", height=4, width=8)
  plot(kidney.parfm)
dev.off()

parfm(Surv(time, status) ~ sex + age, cluster="id", 
      data=kidney, dist="exponential", frailty="ingau")

parfm(Surv(time, status) ~ sex + age, cluster="id", 
      data=kidney, dist="exponential", frailty="possta")

parfm(Surv(time, status) ~ sex + age, cluster="id", 
        data=kidney, dist="exponential", frailty="possta", 
        iniFpar=0.25)

parfm(Surv(time, status) ~ sex + age, cluster="id", 
        data=kidney, dist="exponential", frailty="possta", 
        method="Nelder-Mead")

coxph(Surv(time, status) ~ sex + age +
        frailty(id, distribution="gamma", eps=1e-11),
        outer.max=15, data=kidney)

