##############################################
## NOTES:
##  1.  Amat and Lmat have stepwise functions for sensitivity, because maturity uses a stepwise function of age or length
## CONVERSIONS:
##  1.  Steepness only depends upon SB_a and SR_alpha (not SR_beta)
##  2.  Increasing SR_beta decreases R0
## Updates:
## On 1/27 Cole added some code to test convergence of inner
## loop via line plots of yield curves. Writes this long
## data.frame to file for later plotting
## On 3/2 Cole added some really minor tweaks for plotting.
## On 6/11 Cole added code to run over a combination of M and LMARR values and
## plot a phase plane of elasticities.
## On 6/13 Cole restructured the code to make the minimization loops into a
## function, also reorganized some of the code
##############################################

## setwd("C:/Users/James.Thorson/Desktop/Project_git/Elasticity")
## Paths are now relative so make sure the wd is in the right folder
rm(list=ls())
File <- paste0(getwd(), "/")
library(ggplot2)
library(shape)    # has a better arrows function
library(reshape2)
source("functions.R")                   # core functions for analysis

## Settings for the example life histories
## SpeciesSet_Example = c("Rockfish","Snapper","Sardine")
SpeciesSet_Example= c("long-lived", "medium-lived", "short-lived")
M_Set_Example<- c(.05, .2 ,.6)
Linf_Set_Example<- c(45,90,30)
LMARR_Set_Example<- c(-1,2,1)
## what is steepness?
MLRR_Set_Example= exp(LMARR_Set_Example) / (1 - exp(-M_Set_Example))
round(MLRR_Set_Example / (4 + MLRR_Set_Example),2)

### The full sets across which we loop, used for the example life
### histories and other runs
F_Set_Full = c("Optimum", "Default")
## Parameter set across which to estimate elasticity
ParSet_Full <- c("M_a_All", "M_a_Juv", "M_a_Adult", "W_1", "W_2",
                 "SR_alpha", "SR_beta", "Amat", "Lmat", "Linf", "K", "h",
                 "R0")[-c(5, 7)]
## Sensitivity set is the amount to jitter each parameter
SensSet_Full <-
    c("Big-","Small-","Diff-","Default","Diff+","Small+","Big+")
ScenarioSet_Full = c("Orig","New","Wrong")


## ---------------
## Run the three example life histories
## ---------------
Results.examples <-
    run.scenario(M_Set=M_Set_Example, Linf_Set=Linf_Set_Example,
                 LMARR_Set=LMARR_Set_Example,
                 SensSet=SensSet_Full, F_Set=F_Set_Full,
                 ParSet=ParSet_Full, ScenarioSet=ScenarioSet_Full,
                 SpeciesSet=SpeciesSet_Example, save.curve.data=TRUE,
                 save.lifehist.curves=TRUE)
## The result object now contains the Results array, but also the age and
## yield curves as elements.
names(Results.examples)
## Save the whole run to file
save(Results.examples, file=paste0("Results/Results.examples.RData"))

## ---------------
## Process the results
## ---------------
load("Results/Results.examples.RData")
## Eliminate stuff that crashes the population
Results.examples$Results[1,1,1,7,'New',c('Fmsy','SBmsy','Cmsy')] = 0
## Calculate elasticities and manipulate data for plotting
Diff.examples <- with(Results.examples,
                      log(Results['Optimum',,ParOrder,"Diff+",'New',]) -
                      log(Results['Optimum',,ParOrder,"Diff-",'New',]))
## Create elasticities
Elas.examples = Diff.examples / (2 * Delta)
Elas.table <- rbind(Elas.examples[,,'Fmsy'], Elas.examples[,,'SBmsy'], Elas.examples[,,'Cmsy'])
## Combine the biological parameters together
scalars.long <- melt(lapply(Results.examples$Unfished, function(x) x$scalar.params))
names(scalars.long) <- c("value", "parameter", "species")
scalars.wide <- dcast(scalars.long, parameter~species, value.var="value", fun.aggregate=median)
par.order <- c("Linf", "LMARR", "M", "Selex_A50_times_Amat", "Selex_Slope",
               "W0", "theta1", "theta2", "R0", "Winf", "B", "A", "L0", "K", "Amat", "Lmat", "W_amat",
               "MLRR", "h", "SB_0_hat","SR_alpha", "SR_beta")
##               "Fmsy", "Cmsy", "SBmsy", "Rmsy",
scalars.wide.table <- subset(scalars.wide, parameter %in% par.order)
## Terribly ugly way to make a table. Sort it so the table matches paper easier
scalars.wide.table <- scalars.wide.table[match(par.order, scalars.wide.table$parameter),]
temp <- data.frame(t(Results.examples$Results["Optimum",,1, "Default", "New",]))
temp <- cbind(parameter=row.names(temp), temp)
names(temp) <- names(scalars.wide.table)
scalars.wide.table <- rbind(scalars.wide.table, temp)
row.names(scalars.wide.table) <- scalars.wide.table$parameter
scalars.wide.table$parameter <- NULL

## ---------------
## Write Tables
## ---------------
write.table( format(round(Elas.table,2), digits=2, nsmall=2, scientific=FALSE),
            file=paste0(File,"Results/table2.csv"), sep=",")
with(Results.examples, write.csv("Results/age.curves.csv", x=age.curves))
with(Results.examples, write.csv("Results/yield.curves.csv", x=yield.curves))
write.csv("Results/biology_scalars.csv", x=scalars.wide)
write.csv("Results/Table1.csv", x= scalars.wide.table)

## ---------------
## Make figures and plots
## ---------------
## This file makes some exploratory ggplots for checking convergence and
## biological properties of the life histories, takes ~ 1 minute
#source("make_ggplots.R")
## Make a subset of data for plotting, since no need to plot the Diff ones,
## which are far too small to see the effect
PlotResults <-
    Results.examples$Results[,,,c("Big-","Small-","Default","Small+","Big+"),,]
## Some global settings for tweaking and prettyification
species.ltys <- c(1,3, 1)
border.col <- gray(.5)
species.cols <- c(gray(0), gray(0), gray(.5))
species.pchs <- c(16, 16, 16)
label.col <- gray(.3)
## When run, these functions make the figure given the file inputs
res <- 500
## This file has functions for creating the manuscript figures, but also
## produces some exploratory ggplots
source("figures.R")
figure1("png", width=3, height=5, res=res)
figure2a("Fmsy", "png", width=6.5, height=4.5, res=res)
figure2b("Fmsy", "png", width=6.5, height=4.5, res=res)
figure2c("Fmsy", "png", width=6.5, height=4.5, res=res)
figure3("png", width=6.5, height=4.5, res=res)
figure4("png", width=6.5, height=4.5, res=res)
figure_new("png", width=6.5, height=4.5, res=res)
## End of code
## ------------------------------------------------------------


## Everything below here is old or interim code for exploration. Doesn't
## need to be run for the manuscript results.

## ------------------------------------------------------------
## This is some old code when we tried to look at plotting a contour
## surface. For now we are abandoning this approach but leaving the code
## here in case we want to come back to it at some point..
## ## Make a grid of M and LMARR for contour plotting
## n.grid <- 20                         # length of each dimension
## ## number of example species, i.e. those that are not for the grid and we
## ## want to save more information about for plotting
## M.seq <- seq(.05, 1, len=n.grid)
## LMARR.seq <- seq(-1,2, len=n.grid)
## grid <- expand.grid(M=M.seq, LMARR=LMARR.seq, Linf=60)
## grid$Species <- paste0("grid", row.names(grid))
## ## Turn off the extra data calculations since don't need them and they are
## ## slow, as well as combinations that aren't used.
## Results.grid <-
##     run.scenario(M_Set= grid$M, Linf_Set=grid$Linf, LMARR_Set=grid$LMARR,
##                  SpeciesSet=grid$Species, SensSet=c("Diff+","Default", "Diff-"),
##                  F_Set=F_Set_Full, ParSet=ParSet_Full, ScenarioSet=ScenarioSet_Full,
##                  save.curve.data=FALSE, save.lifehist.curves=FALSE)
## ## Save the whole run to file
## save(Results.grid, file=paste0("Results/Results.grid.RData"))
## load("Results/Results.grid.RData")
## Diff.grid = with(Results.grid, log(Results['Optimum',,ParOrder,"Diff+",'New',]) -
##     log(Results['Optimum',,ParOrder,"Diff-",'New',]))
## Elas.grid = Diff.grid / (2 * Delta)
## Elas.grid.long <- reshape2::melt(Elas.grid)
## names(Elas.grid.long) <- c("Species", "Parameter", "Metric", "Elasticity")
## Elas.grid.long <- merge(Elas.grid.long, grid)


## ## ---------------
## ## Run some explicit examples to look at convergence and behavior
## Results.topleft <-
##     run.scenario(M_Set= .8, Linf_Set=60, LMARR_Set=-1,
##                  SpeciesSet="topleft", SensSet=SensSet_Full,
##                  F_Set=F_Set_Full, ParSet=ParSet_Full, ScenarioSet=ScenarioSet_Full,
##                  save.curve.data=TRUE, save.lifehist.curves=TRUE)
## save(Results.topleft, file=paste0("Results/Results.topleft.RData"))
## ## Plot yield curves
## g <- ggplot(Results.topleft$yield.curves, aes(x=F, y=SSB))
## g <- g +  geom_line(aes(colour=Multiplier, group=Delta)) +
##     geom_point(aes(x=Fmsy, y=SBmsy, colour=Multiplier)) +
##     facet_grid(Species~Parameter, scales="free_y")
## ggsave(filename="Figures/topleft_yield_curves.png", g, width=12, height=6,
##        units="in")
## ## Melt down the age
## topleft.age.long <- melt(Results.topleft$age.curves,
##                        measure.vars=c("B_a_hat", "Catch_a", "W_a_hat", "Surv_a_hat", "BPR_a_hat"))
## g <- ggplot(topleft.age.long, aes(x=age, y=value, colour=Multiplier))
## g <- g +  geom_line(aes(group=ID))  + xlim(0,15)+
##     facet_grid(variable~Parameter, scales="free_y")
## ggsave(filename="Figures/topleft_age_curves_weight.png", g, width=12, height=6,
##        units="in")
## ## Run an example on the crazy ridge
## temp <- Elas.grid.long[which.max(abs(Elas.grid.long$Elasticity)),]
## Results.ridge <-
##     run.scenario(M_Set= temp$M, Linf_Set=temp$Linf, LMARR_Set=temp$LMARR,
##                  SpeciesSet="ridge", SensSet=SensSet_Full,
##                  F_Set=F_Set_Full, ParSet=ParSet_Full, ScenarioSet=ScenarioSet_Full,
##                  save.curve.data=TRUE, save.lifehist.curves=TRUE)
## save(Results.ridge, file=paste0("Results/Results.ridge.RData"))
## ## Plot yield curves
## g <- ggplot(Results.ridge$yield.curves, aes(x=F, y=SSB))
## g <- g +  geom_line(aes(colour=Multiplier, group=Delta)) +
##     geom_point(aes(x=Fmsy, y=SBmsy, colour=Multiplier)) +
##     facet_grid(Species~Parameter, scales="free_y")
## ggsave(filename="Figures/ridge_yield_curves.png", g, width=12, height=6,
##        units="in")
## ## Melt down the age
## ridge.age.long <- melt(Results.ridge$age.curves,
##                        measure.vars=c("B_a_hat", "Catch_a", "W_a_hat", "Surv_a_hat", "BPR_a_hat"))
## g <- ggplot(ridge.age.long, aes(x=age, y=value, colour=Multiplier))
## g <- g +  geom_line(aes(group=ID))  + xlim(0,15)+
##     facet_grid(variable~Parameter, scales="free_y")
## ggsave(filename="Figures/ridge_age_curves_weight.png", g, width=12, height=6,
##        units="in")


## ## RUn a test to see the influence of Linf while holding M and LMARR
## ## constant
## Results.Linf <-
##     run.scenario(M_Set=c(.3, .3), Linf_Set=c(30,60),
##                  LMARR_Set=c(1, 1),
##                  SensSet=SensSet_Full, F_Set=F_Set_Full,
##                  ParSet=ParSet_Full, ScenarioSet=ScenarioSet_Full,
##                  SpeciesSet=c("short", "long"), save.curve.data=TRUE,
##                  save.lifehist.curves=TRUE)
## g <- ggplot(Results.Linf$yield.curves, aes(x=F, y=Catch, colour=Multiplier))
## g <- g +  geom_line(aes(group=Multiplier)) +
##     geom_point(aes(x=Fmsy, y=Cmsy)) +
##     facet_grid(Species~Parameter, scales="free_y")
## ggsave(filename="Figures/Linf_example_F_curves_catch.png",g, width=12, height=6,
##        units="in")
## g <- ggplot(Results.Linf$age.curves, aes(x=age, y=W_a_hat, colour=Multiplier))
## g <- g +  geom_line(aes(group=ID))  + xlim(0,50)+
##     facet_grid(Species~Parameter, scales="free_y")
## ggsave(filename="Figures/Linf_example_age_curves_weight.png", g, width=12, height=6,
##        units="in")
## ## How do elasticities compare?
## Linf_Diff <- with(Results.Linf, log(Results['Optimum',,ParOrder,"Diff+",'New',]) - log(Results['Optimum',,ParOrder,"Diff-",'New',]))
## Linf_Elas <- Linf_Diff / (2 * Delta)
## do.call(rbind, as.list(Linf_Elas))
## Linf_Elas_long <- melt(Linf_Elas)
## names(Linf_Elas_long) <- c("Linf", "variable", "metric", "elasticity")
## Linf_Elas_wide <- dcast(Linf_Elas_long, formula=variable+metric~Linf)
## write.csv(Linf_Elas_wide, file="Results/Linf_example_elasticities.csv")
##
## Write the elasticities to file
## for(SpeciesI in 1:length(SpeciesSet)){
##     Diff = log(Results['Optimum',,ParOrder,"Diff+",'New',]) - log(Results['Optimum',,ParOrder,"Diff-",'New',])
##     ##Response = log(Results['Optimum',,,c("Diff-","Default","Diff+"),'New',]) - log(Results['Optimum',,,c("Diff-","Default","Diff+"),'Orig',])
##     ##Elas = ( log(Response[,,'Diff+',])-log(Response[,,'Diff-',]) ) / 2
##     Elas = Diff / (2 * Delta)
##     Mat = rbind(Elas[,,'Fmsy'], Elas[,,'SBmsy'], Elas[,,'Cmsy'])
##     Elas.long <- reshape2::melt(Elas)
##     names(Elas.long) <- c("Species", "Parameter", "Metric", "Elasticity")
##     Elas.long <- merge(Elas.long, grid)
##     write.csv( Mat, file=paste(File,"Results/Elasticity_Species=",SpeciesSet[SpeciesI],".csv",sep=""))
## }

## max(Elas.grid.long$Elasticity, na.rm=TRUE)
## df.temp <- data.frame(M=M_Set_Example, LMARR=LMARR_Set_Example)
## elas.max <- 20
## Elas.grid.long$elas.trunc <- Elas.grid.long$Elasticity
## Elas.grid.long$elas.trunc <-
##     ifelse(Elas.grid.long$elas.trunc>elas.max, elas.max, Elas.grid.long$elas.trunc)
## Elas.grid.long$elas.trunc <-
##     ifelse(Elas.grid.long$elas.trunc < -elas.max, -elas.max, Elas.grid.long$elas.trunc)
## g <- ggplot(Elas.grid.long, aes(x=LMARR, y=M))+
##     facet_grid(Metric~Parameter) +
##     geom_point(aes(size=sqrt(abs(elas.trunc)), alpha=sqrt(abs(elas.trunc)), colour=elas.trunc>0)) +
##     geom_polygon(data=df.temp, aes(LMARR, M), fill=NA, linetype=1,
##                  colour=gray(.8)) +
##     geom_point(data=df.temp, aes(LMARR, M), shape=1) +
##     scale_size_continuous(limits = c(0,elas.max))
## ggsave("Figures/elasticity_surfaces.png",g,  width=12, height=6)


## End of file
## ------------------------------------------------------------
