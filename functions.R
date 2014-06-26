## testing stuff
## Opt <- optimize( f=ProjectionFn, F=.1, interval=Interval,
##                 OutputType="Catch", M_a=M_a, W_1= W_1,
##                 W_2=W_2, SR_alpha=SR_alpha, SR_beta=SR_beta, Amat=Amat,
##                 Lmat=Lmat, Linf=Linf, K=K, h=h, R0=R0, Selex_a=Selex_a )
## Opt_Orig <- optimize( f=EquilibriumFn, interval=c(0,1), maximum=TRUE,
##                      Output="Catch", M_a=M_a, W_1=W_1, W_2=W_2,
##                      SR_alpha=SR_alpha, SR_beta=SR_beta, Amat=Amat,
##                      Lmat=Lmat, Linf=Linf, K=K, h=h, R0=R0,
##                      Selex_a=Selex_a)

## Global Settings
Selex_A50_times_Amat = 1                                   ################# THIS HAS BEEN CHANGED
Selex_Slope = 10000
## Scaling parameters
W_1 = 0.01 # Fusiform prior given weight in gram and length in cm
W_2 = 3.00 # Necessary for derivations used below (although fusiform prior is 3.04 given weight in gram and length in cm)
## W0 was changed on 6/23/14 and moved inside the run.scenario function
## Need to set a global here since R0 is overwritten in some cases inside
## the loop
R0_global <- 1e9
Interval = log( c(1.01, R0_global*2) )
AgeSet = 0:100
M_Type = c("Constant","Length-Based")[2]
Delta = 1e-6
ParOrder <- c("M_a_All","M_a_Juv","M_a_Adult",  "h","SR_alpha","R0",
              "Linf","K","W_1",  "Lmat","Amat")
## Names for plotting
ParName <- c("All ages", "Juvenile", "Adult", "Steepness", "MRPSB",
             expression(R[0]), expression(L[inf]), "k", expression(theta[1]),
             expression(L[mat]), expression(A[mat]))
Multiplier.values <- data.frame("Big-"=exp(-0.20), "Small-"=exp(-0.1),
                                "Default"=1, "Small+"=exp(0.1), "Big+"=exp(0.20))

LambdaFn <- function(M_a, Mat_a, W_a, SR_alpha){
  Trans = matrix(0, ncol=length(AgeSet), nrow=length(AgeSet))
  for(AgeI in 2:length(AgeSet) ) Trans[AgeI,AgeI-1] = exp(-M_a[AgeI-1])
  for(AgeI in 1:length(AgeSet) ) Trans[1,AgeI] = Mat_a[AgeI] * W_a[AgeI] * SR_alpha
  Eigen = eigen(Trans)
  if( Im(Eigen$value[1])!=0 ){
    Return = NA
  }else{
    Return = Re(Eigen$value[1])
  }
  return( Return )
}
## Projection function
ProjectionFn <- function(F, Ln_Rinput=NA, OutputType="Catch", M_a, W_1,
                         W_2, SR_alpha, SR_beta, Amat, Lmat, Linf, K, h,
                         R0, Selex_a, L0, W0){
    if(!is.na(Lmat)){
        Amat = -log( (Linf-Lmat)/(Linf-L0) ) / K
    }
    ## 95% of maturity occurs within Amat +/- 25% of Amat
    Mat_a_hat = pnorm( (AgeSet-Amat)/(Amat/4)*1.96 )
    L_a_hat = Linf - (Linf - L0) * exp(-K * AgeSet)
    W_a_hat = W_1 * L_a_hat^W_2
    Surv_a_hat = exp(-cumsum((M_a+F*Selex_a)))
    BPR_a_hat = Surv_a_hat * W_a_hat
    SBPR_a_hat = BPR_a_hat * Mat_a_hat
    B_a_hat = exp(Ln_Rinput) * BPR_a_hat
    SB_a_hat = exp(Ln_Rinput) * SBPR_a_hat
    SB_average = sum(SB_a_hat)
    SBPR_0_hat = sum( exp(-cumsum(M_a)) * W_a_hat * Mat_a_hat )
    if(!is.na(h) & !is.na(R0)){
        SB_0_hat = R0 * SBPR_0_hat
        MLRR =  4*h / (1-h)
        SR_alpha = MLRR / sum(SBPR_0_hat)
        SR_beta = ( (SR_alpha * SB_0_hat) / R0 - 1) / SB_0_hat
        SR_curve <- t(sapply(seq(0,SB_0_hat,length=101),
                             FUN=function(Num){c(Num,SR_alpha*Num/(1+SR_beta*Num))}))
    } else {
        SB_0_hat = NULL
        SR_curve = NULL
    }
    R_average = SR_alpha * SB_average / (1 + SR_beta*SB_average)
    CPR_a <- Surv_a_hat * ((Selex_a*F)/(M_a+Selex_a*F)) *
        (1 - exp(-(M_a+F*Selex_a))) * W_a_hat
    C_a = exp(Ln_Rinput) * CPR_a
    Recruitment_Objective = (exp(Ln_Rinput)/R_average - 1)^2
    if(OutputType=="Recruitment_Objective") Return = Recruitment_Objective
    if(OutputType=="Catch") Return = sum(C_a)
    if(OutputType=="SB_average") Return = SB_average
    if(OutputType=="Diag")
        Return = list( 'Recruitment_Objective'=Recruitment_Objective,
        'L_a_hat'=L_a_hat, 'W_a_hat'=W_a_hat, 'Surv_a_hat'=Surv_a_hat,
        'BPR_a_hat'=BPR_a_hat, 'SBPR_a_hat'=SBPR_a_hat, 'B_a_hat'=B_a_hat,
        'Mat_a_hat'=Mat_a_hat, 'SB_a_hat'=SB_a_hat,
        'SB_average'=SB_average, 'SB_0_hat'=SB_0_hat,
        'R_average'=R_average, 'Rinput'=exp(Ln_Rinput), 'Catch_a'=C_a,
        'Catch'=sum(C_a), 'SR_alpha'=SR_alpha, 'SR_beta'=SR_beta,
        'SR_curve'=SR_curve )
    if(OutputType=="All"){
        ## Return all biological and derived params
        Return <-
            list(age.params=list( 'Recruitment_Objective'=Recruitment_Objective,
                 ## derived params
                 'L_a_hat'=L_a_hat, 'W_a_hat'=W_a_hat, 'Surv_a_hat'=Surv_a_hat,
                 'BPR_a_hat'=BPR_a_hat, 'SBPR_a_hat'=SBPR_a_hat, 'B_a_hat'=B_a_hat,
                 'Mat_a_hat'=Mat_a_hat, 'SB_a_hat'=SB_a_hat,
                 'Selex_a'=Selex_a, 'Catch_a'=C_a, 'M_a'=M_a),
             scalar.params=list('SB_average'=SB_average, 'SB_0_hat'=SB_0_hat,
                 'R_average'=R_average, 'Rinput'=exp(Ln_Rinput),
                 'Catch'=sum(C_a), 'SR_alpha'=SR_alpha, 'SR_beta'=SR_beta,
                 'W_1'=W_1, 'W_2'=W_2, 'SR_alpha'=SR_alpha, 'SR_beta'=SR_beta,
                 'Amat'=Amat, 'Lmat'=Lmat, 'Linf'=Linf, 'K'=K, 'h'=h, 'R0'=R0,
                 'L0'=L0, 'W0'=W0))
    }
    return(Return)
}

## M_a=M_a_Input; W_1=W_1_Input; W_2=W_2_Input; SR_alpha=SR_alpha_Input; SR_beta=SR_beta_Input; Amat=Amat_Input; Lmat=Lmat_Input; Linf=Linf_Input; K=K_Input
EquilibriumFn <- function(F, Output="Catch", M_a, W_1, W_2, SR_alpha,
                          SR_beta, Amat, Lmat, Linf, K, h, R0, Selex_a, L0,
                          W0){
    Opt <- optimize( f=ProjectionFn, F=F, interval=Interval,
                    OutputType="Recruitment_Objective", M_a=M_a, W_1= W_1,
                    W_2=W_2, SR_alpha=SR_alpha, SR_beta=SR_beta, Amat=Amat,
                    Lmat=Lmat, Linf=Linf, K=K, h=h, R0=R0, Selex_a=Selex_a,
                    L0=L0, W0=W0)
    ## ProjectionFn(F=F, Ln_Rinput=Interval[1], OutputType="Diag", M_a, W_1, W_2, SR_alpha, SR_beta, Amat, Lmat, Linf, K, h, R0 )
    if(Output=="Catch") Return = ProjectionFn(F=F, Ln_Rinput=Opt$minimum,
                        OutputType="Catch", M_a, W_1, W_2, SR_alpha,
                        SR_beta, Amat, Lmat, Linf, K, h, R0,
                        Selex_a=Selex_a, L0=L0, W0=W0)
    if(Output=="SB_average") Return = ProjectionFn(F=F,
                             Ln_Rinput=Opt$minimum,
                             OutputType="SB_average", M_a, W_1, W_2,
                             SR_alpha, SR_beta, Amat, Lmat, Linf, K, h, R0,
                             Selex_a=Selex_a, L0=L0, W0=W0)
    if(Output=="Diag"){
        Return <- ProjectionFn(F=F, Ln_Rinput=Opt$minimum,
                               OutputType="Diag", M_a, W_1, W_2, SR_alpha,
                               SR_beta, Amat, Lmat, Linf, K, h, R0,
                               Selex_a=Selex_a,L0=L0, W0=W0)
        Return[['Opt']] = Opt
    }
    if(Output=="Recruitment_Objective") Return = Opt$objective
    return(Return)
}

run.scenario <- function(M_Set, Linf_Set, LMARR_Set, SensSet, F_Set,
                         ParSet, ScenarioSet, SpeciesSet,
                         save.curve.data, save.lifehist.curves){
    ## This function runs an optimization scenario across the sets of
    ## inputs given. It's broken down like this since in some cases we want
    ## to run only the elasticities (for the grid) while for the example
    ## life histories we want the full information and age and yield curves
    ## too. The save arguemnts determine whether to run and return those
    ## extra data.
    ##
    ## Results objects
    Results <- array(NA, dim=c(length(F_Set),length(SpeciesSet),length(ParSet),length(SensSet),length(ScenarioSet),3),
                     dimnames=list(F_Set,SpeciesSet,ParSet,SensSet,ScenarioSet,c("Fmsy","SBmsy","Cmsy")))
    ## These lists hold temporary results for the age and yield curves
    age.list <- yields.list <- list()
    age.curves <- yield.curves <- NULL
    ## Thie lifehistory list is used to plot figure 1, which shows the
    ## basic properties of the model given a life history.
    lifehist.curves <- list()
    Unfished <- list()
    k <- 1                              # used in loop below
for(SpeciesI in 1:length(SpeciesSet)){
for(F_I in 1:length(F_Set)){
for(ParI in 1:length(ParSet)){
for(SensI in c(ceiling(length(SensSet)/2),(1:length(SensSet))[-ceiling(length(SensSet)/2)])){
    F_Type = F_Set[F_I]
    M = M_Set[SpeciesI]   # Scales catabolism
    Linf = Linf_Set[SpeciesI]  # Scales anabolism
    LMARR = LMARR_Set[SpeciesI] # Red snapper, clupeidae from Myers et al. 1999
    print(paste(collapse = ", ", c(M, Linf, LMARR, F_Set[F_I], SpeciesSet[SpeciesI],
                                   ParSet[ParI], SensSet[SensI])))
    ## print(paste(collapse = ", ", c(F_I, SpeciesI, ParI, SensI)))
    ## Derived parameters
    B = M / (3^(1/2)/2^(3/2))     # M = (3^(1/2)/2^(3/2))*3 * K = 1.837*K
    K = B / 3
    Winf = W_1 * Linf^W_2
    W0 = 0.1 #  Mass of larva from Andersen et al. 2009
    ## W0=Winf/1e6 we experimented with this parmaeterization but didn't work and took it out (6/23/14) CCM
    L0 = (W0/W_1)^(1/W_2)
    A = B * Winf^(1/3)
    Amat = log((3*K+M)/M) / K   # From Shertzer and Williams 2003
    Lmat = Linf - (Linf-L0)*exp(-K*Amat)
    W_amat = W_1 * Lmat^W_2
    Mat_a = pnorm( (AgeSet-Amat)/(Amat/4)*1.96 )     # 95% of maturity occurs within Amat +/- 25% of Amat
    ## Length at age
    L_a = Linf - (Linf - L0) * exp(-K * AgeSet)
    W_a = W_1 * L_a^W_2
    Lmat = Linf - (Linf - L0) * exp(-K * Amat)
    ## Mortality at age
    if(M_Type=="Length-Based") M_a = (W_a/Winf)^(-1/2) * K
    if(M_Type=="Constant") M_a = rep(M, length=length(AgeSet))
    ## Spawning output at age
    ## R0 can't be in globals since overwritten with NA in the loop below
    R0 = R0_global
    Surv_a = exp(-cumsum(M_a))  # Survival probability through each age
    BPR_a = Surv_a * W_a
    SBPR_a = BPR_a * Mat_a
    B_a = R0 * BPR_a
    SB_a = R0 * SBPR_a
    ## B-H recruitment parameters
    MLRR = exp(LMARR) / (1 - exp(-M))
    h = MLRR / (4 + MLRR)
    SB_0 = sum(SB_a)
    SR_alpha = MLRR / sum(SBPR_a)               # From Myers et al. 1999
    SR_beta = ( (SR_alpha * SB_0) / R0 - 1) / SB_0   # From definition of B-H SRR
    ## Check steepness (Mangel et al. 2013 Eq. 11): sum(SBPR_a)*SR_alpha / (4 + sum(SBPR_a)*SR_alpha)
    ## Check R0 (from B-H definition): SR_alpha * SB_0 / (1 + SR_beta*SB_0)
    ## Check Bmsy/B0 (Mangel et al. 2013 Eq. 16): (sqrt(4*h/(1-h))-1)/((4*h/(1-h))-1)
    ## Fishing params
    Selex_A50 = Selex_A50_times_Amat * Amat
    Selex_a = plogis( Selex_Slope * (AgeSet-Selex_A50))
    ##
    ## The default model behavior is plotted as figure 1, so need to
    ## save that data for the default run (F_I=1) and of course it is
    ## the same across the other dimensions so only save it once.
    if(SensI==1 &  ParI==1 & F_I==1){
        ## Get the biological parameter at unfished levels
        Opt.temp <- optimize( f=ProjectionFn, F=0, interval=Interval,
                             OutputType="Recruitment_Objective", M_a=M_a, W_1= W_1,
                             W_2=W_2, SR_alpha=SR_alpha, SR_beta=SR_beta, Amat=Amat,
                             Lmat=Lmat, Linf=Linf, K=K, h=h, R0=R0, Selex_a=Selex_a,
                             L0=L0, W0=W0)
        temp <-
            ProjectionFn(F=0, Ln_Rinput=Opt.temp$minimum,
                         OutputType="All", M_a, W_1, W_2, SR_alpha,
                         SR_beta, Amat, Lmat, Linf, K, h, R0,
                         Selex_a=Selex_a,L0=L0, W0=W0)
        ## Add in some loose parameters here to the scalars, really ugly
        ## way of making tables
        temp$scalar.params$lambda <- LambdaFn(M_a=M_a, Mat_a= Mat_a, W_a= W_a, SR_alpha=SR_alpha)
        temp$scalar.params$Winf <- Winf
        temp$scalar.params$B <- B
        temp$scalar.params$A <- A
        temp$scalar.params$W_amat <- W_amat
        temp$scalar.params$MLRR <- MLRR
        temp$scalar.params$Linf <- Linf
        temp$scalar.params$LMARR <- LMARR
        temp$scalar.params$M <- M
        temp$scalar.params$Selex_A50_times_Amat <- Selex_A50_times_Amat
        temp$scalar.params$Selex_Slope <- Selex_Slope
        temp$scalar.params$W0 <- W0
        temp$scalar.params$theta1<- W_1
        temp$scalar.params$theta2 <- W_2
        temp$scalar.params$R0 <- R0_global
        Unfished[[SpeciesSet_Example[SpeciesI]]] <- temp
        if(save.lifehist.curves==TRUE){
            SB.temp <- seq(0, SB_0, length=1000)
            message(paste("saving life history data for: ", SpeciesSet[SpeciesI]))
            lifehist.curves[[SpeciesSet[SpeciesI]]] <-
                list(SB=SB.temp, SR_alpha=SR_alpha, SR_beta=SR_beta, AgeSet=AgeSet,
                     Surv_a=Surv_a, BPR_a=BPR_a, SBPR_a=SBPR_a)
        }
    }
    ## Magnitude of change
    Multiplier <- switch(SensSet[SensI], "Big-"=exp(-0.20),
                         "Small-"=exp(-0.1), "Diff-"=exp(-Delta), "Default"=1,
                         "Diff+"=exp(Delta), "Small+"=exp(0.1), "Big+"=exp(0.20))
    ## Changes to inputs
    if(ParSet[ParI]=="M_a_Juv" | ParSet[ParI]=="M_a_Adult" | ParSet[ParI]=="M_a_All"){
        if(ParSet[ParI]=="M_a_All") M_a_Input = c(M_a[1], Multiplier*M_a[-1])
        if(ParSet[ParI]=="M_a_Juv") M_a_Input = c(M_a[1], Multiplier*M_a[c(2:floor(Amat))], M_a[-c(1:floor(Amat))])
        if(ParSet[ParI]=="M_a_Adult") M_a_Input = c(M_a[1:floor(Amat)], M_a[-c(1:floor(Amat))]*Multiplier)
    }else{
        M_a_Input = M_a
    }
    if(ParSet[ParI]=="W_1"){ W_1_Input = W_1 * Multiplier }else{ W_1_Input = W_1 }
    if(ParSet[ParI]=="W_2"){ W_2_Input = W_2 * Multiplier }else{ W_2_Input = W_2 }
    if(ParSet[ParI]=="SR_alpha"){ SR_alpha_Input = SR_alpha * Multiplier }else{ SR_alpha_Input = SR_alpha }
    if(ParSet[ParI]=="SR_beta"){ SR_beta_Input = SR_beta * Multiplier }else{ SR_beta_Input = SR_beta }
    if(ParSet[ParI]=="Amat"){ Amat_Input = Amat * Multiplier }else{ Amat_Input = Amat }
    if(ParSet[ParI]=="Lmat"){
        Lmat_Input = Lmat * Multiplier
        Lmat = Lmat * 1
    }else{
        Lmat = Lmat_Input = NA
    }
    if(ParSet[ParI]=="Linf"){ Linf_Input = Linf * Multiplier }else{ Linf_Input = Linf }
    if(ParSet[ParI]=="K"){ K_Input = K * Multiplier }else{ K_Input = K }
    ## Special case since we have to parameterizations of recruitment
    if(ParSet[ParI]=="h" | ParSet[ParI]=="R0"){
        if(ParSet[ParI]=="h"){
            h_Input = h * Multiplier
            R0_Input = R0
        }
        if(ParSet[ParI]=="R0"){
            h_Input = h
            R0_Input = R0 * Multiplier
        }
    }else{
        h = h_Input = NA
        R0 = R0_Input = NA
    }
    ## Original MSY
    if( (is.na(h_Input) | h_Input<1) ){    # Exclude impossible parameterizations
        if(F_Type=="Optimum" | SensSet[SensI]=="Default")
            Opt_Orig <- optimize( f=EquilibriumFn, interval=c(0,1), maximum=TRUE,
                                 Output="Catch", M_a=M_a, W_1=W_1, W_2=W_2,
                                 SR_alpha=SR_alpha, SR_beta=SR_beta, Amat=Amat,
                                 Lmat=Lmat, Linf=Linf, K=K, h=h, R0=R0,
                                 Selex_a=Selex_a, L0=L0, W0=W0)
        if(F_Type=="Default" & SensSet[SensI]!="Default")
            Opt_Orig = list("maximum"=Results[F_I,SpeciesI,ParI,"Default",'Orig','Fmsy'] )
        Results[F_I,SpeciesI,ParI,SensI,'Orig','Fmsy'] = Opt_Orig$maximum
        Diag_Orig <- EquilibriumFn(F=Opt_Orig$maximum, Output="Diag",
                                   M_a=M_a, W_1=W_1, W_2=W_2,
                                   SR_alpha=SR_alpha, SR_beta=SR_beta,
                                   Amat=Amat, Lmat=Lmat, Linf=Linf, K=K,
                                   h=h, R0=R0,Selex_a=Selex_a, L0=L0, W0=W0 )
        Results[F_I,SpeciesI,ParI,SensI,'Orig','Cmsy'] = Diag_Orig$Catch
        Results[F_I,SpeciesI,ParI,SensI,'Orig','SBmsy'] = Diag_Orig$SB_average
        ## Check if the population crashes
        Diag <- EquilibriumFn(F=0, Output="Diag",
                              M_a=M_a_Input, W_1=W_1_Input,
                              W_2=W_2_Input,
                              SR_alpha=SR_alpha_Input,
                              SR_beta=SR_beta_Input,
                              Amat=Amat_Input,
                              Lmat=Lmat_Input,
                              Linf=Linf_Input, K=K_Input,
                              h=h_Input, R0=R0_Input,Selex_a=Selex_a,
                              L0=L0, W0=W0)
        if(Diag$Recruitment_Objective < 0.01){
            ## New Fmsy
            if(F_Type=="Optimum" | SensSet[SensI]=="Default"){
                Opt_New <- optimize( f=EquilibriumFn,
                                    interval=c(0,Opt_Orig$maximum*(1/4)),
                                    maximum=TRUE, Output="Catch",
                                    M_a=M_a_Input, W_1=W_1_Input,
                                    W_2=W_2_Input, SR_alpha=SR_alpha_Input,
                                    SR_beta=SR_beta_Input, Amat=Amat_Input,
                                    Lmat=Lmat_Input, Linf=Linf_Input,
                                    K=K_Input, h=h_Input,
                                    R0=R0_Input, Selex_a=Selex_a, L0=L0, W0=W0 )
                if(Opt_New$maximum >= Opt_Orig$maximum*(1/4)*0.75)
                    Opt_New <- optimize( f=EquilibriumFn,
                                        interval=c(0,Opt_Orig$maximum*(1)),
                                        maximum=TRUE, Output="Catch", M_a=M_a_Input,
                                        W_1=W_1_Input, W_2=W_2_Input,
                                        SR_alpha=SR_alpha_Input, SR_beta=SR_beta_Input,
                                        Amat=Amat_Input, Lmat=Lmat_Input,
                                        Linf=Linf_Input, K=K_Input, h=h_Input,
                                        R0=R0_Input,Selex_a=Selex_a, L0=L0, W0=W0 )
                if(Opt_New$maximum >= Opt_Orig$maximum*(1)*0.75)
                    Opt_New <- optimize( f=EquilibriumFn,
                                        interval=c(0,Opt_Orig$maximum*(4)),
                                        maximum=TRUE, Output="Catch",
                                        M_a=M_a_Input, W_1=W_1_Input,
                                        W_2=W_2_Input,
                                        SR_alpha=SR_alpha_Input,
                                        SR_beta=SR_beta_Input,
                                        Amat=Amat_Input, Lmat=Lmat_Input,
                                        Linf=Linf_Input, K=K_Input,
                                        h=h_Input,
                                        R0=R0_Input,Selex_a=Selex_a, L0=L0, W0=W0 )
                if(Opt_New$maximum >= Opt_Orig$maximum*(4)*0.75)
                    Opt_New <- optimize( f=EquilibriumFn,
                                        interval=c(0,Opt_Orig$maximum*(16)),
                                        maximum=TRUE,
                                        Output="Catch", M_a=M_a_Input, W_1=W_1_Input,
                                        W_2=W_2_Input,
                                        SR_alpha=SR_alpha_Input,
                                        SR_beta=SR_beta_Input, Amat=Amat_Input,
                                        Lmat=Lmat_Input, Linf=Linf_Input,
                                        K=K_Input, h=h_Input,
                                        R0=R0_Input, Selex_a=Selex_a, L0=L0, W0=W0)
                if(Opt_New$maximum >= Opt_Orig$maximum*(16)*0.75)
                    stop("Check performance of 1-dimensional optimizer")
            }
            if(F_Type=="Default" & SensSet[SensI]!="Default")
                Opt_New <- list(
                    "maximum"=Results[F_I,SpeciesI,ParI,"Default",'New','Fmsy'] )
            Results[F_I,SpeciesI,ParI,SensI,'New','Fmsy'] = Opt_New$maximum
            Diag_New <- EquilibriumFn(F=Opt_New$maximum,
                                      Output="Diag", M_a=M_a_Input, W_1=W_1_Input,
                                      W_2=W_2_Input, SR_alpha=SR_alpha_Input,
                                      SR_beta=SR_beta_Input, Amat=Amat_Input,
                                      Lmat=Lmat_Input, Linf=Linf_Input, K=K_Input,
                                      h=h_Input, R0=R0_Input,Selex_a=Selex_a, L0=L0, W0=W0  )
            Results[F_I,SpeciesI,ParI,SensI,'New','Cmsy'] = Diag_New$Catch
            Results[F_I,SpeciesI,ParI,SensI,'New','SBmsy'] = Diag_New$SB_average
            ## Wrong F
            Results[F_I,SpeciesI,ParI,SensI,'Wrong','Fmsy'] = Opt_Orig$maximum
            Diag_Wrong <- EquilibriumFn(F=Opt_Orig$maximum,
                                        Output="Diag", M_a=M_a_Input, W_1=W_1_Input, W_2=W_2_Input,
                                        SR_alpha=SR_alpha_Input, SR_beta=SR_beta_Input, Amat=Amat_Input,
                                        Lmat=Lmat_Input, Linf=Linf_Input, K=K_Input,
                                        h=h_Input, R0=R0_Input,
                                        Selex_a=Selex_a, L0=L0, W0=W0 )
            Results[F_I,SpeciesI,ParI,SensI,'Wrong','Cmsy'] = Diag_Wrong$Catch
            Results[F_I,SpeciesI,ParI,SensI,'Wrong','SBmsy'] = Diag_Wrong$SB_average
        } else {
            warning(paste("population crashed for: ", paste(collapse = ", ",  c(F_Set[F_I], SpeciesSet[SpeciesI], ParSet[ParI], SensSet[SensI]))))
        }
    }
    ## Test the equilibrium function yield curves, again only for the
    ## example species
    if(save.curve.data==TRUE){
        F.seq <- seq(0, .5, len=50)
        age.temp <- yield.temp <- list()
        ## Loop through range of F values, store the equilibrium yield
        for(i in 1:length(F.seq)){
            test <- EquilibriumFn(F=F.seq[i],
                                  Output="Diag", M_a=M_a_Input, W_1=W_1_Input,
                                  W_2=W_2_Input, SR_alpha=SR_alpha_Input,
                                  SR_beta=SR_beta_Input, Amat=Amat_Input,
                                  Lmat=Lmat_Input, Linf=Linf_Input,
                                  K=K_Input, h=h_Input,
                                  R0=R0_Input, Selex_a=Selex_a, L0=L0, W0=W0)
            SR_curve.temp <- ifelse(is.null(test$SR_curve), NA, test$SR_curve)
            yield.temp[[i]] <- data.frame("Catch"=test$Catch,
                                          "SSB"=test$SB_average,
                                          "SBmsy"=Diag_New$SB_average)
        }
        age.temp <- do.call(rbind, age.temp)
        yield.temp <- do.call(rbind, yield.temp)
        ## Note if h>1 the loop skips, so I need to fill those in
        ## with NA values below, hence the ifelse statement.
        age.list[[k]] <-  as.data.frame(cbind(
            ID=k, F_Type, M=M, M_Type, Multiplier, ParSet[ParI], SensSet[SensI],
            SpeciesSet[SpeciesI], AgeSet, Fmsy=Opt_New$maximum,
            Cmsy=ifelse(is.null(Opt_New$objective), NA, Opt_New$objective),
            data.frame(test[c(2:9,14)]), SR_curve=SR_curve.temp))
        yields.list[[k]] <- data.frame(
            F_Type, M=M, M_Type, Multiplier, ParSet[ParI], SensSet[SensI],
            SpeciesSet[SpeciesI] , F=F.seq, yield.temp, Fmsy=Opt_New$maximum,
            Cmsy=ifelse(is.null(Opt_New$objective), NA, Opt_New$objective))
        k <- k+1
    }
}}}}
    ## End for loops across all scenarios
    ##
    ## Process the results a bit and then return them as a list
    if(save.curve.data){
        age.curves <- as.data.frame(do.call(rbind, age.list))
        names(age.curves)[1:9] <- c("ID", "F_Type", "M","M_Type", "Multiplier","Parameter",
                                    "Delta", "Species", "age")
        age.curves <- droplevels(subset(age.curves, F_Type=="Optimum"))
        yield.curves <- as.data.frame(do.call(rbind, yields.list))
        names(yield.curves) <- c("F_Type", "M","M_Type", "Multiplier","Parameter",
                                 "Delta", "Species", "F", "Catch",  "SSB", "SBmsy",
                                 "Fmsy", "Cmsy")
        yield.curves <- droplevels(subset(yield.curves, F_Type=="Optimum"))
    } else {
        age.curves <- yield.curves <- NULL
    }
    results.list <- list(Results=Results, Unfished=Unfished,
                         age.curves=age.curves,
                         yield.curves=yield.curves,
                         lifehist.curves=lifehist.curves)
    return(results.list)
}

## Old, shouldn't need these anymore. 6/13/2014 CCM
## draw.circle = function(x, y, radius, nv=101, Log=FALSE, ...){
##     for(i in 1:length(radius)){
##         X = ifelse(Log==TRUE,log(x),x) + radius[i] * sin( seq(0,2*pi, length=nv) )
##         Y = ifelse(Log==TRUE,log(y),y) + radius[i] * cos( seq(0,2*pi, length=nv) )
##         if(Log==TRUE) lines( x=exp(X), y=exp(Y), ...)
##         if(Log==FALSE) lines( x=X, y=Y, ...)
##     }
## }
## Optimize = function(f, interval, Eval_Num=51, maximum=FALSE, ...){
##     Eval_Loc = seq(interval[1], interval[2], length=Eval_Num)
##     Eval_Obj = rep(NA, Eval_Num)
##     for(i in 1:length(Eval_Obj)) Eval_Obj[i] = f( Eval_Loc, ... )
## }
