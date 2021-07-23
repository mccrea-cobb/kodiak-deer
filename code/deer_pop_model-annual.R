
sink("./code/jags/deer_pop_model.jags")
cat("
    model { 
    
    #### Priors and constraints #####
    
    # Observation detection rate 
    for (i in 1:nclasses) {
      p[i] ~ dunif(0,1)
    }

    # Initial abundances
    for(i in 1:nclasses){       # For each age/sex class,
        N[i,1] ~ dcat(psi[,i])    # Abundance at time 1 is a categorical dist, drawn from probs in psi
    } #i
    
    # Initial lagged density at time 1
    D.l[1] <- D[1]
    
    
    # Relationships for vital rate parameters
    for (t in 1:(nyears-1)){      
      logit(saf[t]) <- logit.b0.saf + eps[t,1]    # Adult female survival
      logit(sam[t]) <- logit.b0.sam + eps[t,2]    # Adult male survival
      logit(sf[t]) <- logit.b0.sf + eps[t,3]      # Fawn survival
      log(f[t]) <- log.bo.f + eps[t,4]            # Fecundity
    }
    
    
    # Priors for precision of vital rates 
    for (i in 1:4) {zero[i] <- 0}
    for (t in 1:(ti-1)){
      eps[t,1:4] ~ dmnorm(zero[], Omega[,])
    }
    
    # Priors for precision matrix
    Omega[1:4, 1:4] ~ dwish(R[ , ], 5)
    Sigma[1:4, 1:4] <- inverse(Omega[ , ])
    

    # Priors for mean demographic rates
    b0.saf ~ dunif(0, 1)    # Adult female survival
    b0.sam ~ dunif(0, 1)    # Adult male survival
    b0.sf ~ dunif(0, 1)     # Fawn survival
    b0.f ~ dunif(0, 2)      # Fecundity
    
    # Back-transformation of mean demographic rates
    logit.b0.saf <- log(b0.saf / (1 - b0.saf))    # Adult female survival
    logit.b0.sam <- log(b0.sam / (1 - b0.sam))    # Adult male survival
    logit.b0.sf <- log(b0.sf / (1 - b0.sf))       # Fawn survival
    log.b0.f <- log(b0.f)                         # Fecundity


    ##### Likelihoods #####
    
    # System process
    for (t in 2:nyears){  
      N[1,t] ~ dpois((Nbf[t]) * f[t])                   # Newborn fawn abundance (spring) at time t (surviving breeding females*fecundity rate)
      
      N[2,t] ~ dbin(sf[t-1], N[1,t-1])                  # Yearling abundance at time t (fawns at t-1 that survived)
      Nfy[t] <- dbin(0.5, N[2,t])                       # Yearling female abundance at time t (to estimate fawn abundance and adult females)
      Nmy[t] <- 1 - Nfy[t]                              # Yearling male abundance at time t (to estimate adult males)
      
      N[3,t] ~ dbin(saf[t-1], (Nfy[t-1] + N[3,t-1]))    # Adult female abundance at time t (adult female survival at time t-1 * cow and female yearling abundance at time t, adjusted for removals)
      N[4,t] ~ dbin(sam[t-1], (Nmy[t-1] + N[4,t-1]))    # Adult male abundance at time t (adult male survival at time t-1 * bull and male yearling abundance at time t, adjusted for removals)
    }
    
    
    # Observation model
    for (i in 1:nclasses) {             # For each age/sex class
      for(t in 1:nyears) {                # For each year
        for (j in 1:nsurveys) {             # For each survey rep
          C[i,j,t] ~ dbin(p[i], N[i,t])
        } #j
      } #t
    } #i
    
    for (t in (1:nyears)){
      Nbf[t] <- Nfy[t] + N[3,t]  # Number of breeding females (yearlings and adults)
    }
    

    ##### Derived parameters ######
    
    # Total abundance
    for (t in (1:nyears)){
      Ntot[t] <- round(sum(N[1:nclasses,t]))
    }
    
    # Density
    for (t in (1:nyears)){
      D[t] <- Ntot[t]/area[t]
    }

    # Density, lagged 1 year
    for (t in (2:nyears)){
      D.l[t] <- D[t-1]
    }
    
    # Annual herd growth rate (lambda)
     for (t in 1:(nyears-1)){
     lambda[t] <- Ntot[t+1] / (Ntot[t] + 1E-10)   # Add a small values to avoid zeros in the denominator
     }

    # Average adult female survival
    saf.avg <- mean(saf)

    # Average adult male survival
    sam.avg <- mean(sam)

    # Average yearling survival
    sy.avg <- mean(sy)

    # Average fecundity
    f.avg <- mean(f)

    # Average lambda
    lambda.avg <- mean(lambda)
    
    # # Goodness-of-fit statistics
    # for(t in 1:nyears){
    #   gf.f[t] <- pow((pow(C[1,t],0.5) - pow(N[1,t],0.5)),2)   # Fawns
    #   gf.y[t] <- pow((pow(C[2,t],0.5) - pow(N[2,t],0.5)),2)   # Yearlings
    #   gf.af[t] <- pow((pow(C[4,t],0.5) - pow(N[4,t],0.5)),2)  # Adult females
    #   gf.am[t] <- pow((pow(C[5,t],0.5) - pow(N[5,t],0.5)),2)  # Adult males
    # }
    
    }
  ", fill = TRUE)
sink()


require(jagsUI)

model <- "./code/jags/deer_pop_model.jags"

startyear <- 1
endyear <- 10

# Load data
source("./code/functions/import_format_array.R")
C <- import_format_array(filedir = "./data/raw_data/survey/survey_data_reps.csv",
                         herd = "Tomales", 
                         startyear = startyear,
                         endyear = endyear)

# Define parameters to monitor
parameters <- c("Ntot", "N", "D", "lambda",
                "f", "sy", "saf", "sam", 
                "f.avg", "sy.avg", "saf.avg", "sam.avg", "lambda.avg",
                "b0.f", "b0.sy", "b0.saf", "b0.sam", "p"
                )

# MCMC settings
nc <- 3
na <- NULL
ni <- 100000
nb <- 10000
nt <- 10

nclasses <-  nrow(C) # age classes
nyears <- ncol(C[1,,])  # years
nsurveys <- 3  # survey replicates

# Calculate psi (for dcat distribution), based on the max value (count) in each row (plot) for each dimension (year)
maxCyear1 <- apply(C, c(1,3), max, na.rm = T)[,1]   # Max counts/plot in year 1
maxCyear1 <- replace(maxCyear1, maxCyear1==0, 1)  # Replace zeros with 1
psi_range <- max(maxCyear1)*3  # Define the range of values for psi

psi <- matrix(0, psi_range, nclasses)   # An empty matrix to hold psi
for(i in 1:nclasses){                   # Fill that matrix with probabilities summing to one, zeros for anything less than the maxCyear1 value
  psi[1:(maxCyear1[i]-1), i] <- 0
  psi[maxCyear1[i]:psi_range, i] <- 1/((psi_range+1) - maxCyear1[i])
}

save_dir <- "./products/models/rdata/tomales/"

ipm = NULL
start <- Sys.time()

## Initial values
Nst <- apply(C, c(1,3), max, na.rm=TRUE) + 1
Nst[Nst == "-Inf"] <- 1
inits <- function(){
  list(b0.saf = runif(1,0,0.5),
       b0.sam = runif(1,0,0.5),
       b0.sy = runif(1,0,0.5),
       b0.p = runif(1,0,0.5),
       b = rnorm(3, 0.2, 0.5),
       Omega = diag(4),
       N = Nst,
       po = rep(0.8, nclasses)
       # po = 0.5,
       # po = array(.8, c(nclasses, nyears))
       )
}

bugs_data <-  
  list(
    C = C,
    R.c = dat[, "removals_cows"],
    R.b = dat[, "removals_bulls"],
    contra = dat[, "contra"],
    ti = ncol(C[1,,]),
    area = dat[, "herd_range_size"],
    R = diag(c(1,1,1,1), ncol = 4),  
    nyears = ncol(C[1,,]),
    nclasses = nclasses,
    nsurveys = nsurveys,    # survey replicates
    psi = psi
  )

ipm <-
  jags(data = bugs_data, 
       inits = inits, 
       parameters.to.save = parameters, 
       model.file = "code/jags/tomales/null.jags", 
       n.chains = nc,
       n.adapt = na, 
       n.iter = ni, 
       n.burnin = nb, 
       n.thin = nt, 
       parallel = T
       )

save(ipm, file = paste0(save_dir, paste("null"), ".RData"))

end = Sys.time()
end - start
