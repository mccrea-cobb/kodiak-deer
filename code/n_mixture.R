cat(file = "./code/jags/nmix.jags", "
    model {
    
    # Priors
    for (i in 1:nclasses) {
      p[i] ~ dunif(0,1)  # Detection
    }
    

    for(i in 1:nclasses){                   # For each age/sex class,
        N[i,1] ~ dcat(psi[,i])                 # Abundance at time 1 is a categorical dist, drawn from probs in psi
    } #i
     
    # Ecological model
    for (i in 1:nclasses) {  # For each age/sex class
      for(t in 2:nyears) {   # For each year
        N[i,t] ~ dpois(N[i,t-1])
      }
    }
    
    # Observation model
    for (i in 1:nclasses) {             # For each age/sex class
      for(t in 1:nyears) {                # For each year
        for (j in 1:nsurveys) {             # For each survey rep
          C[i,j,t] ~ dbin(p[i], N[i,t])
        } #j
      } #t
    } #i
    
    
    # Derived parameters
  
    # Total elk abundance
    for (t in (1:nyears)){
      Ntot[t] <- round(sum(N[1:5,t]))
    }
    }
    ")

# Load data
source("./code/functions/import_format.R")
dat <- import_format("Tomales", 1978, 2020)

# Load data
source("./code/functions/import_format_array.R")
C <- import_format_array("Tomales", 1978, 2020)

# Define parameters to monitor
parameters <- c("Ntot", "N", "p")

# MCMC settings
nc <- 3
na <- NULL
ni <- 100000
nb <- 10000
nt <- 10

nclasses <-  5 # age classes
nyears <- ncol(C[1,,])  # years
nsurveys <- 3  # survey replicates

# Calculate psi (for dcat distribution), based on the max value (count) in each row (plot) for each dimension (year)
maxCyear1 <- apply(C, c(1,3), max, na.rm = T)[,1]   # Max counts/plot in year 1
maxCyear1 <- replace(maxCyear1, maxCyear1==0, 1)  # Replace zeros with 1
psi_range <- max(maxCyear1)*2  # Define the range of values for psi

psi <- matrix(0, psi_range, nclasses)   # An empty matrix to hold psi
for(i in 1:nclasses){                   # Fill that matrix with probabilities summing to one, zeros for anything less than the maxCyear1 value
  psi[1:(maxCyear1[i]-1), i] <- 0
  psi[maxCyear1[i]:psi_range, i] <- 1/((psi_range+1) - maxCyear1[i])
}

save_dir <- "./products/models/rdata/"

start <- Sys.time()

## Initial values
Nst <- apply(C, c(1,3), max, na.rm=TRUE) + 1
Nst[Nst == "-Inf"] <- 1
inits <- function() list(N = Nst,
                         p = rep(0.5, nclasses)
                         )

bugs_data <-  
  list(
    C = C,
    nyears = ncol(C[1,,]),
    nclasses = nclasses,
    nsurveys = nsurveys,    # survey replicates
    psi = psi
  )

ipm <-
  jags(data = bugs_data, 
       inits = inits, 
       parameters.to.save = parameters, 
       model.file = "code/jags/Nmix.jags", 
       n.chains = nc,
       n.adapt = na, 
       n.iter = ni, 
       n.burnin = nb, 
       n.thin = nt, 
       parallel = T
  )
