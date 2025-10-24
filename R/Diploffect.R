
Diploffect_Zupdater <- setRefClass(
  "Diploffect_Zupdater",
  fields=list(
    hap.names  = "character", # haplotype names
    num.haps   = "integer", # number of haplotypes
    num.dips   = "integer", 
    dip.names  = "character", # diplotype names in prob.format order
    prob.D     = "matrix", # prior diplotype probabilities
    num.obs    = "integer",
    prob.format = "character", # upper.tri or lower.tri
    components = "list", # list(<name>="additive", <name>="dominance", or <name>="diplotype")
    Amat = "matrix" # multiplier matrix converting diplotypes to haplotype counts
    )
)

Diploffect_Zupdater$methods(
  initialize = function(
                       hap.names,
                       prob.D,
                       prob.format,
                       components)
  {
    .self$hap.names <- hap.names
    .self$num.haps <- length(hap.names)
    if (ncol(prob.D)!=num.haps*(num.haps+1)/2){
      stop(num.haps, "haplotypes implies", num.haps*(num.haps-1)/2, "diplotype probabilities,",
           "not", ncol(prob.D))
    }
    .self$prob.D <- prob.D
    if ("upper.tri"!=prob.format & "lower.tri"!=prob.format){
      stop("prob.format must be upper.tri or lower.tri, not '", prob.format, "'")
    }
    .self$num.dips <- ncol(prob.D)
    .self$num.obs  <- nrow(prob.D)
    .self$prob.format <- prob.format
    .self$Amat <- make.dip2hap.matrix(num.haps, prob.format=prob.format, hap.names=hap.names)
    .self$dip.names <- rownames(Amat)
    for (r in names(components)){
      if (!(components[[r]] %in% c("additive", "dominance", "diplotype"))) {
        stop("Unrecognized type ", components[[r]], " for component '", r, "'")
      }
    }
    .self$components <- components
  }
)

Diploffect_Zupdater$methods(
  whichComponents = function(){
    names(.self$components)
  }
)

Diploffect_Zupdater$methods(
  sample.Dip = function(probmat){
    Dip <- matrix(0, nrow=num.obs, ncol=num.dips)
    colnames(Dip) <- dip.names
    rownames(Dip) <- rownames(probmat)
    jj <- 1:num.dips
    for (i in 1:num.obs){
      d <- sample(jj, size=1, prob=probmat[i, ])
      Dip[i, d] <- 1
    }
    Dip
  }
)

Diploffect_Zupdater$methods(
  sample.posterior.Dip = function(yres, Vinv, state){
    # make predicted y for each possible diplotype state
    Dip.reps <- diag(rep(1, num.dips)) # 1 representative of each diplotype state
    Z.reps   <- .self$Dip.to.Z(Dip.reps) # corresponding Z's for each dip
    ypred.reps <- numeric(num.dips) # predicted value based on ranef and Z's for each dip
    for (r in names(components)){
      rpred <- Z.reps[[r]] %*% state$random[[r]]$ranef # contribution from r
      ypred.reps <- ypred.reps + rpred
    }
    # make p.D.given.y for each individual
    prob.D.given.y <- matrix(0, nrow=num.obs, ncol=num.dips)
    for (i in 1:num.obs){
      prob.D.given.y[i, ] <- Vinv[i, i]*(yres[i] - ypred.reps)^2
    }
    prob.D.given.y <- exp( - 0.5 * state$sigma2inv * prob.D.given.y )
    prob.D.given.y <- prob.D.given.y * prob.D
    prob.D.given.y <- prob.D.given.y / rowSums(prob.D.given.y)
    if (FALSE){ # diagnostics
      im <- which.max((rowSums((prob.D-prob.D.given.y)^2)))
      plot(prob.D[im,], prob.D.given.y[im,])
      # add KL divergence, etc?
    }
    # sample a new Dip
    .self$sample.Dip(prob.D.given.y)
  }
)


Diploffect_Zupdater$methods(
  sample.posterior.Z = function(yres, Vinv, state){
    Dip <- sample.posterior.Dip(yres, Vinv, state)
    # TODO: option to save Dip 
    .self$Dip.to.Z(Dip)
  }
)                                   
    
Diploffect_Zupdater$methods(
  sample.prior.Dip = function(){
    .self$sample.Dip(.self$prob.D)
  }
)

Diploffect_Zupdater$methods(
  Dip.to.Z = function(Dip){
    Z.list <- list()
    for (r in names(.self$components)){
      type <- .self$components[[r]]
      if ("additive"==type){
        Z.list[[r]] <- Dip %*% .self$Amat
      } else if ("dominance"==type) {
        Z.list[[r]] <- Dip[, -(1:num.haps)] # indicator for hets
      } else {
        Z.list[[r]] <- Dip # diplotype effect
      }
    }
    Z.list
  }
)

Diploffect_Zupdater$methods(
  sample.prior.Z = function(){
    Dip <- .self$sample.prior.Dip()
    .self$Dip.to.Z(Dip)
  }
)
        




