library(WVmisc)
library(mnormt)
library(coda)
library(Matrix)

source.Rfiles <- function(dir){
  files <- list.files(path=dir, pattern="*.R$", full.names=TRUE)
  for (f in files) source(f)
}

# manually load R scripts
code.dir <- "/Users/wvaldar/Dropbox/wvaldar/code/R/packages/"
gibbs.code.dir <- file.path(code.dir, "/WVGibbsLmm/github/WVGibbsLmm/")
diplo.code.dir <- file.path(code.dir, "Diploffect.Gibbs/github/Diploffect.Gibbs/")
source.Rfiles(file.path(gibbs.code.dir, "R"))
source.Rfiles(file.path(diplo.code.dir, "R"))

## RBKS data
dat <- readRDS("data/diplo_tests/RBKS_diplo_example.RDS")
hap.names <- LETTERS[1:8] # can use any vector of strings

dzu <- Diploffect_Zupdater$new(hap.names=hap.names,
                               prob.D=dat$happy_locusmatrix,
                               prob.format="upper.tri",
                               components=list(hap="additive", dom="dominance"))
Z.list.init <- dzu$sample.prior.Z()
#Z.list.init <- dzu$expect.prior.Z()

X <- model.matrix(y ~ expr.pc1 + SexM + DietLOW.FAT, data=dat$data)

# instantiate ... allow all effects to be sampled in a block (this is the default anyway)
g <- lmmgibbs(y=dat$data$y, X=X,
              prior.beta = normal.prior(mean=rep(0, ncol(X)), var=10e3),
              prior.sigma2inv = gamma.prior(shape=1, rate=1),
              prior.tau2inv = gamma.prior(shape=1, rate=1),
              random=list(
                kinship = normal.ranef(Z=diag(rep(1,nrow(X))),
                                       corrmat=dat$K,
                                       stz="project"),
                hap = normal.ranef(Z=Z.list.init$hap, stz="project"), # "project" or "basis"
                dom = normal.ranef(Z=Z.list.init$dom, stz="project") # "project" or "basis"
              ),
              Z.updater=dzu,
              Z.update.every=1,
              Z.update.from="posterior", # "prior" for MI or "posterior" for full Gibbs sampling
              )

date()
out <- g$sample(wanted=c("fixed", "random.hap", "random.dom", "sigma2inv", "tau2inv"),
                n.iter=1e4, verbose.at=100, thin=10)
date()

# example posterior contrast: probability of hap.G < hap.H is
# mean((out[, "random.hap.ranef.G"] - out[, "random.hap.ranef.H"])<0)


# to clean up plot labels
rm.label.junk <- function(s){
  s <- sub("random.", "", s)
  s <- sub("ranef.", "", s)
  s <- sub("fixed.", "", s)
  s
}
  
pdf(file="hpd.pdf", height=10, width=6)
# fixed effects
wanted <- grep("fixed.", colnames(out), value=TRUE)
plot.hpd(out, wanted=wanted, names=rm.label.junk(wanted))
abline(v=0, col=gray(0.8))
# haplotypes
wanted <- setdiff(grep("hap.", colnames(out), value=TRUE), "random.hap.tau2inv")
plot.hpd(out, wanted=wanted, names=rm.label.junk(wanted))
abline(v=0, col=gray(0.8))
# dominance deviations
wanted <- setdiff(grep("dom.", colnames(out), value=TRUE), grep("inv", colnames(out), value=TRUE))
plot.hpd(out, wanted=wanted, names=rm.label.junk(wanted))
abline(v=0, col=gray(0.8))
dev.off()





## CHOL data
# load data
load(file.path(diplo.code.dir, "/data/diplo_tests/CHOL_thu_ex.RData"))
kk <- match(data$SUBJECT.NAME, colnames(K))
K <- as.matrix(K)
K <- K[kk, kk]
n <- nrow(data)
all(data$SUBJECT.NAME==colnames(K)) # sanity check
all(data$SUBJECT.NAME==rownames(locusmatrix)) # sanity check
hap.names <- LETTERS[1:8] # can use any vector of strings

# Make WVGibbsLmm add-on for updating diplotype Z matrices
# TODO: test for upper tri vs lower tri
dzu <- Diploffect_Zupdater$new(hap.names=hap.names, prob.D=locusmatrix,
                              prob.format="upper.tri",
                              components=list(hap="additive", dom="dominance"))
Z.list.init <- dzu$sample.prior.Z()

X <- model.matrix(rint.CHOL ~ Sex + Diet, data=data)

# instantiate ... allow all effects to be sampled in a block (this is the default anyway)
g <- lmmgibbs(y=data$rint.CHOL, X=X,
              prior.beta = normal.prior(mean=rep(0, ncol(X)), var=10e3),
              prior.sigma2inv = gamma.prior(shape=1, rate=1),
              prior.tau2inv = gamma.prior(shape=1, rate=1),
              random=list(
                kinship = normal.ranef(Z=diag(rep(1,n)), corrmat=K, stz="project"),
                hap = normal.ranef(Z=Z.list.init$hap, stz="basis"),
                dom = normal.ranef(Z=Z.list.init$dom, stz="basis")
              ),
              Z.updater=dzu,
              Z.update.every=5,
              )

date()
out <- g$sample(wanted=c("fixed", "random.hap", "random.dom", "sigma2inv", "tau2inv"),
                n.iter=10, verbose.at=1, thin=1)
date()

# to clean up plot labels
rm.label.junk <- function(s){
  s <- sub("random.", "", s)
  s <- sub("ranef.", "", s)
  s <- sub("fixed.", "", s)
  s
}
  
pdf(file="hpd.pdf", height=10, width=6)
# fixed effects
wanted <- grep("fixed.", colnames(out), value=TRUE)
plot.hpd(out, wanted=wanted, names=rm.label.junk(wanted))
abline(v=0, col=gray(0.8))
# haplotypes
wanted <- setdiff(grep("hap.", colnames(out), value=TRUE), "random.hap.tau2inv")
plot.hpd(out, wanted=wanted, names=rm.label.junk(wanted))
abline(v=0, col=gray(0.8))
# dominance deviations
wanted <- setdiff(grep("dom.", colnames(out), value=TRUE), grep("inv", colnames(out), value=TRUE))
plot.hpd(out, wanted=wanted, names=rm.label.junk(wanted))
abline(v=0, col=gray(0.8))
dev.off()





