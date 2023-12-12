## Install packages (from CRAN)
install.packages("qtlpoly")
install.packages("mappoly")

## ## (Optional) Install QTLpoly's development version
## install.packages(devtools)
## devtools::install_github("gabrielgesteira/QTLpoly")

## Load packages
library(qtlpoly)
library(mappoly)

## Calculate genotype probabilities
genoprob4x = lapply(maps4x, calc_genoprob)

## Read phenotypic data
## Make sure that:
## - individual names match the ones in the map
## - individual names read as row names
## - one column for each phenotype
## Missing data is allowed - QTLpoly will handle that
phenotypes = read.table("pheno.txt", header=T, sep=" ", row.names = 1)

## Check individual names
dimnames(genoprob4x[[1]]$probs)[[3]]
rownames(phenotypes)
all(rownames(phenotypes) %in% dimnames(genoprob4x[[1]]$probs)[[3]])

## Read data to a QTLpoly object
dat = read_data(ploidy = 4, geno.prob = genoprob4x, pheno = phenotypes, step = 1)
print(dat, detailed = TRUE)

## Score-based resampling method
data.sim = simulate_qtl(data = dat, mu = 0, h2.qtl = NULL, var.error = 1, n.sim = 10, missing = TRUE, seed = 123)
score.null = null_model(data = data.sim$results, n.clusters = 1, plot = NULL)
## score.null = null_model2(data = data.sim$results, n.clusters = 12, plot = NULL)
min.pvl = unlist(lapply(score.null$results, function(x) return(x$pval[which.max(x$stat)])))
quantile(sort(min.pvl), c(0.2,0.05))
## Setting thresholds
sig.fwd = quantile(sort(min.pvl), 0.2)
sig.bwd = quantile(sort(min.pvl), 0.05)
## sig.fwd = 0.0011493379; sig.bwd = 0.0002284465  # 20% and 5% genome-wide significance level (1000 samples)

## Manual QTL detection
null.mod = null_model(data = dat, pheno.col = 1, n.clusters = 1)
## null.mod = null_model2(data = dat, pheno.col = 1, n.clusters = 12)
print(null.mod)
search.mod = search_qtl(data = dat, model = null.mod, w.size = 15, sig.fwd = 0.01, n.clusters = 12)
print(search.mod)
optimize.mod = optimize_qtl(data = dat, model = search.mod, sig.bwd = 1e-04, n.clusters = 12)
print(optimize.mod)
search.mod2 = search_qtl(data = dat, model = optimize.mod, sig.fwd = 1e-04, n.clusters = 12)
print(search.mod2)
profile.mod = profile_qtl(data = dat, model = optimize.mod, d.sint = 1.5, polygenes = FALSE, n.clusters = 12)
print(profile.mod)
print(profile.mod, sint = "lower")
print(profile.mod, sint = "upper")

## Automatic QTL detection
remim.mod = remim(data = dat, pheno.col = 1, w.size = 15, sig.fwd = sig.fwd, sig.bwd = sig.bwd, d.sint = 1.5, n.clusters = 1)
remim.mod = remim2(data = dat, pheno.col = 1, w.size = 15, sig.fwd = sig.fwd, sig.bwd = sig.bwd, d.sint = 1.5, n.clusters = 12)
print(remim.mod)
print(remim.mod, sint = "lower")
print(remim.mod, sint = "upper")
plot_profile(data = dat, model = remim.mod, grid = TRUE)
plot_sint(data = dat, model = remim.mod)

## Fit final QTL model
fitted.mod = fit_model(data = dat, model = remim.mod)
summary(fitted.mod)
plot_qtl(data = dat, model = remim.mod, fitted = fitted.mod)

## Estimate allele effects
est.effects = qtl_effects(ploidy = 4, fitted = fitted.mod)
plot(est.effects, p1 = "Atlantic", p2 = "B1829-5")

## Predict QTL-based breeding values
y.hat = breeding_values(data = dat, fitted = fitted.mod)
plot(y.hat)
rownames(y.hat$results$FM07$y.hat)[which.max(y.hat$results$FM07$y.hat)] # Best individual

## Run FEIM for comparison
perm = permutations(data = dat, n.sim = 10, n.clusters = 1)
print(perm)
(sig.lod = perm$sig.lod$`0.95`)
## sig.lod = c(5.68, 5.78, 5.6)  # 5% genome-wide significance level
feim.mod = feim(data = dat, w.size = 15, sig.lod = sig.lod)
print(feim.mod)
plot_profile(data = dat, model = feim.mod, grid = TRUE)

## Exporting results to VIEWpoly
save(maps4x, file = "mappoly.maps.RData")
save(dat, file = "qtlpoly.data.RData")
save(remim.mod, file = "qtlpoly.remim.mod.RData")
save(fitted.mod, file = "qtlpoly.fitted.mod.RData")
save(est.effects, file = "qtlpoly.est.effects.RData")
