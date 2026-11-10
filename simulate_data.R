
data_dir <- '~/indepth_hackathon/data'
list.files(data_dir)

bmi <- haven::read_sav(file.path(data_dir, 'GenR_sdsfbmi_for_hackathon.sav'))
  # dplyr::mutate(IDC = as.numeric(IDC))

summary(bmi[, grepl('bmi', names(bmi))])

adhd <- read.csv(file.path(data_dir, 'ADHD_df.csv')) # ignore z scores 
summary(adhd)

cell <- read.csv(file.path(data_dir, 'birth_450k_unilife_observed_release4_zscored.csv'))
summary(cell)

data <- merge(bmi, adhd, by = 'IDC')

N = nrow(adhd)

set.seed(73)
adhd$pheno_null_norm <- rnorm(n=N, mean = 0, sd = 1)
adhd$pheno_null_count <- rnbinom(n=N, size = 1, mu=3) # rpois lambda=1 not dispersed enough
adhd$pheno_null_binar <- rbinom(n=N, size = 1, prob = 0.15)

summary(adhd)
adhd$sum_att_9m_z <- NULL
adhd$Attention_Deficit__Hyperactivity_Problems_TScore_9m <- NULL

saveRDS(adhd, 'sim_fullcohort.rds')

# ==============================================================================

selection_file <- haven::read_sav('~/methylCHART/data/DNAm_selection_file.sav') |>
  dplyr::mutate(
    Period = factor(Period, levels = c("Birth", "Age5", "Age9", "Age13", "Age17")),
    Array = as.factor(Batch),
    Batch = interaction(Period, Array, drop = TRUE))

table(selection_file$Batch)

birth450k <- selection_file |> dplyr::filter(Batch == 'Birth.450k')
age9y450k <- selection_file |> dplyr::filter(Batch == 'Age9.450k')

table(birth450k$Bridge)

# waldo::compare(birth450k$Sample_ID, cell$Sample_ID)
# setdiff(cell$Sample_ID, birth450k$Sample_ID)
samples_dropped <- setdiff(birth450k$Sample_ID, cell$Sample_ID)
children_dropped <- setdiff(birth450k$IDC, cell$IDC)

length(samples_dropped)

bridges <- birth450k$Sample_ID[birth450k$Bridge=='Bridge_Sample']

intersect(samples_dropped, bridges) # 56 samples dropped and none are not bridges?? why?

# Simulate some data ===========================================================

set.seed(73)

genr_dnam <- bigmemory::attach.big.matrix(
  file.path('~/methylCHART/data/genr_dnam.desc'), readonly = TRUE, lockfile = TRUE)


make_association <- function(subset, name='450K_Birth') {
  
  dnam_sub <- genr_dnam[, 1:10]
  dim(dnam_sub)
  dnam_sub <- dnam_sub[subset$Sample_ID, ]
  dim(dnam_sub)
  
  N = nrow(subset)
  
  betas <- c(seq(-0.9, -0.1, by=0.2), seq(0.1, 0.9, by=0.2))
  # eta <- as.vector(dnam_sub %*% betas)
  
  S <- cov(dnam_sub)
  beta_joint <- solve(S, diag(S) * betas)  # joint coefs that give marginal slopes = betas
  etaj <- as.vector(dnam_sub %*% beta_joint)
  
  subset$pheno_10hits_norm <- etaj + rnorm(n=N, mean=0, sd=0.3)
  subset$pheno_10hits_count <- rnbinom(n=N, size = 2, mu = exp(1+0.3*etaj))
  subset$pheno_10hits_binar <- rbinom(n=N, size = 1, prob = plogis(-0.5 + etaj))
  
  names(betas) <- names(beta_joint)
  sink(paste0('sim_beta_10hits',name,'.txt'))
  print(betas)
  print(beta_joint)
  sink()
  
  saveRDS(subset, file=paste0('sim_data_',name,'.rds'))
  
  return(subset)
}

d1 <- make_association(birth450k, '450K_Birth')
d2 <- make_association(age9y450k, '450K_Age9y')

# check marginals 
# dnam_sub <- genr_dnam[, 1:10]
# dim(dnam_sub)
# dnam_sub <- dnam_sub[d1$Sample_ID, ]
# 
# d2 <- cbind(d1, dnam_sub)
# 
# for (cg in colnames(dnam_sub)) {
#   fit = lm(as.formula(paste0('pheno_10hits_norm ~ ', cg)), data=d2)
#   
#   message(cg, ' : ', round(fit$coef[2], 2))
# }

N = nrow(birth450k)

birth450k$pheno_null_norm <- rnorm(n=N, mean = 0, sd = 1)
birth450k$pheno_null_count <- rnbinom(n=N, size = 1, mu=3) # rpois lambda=1 not dispersed enough
birth450k$pheno_null_binar <- rbinom(n=N, size = 1, prob = 0.15)

hist(birth450k$pheno_null_norm)
hist(birth450k$pheno_null_count)
table(birth450k$pheno_null_binar)

# Now some associated phenos 


# R2 <- 0.5
# sigma <- sqrt(var(eta) * (1 - R2) / R2)
# y_norm <- eta + rnorm(N, 0, sigma)

# Counts: drop the 0.3 shrinkage (or rescale eta), less overdispersion
# eta_s <- eta / sd(eta)
# y_count <- rnbinom(N, size = 5, mu = exp(1 + 0.5 * eta_s))

# Binary
# y_binar <- rbinom(N, 1, plogis(-0.5 + 1 * eta_s))

# cor(y_norm, X, method = "spearman")

