# SCLC survival modelling and external validation
# Configure file paths with environment variables before running.

required_packages <- c(
  "MASS",
  "reshape",
  "data.table",
  "survival",
  "survminer",
  "survIDINRI",
  "survivalsvm",
  "haven",
  "readxl",
  "tableone",
  "lubridate",
  "grid",
  "gridExtra",
  "RColorBrewer",
  "ggExtra",
  "ggsignif",
  "ggpubr",
  "ggprism",
  "ggDCA",
  "rms",
  "nomogramFormula",
  "timeROC",
  "recipes",
  "broom",
  "stringr",
  "tibble",
  "tidyr",
  "dplyr",
  "ggplot2",
  "tidyverse",
  "mlr3verse",
  "mlr3proba",
  "mlr3extralearners"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

development_data_path <- Sys.getenv("MOLAN_DEVELOPMENT_DATA", unset = file.path("data", "development_cohort.sav"))
external_data_path <- Sys.getenv("MOLAN_EXTERNAL_DATA", unset = file.path("data", "external_validation_cohort.sav"))
output_dir <- Sys.getenv("MOLAN_OUTPUT_DIR", unset = "results")
model_dir <- file.path(output_dir, "models")
figure_dir <- file.path(output_dir, "figures")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(development_data_path)) stop("Development data not found: ", development_data_path)
if (!file.exists(external_data_path)) stop("External validation data not found: ", external_data_path)



development_data <- haven::read_sav(development_data_path)
external_data <- haven::read_sav(external_data_path)
set.seed(347)
split <- rsample::initial_split(development_data, prop = 0.7, strata = "status_os")
train_data <- rsample::training(split)
validation_data   <- rsample::testing(split)
Charcutpoint <- surv_cutpoint(train_data, time = "OS", event = "status_os",
                              variables = "ITHScore")
plot(Charcutpoint,labels=c("High Risk","Low Risk"),palette = c("#5396E6","#DB0138"))

median_score <- as.numeric(Charcutpoint["ITHScore"][[1]]["estimate"][[1]][[1]])
train_data$ITHScoreG <- ifelse(train_data$ITHScore > median_score, "High", "Low")
external_data$ITHScoreG <- ifelse(external_data$ITHScore > median_score, "High", "Low")
validation_data$ITHScoreG <- ifelse(validation_data$ITHScore > median_score, "High", "Low")

p_train <- coxph(Surv(OS, status_os) ~ ITHScoreG, data = train_data) %>% summary() %>% coefficients %>% .[,"Pr(>|z|)"]
p_external_cohort <- coxph(Surv(OS, status_os) ~ ITHScoreG, data = external_data) %>% summary() %>% coefficients %>% .[,"Pr(>|z|)"]
survfit(Surv(OS,status_os)~ITHScoreG, data = train_data) %>%
  ggsurvplot(
    conf.int = TRUE,
    xlab = "Time (months)",
    ylab = "OS",
    title = "Train Cohort-ITHScoreG",
    pval = TRUE,
    risk.table = TRUE,
    xlim = c(0,60),
    break.time.by = 12,
    censor.size = 7,
    palette = c("#3c5769","#b40000")
  )

survfit(Surv(OS,status_os)~ITHScoreG, data = external_data) %>%
  ggsurvplot(
    conf.int = TRUE,
    xlab = "Time (months)",
    ylab = "OS",
    title = "External Cohort-ITHScoreG",
    pval = TRUE,
    risk.table = TRUE,
    xlim = c(0,60),
    break.time.by = 12,
    censor.size = 7,
    palette = c("#3c5769","#b40000")
  )
survfit(Surv(OS,status_os)~ITHScoreG, data = validation_data) %>%
  ggsurvplot(
    conf.int = TRUE,
    xlab = "Time (months)",
    ylab = "OS",
    title = "Validation Cohort-ITHScoreG",
    pval = TRUE,
    risk.table = TRUE,
    xlim = c(0,60),
    break.time.by = 12,
    censor.size = 7,
    palette = c("#3c5769","#b40000")
  )

UnicoX <- function(x, data){
  FML <- as.formula(paste0("Surv(OS, status_os) ~ ", x))
  coxph_model <- survival::coxph(FML, data = data)
  diff_data <- summary(coxph_model)
  HR <- round(diff_data$coefficients[,2], 3)
  HR_L <- round(diff_data$conf.int[,3], 3)
  HR_U <- round(diff_data$conf.int[,4], 3)
  Unicox <- data.frame(
    "Characteristics" = rownames(diff_data$coefficients),
    "HR(95%CI)" = paste(HR, "(", HR_L, "-", HR_U, ")"),
    "P value" = round(diff_data$coefficients[,"Pr(>|z|)"], 3),
    "var_name" = x
  )
  return(Unicox)
}
uni_names <- colnames(train_data)[c(2:9,29)]
uni_univar <- lapply(uni_names, UnicoX, data = train_data)
uni_univar <- do.call(rbind, uni_univar)
sig_vars <- unique(uni_univar[uni_univar$P.value<0.05,"var_name"])

traindata <- train_data[,c(sig_vars,"OS","status_os")]
external_model_data <- external_data[,c(sig_vars,"OS","status_os")]
validation_data <- validation_data[,c(sig_vars,"OS","status_os")]

cat("\n========== 开始训练 Cox 模型 ==========\n")

datarecipe_coxph <- recipe(OS + status_os ~ ., traindata) %>%
  prep()
datarecipe_coxph

traindata2_coxph <- bake(datarecipe_coxph, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
external_data2_coxph <- bake(datarecipe_coxph, new_data = external_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_coxph <- as_task_surv(
  traindata2_coxph, 
  time = "OS",
  event = "status_os", 
  type = "right"
)
task_external_coxph <- as_task_surv(
  external_data2_coxph, 
  time = "OS",
  event = "status_os", 
  type = "right"
)

learner_coxph <- lrn("surv.coxph")
learner_coxph

set.seed(seednum)
learner_coxph$train(task_train_coxph)
learner_coxph$model
summary(learner_coxph$model)

predtrain_coxph <- learner_coxph$predict(task_train_coxph)
predtrain_coxph$score(msrs(c("surv.cindex")))

predprobtrain_coxph <- 
  predtrain_coxph$distr[1:nrow(traindata2_coxph)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "coxph",
         dataset = "train",
         time = traindata2_coxph$OS,
         status_os = traindata2_coxph$status_os)

evaltrain_coxph <- eval4sa(
  predprob = predprobtrain_coxph,
  preddata = traindata2_coxph,
  etime = "OS",
  estatus = "status_os",
  model = "coxph",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predexternal_coxph <- learner_coxph$predict(task_external_coxph)
predprobexternal_coxph <- 
  predexternal_coxph$distr[1:nrow(external_data2_coxph)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "coxph",
         dataset = "external",
         time = external_data2_coxph$OS,
         status_os = external_data2_coxph$status_os)

evalexternal_coxph <- eval4sa(
  predprob = predprobexternal_coxph,
  preddata = external_data2_coxph,
  etime = "OS",
  estatus = "status_os",
  model = "coxph",
  dataset = "external",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)
evaltrain_coxph$auc
evalexternal_coxph$auc
save(predtrain_coxph, predprobtrain_coxph, evaltrain_coxph,
     predexternal_coxph, predprobexternal_coxph, evalexternal_coxph,
     file = file.path(model_dir, "coxph.RData"))

cat("========== Cox 模型训练完成 ==========\n\n")

cat("\n========== 开始训练 Lasso 模型 ==========\n")

datarecipe_lasso <- recipe(OS + status_os ~ ., traindata) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep()

traindata2_lasso <- bake(datarecipe_lasso, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
external_data2_lasso <- bake(datarecipe_lasso, new_data = external_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_lasso <- as_task_surv(
  traindata2_lasso, 
  time = "OS",
  event = "status_os", 
  type = "right"
)
task_external_lasso <- as_task_surv(
  external_data2_lasso, 
  time = "OS",
  event = "status_os", 
  type = "right"
)

learner_lasso_0 <- lrn(
  "surv.cv_glmnet",
  type.measure = "C",
  s = "lambda.1se"
)
learner_lasso_1 <- ppl(
  "distrcompositor",
  learner = learner_lasso_0,
  estimator = "kaplan",
  form = "ph"
)
learner_lasso <- as_learner(learner_lasso_1)

set.seed(seednum)
learner_lasso$train(task_train_lasso)
coef(learner_lasso$model$surv.cv_glmnet$model)

predtrain_lasso <- learner_lasso$predict(task_train_lasso)
predprobtrain_lasso <- 
  predtrain_lasso$distr[1:nrow(traindata2_lasso)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "lasso",
         dataset = "train",
         time = traindata2_lasso$OS,
         status = traindata2_lasso$status_os)

evaltrain_lasso <- eval4sa(
  predprob = predprobtrain_lasso,
  preddata = traindata2_lasso,
  etime = "OS",
  estatus = "status_os",
  model = "lasso",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predexternal_lasso <- learner_lasso$predict(task_external_lasso)
predprobexternal_lasso <- 
  predexternal_lasso$distr[1:nrow(external_data2_lasso)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "lasso",
         dataset = "external",
         time = external_data2_lasso$OS,
         status = external_data2_lasso$status_os)

evalexternal_lasso <- eval4sa(
  predprob = predprobexternal_lasso,
  preddata = external_data2_lasso,
  etime = "OS",
  estatus = "status_os",
  model = "lasso",
  dataset = "external",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)
evaltrain_lasso$auc
evalexternal_lasso$auc
save(predtrain_lasso, predprobtrain_lasso, evaltrain_lasso,
     predexternal_lasso, predprobexternal_lasso, evalexternal_lasso,
     file = file.path(model_dir, "lasso.RData"))

cat("========== Lasso 模型训练完成 ==========\n\n")

cat("\n========== 开始训练 决策树 模型 ==========\n")

datarecipe_rpart <- recipe(OS + status_os ~ ., traindata) %>%
  prep()

traindata2_rpart <- bake(datarecipe_rpart, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
external_data2_rpart <- bake(datarecipe_rpart, new_data = external_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_rpart <- as_task_surv(
  traindata2_rpart, 
  time = "OS",
  event = "status_os", 
  type = "right"
)
task_external_rpart <- as_task_surv(
  external_data2_rpart, 
  time = "OS",
  event = "status_os", 
  type = "right"
)

learner_rpart_0 <- ppl(
  "distrcompositor",
  learner = lrn(
    "surv.rpart",
    cp = to_tune(0.001, 0.5),
    minbucket = to_tune(5, 9)
  ),
  estimator = "kaplan",
  form = "ph"
) %>%
  as_learner()

learner_rpart <- auto_tuner(
  tuner = tnr(
    "grid_search", 
    param_resolutions = c(surv.rpart.cp = 5, 
                          surv.rpart.minbucket = 3), 
    batch_size = 4
  ),
  learner = learner_rpart_0,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none")
)

future::plan(future::sequential)

set.seed(seednum)
learner_rpart$train(task_train_rpart)
learner_rpart$tuning_result

predtrain_rpart <- learner_rpart$predict(task_train_rpart)
predprobtrain_rpart <- 
  predtrain_rpart$distr[1:nrow(traindata2_rpart)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "rpart",
         dataset = "train",
         time = traindata2_rpart$OS,
         status = traindata2_rpart$status_os)
colnames(predprobtrain_rpart)[1:length(itps)] <- as.character(itps)

evaltrain_rpart <- eval4sa(
  predprob = predprobtrain_rpart,
  preddata = traindata2_rpart,
  etime = "OS",
  estatus = "status_os",
  model = "rpart",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predexternal_rpart <- learner_rpart$predict(task_external_rpart)
predprobexternal_rpart <- 
  predexternal_rpart$distr[1:nrow(external_data2_rpart)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "rpart",
         dataset = "external",
         time = external_data2_rpart$OS,
         status = external_data2_rpart$status_os)
colnames(predprobexternal_rpart)[1:length(itps)] <- as.character(itps)

evalexternal_rpart <- eval4sa(
  predprob = predprobexternal_rpart,
  preddata = external_data2_rpart,
  etime = "OS",
  estatus = "status_os",
  model = "rpart",
  dataset = "external",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)
evaltrain_rpart$auc
evalexternal_rpart$auc
save(predtrain_rpart, predprobtrain_rpart, evaltrain_rpart,
     predexternal_rpart, predprobexternal_rpart, evalexternal_rpart,
     file = file.path(model_dir, "ctree.RData"))

cat("========== 决策树 模型训练完成 ==========\n\n")

cat("\n========== 开始训练 随机森林 模型 ==========\n")

datarecipe_rsf <- recipe(OS + status_os ~ ., traindata) %>%
  prep()

traindata2_rsf <- bake(datarecipe_rsf, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
external_data2_rsf <- bake(datarecipe_rsf, new_data = external_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_rsf <- as_task_surv(
  traindata2_rsf, 
  time = "OS",
  event = "status_os", 
  type = "right"
)
task_external_rsf <- as_task_surv(
  external_data2_rsf, 
  time = "OS",
  event = "status_os", 
  type = "right"
)

learner_rsf_0 <- lrn(
  "surv.rfsrc",
  ntree = to_tune(200, 500),
  mtry = to_tune(3, 5),
  nodesize = to_tune(15, 21),
  predict_type = "distr"
)

learner_rsf <- auto_tuner(
  tuner = tnr("grid_search", resolution = 3, batch_size = 3),
  learner = learner_rsf_0,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none")
)

future::plan(future::sequential)

set.seed(seednum)
learner_rsf$train(task_train_rsf)
learner_rsf$tuning_result

predtrain_rsf <- learner_rsf$predict(task_train_rsf)
predprobtrain_rsf <- 
  predtrain_rsf$distr[1:nrow(traindata2_rsf)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "rsf",
         dataset = "train",
         time = traindata2_rsf$OS,
         status = traindata2_rsf$status_os)

evaltrain_rsf <- eval4sa(
  predprob = predprobtrain_rsf,
  preddata = traindata2_rsf,
  etime = "OS",
  estatus = "status_os",
  model = "rsf",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predexternal_rsf <- learner_rsf$predict(task_external_rsf)
predprobexternal_rsf <- 
  predexternal_rsf$distr[1:nrow(external_data2_rsf)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "rsf",
         dataset = "external",
         time = external_data2_rsf$OS,
         status = external_data2_rsf$status_os)

evalexternal_rsf <- eval4sa(
  predprob = predprobexternal_rsf,
  preddata = external_data2_rsf,
  etime = "OS",
  estatus = "status_os",
  model = "rsf",
  dataset = "external",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)
save(predtrain_rsf, predprobtrain_rsf, evaltrain_rsf,
     predexternal_rsf, predprobexternal_rsf, evalexternal_rsf,
     file = file.path(model_dir, "rsf.RData"))

cat("========== 随机森林 模型训练完成 ==========\n\n")

cat("\n========== 开始训练 GBM 模型 ==========\n")

datarecipe_gbm <- recipe(OS + status_os ~ ., traindata) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep()

traindata2_gbm <- bake(datarecipe_gbm, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
external_data2_gbm <- bake(datarecipe_gbm, new_data = external_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_gbm <- as_task_surv(
  traindata2_gbm, 
  time = "OS",
  event = "status_os", 
  type = "right"
)
task_external_gbm <- as_task_surv(
  external_data2_gbm, 
  time = "OS",
  event = "status_os", 
  type = "right"
)

learner_gbm_0 <- lrn(
  "surv.gbm",
  n.trees = to_tune(100L, 300L),
  interaction.depth = to_tune(1L, 3L),
  n.minobsinnode = to_tune(5L, 15L),
  shrinkage = to_tune(0.01, 0.05),
  
  bag.fraction = 1,  # 关键：每棵树使用全部训练样本
  cv.folds = 0L,     # 内部不再做CV，外层mlr3已经做3折
  n.cores = 1L,
  keep.data = FALSE,
  verbose = FALSE
)

learner_gbm_1 <- ppl(
  "distrcompositor",
  learner = learner_gbm_0,
  estimator = "kaplan",
  form = "ph"
)

learner_gbm_2 <- as_learner(learner_gbm_1)
learner_gbm <- auto_tuner(
  tuner = tnr(
    "grid_search",
    resolution = 2,
    batch_size = 1
  ),
  learner = learner_gbm_2,
  resampling = rsmp("cv", folds = 3),
  measure = msr("surv.cindex"),
  terminator = trm("none")
)

future::plan(future::sequential)

set.seed(seednum)

system.time({
  learner_gbm$train(task_train_gbm)
})

learner_gbm$tuning_result
learner_gbm$archive$data

predtrain_gbm <- learner_gbm$predict(task_train_gbm)
predprobtrain_gbm <- 
  predtrain_gbm$distr[1:nrow(traindata2_gbm)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "gbm",
         dataset = "train",
         time = traindata2_gbm$OS,
         status = traindata2_gbm$status_os)

evaltrain_gbm <- eval4sa(
  predprob = predprobtrain_gbm,
  preddata = traindata2_gbm,
  etime = "OS",
  estatus = "status_os",
  model = "gbm",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predexternal_gbm <- learner_gbm$predict(task_external_gbm)
predprobexternal_gbm <- 
  predexternal_gbm$distr[1:nrow(external_data2_gbm)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "gbm",
         dataset = "external",
         time = external_data2_gbm$OS,
         status = external_data2_gbm$status_os)

evalexternal_gbm <- eval4sa(
  predprob = predprobexternal_gbm,
  preddata = external_data2_gbm,
  etime = "OS",
  estatus = "status_os",
  model = "gbm",
  dataset = "external",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

evaltrain_rsf$auc
evalexternal_rsf$auc

save(predtrain_gbm, predprobtrain_gbm, evaltrain_gbm,
     predexternal_gbm, predprobexternal_gbm, evalexternal_gbm,
     file = file.path(model_dir, "gbm.RData"))

cat("========== GBM 模型训练完成 ==========\n\n")

cat("\n========== 开始训练 SVM 模型 ==========\n")

datarecipe_svm <- recipe(OS + status_os ~ ., traindata) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_normalize(all_predictors()) %>%
  prep()

traindata2_svm <- bake(datarecipe_svm, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
external_data2_svm <- bake(datarecipe_svm, new_data = external_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_svm <- as_task_surv(
  traindata2_svm, 
  time = "OS",
  event = "status_os", 
  type = "right"
)
task_external_svm <- as_task_surv(
  external_data2_svm, 
  time = "OS",
  event = "status_os", 
  type = "right"
)

learner_svm_0 <- lrn(
  "surv.svm",
  type = "vanbelle1", 
  diff.meth = "makediff1",
  kernel = "rbf_kernel",
  opt.meth = "ipop",
  gamma.mu = to_tune(p_dbl(lower = 1e-2, upper = 1))
)

learner_svm_1 <- ppl(
  "distrcompositor",
  learner = learner_svm_0,
  estimator = "kaplan",
  form = "ph"
)

learner_svm_2 <- as_learner(learner_svm_1)

learner_svm <- auto_tuner(
  tuner = tnr("random_search", batch_size = 4),
  learner = learner_svm_2,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("evals", n_evals = 10)
)

future::plan(future::sequential)

set.seed(seednum)
learner_svm$train(task_train_svm)
learner_svm$tuning_result

predtrain_svm <- learner_svm$predict(task_train_svm)
predprobtrain_svm <- 
  predtrain_svm$distr[1:nrow(traindata2_svm)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "svm",
         dataset = "train",
         time = traindata2_svm$OS,
         status = traindata2_svm$status_os)

evaltrain_svm <- eval4sa(
  predprob = predprobtrain_svm,
  preddata = traindata2_svm,
  etime = "OS",
  estatus = "status_os",
  model = "svm",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predexternal_svm <- learner_svm$predict(task_external_svm)
predprobexternal_svm <- 
  predexternal_svm$distr[1:nrow(external_data2_svm)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "svm",
         dataset = "external",
         time = external_data2_svm$OS,
         status = external_data2_svm$status_os)

evalexternal_svm <- eval4sa(
  predprob = predprobexternal_svm,
  preddata = external_data2_svm,
  etime = "OS",
  estatus = "status_os",
  model = "svm",
  dataset = "external",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)
evaltrain_svm$auc
evalexternal_svm$auc
save(predtrain_svm, predprobtrain_svm, evaltrain_svm,
     predexternal_svm, predprobexternal_svm, evalexternal_svm,
     file = file.path(model_dir, "svm.RData"))

cat("========== SVM 模型训练完成 ==========\n\n")

cat("\n========== 开始训练 XGBoost 模型 ==========\n")

datarecipe_xgboost <- recipe(OS + status_os ~ ., traindata) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep()

traindata2_xgboost <- bake(datarecipe_xgboost, new_data = NULL) %>%
  dplyr::select(OS, status_os, everything())
external_data2_xgboost <- bake(datarecipe_xgboost, new_data = external_model_data) %>%
  dplyr::select(OS, status_os, everything())

task_train_xgboost <- as_task_surv(
  traindata2_xgboost, 
  time = "OS",
  event = "status_os", 
  type = "right"
)
task_external_xgboost <- as_task_surv(
  external_data2_xgboost, 
  time = "OS",
  event = "status_os", 
  type = "right"
)

learner_xgboost_0 <- lrn(
  "surv.xgboost",
  nrounds = to_tune(100, 500), 
  max_depth = to_tune(1, 5),
  eta = to_tune(1e-4, 1)
)

learner_xgboost_1 <- ppl(
  "distrcompositor",
  learner = learner_xgboost_0,
  estimator = "kaplan",
  form = "ph"
)

learner_xgboost_2 <- as_learner(learner_xgboost_1)

learner_xgboost <- auto_tuner(
  tuner = tnr("random_search", batch_size = 4),
  learner = learner_xgboost_2,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("evals", n_evals = 40)
)

future::plan(future::sequential)

set.seed(seednum)
learner_xgboost$train(task_train_xgboost)
learner_xgboost$tuning_result

predtrain_xgboost <- learner_xgboost$predict(task_train_xgboost)
predprobtrain_xgboost <- 
  predtrain_xgboost$distr[1:nrow(traindata2_xgboost)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "xgboost",
         dataset = "train",
         time = traindata2_xgboost$OS,
         status = traindata2_xgboost$status_os)

evaltrain_xgboost <- eval4sa(
  predprob = predprobtrain_xgboost,
  preddata = traindata2_xgboost,
  etime = "OS",
  estatus = "status_os",
  model = "xgboost",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)

predexternal_xgboost <- learner_xgboost$predict(task_external_xgboost)
predprobexternal_xgboost <- 
  predexternal_xgboost$distr[1:nrow(external_data2_xgboost)]$survival(itps) %>%
  t() %>%
  as.data.frame() %>%
  mutate(model = "xgboost",
         dataset = "external",
         time = external_data2_xgboost$OS,
         status = external_data2_xgboost$status_os)

evalexternal_xgboost <- eval4sa(
  predprob = predprobexternal_xgboost,
  preddata = external_data2_xgboost,
  etime = "OS",
  estatus = "status_os",
  model = "xgboost",
  dataset = "external",
  timepoints = itps,
  plotcalimethod = "quantile",
  bw4nne = NULL,
  q4quantile = 1
)
evaltrain_xgboost$auc
evalexternal_xgboost$auc
save(
  predtrain_xgboost, predprobtrain_xgboost, evaltrain_xgboost,
  predexternal_xgboost, predprobexternal_xgboost, evalexternal_xgboost,
  file = file.path(model_dir, "xgboost.RData"))

external_validation_dir <- file.path(output_dir, "external_validation")
dir.create(external_validation_dir, recursive = TRUE, showWarnings = FALSE)

evaluate_external_validation <- function(
    learner, recipe_object, model_id, raw_validation_data,
    timepoints = itps
) {
  baked_data <- bake(recipe_object, new_data = raw_validation_data) %>%
    dplyr::select(OS, status_os, everything())
  
  finite_followup <- baked_data$OS[is.finite(baked_data$OS)]
  if (length(finite_followup) == 0L) {
    stop("External validation has no finite OS follow-up values.")
  }
  max_followup <- max(finite_followup)
  evaluable_timepoints <- sort(unique(
    timepoints[is.finite(timepoints) & timepoints < max_followup]
  ))
  if (length(evaluable_timepoints) == 0L) {
    stop(
      "No requested evaluation time is below the external cohort maximum follow-up (",
      round(max_followup, 3), " months)."
    )
  }
  omitted_timepoints <- setdiff(timepoints, evaluable_timepoints)
  if (length(omitted_timepoints) > 0L) {
    message(
      model_id, ": external follow-up max = ", round(max_followup, 3),
      " months; evaluating ", paste(evaluable_timepoints, collapse = ", "),
      " months and omitting ", paste(omitted_timepoints, collapse = ", "),
      " months."
    )
  }
  
  validation_task <- as_task_surv(
    baked_data,
    time = "OS",
    event = "status_os",
    type = "right"
  )
  prediction <- learner$predict(validation_task)
  probability <- prediction$distr[seq_len(nrow(baked_data))]$survival(evaluable_timepoints) %>%
    t() %>%
    as.data.frame() %>%
    mutate(
      model = model_id,
      dataset = "external",
      time = baked_data$OS,
      status_os = baked_data$status_os
    )
  colnames(probability)[seq_along(evaluable_timepoints)] <- as.character(evaluable_timepoints)
  
  evaluation <- eval4sa(
    predprob = probability,
    preddata = baked_data,
    etime = "OS",
    estatus = "status_os",
    model = model_id,
    dataset = "external",
    timepoints = evaluable_timepoints,
    plotcalimethod = "quantile",
    bw4nne = NULL,
    q4quantile = 3
  )
  
  list(
    data = baked_data,
    task = validation_task,
    prediction = prediction,
    probability = probability,
    evaluation = evaluation,
    requested_timepoints = timepoints,
    evaluable_timepoints = evaluable_timepoints,
    omitted_timepoints = omitted_timepoints,
    max_followup = max_followup
  )
}

external_validation_results <- list(
  coxph = evaluate_external_validation(
    learner_coxph, datarecipe_coxph, "coxph", validation_data
  ),
  lasso = evaluate_external_validation(
    learner_lasso, datarecipe_lasso, "lasso", validation_data
  ),
  rpart = evaluate_external_validation(
    learner_rpart, datarecipe_rpart, "rpart", validation_data
  ),
  rsf = evaluate_external_validation(
    learner_rsf, datarecipe_rsf, "rsf", validation_data
  ),
  gbm = evaluate_external_validation(
    learner_gbm, datarecipe_gbm, "gbm", validation_data
  ),
  svm = evaluate_external_validation(
    learner_svm, datarecipe_svm, "svm", validation_data
  ),
  xgboost = evaluate_external_validation(
    learner_xgboost, datarecipe_xgboost, "xgboost", validation_data
  )
)

external_validation_data2_coxph <- external_validation_results$coxph$data
external_validation_data2_gbm <- external_validation_results$gbm$data
external_validation_data2_xgboost <- external_validation_results$xgboost$data
predexternal_validation_coxph <- external_validation_results$coxph$prediction
predexternal_validation_gbm <- external_validation_results$gbm$prediction
predexternal_validation_xgboost <- external_validation_results$xgboost$prediction
evalexternal_validation_coxph <- external_validation_results$coxph$evaluation
evalexternal_validation_lasso <- external_validation_results$lasso$evaluation
evalexternal_validation_rpart <- external_validation_results$rpart$evaluation
evalexternal_validation_rsf <- external_validation_results$rsf$evaluation
evalexternal_validation_gbm <- external_validation_results$gbm$evaluation
evalexternal_validation_svm <- external_validation_results$svm$evaluation
evalexternal_validation_xgboost <- external_validation_results$xgboost$evaluation

external_auc <- bind_rows(lapply(
  external_validation_results,
  function(x) x$evaluation$auc
))
write.csv(
  external_auc,
  file.path(external_validation_dir, "external_validation_AUC.csv"),
  row.names = FALSE
)

for (model_id in names(external_validation_results)) {
  model_result <- external_validation_results[[model_id]]
  write.csv(
    model_result$probability,
    file.path(external_validation_dir, paste0(model_id, "_survival_probability.csv")),
    row.names = FALSE
  )
  write.csv(
    model_result$evaluation$brierscore,
    file.path(external_validation_dir, paste0(model_id, "_brier.csv")),
    row.names = FALSE
  )
  ggsave(
    file.path(external_validation_dir, paste0(model_id, "_ROC.pdf")),
    model_result$evaluation$rocplot,
    width = 7,
    height = 6
  )
  ggsave(
    file.path(external_validation_dir, paste0(model_id, "_calibration.pdf")),
    model_result$evaluation$calibrationplot,
    width = 7,
    height = 6
  )
}

saveRDS(
  external_validation_results,
  file.path(external_validation_dir, "external_validation_all_models.rds")
)

explain_dir <- file.path(output_dir, "explainability")
dir.create(explain_dir, recursive = TRUE, showWarnings = FALSE)

cox_summary <- summary(learner_coxph$model)
cox_explanation <- data.frame(
  variable = rownames(cox_summary$coefficients),
  beta = cox_summary$coefficients[, "coef"],
  HR = cox_summary$conf.int[, "exp(coef)"],
  CI_lower = cox_summary$conf.int[, "lower .95"],
  CI_upper = cox_summary$conf.int[, "upper .95"],
  p_value = cox_summary$coefficients[, "Pr(>|z|)"],
  row.names = NULL,
  check.names = FALSE
) %>%
  arrange(p_value)
write.csv(
  cox_explanation,
  file.path(explain_dir, "coxph_HR_95CI.csv"),
  row.names = FALSE
)

permutation_importance_surv <- function(
    learner, data, time_col = "OS", event_col = "status_os",
    model_name, dataset_name, n_repeats = 20L, seed = 26L
) {
  feature_names <- setdiff(names(data), c(time_col, event_col))
  base_task <- as_task_surv(
    data, time = time_col, event = event_col, type = "right"
  )
  cindex_measure <- msr("surv.cindex")
  baseline_cindex <- unname(
    learner$predict(base_task)$score(cindex_measure)[[1L]]
  )
  
  set.seed(seed)
  result <- lapply(feature_names, function(feature_name) {
    permuted_cindex <- replicate(n_repeats, {
      permuted_data <- data
      permuted_data[[feature_name]] <- sample(permuted_data[[feature_name]])
      permuted_task <- as_task_surv(
        permuted_data, time = time_col, event = event_col, type = "right"
      )
      unname(learner$predict(permuted_task)$score(cindex_measure)[[1L]])
    })
    
    data.frame(
      model = model_name,
      dataset = dataset_name,
      variable = feature_name,
      baseline_cindex = baseline_cindex,
      permuted_cindex = mean(permuted_cindex, na.rm = TRUE),
      importance = baseline_cindex - mean(permuted_cindex, na.rm = TRUE),
      importance_sd = stats::sd(
        baseline_cindex - permuted_cindex, na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(result) %>% arrange(desc(importance))
}

permutation_importance_external <- bind_rows(
  permutation_importance_surv(
    learner_coxph, external_validation_data2_coxph,
    model_name = "CoxPH", dataset_name = "external", seed = seednum
  ),
  permutation_importance_surv(
    learner_gbm, external_validation_data2_gbm,
    model_name = "GBM", dataset_name = "external", seed = seednum
  ),
  permutation_importance_surv(
    learner_xgboost, external_validation_data2_xgboost,
    model_name = "XGBoost", dataset_name = "external", seed = seednum
  )
)

write.csv(
  permutation_importance_external,
  file.path(explain_dir, "permutation_importance_external.csv"),
  row.names = FALSE
)

p_permutation_importance <- permutation_importance_external %>%
  mutate(
    variable = forcats::fct_reorder(variable, importance),
    model = factor(model, levels = c("CoxPH", "GBM", "XGBoost"))
  ) %>%
  ggplot(aes(x = importance, y = variable, color = model)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  geom_errorbar(
    aes(
      xmin = importance - importance_sd,
      xmax = importance + importance_sd
    ),
    position = position_dodge(width = 0.6),
    width = 0
  ) +
  labs(
    x = "Permutation importance (decrease in external C-index)",
    y = NULL,
    color = "Model",
    title = "Model-agnostic variable importance"
  ) +
  theme_bw()

ggsave(
  file.path(explain_dir, "permutation_importance_external.pdf"),
  p_permutation_importance,
  width = 8,
  height = max(
    4,
    0.35 * length(unique(permutation_importance_external$variable)) + 2
  )
)

save(
  cox_explanation,
  permutation_importance_external,
  p_permutation_importance,
  file = file.path(explain_dir, "model_explainability.RData")
)

  shap_packages <- c("kernelshap", "shapviz", "forcats", "ragg", "svglite")
missing_shap_packages <- shap_packages[
  !vapply(shap_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_shap_packages) > 0L) {
  warning(
    "SHAP analysis skipped. Install required R packages and rerun: ",
    paste(missing_shap_packages, collapse = ", "),
      "\ninstall.packages(c('kernelshap', 'shapviz', 'forcats', 'ragg', 'svglite'))"
  )
} else {
  shap_dir <- file.path(explain_dir, "SHAP")
  dir.create(shap_dir, recursive = TRUE, showWarnings = FALSE)
  
  shap_palette <- c(low = "#4B146E", mid = "#C43C63", high = "#FDB515")
  shap_model_palette <- c(CoxPH = "#3B6FB6", GBM = "#159D89", XGBoost = "#D67A2C")
  
  theme_shap_publication <- function(base_size = 8.5) {
    ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
      ggplot2::theme(
        axis.line = ggplot2::element_line(linewidth = 0.45, colour = "black"),
        axis.ticks = ggplot2::element_line(linewidth = 0.4, colour = "black"),
        axis.text = ggplot2::element_text(colour = "black"),
        axis.title = ggplot2::element_text(colour = "black"),
        legend.title = ggplot2::element_text(face = "bold"),
        legend.key.height = grid::unit(14, "mm"),
        plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
        plot.subtitle = ggplot2::element_text(colour = "#555555"),
        plot.margin = ggplot2::margin(7, 9, 7, 7)
      )
  }
  
  pretty_feature_label <- function(x) {
    x <- gsub("[._]+", " ", x)
    x <- gsub("\\s+", " ", x)
    trimws(x)
  }
  
  canonical_feature_label <- function(x) {
    sub(" (Yes|No|Low|High)$", "", pretty_feature_label(x))
  }
  
  save_shap_plot <- function(plot_object, filename, width = 8, height = 6) {
    filename_base <- tools::file_path_sans_ext(filename)
    tryCatch(
      {
        for (extension in c("pdf", "svg", "tiff", "png")) {
          output_file <- paste0(filename_base, ".", extension)
          save_args <- list(
            filename = output_file, plot = plot_object,
            width = width, height = height, units = "in",
            limitsize = FALSE, bg = "white"
          )
          if (extension %in% c("tiff", "png")) save_args$dpi <- 600
          if (extension == "tiff") save_args$compression <- "lzw"
          do.call(ggplot2::ggsave, save_args)
        }
        TRUE
      },
      error = function(e) {
        warning("Could not save SHAP plot ", filename, ": ", conditionMessage(e))
        FALSE
      }
    )
  }
  
  safe_filename <- function(x) {
    x <- gsub("[^[:alnum:]_.-]+", "_", x)
    substr(x, 1L, 120L)
  }
  
  mlr3_crank_predict <- function(object, X) {
    as.numeric(object$predict_newdata(as.data.frame(X))$crank)
  }
  
  calculate_and_plot_shap <- function(
    learner, train_data, external_data, model_name,
    time_col = "OS", event_col = "status_os",
    background_n = 50L, explain_n = 100L, seed = 26L
  ) {
    model_dir <- file.path(shap_dir, model_name)
    dependence_dir <- file.path(model_dir, "dependence")
    local_dir <- file.path(model_dir, "local")
    dir.create(dependence_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
    
    train_x <- train_data[, setdiff(names(train_data), c(time_col, event_col)), drop = FALSE]
    test_x <- external_data[, setdiff(names(external_data), c(time_col, event_col)), drop = FALSE]
    
    set.seed(seed)
    background_ids <- sample(
      seq_len(nrow(train_x)),
      size = min(background_n, nrow(train_x))
    )
    explain_ids <- if (nrow(test_x) <= explain_n) {
      seq_len(nrow(test_x))
    } else {
      sort(sample(seq_len(nrow(test_x)), explain_n))
    }
    background_x <- train_x[background_ids, , drop = FALSE]
    explain_x <- test_x[explain_ids, , drop = FALSE]
    
    set.seed(seed)
    shap_kernel <- kernelshap::kernelshap(
      object = learner,
      X = explain_x,
      bg_X = background_x,
      pred_fun = mlr3_crank_predict,
      verbose = FALSE
    )
    shap_object <- shapviz::shapviz(shap_kernel)
    
    shap_values <- as.data.frame(shap_object$S)
    shap_values <- data.frame(
      explained_row = explain_ids,
      shap_values,
      check.names = FALSE
    )
    write.csv(
      shap_values,
      file.path(model_dir, "SHAP_values.csv"),
      row.names = FALSE
    )
    
    mean_abs_shap <- data.frame(
      model = model_name,
      variable = colnames(shap_object$S),
      mean_abs_shap = colMeans(abs(shap_object$S), na.rm = TRUE),
      row.names = NULL
    ) %>%
      mutate(
        importance_share = mean_abs_shap / sum(mean_abs_shap, na.rm = TRUE)
      ) %>%
      arrange(desc(mean_abs_shap))
    write.csv(
      mean_abs_shap,
      file.path(model_dir, "SHAP_mean_absolute_importance.csv"),
      row.names = FALSE
    )
    
    feature_order <- mean_abs_shap$variable
    bar_data <- mean_abs_shap %>%
      mutate(
        variable_label = pretty_feature_label(variable),
        variable_label = factor(variable_label, levels = rev(pretty_feature_label(feature_order)))
      )
    p_bar <- ggplot2::ggplot(
      bar_data,
      ggplot2::aes(x = mean_abs_shap, y = variable_label)
    ) +
      ggplot2::geom_col(width = 0.68, fill = shap_model_palette[[model_name]]) +
      ggplot2::geom_text(
        ggplot2::aes(label = formatC(mean_abs_shap, digits = 3, format = "f")),
        hjust = 1.08, colour = "white", size = 2.8, family = "Arial"
      ) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
      ggplot2::labs(
        x = expression("Mean |SHAP value|"), y = NULL,
        title = paste0(model_name, " — global feature importance"),
        subtitle = "External validation cohort"
      ) +
      theme_shap_publication()
    save_shap_plot(
      p_bar,
      file.path(model_dir, "SHAP_importance_bar.pdf"),
      width = 6.2,
      height = max(4, 0.27 * ncol(explain_x) + 1.8)
    )
    
    set.seed(seed)
    beeswarm_display_order <- rev(feature_order)
    beeswarm_rows <- lapply(seq_along(beeswarm_display_order), function(i) {
      feature_name <- beeswarm_display_order[[i]]
      raw_value <- explain_x[[feature_name]]
      numeric_value <- if (is.numeric(raw_value)) {
        as.numeric(raw_value)
      } else {
        as.numeric(factor(raw_value, levels = unique(raw_value)))
      }
      value_range <- range(numeric_value, na.rm = TRUE)
      normalized_value <- if (!all(is.finite(value_range)) || diff(value_range) == 0) {
        rep(0.5, length(numeric_value))
      } else {
        (numeric_value - value_range[1]) / diff(value_range)
      }
      data.frame(
        variable = feature_name,
        variable_label = pretty_feature_label(feature_name),
        shap_value = as.numeric(shap_object$S[, feature_name]),
        feature_value = normalized_value,
        y = i + stats::runif(length(numeric_value), -0.22, 0.22),
        stringsAsFactors = FALSE
      )
    })
    beeswarm_data <- dplyr::bind_rows(beeswarm_rows)
    beeswarm_data$variable_label <- factor(
      beeswarm_data$variable_label,
      levels = pretty_feature_label(beeswarm_display_order)
    )
    beeswarm_data$y <- as.numeric(beeswarm_data$variable_label) +
      (beeswarm_data$y - round(beeswarm_data$y))
    
    p_beeswarm <- ggplot2::ggplot(
      beeswarm_data,
      ggplot2::aes(x = shap_value, y = y, colour = feature_value)
    ) +
      ggplot2::geom_vline(xintercept = 0, colour = "#A9A9A9", linewidth = 0.55) +
      ggplot2::geom_point(size = 1.55, alpha = 0.82, stroke = 0) +
      ggplot2::scale_y_continuous(
        breaks = seq_along(beeswarm_display_order),
        labels = pretty_feature_label(beeswarm_display_order),
        expand = ggplot2::expansion(add = 0.6)
      ) +
      ggplot2::scale_colour_gradientn(
        colours = unname(shap_palette), limits = c(0, 1),
        breaks = c(0, 1), labels = c("Low", "High"),
        name = "Feature value"
      ) +
      ggplot2::labs(
        x = "SHAP value (impact on predicted risk)", y = NULL,
        title = paste0(model_name, " — SHAP summary"),
        subtitle = "Positive values increase predicted risk; negative values decrease it"
      ) +
      theme_shap_publication() +
      ggplot2::theme(legend.position = "right")
    save_shap_plot(
      p_beeswarm,
      file.path(model_dir, "SHAP_beeswarm.pdf"),
      width = 7.2,
      height = max(4.2, 0.28 * ncol(explain_x) + 2)
    )
    
    for (feature_name in colnames(explain_x)) {
      feature_file <- safe_filename(feature_name)
      p_dependence <- shapviz::sv_dependence(
        shap_object,
        v = feature_name,
        color_var = NULL
      ) +
        ggplot2::geom_hline(yintercept = 0, colour = "#A9A9A9", linewidth = 0.45) +
        ggplot2::labs(
          title = paste0(model_name, " — ", pretty_feature_label(feature_name)),
          subtitle = "SHAP dependence in the external validation cohort",
          x = pretty_feature_label(feature_name), y = "SHAP value"
        ) +
        theme_shap_publication()
      save_shap_plot(
        p_dependence,
        file.path(dependence_dir, paste0(feature_file, "_dependence.pdf"))
      )
      
      if (ncol(explain_x) > 1L) {
        p_dependence_interaction <- tryCatch(
          shapviz::sv_dependence(
            shap_object,
            v = feature_name,
            color_var = "auto"
          ) +
            ggplot2::geom_hline(yintercept = 0, colour = "#A9A9A9", linewidth = 0.45) +
            ggplot2::labs(
              title = paste0(model_name, " — ", pretty_feature_label(feature_name)),
              subtitle = "SHAP dependence; colour shows the strongest interaction",
              x = pretty_feature_label(feature_name), y = "SHAP value",
              colour = "Interaction value"
            ) +
            ggplot2::scale_colour_gradientn(colours = unname(shap_palette)) +
            theme_shap_publication(),
          error = function(e) NULL
        )
        if (!is.null(p_dependence_interaction)) {
          save_shap_plot(
            p_dependence_interaction,
            file.path(
              dependence_dir,
              paste0(feature_file, "_dependence_interaction.pdf")
            )
          )
        }
      }
    }
    
    explained_risk <- mlr3_crank_predict(learner, explain_x)
    local_rows <- c(
      low_risk = which.min(explained_risk),
      median_risk = which.min(abs(explained_risk - stats::median(explained_risk))),
      high_risk = which.max(explained_risk)
    )
    
    local_index <- data.frame(
      model = model_name,
      explanation = names(local_rows),
      shap_row = as.integer(local_rows),
      external_row = explain_ids[as.integer(local_rows)],
      predicted_crank = explained_risk[as.integer(local_rows)],
      row.names = NULL
    )
    write.csv(
      local_index,
      file.path(local_dir, "local_explanation_rows.csv"),
      row.names = FALSE
    )
    
    for (local_name in names(local_rows)) {
      row_id <- as.integer(local_rows[[local_name]])
      p_waterfall <- shapviz::sv_waterfall(
        shap_object,
        row_id = row_id,
        max_display = ncol(explain_x)
      ) +
        ggplot2::labs(
          title = paste0(model_name, " — ", gsub("_", " ", local_name), " patient"),
          subtitle = "Local SHAP explanation of predicted risk"
        ) +
        theme_shap_publication()
      save_shap_plot(
        p_waterfall,
        file.path(local_dir, paste0(local_name, "_SHAP_waterfall.pdf"))
      )
      
      p_force <- tryCatch(
        shapviz::sv_force(
          shap_object,
          row_id = row_id,
          max_display = ncol(explain_x)
        ) +
          ggplot2::labs(
            title = paste0(model_name, " — ", gsub("_", " ", local_name), " patient"),
            subtitle = "Yellow increases risk; purple decreases risk"
          ) +
          theme_shap_publication(),
        error = function(e) NULL
      )
      if (!is.null(p_force)) {
        save_shap_plot(
          p_force,
          file.path(local_dir, paste0(local_name, "_SHAP_force.pdf")),
          width = 10,
          height = 4
        )
      }
    }
    
    saveRDS(shap_kernel, file.path(model_dir, "kernelshap_result.rds"))
    saveRDS(shap_object, file.path(model_dir, "shapviz_object.rds"))
    
    list(
      kernelshap = shap_kernel,
      shapviz = shap_object,
      importance = mean_abs_shap,
      local_index = local_index
    )
  }
  
  shap_results <- list(
    CoxPH = calculate_and_plot_shap(
      learner_coxph, traindata2_coxph, external_data2_coxph,
      model_name = "CoxPH", seed = seednum
    ),
    GBM = calculate_and_plot_shap(
      learner_gbm, traindata2_gbm, external_data2_gbm,
      model_name = "GBM", seed = seednum
    ),
    XGBoost = calculate_and_plot_shap(
      learner_xgboost, traindata2_xgboost, external_data2_xgboost,
      model_name = "XGBoost", seed = seednum
    )
  )
  
  shap_importance_all_models <- bind_rows(
    lapply(shap_results, `[[`, "importance")
  )
  write.csv(
    shap_importance_all_models,
    file.path(shap_dir, "SHAP_importance_all_models.csv"),
    row.names = FALSE
  )
  
  p_shap_comparison <- shap_importance_all_models %>%
    mutate(
      model = factor(model, levels = c("CoxPH", "GBM", "XGBoost")),
      variable = canonical_feature_label(variable),
      variable = forcats::fct_reorder(variable, importance_share)
    ) %>%
    ggplot(aes(x = importance_share, y = variable, color = model)) +
    geom_line(
      aes(group = model),
      position = position_dodge(width = 0.55),
      linewidth = 0.45, alpha = 0.45
    ) +
    geom_point(position = position_dodge(width = 0.55), size = 2.3) +
    scale_color_manual(values = shap_model_palette) +
    labs(
      x = "Within-model share of mean absolute SHAP",
      y = NULL,
      color = "Model",
      title = "SHAP importance comparison",
      subtitle = "Importance is normalized within each model"
    ) +
    theme_shap_publication() +
    theme(legend.position = "top")
  save_shap_plot(
    p_shap_comparison,
    file.path(shap_dir, "SHAP_importance_all_models.pdf"),
    height = max(
      4,
      0.35 * length(unique(shap_importance_all_models$variable)) + 2
    )
  )
  
  save(
    shap_results,
    shap_importance_all_models,
    p_shap_comparison,
    file = file.path(shap_dir, "SHAP_all_results.RData")
  )
}

cat("========== XGBoost 模型训练完成 ==========\n\n")

evaltrain_coxph$auc
evalexternal_coxph$auc
evaltrain_gbm$auc
evalexternal_gbm$auc
evaltrain_lasso$auc
evalexternal_lasso$auc
evaltrain_rpart$auc
evalexternal_rpart$auc
evaltrain_rsf$auc
evalexternal_rsf$auc
evaltrain_svm$auc
evalexternal_svm$auc
evaltrain_xgboost$auc
evalexternal_xgboost$auc

all_auc_list <- list(
  evaltrain_coxph$auc,
  evalexternal_coxph$auc,
  evaltrain_gbm$auc,
  evalexternal_gbm$auc,
  evaltrain_lasso$auc,
  evalexternal_lasso$auc,
  evaltrain_rpart$auc,
  evalexternal_rpart$auc,
  evaltrain_rsf$auc,
  evalexternal_rsf$auc,
  evaltrain_svm$auc,
  evalexternal_svm$auc,
  evaltrain_xgboost$auc,
  evalexternal_xgboost$auc,
  evalexternal_validation_coxph$auc,
  evalexternal_validation_gbm$auc,
  evalexternal_validation_lasso$auc,
  evalexternal_validation_rpart$auc,
  evalexternal_validation_rsf$auc,
  evalexternal_validation_svm$auc,
  evalexternal_validation_xgboost$auc
)
all_auc_df <- bind_rows(all_auc_list)

auc_wide_df <- copy(all_auc_df)
auc_wide_df[, col_name := paste0(times, "month", dataset)]
auc_wide_df <- dcast(
  auc_wide_df[, .(model, col_name, AUC)],
  model ~ col_name,
  value.var = "AUC"
)
auc_columns <- c(
  "model",
  "12monthtrain", "24monthtrain", "36monthtrain", "48monthtrain", "60monthtrain",
  "12monthtest", "24monthtest", "36monthtest","48monthtest", "60monthtest",
  "12monthexternal", "24monthexternal", "36monthexternal","48monthexternal", "60monthexternal"
)
for (missing_column in setdiff(auc_columns, names(auc_wide_df))) {
  auc_wide_df[, (missing_column) := NA_real_]
}
auc_wide_df <- auc_wide_df[, ..auc_columns]
numeric_auc_columns <- setdiff(names(auc_wide_df), "model")
auc_wide_df[, (numeric_auc_columns) := lapply(
  .SD, round, digits = 3
), .SDcols = numeric_auc_columns]

evaltrain_gbm$auc
evaltrain_gbm$roc
evaltrain_gbm$rocplot
evaltrain_gbm$brierscore
evaltrain_gbm$brierscoretest
evaltrain_gbm$calibration
evaltrain_gbm$calibrationplot
evalexternal_gbm$calibrationplot
evaltrain_coxph$calibrationplot
evalexternal_coxph$calibrationplot

traindata$risk_coxph <- predtrain_coxph$crank
traindata$risk_gbm <- predtrain_gbm$crank
traindata$risk_lasso <- predtrain_lasso$crank
traindata$risk_rpart <- predtrain_rpart$crank
traindata$risk_rsf <- predtrain_rsf$crank
traindata$risk_svm <- predtrain_svm$crank
traindata$risk_xgboost <- predtrain_xgboost$crank

external_model_data$risk_coxph <- predexternal_coxph$crank
external_model_data$risk_gbm <- predexternal_gbm$crank
external_model_data$risk_lasso <- predexternal_lasso$crank
external_model_data$risk_rpart <- predexternal_rpart$crank
external_model_data$risk_rsf <- predexternal_rsf$crank
external_model_data$risk_svm <- predexternal_svm$crank
external_model_data$risk_xgboost <- predexternal_xgboost$crank

validation_data$risk_coxph <- external_validation_results$coxph$prediction$crank
validation_data$risk_gbm <- external_validation_results$gbm$prediction$crank
validation_data$risk_lasso <- external_validation_results$lasso$prediction$crank
validation_data$risk_rpart <- external_validation_results$rpart$prediction$crank
validation_data$risk_rsf <- external_validation_results$rsf$prediction$crank
validation_data$risk_svm <- external_validation_results$svm$prediction$crank
validation_data$risk_xgboost <- external_validation_results$xgboost$prediction$crank

time_points <- c(6,12,18, 24,30, 36,42, 48, 54,60)
lancet_colors <- c(
  "#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", 
  "#FDAF91", "#AD002A"
)

plot_roc_table <- function(data, data_name, requested_time_points = time_points) {
  max_followup <- max(data$OS[is.finite(data$OS)], na.rm = TRUE)
  plot_time_points <- requested_time_points[requested_time_points < max_followup]
  if (length(plot_time_points) == 0L) {
    stop("No ROC evaluation time is below maximal follow-up for ", data_name, ".")
  }
  omitted_time_points <- setdiff(requested_time_points, plot_time_points)
  if (length(omitted_time_points) > 0L) {
    message(
      data_name, ": ROC omitted non-estimable time points: ",
      paste(omitted_time_points, collapse = ", "), " months."
    )
  }
  model_list <- list(
    CoxPH = data$risk_coxph,
    GBM = data$risk_gbm,
    LASSO = data$risk_lasso,
    RPart = data$risk_rpart,
    RSF = data$risk_rsf,
    SVM = data$risk_svm,
    XGBoost = data$risk_xgboost
  )
  
  auc_ci_list <- lapply(names(model_list), function(model_name) {
    roc_obj <- timeROC(
      T = data$OS, delta = data$status_os,
      marker = model_list[[model_name]], cause = 1,
      times = plot_time_points, iid = TRUE
    )
    
    ci_obj <- confint(roc_obj)
    auc_vals <- roc_obj$AUC
    ci_vals <- if (!is.null(ci_obj$CI_AUC)) ci_obj$CI_AUC else ci_obj$CB_AUC
    
    df_temp <- data.frame(
      Time = plot_time_points,
      Model = model_name,
      AUC_raw = auc_vals,
      AUC_lower = round(ci_vals[, "2.5%"] / 100, 3),
      AUC_upper = round(ci_vals[, "97.5%"] / 100, 3)
    ) %>%
      mutate(
        AUC = round(AUC_raw, 3),
        `AUC (95%CI)` = paste0(AUC, " (", AUC_lower, "–", AUC_upper, ")"),
        Time_label = paste0(Time, " months")
      )
    return(df_temp)
  })
  
  auc_ci_combined <- do.call(rbind, auc_ci_list)
  auc_ci_results_plot <- auc_ci_combined %>% dplyr::select(Model, Time, AUC)
  auc_ci_results_table <- auc_ci_combined %>% 
    dplyr::select(Model, Time_label, `AUC (95%CI)`) %>%
    pivot_wider(names_from = Time_label, values_from = `AUC (95%CI)`)
  
  roc_plot <- ggplot(auc_ci_results_plot, aes(x = Time, y = AUC, color = Model)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(values = lancet_colors) +
    scale_y_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, 0.1), expand = c(0, 0)) +
    scale_x_continuous(
      limits = range(plot_time_points) + c(-1, 1),
      breaks = plot_time_points,
      expand = c(0, 0)
    ) +
    theme_bw() +
    theme(
      text = element_text(family = "sans", size = 10),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.direction = "horizontal",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      panel.background = element_rect(fill = "white"),
      plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm")
    ) +
    labs(
      title = paste0("Time-dependent AUC of Multiple Prognostic Models (", data_name, ")"),
      x = "Time (months)",
      y = "Area Under the ROC Curve (AUC)"
    )
  
  lancet_table_theme <- ttheme_minimal(
    core = list(fg_params = list(fontsize = 9, fontfamily = "sans")),
    colhead = list(fg_params = list(fontsize = 10, fontface = "bold"), bg_params = list(fill = "#f0f0f0")),
    rowhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
  )
  auc_table <- tableGrob(
    auc_ci_results_table,
    theme = lancet_table_theme,
    rows = NULL,
    cols = c("Model", paste0(plot_time_points, " months"))
  )
  auc_table$widths <- unit(rep(1, ncol(auc_table)), "null")
  
  final_plot <- grid.arrange(roc_plot, auc_table, ncol = 1, heights = c(4, 2))
  return(final_plot)
}

p_train <- plot_roc_table(traindata, "Training Set")
print(p_train)

p_external_cohort <- plot_roc_table(external_model_data, "Validation Set")
print(p_external_cohort)

p_external <- plot_roc_table(validation_data, "External Validation Set")
print(p_external)

ggsave(file.path(figure_dir, "TRAIN_ROC.pdf"), p_train, width=18, height=10, dpi=300)
ggsave(file.path(figure_dir, "Test_ROC.pdf"), p_external_cohort, width=18, height=10, dpi=300)
ggsave(
  file.path(figure_dir, "External_Validation_ROC.pdf"),
  p_external,
  width = 18,
  height = 10,
  dpi = 300
)

lancet_colors <- c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A")
model_names <- c("CoxPH", "GBM", "LASSO", "RPart", "RSF", "SVM", "XGBoost")
risk_cols <- c("risk_coxph", "risk_gbm", "risk_lasso", "risk_rpart", "risk_rsf", "risk_svm", "risk_xgboost")

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo1 <- cph(Surv(OS, status_os) ~ risk_coxph, x=T,y=T, data=traindata, surv=T, time.inc=12)
cal1 <- calibrate(trainnomo1, cmethod='KM', method='boot', u=12, m=300, B=1000)
plot(cal1, lwd=2, lty=1, errbar.col=lancet_colors[1],
     xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of 12 months MACE",
     ylab="Actual 12 months MACE (proportion)",
     col=lancet_colors[1])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo2 <- cph(Surv(OS, status_os) ~ risk_gbm, x=T,y=T, data=traindata, surv=T, time.inc=12)
cal2 <- calibrate(trainnomo2, cmethod='KM', method='boot', u=12, m=300, B=1000)
plot(cal2, add=T, lwd=2, lty=1, errbar.col=lancet_colors[2], col=lancet_colors[2])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo3 <- cph(Surv(OS, status_os) ~ risk_lasso, x=T,y=T, data=traindata, surv=T, time.inc=12)
cal3 <- calibrate(trainnomo3, cmethod='KM', method='boot', u=12, m=300, B=1000)
plot(cal3, add=T, lwd=2, lty=1, errbar.col=lancet_colors[3], col=lancet_colors[3])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo4 <- cph(Surv(OS, status_os) ~ risk_rpart, x=T,y=T, data=traindata, surv=T, time.inc=12)
cal4 <- calibrate(trainnomo4, cmethod='KM', method='boot', u=12, m=300, B=1000)
plot(cal4, add=T, lwd=2, lty=1, errbar.col=lancet_colors[4], col=lancet_colors[4])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo5 <- cph(Surv(OS, status_os) ~ risk_rsf, x=T,y=T, data=traindata, surv=T, time.inc=12)
cal5 <- calibrate(trainnomo5, cmethod='KM', method='boot', u=12, m=300, B=1000)
plot(cal5, add=T, lwd=2, lty=1, errbar.col=lancet_colors[5], col=lancet_colors[5])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo6 <- cph(Surv(OS, status_os) ~ risk_svm, x=T,y=T, data=traindata, surv=T, time.inc=12)
cal6 <- calibrate(trainnomo6, cmethod='KM', method='boot', u=12, m=300, B=1000)
plot(cal6, add=T, lwd=2, lty=1, errbar.col=lancet_colors[6], col=lancet_colors[6])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo7 <- cph(Surv(OS, status_os) ~ risk_xgboost, x=T,y=T, data=traindata, surv=T, time.inc=12)
cal7 <- calibrate(trainnomo7, cmethod='KM', method='boot', u=12, m=300, B=1000)
plot(cal7, add=T, lwd=2, lty=1, errbar.col=lancet_colors[7], col=lancet_colors[7])

legend(0.6,0.4,
       legend = model_names,
       x.intersp=1, y.intersp=1, lty=1, lwd=3, col=lancet_colors, bty="n", seg.len=1, cex=0.8)

legend("topleft",
       legend = c("Training Cohort", "MACE"),
       x.intersp=1, y.intersp=0.8, lty=1, lwd=2, col="black", bty="n", seg.len=1, cex=0.8)

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo1_36 <- cph(Surv(OS, status_os) ~ risk_coxph, x=T,y=T, data=traindata, surv=T, time.inc=36)
cal1_36 <- calibrate(trainnomo1_36, cmethod='KM', method='boot', u=36, m=300, B=1000)
plot(cal1_36, lwd=2, lty=1, errbar.col=lancet_colors[1],
     xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of 36 months MACE",
     ylab="Actual 36 months MACE (proportion)",
     col=lancet_colors[1])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo2_36 <- cph(Surv(OS, status_os) ~ risk_gbm, x=T,y=T, data=traindata, surv=T, time.inc=36)
cal2_36 <- calibrate(trainnomo2_36, cmethod='KM', method='boot', u=36, m=300, B=1000)
plot(cal2_36, add=T, lwd=2, lty=1, errbar.col=lancet_colors[2], col=lancet_colors[2])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo3_36 <- cph(Surv(OS, status_os) ~ risk_lasso, x=T,y=T, data=traindata, surv=T, time.inc=36)
cal3_36 <- calibrate(trainnomo3_36, cmethod='KM', method='boot', u=36, m=300, B=1000)
plot(cal3_36, add=T, lwd=2, lty=1, errbar.col=lancet_colors[3], col=lancet_colors[3])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo4_36 <- cph(Surv(OS, status_os) ~ risk_rpart, x=T,y=T, data=traindata, surv=T, time.inc=36)
cal4_36 <- calibrate(trainnomo4_36, cmethod='KM', method='boot', u=36, m=300, B=1000)
plot(cal4_36, add=T, lwd=2, lty=1, errbar.col=lancet_colors[4], col=lancet_colors[4])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo5_36 <- cph(Surv(OS, status_os) ~ risk_rsf, x=T,y=T, data=traindata, surv=T, time.inc=36)
cal5_36 <- calibrate(trainnomo5_36, cmethod='KM', method='boot', u=36, m=300, B=1000)
plot(cal5_36, add=T, lwd=2, lty=1, errbar.col=lancet_colors[5], col=lancet_colors[5])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo6_36 <- cph(Surv(OS, status_os) ~ risk_svm, x=T,y=T, data=traindata, surv=T, time.inc=36)
cal6_36 <- calibrate(trainnomo6_36, cmethod='KM', method='boot', u=36, m=300, B=1000)
plot(cal6_36, add=T, lwd=2, lty=1, errbar.col=lancet_colors[6], col=lancet_colors[6])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo7_36 <- cph(Surv(OS, status_os) ~ risk_xgboost, x=T,y=T, data=traindata, surv=T, time.inc=36)
cal7_36 <- calibrate(trainnomo7_36, cmethod='KM', method='boot', u=36, m=300, B=1000)
plot(cal7_36, add=T, lwd=2, lty=1, errbar.col=lancet_colors[7], col=lancet_colors[7])

legend(0.6,0.4, legend=model_names, x.intersp=1,y.intersp=1,lty=1,lwd=3,col=lancet_colors,bty="n",seg.len=1,cex=0.8)
legend("topleft", c("Training Cohort","MACE"), x.intersp=1,y.intersp=0.8,lty=1,lwd=2,col="black",bty="n",seg.len=1,cex=0.8)

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo1_60 <- cph(Surv(OS, status_os) ~ risk_coxph, x=T,y=T, data=traindata, surv=T, time.inc=60)
cal1_60 <- calibrate(trainnomo1_60, cmethod='KM', method='boot', u=60, m=300, B=1000)
plot(cal1_60, lwd=2, lty=1, errbar.col=lancet_colors[1],
     xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of 60 months MACE",
     ylab="Actual 60 months MACE (proportion)",
     col=lancet_colors[1])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo2_60 <- cph(Surv(OS, status_os) ~ risk_gbm, x=T,y=T, data=traindata, surv=T, time.inc=60)
cal2_60 <- calibrate(trainnomo2_60, cmethod='KM', method='boot', u=60, m=300, B=1000)
plot(cal2_60, add=T, lwd=2, lty=1, errbar.col=lancet_colors[2], col=lancet_colors[2])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo3_60 <- cph(Surv(OS, status_os) ~ risk_lasso, x=T,y=T, data=traindata, surv=T, time.inc=60)
cal3_60 <- calibrate(trainnomo3_60, cmethod='KM', method='boot', u=60, m=300, B=1000)
plot(cal3_60, add=T, lwd=2, lty=1, errbar.col=lancet_colors[3], col=lancet_colors[3])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo4_60 <- cph(Surv(OS, status_os) ~ risk_rpart, x=T,y=T, data=traindata, surv=T, time.inc=60)
cal4_60 <- calibrate(trainnomo4_60, cmethod='KM', method='boot', u=60, m=300, B=1000)
plot(cal4_60, add=T, lwd=2, lty=1, errbar.col=lancet_colors[4], col=lancet_colors[4])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo5_60 <- cph(Surv(OS, status_os) ~ risk_rsf, x=T,y=T, data=traindata, surv=T, time.inc=60)
cal5_60 <- calibrate(trainnomo5_60, cmethod='KM', method='boot', u=60, m=300, B=1000)
plot(cal5_60, add=T, lwd=2, lty=1, errbar.col=lancet_colors[5], col=lancet_colors[5])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo6_60 <- cph(Surv(OS, status_os) ~ risk_svm, x=T,y=T, data=traindata, surv=T, time.inc=60)
cal6_60 <- calibrate(trainnomo6_60, cmethod='KM', method='boot', u=60, m=300, B=1000)
plot(cal6_60, add=T, lwd=2, lty=1, errbar.col=lancet_colors[6], col=lancet_colors[6])

dd <- datadist(traindata)
options(datadist = 'dd')
trainnomo7_60 <- cph(Surv(OS, status_os) ~ risk_xgboost, x=T,y=T, data=traindata, surv=T, time.inc=60)
cal7_60 <- calibrate(trainnomo7_60, cmethod='KM', method='boot', u=60, m=300, B=1000)
plot(cal7_60, add=T, lwd=2, lty=1, errbar.col=lancet_colors[7], col=lancet_colors[7])

legend(0.6,0.4, legend=model_names, x.intersp=1,y.intersp=1,lty=1,lwd=3,col=lancet_colors,bty="n",seg.len=1,cex=0.8)
legend("topleft", c("Training Cohort","MACE"), x.intersp=1,y.intersp=0.8,lty=1,lwd=2,col="black",bty="n",seg.len=1,cex=0.8)

lancet_colors <- c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A")
model_names <- c("CoxPH", "GBM", "LASSO", "RPart", "RSF", "SVM", "XGBoost")
risk_cols <- c("risk_coxph", "risk_gbm", "risk_lasso", "risk_rpart", "risk_rsf", "risk_svm", "risk_xgboost")

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo1 <- cph(Surv(OS, status_os) ~ risk_coxph, x=T,y=T, data=external_model_data, surv=T, time.inc=12)
cal1 <- calibrate(testnomo1, cmethod='KM', method='boot', u=12, m=100, B=1000)
plot(cal1, lwd=2, lty=1, errbar.col=lancet_colors[1],
     xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of 12 months MACE",
     ylab="Actual 12 months MACE (proportion)",
     col=lancet_colors[1])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo2 <- cph(Surv(OS, status_os) ~ risk_gbm, x=T,y=T, data=external_model_data, surv=T, time.inc=12)
cal2 <- calibrate(testnomo2, cmethod='KM', method='boot', u=12, m=100, B=1000)
plot(cal2, add=T, lwd=2, lty=1, errbar.col=lancet_colors[2], col=lancet_colors[2])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo3 <- cph(Surv(OS, status_os) ~ risk_lasso, x=T,y=T, data=external_model_data, surv=T, time.inc=12)
cal3 <- calibrate(testnomo3, cmethod='KM', method='boot', u=12, m=100, B=1000)
plot(cal3, add=T, lwd=2, lty=1, errbar.col=lancet_colors[3], col=lancet_colors[3])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo4 <- cph(Surv(OS, status_os) ~ risk_rpart, x=T,y=T, data=external_model_data, surv=T, time.inc=12)
cal4 <- calibrate(testnomo4, cmethod='KM', method='boot', u=12, m=100, B=1000)
plot(cal4, add=T, lwd=2, lty=1, errbar.col=lancet_colors[4], col=lancet_colors[4])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo5 <- cph(Surv(OS, status_os) ~ risk_rsf, x=T,y=T, data=external_model_data, surv=T, time.inc=12)
cal5 <- calibrate(testnomo5, cmethod='KM', method='boot', u=12, m=100, B=1000)
plot(cal5, add=T, lwd=2, lty=1, errbar.col=lancet_colors[5], col=lancet_colors[5])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo6 <- cph(Surv(OS, status_os) ~ risk_svm, x=T,y=T, data=external_model_data, surv=T, time.inc=12)
cal6 <- calibrate(testnomo6, cmethod='KM', method='boot', u=12, m=130, B=1000)
plot(cal6, add=T, lwd=2, lty=1, errbar.col=lancet_colors[6], col=lancet_colors[6])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo7 <- cph(Surv(OS, status_os) ~ risk_xgboost, x=T,y=T, data=external_model_data, surv=T, time.inc=12)
cal7 <- calibrate(testnomo7, cmethod='KM', method='boot', u=12, m=100, B=1000)
plot(cal7, add=T, lwd=2, lty=1, errbar.col=lancet_colors[7], col=lancet_colors[7])

legend(0.6,0.4,
       legend = model_names,
       x.intersp=1, y.intersp=1, lty=1, lwd=3, col=lancet_colors, bty="n", seg.len=1, cex=0.8)

legend("topleft",
       legend = c("Validation Set", "MACE"),
       x.intersp=1, y.intersp=0.8, lty=1, lwd=2, col="black", bty="n", seg.len=1, cex=0.8)

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo1_36 <- cph(Surv(OS, status_os) ~ risk_coxph, x=T,y=T, data=external_model_data, surv=T, time.inc=36)
cal1_36 <- calibrate(testnomo1_36, cmethod='KM', method='boot', u=36, m=100, B=1000)
plot(cal1_36, lwd=2, lty=1, errbar.col=lancet_colors[1],
     xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of 36 months MACE",
     ylab="Actual 36 months MACE (proportion)",
     col=lancet_colors[1])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo2_36 <- cph(Surv(OS, status_os) ~ risk_gbm, x=T,y=T, data=external_model_data, surv=T, time.inc=36)
cal2_36 <- calibrate(testnomo2_36, cmethod='KM', method='boot', u=36, m=100, B=1000)
plot(cal2_36, add=T, lwd=2, lty=1, errbar.col=lancet_colors[2], col=lancet_colors[2])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo3_36 <- cph(Surv(OS, status_os) ~ risk_lasso, x=T,y=T, data=external_model_data, surv=T, time.inc=36)
cal3_36 <- calibrate(testnomo3_36, cmethod='KM', method='boot', u=36, m=100, B=1000)
plot(cal3_36, add=T, lwd=2, lty=1, errbar.col=lancet_colors[3], col=lancet_colors[3])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo4_36 <- cph(Surv(OS, status_os) ~ risk_rpart, x=T,y=T, data=external_model_data, surv=T, time.inc=36)
cal4_36 <- calibrate(testnomo4_36, cmethod='KM', method='boot', u=36, m=100, B=1000)
plot(cal4_36, add=T, lwd=2, lty=1, errbar.col=lancet_colors[4], col=lancet_colors[4])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo5_36 <- cph(Surv(OS, status_os) ~ risk_rsf, x=T,y=T, data=external_model_data, surv=T, time.inc=36)
cal5_36 <- calibrate(testnomo5_36, cmethod='KM', method='boot', u=36, m=100, B=1000)
plot(cal5_36, add=T, lwd=2, lty=1, errbar.col=lancet_colors[5], col=lancet_colors[5])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo6_36 <- cph(Surv(OS, status_os) ~ risk_svm, x=T,y=T, data=external_model_data, surv=T, time.inc=36)
cal6_36 <- calibrate(testnomo6_36, cmethod='KM', method='boot', u=36, m=130, B=1000)
plot(cal6_36, add=T, lwd=2, lty=1, errbar.col=lancet_colors[6], col=lancet_colors[6])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo7_36 <- cph(Surv(OS, status_os) ~ risk_xgboost, x=T,y=T, data=external_model_data, surv=T, time.inc=36)
cal7_36 <- calibrate(testnomo7_36, cmethod='KM', method='boot', u=36, m=100, B=1000)
plot(cal7_36, add=T, lwd=2, lty=1, errbar.col=lancet_colors[7], col=lancet_colors[7])

legend(0.6,0.4, legend=model_names, x.intersp=1,y.intersp=1,lty=1,lwd=3,col=lancet_colors,bty="n",seg.len=1,cex=0.8)
legend("topleft", c("Validation Set","MACE"), x.intersp=1,y.intersp=0.8,lty=1,lwd=2,col="black",bty="n",seg.len=1,cex=0.8)

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo1_60 <- cph(Surv(OS, status_os) ~ risk_coxph, x=T,y=T, data=external_model_data, surv=T, time.inc=60)
cal1_60 <- calibrate(testnomo1_60, cmethod='KM', method='boot', u=60, m=100, B=1000)
plot(cal1_60, lwd=2, lty=1, errbar.col=lancet_colors[1],
     xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of 60 months MACE",
     ylab="Actual 60 months MACE (proportion)",
     col=lancet_colors[1])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo2_60 <- cph(Surv(OS, status_os) ~ risk_gbm, x=T,y=T, data=external_model_data, surv=T, time.inc=60)
cal2_60 <- calibrate(testnomo2_60, cmethod='KM', method='boot', u=60, m=100, B=1000)
plot(cal2_60, add=T, lwd=2, lty=1, errbar.col=lancet_colors[2], col=lancet_colors[2])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo3_60 <- cph(Surv(OS, status_os) ~ risk_lasso, x=T,y=T, data=external_model_data, surv=T, time.inc=60)
cal3_60 <- calibrate(testnomo3_60, cmethod='KM', method='boot', u=60, m=100, B=1000)
plot(cal3_60, add=T, lwd=2, lty=1, errbar.col=lancet_colors[3], col=lancet_colors[3])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo4_60 <- cph(Surv(OS, status_os) ~ risk_rpart, x=T,y=T, data=external_model_data, surv=T, time.inc=60)
cal4_60 <- calibrate(testnomo4_60, cmethod='KM', method='boot', u=60, m=100, B=1000)
plot(cal4_60, add=T, lwd=2, lty=1, errbar.col=lancet_colors[4], col=lancet_colors[4])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo5_60 <- cph(Surv(OS, status_os) ~ risk_rsf, x=T,y=T, data=external_model_data, surv=T, time.inc=60)
cal5_60 <- calibrate(testnomo5_60, cmethod='KM', method='boot', u=60, m=100, B=1000)
plot(cal5_60, add=T, lwd=2, lty=1, errbar.col=lancet_colors[5], col=lancet_colors[5])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo6_60 <- cph(Surv(OS, status_os) ~ risk_svm, x=T,y=T, data=external_model_data, surv=T, time.inc=60)
cal6_60 <- calibrate(testnomo6_60, cmethod='KM', method='boot', u=60, m=300, B=1000)
plot(cal6_60, add=T, lwd=2, lty=1, errbar.col=lancet_colors[6], col=lancet_colors[6])

dd <- datadist(external_model_data)
options(datadist = 'dd')
testnomo7_60 <- cph(Surv(OS, status_os) ~ risk_xgboost, x=T,y=T, data=external_model_data, surv=T, time.inc=60)
cal7_60 <- calibrate(testnomo7_60, cmethod='KM', method='boot', u=60, m=150, B=1000)
plot(cal7_60, add=T, lwd=2, lty=1, errbar.col=lancet_colors[7], col=lancet_colors[7])

legend(0.6,0.4, legend=model_names, x.intersp=1,y.intersp=1,lty=1,lwd=3,col=lancet_colors,bty="n",seg.len=1,cex=0.8)
legend("topleft", c("Validation Set","MACE"), x.intersp=1,y.intersp=0.8,lty=1,lwd=2,col="black",bty="n",seg.len=1,cex=0.8)

model_names <- c("CoxPH", "GBM", "LASSO", "RPart", "RSF", "SVM", "XGBoost")
risk_vars   <- c("risk_coxph","risk_gbm","risk_lasso","risk_rpart","risk_rsf","risk_svm","risk_xgboost")
colors      <- c("#00468B","#ED0000","#42B540","#0099B4","#925E9F","#FDAF91","#AD002A","black","grey")

dca_train_12 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=traindata),
  model.names = model_names, times=12
)
p_train_12 <- ggplot(dca_train_12) + 
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Training Set - 12 Months")

dca_train_36 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=traindata),
  model.names = model_names, times=36
)
p_train_36 <- ggplot(dca_train_36) + 
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Training Set - 36 Months")

dca_train_60 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=traindata),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=traindata),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=traindata),
  model.names = model_names, times=60
)
p_train_60 <- ggplot(dca_train_60) + 
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Training Set - 60 Months")

dca_test_12 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=external_model_data),
  model.names = model_names, times=12
)
p_test_12 <- ggplot(dca_test_12) + 
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Test Set - 12 Months")

dca_test_36 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=external_model_data),
  model.names = model_names, times=36
)
p_test_36 <- ggplot(dca_test_36) + 
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Test Set - 36 Months")

dca_test_60 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=external_model_data),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=external_model_data),
  model.names = model_names, times=60
)
p_test_60 <- ggplot(dca_test_60) + 
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Test Set - 60 Months")

dca_val_12 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=validation_data),
  model.names = model_names, times=12
)
p_val_12 <- ggplot(dca_val_12) + 
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Validation Set - 12 Months")

dca_val_36 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=validation_data),
  model.names = model_names, times=36
)
p_val_36 <- ggplot(dca_val_36) + 
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Validation Set - 36 Months")

dca_val_60 <- dca(
  coxph(Surv(OS, status_os) ~ risk_coxph,    data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_gbm,      data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_lasso,    data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_rpart,    data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_rsf,      data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_svm,      data=validation_data),
  coxph(Surv(OS, status_os) ~ risk_xgboost,  data=validation_data),
  model.names = model_names, times=60
)
p_val_60 <- ggplot(dca_val_60) + 
  scale_color_manual(values=colors) + theme_bw() +
  ggtitle("Validation Set - 60 Months")

p_train_12
p_train_36
p_train_60

p_test_12
p_test_36
p_test_60

allplots <- list(
  p_train_12, p_train_36, p_train_60,
  p_test_12,  p_test_36,  p_test_60,
  p_val_12,  p_val_36,  p_val_60
)

ok <- sapply(allplots, inherits, "ggplot")
if(all(ok)){
  combined <- grid.arrange(
    p_train_12, p_train_36, p_train_60,
    p_test_12,  p_test_36,  p_test_60,
    p_val_12,  p_val_36,  p_val_60,
    ncol=3, nrow=3,
    top="Decision Curve Analysis (7 models)"
  )
  ggsave(file.path(figure_dir, "DCA_6plots.pdf"), combined, width=18, height=10, dpi=300)
} else {
  cat("失败的图：", which(!ok), "\n")
}

Charcutpoint <- surv_cutpoint(traindata, time = "OS", event = "status_os",
                              variables = "risk_gbm")
plot(Charcutpoint)
nomocutoff <- as.numeric(Charcutpoint["risk_gbm"][[1]]["estimate"][[1]][[1]])

traindata$xgboost_group <- ifelse(traindata$risk_gbm>=nomocutoff,"High","Low")
external_model_data$xgboost_group <- ifelse(external_model_data$risk_gbm>=nomocutoff,"High","Low")
validation_data$xgboost_group <- ifelse(validation_data$risk_gbm>=nomocutoff,"High","Low")

fit <- survfit(Surv(traindata$OS,traindata$status_os)~xgboost_group,data = traindata)
fit %>% summary(times = c(12,24,36,48,60))
ggsurvplot(fit,linetype = c(1,1),
           risk.table=T,
           data = traindata,
           title="Training Cohort",
           ylab="Mace (percentage)",xlab = "Time (Months)",
           pval = T,
           xlim = c(0,60),
           break.time.by = 12,
           palette = c("#3c5769","#b40000"))
fit <- survfit(Surv(external_model_data$OS,external_model_data$status_os)~xgboost_group,data = external_model_data)
fit %>% summary(times = c(12,24,36,48,60))
ggsurvplot(fit,linetype = c(1,1),
           risk.table=T,
           data = external_model_data,
           title="Validation Cohort",
           ylab="Mace (percentage)",xlab = "Time (Months)",
           pval = T,
           xlim = c(0,60),
           break.time.by = 12,
           palette = c("#3c5769","#b40000"))
fit <- survfit(Surv(validation_data$OS,validation_data$status_os)~xgboost_group,data = validation_data)
fit %>% summary(times = c(12,24,36,48,60))
ggsurvplot(fit,linetype = c(1,1),
           risk.table=T,
           data = validation_data,
           title="Validation Cohort",
           ylab="Mace (percentage)",xlab = "Time (Months)",
           pval = T,
           xlim = c(0,60),
           break.time.by = 12,
           palette = c("#3c5769","#b40000"))

dd <- datadist(train_data)
options(datadist = "dd")

FML <- as.formula(paste0("Surv(OS, status_os) ~ ", paste(sig_vars, collapse = " + ")))
coxph_model <- rms::cph(FML, data = train_data, x = TRUE, y = TRUE, surv = TRUE)
stepwise_model <- stats::step(coxph_model, trace = 0)
nobs <- NROW(train_data)
surv <- Survival(stepwise_model)
surv1 <- function(x) surv(1 * 12, lp = x)
surv2 <- function(x) surv(2 * 12, lp = x)
surv3 <- function(x) surv(3 * 12, lp = x)
nom <- nomogram(stepwise_model, fun = list(surv1, surv2, surv3), 
                fun.at = c(0.9, 0.8, 0.7,0.6,0.5,0.4,0.3,0.2,0.1), funlabel = c("1 year", "2 years", "3 years"))
pdf(file.path(figure_dir, "nomplot.pdf"),
    width = 18,
    height = 9)

plot(nom,
     lplabel = "Linear Predictor",
     main = "Nomogram for BMFS in SCLC Patients without PCI")

dev.off() # 必须关闭设备！！
rcorrcens(Surv(OS, status_os)~predict(stepwise_model),data = train_data)
train_data$nomosco <- TotalPoints.rms(rd = train_data, fit = stepwise_model,
                                      nom = nom)[["total points"]]
external_data$nomosco <- TotalPoints.rms(rd = external_data, fit = stepwise_model,
                                     nom = nom)[["total points"]]
validation_data$nomosco <- TotalPoints.rms(rd = validation_data, fit = stepwise_model,
                                           nom = nom)[["total points"]]
roc_train <- timeROC::timeROC(T = train_data$OS, 
                              delta = train_data$status_os, 
                              marker = train_data$nomosco, 
                              cause = 1, 
                              weighting = "marginal", 
                              times = c(12, 24, 36, 48,60), 
                              iid = TRUE)
roc_external <- timeROC::timeROC(T = external_data$OS,
                             delta = external_data$status_os, 
                             marker = external_data$nomosco, 
                             cause = 1, 
                             weighting = "marginal", 
                             times = c(12, 24, 36, 48,60), 
                             iid = TRUE)
roc_validation <- timeROC::timeROC(T = validation_data$OS,
                                   delta = validation_data$status_os, 
                                   marker = validation_data$nomosco, 
                                   cause = 1, 
                                   weighting = "marginal", 
                                   times = c(12, 24, 36, 48,60), 
                                   iid = TRUE)
auc_train_mat <- as.numeric(roc_train$AUC) 
auc_external_mat <- as.numeric(roc_external$AUC)
auc_validation_mat <- as.numeric(roc_validation$AUC)
times <- c(6,12, 18,24, 30,36,42, 48,54, 60)

train_auc <- auc_train_mat
external_auc_nomogram <- auc_external_mat
validation_auc <- auc_validation_mat

plot_df <- data.frame(
  Time = rep(times, 3),
  AUC = c(train_auc, external_auc_nomogram,validation_auc),
  Dataset = rep(c("Training Set", "External Validation Set","Validation Set"), each = 5)
)

train_labels <- paste0(
  "Training ", times, "M: ", 
  round(roc_train$AUC, 3), 
  " (95%CI: ", 
  round(confint(roc_train)[[1]][,1], 3), "-", round(confint(roc_train)[[1]][,2], 3), ")"
)

external_labels_nomogram <- paste0(
  "External ", times, "M: ", 
  round(roc_external$AUC, 3), 
  " (95%CI: ", 
  round(confint(roc_external)[[1]][,1], 3), "-", round(confint(roc_external)[[1]][,2], 3), ")"
)

validation_labels <- paste0(
  "Validation ", times, "M: ", 
  round(roc_validation$AUC, 3), 
  " (95%CI: ", 
  round(confint(roc_validation)[[1]][,1], 3), "-", round(confint(roc_validation)[[1]][,2], 3), ")"
)

all_labels <- c(train_labels, external_labels_nomogram,validation_labels)

p <- ggplot(plot_df, aes(x = Time, y = AUC, color = Dataset, group = Dataset)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  labs(
    x = "Time (Months)",
    y = "Time-Dependent AUC"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 8),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = "right", # 调整位置避免重叠
    legend.text = element_text(size = 6), # 缩小字体适配长标签
    panel.grid = element_blank(),
    legend.key = element_blank() # 去除图例背景
  ) +
  scale_color_manual(
    values = c("Training Set" = "#DC143C", "External Validation Set" = "#FFA500", "Validation Set" = "#1E90FF"),
    labels = list(
      paste(train_labels, collapse = "\n"),  # 训练集所有时间点合并为一个图例项
      paste(external_labels_nomogram, collapse = "\n"),  # 训练集所有时间点合并为一个图例项
      paste(validation_labels, collapse = "\n")     # 测试集所有时间点合并为一个图例项
    )
  ) +
  scale_x_continuous(breaks = times) +
  ylim(0.3, 1.0) +
  guides(
    color = guide_legend(
      label.position = "right",
      label.hjust = 0,
      keywidth = 1,
      keyheight = 5
    )
  )
p

auc.folds <- c("T12","T24","T36","T48","T60","V12","V24","V36","V48","V60")
auc.folds.data <- data.frame(auc.folds)
auc.folds.data$folds1 <- NA
auc.folds.data$folds2 <- NA
auc.folds.data$folds3 <- NA
auc.folds.data$folds4 <- NA
auc.folds.data$folds5 <- NA
auc.folds.data <- data.frame(t(auc.folds.data))
colnames(auc.folds.data) <- auc.folds.data[1,]
auc.folds.data <- auc.folds.data[-1,]
set.seed(21)
K_folds=5
folds.data <- train_data
folds.data$random_num <- sample(1:nrow(folds.data),nrow(folds.data),replace=F)
folds.data$group <- folds.data$random_num %% K_folds+1
model_roc <- data.frame()
for (i in min(folds.data$group):max(folds.data$group)) {
  folds.data.validation <- subset(folds.data,folds.data$group==i)
  folds.data.train <- subset(folds.data,folds.data$group!=i)
  y.train.roc <- timeROC(T = folds.data.train$OS,
                         delta = folds.data.train$status_os,
                         marker = folds.data.train$nomosco,
                         cause = 1,
                         times = c(12,24,36,48,60),
                         iid = T)
  y.validation.roc <- timeROC(T = folds.data.validation$OS,
                              delta = folds.data.validation$status_os,
                              marker = folds.data.validation$nomosco,
                              cause = 1,
                              times = c(12,24,36,48,60),
                              iid = T)
  auc.folds.data[i,1] <- round(y.train.roc$AUC[[1]],3)
  auc.folds.data[i,2] <- round(y.train.roc$AUC[[2]],3)
  auc.folds.data[i,3] <- round(y.train.roc$AUC[[3]],3)
  auc.folds.data[i,4] <- round(y.train.roc$AUC[[4]],3)
  auc.folds.data[i,5] <- round(y.train.roc$AUC[[5]],3)
  auc.folds.data[i,6] <- round(y.validation.roc$AUC[[1]],3)
  auc.folds.data[i,7] <- round(y.validation.roc$AUC[[2]],3)
  auc.folds.data[i,8] <- round(y.validation.roc$AUC[[3]],3)
  auc.folds.data[i,9] <- round(y.validation.roc$AUC[[4]],3)
  auc.folds.data[i,10] <- round(y.validation.roc$AUC[[5]],3)
}
model_roc <- data.frame("T12"=as.numeric(auc.folds.data[,1]),
                        "T24"=as.numeric(auc.folds.data[,2]),
                        "T36"=as.numeric(auc.folds.data[,3]),
                        "T48"=as.numeric(auc.folds.data[,4]),
                        "T60"=as.numeric(auc.folds.data[,5]),
                        "V12"=as.numeric(auc.folds.data[,6]),
                        "V24"=as.numeric(auc.folds.data[,7]),
                        "V36"=as.numeric(auc.folds.data[,8]),
                        "V48"=as.numeric(auc.folds.data[,9]),
                        "V60"=as.numeric(auc.folds.data[,10]))
summary(model_roc)
train.model.roc <- model_roc[,1:5]
names(train.model.roc) <- c("12Months","24Months","36Months","48Months","60Months")
train.model.roc$group <- "Train"
validation_model_roc <- model_roc[,6:10]
names(validation_model_roc) <- c("12Months","24Months","36Months","48Months","60Months")
validation_model_roc$group <- "Validation"
model_roc_plot <- rbind.data.frame(train.model.roc,validation_model_roc)
model_roc_plot <- melt(model_roc_plot)
p1<-ggplot(model_roc_plot,aes(group,value,color=group,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge")+ 
  stat_summary(geom = "errorbar",fun.data = 'mean_sd', width = 0.3)+
  labs(x="5 Cross Validation AUC",y=NULL)+
  theme_prism(palette = "candy_bright",
              base_fontface = "bold", 
              base_family = "Arial", 
              base_size = 16,  
              base_line_size = 0.8, 
              axis_text_angle = 45)+ 
  scale_fill_prism(palette = "candy_bright")
p2<-p1+facet_grid(~variable,scales = 'free')
p2
p3 <- p2+geom_point(data=model_roc_plot,aes(group,value),size=2,pch=20,color="black") 
p3

ggsave(file.path(figure_dir, "crossvalidation.pdf"),
       plot = p3,
       width = 18,
       height = 9)

dd<-datadist(train_data)
options(datadist = 'dd')
trainnomo<-cph(Surv(OS,status_os)~nomosco,x=T,y=T,data = train_data,surv = T,time.inc = 12)
cal_t_24 <- calibrate(trainnomo, cmethod='KM', method='boot', u=12, m=60, B=1000)
plot(cal_t_24,lwd=2,lty=1,errbar.col=NA,
     xlim=c(0,1),ylim=c(0,1),
     xlab="预测概率",
     ylab="实际发生概率",
     col="#3B83C0")
trainnomo<-cph(Surv(OS,status_os)~nomosco,x=T,y=T,data = train_data,surv = T,time.inc = 36)
cal_t_12 <- calibrate(trainnomo, cmethod='KM', method='boot', u=36, m=60, B=1000)
plot(cal_t_12,add=T,lwd=2,lty=1,xlim=c(0,1),ylim=c(0,1),errbar.col=NA,col="#EE362B")
trainnomo<-cph(Surv(OS,status_os)~nomosco,x=T,y=T,data = train_data,surv = T,time.inc = 60)
cal_t_12 <- calibrate(trainnomo, cmethod='KM', method='boot', u=60, m=60, B=1000)
plot(cal_t_12,add=T,lwd=2,lty=1,xlim=c(0,1),ylim=c(0,1),errbar.col=NA,col="#F8CB75")
legend(0.6,0.4,c(paste("列线图（1年）"),
                 paste("列线图（3年）"),
                 paste("列线图（5年）")),
       x.intersp = 1,
       y.intersp = 1,
       lty = 5,lwd = 3,col = c("#3B83C0","#EE362B","#F8CB75"),bty = "n",
       seg.len = 1,cex = 0.9)
legend("topleft",c(paste("训练集")),x.intersp = 1,
       y.intersp = 0.8,lty = 1,lwd = 2,col = "black",
       bty = "n",seg.len = 1,cex = 0.8)

dd<-datadist(external_data)
options(datadist = 'dd')
trainnomo<-cph(Surv(OS,status_os)~nomosco,x=T,y=T,data = external_data,surv = T,time.inc = 12)
cal_t_24 <- calibrate(trainnomo, cmethod='KM', method='boot', u=12, m=27, B=1000)
plot(cal_t_24,lwd=2,lty=1,errbar.col=NA,
     xlim=c(0,1),ylim=c(0,1),
     xlab="预测概率",
     ylab="实际发生概率",
     col="#3B83C0")
trainnomo<-cph(Surv(OS,status_os)~nomosco,x=T,y=T,data = external_data,surv = T,time.inc = 36)
cal_t_12 <- calibrate(trainnomo, cmethod='KM', method='boot', u=36, m=27, B=1000)
plot(cal_t_12,add=T,lwd=2,lty=1,errbar.col=NA,col="#EE362B")
trainnomo<-cph(Surv(OS,status_os)~nomosco,x=T,y=T,data = external_data,surv = T,time.inc = 54)
cal_t_12 <- calibrate(trainnomo, cmethod='KM', method='boot', u=54, m=27, B=1000)
plot(cal_t_12,add=T,lwd=2,lty=1,errbar.col=NA,col="#F8CB75")
legend(0.6,0.4,c(paste("列线图（1年）"),
                 paste("列线图（3年）"),
                 paste("列线图（5年）")),
       x.intersp = 1,
       y.intersp = 1,
       lty = 5,lwd = 3,col = c("#3B83C0","#EE362B","#F8CB75"),bty = "n",
       seg.len = 1,cex = 0.9)
legend("topleft",c(paste("验证集")),x.intersp = 1,
       y.intersp = 0.8,lty = 1,lwd = 2,col = "black",
       bty = "n",seg.len = 1,cex = 0.8)

dd<-datadist(validation_data)
options(datadist = 'dd')
trainnomo<-cph(Surv(OS,status_os)~nomosco,x=T,y=T,data = validation_data,surv = T,time.inc = 12)
cal_t_24 <- calibrate(trainnomo, cmethod='KM', method='boot', u=12, m=27, B=1000)
plot(cal_t_24,lwd=2,lty=1,errbar.col=NA,
     xlim=c(0,1),ylim=c(0,1),
     xlab="预测概率",
     ylab="实际发生概率",
     col="#3B83C0")
trainnomo<-cph(Surv(OS,status_os)~nomosco,x=T,y=T,data = validation_data,surv = T,time.inc = 36)
cal_t_12 <- calibrate(trainnomo, cmethod='KM', method='boot', u=36, m=27, B=1000)
plot(cal_t_12,add=T,lwd=2,lty=1,errbar.col=NA,col="#EE362B")
trainnomo<-cph(Surv(OS,status_os)~nomosco,x=T,y=T,data = validation_data,surv = T,time.inc = 54)
cal_t_12 <- calibrate(trainnomo, cmethod='KM', method='boot', u=54, m=27, B=1000)
plot(cal_t_12,add=T,lwd=2,lty=1,errbar.col=NA,col="#F8CB75")

legend(0.6,0.4,c(paste("列线图（1年）"),
                 paste("列线图（3年）"),
                 paste("列线图（5年）")),
       x.intersp = 1,
       y.intersp = 1,
       lty = 5,lwd = 3,col = c("#3B83C0","#EE362B","#F8CB75"),bty = "n",
       seg.len = 1,cex = 0.9)
legend("topleft",c(paste("Validation")),x.intersp = 1,
       y.intersp = 0.8,lty = 1,lwd = 2,col = "black",
       bty = "n",seg.len = 1,cex = 0.8)

Charcutpoint <- surv_cutpoint(train_data, time = "OS", event = "status_os",
                              variables = "nomosco")
plot(Charcutpoint)
median_score <- as.numeric(Charcutpoint["nomosco"][[1]]["estimate"][[1]][[1]])
train_data$ITHScoreGG <- ifelse(train_data$nomosco >= median_score, "High", "Low")
external_data$ITHScoreGG <- ifelse(external_data$nomosco >= median_score, "High", "Low")
validation_data$ITHScoreGG <- ifelse(validation_data$nomosco >= median_score, "High", "Low")
survfit(Surv(OS,status_os)~ITHScoreGG, data = train_data) %>%
  ggsurvplot(
    conf.int = TRUE,
    xlab = "Time (months)",
    ylab = "OS",
    title = "Train Cohort-Nomogram",
    pval = TRUE,
    risk.table = TRUE,
    xlim = c(0,60),
    break.time.by = 12,
    censor.size = 7,
    palette = c("#3c5769","#b40000")
  )
survfit(Surv(OS,status_os)~ITHScoreGG, data = external_data) %>%
  ggsurvplot(
    conf.int = TRUE,
    xlab = "Time (months)",
    ylab = "OS",
    title = "External Cohort-Nomogram",
    pval = TRUE,
    risk.table = TRUE,
    xlim = c(0,60),
    break.time.by = 12,
    censor.size = 7,
    palette = c("#3c5769","#b40000")
  )
survfit(Surv(OS,status_os)~ITHScoreGG, data = validation_data) %>%
  ggsurvplot(
    conf.int = TRUE,
    xlab = "Time (months)",
    ylab = "OS",
    title = "Validation Cohort-Nomogram",
    pval = TRUE,
    risk.table = TRUE,
    xlim = c(0,60),
    break.time.by = 12,
    censor.size = 7,
    palette = c("#3c5769","#b40000")
  )



