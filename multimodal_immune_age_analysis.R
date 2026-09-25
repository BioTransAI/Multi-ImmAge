# Multimodal immune age in UK Biobank
# Integrated analysis workflow
#
# Replace the <...> input/output path placeholders and run from top to bottom.
# All model-fitting functions and clinical analyses are defined in this file.
#
# Analysis sequence:
#   1. Model configuration and nested stacking
#   2. Held-out model comparison and modality contributions
#   3. Cross-modal immune-aging concordance
#   4. Numeric phenotype associations
#   5. Disease, disease burden, and all-cause mortality
#   6. Incident-disease multiplicity and cross-clock overlap
#   7. Mortality dose response, absolute risk, and prognostic performance
#
# Statistical methods and parameter values follow the supplied analysis.
# Inputs must be preprocessed and aligned by participant identifier.

# 1. Model configuration and nested stacking ----

library(mlr3verse)
library(mlr3learners)
library(mlr3extralearners)
library(data.table)
library(parallel)
setwd("<OUTPUT_DIRECTORY>")

# 1.1 Model settings and preprocessed inputs ----

# all_data: participant IDs in row names, numeric age in the first column
# (label), and finite numeric features. immlist maps modalities to feature names.
# S sets the seed; KOUT/KIN set fold counts; PTRAIN sets the development fraction
# NC sets worker count; TLIM sets per-fit time limits; KEEP retains fitted models.
# L1/L2 define candidate learners; FUS defines the fixed fusion methods.
S <- 666L
KOUT <- 5L
KIN <- 5L
NC <- 15L
PTRAIN <- .8
KEEP <- FALSE
TLIM <- 1200
all_data <- readRDS("<ANALYSIS_DATA_RDS>")
immune_features <- readRDS("<IMMUNE_FEATURE_LIST_RDS>")
immlist <- immune_features[c("immune_hematology", "immune_metabolite", "immune_organ_hpa_malmstrom")]
names(immlist) <- c("immune_hematology", "immune_metabolomics", "immune_proteomics")
str(immlist)
L1 <- list(Ridge = lrn("regr.glmnet", alpha = 0), LinearSVR = lrn("regr.liblinear"), LightGBM = lrn("regr.lightgbm"),
  FNN = lrn("regr.fnn"), LASSO = lrn("regr.glmnet", alpha = 1))
L2 <- list(Ridge = lrn("regr.glmnet", alpha = 0), LM = lrn("regr.lm"), NotNegRidge = lrn("regr.glmnet", alpha =
  0, lower.limits = 0))
FUS <- c("EqualWeight", "MAEWeight")
t1 = Sys.time()

# 1.2 Validate learner specifications, inputs, and module feature mappings ----

normL <- function(x, s) {
  if (!is.list(x) || !length(x) || is.null(names(x)) || anyNA(names(x)) || any(!nzchar(trimws(names(x)))) ||
    anyDuplicated(names(x))) stop(s, " must be a uniquely named non-empty list")
  z <- lapply(x, \(q) tryCatch({
    if (is.character(q)) {
      if (length(q) != 1L) stop("learner ID must have length 1")
      q <- lrn(if (startsWith(q, "regr.")) q else paste0("regr.", q))
    } else if (is.list(q) && !inherits(q, "Learner")) {
      id <- q$id
      q$id <- NULL
      if (!is.character(id) || length(id) != 1L) stop("learner specification needs one id")
      id <- if (startsWith(id, "regr.")) id else paste0("regr.", id)
      q <- do.call(lrn, c(list(id), q))
    }
    if (!inherits(q, "Learner") || q$task_type != "regr") stop("not a regression learner")
    q
  }, error = identity))
  ok <- !vapply(z, inherits, logical(1), "error")
  st <- data.table(name = names(x), registered = ok, registry_error = vapply(z, \(q) if (inherits(q, "error")) conditionMessage(q) else "",
    ""))
  if (any(!ok)) warning("Skipping unregistered ", s, ": ", paste(st$name[!ok], collapse = ", "))
  z <- z[ok]
  if (!length(z)) stop("No registered learner in ", s)
  list(learners = z, status = st)
}

cap <- function(z) {
  fb <- lrn("regr.featureless")
  if ("when" %in% names(formals(z$encapsulate))) z$encapsulate("callr", fallback = fb, when = \(cond, stage) FALSE) else z$encapsulate("callr",
    fallback = fb)
  z$timeout <- c(train = TLIM, predict = TLIM)
  z
}

a1 <- normL(L1, "L1")
a2 <- normL(L2, "L2")
L1 <- lapply(a1$learners, cap)
L2 <- lapply(a2$learners, cap)
l1_registry_status <- a1$status
l2_registry_status <- a2$status
setnames(l1_registry_status, "name", "stage1")
setnames(l2_registry_status, "name", "stage2")
rm(a1, a2)
if (any(grepl("__", c(names(L1), names(L2)), fixed = TRUE)) || any(FUS %in% names(L2))) stop("Learner names cannot contain '__' or duplicate mandatory fusion names")
stopifnot(KOUT >= 2L, KIN >= 2L, NC >= 1L, PTRAIN > 0, PTRAIN < 1, TLIM > 0, names(all_data)[1L] == "label",
  is.numeric(all_data$label), !anyNA(names(all_data)), !anyDuplicated(names(all_data)), is.list(immlist),
  length(immlist) > 0L)
if (is.null(names(immlist)) || anyNA(names(immlist)) || any(!nzchar(trimws(names(immlist)))) || anyDuplicated(names(immlist))) stop("Module names are missing, blank or duplicated")
eid <- rownames(all_data)
if (is.null(eid) || anyNA(eid) || any(!nzchar(eid)) || anyDuplicated(eid)) stop("EID is missing, blank or duplicated")
ip <- rownames(installed.packages())
if (!"callr" %in% ip) stop("Missing mandatory package: callr")
l1miss <- setNames(lapply(L1, \(z) setdiff(z$packages, ip)), names(L1))
l1_package_status <- data.table(stage1 = names(L1), available = !lengths(l1miss), missing = vapply(l1miss,
  paste, "", collapse = ", "))
badL1 <- l1_package_status$stage1[!l1_package_status$available]
if (length(badL1)) warning("Skipping unavailable L1: ", paste(badL1, collapse = ", "))
L1 <- L1[l1_package_status$stage1[l1_package_status$available]]
if (!length(L1)) stop("No available L1 learner")
l2miss <- setNames(lapply(L2, \(z) setdiff(z$packages, ip)), names(L2))
l2_package_status <- data.table(stage2 = names(L2), available = !lengths(l2miss), missing = vapply(l2miss,
  paste, "", collapse = ", "))
badL2 <- l2_package_status$stage2[!l2_package_status$available]
if (length(badL2)) warning("Skipping unavailable L2: ", paste(badL2, collapse = ", "))
L2 <- L2[l2_package_status$stage2[l2_package_status$available]]
ann <- lapply(immlist, \(z) unique(as.character(z)))
mods0 <- lapply(ann, intersect, names(all_data)[ - 1L])
unmapped <- Map(setdiff, ann, mods0)
if (any(!lengths(mods0))) stop("No mapped feature: ", paste(names(mods0)[!lengths(mods0)], collapse = ", "))
used <- unique(unlist(mods0, use.names = FALSE))
fi <- sprintf("F%05d", seq_along(used))
mi <- sprintf("M%03d", seq_along(mods0))
f2i <- setNames(fi, used)
mods <- lapply(mods0, \(z) unname(f2i[z]))
names(mods) <- mi
dat <- all_data[, c("label", used), drop = FALSE]
names(dat) <- c("label", fi)
if (!all(vapply(dat, is.numeric, logical(1))) || !all(vapply(dat, \(z) all(is.finite(z)), logical(1)))) stop("label/features must be finite numeric values")
constant <- used[vapply(dat[ - 1L], \(z) length(unique(z)) < 2L, logical(1))]
if (length(constant)) warning("Constant features retained because feature selection is disabled: ", paste(constant,
  collapse = ", "))
feature_name_map <- data.table(original = used, internal = fi)
module_name_map <- data.table(original = names(mods0), internal = mi)
module_feature_map <- rbindlist(Map(\(m, o, i) data.table(module = m, original = o, internal = i), names(mods0),
  mods0, mods))
reserved <- c("eid", "label", "fold", "prediction", "immune_age", "age_gap_raw", "age_acceleration")
if (any(names(mods0) %in% reserved)) stop("Module names conflict with output columns: ", paste(intersect(names(mods0),
  reserved), collapse = ", "))
mapping <- data.table(module = names(mods0), module_internal = mi, annotated = lengths(ann), mapped = lengths(mods0),
  excluded_unmapped = lengths(unmapped))
y <- dat$label
n <- length(y)
tasks <- Map(\(z, i) TaskRegr$new(paste0("module", i), as.data.table(dat[, c("label", z), drop = FALSE]),
  target = "label"), mods, seq_along(mods))
names(tasks) <- names(mods)
stopifnot(all(vapply(tasks, \(z) identical(z$row_ids, seq_len(n)), logical(1))))

# 1.3 Define age-stratified folds, performance metrics, and learner checks ----

mkfold <- function(a, k, seed, b = 10L) {
  if (k > length(a)) stop("folds > samples")
  set.seed(seed)
  q <- unique(quantile(a, seq(0, 1, length.out = b + 1L), type = 8))
  g <- if (length(q) > 2L) cut(a, q, include.lowest = TRUE) else factor(rep(1L, length(a)))
  f <- integer(length(a))
  for (ii in split(seq_along(a), g, drop = TRUE)) {
    v <- rep_len(seq_len(k), length(ii))
    f[ii] <- v[sample.int(length(ii))]
  }
  if (any(tabulate(f, k) == 0L)) {
    v <- rep_len(seq_len(k), length(a))
    f <- v[sample.int(length(a))]
  }
  f
}

mksplit <- function(a, p, seed, b = 10L) {
  set.seed(seed)
  q <- unique(quantile(a, seq(0, 1, length.out = b + 1L), type = 8))
  g <- if (length(q) > 2L) cut(a, q, include.lowest = TRUE) else factor(rep(1L, length(a)))
  sort(unlist(lapply(split(seq_along(a), g, drop = TRUE), \(ii) ii[sample.int(length(ii), max(1L, floor(p * length(ii))))]),
    use.names = FALSE))
}

met <- function(a, p) {
  s <- sum((a - mean(a)) ^ 2)
  data.table(N = length(a), MAE = mean(abs(a - p)), PCC = if (length(a) > 1L && sd(a) > 0 && sd(p) > 0) cor(a,
    p) else NA_real_, R2 = if (s > 0) 1 - sum((a - p) ^ 2) / s else NA_real_, RMSE = sqrt(mean((a - p) ^ 2)))
}

restoreM <- function(M) {
  j <- match(colnames(M), module_name_map$internal)
  if (anyNA(j)) stop("Unknown internal module name")
  colnames(M) <- module_name_map$original[j]
  M
}

restoreS <- function(D) {
  D <- copy(D)
  D$module <- module_name_map$original[match(D$module_internal, module_name_map$internal)]
  if (anyNA(D$module)) stop("Unknown module in stage-1 results")
  setcolorder(D, c("module", setdiff(names(D), "module")))
  D
}

restoreW <- function(W) {
  if (is.null(W) || !nrow(W)) return(data.table())
  W <- copy(W)
  W$module <- module_name_map$original[match(W$module_internal, module_name_map$internal)]
  if (anyNA(W$module)) stop("Unknown module in weights")
  setcolorder(W, c("module", setdiff(names(W), "module")))
  W
}

mperf <- function(a, M, set, f = NULL, win) {
  one <- function(ii, v) {
    z <- rbindlist(lapply(seq_len(ncol(M)), \(j) cbind(data.table(set = set, module = colnames(M)[j], fold =
      v), met(a[ii], M[ii, j]))))
    if (is.null(f)) z[, stage1 := win$stage1[match(module, win$module)]] else if (is.na(v)) z[, stage1 := "fold_specific_inner_CV"] else z[,
      stage1 := win$stage1[match(paste(module, v), paste(win$module, win$fold))]]
    if (anyNA(z$stage1)) stop("Missing stage-1 winner metadata")
    setcolorder(z, c("set", "module", "stage1", "fold", setdiff(names(z), c("set", "module", "stage1", "fold"))))
    z
  }
  if (is.null(f)) one(seq_along(a), NA_integer_) else rbindlist(c(lapply(sort(unique(f)), \(v) one(which(f == v),
    v)), list(one(seq_along(a), NA_integer_))))
}

MAP <- function(x, f, ...) {
  if (NC == 1L || .Platform$OS.type == "windows") lapply(x, f, ...) else mclapply(x, f, ..., mc.cores = min(NC,
    length(x)), mc.set.seed = TRUE, mc.preschedule = FALSE)
}

badfit <- function(l) {
  e <- tryCatch(l$errors, error = \(e) NULL)
  if ((is.numeric(e) && any(e > 0)) || (!is.null(e) && !is.numeric(e) && length(e) > 0L)) return(TRUE)
  z <- tryCatch(l$log, error = \(e) NULL)
  if (is.null(z) || !NROW(z)) return(FALSE)
  if ("class" %in% names(z) && any(tolower(as.character(z$class)) == "error")) return(TRUE)
  x <- if ("msg" %in% names(z)) z$msg else if ("condition" %in% names(z)) vapply(z$condition, \(q) if (inherits(q,
    "condition")) conditionMessage(q) else paste(as.character(q), collapse = " "), "") else unlist(z, use.names =
    FALSE)
  any(grepl("using fallback learner", paste(x, collapse = " "), ignore.case = TRUE))
}

chk <- function(l) {
  if (is.null(l$model) || badfit(l)) stop("Learner failed or timed out: ", l$id, call. = FALSE)
  invisible(l)
}
trn <- function(l, t, id = NULL) {
  l$train(t, row_ids = id)
  chk(l)
}

pp <- function(l, t, id) {
  q <- l$predict(t, row_ids = id)
  chk(l)
  j <- match(id, q$row_ids)
  if (anyNA(j) || anyDuplicated(q$row_ids)) stop("Invalid/misaligned stage-1 row_ids")
  p <- q$response[j]
  if (length(p) != length(id) || any(!is.finite(p))) stop("Invalid stage-1 predictions")
  p
}

# 1.4 Select the stage-1 learner separately for each module ----

oof1 <- function(task, tr, base, inner, seed) {
  p <- numeric(length(tr))
  for (v in seq_len(KIN)) {
    va <- which(inner == v)
    l <- base$clone(deep = TRUE)
    l$reset()
    set.seed(seed + v)
    trn(l, task, tr[inner != v])
    p[va] <- pp(l, task, tr[va])
  }
  if (any(!is.finite(p))) stop("Invalid inner OOF predictions")
  p
}

best1 <- function(task, tr, te, A, inner, seed, keep = FALSE) {
  z <- setNames(lapply(seq_along(A), \(i) {
    nm <- names(A)[i]
    t0 <- proc.time()[3L]
    q <- tryCatch(oof1(task, tr, A[[i]], inner, seed + 100L * i), error = identity)
    sec <- unname(proc.time()[3L] - t0)
    if (inherits(q, "error")) list(pred = NULL, row = data.table(stage1 = nm, ok = FALSE, seconds = sec, message =
      conditionMessage(q), N = length(tr), MAE = NA_real_, PCC = NA_real_, R2 = NA_real_, RMSE = NA_real_)) else list(pred =
      q, row = cbind(data.table(stage1 = nm, ok = TRUE, seconds = sec, message = ""), met(y[tr], q)))
  }), names(A))
  sc <- rbindlist(lapply(z, \(q) q$row), fill = TRUE)
  sc$refit_ok <- NA
  sc$refit_seconds <- NA_real_
  sc$selected <- FALSE
  ok <- which(sc$ok & is.finite(sc$MAE))
  if (!length(ok)) stop("All L1 learners failed for ", task$id)
  ord <- ok[order(sc$MAE[ok], sc$RMSE[ok], - sc$PCC[ok], - sc$R2[ok], na.last = TRUE)]
  win <- NULL
  qwin <- NULL
  for (k in ord) {
    nm <- sc$stage1[k]
    t0 <- proc.time()[3L]
    q <- tryCatch({
      l <- A[[nm]]$clone(deep = TRUE)
      l$reset()
      set.seed(seed + 90000L + k)
      trn(l, task, tr)
      list(pred = pp(l, task, te), model = if (keep) l else NULL)
    }, error = identity)
    sc$refit_seconds[k] <- unname(proc.time()[3L] - t0)
    if (inherits(q, "error")) {
      sc$refit_ok[k] <- FALSE
      sc$message[k] <- paste0("Refit: ", conditionMessage(q))
    } else {
      sc$refit_ok[k] <- TRUE
      win <- nm
      qwin <- q
      break
    }
  }
  if (is.null(win)) stop("All ranked L1 refits failed for ", task$id)
  sc$selected <- sc$stage1 == win
  list(train = z[[win]]$pred, test = qwin$pred, winner = win, scores = sc, model = qwin$model)
}

fitmods <- function(tr, te, A = L1, inner = mkfold(y[tr], KIN, 1L), seed = 1L, keep = FALSE) {
  z <- MAP(seq_along(tasks), \(m) tryCatch(best1(tasks[[m]], tr, te, A, inner, seed + 1000L * m, keep), error =
    identity))
  okz <- vapply(z, \(q) is.list(q) && !inherits(q, "condition") && all(c("train", "test", "scores") %in% names(q)),
    logical(1))
  if (any(!okz)) {
    k <- which(!okz)[1L]
    msg <- if (inherits(z[[k]], "condition")) conditionMessage(z[[k]]) else paste(as.character(z[[k]]), collapse =
      " ")
    stop("Module ", names(tasks)[k], ": ", msg)
  }
  z0 <- do.call(cbind, lapply(z, \(q) q$train))
  z1 <- do.call(cbind, lapply(z, \(q) q$test))
  colnames(z0) <- colnames(z1) <- names(tasks)
  sc <- rbindlist(Map(\(q, m) {
    x <- copy(q$scores)
    x$module_internal <- m
    x
  }, z, names(tasks)), fill = TRUE)
  setcolorder(sc, c("module_internal", setdiff(names(sc), "module_internal")))
  list(train = z0, test = z1, scores = sc, selection = sc[sc$selected == TRUE], models = if (keep) setNames(lapply(z,
    \(q) q$model), names(tasks)) else NULL)
}

# 1.5 Integrate module predictions using fusion rules or stage-2 learners ----

fit2 <- function(z0, z1, a, base, seed, keep = FALSE) {
  set.seed(seed)
  t <- TaskRegr$new("stage2", as.data.table(cbind(label = a, z0)), target = "label")
  l <- base$clone(deep = TRUE)
  l$reset()
  trn(l, t)
  nd <- as.data.table(z1)
  setcolorder(nd, t$feature_names)
  q <- l$predict_newdata(nd, task = t)
  chk(l)
  j <- match(seq_len(nrow(z1)), q$row_ids)
  if (anyNA(j) || anyDuplicated(q$row_ids)) stop("Invalid/misaligned stage-2 row_ids")
  p <- q$response[j]
  if (length(p) != nrow(z1) || any(!is.finite(p))) stop("Invalid stage-2 predictions")
  list(pred = p, model = if (keep) l else NULL)
}

fuse2 <- function(z0, z1, a, method) {
  if (!identical(colnames(z0), colnames(z1)) || any(!is.finite(z0)) || any(!is.finite(z1))) stop("Invalid/misaligned module predictions")
  e <- colMeans(abs(sweep(z0, 1L, a, "-")))
  if (method == "EqualWeight") w <- rep(1 / ncol(z0), ncol(z0)) else {
    zero <- e <= sqrt(.Machine$double.eps)
    if (any(zero)) w <- zero / sum(zero) else {
      w <- 1 / (e + sqrt(.Machine$double.eps) * max(1, median(e)))
      w <- w / sum(w)
    }
  }
  names(w) <- colnames(z0)
  list(pred = drop(z1[, names(w), drop = FALSE] %*% w), weights = data.table(module_internal = names(w), module_MAE =
    e[names(w)], weight = as.numeric(w)))
}

run2 <- function(z0, z1, a, method, B, seed, keep = FALSE) {
  if (method %in% FUS) {
    q <- fuse2(z0, z1, a, method)
    q$model <- if (keep) list(method = method, weights = q$weights) else NULL
    q
  } else {
    if (is.null(B[[method]])) stop("Unknown stage-2 method: ", method)
    q <- fit2(z0, z1, a, B[[method]], seed, keep)
    q$weights <- NULL
    q
  }
}

once <- function(tr, te, method, seed, keep = FALSE, A = L1, B = L2, inner = mkfold(y[tr], KIN, seed)) {
  z <- fitmods(tr, te, A, inner, seed, keep)
  q <- run2(z$train, z$test, y[tr], method, B, seed + 50000L, keep)
  list(pred = q$pred, module = z$test, weights = q$weights, stage1_scores = z$scores, stage1_selection = z$selection,
    models = if (keep) list(stage1 = z$models, stage2 = q$model) else NULL)
}

# 1.6 Evaluate complete pipelines with nested cross-validation ----

cvstack <- function(ids, A = L1, B = L2, seed = S, keep = FALSE, methods = c(FUS, names(B))) {
  methods <- unique(methods)
  badm <- setdiff(methods, c(FUS, names(B)))
  if (length(badm)) stop("Unknown stage-2 method: ", paste(badm, collapse = ", "))
  f <- mkfold(y[ids], KOUT, seed)
  g <- data.table(stage2 = methods, pipeline = paste0("ModuleSpecific__", methods))
  P <- matrix(NA_real_, length(ids), nrow(g), dimnames = list(NULL, g$pipeline))
  M <- matrix(NA_real_, length(ids), length(tasks), dimnames = list(NULL, names(tasks)))
  S1 <- W <- E <- list()
  fm <- if (keep) vector("list", KOUT) else NULL
  for (v in seq_len(KOUT)) {
    iv <- which(f == v)
    tr <- ids[f != v]
    te <- ids[iv]
    inn <- mkfold(y[tr], KIN, seed + 100L * v)
    z <- fitmods(tr, te, A, inn, seed + 10000L * v, keep)
    M[iv, ] <- z$test
    sx <- copy(z$scores)
    sx$fold <- v
    S1[[length(S1) + 1L]] <- sx
    if (keep) fm[[v]] <- list(stage1 = z$models, stage2 = list())
    for (j in seq_along(methods)) {
      bn <- methods[j]
      q <- tryCatch(run2(z$train, z$test, y[tr], bn, B, seed + 50000L * v + j, keep), error = identity)
      if (inherits(q, "error")) {
        E[[length(E) + 1L]] <- data.table(fold = v, stage2 = bn, message = conditionMessage(q))
        next
      }
      P[iv, g$pipeline[j]] <- q$pred
      if (!is.null(q$weights)) W[[length(W) + 1L]] <- cbind(data.table(fold = v, stage2 = bn), q$weights)
      if (keep) fm[[v]]$stage2[[bn]] <- q$model
    }
  }
  if (!all(is.finite(M))) stop("Nested CV has invalid module predictions")
  ok <- vapply(seq_len(ncol(P)), \(j) all(is.finite(P[, j])), logical(1))
  if (!any(ok)) stop("All stage-2 pipelines failed")
  P <- P[, ok, drop = FALSE]
  g <- g[match(colnames(P), g$pipeline)]
  W <- if (length(W)) rbindlist(W, fill = TRUE) else data.table()
  E <- if (length(E)) rbindlist(E, fill = TRUE) else data.table()
  if (nrow(W)) W <- W[W$stage2 %in% g$stage2]
  sc <- rbindlist(S1, fill = TRUE)
  list(pred = P, module = M, fold = f, grid = g, weights = W, stage2_errors = E, stage1_scores = sc, stage1_selection =
    sc[sc$selected == TRUE], models = fm)
}

score <- function(a, P, g) rbindlist(lapply(seq_len(nrow(g)), \(j) cbind(g[as.integer(j)], met(a, P[, j]))))
scorefold <- function(a, P, g, f) rbindlist(lapply(sort(unique(f)), \(v) cbind(data.table(fold = v), score(a[f == v],
  P[f == v, , drop = FALSE], g))))

# 1.7 Select the pipeline in the development set and evaluate the held-out set ----

ids <- seq_len(n)
tr <- mksplit(y, PTRAIN, S)
te <- setdiff(ids, tr)
if (!length(te) || length(tr) < KOUT || floor(length(tr) * (KOUT - 1) / KOUT) < KIN) stop("Invalid split/fold settings")
t0 <- proc.time()[3L]
sel <- cvstack(tr, L1, L2, S, KEEP)
selection_seconds <- unname(proc.time()[3L] - t0)
leaderboard <- score(y[tr], sel$pred, sel$grid)
setorder(leaderboard, MAE, RMSE, - PCC, - R2)
leaderboard$role <- "selection_secondary"
BEST2 <- leaderboard$stage2[1L]
key <- leaderboard$pipeline[1L]
fold_performance <- scorefold(y[tr], sel$pred, sel$grid, sel$fold)
fold_performance$role <- "selection_secondary"
testfit <- once(tr, te, BEST2, S + 9000L, KEEP, L1, L2)
split_performance <- rbindlist(list(cbind(data.table(set = "train_selection_oof", role = "selection_secondary"),
  met(y[tr], sel$pred[, key])), cbind(data.table(set = "heldout_test", role = "primary"), met(y[te], testfit$pred))))
train_prediction <- data.frame(eid = eid[tr], label = y[tr], fold = sel$fold, prediction = sel$pred[, key])
test_prediction <- data.frame(eid = eid[te], label = y[te], prediction = testfit$pred)
trM <- restoreM(sel$module)
teM <- restoreM(testfit$module)
train_module_prediction <- data.frame(eid = eid[tr], label = y[tr], fold = sel$fold, trM, check.names = FALSE)
test_module_prediction <- data.frame(eid = eid[te], label = y[te], teM, check.names = FALSE)
train_stage1_scores <- restoreS(sel$stage1_scores)
train_stage1_winners <- train_stage1_scores[train_stage1_scores$selected == TRUE]
test_stage1_scores <- restoreS(testfit$stage1_scores)
test_stage1_winners <- test_stage1_scores[test_stage1_scores$selected == TRUE]
train_module_performance <- mperf(y[tr], trM, "train_selection_oof", sel$fold, train_stage1_winners)
train_module_performance$role <- "selection_secondary"
test_module_performance <- mperf(y[te], teM, "heldout_test", NULL, test_stage1_winners)
test_module_performance$role <- "heldout_descriptive_secondary"
train_stage1_frequency <- train_stage1_winners[, .(selected_outer_folds = .N), by = c("module", "stage1")]
selection_fusion_weights_all <- restoreW(sel$weights)
selection_fusion_weights <- restoreW(if (nrow(sel$weights)) sel$weights[sel$weights$stage2 == BEST2] else data.table())
test_fusion_weights <- restoreW(if (is.null(testfit$weights)) data.table() else cbind(data.table(stage2 =
  BEST2), testfit$weights))

# 1.8 Generate full-sample OOF predictions with the selected stage-2 method ----

full <- cvstack(ids, L1, L2, S + 20000L, KEEP, methods = BEST2)
stopifnot(identical(full$grid$pipeline, paste0("ModuleSpecific__", BEST2)))
fp <- full$pred[, 1L]
fmout <- restoreM(full$module)
full_performance <- rbindlist(list(scorefold(y, full$pred, full$grid, full$fold), cbind(data.table(fold =
  NA_integer_), score(y, full$pred, full$grid))), fill = TRUE)
full_performance$role <- "post_selection_secondary"
full_stage1_scores <- restoreS(full$stage1_scores)
full_stage1_winners <- full_stage1_scores[full_stage1_scores$selected == TRUE]
full_stage1_frequency <- full_stage1_winners[, .(selected_outer_folds = .N), by = c("module", "stage1")]
module_performance <- mperf(y, fmout, "full_nested_oof", full$fold, full_stage1_winners)
module_performance$role <- "post_selection_secondary"
full_module_prediction <- data.frame(eid = eid, label = y, fold = full$fold, fmout, check.names = FALSE)
qmae <- vapply(seq_len(ncol(teM)), \(j) mean(abs(y[te] - teM[, j])), 0)
stopifnot(!any(vapply(list(train_module_performance, test_module_performance, module_performance), \(z) "stage2" %in%
  names(z), TRUE)), !anyNA(test_module_performance$stage1), isTRUE(all.equal(unname(qmae), test_module_performance$MAE[match(colnames(teM),
  test_module_performance$module)], tolerance = 1e-12)))
rm(qmae)
full_fusion_weights <- restoreW(full$weights)
full_oof <- data.frame(eid = eid, label = y, fold = full$fold, immune_age = fp, age_gap_raw = fp - y, age_acceleration =
  resid(lm(fp ~ y)), fmout, check.names = FALSE)

# 1.9 Assemble and export model results and analysis provenance ----

# The held-out test is the primary evaluation. Full-cohort OOF results use the
# selected stage-2 method and are secondary, post-selection results.
# age_acceleration stores raw age residuals here; association modules standardize
# their residuals after pooling all OOF predictions.
desc <- function(x) lapply(x, \(z) list(id = z$id, params = z$param_set$values, packages = z$packages, timeout =
  z$timeout, encapsulation = z$encapsulation))
result <- list(best = list(stage1 = "module_specific_inner_CV", stage1_by_module = test_stage1_winners, stage2 =
  BEST2), mapping = mapping, excluded_unmapped = unmapped, feature_name_map = feature_name_map, module_name_map =
  module_name_map, module_feature_map = module_feature_map, module_features = mods0, registry_status = list(stage1 =
  l1_registry_status, stage2 = l2_registry_status), package_status = list(stage1 = l1_package_status, stage2 =
  l2_package_status), split = list(train_eid = eid[tr], test_eid = eid[te], selection_seconds = selection_seconds,
  stage2_errors = sel$stage2_errors, leaderboard = leaderboard, fold_performance = fold_performance, performance =
  split_performance, stage1_candidate_scores = train_stage1_scores, stage1_winners = train_stage1_winners,
  stage1_winner_frequency = train_stage1_frequency, test_stage1_candidate_scores = test_stage1_scores, test_stage1_winners =
  test_stage1_winners, train_prediction = train_prediction, test_prediction = test_prediction, train_module_prediction =
  train_module_prediction, train_module_performance = train_module_performance, test_module_prediction = test_module_prediction,
  test_module_performance = test_module_performance, fusion_weights_all = selection_fusion_weights_all, fusion_weights =
  selection_fusion_weights, test_fusion_weights = test_fusion_weights, test_models = if (KEEP) testfit$models else NULL),
  full_cv = list(oof = full_oof, performance = full_performance, module_prediction = full_module_prediction,
  module_performance = module_performance, stage1_candidate_scores = full_stage1_scores, stage1_winners =
  full_stage1_winners, stage1_winner_frequency = full_stage1_frequency, fusion_weights = full_fusion_weights,
  stage2_errors = full$stage2_errors, module_prediction_cor = cor(fmout), fold_models = full$models), config =
  list(seed = S, train_fraction = PTRAIN, outer_folds = KOUT, inner_folds = KIN, cores = NC, per_fit_timeout_seconds =
  TLIM, stage1_selection = "module-specific inner-OOF MAE", stage1_tiebreak = c("RMSE", "PCC(desc)", "R2(desc)"),
  stage2_selection = "pooled training outer-OOF MAE", metrics = c("MAE", "PCC", "R2", "RMSE"), mandatory_fusions =
  FUS, mae_weighting = "normalized inverse inner-OOF module MAE", external_feature_selection = "none", unmapped_features =
  "excluded", imputation = "none", test_role = "primary", full_cv_role = "post-selection secondary"), learners =
  list(stage1 = desc(L1), stage2 = desc(L2)), session = sessionInfo())
print(mapping)
print(l1_registry_status)
print(l1_package_status)
print(l2_registry_status)
print(l2_package_status)
print(leaderboard)
print(test_stage1_winners[, c("module", "stage1", "MAE", "PCC", "R2", "RMSE"), with = FALSE])
print(split_performance)
print(full_performance[is.na(full_performance$fold)])
saveRDS(result, "immune_age_nested_stacking_result.rds")
t2 = Sys.time()
t2 - t1

# 2. Held-out model comparison and modality contributions ----

if (!exists("result")) result <- readRDS("immune_age_nested_stacking_result.rds")
mods <- c("immune_hematology", "immune_metabolomics", "immune_proteomics")
a <- as.data.table(result$split$test_prediction)
m <- as.data.table(result$split$test_module_prediction)
stopifnot(setequal(a$eid, m$eid), !anyDuplicated(a$eid), !anyDuplicated(m$eid))
d <- merge(a[, .(eid, label, Multimodal = prediction)], m[, c("eid", "label", mods), with = FALSE], by = c("eid",
  "label"))
stopifnot(nrow(d) == nrow(a), all(is.finite(as.matrix(d[, c("label", "Multimodal", mods), with = FALSE]))))
y <- d$label
P <- as.matrix(d[, c("Multimodal", mods), with = FALSE])

# 2.1 Paired bootstrap comparison of held-out MAE ----

E <- abs(sweep(P, 1, y, "-"))
mae <- colMeans(E)
set.seed(666)
B <- 2000L
bt <- replicate(B, colMeans(E[sample.int(nrow(E), nrow(E), TRUE), , drop = FALSE]))
bd <- sweep(bt[ - 1, , drop = FALSE], 2, bt[1, ], \(s, m) m - s)
g <- data.table(module = mods, d = mae[1] - mae[ - 1], lo = apply(bd, 1, quantile, .025, names = FALSE), hi =
  apply(bd, 1, quantile, .975, names = FALSE))
best <- as.data.table(result$split$train_module_performance)[is.na(fold)][order(MAE), module][1]

# 2.2 Linear fusion and centered modality contributions ----

stopifnot(result$best$stage2 %in% c("EqualWeight", "MAEWeight", "Ridge", "LM", "NotNegRidge"))
X <- P[, mods, drop = FALSE]
Z <- cbind(Intercept = 1, X)
f <- lm.fit(Z, d$Multimodal)
names(f$coefficients) <- colnames(Z)
b <- f$coefficients[mods]
if (any(!is.finite(f$coefficients)) || f$rank < ncol(Z)) stop("Non-identifiable contributions")
pr <- drop(Z %*% f$coefficients)
er <- max(abs(pr - d$Multimodal))
if (er > 1e-6 * max(1, max(abs(d$Multimodal)))) stop("Stage-2 reconstruction failed")
C <- sweep(sweep(X, 2, colMeans(X), "-"), 2, b, "*")
u <- colMeans(abs(C))
if (!is.finite(sum(u)) || sum(u) <= 0) stop("Zero contribution")
co <- data.table(module = mods, beta = unname(b), share = 100 * u / sum(u))
co[, s := sign(beta) * share]

# 2.3 Predicted-age and age-adjusted residual correlations ----

Q <- P[, mods, drop = FALSE]
R <- list(cor(Q), cor(apply(Q, 2, \(z) resid(lm(z ~ y)))))

# 3. Cross-modal immune-aging concordance ----

R <- if (exists("result", inherits = TRUE)) get("result", inherits = TRUE) else readRDS("immune_age_nested_stacking_result.rds")
m <- c("immune_hematology", "immune_metabolomics", "immune_proteomics")
d <- as.data.table(copy(R$full_cv$oof))
stopifnot(all(c("eid", "label", "fold", "immune_age", m) %in% names(d)), !anyDuplicated(d$eid))
d <- d[, c("eid", "label", "fold", "immune_age", m), with = FALSE]
setnames(d, c("label", "immune_age", m), c("age", "predOverall", "predH", "predM", "predP"))
rr <- function(p, a) lm.fit(cbind(1, a), p)$residuals
pr <- c("predOverall", "predH", "predM", "predP")
ac <- c("accOverall", "accH", "accM", "accP")
zv <- c("zOverall", "zH", "zM", "zP")
d[, (ac) := lapply(.SD, rr, a = age), .SDcols = pr]
d[, (zv) := lapply(.SD, \(x) as.numeric(scale(x))), .SDcols = ac]
stopifnot(all(is.finite(unlist(d[, ..zv]))))
d[, `:=`(sharedAgingDirection = (zH + zM + zP) / 3, crossModalDiscordance = sqrt(((zH - zM) ^ 2 + (zH - zP) ^ 2 +
  (zM - zP) ^ 2) / 6))][, directionalConcordance := sharedAgingDirection / (1 + crossModalDiscordance)]
lev <- c("No accelerated modality", "Single-modality acceleration", "Dual-modality acceleration", "Three-modality acceleration")
d[, nAccelerated := rowSums(cbind(zH, zM, zP) >= 0)][, Pattern := factor(lev[nAccelerated + 1L], lev)]
a <- d[, .(N = .N), by = Pattern][order(Pattern)][, pct := N / sum(N)]
nm <- c(zOverall = "Integrated multimodal immune age", zH = "Hematologic immune age", zM = "Metabolomic immune age",
  zP = "Proteomic immune age")
hc <- melt(d, id.vars = "directionalConcordance", measure.vars = names(nm), variable.name = "Clock", value.name =
  "z", variable.factor = FALSE)
hc[, Clock := factor(Clock, names(nm), nm)]
cr <- hc[, {
  q <- cor.test(z, directionalConcordance, method = "pearson")
  .(N = .N, r = unname(q$estimate), lo = q$conf.int[1], hi = q$conf.int[2])
}, by = Clock]
print(a[, .(Pattern, N, Percent = round(100 * pct, 1))])
print(cr[, .(Clock, N, r, lo, hi)])

# 4. Numeric phenotype associations ----

# 4.1 Phenotype inputs and covariates ----

# Blood-collection records support the original positional covariate alignment.
blood = readRDS("<BLOOD_COLLECTION_RDS>")
cov = readRDS("<COVARIATES_RDS>")
age = readRDS("<CHRONOLOGICAL_AGE_RDS>")
com = intersect(rownames(cov), rownames(age))
cov = cbind(age = age[com, ]$label, cov[com, ])
cov = cov[, c("age", "sex", "assessment_centre", "ethnicity_5cat", "town_index", "smoking", "alcohol", "physical_activity",
  "bmi")]
x <- as.character(cov$assessment_centre)
cov$assessment_centre <- factor(fcase(x %in% c("11004", "11005"), "Scotland", x %in% c("11003", "11022", "11023"),
  "Wales", x %in% c("11009", "11017", "11027", "11010", "11014"), "NE/Yorkshire", x %in% c("11008", "11001",
  "11016", "10003", "11024", "11025"), "North West", x %in% c("11021", "11013", "11006"), "Midlands", x %in%
  c("11011", "11028", "11002", "11007", "11026"), "South", x %in% c("11012", "11020", "11018"), "London"))
phe = readRDS("<PHENOTYPE_DATA_RDS>")
if (getRversion() < "4.1.0") stop("R >= 4.1 is required")
library(data.table)
library(parallel)

# 4.2 Phenotype analysis settings ----

N_WORKERS <- 5L
MIN_N <- 100L
ALPHA <- .05
COV_ROW_ALIGNED_WITH_BLOOD <- TRUE
SPECIAL_MISSING_ALREADY_NA <- TRUE
OUT <- "immune_age_numeric_phenotype"
MOD <- c(hematologic = "immune_hematology", metabolomic = "immune_metabolomics", proteomic = "immune_proteomics")
AA <- c("AA_integrated", "AA_hematologic", "AA_metabolomic", "AA_proteomic")
CLOCK <- setNames(c("Integrated multimodal", "Hematologic", "Metabolomic", "Proteomic"), AA)
COVARS <- c("age", "sex", "assessment_centre", "ethnicity_5cat", "town_index", "smoking", "alcohol", "physical_activity",
  "bmi")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

# 4.3 Baseline phenotype checks ----

# Participant IDs are stored in row names; character columns are excluded.
stopifnot(inherits(phe, "data.frame"), inherits(cov, "data.frame"), is.list(result))
rn <- rownames(phe)
if (is.null(rn) || length(rn) != nrow(phe) || identical(rn, as.character(seq_len(nrow(phe))))) stop("phe row names are not valid participant IDs")
phe_id <- as.character(rn)
if (anyNA(phe_id) || any(!nzchar(phe_id)) || anyDuplicated(phe_id)) stop("phe row-name IDs contain missing, blank, or duplicate values")
ph <- as.data.table(phe)
fid <- function(x) {
  z <- tolower(trimws(x))
  z <- sub("^f\\.([0-9]+)\\..*$", "\\1", z, perl = TRUE)
  z <- sub("^p", "", z)
  z <- sub("_i[0-9]+(?:_a[0-9]+)?$", "", z, perl = TRUE)
  z <- sub("-[0-9]+\\.[0-9]+$", "", z)
  z[!grepl("^[0-9]+$", z)] <- NA_character_
  z
}
src <- names(ph)
is_num <- vapply(src, \(v) is.numeric(ph[[v]]), logical(1))
is_chr <- vapply(src, \(v) is.character(ph[[v]]), logical(1))
if (any(!is_num & !is_chr)) stop("Phenotype columns are neither numeric nor character: ", paste(src[!is_num &
  !is_chr], collapse = ", "))
field_map <- data.table(source_col = src, field_id = fid(src), R_class = vapply(src, \(v) paste(class(ph[[v]]),
  collapse = "/"), character(1)), is_numeric = is_num, is_character = is_chr)
bad <- field_map[is_numeric == TRUE & is.na(field_id)]
if (nrow(bad)) stop("Unrecognized field IDs in numeric columns: ", paste(bad$source_col, collapse = ", "))
field_map[, selection := fifelse(is_character == TRUE, "excluded_character", "analysed_numeric")]
cand <- field_map[selection == "analysed_numeric"]
du <- cand[duplicated(field_id) | duplicated(field_id, fromLast = TRUE)]
if (nrow(du)) stop("Multiple numeric columns map to the same field ID: ", paste(du[, paste0(field_id, "=",
  paste(source_col, collapse = "/")), by = field_id]$V1, collapse = ", "))
num <- cand$source_col
if (!length(num)) stop("No numeric phenotypes are available in phe")
if (!isTRUE(SPECIAL_MISSING_ALREADY_NA)) stop("Convert field-specific UKB missing-value codes to NA before analysis; do not remove all negative values")

# 4.4 Standardized age accelerations from pooled OOF predictions ----

o <- copy(as.data.table(result$full_cv$oof))
need <- c("eid", "label", "immune_age", unname(MOD))
if (!all(need %chin% names(o))) stop("result$full_cv$oof is missing: ", paste(setdiff(need, names(o)), collapse =
  ", "))
if (any(!vapply(need[ - 1L], \(v) all(is.finite(o[[v]])), logical(1)))) stop("OOF age labels or predictions contain NA/Inf")
z <- function(x) {
  s <- sd(x)
  if (any(!is.finite(x)) || !is.finite(s) || s <= 0) stop("Age acceleration cannot be standardized")
  as.numeric((x - mean(x)) / s)
}
rz <- function(p, a) z(lm.fit(cbind(1, a), p)$residuals)
aa <- o[, .(eid = as.character(eid), clock_age = as.numeric(label), AA_integrated = rz(immune_age, label))]
aa[, (AA[ - 1L]) := lapply(unname(MOD), \(v) rz(o[[v]], o$label))]
if (anyNA(aa) || anyDuplicated(aa$eid)) stop("OOF age data contain missing or invalid values or duplicate IDs")

# 4.5 Participant alignment and covariate design ----

cv <- copy(as.data.table(cov))
if (!all(COVARS %chin% names(cv))) stop("cov is missing: ", paste(setdiff(COVARS, names(cv)), collapse = ", "))
fac <- COVARS[vapply(COVARS, \(v) is.factor(cv[[v]]), logical(1))]
lev <- setNames(lapply(fac, \(v) levels(cv[[v]])), fac)
if (!"eid" %chin% names(cv)) {
  if (!COV_ROW_ALIGNED_WITH_BLOOD || !exists("blood", inherits = FALSE) || nrow(cv) != nrow(blood)) stop("cov has no eid; positional alignment requires confirmed row correspondence and equal row counts with blood")
  cv[, eid := as.character(blood$eid)]
} else cv[, eid := as.character(eid)]
if (anyNA(cv$eid) || anyDuplicated(cv$eid)) stop("cov$eid contains missing or duplicate IDs")
cv <- cv[, c("eid", COVARS), with = FALSE]
base <- merge(aa, cv, by = "eid", sort = FALSE)
if (nrow(base) != nrow(aa) || anyNA(base[, c("eid", "clock_age", AA, COVARS), with = FALSE])) stop("OOF and cov are not fully matched, or covariates contain missing values")
for (v in fac) set(base, j = v, value = droplevels(factor(as.character(base[[v]]), levels = lev[[v]])))
age_audit <- base[, .(N = .N, r = cor(clock_age, age), Q95 = unname(quantile(abs(clock_age - age), .95)),
  MAX = max(abs(clock_age - age)))]
if (!is.finite(age_audit$r) || age_audit$r < .99 || age_audit$Q95 > 2) stop("OOF labels and cov$age do not match; check cov/blood row order and age definitions")
ix <- match(base$eid, phe_id)
if (anyNA(ix)) {
  warning(sum(is.na(ix)), " OOF participants are absent from phe and excluded consistently across phenotypes")
  keep <- which(!is.na(ix))
  base <- base[keep]
  ix <- ix[keep]
}
if (!nrow(base)) stop("OOF and phe have no participants in common")
P <- ph[ix, ..num]
if (!identical(base$eid, phe_id[ix])) stop("Participant-ID alignment failed between phe and the analysis cohort")
COV_USE <- COVARS[vapply(COVARS, \(v) uniqueN(base[[v]]) > 1L, logical(1))]
COV_DROPPED <- setdiff(COVARS, COV_USE)
if (length(COV_DROPPED)) warning("Constant covariates cannot be adjusted for and have been recorded: ", paste(COV_DROPPED,
  collapse = ", "))
Z <- model.matrix(reformulate(COV_USE), data = base)
AV <- as.matrix(base[, ..AA])
storage.mode(AV) <- "double"

# 4.6 Phenotype-specific quality checks and eligibility ----

qc <- rbindlist(lapply(seq_along(num), \(i) {
  y <- P[[i]]
  ok <- is.finite(y)
  n <- sum(ok)
  s <- if (n > 1L) sd(y[ok]) else NA_real_
  u <- uniqueN(y[ok])
  data.table(source_col = num[i], N_mother = nrow(P), N = n, N_missing = sum(is.na(y)), N_nonfinite = sum(!is.na(y) &
    !is.finite(y)), missing_rate = 1 - n / nrow(P), n_unique = u, outcome_mean = if (n) mean(y[ok]) else NA_real_,
    outcome_sd = s, outcome_kind = fifelse(u == 2L, "binary_numeric", fifelse(u <= 10L, "discrete_numeric",
    "continuous_numeric")), eligible = n >= MIN_N & u >= 2L & is.finite(s) & s > 0)
}))
qc <- merge(qc, field_map[, .(source_col, field_id, R_class)], by = "source_col", sort = FALSE)
qc[, failure_reason := fcase(eligible == TRUE, "", N < MIN_N, paste0("N<", MIN_N), n_unique < 2L, "fewer than 2 unique values",
  !is.finite(outcome_sd) | outcome_sd <= 0, "zero/non-finite SD", default = "not estimable")]
todo <- qc[eligible == TRUE]
M <- nrow(todo)
if (!M) stop("No numeric phenotypes meet the minimum estimability criteria")

# 4.7 Separate OLS associations with HC1 robust standard errors ----

safe <- function(expr) {
  w <- character()
  v <- withCallingHandlers(tryCatch(expr, error = identity), warning = \(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  })
  list(value = v, warning = paste(unique(w[nzchar(w)]), collapse = " | "))
}
fail4 <- function(q, status, reason, w = "", rank = NA_integer_, rd = NA, alias = "") data.table(field_id =
  q$field_id, source_col = q$source_col, outcome_kind = q$outcome_kind, N = q$N, N_missing = q$N_missing,
  N_nonfinite = q$N_nonfinite, N_excluded_y = q$N_mother - q$N, missing_rate = q$missing_rate, n_unique =
  q$n_unique, outcome_mean = q$outcome_mean, outcome_sd = q$outcome_sd, exposure = AA, clock = unname(CLOCK[AA]),
  beta = NA_real_, SE_model = NA_real_, P_model = NA_real_, SE_HC1 = NA_real_, LCL = NA_real_, UCL = NA_real_,
  statistic = NA_real_, df = NA_real_, partial_R2 = NA_real_, P = NA_real_, covariate_rank = rank, rank_deficient =
  rd, aliased_covariate_columns = alias, status = status, warning = w, failure_reason = reason)
scan <- function(i) {
  setDTthreads(1L)
  q <- todo[i]
  y <- P[[q$source_col]]
  ii <- which(is.finite(y))
  ys <- q$outcome_sd
  X <- Z[ii, , drop = FALSE]
  Y <- cbind((y[ii] - q$outcome_mean) / ys, AV[ii, , drop = FALSE])
  a <- safe(lm.fit(X, Y, singular.ok = TRUE))
  if (inherits(a$value, "error")) return(fail4(q, "error", conditionMessage(a$value), a$warning))
  f <- a$value
  R <- f$residuals
  rk <- f$rank
  rd <- rk < ncol(X)
  alias <- if (rd) paste(colnames(X)[f$qr$pivot[seq.int(rk + 1L, ncol(X))]], collapse = " | ") else ""
  wx <- c(a$warning, if (rd) "covariate design matrix is rank deficient")
  w <- paste(unique(wx[nzchar(wx)]), collapse = " | ")
  n <- length(ii)
  dfs <- n - rk - 1L
  rbindlist(lapply(seq_along(AA), \(j) {
    rx <- R[, j + 1L]
    ry <- R[, 1L]
    sxx <- sum(rx ^ 2)
    if (!is.finite(sxx) || sxx <= sqrt(.Machine$double.eps) || dfs <= 0) return(fail4(q, "not_estimable",
      "exposure has no residual variation or insufficient df", w, rk, rd, alias)[exposure == AA[j]])
    b <- sum(rx * ry) / sxx
    e <- ry - b * rx
    sm <- sqrt(sum(e ^ 2) / dfs / sxx)
    sr <- sqrt((n / dfs) * sum((rx * e) ^ 2) / sxx ^ 2)
    tm <- b / sm
    tt <- b / sr
    pm <- 2 * pt(abs(tm), dfs, lower.tail = FALSE)
    pv <- 2 * pt(abs(tt), dfs, lower.tail = FALSE)
    qq <- qt(1 - ALPHA / 2, dfs)
    pr <- 1 - sum(e ^ 2) / sum(ry ^ 2)
    ok <- all(is.finite(c(b, sm, sr, tm, tt, pm, pv, pr)))
    if (!ok) return(fail4(q, "not_estimable", "non-finite target estimate", w, rk, rd, alias)[exposure == AA[j]])
    data.table(field_id = q$field_id, source_col = q$source_col, outcome_kind = q$outcome_kind, N = n, N_missing =
      q$N_missing, N_nonfinite = q$N_nonfinite, N_excluded_y = q$N_mother - n, missing_rate = q$missing_rate,
      n_unique = q$n_unique, outcome_mean = q$outcome_mean, outcome_sd = ys, exposure = AA[j], clock = unname(CLOCK[AA[j]]),
      beta = b, SE_model = sm, P_model = pm, SE_HC1 = sr, LCL = b - qq * sr, UCL = b + qq * sr, statistic =
      tt, df = dfs, partial_R2 = max(0, min(1, pr)), P = pv, covariate_rank = rk, rank_deficient = rd, aliased_covariate_columns =
      alias, status = fifelse(nzchar(w), "ok_with_warning", "ok"), warning = w, failure_reason = "")
  }))
}

# 4.8 Parallel estimation and multiple-testing correction ----

# Retain warnings and failed phenotype-clock rows.
nw <- as.integer(max(1L, min(N_WORKERS, nrow(todo), ifelse(is.na(detectCores(logical = FALSE)), N_WORKERS,
  detectCores(logical = FALSE)))))
message("Numeric phenotype scans: ", M, " phenotypes × 4 clocks; ", nw, " fork workers")
ans <- mclapply(seq_len(M), \(i) tryCatch(scan(i), error = \(e) fail4(todo[i], "error", conditionMessage(e))),
  mc.cores = nw, mc.preschedule = TRUE, mc.set.seed = FALSE)
if (any(!vapply(ans, \(x) inherits(x, "data.frame") && nrow(x) == 4L, logical(1)))) stop("At least one parallel worker failed to return four result rows")
unfit <- qc[eligible == FALSE]
unfit <- if (nrow(unfit)) rbindlist(lapply(seq_len(nrow(unfit)), \(i) fail4(unfit[i], "not_estimable", unfit$failure_reason[i]))) else data.table()
association <- rbindlist(c(ans, list(unfit)), fill = TRUE)
N_TESTS_PLANNED <- 4L * M
association[, P_bonferroni := pmin(1, P * N_TESTS_PLANNED)][, significant_bonferroni := is.finite(P_bonferroni) &
  P_bonferroni < ALPHA][, sort_missing := is.na(P)]
setorder(association, sort_missing, P_bonferroni, P, field_id, exposure)
association[, sort_missing := NULL]

# The supplied analysis ends with the association and qc objects in memory.
# Add project-specific result export here if required.

# 5. Disease, disease burden, and all-cause mortality ----

# 5.1 Clinical inputs and covariates ----

# Use the full-cohort OOF predictions stored in result by module 1.
# Hospital events require eid, icd3, and event_date; blood collection requires
# eid and blood_date; death records require eid and death_date.
# Retain verified covariate row alignment and the original factor coding.
dis = readRDS("<DISEASE_AND_DEATH_RDS>")$hosp_events
death = readRDS("<DISEASE_AND_DEATH_RDS>")$death_data
blood = readRDS("<BLOOD_COLLECTION_RDS>")
cov = readRDS("<COVARIATES_RDS>")
age = readRDS("<CHRONOLOGICAL_AGE_RDS>")
com = intersect(rownames(cov), rownames(age))
cov = cbind(age = age[com, ]$label, cov[com, ])
cov = cov[, c("age", "sex", "assessment_centre", "ethnicity_5cat", "town_index", "smoking", "alcohol", "physical_activity",
  "bmi")]
x <- as.character(cov$assessment_centre)
cov$assessment_centre <- factor(fcase(x %in% c("11004", "11005"), "Scotland", x %in% c("11003", "11022", "11023"),
  "Wales", x %in% c("11009", "11017", "11027", "11010", "11014"), "NE/Yorkshire", x %in% c("11008", "11001",
  "11016", "10003", "11024", "11025"), "North West", x %in% c("11021", "11013", "11006"), "Midlands", x %in%
  c("11011", "11028", "11002", "11007", "11026"), "South", x %in% c("11012", "11020", "11018"), "London"))

# 5.2 Clinical model settings and follow-up dates ----

if (getRversion() < "4.1.0") stop("R >= 4.1 is required.")
library(data.table)
library(survival)
library(parallel)
if (!requireNamespace("MASS", quietly = TRUE)) stop("Install the MASS package before running this script.")

# Censoring dates are retained from the supplied analysis. Use dates appropriate
# to the authorized data release if its administrative follow-up ends earlier.
MIN_CASE <- 100L
MIN_CONTROL <- 100L
N_WORKERS <- 5L
COV_ROW_ALIGNED_WITH_BLOOD <- TRUE
HES_END <- setNames(as.IDate(c("2023-03-31", "2022-08-31", "2022-05-31")), c("England", "Scotland", "Wales"))
DEATH_END <- setNames(as.IDate(c("2024-08-31", "2024-11-30", "2024-08-31")), c("England", "Scotland", "Wales"))
KEEP_ICD_CHAPTERS <- c(LETTERS[1:17], "U")
MOD <- c(hematologic = "immune_hematology", metabolomic = "immune_metabolomics", proteomic = "immune_proteomics")
AA <- c("AA_integrated", "AA_hematologic", "AA_metabolomic", "AA_proteomic")
COVARS <- c("age", "sex", "assessment_centre", "ethnicity_5cat", "town_index", "smoking", "alcohol", "physical_activity",
  "bmi")

# 5.3 Four age accelerations from full-cohort nested-CV OOF predictions ----

# Pool predictions across all folds, regress predicted age on chronological age,
# and standardize the residuals across the pooled sample.
stopifnot(inherits(dis, "data.frame"), inherits(blood, "data.frame"), inherits(cov, "data.frame"), inherits(death,
  "data.frame"), is.list(result))
o <- copy(as.data.table(result$full_cv$oof))
need <- c("eid", "label", "immune_age", unname(MOD))
if (!all(need %chin% names(o))) stop("Missing columns in result$full_cv$oof: ", paste(setdiff(need, names(o)),
  collapse = ", "))
z <- function(x) {
  s <- sd(x)
  if (any(!is.finite(x)) || !is.finite(s) || s == 0) stop("Age acceleration cannot be standardized.")
  as.numeric((x - mean(x)) / s)
}
rz <- function(p, a) z(lm.fit(cbind(1, a), p)$residuals)
aa <- o[, .(eid = as.integer(eid), clock_age = label, AA_integrated = rz(immune_age, label))]
aa[, (AA[ - 1L]) := lapply(unname(MOD), \(v) rz(o[[v]], o$label))]
if (anyNA(aa) || anyDuplicated(aa$eid)) stop("OOF age data contain invalid values or duplicate participant IDs.")

# 5.4 Merge blood collection dates, covariates, and death records ----

# Restore factor coding only for covariates that were originally factors.
fac <- intersect(COVARS, names(cov)[vapply(cov, is.factor, logical(1))])
lev <- setNames(lapply(fac, \(v) levels(cov[[v]])), fac)
restore_fac <- function(d, drop = FALSE) {
  for (v in intersect(fac, names(d))) set(d, j = v, value = factor(as.character(d[[v]]), levels = lev[[v]]))
  if (drop && length(fac)) d[, (fac) := lapply(.SD, droplevels), .SDcols = fac]
  d[]
}
bl <- as.data.table(blood)[, .(eid = as.integer(eid), blood_date = as.IDate(blood_date))]
cv <- copy(as.data.table(cov))
if (!"eid" %chin% names(cv)) {
  if (!COV_ROW_ALIGNED_WITH_BLOOD || nrow(cv) != nrow(bl)) stop("cov has no eid column, and row alignment with blood cannot be confirmed.")
  cv[, eid := bl$eid]
}
if (!all(COVARS %chin% names(cv))) stop("Missing columns in cov: ", paste(setdiff(COVARS, names(cv)), collapse =
  ", "))
cv[, eid := as.integer(eid)]
cv <- cv[, c("eid", COVARS), with = FALSE]
dd0 <- copy(as.data.table(death))
if (!all(c("eid", "death_date") %chin% names(dd0))) stop("death must contain eid and death_date.")
dc <- grep("^40000-", names(dd0), value = TRUE)
if (length(dc)) {
  chk <- do.call(pmin, c(lapply(dd0[, ..dc], \(x) as.numeric(as.IDate(x))), list(na.rm = TRUE)))
  chk[!is.finite(chk)] <- NA_real_
  giv <- as.numeric(as.IDate(dd0$death_date))
  if (any(xor(is.na(chk), is.na(giv)) | (!is.na(chk) & chk != giv))) stop("death_date is not the earliest date across instances of UKB Field 40000.")
}
dd <- dd0[, .(eid = as.integer(eid), death_date = as.IDate(death_date))]
if (anyNA(bl$eid) || anyNA(bl$blood_date) || anyNA(cv$eid) || anyNA(dd$eid) || anyDuplicated(bl$eid) || anyDuplicated(cv$eid) ||
  anyDuplicated(dd$eid)) stop("blood, cov, or death contain invalid or duplicate participant IDs, or invalid blood collection dates.")
base <- Reduce(\(x, y) merge(x, y, by = "eid", sort = FALSE), list(aa, bl, cv))
base <- merge(base, dd, by = "eid", all.x = TRUE, sort = FALSE)
restore_fac(base)
if (nrow(base) != nrow(aa) || anyNA(base[, c("eid", "clock_age", "blood_date", AA, COVARS), with = FALSE])) stop("Merging OOF data, blood, cov, and death failed, or required variables are missing.")
sl <- unname(lm.fit(cbind(1, base$age), base$clock_age)$coefficients[2])
age_audit <- base[, .(N = .N, r = cor(clock_age, age), slope = sl, mean_diff = mean(clock_age - age), median_abs_diff =
  median(abs(clock_age - age)), Q95 = unname(quantile(abs(clock_age - age), .95)), MAX = max(abs(clock_age - age)))]
print(age_audit)
if (!is.finite(age_audit$r) || age_audit$r < .99 || !is.finite(age_audit$slope) || age_audit$slope < .9 ||
  age_audit$slope > 1.1) stop("cov and blood are unlikely to be aligned by row; add verified participant IDs to cov before merging.")
if (age_audit$Q95 > 2) stop("OOF label and cov$age have inconsistent age reference times; check their sources.")
base[, age_cov_original := age][, age := clock_age][, clock_age := NULL]
ct <- tolower(trimws(as.character(base$assessment_centre)))
if (anyNA(ct) || any(!nzchar(ct))) stop("assessment_centre contains missing or empty values.")
base[, region := fcase(grepl("scot|edin|glasg", ct), "Scotland", grepl("wales|welsh|cardiff|swansea|wrexham",
  ct), "Wales", default = "England")]
base[, `:=`(disease_admin_end = pmin(HES_END[region], DEATH_END[region]), death_admin_end = DEATH_END[region])]
if (anyNA(base$disease_admin_end) || anyNA(base$death_admin_end) || any(base$blood_date >= base$disease_admin_end) ||
  any(base$blood_date >= base$death_admin_end) || any(base$death_date <= base$blood_date, na.rm = TRUE)) stop("Invalid regional censoring dates, death dates, or blood collection dates.")
base[, `:=`(disease_censor_date = disease_admin_end, mortality_event = as.integer(!is.na(death_date) & death_date <= death_admin_end),
  mortality_end_date = death_admin_end)][!is.na(death_date) & death_date < disease_admin_end, disease_censor_date := death_date][mortality_event == 1L,
  mortality_end_date := death_date]
base[, `:=`(disease_followup_years = as.numeric(disease_censor_date - blood_date) / 365.25, mortality_followup_years =
  as.numeric(mortality_end_date - blood_date) / 365.25)]
if (any(!is.finite(base$disease_followup_years) | base$disease_followup_years <= 0) || any(!is.finite(base$mortality_followup_years) |
  base$mortality_followup_years <= 0)) stop("Nonpositive follow-up times were detected.")

# 5.5 First recorded date for each three-character ICD-10 code ----

# Diagnosis on the blood collection date is prevalent. A disease event on the
# date of death is counted as a disease event.
d0 <- as.data.table(dis)
if (!all(c("eid", "icd3", "event_date") %chin% names(d0))) stop("dis must contain eid, icd3, and event_date.")
dx <- d0[eid %in% base$eid & !is.na(icd3) & !is.na(event_date), .(eid = as.integer(eid), icd3 = toupper(trimws(icd3)),
  event_date = as.IDate(event_date))]
dx <- dx[grepl("^[A-Z][0-9]{2}$", icd3) & substr(icd3, 1, 1) %chin% KEEP_ICD_CHAPTERS][base[, .(eid, blood_date,
  disease_admin_end, disease_censor_date)], on = "eid", nomatch = 0][event_date <= disease_admin_end]
first <- dx[, .(first_date = min(event_date), blood_date = blood_date[1L], disease_censor_date = disease_censor_date[1L]),
  by = .(eid, icd3)]
prev <- first[first_date <= blood_date, .(eid, icd3, first_date)]
inc <- first[first_date > blood_date & first_date <= disease_censor_date, .(eid, icd3, first_date)]
pc <- prev[, .(cases = uniqueN(eid)), by = icd3][, `:=`(risk_n = nrow(base), controls = nrow(base) - cases)][cases >= MIN_CASE &
  controls >= MIN_CONTROL]
px <- prev[, .(prevalent_excluded = uniqueN(eid)), by = icd3]
ic <- merge(inc[, .(events = uniqueN(eid)), by = icd3], px, by = "icd3", all.x = TRUE, sort = FALSE)
ic[is.na(prevalent_excluded), prevalent_excluded := 0L][, `:=`(risk_n = nrow(base) - prevalent_excluded, censored =
  nrow(base) - prevalent_excluded - events)]
ic <- ic[events >= MIN_CASE & censored >= MIN_CONTROL]
setorder(pc, - cases)
setorder(ic, - events)
if (!nrow(pc) && !nrow(ic)) stop("No three-character ICD-10 endpoint meets the case and noncase count thresholds.")

# 5.6 Participant-level disease burden ----

# Summarize baseline prevalence, incident diseases, and observed follow-up.
ps <- prev[, .(n_prevalent = uniqueN(icd3)), by = eid]
ins <- inc[, .(n_incident = uniqueN(icd3), first_incident = min(first_date)), by = eid]
base <- merge(merge(base, ps, by = "eid", all.x = TRUE, sort = FALSE), ins, by = "eid", all.x = TRUE, sort =
  FALSE)
restore_fac(base)
base[is.na(n_prevalent), n_prevalent := 0L][is.na(n_incident), n_incident := 0L][, `:=`(any_prevalent = as.integer(n_prevalent > 0L),
  any_incident = as.integer(n_incident > 0L), n_total_observed = n_prevalent + n_incident)]

# 5.7 Disease associations with each age acceleration modeled separately ----

# Prevalent disease: logistic regression. Incident disease: cause-specific Cox.
prep <- function(d) {
  d <- copy(d)
  restore_fac(d, TRUE)
  d
}
safe <- function(expr) {
  w <- character()
  f <- withCallingHandlers(tryCatch(expr, error = identity), warning = \(q) {
    w <<- c(w, conditionMessage(q))
    invokeRestart("muffleWarning")
  })
  list(f = f, w = paste(unique(w), collapse = " | "))
}
adjvars <- function(d) COVARS[vapply(d[, ..COVARS], \(v) uniqueN(v) > 1L, logical(1))]
fail <- function(N, events = NA_integer_, w = "") data.table(N = N, events = events, beta = NA_real_, SE =
  NA_real_, effect = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_, PH_p = NA_real_, PH_global_p =
  NA_real_, concordance = NA_real_, converged = FALSE, warning = w)
fit_glm <- function(d, x, type = c("binary", "count"), offset = FALSE) {
  type <- match.arg(type)
  rhs <- c(x, adjvars(d), if (offset) "offset(log(disease_followup_years))")
  fm <- as.formula(paste("outcome~", paste(rhs, collapse = "+")))
  q <- safe(if (type == "binary") glm(fm, data = d, family = binomial(), control = glm.control(maxit = 100)) else MASS::glm.nb(fm,
    data = d, control = glm.control(maxit = 100)))
  f <- q$f
  if (inherits(f, "error") || !x %chin% rownames(coef(summary(f)))) return(fail(nrow(d), if (type == "binary") sum(d$outcome == 1L) else sum(d$outcome),
    paste(q$w, if (inherits(f, "error")) conditionMessage(f) else "")))
  a <- coef(summary(f))[x, ]
  b <- unname(a[1])
  s <- unname(a[2])
  pv <- if ("Pr(>|z|)" %in% names(a)) unname(a["Pr(>|z|)"]) else unname(a[4])
  ef <- exp(b)
  lo <- exp(b - 1.96 * s)
  hi <- exp(b + 1.96 * s)
  ok <- isTRUE(f$converged) && all(is.finite(c(b, s, ef, lo, hi, pv)))
  if (!ok) return(fail(nrow(d), if (type == "binary") sum(d$outcome == 1L) else sum(d$outcome), q$w))
  data.table(N = nobs(f), events = if (type == "binary") sum(d$outcome == 1L) else sum(d$outcome), beta =
    b, SE = s, effect = ef, LCL = lo, UCL = hi, p = pv, PH_p = NA_real_, PH_global_p = NA_real_, concordance =
    NA_real_, converged = TRUE, warning = q$w)
}
fit_cox <- function(d, x) {
  fm <- as.formula(paste("Surv(time,event)~", paste(c(x, adjvars(d)), collapse = "+")))
  q <- safe(coxph(fm, data = d, ties = "efron", x = TRUE, model = FALSE, y = TRUE))
  f <- q$f
  if (inherits(f, "error") || !x %chin% rownames(coef(summary(f)))) return(fail(nrow(d), sum(d$event), paste(q$w,
    if (inherits(f, "error")) conditionMessage(f) else "")))
  sm <- summary(f)
  a <- sm$coefficients[x, ]
  ci <- sm$conf.int[x, ]
  ok <- all(is.finite(c(a[c("coef", "se(coef)", "exp(coef)", "Pr(>|z|)")], ci[c("lower .95", "upper .95")])))
  if (!ok) return(fail(nrow(d), sum(d$event), q$w))
  zp <- tryCatch(cox.zph(f, transform = "km")$table, error = \(e) NULL)
  ph <- if (is.null(zp) || !x %in% rownames(zp)) NA_real_ else zp[x, "p"]
  pg <- if (is.null(zp) || !"GLOBAL" %in% rownames(zp)) NA_real_ else zp["GLOBAL", "p"]
  data.table(N = f$n, events = f$nevent, beta = unname(a["coef"]), SE = unname(a["se(coef)"]), effect = unname(a["exp(coef)"]),
    LCL = unname(ci["lower .95"]), UCL = unname(ci["upper .95"]), p = unname(a["Pr(>|z|)"]), PH_p = ph, PH_global_p =
    pg, concordance = unname(sm$concordance[1]), converged = TRUE, warning = q$w)
}
pids <- split(prev$eid, prev$icd3)
iev <- split(inc[, .(eid, first_date)], inc$icd3)
ids <- function(L, k) if (is.null(L[[k]])) integer() else L[[k]]
scan_prev <- function(code) {
  d <- prep(base)
  d[, outcome := as.integer(eid %in% ids(pids, code))]
  rbindlist(lapply(AA, \(x) cbind(data.table(analysis = "prevalent", icd3 = code, exposure = x, measure =
    "OR"), fit_glm(d, x, "binary"))))
}
scan_inc <- function(code) {
  d <- prep(base[!eid %in% ids(pids, code)])
  ev <- iev[[code]]
  if (is.null(ev) || !nrow(ev) || anyDuplicated(ev$eid)) stop("Invalid incident event table: ", code)
  d[, `:=`(event = 0L, end_date = disease_censor_date)]
  j <- match(ev$eid, d$eid)
  if (anyNA(j)) stop("Participant ID matching failed for incident events: ", code)
  set(d, j, "event", 1L)
  set(d, j, "end_date", ev$first_date)
  d[, time := as.numeric(end_date - blood_date) / 365.25]
  if (any(!is.finite(d$time) | d$time <= 0) || sum(d$event) != nrow(ev)) stop("Cox risk-set construction failed: ",
    code)
  nd <- sum(d$event == 0L & !is.na(d$death_date) & d$death_date <= d$disease_admin_end)
  tie <- sum(!is.na(d$death_date[j]) & ev$first_date == d$death_date[j])
  py <- sum(d$time)
  rbindlist(lapply(AA, \(x) cbind(data.table(analysis = "incident_cause_specific_cox", icd3 = code, exposure =
    x, measure = "HR", person_years = py, deaths_before_event = nd, same_day_disease_death = tie), fit_cox(d,
    x))))
}
jobs <- rbindlist(list(data.table(kind = "prevalent", code = pc$icd3), data.table(kind = "incident", code =
  ic$icd3)))
stopifnot(nrow(jobs) > 0L, N_WORKERS >= 1L)
ncpu <- detectCores(logical = FALSE)
if (is.na(ncpu)) ncpu <- N_WORKERS
nw <- as.integer(max(1L, min(N_WORKERS, ncpu, nrow(jobs))))
message("Parallel disease scans: ", nrow(jobs), " endpoints × 4 clocks, ", nw, " workers")
run_job <- function(i) {
  old <- getDTthreads()
  on.exit(setDTthreads(old), add = TRUE)
  setDTthreads(1L)
  tryCatch(list(ok = TRUE, value = if (jobs$kind[i] == "prevalent") scan_prev(jobs$code[i]) else scan_inc(jobs$code[i])),
    error = \(e) list(ok = FALSE, error = conditionMessage(e)))
}
ans <- mclapply(seq_len(nrow(jobs)), run_job, mc.cores = nw, mc.preschedule = FALSE, mc.set.seed = FALSE)
ok <- vapply(ans, \(x) is.list(x) && isTRUE(x$ok) && inherits(x$value, "data.frame") && nrow(x$value) == length(AA),
  logical(1))
if (any(!ok)) {
  msg <- vapply(which(!ok), \(i) {
    x <- ans[[i]]
    e <- if (inherits(x, "try-error")) as.character(x) else if (is.list(x) && !is.null(x$error)) x$error else "worker crashed/invalid result"
    paste0(jobs$kind[i], " ", jobs$code[i], ": ", e)
  }, character(1))
  stop("Parallel disease models failed:\n", paste(msg, collapse = "\n"), call. = FALSE)
}
assoc_disease <- rbindlist(lapply(ans, `[[`, "value"), fill = TRUE)
assoc_disease[, FDR_within_clock := p.adjust(p, "BH"), by = .(analysis, exposure)][, FDR_all_clocks := p.adjust(p,
  "BH"), by = analysis][, PH_flag := analysis == "incident_cause_specific_cox" & !is.na(PH_p) & PH_p < .01]
setorder(assoc_disease, analysis, exposure, FDR_within_clock, p)

# 5.8 Disease presence and counts ----

# Analyze baseline burden separately from incident burden during follow-up.
burden <- list()
d <- prep(base)
d[, outcome := any_prevalent]
burden[[1]] <- rbindlist(lapply(AA, \(x) cbind(data.table(analysis = "any_prevalent", exposure = x, measure =
  "OR"), fit_glm(d, x, "binary"))))
d <- prep(base)
d[, outcome := n_prevalent]
burden[[2]] <- rbindlist(lapply(AA, \(x) cbind(data.table(analysis = "number_prevalent_ICD3", exposure = x,
  measure = "IRR"), fit_glm(d, x, "count"))))
d <- prep(base)
d[, `:=`(event = as.integer(!is.na(first_incident)), end_date = disease_censor_date)][event == 1L, end_date := first_incident][,
  time := as.numeric(end_date - blood_date) / 365.25]
if (any(!is.finite(d$time) | d$time <= 0)) stop("Invalid follow-up time for the first incident disease analysis.")
burden[[3]] <- rbindlist(lapply(AA, \(x) cbind(data.table(analysis = "time_to_first_new_ICD3", exposure =
  x, measure = "HR"), fit_cox(d, x))))
d <- prep(base)
d[, outcome := n_incident]
burden[[4]] <- rbindlist(lapply(AA, \(x) cbind(data.table(analysis = "rate_of_new_ICD3", exposure = x, measure =
  "IRR"), fit_glm(d, x, "count", TRUE))))
assoc_burden <- rbindlist(burden, fill = TRUE)
assoc_burden[, FDR_within_clock := p.adjust(p, "BH"), by = exposure][, FDR_all := p.adjust(p, "BH")][, PH_flag := grepl("time_to_",
  analysis) & !is.na(PH_p) & PH_p < .01]

# 5.9 All-cause mortality ----

# Use death-registry follow-up and fit a separate Cox model for each standardized
# age acceleration.
dm <- prep(base)
dm[, `:=`(event = mortality_event, time = mortality_followup_years)]
if (sum(dm$event) < MIN_CASE || sum(dm$event == 0L) < MIN_CONTROL) stop("Insufficient events or nonevents for the all-cause mortality analysis.")
sf <- survfit(Surv(time, 1L - event) ~ 1, data = dm)
st <- summary(sf)$table
rkmed <- unname(if (is.matrix(st)) st[1L, "median"] else st["median"])
mortality_summary <- dm[, .(N = .N, deaths = sum(event), alive_or_censored = sum(event == 0L), post_censor_deaths =
  sum(!is.na(death_date) & death_date > death_admin_end), person_years = sum(time), median_observed_followup_years =
  median(time), reverse_KM_median_followup_years = rkmed)]
mortality_by_region <- dm[, .(N = .N, deaths = sum(event), alive_or_censored = sum(event == 0L), post_censor_deaths =
  sum(!is.na(death_date) & death_date > death_admin_end), person_years = sum(time)), by = .(region, death_admin_end)]
assoc_mortality <- rbindlist(lapply(AA, \(x) cbind(data.table(analysis = "all_cause_mortality", exposure =
  x, measure = "HR", person_years = mortality_summary$person_years, reverse_KM_median_followup_years = mortality_summary$reverse_KM_median_followup_years),
  fit_cox(dm, x))))
assoc_mortality[, `:=`(FDR_BH = p.adjust(p, "BH"), PH_flag = !is.na(PH_p) & PH_p < .01)]
assoc_mortality <- assoc_mortality[, .(analysis, exposure, measure, N, deaths = events, person_years, reverse_KM_median_followup_years,
  beta, SE, HR = effect, HR_LCL = LCL, HR_UCL = UCL, P = p, FDR_BH, PH_P = PH_p, PH_global_P = PH_global_p,
  PH_flag, concordance, converged, warning)]
setorder(assoc_mortality, FDR_BH, P)

# 5.10 Save disease and mortality results ----

# Retain source-specific censoring dates and use long-form disease event tables.
disease_result <- list(definition = list(time_zero = "blood_date", dis_source = "inpatient ICD-10 fields 41270/41280",
  AA_scale = "HR/OR/IRR per 1-SD globally age-adjusted OOF acceleration", incident_model = "cause-specific Cox; death censored",
  mortality_model = "all-cause Cox", warning_policy = "warnings retained; target AA estimate accepted only when beta, SE, effect, CI and P are finite",
  parallel = list(method = "parallel::mclapply", requested_workers = N_WORKERS, used_workers = nw, mc.preschedule =
  FALSE), same_day_rule = "disease/death event takes priority over censoring", loss_to_follow_up = "Field 191 not supplied; no loss-date censoring applied",
  min_cases = MIN_CASE, min_controls = MIN_CONTROL, HES_end = HES_END, death_registry_end = DEATH_END, censoring_source =
  "https://biobank.ndph.ox.ac.uk/ukb/exinfo.cgi?src=Data_providers_and_dates", dates_checked = "2026-09-01",
  ICD_chapters = KEEP_ICD_CHAPTERS), age_alignment = age_audit, centre_region_map = unique(base[, .(assessment_centre,
  region, disease_admin_end, death_admin_end)]), participant_data = base, event_tables = list(prevalent =
  prev, incident_observed = inc), case_counts = list(prevalent = pc, incident = ic), mortality_summary = list(overall =
  mortality_summary, by_region = mortality_by_region), association = list(disease_specific = assoc_disease,
  burden = assoc_burden, all_cause_mortality = assoc_mortality), session = sessionInfo())
saveRDS(disease_result, "immune_age_disease_cox_result.rds")
fwrite(assoc_disease, "disease_specific_associations_cox.csv")
fwrite(assoc_burden, "disease_burden_associations.csv")
fwrite(assoc_mortality, "all_cause_mortality_associations.csv")
fwrite(mortality_summary, "all_cause_mortality_followup_summary.csv")
fwrite(mortality_by_region, "all_cause_mortality_followup_by_region.csv")
fwrite(pc, "prevalent_case_counts.csv")
fwrite(ic, "incident_event_counts.csv")
print(base[, .(N = .N, deaths_for_mortality = sum(mortality_event), any_prevalent = sum(any_prevalent), any_incident =
  sum(any_incident), median_disease_followup = median(disease_followup_years), median_mortality_followup =
  median(mortality_followup_years))])
print(pc)
print(ic)
print(assoc_burden)
print(mortality_summary)
print(assoc_mortality)

# 6. Incident-disease multiplicity and cross-clock overlap ----

ALPHA <- .05
aa <- c("AA_integrated", "AA_hematologic", "AA_metabolomic", "AA_proteomic")
id <- c("Integrated", "Hematologic", "Metabolomic", "Proteomic")
nm <- c("Integrated multimodal", "Hematologic", "Metabolomic", "Proteomic")
full <- paste0(nm, " immune age")
z <- copy(as.data.table(disease_result$association$disease_specific))
need <- c("analysis", "icd3", "exposure", "measure", "events", "beta", "SE", "effect", "LCL", "UCL", "p",
  "PH_p", "converged", "warning")
if (length(miss <- setdiff(need, names(z)))) stop("Results are missing these variables: ", paste(miss, collapse =
  ", "))
z <- z[analysis == "incident_cause_specific_cox" & measure == "HR" & exposure %chin% aa]
if (!nrow(z)) stop("No incident cause-specific Cox results were found")
z[, icd3 := toupper(trimws(as.character(icd3)))]
if (anyDuplicated(z[, .(icd3, exposure)])) stop("Duplicate disease-clock result rows were found")
if (nrow(x <- z[, .N, icd3][N != 4L])) stop("These diseases do not have four expected model rows: ", paste(x$icd3,
  collapse = ", "))
planned <- NULL
if (!is.null(disease_result$case_counts$incident)) {
  pc <- as.data.table(disease_result$case_counts$incident)
  if (!"icd3" %chin% names(pc)) stop("case_counts$incident is missing icd3")
  planned <- unique(toupper(trimws(as.character(pc$icd3))))
  planned <- planned[grepl("^[A-Z][0-9]{2}$", planned)]
  if (!setequal(unique(z$icd3), planned)) stop("Incident Cox results do not match the prespecified case_counts$incident disease universe; the Bonferroni denominator must not be reduced")
}
M0 <- if (length(planned)) length(planned) else uniqueN(z$icd3)
if (!length(planned)) warning("case_counts$incident is unavailable; the Bonferroni denominator is derived from the complete incident Cox results")
BONF_FAMILY <- 4L * M0
BONF_THRESHOLD <- ALPHA / BONF_FAMILY
z[, cid := id[match(exposure, aa)]]
z[, valid := converged %in% TRUE & is.finite(beta) & is.finite(SE) & SE > 0 & is.finite(effect) & effect > 0 &
  is.finite(LCL) & LCL > 0 & is.finite(UCL) & UCL > 0 & is.finite(p) & between(p, 0, 1) & LCL <= effect &
  effect <= UCL]
z[, P_bonferroni := NA_real_]
z[valid == TRUE, P_bonferroni := pmin(1, p * BONF_FAMILY)]
z[, failure_reason := fcase(valid == TRUE, NA_character_, !(converged %in% TRUE), "Model did not converge",
  !is.finite(beta) | !is.finite(SE) | !is.finite(effect) | !is.finite(LCL) | !is.finite(UCL), "Non-finite estimate or CI",
  !is.finite(p), "Missing P value", default = "Model unavailable")]
endpoint_audit <- z[, .(valid_clocks = sum(valid == TRUE)), icd3]
excluded <- endpoint_audit[valid_clocks < 4L]
U <- endpoint_audit[valid_clocks == 4L, icd3]
if (!length(U)) stop("No incident diseases have valid estimates for all four clocks")
if (nrow(excluded)) warning(nrow(excluded), " diseases were excluded from cross-clock comparisons because at least one clock lacked a valid estimate; see excluded")
ch <- c("Infectious diseases", "Neoplasms", "Blood and immune disorders", "Endocrine and metabolic", "Mental and behavioural",
  "Nervous system", "Eye", "Ear", "Circulatory", "Respiratory", "Digestive", "Skin", "Musculoskeletal", "Genitourinary",
  "Pregnancy", "Perinatal", "Congenital", "Special purposes")
z[, `:=`(letter = substr(icd3, 1, 1), number = suppressWarnings(as.integer(substr(icd3, 2, 3))))]
z[, chapter := fcase(letter %chin% c("A", "B"), ch[1], letter == "C" | (letter == "D" & number <= 48), ch[2],
  letter == "D" & number >= 50, ch[3], letter == "E", ch[4], letter == "F", ch[5], letter == "G", ch[6], letter == "H" &
  number <= 59, ch[7], letter == "H" & number >= 60, ch[8], letter == "I", ch[9], letter == "J", ch[10], letter == "K",
  ch[11], letter == "L", ch[12], letter == "M", ch[13], letter == "N", ch[14], letter == "O", ch[15], letter == "P",
  ch[16], letter == "Q", ch[17], letter == "U", ch[18], default = NA_character_)]
if (z[is.na(chapter), .N]) stop("Unable to map ICD3 chapters for: ", paste(unique(z[is.na(chapter), icd3]),
  collapse = ", "))
z[, sig := valid == TRUE & is.finite(P_bonferroni) & P_bonferroni < ALPHA]
d <- copy(z[icd3 %chin% U])

# 6.1 Recovery of single-clock significant disease sets ----

# Restrict comparisons to endpoints with valid estimates for all four clocks.
ist <- d[cid == "Integrated", .(icd3, integrated_sig = sig)]
ss <- merge(d[cid != "Integrated" & sig == TRUE, .(icd3, cid)], ist, "icd3")
st <- c("Recovered by integrated clock", "Not recovered")
ss[, status := fifelse(integrated_sig == TRUE, st[1], st[2])]
den <- merge(data.table(cid = id[ - 1]), ss[, .(den = .N), cid], "cid", all.x = TRUE)[is.na(den), den := 0L]
cap <- merge(CJ(cid = id[ - 1], status = st, unique = TRUE), ss[, .N, .(cid, status)], c("cid", "status"),
  all.x = TRUE)
cap <- merge(cap, den, "cid", all.x = TRUE)[is.na(N), N := 0L]
cap[, `:=`(status = factor(status, st), pct = fifelse(den > 0, N / den, NA_real_))]
setorder(cap, cid, status)

# 6.2 Shared and integrated-only disease sets ----

mem <- d[, .(all_four = all(sig == TRUE), integrated_only = sig[cid == "Integrated"] & all(!sig[cid != "Integrated"])),
  .(icd3, chapter)]
all4 <- mem[all_four == TRUE, icd3]
only <- mem[integrated_only == TRUE, icd3]
if (!length(all4)) stop("No diseases are globally Bonferroni-significant for all four clocks")
if (!length(only)) stop("No diseases are globally Bonferroni-significant only for the integrated clock")

# 7. Mortality dose response, absolute risk, and prognostic performance ----

RESULT_RDS <- "immune_age_disease_cox_result.rds"
S <- 666L
B_BOOT <- 200L
NC <- 5L
STD_N <- 20000L
EVAL_TIMES <- c(5, 10, 15)
TAU <- 10
CACHE <- "mortality_bootstrap.rds"
pk <- c("data.table", "survival", "timeROC", "survAUC")
mis <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
if (length(mis)) stop("Install the required packages: ", paste(mis, collapse = ", "))
invisible(lapply(pk, library, character.only = TRUE))
if (!exists("disease_result", inherits = FALSE)) {
  if (!file.exists(RESULT_RDS)) stop("File not found: ", RESULT_RDS)
  disease_result <- readRDS(RESULT_RDS)
}
AA <- c("Integrated multimodal" = "AA_integrated", "Hematologic" = "AA_hematologic", "Metabolomic" = "AA_metabolomic",
  "Proteomic" = "AA_proteomic")
COV0 <- c("age", "sex", "assessment_centre", "ethnicity_5cat", "town_index", "smoking", "alcohol", "physical_activity",
  "bmi")

# 7.1 Mortality follow-up and baseline Cox models ----

# Use participant data and OOF accelerations from module 5.
D <- copy(as.data.table(disease_result$participant_data))
need <- c("mortality_followup_years", "mortality_event", unname(AA), COV0)
if (!all(need %chin% names(D))) stop("participant_data is missing: ", paste(setdiff(need, names(D)), collapse =
  ", "))
D <- D[, ..need]
setnames(D, c("mortality_followup_years", "mortality_event"), c("time", "event"))
if (anyNA(D) || any(!is.finite(D$time)) || any(D$time <= 0) || !all(D$event %in% c(0L, 1L))) stop("Mortality data contain missing values or invalid follow-up times or outcomes")
D <- as.data.table(droplevels(as.data.frame(D)))
COV <- COV0[vapply(D[, ..COV0], uniqueN, integer(1)) > 1L]
if (any(EVAL_TIMES <= 0) || anyDuplicated(EVAL_TIMES) || is.unsorted(EVAL_TIMES) || !TAU %in% EVAL_TIMES) stop("EVAL_TIMES must be increasing and include TAU")
Ghat <- \(d, t) summary(survfit(Surv(time, 1 - event) ~ 1, data = d), times = t, extend = TRUE)$surv
badT <- EVAL_TIMES[EVAL_TIMES >= max(D$time) | vapply(EVAL_TIMES, \(t) sum(D$event == 1L & D$time <= t) < 50L ||
  sum(D$time > t) < 100L, logical(1)) | Ghat(D, EVAL_TIMES) < .05]
if (length(badT)) stop("These horizons lack sufficient events, participants at risk, or IPCW support; revise EVAL_TIMES: ",
  paste(badT, collapse = ", "))
RISK_TIMES <- seq(0, max(EVAL_TIMES), length.out = 121)
dc <- parallel::detectCores(logical = FALSE)
if (is.na(dc)) dc <- NC
NC <- max(1L, min(as.integer(NC), dc, B_BOOT))
if (B_BOOT < 50L) stop("B_BOOT must be at least 50")
set.seed(S)
STD_I <- sample.int(nrow(D), min(as.integer(STD_N), nrow(D)))
STD <- D[STD_I]
rhs <- paste(COV, collapse = "+")
terms <- c("Covariates" = "", AA)
FML <- lapply(terms, \(x) as.formula(paste0("Surv(time,event)~", paste(c(if (nzchar(x)) x, COV), collapse =
  "+"))))
FIT <- lapply(FML, \(f) coxph(f, D, ties = "efron", x = TRUE, y = TRUE))
pval <- \(p) format.pval(p, digits = 2, eps = 1e-4)
lrt <- \(big, small) {
  a <- logLik(big)
  b <- logLik(small)
  pchisq(2 * (as.numeric(a) - as.numeric(b)), attr(a, "df") - attr(b, "df"), lower.tail = FALSE)
}

# 7.2 Linear Cox associations and four-clock Holm correction ----

M <- copy(as.data.table(disease_result$association$all_cause_mortality))
req <- c("exposure", "HR", "HR_LCL", "HR_UCL", "P", "PH_P", "converged")
if (!all(req %chin% names(M))) stop("all_cause_mortality is missing: ", paste(setdiff(req, names(M)), collapse =
  ", "))
M <- M[exposure %chin% unname(AA)]
if (nrow(M) != 4L || anyDuplicated(M$exposure) || any(is.na(M$converged) | M$converged != TRUE) || any(!is.finite(M[,
  unlist(.SD), .SDcols = c("HR", "HR_LCL", "HR_UCL", "P")]))) stop("The four mortality Cox results are incomplete")
M[, `:=`(clock = names(AA)[match(exposure, AA)], P_Holm = p.adjust(P, "holm"))][, clock := factor(clock, levels =
  rev(names(AA)))]
if (any(M$PH_P < .05, na.rm = TRUE)) warning("At least one age-acceleration PH P value is below 0.05; inspect Schoenfeld residuals and consider time-varying hazard ratios before interpretation")

# 7.3 Spline dose responses and likelihood-ratio tests ----

# Knots: 5th/35th/65th/95th percentiles; reference = 0 SD.
rcs_formula <- function(v, k) as.formula(sprintf("Surv(time,event)~splines::ns(%s,knots=c(%.17g,%.17g),Boundary.knots=c(%.17g,%.17g))+%s",
  v, k[2], k[3], k[1], k[4], rhs))
mm <- function(f, d) {
  x <- model.matrix(f, data = d)
  x[, names(coef(f)), drop = FALSE]
}
rcs_one <- function(nm, v) {
  k <- unname(quantile(D[[v]], c(.05, .35, .65, .95)))
  fs <- coxph(rcs_formula(v, k), D, ties = "efron", x = TRUE, y = TRUE)
  g <- seq(quantile(D[[v]], .01), quantile(D[[v]], .99), length.out = 181)
  nd <- D[rep(1L, length(g))]
  set(nd, j = v, value = g)
  n0 <- copy(nd)
  set(n0, j = v, value = 0)
  dx <- mm(fs, nd) - mm(fs, n0)
  eta <- drop(dx %*% coef(fs))
  se <- sqrt(pmax(0, rowSums((dx %*% vcov(fs)) * dx)))
  list(curve = data.table(clock = nm, x = g, HR = exp(eta), lo = exp(eta - 1.96 * se), hi = exp(eta + 1.96 * se)),
    ann = data.table(clock = nm, label = sprintf("P-overall %s\nP-nonlinear %s", pval(lrt(fs, FIT[["Covariates"]])),
    pval(lrt(fs, FIT[[nm]])))), fit = fs)
}
SP <- Map(rcs_one, names(AA), unname(AA))
SC <- rbindlist(lapply(SP, `[[`, "curve"))
SA <- rbindlist(lapply(SP, `[[`, "ann"))
SC[, clock := factor(clock, levels = names(AA))]
SA[, clock := factor(clock, levels = names(AA))]

# 7.4 Population-standardized cumulative mortality risk ----

std_risk <- function(f, pop, times) {
  bh <- basehaz(f, centered = FALSE)
  H <- approx(c(0, bh$time), c(0, bh$hazard), times, method = "constant", f = 0, rule = 2)$y
  if (any(!is.finite(H))) stop("Invalid baseline cumulative hazard")
  rbindlist(lapply(c( - 1, 0, 1), \(a) {
    nd <- copy(pop)
    set(nd, j = "AA_integrated", value = a)
    lp <- as.numeric(predict(f, newdata = nd, type = "lp", reference = "zero"))
    ee <- exp(lp)
    if (any(!is.finite(ee))) stop("Invalid standardized predictions")
    z <- vapply(H, \(h) mean( - expm1( - h * ee)), numeric(1))
    if (any(!is.finite(z))) stop("Invalid standardized absolute risk")
    data.table(scenario = factor(c("Lower (−1 SD)", "Average (0 SD)", "Higher (+1 SD)")[a + 2L], levels =
      c("Lower (−1 SD)", "Average (0 SD)", "Higher (+1 SD)")), time = times, risk = z)
  }))
}
KINT <- unname(quantile(D$AA_integrated, c(.05, .35, .65, .95)))
FINTS <- rcs_formula("AA_integrated", KINT)
C0 <- std_risk(coxph(FINTS, D, ties = "efron", x = TRUE, y = TRUE), D, RISK_TIMES)

# 7.5 Paired bootstrap and out-of-bag discrimination ----

boot_one <- function(b) {
  old <- getDTthreads()
  on.exit(setDTthreads(old), add = TRUE)
  setDTthreads(1L)
  set.seed(S + b)
  ib <- sample.int(nrow(D), nrow(D), replace = TRUE)
  ob <- setdiff(seq_len(nrow(D)), unique(ib))
  tr <- D[ib]
  te <- D[ob]
  if (any(vapply(EVAL_TIMES, \(t) sum(te$event == 1L & te$time <= t) < 20L || sum(te$time > t) < 50L, logical(1))) ||
    any(Ghat(te, EVAL_TIMES) < .05)) stop("Insufficient OOB events, participants at risk, or IPCW support")
  ff <- lapply(FML, \(f) coxph(f, tr, ties = "efron", x = TRUE, y = TRUE))
  lp <- lapply(ff, \(f) as.numeric(predict(f, newdata = te, type = "lp", reference = "zero")))
  if (any(!is.finite(unlist(lp)))) stop("Invalid OOB linear predictions")
  au <- rbindlist(Map(\(x, nm) {
    q <- timeROC::timeROC(T = te$time, delta = te$event, marker = x, cause = 1, weighting = "marginal", times =
      EVAL_TIMES, ROC = FALSE, iid = FALSE)
    j <- match(EVAL_TIMES, q$times)
    if (anyNA(j)) stop("timeROC did not return the requested horizons")
    data.table(model = nm, time = EVAL_TIMES, AUC = q$AUC[j])
  }, lp, names(lp)))
  uc <- data.table(model = names(lp), UnoC = vapply(lp, \(x) survAUC::UnoC(Surv(tr$time, tr$event), Surv(te$time,
    te$event), x, time = TAU), numeric(1)))
  fs <- coxph(FINTS, tr, ties = "efron", x = TRUE, y = TRUE)
  list(auc = au, uno = uc, risk = std_risk(fs, STD, RISK_TIMES))
}
valid_boot <- function(x) tryCatch({
  m <- unique(x$auc$model)
  is.list(x) && all(c("auc", "uno", "risk") %in% names(x)) && nrow(x$auc) == length(m) * length(EVAL_TIMES) && nrow(x$uno) == length(m) && nrow(x$risk) == 3L * length(RISK_TIMES) && !anyDuplicated(x$auc[,
    .(model, time)]) && all(names(FML) %in% m) && setequal(m, x$uno$model) && all(is.finite(x$auc$AUC)) && all(is.finite(x$uno$UnoC)) && all(is.finite(x$risk$risk))
}, error = \(e) FALSE)
fp <- c(sum(D$time), sum(D$time ^ 2), sum(D$time * D$event), vapply(unname(AA), \(v) sum(D[[v]]), numeric(1)),
  vapply(unname(AA), \(v) sum(D[[v]] ^ 2), numeric(1)))
key <- list(method = "paired_OOB_bootstrap_v3", N = nrow(D), deaths = sum(D$event), B = B_BOOT, times = EVAL_TIMES,
  risk_times = RISK_TIMES, tau = TAU, seed = S, std_n = nrow(STD), std_i = STD_I, knots = KINT, cov = COV,
  models = names(FML), fingerprint = fp)
BO <- NULL
FAIL <- character()
if (file.exists(CACHE)) {
  cc <- readRDS(CACHE)
  kk <- setdiff(names(key), "models")
  if (is.list(cc) && is.list(cc$key) && identical(cc$key[kk], key[kk]) && all(names(FML) %in% cc$key$models) &&
    length(cc$boot) && all(vapply(cc$boot, valid_boot, logical(1)))) {
    BO <- cc$boot
    FAIL <- cc$failures
  }
}
if (is.null(BO)) {
  message("Paired OOB bootstrap: ", B_BOOT, " resamples, ", NC, " Linux fork workers")
  zz <- parallel::mclapply(seq_len(B_BOOT), \(b) tryCatch(boot_one(b), error = \(e) list(error = conditionMessage(e))),
    mc.cores = NC, mc.preschedule = FALSE, mc.set.seed = FALSE)
  ok <- vapply(zz, valid_boot, logical(1))
  FAIL <- vapply(zz[!ok], \(x) if (is.list(x) && "error" %in% names(x)) x$error else "Invalid bootstrap output",
    character(1))
  if (sum(ok) < ceiling(.9 * B_BOOT)) stop("Fewer than 90% of bootstrap resamples are valid; first error: ",
    FAIL[1L])
  if (length(FAIL)) warning(length(FAIL), " bootstrap resamples failed; at least 90% succeeded, and errors were saved")
  BO <- zz[ok]
  saveRDS(list(key = key, boot = BO, failures = FAIL), CACHE)
}
message("Valid bootstrap resamples: ", length(BO), "/", B_BOOT)

# 7.6 Bootstrap confidence intervals for absolute risk ----

BR <- rbindlist(lapply(seq_along(BO), \(i) cbind(b = i, BO[[i]]$risk)))
BC <- BR[, .(lo = quantile(risk, .025), hi = quantile(risk, .975)), by = .(scenario, time)]
CC <- merge(C0, BC, by = c("scenario", "time"), sort = FALSE)

# 7.7 Time-dependent AUC gains over the covariate-only model ----

AU <- rbindlist(lapply(seq_along(BO), \(i) cbind(b = i, BO[[i]]$auc)))
AD <- merge(AU[model %chin% names(AA)], AU[model == "Covariates", .(b, time, base = AUC)], by = c("b", "time"))
AD[, gain := 100 * (AUC - base)]
DS <- AD[, .(gain = mean(gain), lo = quantile(gain, .025), hi = quantile(gain, .975)), by = .(model, time)][,
  model := factor(model, levels = names(AA))]

# 7.8 Discrimination and paired comparisons at the selected horizon ----

UC <- rbindlist(lapply(seq_along(BO), \(i) cbind(b = i, BO[[i]]$uno)))
EA <- AU[model %chin% names(AA) & time == TAU, .(value = mean(AUC), lo = quantile(AUC, .025), hi = quantile(AUC,
  .975)), by = model][, metric := sprintf("%g-year AUC", TAU)]
EU <- UC[model %chin% names(AA), .(value = mean(UnoC), lo = quantile(UnoC, .025), hi = quantile(UnoC, .975)),
  by = model][, metric := sprintf("Uno C (≤%g years)", TAU)]
EE <- rbind(EA, EU)[, model := factor(model, levels = rev(names(AA)))]
int <- names(AA)[1L]
sing <- names(AA)[ - 1L]
ed1 <- merge(AU[model == int & time == TAU, .(b, int = AUC)], AU[model %chin% sing & time == TAU], by = "b")[,
  .(b, comparison = model, metric = sprintf("%g-year AUC", TAU), difference = int - AUC)]
ed2 <- merge(UC[model == int, .(b, int = UnoC)], UC[model %chin% sing], by = "b")[, .(b, comparison = model,
  metric = sprintf("Uno C (≤%g years)", TAU), difference = int - UnoC)]
EDIFF <- rbind(ed1, ed2)[, .(difference = mean(difference), lo = quantile(difference, .025), hi = quantile(difference,
  .975)), by = .(comparison, metric)]

# The risk and performance summaries remain in CC, DS, EE, and EDIFF.
# Add project-specific summary export here if required.
