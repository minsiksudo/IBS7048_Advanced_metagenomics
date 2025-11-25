# SIAMCAT

if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
BiocManager::install("SIAMCAT")

library("SIAMCAT")

data("feat_crc_zeller", package="SIAMCAT")
data("meta_crc_zeller", package="SIAMCAT")

feat.crc.zeller[1:3, 1:3]


dim(feat.crc.zeller)

head(meta.crc.zeller)

#Creation of labels
label.crc.zeller <- create.label(meta=meta.crc.zeller,
                                 label='Group', case='CRC')

#SIAMCAT object creation
sc.obj <- siamcat(feat=feat.crc.zeller,
                  label=label.crc.zeller,
                  meta=meta.crc.zeller)

show(sc.obj)

sc.obj <- filter.features(sc.obj,
                          filter.method = 'abundance',
                          cutoff = 0.001)

sc.obj <- check.associations(sc.obj, log.n0 = 1e-06, alpha = 0.05)
association.plot(sc.obj, sort.by = 'fc', 
                 panels = c('fc', 'prevalence', 'auroc'))

check.confounders(sc.obj, fn.plot = 'confounder_plots.pdf',
                  meta.in = NULL, feature.type = 'filtered')


sc.obj <- normalize.features(sc.obj, norm.method = "log.unit",
                             norm.param = list(log.n0 = 1e-06, n.p = 2,norm.margin = 1))

sc.obj <-  create.data.split(sc.obj, num.folds = 5, num.resample = 2)


sc.obj <- train.model(sc.obj, method = "lasso")

model_type(sc.obj)

models <- models(sc.obj)
models[[1]]$model

sc.obj <- make.predictions(sc.obj)
pred_matrix <- pred_matrix(sc.obj)

sc.obj <-  evaluate.predictions(sc.obj)
model.evaluation.plot(sc.obj)

model.interpretation.plot(sc.obj, fn.plot = 'interpretation.pdf',
                          consens.thres = 0.5, limits = c(-3, 3), heatmap.type = 'zscore')

