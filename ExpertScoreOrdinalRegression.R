library(readxl)
library(Rcpp)
library(ordinal)
library(Metrics)
library(ggplot2)

#install.packages(c("readxl", "Rcpp", "ordinal", "Metrics"))

files <- list.files("W:/PYWELL_SHARED/Pywell Projects/BRC/Oli Pescott/235 1km SDMs/validationResults/",
                    pattern = "_FINAL.x",
                    full.names = T)

dat <- lapply(files,
              read_excel)

dat <- data.frame(do.call("rbind", dat))

head(dat)

dat <- dat[, c("charlieName", "monads", "q3", "q6")]

colnames(dat)[1] <- "species"

dat <- dat[-1,]

dat <- dat[-which(dat$q6 == "-"), ]

dat$q6 <- as.factor(dat$q6)

files <- list.files("W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/SDMValidation/averageSkill/",
                    full.names = T)

skill <- lapply(files,
                read.csv)

skill <- do.call("rbind", skill)

dat <- merge(dat, skill, by = "species", all.y = FALSE)

files <- list.files("W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/uncertaintySums/",
                          full.names = T)

uncertainty <- lapply(files,
                      read.csv)

uncertainty <- do.call("rbind", uncertainty)

dat <- merge(dat, uncertainty, by = "species", all.y = FALSE)

files <- list.files("W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/TSDA/specialisationIndices/",
                      full.names = T)

special <- lapply(files,
                  read.csv)

special <- do.call("rbind", special)

colnames(special)[1] <- "species"

dat <- merge(dat, special, by = "species", all.y = FALSE)

usdm::vif(dat[,c("monads", "auc", "sumSD", "sumCV", "SSI")])

ff <- q6 ~ auc + sumSD + monads

fr <- q6 ~ auc + sumSD + monads + (1|group)

fitAndEvaluate <- function(formula, mixedMod) {

  if (mixedMod == TRUE) {

    mod <- ordinal::clmm(formula,
                        data = dat)

  } else {

    mod <- ordinal::clm(formula,
                        data = dat)

  }


  print(summary(mod))

  #if (formula == ff) fit <- ordinal::nominal_test(mod)

  preds <- predict(mod,
                   newdat = dat,
                   type = "class")


  predObs <- cbind(as.numeric(dat$q6), as.numeric(preds$fit))

  cor <- cor.test(as.numeric(dat$q6),
                  as.numeric(preds$fit),
                  method = "spearman")

  mae <- Metrics::mae(predObs[,1], predObs[,2])

  out <- data.frame(formula = as.character(formula),
                    cor = cor[4],
                    mae = mae)

}

evalFF <- fitAndEvaluate(ff, FALSE)

evalFR <- fitAndEvaluate(fr, TRUE)

ggplot(data = dat, aes(x = as.numeric(q6))) +
  geom_histogram() +
  facet_wrap(~group) +
  theme_linedraw() +
  xlab("Accuracy score")
