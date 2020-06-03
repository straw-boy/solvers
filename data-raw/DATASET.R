# code to prepare `DATASET` dataset goes here

# bodyfat (gaussian) ------------------------------------------------------
library(Matrix)

temp_file <- tempfile()

download.file(
  "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression/bodyfat",
  temp_file
)

tmp <- e1071::read.matrix.csr(temp_file, fac = FALSE)
unlink(temp_file)

data <- as(tmp$x,"sparseMatrix")
data2 <- as.matrix(data)

colnames(data2) <- c("siri_1956",
                     "age",
                     "weight",
                     "height",
                     "neck",
                     "chest",
                     "abdomen",
                     "hip",
                     "thigh",
                     "knee",
                     "ankle",
                     "biceps",
                     "foream",
                     "wrist")

# use the Siri 1956 equation as response
bodyfat <- list(x = data2[, -1], y = data2[, 1])

usethis::use_data(bodyfat, overwrite = TRUE)

# abalone (poisson) -------------------------------------------------------

library(SparseM)
library(Matrix)
library(caret)

temp_file <- tempfile()

download.file(
  "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression/abalone",
  temp_file
)

tmp <- e1071::read.matrix.csr(temp_file, fac = FALSE)
unlink(temp_file)
tmp_x <- as.data.frame(as.matrix(tmp$x))
tmp_x$V1 <- as.factor(tmp_x$V1)

tmp_x <- as.data.frame(model.matrix(~ ., tmp_x))[, -1]
colnames(tmp_x) <- c("sex",
                     "infant",
                     "length",
                     "diameter",
                     "height",
                     "weight_whole",
                     "weight_shucked",
                     "weight_viscera",
                     "weight_shell")

# randomly select a subset of the rows
part <- caret::createDataPartition(tmp$y, p = 0.05, list = FALSE)

abalone <- list(x = tmp_x[part, ], y = tmp$y[part])

usethis::use_data(abalone, overwrite = TRUE)

# heart (binomial) --------------------------------------------------------

library(SparseM)
library(Matrix)
library(fastDummies)
library(tidyverse)

temp_file <- tempfile()

download.file(
  "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/heart",
  temp_file
)

tmp <- e1071::read.matrix.csr(temp_file, fac = FALSE)
x <- as.data.frame(as.matrix(tmp$x))
colnames(x) <- c("age",
                 "sex",
                 "chest_pain",
                 "bp",
                 "chol",
                 "glucose",
                 "ecg",
                 "hr",
                 "angina",
                 "old_peak",
                 "slope",
                 "vessels",
                 "thal")

x2 <- x %>%
  mutate(sex = factor(sex, labels = c("male", "female")),
         cp = factor(chest_pain,
                     levels = c(4, 1, 2, 3),
                     labels = c("asymtompatic", "typical", "atypical", "nonanginal")),
         ecg = factor(ecg, labels = c("normal", "abnormal", "estes")),
         angina = as.factor(angina),
         glucose = factor(glucose, labels = c("low", "high")),
         slope = factor(slope,
                        levels = c(1, 2, 3),
                        labels = c("upsloping", "flat", "downsloping")),
         thal = factor(thal, labels = c("normal", "fixed", "reversible")))

x3 <- model.matrix(~ ., x2) %>%
  as.data.frame() %>%
  select(age,
         bp,
         chol,
         hr,
         old_peak,
         vessels,
         sex = sexfemale,
         angina = angina1,
         glucose_high = glucosehigh,
         cp_typical = cptypical,
         cp_atypical = cpatypical,
         cp_nonanginal = cpnonanginal,
         ecg_abnormal = ecgabnormal,
         ecg_estes = ecgestes,
         slope_flat = slopeflat,
         slope_downsloping = slopedownsloping,
         thal_fixed = thalfixed,
         thal_reversible = thalreversible)

x4 <- Matrix::Matrix(as.matrix(x3), sparse = TRUE)

# response
y <- factor(tmp$y, labels = c("absence", "presence"))

heart <- list(x = x4, y = y)

usethis::use_data(heart, overwrite = TRUE)

# wine (multiclass) -------------------------------------------------------

library(Matrix)
library(SparseM)

temp_file <- tempfile(fileext = ".csv")

download.file(
  "https://raw.githubusercontent.com/hadley/rminds/master/1-data/wine.csv",
  temp_file
)

tmp <- read.csv(temp_file)

x <- as.matrix(tmp[, -1])
y <- as.factor(tmp[, 1])

wine <- list(x = x, y = y)

usethis::use_data(wine, overwrite = TRUE)


# Student (multi-task) ----------------------------------------------------

tmp_file <- tempfile()
tmp_dir <- tempdir()

download.file(
  "https://archive.ics.uci.edu/ml/machine-learning-databases/00320/student.zip",
  tmp_file
)

unzip(tmp_file, exdir = tmp_dir)

d1 <-
  read.table(file.path(tmp_dir, "student-mat.csv"), sep = ";", header = TRUE)
d2 <-
  read.table(file.path(tmp_dir, "student-por.csv"), sep = ";", header = TRUE)
d3 <- merge(d1, d2, by = c("school",
                           "sex",
                           "age",
                           "address",
                           "famsize",
                           "Pstatus",
                           "Medu",
                           "Fedu",
                           "Mjob",
                           "Fjob",
                           "reason",
                           "nursery",
                           "internet"),
            suffixes = c("_math", "_port"))
y <- with(d3, cbind(G3_math, G3_port))
x1 <- subset(d3,
             select = c(
               "school",
               "sex",
               "age",
               "address",
               "famsize",
               "Pstatus",
               "Medu",
               "Fedu",
               "Mjob",
               "Fjob",
               "reason",
               "nursery",
               "internet"
             )
)
x1$famsize <- relevel(x1$famsize, ref = "LE3")
x2 <- model.matrix(~ ., x1)[, -1]

colnames(x2) <- c("school_ms",
                  "sex",
                  "age",
                  "urban",
                  "large_family",
                  "cohabitation",
                  "Medu",
                  "Fedu",
                  "Mjob_health",
                  "Mjob_other",
                  "Mjob_services",
                  "Mjob_teacher",
                  "Fjob_health",
                  "Fjob_other",
                  "Fjob_services",
                  "Fjob_teacher",
                  "reason_home",
                  "reason_other",
                  "reason_rep",
                  "nusery",
                  "internet")
colnames(y) <- c("math", "portugese")

student <- list(x = x2, y = y)

usethis::use_data(student, overwrite = TRUE)

unlink(tmp_file)


# # Covtype (binary) ----------------------------------------------------

library(Matrix)
tmp_file <- tempfile()

download.file(
  "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/covtype.libsvm.binary.bz2",
  tmp_file
)

tmp <- e1071::read.matrix.csr(tmp_file, fac = FALSE)

unlink(tmp_file)

x <- as(tmp$x,"sparseMatrix")
x2 <- as.matrix(x)

y <- as.matrix(tmp$y)
binary_covtype <- list(x = x2, y = y)

usethis::use_data(binary_covtype, overwrite = TRUE)

# # Golub ----------------------------------------------------

tmp <- tempfile()

download.file(
  "https://statweb.stanford.edu/~tibs/strong/ltrain.x",
  tmp
)

x <- t(matrix(scan(tmp), nrow = 7129, byrow = TRUE))
y <- double(38)
y[1:27] <- 0
y[28:38] <- 1

golub <- list(x = x, y = y)

usethis::use_data(golub, overwrite = TRUE)

unlink(tmp)

# # Gisette ----------------------------------------------------


library(Matrix)

tmp_data <- tempfile()
tmp_labels <- tempfile()

download.file(
  "https://statweb.stanford.edu/~tibs/strong/realdata/gisette_train.data",
  tmp_data
)

download.file(
  "https://statweb.stanford.edu/~tibs/strong/realdata/gisette_train.labels",
  tmp_labels
)

d <- scan(tmp_data, sep = " ")
x <- matrix(d, nrow = 6000, byrow = TRUE)[, 1:5000]
x <- as(x, "sparseMatrix")

dy <- scan(tmp_labels, sep = " ")
y <- (dy + 1)/2

gisette <- list(x = x, y = y)

usethis::use_data(gisette, overwrite = TRUE)

# # Dorothea ----------------------------------------------------

library(Matrix)

tmp_data <- tempfile()
tmp_labels <- tempfile()

download.file(
  "https://statweb.stanford.edu/~tibs/strong/realdata/dorothea_train.data",
  tmp_data
)

download.file(
  "https://statweb.stanford.edu/~tibs/strong/realdata/dorothea_train.labels",
  tmp_labels
)

d <- readLines(tmp_data, n = 800)

makevec <- function(d) {
  j <- 1
  j0 <- 1
  nn <- nchar(d)
  xx <- NULL
  while (j < nn) {
    while (substring(d, j, j) != " ") {
      j <- j + 1
    }
    xx <- c(xx, as.numeric(substring(d, j0, j - 1)))
    j0 <- j + 1
    j <- j + 1
  }
  xx
}

x <- matrix(0,
            nrow = 800,
            ncol = 100000,
            byrow = TRUE)

for (i in 1:800) {
  o <- makevec(d[i])
  x[i, o] <- 1
}

oo <- colSums(x)
x <- x[, oo > 0]

dy <- scan(tmp_labels, sep = " ")
y <- (dy + 1) / 2

x <- scale(x, TRUE, FALSE)

dorothea <- list(x = x, y = y)

usethis::use_data(dorothea, overwrite = TRUE)


# # cpusmall ----------------------------------------------------


library(e1071)
library(SparseM)

temp_file <- tempfile()

download.file(
  "https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/regression/cpusmall",
  temp_file
)

tmp <- e1071::read.matrix.csr(temp_file, fac = FALSE)

unlink(temp_file)

x <- as.matrix(tmp$x)
y <- tmp$y

cpusmall <- list(x = x, y = y)

usethis::use_data(cpusmall, overwrite = TRUE)