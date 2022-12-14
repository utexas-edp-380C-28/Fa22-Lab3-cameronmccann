#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 3: An Overview of Linear Regression 
#        and Data Generation
#
# Name:  Cameron McCann
#---------------------------------------------------#

## SETUP -------------------------------------------
# Source rmvnorm()
source("scripts/rmvnorm.R")

## Begin Code --------------------------------------

## 1.a
# Q1) {Set parameters}
p_x <- list(mu = -2,
            sigma = 3)

p_y <- list(beta_0 = 12, 
            beta_1 = 4, 
            sigma = 7)

# Q2) 
generate_X_Y <- function(n, p_x, p_y) {
  
  X <- rnorm(n = n, 
             mean = p_x$mu, 
             sd = p_x$sigma)
  
  Y <- rnorm(n = n, 
        mean = p_y$beta_0 + p_y$beta_1 * X, 
        sd = p_y$sigma) 
  
  # output should be n x 2
  cbind(x = X, y = Y)
}

# Q3) 
analyze <- function(data_set) {
  
  # Run linear regression 
  data_set_lm <- lm(y ~ x, data_set)
  
  #
  means <- sapply(data_set_lm$model, mean)
  names(means) <- c("m_y", "m_x")

  #SDs
  sds <- sapply(data_set_lm$model, sd)
  names(sds) <- c("s_y", "s_x")

  #cor
  cor <- cor(data_set_lm$model)[2, 1]
  names(cor) <- c("cor")

  #b0 & b1
  betas <- c(data_set_lm$coefficients[1], data_set_lm$coefficients[2])
  names(betas) <- c("b0", "b1")

  #s_e
  s_e <- sd(data_set_lm$residuals)
  names(s_e) <- c("s_e")

  c(means, sds, cor, betas, s_e)
  
  # function 
    #input = data_set; 
    #output = vector with following labels c(m_x, m_y, s_x, s_y, cor, b0, b1, s_e)
}

# Q4, Q5, & Q6) 
set.seed(17290)

sim_data <- replicate(n = 500, generate_X_Y(n = 100, 
                                            p_x = p_x,
                                            p_y = p_y))

# pre-fill results dataset
## {can do this with: results <- apply(rep_data, 3, analyze)}
results <- data.frame(matrix(1, 500, 8))

for(i in 1:dim(sim_data)[3]) {
  
  results[i, ] <- analyze(as.data.frame(sim_data[, , i]))
  
  # add column names to results table
  colnames(results) <- c("m_y", "m_x", "s_y", "s_x", "cor", "b0", "b1", "s_e")
}


# Answer
list(mean = apply(results, 2, mean), 
     sd = apply(results, 2, sd))

# $mean
# m_y        m_x        s_y        s_x        cor         b0         b1 
# 3.9810216 -2.0061317 13.8596003  3.0034995  0.8631812 11.9780298  3.9873856 
# s_e 
# 6.9396095 
# 
# $sd
# m_y        m_x        s_y        s_x        cor         b0         b1 
# 1.41924855 0.31115637 1.03767609 0.20911389 0.02432744 0.86119359 0.24140973 
# s_e 
# 0.48821891 


## 2.a

# OLS function
OLSfunc <- function(y, X) {
  
  # Drop 1st column if all 1s
  ifelse(sum(X[, 1]) == nrow(X), 
         X,
         X <- cbind(rep(1, nrow(X)), X)
  )
  
  # Estimating Bs
  Betas <- solve(crossprod(X)) %*% crossprod(X, y) 
  
  # Hat matrix
  # H_Matrix <- X %*% solve(t(X) %*% X) %*% t(X)
  
  # Residuals
  resid <- y - (X %*% Betas)
  
  # SD(e)
  SD_e <- sqrt( (t(resid) %*% resid) / (nrow(X) - ncol(X)) )
  
  # R^2 
  R2 <- (t(Betas) %*% cov(X) %*% Betas) / 
    ((t(Betas) %*% cov(X) %*% Betas) + SD_e^2)
  
  # SE for Betas
  SE_b <- sqrt( SD_e^2 %*% diag(solve(crossprod(X))) )
  
  # t values
  t_vals <- Betas / t(SE_b)
  
  # DF residuals
  dfs_residuals <- (nrow(X) - nrow(Betas))
  
  # Create output
  output <-   data.frame(Estimate = c(Betas, SD_e, R2),
                         SE = c(SE_b, "NA", "NA"),
                         t_value = c(t_vals, "NA", "NA"),
                         p = c(2 * pt(abs(t_vals), dfs_residuals, lower.tail = FALSE), "NA", "NA"))
  
  # row names 
  rownames(output) <- c(paste0(rep("b", nrow(Betas)), 0:(nrow(Betas)-1)), 
                        "SD(e)", "R2")
  # column names 
  colnames(output) <- c("Estimate", "SE", "t value", "p")
  
  # output
  return(output)
}

# Answer 
summary(lm(mpg ~ wt + cyl + gear, data = mtcars))

OLSfunc(as.matrix(mtcars$mpg), as.matrix(mtcars[, c("wt", "cyl", "gear")]))

# Estimate                SE            t value                    p
# b0    42.3863641  4.37899519903836   9.67947260714761  1.9659287164877e-10
# b1    -3.3920819 0.820802491744326  -4.13264087599804 0.000294156959036295
# b2    -1.5280010 0.419753253251725  -3.64023630377079  0.00109260914427508
# b3    -0.5228629 0.778870339652137 -0.671309355170919    0.507524370714282
# SD(e)  2.5921848                NA                 NA                   NA
# R2     0.8182681                NA                 NA                   NA


## 3.a
set.seed(21389)
# Predefine parameters 
param_list <- list(Rho = 0.3,
                   t_mu = c(5, 10),
                   t_sigma = c(1, 2))

# Function that translates parameter list into rmvnorm() inputs 
param2rmvnorm_list <- function(param_list) {
  
  # Create R_x & change diag
  R_x <- matrix(param_list$Rho,
                dim(diag(c(param_list$t_sigma)))[1],
                dim(diag(c(param_list$t_sigma)))[2])
  diag(R_x) <- 1
  
  # Create list
  rmvnorm_list <- list(
    ## Mean vector
    mu = matrix(param_list$t_mu, 
                ncol = 1), 
    ## Sigma matrix 
    Sigma = diag(c(param_list$t_sigma)) %*% 
      R_x %*% 
      diag(c(param_list$t_sigma))
  )
  
  return(rmvnorm_list)
}

# Change parameter list to rmvnorm inputs
pred_rmvnorm_list <- param2rmvnorm_list(param_list)

# Generate predictor data 
pred_data <- rmvnorm(n = 100000, 
        mu = pred_rmvnorm_list$mu,
        Sigma = pred_rmvnorm_list$Sigma)


## 3.b
set.seed(23921)
# Predefine outcome parameters
outcome_param_list <- list(beta1 = 1,
                           beta2 = 1, 
                           R2 = 0.6,
                           mu_y = 10, 
                           sd_y = 5)

# Betas
Betas_3b <- matrix(as.double(
  c(outcome_param_list[startsWith(names(outcome_param_list), 
                                  prefix = "beta")]))
)


# Compute SD error for Y
SD_e <- sqrt(
  t( Betas_3b ) %*% 
    pred_rmvnorm_list$Sigma %*% 
    Betas_3b %*% 
    ( (1/outcome_param_list$R2) - 1)
)

# Find & set intercept
intercept_3b <- outcome_param_list$mu_y - 
  t( pred_rmvnorm_list$mu ) %*% Betas_3b # intercept = -5

# Generate Y data
outcome_data <- rnorm(n = 100000, 
      mean = as.vector(intercept_3b) + pred_data %*% Betas_3b,
      sd = SD_e)

# Answer 

## Analyze with OLS function 
(ols_output_3b <- OLSfunc(outcome_data, pred_data))

## Difference in Pop & Estimates
Diff_Tbl_3b <- as.data.frame(cbind(
  Pop = c(intercept_3b, 
          outcome_param_list[startsWith(names(outcome_param_list),
                                        prefix = "beta")],
          SD_e,
          outcome_param_list$R2),
  Est = c(ols_output_3b$Estimate)))

### Add rownames 
rownames(Diff_Tbl_3b) <- c(paste0("beta", 0:(nrow(Diff_Tbl_3b)-3)),
                           "SD(e)", 
                           "R2")

### Add difference column 
Diff_Tbl_3b$Diff <- c(unlist(Diff_Tbl_3b$Pop) - unlist(Diff_Tbl_3b$Est))

Diff_Tbl_3b

# Estimate                  SE           t value  p
# b0    -4.9370600  0.0404634289861931 -122.012892517685  0
# b1     0.9891621 0.00677154296619353  146.076328907351  0
# b2     0.9997359 0.00337368817415204  296.333224140717  0
# SD(e)  2.0361731                  NA                NA NA
# R2     0.5978875                  NA                NA NA

# Pop       Est          Diff
# beta0      -5  -4.93706 -0.0629399882
# beta1       1 0.9891621  0.0108378625
# beta2       1 0.9997359  0.0002641061
# SD(e) 2.03306  2.036173 -0.0031129887
# R2        0.6 0.5978875  0.0021125328


## 3.c
set.seed(123782)
# Predefine outcome parameters
outcome2_param_list <- list(Rho_y1 = 0.3,
                            Rho_y2 = -0.4,
                            mu_y = 10,
                            SD_y = 5)

# Covariances (Sigma_Xy)
## Change corr to covariances by multiplying by the product of x & y SD
Sigma_Xy_3c <- c(outcome2_param_list$Rho_y1 * 
                   (outcome2_param_list$SD_y * param_list$t_sigma[1]), 
                 
                 
                 outcome2_param_list$Rho_y2 * 
                   (outcome2_param_list$SD_y * param_list$t_sigma[2])
)

## Compute betas
Betas_3c <- solve(param2rmvnorm_list(param_list)$Sigma) %*% Sigma_Xy_3c

# Residual variance 
var_e <- outcome2_param_list$SD_y^2 - t(Betas_3c) %*% Sigma_Xy_3c

# Find & set intercept
intercept_3c <- outcome2_param_list$mu_y - 
  t( pred_rmvnorm_list$mu ) %*% Betas_3c 

# R^2
R2 <- (t(Betas_3c) %*% param2rmvnorm_list(param_list)$Sigma %*% Betas_3c) / 
  ((t(Betas_3c) %*% param2rmvnorm_list(param_list)$Sigma %*% Betas_3c) + var_e)

# Generate Y data
outcome2_data <- rnorm(n = 100000, 
      mean = as.vector(intercept_3c) + pred_data %*% Betas_3c,
      sd = sqrt(var_e) )

# Answer 

## Analyze with OLS function 
(ols_output_3c <- OLSfunc(outcome2_data, pred_data))

## Difference in Pop & Estimates
Diff_Tbl_3c <- data.frame(Pop = c(intercept_3c, 
                                  Betas_3c,
                                  sqrt(var_e),
                                  R2),
                          Est = ols_output_3c$Estimate)

### Add difference column 
Diff_Tbl_3c$Diff <- Diff_Tbl_3c$Pop - Diff_Tbl_3c$Est

### Add rownames 
rownames(Diff_Tbl_3c) <- c(paste0("b", 0:2),
                             "SD(e)", 
                             "R2")

Diff_Tbl_3c

# Estimate                  SE           t value  p
# b0    11.7786196   0.079813607389889  147.576585329678  0
# b1     2.3182149  0.0133567837741075  173.560861847616  0
# b2    -1.3368715 0.00665455771134872 -200.895625192709  0
# SD(e)  4.0163259                  NA                NA NA
# R2     0.3524265                  NA                NA NA

#       Pop        Est         Diff
# b0    11.9230769 11.7786196  0.144457282
# b1     2.3076923  2.3182149 -0.010522596
# b2    -1.3461538 -1.3368715 -0.009282314
# SD(e)  4.0191848  4.0163259  0.002858840
# R2     0.3538462  0.3524265  0.001419650


## 3.d

### 1)

#### Method 1 function 
generate_Y_MultiReg_method1 <- function(n, p_x, p_y) {
  
  # Generate predictors
  ## Change parameter list to rmvnorm inputs
  p_x_list <- param2rmvnorm_list(p_x)
  
  ## Generate predictor data 
  pred_data <- rmvnorm(n = n, 
                       mu = p_x_list$mu,
                       Sigma = p_x_list$Sigma)
  
  
  # Generate outcome data
  ## Betas
  Betas <- matrix(as.double(
    c(p_y[startsWith(names(p_y), 
                                    prefix = "beta")])
  ))
  
  ## Compute SD error for Y
  SD_e <- sqrt(
    t( Betas ) %*% 
      p_x_list$Sigma %*% 
      Betas %*% 
      ( (1/p_y$R2) - 1)
  )
  
  ## Find & set intercept
  intercept <- p_y$mu_y - 
    t( p_x_list$mu ) %*% Betas 
  
  ## Generate Y data
  outcome_data <- rnorm(n = n, 
                        mean = as.vector(intercept) + pred_data %*% Betas,
                        sd = SD_e)
  
  # Create output list object 
  output_data <- list(y = outcome_data,
                      X = pred_data)
  
  # Add population parameter attributes to output 
  attr(output_data, "b0") <- intercept
  attr(output_data, "b1") <- Betas[1]
  attr(output_data, "b2") <- Betas[2]
  attr(output_data, "b3") <- Betas[3]
  attr(output_data, "b4") <- Betas[4]
  attr(output_data, "b5") <- Betas[5]
  attr(output_data, "SD(e)") <- SD_e
  attr(output_data, "R2") <- p_y$R2
  
  # ## Loop relevant population parameters for predictors
  # for (i in seq_along(p_x)) {
  #   attr(output_data, names(p_x[i])) <- p_x[[i]]
  # }
  # 
  # ## Loop relevant population parameters for outcome
  # for (i in seq_along(p_y)) {
  #   attr(output_data, names(p_y[i])) <- p_y[[i]]
  # }

  
  return(output_data)
  
}

#### Method 2 function 
generate_Y_MultiReg_method2 <- function(n, p_x, p_y) {
  
  # Generate predictors
  ## Change parameter list to rmvnorm inputs
  p_x_list <- param2rmvnorm_list(p_x)
  
  ## Generate predictor data 
  pred_data <- rmvnorm(n = n, 
                       mu = p_x_list$mu,
                       Sigma = p_x_list$Sigma)
  
  
  # Generate outcome data
  ## Change corr to covariances by multiplying by the product of x & y SD
  ### populate vector 
  Sigma_Xy <- rep(1, sum(startsWith(names(p_y), prefix = "Rho_y")))
  
  ### loop to convert to cov mat
  for (i in 1:sum(startsWith(names(p_y), prefix = "Rho_y"))) {
    
    Sigma_Xy[i] <- ( p_y[[paste0("Rho_y", i)]] * 
                       (p_y$sd_y * p_x$t_sigma[i]) )
  }
  
  ## Compute betas
  Betas <- solve(param2rmvnorm_list(p_x)$Sigma) %*% Sigma_Xy
  
  # Residual variance 
  var_e <- p_y$sd_y^2 - t(Betas) %*% Sigma_Xy
  
  # Find & set intercept
  intercept <- p_y$mu_y - t(p_x$t_mu) %*% Betas
  
  # R^2
  R2 <- ((t(Betas) %*% p_x_list$Sigma %*% Betas) / 
           ((t(Betas) %*% p_x_list$Sigma %*% Betas) + var_e))
  
  # # R^2 
  # R2 <- (t(Betas) %*% cov(X) %*% Betas) / 
  #   ((t(Betas) %*% cov(X) %*% Betas) + SD_e^2)
  
  # Generate Y data
  outcome_data <- rnorm(n = n,
                        mean = as.vector(intercept) + pred_data %*% Betas,
                        sd = sqrt(var_e))
  
  # Create output list object 
  output_data <- list(y = outcome_data,
                      X = pred_data)
  
  # Add population parameter attributes to output 
  
  # Add population parameter attributes to output 
  attr(output_data, "b0") <- intercept
  attr(output_data, "b1") <- Betas[1]
  attr(output_data, "b2") <- Betas[2]
  attr(output_data, "b3") <- Betas[3]
  attr(output_data, "b4") <- Betas[4]
  attr(output_data, "b5") <- Betas[5]
  attr(output_data, "SD(e)") <- sqrt(var_e)
  attr(output_data, "R2") <- R2
  
  # ## Loop relevant population parameters for predictors
  # for (i in seq_along(p_x)) {
  #   attr(output_data, names(p_x[i])) <- p_x[[i]]
  # }
  # 
  # ## Loop relevant population parameters for outcome
  # for (i in seq_along(p_y)) {
  #   attr(output_data, names(p_y[i])) <- p_y[[i]]
  # }
  
  
  return(output_data)
}


### 2)
set.seed(6972)

# Set parameters
p_x <- list(Rho = 0.15,
            #t_mu = c(2:6),
            t_mu = rep(0, 5),
            t_sigma = sqrt(c(1:5)))

p_y <- list(beta1 = 1,
            beta2 = 1,
            beta3 = 1,
            beta4 = 1,
            beta5 = 1,
            R2 = 0.5,
            mu_y = 25, 
            sd_y = 5)

# Generate Data
method1_output <- generate_Y_MultiReg_method1(n = 100000,
                                              p_x = p_x,
                                              p_y = p_y)

# Answer 
## Analyze with OLS function 
(ols_output_3d_2 <- OLSfunc(method1_output$y, method1_output$X))

## Difference table
### Create initial table
Diff_Tbl_3d_2 <- data.frame(
  Pop = c(attr(method1_output, "b0"),
          attr(method1_output, "b1"),
          attr(method1_output, "b2"),
          attr(method1_output, "b3"),
          attr(method1_output, "b4"),
          attr(method1_output, "b5"),
          attr(method1_output, "SD(e)"),
          attr(method1_output, "R2")),
  Est = ols_output_3d_2[c("b0", "b1", "b2", "b3", "b4", "b5", "SD(e)", "R2"), 
                        "Estimate"])

### Add difference column 
Diff_Tbl_3d_2$Diff <- Diff_Tbl_3d_2$Pop - Diff_Tbl_3d_2$Est

### Add rownames 
rownames(Diff_Tbl_3d_2) <- c(paste0("beta", 0:5),
                           "SD(e)", 
                           "R2")

Diff_Tbl_3d_2

# Estimate                  SE          t value                     p
# b0    5.1263156  0.0807062930012819 63.5181640342362                     0
# b1    0.9651298  0.0291283319587828  33.133712509388 1.94652030530563e-239
# b2    1.0262173   0.014582468605824 70.3733585805874                     0
# b3    0.9872700  0.0097345431646502 101.419247437614                     0
# b4    0.9862443 0.00729792681138579 135.140339809386                     0
# b5    1.0042094  0.0058721056923955 171.013509573655                     0
# SD(e) 8.9583463                  NA               NA                    NA
# R2    0.4982514                  NA               NA                    NA

# Pop        Est         Diff
# beta0 25.000000 25.0219373 -0.021937347
# beta1  1.000000  0.9812441  0.018755865
# beta2  1.000000  1.0199428 -0.019942757
# beta3  1.000000  0.9881404  0.011859598
# beta4  1.000000  0.9852023  0.014797719
# beta5  1.000000  1.0050628 -0.005062766
# SD(e)  4.825922  4.8184817  0.007440416
# R2     0.500000  0.4986015  0.001398544


### 3) 
set.seed(1237)

# Set parameters
p_x <- list(Rho = 0.15,
            t_mu = c(2:6),
            t_sigma = sqrt(c(1:5)))

p_y <- list(Rho_y1 = -0.15,
            Rho_y2 = -0.50,
            Rho_y3 = 0.15,
            Rho_y4 = 0.30,
            Rho_y5 = 0.20,
            mu_y = 10, 
            sd_y = 4)

# Generate Data
method2_output <- generate_Y_MultiReg_method2(n = 100000,
                            p_x = p_x,
                            p_y = p_y)

# Answer 
## Analyze with OLS function 
(ols_output_3d_3 <- OLSfunc(method2_output$y, method2_output$X))

## Difference table
### Create initial table
Diff_Tbl_3d_3 <- data.frame(
  Pop = c(attr(method2_output, "b0"),
          attr(method2_output, "b1"),
          attr(method2_output, "b2"),
          attr(method2_output, "b3"),
          attr(method2_output, "b4"),
          attr(method2_output, "b5"),
          attr(method2_output, "SD(e)"),
          attr(method2_output, "R2")),
  Est = ols_output_3d_3[c("b0", "b1", "b2", "b3", "b4", "b5", "SD(e)", "R2"), 
                        "Estimate"])

### Add difference column 
Diff_Tbl_3d_3$Diff <- Diff_Tbl_3d_3$Pop - Diff_Tbl_3d_3$Est

### Add rownames 
rownames(Diff_Tbl_3d_3) <- c(paste0("b", 0:5),
                             "SD(e)", 
                             "R2")

Diff_Tbl_3d_3

# Estimate                  SE           t value  p
# b0    11.0934528   0.025436973611998  436.115278599631  0
# b1    -0.7012153 0.00922298668829625 -76.0290890898014  0
# b2    -1.1709826 0.00460317900528182 -254.385638950379  0
# b3     0.2310588 0.00307372349751943  75.1722941986507  0
# b4     0.3498870 0.00230840578427237  151.570837049004  0
# b5     0.1890791 0.00185320994975727  102.027907379875  0
# SD(e)  2.8303933                  NA                NA NA
# R2     0.4974886                  NA                NA NA

#         Pop        Est          Diff
# b0     8.7180880  8.7381916 -2.010364e-02
# b1    -0.7058824 -0.7060328  1.504569e-04
# b2    -1.6637807 -1.6624466 -1.334024e-03
# b3     0.4075414  0.4039440  3.597378e-03
# b4     0.7058824  0.7032716  2.610748e-03
# b5     0.4209069  0.4222341 -1.327148e-03
# SD(e)  2.8284271  2.8211829  7.244233e-03
# R2     0.5000000  0.4999777  2.234008e-05
