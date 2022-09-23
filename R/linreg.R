#' linreg class
#'
#' creates multiple linear regression model
#' @field regcoeff regression coefficients
#' @field fitted_value for fitted values
#' @field df_name for input dataset name
#' @field residu for residuals
#' @field degree_of_freedom for degree of freedom
#' @field residual_variance for residual variance
#' @field var_reg_coeff for variance of regression coefficients
#' @field t_value for t value
#' @field pt_value for p value
#' @return none
#' @importFrom stats model.matrix pt
#' @import methods
#' @exportClass linreg
#' @export linreg
#' @export
#'
linreg <- setRefClass("linreg",
                        fields = list(X = "matrix",
                                      y = "numeric",
                                      formula = "formula",
                                      data = "data.frame",
                                      regcoeff = "matrix",
                                      fitted_value = "matrix",
                                      df_name = "character",
                                      residu = "matrix",
                                      degree_of_freedom = "numeric",
                                      residual_variance = "matrix",
                                      var_reg_coeff = "matrix",
                                      t_value = "matrix",
                                      pt_value = "matrix"
                        ),
                        methods = list(
                          initialize = function(formula, data){
                              formula <<- formula
                              data <<- data
                              X <<- model.matrix(formula, data)
                              y <<- unlist(data[all.vars(formula)[1]])
                              df_name <<- deparse(substitute(data))
                              regcoeff <<- ((solve(t(X) %*% X)) %*% t(X)) %*% y
                              fitted_value <<- X %*% regcoeff
                              residu <<- y - fitted_value
                              degree_of_freedom <<- nrow(X) - ncol(X)
                              residual_variance <<- t(residu) %*% residu / degree_of_freedom
                              var_reg_coeff <<- residual_variance[1] * solve(t(X) %*% X)
                              t_value <<- regcoeff / sqrt(diag(var_reg_coeff))
                              pt_value <<- 2 * (pt(abs(t_value), degree_of_freedom, lower.tail = FALSE))
                          },
                          print = function(){
                            form <- deparse(formula)
                            cat(paste("Call:\nlinreg(formula = ", form, ", data = ",df_name, ")\n\n", sep = ""))
                            cat("Coefficients:\n")
                            coeff <- as.vector(regcoeff)
                            names(coeff) <- colnames(X)
                            print.default(coeff)
                          },
                          plot = function(){
                            cat("The addition of a and b : ",a+b)
                          },
                          resid = function(){
                            print.default(as.vector(residu))
                          },
                          coef = function(){
                            cat("Coefficients:\n")
                            coeff <- as.vector(regcoeff)
                            names(coeff) <- colnames(X)
                            print.default(coeff)
                          },
                          pred = function(){
                            print.default(as.vector(fitted_value))
                          },
                          summary = function(){
                            signif_codes <- vector()
                            for(i in 1:length(pt_value)){
                              if(pt_value[i] < 0.001){
                                signif_codes[i] <- "***"
                                }
                              else if(pt_value[i] < 0.01){
                                signif_codes[i] <- "**"
                                }
                              else if(pt_value[i] < 0.05){
                                signif_codes[i] <- "*"
                                }
                              else if(pt_value[i] < 0.1){
                                signif_codes[i] <- "."
                                }
                            }
                            df<- data.frame(as.vector(regcoeff), as.vector(diag(sqrt(abs(var_reg_coeff)))), t_value, pt_value, signif_codes)
                            names(df) <- c("Estimate", "Std. Error", "t value", "p value", "")
                            cat("Coefficients:", "\n")
                            print.data.frame(df)
                            cat("Residual standard error:", as.vector(sqrt(residual_variance)[1]), "on", degree_of_freedom, "degrees of freedom")
                          }
                        )
  )


