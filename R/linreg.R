#' linreg class
#'
#' creates multiple linear regression model
#' @field regcoeff regression coefficients
#' @field fitted_value for fitted values
#' @field df_name for input dataset name
#' @field residu for residuals
#' @field st_residu for standardized residuals
#' @field degree_of_freedom for degree of freedom
#' @field residual_variance for residual variance
#' @field var_reg_coeff for variance of regression coefficients
#' @field t_value for t value
#' @field pt_value for p value
#' @return none
#' @importFrom stats model.matrix pt
#' @importFrom ggplot2 ggplot aes
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
                                      st_residu="numeric",
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
                              st_residu <<- sqrt(abs(as.numeric(residu)/sd(as.numeric(residu))))
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
                            df_1<- data.frame (residu  =as.vector(residu),
                                               fitted_value =as.vector(fitted_value)
                            )
                            df_2<- data.frame (st_residu  =as.vector(st_residu),
                                               fitted_value =as.vector(fitted_value)
                            )
                            
                            residual_fitted <- ggplot(df_1,aes(x=fitted_value,y=residu))+
                              
                              geom_point(shape=1,size=4,stroke = 1)+
                              #geom_text(nudge_x = 0.21, nudge_y = 0.01,check_overlap = T)+
                              
                              geom_hline(yintercept = 0, linetype='dotted',colour = 'grey') +
                              
                              stat_summary(fun=median, geom="line", colour = "red")+
                              xlab("Fitted values\nlm(Petal.Length ~ Species)")+
                              ylab("Residuals")+
                              ggtitle("Residuals vs Fitted")+
                              theme(plot.title = element_text(hjust = 0.5),
                                    axis.line = element_line(colour = "black"),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.background = element_blank(),
                                    panel.border = element_rect(fill='transparent')
                                    
                              )
                            
                            standresidual_fitted <- ggplot(df_2,aes(x=fitted_value,y=st_residu))+
                              
                              geom_point(shape=1,size=4,stroke = 1)+
                              #geom_text(nudge_x = 0.21, nudge_y = 0.01,check_overlap = T)+
                              
                              geom_hline(yintercept = 0, linetype='dotted',colour = 'grey') +
                              
                              stat_summary(fun=median, geom="line", colour = "red")+
                              xlab("Fitted values\nlm(Petal.Length ~ Species)")+
                              ylab("Standardized Residuals")+
                              ggtitle("Scale - Location")+
                              theme(plot.title = element_text(hjust = 0.5),
                                    axis.line = element_line(colour = "black"),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.background = element_blank(),
                                    panel.border = element_rect(fill='transparent')
                                    
                              )
                            plots <- list(residual_fitted,standresidual_fitted)
                            plots
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
