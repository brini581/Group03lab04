linreg1 <- function(formula, data){
  X <- model.matrix(formula, data)
  y <- unlist(data[all.vars(formula)[1]])
  regcoeff <- ((solve(t(X)%*% X)) %*% t(X))%*%y
  fitted_value <- X %*% regcoeff
  residu <- y - fitted_value
  degree_of_freedom <- nrow(X) - ncol(X)
  residual_variance <- (t(residu) %*% residu)/ degree_of_freedom
  var_reg_coeff <- as.numeric(residual_variance) * as.matrix(solve(t(X)%*% X))
  t_value <- (regcoeff / sqrt(var_reg_coeff[2,2])) [2]
  pt_value<- pt(t_value, degree_of_freedom, lower.tail = FALSE)
  obj <- linreg$new(X=X,y=y,regcoeff=regcoeff,fitted_value=fitted_value,residu=residu,degree_of_freedom=degree_of_freedom,residual_variance=residual_variance
                    ,var_reg_coeff=var_reg_coeff,t_value=t_value,pt_value=pt_value)
  return(obj)
}




linreg <- setRefClass("linreg",
                    fields = list(X="matrix",
                                  y="numeric",
                                  regcoeff="matrix",
                                  fitted_value="matrix",
                                  residu="matrix",
                                  degree_of_freedom="integer",
                                  residual_variance="matrix",
                                  var_reg_coeff="matrix",
                                  t_value="numeric",
                                  pt_value="numeric"
                                  ),
                    methods = list(
                      print=function(){
                        cat("The arguments are ",a)
                      },
                      plot=function(){
                        cat("The addition of a and b : ",a+b)
                      },
                      resid=function(){
                        cat(residu)
                      },
                      coef=function(){
                        cat(regcoeff)
                      },
                      pred=function(){
                        cat(fitted_values)
                      },
                      summary=function(){
                        cat("Coefficients:","\n")
                        cat("\t","\t","Estimate Std. Error t value Pr(>|t|)","\n")
                        cat("(Intercept)")
                        cat(X)
                      }
                    )
)

