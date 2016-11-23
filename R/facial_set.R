facial_set <-
function (data, formula) { 

  stopifnot(is.data.frame(data))

  formula <- as.formula(formula)
  freqColName <- as.character(attributes(terms(formula))$variables[2])
  y <- data[,which(colnames(data) == freqColName)]

  if (sum(y>0) < 2)
    stop ("Error: Contingency table must have at least 2 non-zero cells")

  #lpSolveAPI <- require(lpSolveAPI)
  #if (!lpSolveAPI) 
  #  stop ("exact=T needs lpSolveAPI package")

  value <- list ()

  X <- model.matrix(formula, data = data)
  
  z <- array (0, c(length(y)))
  z[y > 0] <- 1
  t <- t(X) %*% z
  Xt <- t(X)

  A <- array (1, c(length(y)))
  A <- as.vector(A)
  A[y > 0] <- 0

  value <- list ()

  if(sum(A)==0) {

    value$formula <- formula
    value$model.dimension <- dim(X)[2]     
    value$status <- "A started empty"
    value$iterations <- 0   
    value$face.dimension <- value$model.dimension
   
    fit <- glm (formula, data = data, family = poisson())
    value$maxloglik <- -sum(new_y) + sum(new_y * log(fit$fitted.values))
    
    facial_set <- z
    value$face <- cbind(data,facial_set)
   
    return (value)
  } 
    
  iter <- 1

  lp <- make.lp (0, dim(Xt)[2])

  for (i in 1:dim(Xt)[1])
    add.constraint(lp, Xt[i,], "=", t[i]) 

  while(T) {

    set.objfn(lp, -A)
    solve(lp)
    A[get.variables(lp) > 0] <- 0
      
    if(sum(A)==0) {

      value$formula <- formula
      value$model.dimension <- dim(X)[2]
      value$status <- "A became empty"
      value$iterations <- iter
      value$face.dimension <- value$model.dimension
      
      facial_set <- 1 - A
      value$face <- cbind(data,facial_set)

      new_data <- data[facial_set == 1,,drop = F]  
      new_y <- y[facial_set == 1]
      fit <- glm (formula, data = new_data, family = poisson())
      value$face.dimension <- fit$rank
      value$maxloglik <- -sum(new_y) + sum(new_y * log(fit$fitted.values))
      
      return(value)  

    }
      
    if(get.objective(lp) == 0) {

      value$formula <- formula
      value$model.dimension <- dim(X)[2]
      value$status <- "Optimal objective value 0"
      value$iterations <- iter
      
      facial_set <- 1 - A
      value$face <- cbind(data, facial_set)

      new_data <- data[facial_set == 1,,drop = F]  
      new_y <- y[facial_set == 1]

      fit <- glm (formula, data = new_data, family = poisson())
      value$face.dimension <- fit$rank
      value$maxloglik <- -sum(new_y) + sum(new_y * log(fit$fitted.values)) 
        
      return(value)
     
    }

    iter <- iter + 1   

  }
}
