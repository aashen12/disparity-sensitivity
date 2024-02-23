# library(quantreg)

estimate_qbounds <- function(weights, Qhat, Y, lambda, tau0){
    #Quantile regression here:
    df_lb = data.frame(Y, Qhat=Qhat[,1], weights = weights)
    df_ub = data.frame(Y, Qhat=Qhat[,2], weights = weights)

    model_quantile.lb = quantreg::rq(Y~ Qhat, tau = 1-tau0, weights = weights, data=df_lb)
    model_quantile.ub = quantreg::rq(Y ~ Qhat, tau = tau0, weights= weights, data = df_ub)

    cutoff.ub = predict(model_quantile.ub, df_ub)
    w_plus = weights*lambda^(sign(Y - cutoff.ub))
    Y_max = (sum((Y - cutoff.ub)*w_plus) + sum(cutoff.ub*weights))/sum(weights)

    cutoff.lb = predict(model_quantile.lb, df_lb)
    w_minus = weights*lambda^(sign(cutoff.lb-Y))
    Y_min = (sum((Y - cutoff.lb)*w_minus) + sum(cutoff.lb*weights))/sum(weights)
    return(c(Y_min, Y_max))
}



run_bootstrap_iter <- function(df_sample, lambda, bootstrap = TRUE){
    
    #Re-Sample:
    if (bootstrap == TRUE) {
        df_sample_bs = df_sample[sample(1:nrow(df_sample), replace=TRUE),]
    } else {
        df_sample_bs = df_sample
    }
    
    #df_sample_bs <- data.frame(Y,Z,X)
    Y = df_sample_bs$Y
    Z = df_sample_bs$Z
    #Estimate weights:
    #model_ps = WeightIt::weightit(Z~ .-Y, data = df_sample_bs, method='ebal', estimand="ATT")
    #weights = model_ps$weights
    model_ps = glm(Z~ .-Y, data = df_sample_bs, family = binomial)
    probs <- predict(model_ps, type = "response")
    weights = probs/(1-probs)

    #Estimate quantiles:
    tau0 = lambda/(1+lambda)
    #Estimate quantiles: (Qhat))
    rq_model = rq(Y~.-Z, tau=c(1-tau0, tau0), data= df_sample_bs[Z==0,])
    Qhat = predict(rq_model, df_sample_bs)
    if(lambda==1){
        Qhat = cbind(Qhat, Qhat)
    }

    #Estimate
    qbounds0 = estimate_qbounds(weights=weights[Z==0], Qhat[Z==0,], Y=Y[Z==0], lambda, tau0)


    return(mean(Y[Z==1]) - c(qbounds0[2], qbounds0[1]))
}

obtainCIsDG <- function(df_sample, lambda, nboot = 1000) {
    point_est <- data.frame(min = rep(0, nboot),
                            max = rep(0, nboot))
    
    for (i in 1:nboot) {
        point_est[i,] <- 
        tryCatch(
            {
                run_bootstrap_iter(df_sample = df_sample, lambda = lambda)
            },
            error=function(cond) {
                message(cond)
                # Choose a return value in case of error
                return(NA)
            }
        ) 
        
    }
    
    return(c(quantile(point_est[,1], probs = 0.025, na.rm = TRUE), quantile(point_est[,2], probs = 0.975, na.rm = TRUE)))
}
