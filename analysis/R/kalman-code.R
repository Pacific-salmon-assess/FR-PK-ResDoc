# R version of Kalman filter code for Ricker and Larkin model
# Authors: Brice MacGregor, Brian Pyper, Brigitte Dorner,
# School of Resource and Environmental Management
# Simon Fraser University
# Burnaby, British Columbia
# Canada V5A 1S6
#
# You are free to use this code for academic purposes. 
# For future inquiries, please contact: bdorner@driftwoodcove.ca

#-----------------------------------------------------------------------------------
"run.kalman.Ricker" <- function(data, initials = list(mean.a = 1, var.a = 1, b = -1, ln.sig.e = -1, ln.sig.w = -1, Ts = 1), referenceYear = 1950, period=1)
{
    # Purpose: A script for running a Kalman filter (random walk) analysis on the Ricker 'a' parameter for a single salmon stock
    # Assumes the Spawner and Recruits-per-Spawner data data$S and data$lnRS are aligned by brood year and that no brood years are skipped (though missing data for some years are ok
    # and should be represented by NAs)
    # For analysis of even and odd-year lines (pinks), or dominant and sub-dominant lines (sockeye), provide a referenceYear that is an "on" year, i.e.,
    # is to be included in the analysis, and the period to be used (e.g., period = 2 means every second year of data will be included)
	
	print(data)
	row.names(data) <- data$BY
	# remove NAs from beginning and end of record
	first.year <- min(data$BY)
	while (is.na(data[as.character(first.year), "lnRS"])) first.year <- first.year + 1
	last.year <- max(data$BY)
	while (is.na(data[as.character(last.year), "lnRS"])) last.year <- last.year - 1
	
	# only use years of the relevant line (the default period of 1 includes all years)
	valid.years <- first.year:last.year
	if (period > 1)
    valid.years <- valid.years[(valid.years %% period) == (referenceYear %% period)]
	
	valid.years <- as.character(valid.years)
	result <- kf.rw.Ricker(initials, data[valid.years, "S"], data[valid.years, "lnRS"])
	data$lnRS.smooth <- data$a.filtered <- data$a.smooth <- data$a.smooth.var <- rep(NA, nrow(data))
	data[valid.years, "lnRS.smooth"] <- result$smoothe.y
	data[valid.years, "a.filtered"] <- result$post.mean.a
	data[valid.years, "a.smooth"] <- result$smoothe.mean.a
	data[valid.years, "a.smooth.var"] <- result$smoothe.var.a
    
	output.summary.stats <- list(b = result$b, sig.e = result$sig.e, sig.w = result$sig.w, N.tot = result$N.tot, N.cond = result$N.cond,
    Param = result$Param, AICc = result$AICc, Report = result$Report)
	
	list(df = data, summary = output.summary.stats)
}


#-----------------------------------------------------------------------------------
"run.kalman.Larkin" <- function(data, initials = list(mean.a = 1, var.a = 1, b = -1, b1 = 0, b2 = 0, b3 = 0, ln.sig.e = -1, ln.sig.w = -1, Ts = 1), referenceYear = 1950, period=1)
{
    # Purpose: A script for running a Kalman filter (random walk) analysis on the Larkin 'a' parameter for a single salmon stock
    # Assumes the Spawner and Recruits-per-Spawner data data$S and data$lnRS are aligned by brood year and that no brood years are skipped (though missing data for some years are ok
    # and should be represented by NAs)
    # Since this function fits a Larkin model, spawner counts for year t-1 (S.t1), year t-2 (S.t2), and year t-3 (S.t3) must also be given as part of the input data.
    # For analysis of dominant and sub-dominant lines, provide a referenceYear that is an "on" year, i.e.,
    # is to be included in the analysis, and the period to be used (e.g., period = 2 means every second year of data will be included)
	
	row.names(data) <- data$BY
	# remove NAs from beginning and end of record
	first.year <- min(data$BY)
	while (is.na(data[as.character(first.year), "lnRS"]) || any(is.na(data[as.character(first.year), c("S.t1", "S.t2", "S.t3")])))
    first.year <- first.year + 1
	last.year <- max(data$BY)
	while (is.na(data[as.character(last.year), "lnRS"]) || any(is.na(data[as.character(last.year), c("S.t1", "S.t2", "S.t3")])))
    last.year <- last.year - 1
	
	# only use years of the relevant line (the default period of 1 includes all years)
	valid.years <- first.year:last.year
	if (period > 1)
	valid.years <- valid.years[(valid.years %% period) == (referenceYear %% period)]
	
	valid.years <- as.character(valid.years)
	result <- kf.rw.Larkin(initials, data[valid.years, c("S", "S.t1", "S.t2", "S.t3")], data[valid.years, "lnRS"])
	data$lnRS.smooth <- data$a.filtered <- data$a.smooth <- data$a.smooth.var <- rep(NA, nrow(data))
	data[valid.years, "lnRS.smooth"] <- result$smoothe.y
	data[valid.years, "a.filtered"] <- result$post.mean.a
	data[valid.years, "a.smooth"] <- result$smoothe.mean.a
	data[valid.years, "a.smooth.var"] <- result$smoothe.var.a
	output.summary.stats <- list(b = result$b, b1 = result$b1, b2 = result$b2, b3 = result$b3, sig.e = result$sig.e, sig.w = result$sig.w,
    N.tot = result$N.tot, N.cond = result$N.cond, Param = result$Param, AICc = result$AICc, Report = result$Report)
	
	list(df = data, summary = output.summary.stats)
}

#-----------------------------------------------------------------------------------
"kf.rw.Ricker" <- function(initial, x, y)
{
    # Written for S-Plus 2000, Professional Release 2
    # Source of this code:
    # School of Resource and Environmental Management
    # Simon Fraser University
    # Burnaby, British Columbia
    # Canada V5A 1S6
    #
    # For future inquiries, contact:
    # Randall M. Peterman at the above School, either by
    #       e-mail at rmpasst@sfu.ca or by phone at 604 291-4683
    #
    # Authors of this code:
    # Brice MacGregor
    # February 14, 2002
    # Modified by Brian Pyper Oct 10, 2002
    # Modified by Brigitte Dorner Oct 29, 2003
	
    # Purpose:
    # ========
    # Finds ML estimates of kalman filter.  KF model assumes the following:
    # - Simple linear regression
    # - Time-varying intercept only
    # - Intercept follows a random walk
    # - Fits 3 PARAMETERS (b, sig.e, sig.w)
	
    # Arguments:
    # ==========
    # initial  LIST with initial parameter values and starting conditions for estimation procedure.
    #    Names for the elements of initial must be:
    #
    #     mean.a  Initial value for mean of intercept in recursive calculations
    #     var.a   Initial value for variance of intercept in recursive calculations
    #     b     Starting value of slope for ML estimation
    #     ln.sig.e, ln.sig.w Starting values for natural logarithms of error terms in
    #          observation and system equations
    #     Ts     Number of observations at start of data set to omit for
    #          calculation of variance in observation equation and concentrated
    #          likelihood function.
    # x, y   Data for the observation equation
	
    #   fit using general-purpose optimizer
    # fit <- optim(par=c(initial$b, initial$ln.sig.e, initial$ln.sig.w), fn=kalman.rw.fit.Ricker, gr=NULL,
    #	method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), lower = -Inf, upper = Inf,
    #   control = list(), hessian = FALSE,
    #	initial$mean.a, initial$var.a, x, y, initial$Ts)
    #cat("ML estimation for Ricker 'b' parameter and error terms using optim:")
    #print(fit)
	
    # fit using nlminb
	fit <- nlminb(start=c(initial$b, initial$ln.sig.e, initial$ln.sig.w), objective=kalman.rw.fit.Ricker,
    gradient = NULL, hessian = NULL, scale = 1, control = list(), lower = -Inf, upper = Inf,
    initial$mean.a, initial$var.a, x, y, initial$Ts)
	
	if (fit$convergence != 0)
	{
		print(fit$message)
		warning("ML estimation for Ricker 'b' parameter and error terms failed to converge!")
	}
	
    # Perform recursive calculations for KF with ML estimates:
	out <- kalman.rw.Ricker(initial$mean.a, initial$var.a, fit$par[1], fit$par[2], fit$par[3], x, y, initial$Ts)
	N <- length(x) - sum(is.na(x))
	param <- 3
	AICc <- 2 * out$cum.neg.log.lik[1] + 2 * param * ((N - initial$Ts)/(N - initial$Ts - param - 1))
	out$N.tot <- N
	out$N.cond <- initial$Ts
	out$Param <- param
	out$AICc <- AICc
	out$Report <- fit
	out
}


#-----------------------------------------------------------------------------------
"kf.rw.Larkin" <- function(initial, X, y)
{
    # kf.rw adapted for Larkin model
    # Purpose:
    # ========
    # Finds ML estimates of kalman filter.  KF model assumes the following:
    # - Simple linear regression
    # - Time-varying intercept only
    # - Intercept follows a random walk
    # - Fits 6 PARAMETERS (b, b1, b2, b3, sig.e, sig.w)
	
    # Arguments:
    # ==========
    # initial  LIST with initial parameter values and starting conditions for estimation procedure.
    #    Names for the elements of initial must be:
    #
    #     mean.a  Initial value for mean of intercept in recursive calculations
    #     var.a   Initial value for variance of intercept in recursive calculations
    #     b      Starting value for ML estimation
    #     b1     Starting value for ML estimation
    #     b2     Starting value for ML estimation
    #     b3     Starting value for ML estimation
    #     ln.sig.e, ln.sig.w Starting values for natural logarithms of error terms in
    #          observation and system equations
    #     Ts     Number of observations at start of data set to omit for
    #          calculation of variance in observation equation and concentrated
    #          likelihood function.
    # X, y   Data for the observation equation
	
	
    # fit using nlminb
	fit <- nlminb(start=c(initial$b, initial$b1, initial$b2, initial$b3, initial$ln.sig.e, initial$ln.sig.w), objective=kalman.rw.fit.Larkin,
    gradient = NULL, hessian = NULL, scale = 1, control = list(), lower = -Inf, upper = Inf,
    initial$mean.a, initial$var.a, X, y, initial$Ts)
	
	if (fit$convergence != 0)
	{
		print(fit$message)
		warning("ML estimation for Larkin 'b' parameters and error terms failed to converge!")
	}
	
    # Perform recursive calculations for KF with ML estimates:
	out <- kalman.rw.Larkin(initial$mean.a, initial$var.a, fit$par[1], fit$par[2], fit$par[3], fit$par[4], fit$par[5], fit$par[6], X, y, initial$Ts)
	N <- length(X[, 1]) - sum(is.na(X[, 1]))
	param <- 6
	AICc <- 2 * out$cum.neg.log.lik[1] + 2 * param * ((N - initial$Ts)/(N - initial$Ts - param - 1))
	out$N.tot <- N
	out$N.cond <- initial$Ts
	out$Param <- param
	out$AICc <- AICc
	out$Report <- fit
	return(out)
}


#-------------------------------------------------------------------------------------------------------------
"kalman.rw.fit.Ricker" <- function(optim.vars, init.mean.a, init.var.a, x, y, Ts)
# a little helper function for ML fitting
# we need this in R because the optimizer functions are different from those in S
{
    # run Kalman filter and return cumulative log likelihood ...
	kalman.rw.Ricker(init.mean.a, init.var.a, optim.vars[1], optim.vars[2], optim.vars[3], x, y, Ts)$cum.neg.log.lik
}

#-------------------------------------------------------------------------------------------------------------
"kalman.rw.fit.Larkin" <- function(optim.vars, init.mean.a, init.var.a, X, y, Ts)
# a little helper function for ML fitting
# we need this in R because the optimizer functions are different from those in S
{
    # run Kalman filter and return cumulative log likelihood ...
	kalman.rw.Larkin(init.mean.a, init.var.a, optim.vars[1], optim.vars[2], optim.vars[3], optim.vars[4], optim.vars[5], optim.vars[6], X, y, Ts)$cum.neg.log.lik
}


#-----------------------------------------------------------------------------------
"kalman.rw.Ricker" <- function(init.mean.a, init.var.a, b, ln.sig.e, ln.sig.w, x, y, Ts = 0)
{
    #
    # Written for S-Plus 2000, Professional Release 2
    # Source of this code:
    # School of Resource and Environmental Management
    # Simon Fraser University
    # Burnaby, British Columbia
    # Canada V5A 1S6
    #
    # For future inquiries, contact:
    # Randall M. Peterman at the above School, either by
    #       e-mail at rmpasst@sfu.ca or by phone at 604 291-4683
    #
    # Authors of this code:
    # Brice MacGregor
    # February 14, 2001
    # Modified by Brian to handle missing data by using one-step-ahead forecasts, Oct. 10, 2002
	
    # Version Features:
    # - Simple linear regression
    # - Time-varying intercept only
    # - Intercept follows a random walk
	
    # Purpose:
    # ========
    # Performs recursive calculations for simple Kalman filter model with a time-varying
    # intercept parameter which follows a random walk.  Calculates the concentrated
    # likelihood function given the data, starting values for the mean and variance of
    # the intercept, and parameter values.  This function is used with the S-Plus function
    # "ms" to find the maximum likelihood estimates for the parameters.
    # Observation Equation:  y(t) = a(t) + b*x(t) + e(t)
    # System Equation:   a(t) = a(t-1) + w(t)
    # where:  v(t)~N(0,sig.e^2), w(t)~N(0,sig.w^2)
	
    # Arguments:
    # ========================
    # init.mean.a  Starting mean for intercept
    # init.var.a  Starting variance for intercept
    # b     Slope parameter
    # ln.sig.e   Natural log of the standard deviation of observation error*
    # ln.sig.w   Natural log of the standard deviation of system error*
    # x     Independent variable in obs. equation
    # y     Dependent variable in obs. equation
    # Ts    Number of years to omit when calculating the concentrated likelihood
    #     for the data set. See Visser and Molenaar (1988).  Default is zero.
    # *The natural log for the standard deviations of noise terms are used as inputs rather
    # than the straight standard deviations to ensure that the maximum likelihood procedure
    # only returns values that are greater than or equal to 1.  The function returns
    # the straight standard deviations in the output.
	
	
    # ignore years where no y value exists
	x[is.na(y)] <- NA
	
    # Calculate standard deviations for noise terms
	sig.e <- exp(ln.sig.e)
	sig.w <- exp(ln.sig.w)
	
    # Length of time series
	Tmax <- length(x)
	
    # Create vectors to store values calculated each year
	prior.mean.a <- rep(NA, Tmax)	# Prior mean of intercept (a)
	prior.var.a <- rep(NA, Tmax)	# Prior variance of intercept (a)
	y.hat <- rep(NA, Tmax)			# Predicted value of y(t) given y(t-1)
	f <- rep(NA, Tmax)				# Prediction variance
	v <- rep(NA, Tmax)				# Prediction error
	post.mean.a <- rep(NA, Tmax)	# Posterior mean of intercept (a)
	post.var.a <- rep(NA, Tmax)	# Posterior variance of intercept (a)
	filter.y <- rep(NA, Tmax)		# Filtered value for y
	neg.log.like <- rep(NA, Tmax)	# Negative log-likelihood - MIN to get ML estimates
	p.star <- rep(NA, Tmax)			# Used in smoothing
	smoothe.mean.a <- rep(NA, Tmax)	# Smoothed mean of intercept (a)
	smoothe.var.a <- rep(NA, Tmax)	# Smoothed variance of intercept (a)
	smoothe.y <- rep(NA, Tmax)		# Smoothed y
	
    # Start loop over time for recursive calculations:
	
    #Brian $$$$$$$
	cum.neg.log.lik <- 0
	for(t in 1:Tmax)
	{
        # Step 1: Calculate prior mean and variance of intercept (a)
        #   If t=1, then initial values are used as posteriors from previous period
        #   Else, posteriors from previous period are used
		
		if(t == 1)
		{
			prior.mean.a[t] <- init.mean.a
			prior.var.a[t] <- init.var.a
		}
		else
		{
			prior.mean.a[t] <- post.mean.a[t - 1]
			prior.var.a[t] <- post.var.a[t - 1] + sig.w^2
		}
		
		if(is.na(x[t]) == T)
		{
            # Step 2: Predict next value for a[t]
            #y.hat[t] <- prior.mean.a[t] + b * x[t]
            #v[t] <- y[t] - y.hat[t]
            #f[t] <- prior.var.a[t] + sig.e^2
			
            # Step 3: Generate posterior distribution for intercept (a):
			post.mean.a[t] <- prior.mean.a[t]
			post.var.a[t] <- prior.var.a[t]
            #filter.y[t] <- post.mean.a[t] + b * x[t]
			
            # Step 4: Calculate the concentrated likelihood function:
			neg.log.like[t] <- NA
			
		}
		
		else
		{
            # Step 2: Generate predicted value for y(t) given y(t-1) and error
			y.hat[t] <- prior.mean.a[t] + b * x[t]
			v[t] <- y[t] - y.hat[t]
			f[t] <- prior.var.a[t] + sig.e^2
			
            # Step 3: Generate posterior distribution for intercept (a):
			post.mean.a[t] <- prior.mean.a[t] + (prior.var.a[t] * (v[t]/f[t]))
			post.var.a[t] <- prior.var.a[t] - (prior.var.a[t]^2/f[t])
			filter.y[t] <- post.mean.a[t] + b * x[t]
			neg.log.like[t] <- (log(f[t]) + (v[t]^2/f[t]))/2
			
            #		cat("prior C: ", exp(- prior.mean.a[t]), " Ni: ", exp(x[t]), " b: ", b, " y.hat: ", y.hat[t], " y: ", y[t], " post C: ", exp(- post.mean.a[t]), "\n")
		}
		
	} # End loop over time
	
    # Step 5: Calculate cumulative value for concentrated negative log-likelihood
	cum.neg.log.lik <- sum(neg.log.like[(Ts+1):Tmax], na.rm=T)
    # Step 6: Smoothing of kalman filter estimates for time-varying intercept
    # Start loop over time (NB: Calculations start with last values first)
	for(t in Tmax:1)
	{
		if(t == Tmax)
		{
			p.star[t] <- NA
			smoothe.mean.a[t] <- post.mean.a[t]
			smoothe.var.a[t] <- post.var.a[t]
		}
		
		else
		{
			p.star[t] <- post.var.a[t]/prior.var.a[t + 1]
			smoothe.mean.a[t] <- post.mean.a[t] + p.star[t] * (smoothe.mean.a[t + 1] - prior.mean.a[t + 1])
			smoothe.var.a[t] <- post.var.a[t] + p.star[t]^2 * (smoothe.var.a[t + 1] - prior.var.a[t + 1])
		}
		
		smoothe.y[t] <- smoothe.mean.a[t] + b * x[t]
		
	} # End loop over time
	
    # Create a list to store output
    # =============================
	
    # Lines to put output in appropriate format
	init.mean.a <- as.vector(init.mean.a)
	init.var.a <- as.vector(init.var.a)
	b <- as.vector(b)
	sig.e <- as.vector(sig.e)
	sig.w <- as.vector(sig.w)
	out <- list(x = x, y = y, prior.mean.a = prior.mean.a,
    prior.var.a = prior.var.a, y.hat = y.hat, f = f, v = v,
    post.mean.a = post.mean.a, post.var.a = post.var.a,
    filter.y = filter.y, neg.log.like = neg.log.like,
    p.star = p.star, smoothe.mean.a = smoothe.mean.a,
    smoothe.var.a = smoothe.var.a, smoothe.y = smoothe.y,
    cum.neg.log.lik = cum.neg.log.lik, init.mean.a =
    init.mean.a, init.var.a = init.var.a, a.bar = NA, b = b,
    sig.e = sig.e, sig.w = sig.w, rho = NA)
	out
}


#-----------------------------------------------------------------------------------
"kalman.rw.Larkin" <- function(init.mean.a, init.var.a, b, b1, b2, b3, ln.sig.e, ln.sig.w, X, y, Ts = 0)
{
    # kalman.rw.Ricker adapted for Larkin model
    # Version Features:
    # - Simple linear regression
    # - Time-varying intercept only
    # - Intercept follows a random walk
	
    # Purpose:
    # ========
    # Performs recursive calculations for simple Kalman filter model with a time-varying
    # intercept parameter which follows a random walk.  Calculates the concentrated
    # likelihood function given the data, starting values for the mean and variance of
    # the intercept, and parameter values.  This function is used to find the maximum
    # likelihood estimates for the parameters.
    # Observation Equation:  y(t) = a(t) + b*X[, 1](t) + b1*X[, 2](t) + b2*X[, 3](t) + b3*X[, 4](t) + e(t)
    # System Equation:   a(t) = a(t-1) + w(t)
    # where:  v(t)~N(0,sig.e^2), w(t)~N(0,sig.w^2)
	
    # Arguments:
    # ========================
    # init.mean.a  Starting mean for intercept
    # init.var.a  Starting variance for intercept
    # b, b1, b2, b3  Static Larkin parameters
    # ln.sig.e   Natural log of the standard deviation of observation error*
    # ln.sig.w   Natural log of the standard deviation of system error*
    # x     Independent variable in obs. equation
    # y     Dependent variable in obs. equation
    # Ts    Number of years to omit when calculating the concentrated likelihood
    #     for the data set. See Visser and Molenaar (1988).  Default is zero.
    # *The natural log for the standard deviations of noise terms are used as inputs rather
    # than the straight standard deviations to ensure that the maximum likelihood procedure
    # only returns values that are greater than or equal to 1.  The function returns
    # the straight standard deviations in the output.
	
	
    # ignore years where no y value exists
	X[is.na(y), ] <- rep(NA, 4)
	
    # Calculate standard deviations for noise terms
	sig.e <- exp(ln.sig.e)
	sig.w <- exp(ln.sig.w)
	
    # Length of time series
	Tmax <- length(y)
	
    # Create vectors to store values calculated each year
	prior.mean.a <- rep(NA, Tmax)	# Prior mean of intercept (a)
	prior.var.a <- rep(NA, Tmax)	# Prior variance of intercept (a)
	y.hat <- rep(NA, Tmax)			# Predicted value of y(t) given y(t-1)
	f <- rep(NA, Tmax)				# Prediction variance
	v <- rep(NA, Tmax)				# Prediction error
	post.mean.a <- rep(NA, Tmax)	# Posterior mean of intercept (a)
	post.var.a <- rep(NA, Tmax)	# Posterior variance of intercept (a)
	filter.y <- rep(NA, Tmax)		# Filtered value for y
	neg.log.like <- rep(NA, Tmax)	# Negative log-likelihood - MIN to get ML estimates
	p.star <- rep(NA, Tmax)			# Used in smoothing
	smoothe.mean.a <- rep(NA, Tmax)	# Smoothed mean of intercept (a)
	smoothe.var.a <- rep(NA, Tmax)	# Smoothed variance of intercept (a)
	smoothe.y <- rep(NA, Tmax)		# Smoothed y
	
    # Start loop over time for recursive calculations:
	
    #Brian $$$$$$$
	cum.neg.log.lik <- 0
	for(t in 1:Tmax)
	{
        # Step 1: Calculate prior mean and variance of intercept (a)
        #   If t=1, then initial values are used as posteriors from previous period
        #   Else, posteriors from previous period are used
		
		if(t == 1)
		{
			prior.mean.a[t] <- init.mean.a
			prior.var.a[t] <- init.var.a
		}
		else
		{
			prior.mean.a[t] <- post.mean.a[t - 1]
			prior.var.a[t] <- post.var.a[t - 1] + sig.w^2
		}
		
		if(any(is.na(X[t, ])) == T)
		{
            # Skip step 2 and go straight to step 3: Generate posterior distribution for intercept (a):
			post.mean.a[t] <- prior.mean.a[t]
			post.var.a[t] <- prior.var.a[t]
			
            # Step 4: Calculate the concentrated likelihood function:
			neg.log.like[t] <- NA
		}
		
		else
		{
            # Step 2: Generate predicted value for y(t) given y(t-1) and error
			y.hat[t] <- prior.mean.a[t] + b * X[t, 1] + b1 * X[t, 2] + b2 * X[t, 3] + b3 * X[t, 4]
			v[t] <- y[t] - y.hat[t]
			f[t] <- prior.var.a[t] + sig.e^2
			
            # Step 3: Generate posterior distribution for intercept (a):
			post.mean.a[t] <- prior.mean.a[t] + (prior.var.a[t] * (v[t]/f[t]))
			post.var.a[t] <- prior.var.a[t] - (prior.var.a[t]^2/f[t])
			filter.y[t] <- post.mean.a[t] + b * X[t, 1] + b1 * X[t, 2] + b2 * X[t, 3] + b3 * X[t, 4]
			neg.log.like[t] <- (log(f[t]) + (v[t]^2/f[t]))/2
			
            #	cat("prior C: ", exp(- prior.mean.a[t]), " Ni: ", exp(x[t]), " b: ", b, " y.hat: ", y.hat[t], " y: ", y[t], " post C: ", exp(- post.mean.a[t]), "\n")
		}
		
	} # End loop over time
	
    # Step 5: Calculate cumulative value for concentrated negative log-likelihood
	cum.neg.log.lik <- sum(neg.log.like[(Ts+1):Tmax], na.rm=T)
    # Step 6: Smoothing of kalman filter estimates for time-varying intercept
    # Start loop over time (NB: Calculations start with last values first)
	for(t in Tmax:1)
	{
		if(t == Tmax)
		{
			p.star[t] <- NA
			smoothe.mean.a[t] <- post.mean.a[t]
			smoothe.var.a[t] <- post.var.a[t]
		}
		
		else
		{
			p.star[t] <- post.var.a[t]/prior.var.a[t + 1]
			smoothe.mean.a[t] <- post.mean.a[t] + p.star[t] * (smoothe.mean.a[t + 1] - prior.mean.a[t + 1])
			smoothe.var.a[t] <- post.var.a[t] + p.star[t]^2 * (smoothe.var.a[t + 1] - prior.var.a[t + 1])
		}
		
		smoothe.y[t] <- smoothe.mean.a[t] +  b * X[t, 1] + b1 * X[t, 2] + b2 * X[t, 3] + b3 * X[t, 4]
        
		
	} # End loop over time
	
    # Create a list to store output
    # =============================
	
    # Lines to put output in appropriate format
	init.mean.a <- as.vector(init.mean.a)
	init.var.a <- as.vector(init.var.a)
	b <- as.vector(b)
	b1 <- as.vector(b1)
	b2 <- as.vector(b2)
	b3 <- as.vector(b3)
	sig.e <- as.vector(sig.e)
	sig.w <- as.vector(sig.w)
	out <- list(X = X, y = y, prior.mean.a = prior.mean.a, 
    prior.var.a = prior.var.a, y.hat = y.hat, f = f, v = v, 
    post.mean.a = post.mean.a, post.var.a = post.var.a, 
    filter.y = filter.y, neg.log.like = neg.log.like, 
    p.star = p.star, smoothe.mean.a = smoothe.mean.a, 
    smoothe.var.a = smoothe.var.a, smoothe.y = smoothe.y, 
    cum.neg.log.lik = cum.neg.log.lik, init.mean.a = 
    init.mean.a, init.var.a = init.var.a, a.bar = NA, b = b, b1 = b1, b2 = b2, b3 = b3,
    sig.e = sig.e, sig.w = sig.w, rho = NA)
	out
}




