# Sgen function
get_Smsy <- function(a, b){
  Smsy <- (1 - lambert_W0(exp(1 - a))) / b
  return(Smsy)
}

get_Sgen <- function(a, b, int_lower, int_upper, Smsy){
  fun_Sgen <- function(Sgen, a, b, Smsy) {Sgen * a * exp(-b * Sgen) - Smsy}
  Sgen <- uniroot(fun_Sgen, interval=c(int_lower, int_upper), a=a, b=b, Smsy=Smsy)$root
  return(Sgen)
}

quantile_df <- function(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)){ #where probs is a vector of percentiles (0-100)
  tibble(val = quantile(x, probs, na.rm = TRUE))
}


# Several HCRs----------------------------------------------------------------------------
current_HCR <- function(R, OU){
  if(R <= 7.059){C <- ((R*(0.15/7.059))*R)*OU
  C <- ifelse(C<0, 0, C)
  S <- R-C}
  if(7.059 < R & R <= 20){C <- (R-6)*OU
  S <- R-C}
  if(R > 20){C <- (R*0.7)*OU
  C <- ifelse(C>R, R-.0001, C) #leave 100 spawners if catch OU makes C>R
  S <- R-C}

  U <- C/R

  return(c(S,C,U))
}

alt_HCR <- function(R, OU){
  if(R <= 6){C <- 0
  S <- R}
  if(R > 6 & R <= 15){C <- (R-6)*OU
  S <- R-C}
  if(R > 15){C <- (R*0.6)*OU
  C <- ifelse(C>R, R-.0001, C) #leave 100 spawners if catch OU makes C>R
  S <- R-C}

  U <- C/R

  return(c(S,C,U))
}

BB_HCR <- function(R, OU, Sgen, Smsy, Umsy){
  if(R <= Sgen){C <- 0
  S <- R} #might need to add protection here if forecast error makes R<0
  if(R > Sgen & R < Smsy){ C <- (R*(Umsy/Sgen)*R)
  S <- R-C}
  if(R > Smsy){C <- R*Umsy
  C <- ifelse(C>R, R-.0001, C)
  S <- R-C}

  U <- C/R

  return(c(S,C,U))
}