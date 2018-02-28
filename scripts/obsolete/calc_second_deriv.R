# Function from A. Raiho to do the sum of squares of the second deriviative (our stability index)

calc.second.deriv <- function(dataset,h=1, H=1){
  # biomassCI = dataset that you want to calculate the second deriv on
  # h = timestep for the derivative; default = 1
  # H = real time step (e.g. annual, 100-yr bins)
  
  T.tot <- length(dataset) - h
  
  t1 <- h + 1
  
  second.deriv <- sum(((dataset[(t1:T.tot)+h]-2*dataset[(t1:T.tot)]+dataset[(t1:T.tot)-h])/((h*H)^2))^2)
  
  return(second.deriv)
  
}


