log2.scale <- function(x) {
  if (min(x) > 0) {
    
  }
  else {
    mindata.norm = abs(min(x)) + .001
    data.norm = x + mindata.norm
  }
  data.log <- t(apply(data.norm, 1, log2))
}
