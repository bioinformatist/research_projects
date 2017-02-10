NormalizeConstantAsCol <- function(col, scaling.constant = 500) {
  scaling.factor <- scaling.constant/mean(col, trim = 0.02)
  return(scaling.factor * col)
}

NormalizeConstant <- function(eData, scaling.constant = 500) {
  # Normalize the expression matrix. Only used for data.frame objects.
  #
  # Args:
  #   eData: Matrix that contains expression intensity.
  #   sc: A constant as reference of normalization. Default: 500.
  #
  # Returns:
  #   The normalized expression matrix.
  return(apply(eData, 2, NormalizeConstantAsCol, scaling.constant = scaling.constant
  ))
}
