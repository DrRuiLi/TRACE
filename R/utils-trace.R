str_digit <- function(x, digit = 2) {
  sprintf(paste0("%.", digit, "f"), x)
}

num2percent <- function(x, digit = 2) {
  paste0(str_digit(x * 100, digit), "%")
}
