H_md.Y <- function(freq,sex,age) {
  a <- subset(alpha, Sex == sex & Freq == freq, select = c('Alpha'))
  as.numeric(a * (age - 18)^2)
}

Z_score <- function(threshold, freq, sex, age){
  theH <- H_md.Y(freq, sex, age)
  if (threshold > theH){
    theB <- subset(bubl, Freq == freq & Sex == sex, select=c('Bu'))
    theS <- theB + 0.445 * theH
    theZ <- (threshold - theH) / theS
  } else {
    theB <- subset(bubl, Freq == freq & Sex == sex, select=c('Bl'))
    theS <- theB + 0.356 * theH
    theZ <- (threshold - theH) / theS
  }
  as.numeric(theZ)
}

Z248 <- function(T2, T4, T8, sex, age){
  Z2 <- Z_score(T2, 2000, sex, age)
  Z4 <- Z_score(T4, 4000, sex, age)
  Z8 <- Z_score(T8, 8000, sex, age)
  the248 <- c(Z2,Z4,Z8)
  c(mean(the248),sd(the248))
}
