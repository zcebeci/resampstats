# Fonksiyon 3.1: Açýklayýcý istatistikleri hesaplama fonksiyonu 
calc.desc <- function(x, indices, gama=0.1){
   if(missing(indices)) indices <- 1:length(x)
   x <- x[indices]
   xmean <- mean(x)               #Aritmetik ortalama
   xmed <- median(x)              #Ortanca
   xtrmean <- mean(x, trim=gama)  #Budanmýþ ortalama
   xq1 <- unname(quantile(x))[2]  #Q1
   xq3 <- unname(quantile(x))[4]  #Q3
   xmin <- min(x)                 #Minimum
   xmax <- max(x)                 #Maksimum
   xiqr <- IQR(x)                 #Çeyreklikler arasý geniþlik
   xvar <- var(x)                 #Varyans
   xsd <- sd(x)                   #Std. sapma
   xse <- xsd/sqrt(length(x))     #Std. Hata
   xcv <- xsd / xmean * 100       #Varyasyon katsayýsý (%)
   descstats <- list(xmean, xmed, xtrmean, xmin, xmax, 
     xq1, xq3, xiqr, xvar, xsd, xse, xcv)
   names(descstats)<- c("mean", "med", "trmean", 
     "min", "max", "Q1", "Q3", "IQR", "var", "sd", "se", "cv")
   class(descstats) <- "desc"
   return(descstats)
}
