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

# Fonksiyon 3.2: Açýklayýcý istatistikleri görüntüleme fonksiyonu 
print.desc <- function(x){
    cat("Ortalama : ", x$mean, "\n")
    cat("Ortanca  : ", x$med, "\n")
    cat("Bud.Ort. : ", x$trmean, "\n")
    cat("Min      : ", x$min, "\n")
    cat("Max      : ", x$max, "\n")
    cat("Q1       : ", x$Q1, "\n")
    cat("Q3       : ", x$Q3, "\n")
    cat("IQR      : ", x$IQR, "\n")
    cat("Varyans  : ", x$var, "\n")
    cat("Std.sapma: ", x$sd, "\n")
    cat("Std.Hata : ", x$se, "\n")
    cat("Var.Kat.%: ", x$cv, "\n")
 }

# Fonksiyon 3.3: Açýklayýcý istatistik grafikleri
# Baðýmlýlýk – Paket: vioplot
univar.plot <- function(x){
    par(mfrow=c(3,2))
    hist(x, col="gray90", prob=TRUE, ylab="Yoðunluk",
      main="Histogram")
    lines(density(x), col=4)
    if(!require(vioplot, quietly=TRUE)){ 
       boxplot(x, col="gray90", bg=3, xlab="x",
         main="Kutu-býyýk grafiði")
    }else{
       vioplot::vioplot(x, col="gray90", bg=3,
         xlab="", main="Keman grafiði")
    }
     plot(x, col=1, bg="gray90", xlab="", main="Serpilme grafiði")
    abline(h=median(x), col=2)
    abline(h=mean(x), col=4)
    stripchart(x, col=1, bg="gray90", pch=21, cex=1, 
      xlab="x", main="Þerit grafik")
    plot.ecdf(x, col=1, bg="gray90", main="ECDF")
    curve(pnorm(x,mean(x), sd(x)), col=4, add=TRUE)
    qqnorm(x, datax=TRUE, pch=21, col=1, bg="gray90",
       main="Normal Q-Q grafiði")
    qqline(x, datax=TRUE, col=4,lwd=1) # Q-Q doðrusu
 }



