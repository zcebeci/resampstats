#### BÖLÜM 2 ###########################################################

# Fonksiyon 2.1: Ortalama hesaplama fonksiyonu
ortalama <- function(x, yontem) {
   if(yontem==1){
      xort <- mean(x)
   }else if(yontem==2){
      xort <- median(x)
   }else if(yontem==3){
      xort <- mean(x, trim=0.1)
   }else{
     xort <- NA
     stop("Yanlýþ yöntem seçimi")
     }
     return(xort)
}

#### BÖLÜM 3 ###########################################################
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

# BÖLÜM 4 ####################################################################

# Fonksiyon 4.1: Ortalama hesaplama fonksiyonu
 calc.mean <- function(x, indices){
   if(is.vector(x)) x <- as.data.frame(x)
   if(missing(indices)) indices <- 1:nrow(x)
   xmean <- mean(x[indices,1])
   return(xmean)
 }

# Fonksiyon 4.2a: Parametrik olmayan bootstrap fonksiyonu 1
 bs.nonpar1 <- function(dset, statistic, R=2000){
   if(!is.data.frame(dset) | !is.matrix(dset)) 
     dset <- as.data.frame(dset)
   n <- nrow(dset)  # Örneklem büyüklüðü
   # Orijinal örneklem kestirimleri
   tetahat <- statistic(dset, indices=1:n) 
   nstats <- length(tetahat) # Ýstatistik sayýsý
   # Bootstrap kestirimleri matrisi
   tetastar <- matrix(NA, nrow=R, ncol=nstats)
   for (i in 1:R){
     indices <- sample(1:n, replace=TRUE)
     tetastar[i,] <- statistic(dset,indices=indices)
   }
   return(list(t0=tetahat, t=tetastar, data=dset,
     statistic=statistic, R=R, call=match.call()))
 }

# Fonksiyon 4.2b: Parametrik olmayan bootstrap fonksiyonu 2
 bs.nonpar2 <- function(dset, statistic, R=2000){
   if(!is.data.frame(dset) | !is.matrix(dset)) 
     dset <- as.data.frame(dset)
   n <- nrow(dset)  # Örneklem büyüklüðü
   tetahat <- statistic(dset)   # Orijinal örneklem kestirimleri
   nstats <- length(tetahat)    # Ýstatistik sayýsý
   # Bootstrap yeniden örneklem indisleri
   indices <- sample(1:n, n, replace=TRUE)
   # Bootstrap yeniden örneklemler matrisi
   bssamps <- array(dset[indices,], c(R, n))
   # Bootstrap kestirimleri
   tetastar <- matrix(apply(bssamps, 1, statistic), 
     ncol=nstats, nrow=R, byrow=TRUE)
   return(list(t0=tetahat, t=tetastar, data=dset,
     statistic=statistic, R=R, call=match.call()))
 }

# Fonksiyon 4.3a : Çifte bootstrap fonksiyonu 1
 dbs.nonpar1 <- function(dset, statistic, double=FALSE, 
   R=2000, R2=NULL){
   if(!is.data.frame(dset) | !is.matrix(dset)) 
     dset <- as.data.frame(dset)
   n <- nrow(dset)  # Örneklem büyüklüðü
 # Orijinal örneklem kestirimleri
   tetahat <- statistic(dset, indices=1:n) 
   nstats <- length(tetahat) # Ýstatistik sayýsý
 # Bootstrap kestirimleri matrisi
   tetastar <- matrix(NA, nrow=R, ncol=nstats)  
 # Bootstrap t deðerleri matrisi
   tvals <- matrix(NA, nrow=R, ncol=nstats)  
   for (i in 1:R){
     indices <- sample(1:n, replace=TRUE)
     tetastar[i,] <- statistic(dset,indices=indices)
     if(double){
     ## Ýç bootstrap baþý
        if(is.null(R2) | !is.numeric(R2)) R2 <- floor(R*0.1)
        dset2 <- dset[indices,]
      #Ýç bootstrap kestirimleri 
        tetastar2 <- matrix(NA, nrow=R2, ncol=nstats)  
        for(j in 1:R2){
          indices <- sample(1:n, replace=TRUE)
          tetastar2[j,] <- statistic(dset2,indices=indices)
        }
        tvals[i,] <- (tetastar[i,]-tetahat)/
          apply(tetastar2, 2, sd)
        # tvals[which(tvals==Inf | tvals==-Inf)] <- 0
     } else{tvals <- NULL}
     ## Ýç bootstrap sonu
   }
   return(list(t0=tetahat, t=tetastar, tvals=tvals, data=dset, 
    statistic=statistic, double=double, R=R, R2=R2,
     call=match.call()))
 }

# Fonksiyon 4.3b : Çifte bootstrap fonksiyonu 2
 dbs.nonpar2 <- function(dset, statistic, double=FALSE,
    R=2000, R2=NULL){
 ## Ýç bootstrap fonksiyonu 
    bs.inner <- function(dset, statistic, R=200){
       if(!is.data.frame(dset) | !is.matrix(dset)) 
          dset <- as.data.frame(dset)
       n <- nrow(dset)  # Örneklem büyüklüðü
   # Orijinal örneklem kestirimleri
       tetahat <- statistic(dset, indices=1:n) 
       nstats <- length(tetahat) # Ýstatistik sayýsý
       tetastar2 <- matrix(NA, nrow=R, ncol=nstats)  
   # Ýç bootstrap kestirimleri 
       for(j in 1:R){
          indices <- sample(1:n, replace=TRUE)
          tetastar2[j,] <- statistic(dset,indices=indices)
       }
       se <- apply(tetastar2, 2, sd)
       return(se)
     }
 ## Dýþ bootstrap 
     if(!is.data.frame(dset) | !is.matrix(dset)) 
       dset <- as.data.frame(dset)
     n <- nrow(dset)  # Örneklem büyüklüðü
    # Orijinal örneklem kestirimleri
     tetahat <- statistic(dset, indices=1:n) 
     nstats <- length(tetahat) # Ýstatistik sayýsý
    # Bootstrap kestirimleri matrisi
     tetastar <- matrix(NA, nrow=R, ncol=nstats)  
    # Bootstrap t deðerleri matrisi
     tvals <- matrix(NA, nrow=R, ncol=nstats)  
     for(i in 1:R){
        indices <- sample(1:n, replace=TRUE)
        tetastar[i,] <- statistic(dset,indices=indices)
        if(double){
           if(is.null(R2) | !is.numeric(R2)) 
              R2 <- floor(R*0.1)
           se <- bs.inner(dset[indices,], statistic, R=R2)
           tvals[i,] <- (tetastar[i,]-tetahat)/se
           # tvals[which(tvals==Inf | tvals==-Inf)] <- 0
        }else{
           tvals <- NULL
           R2 <- NULL
        }
     }
     return(list(t0=tetahat, t=tetastar, tvals=tvals, data=dset,
       statistic=statistic, double=double, R=R, R2=R2,
       call=match.call()))
 }

# Fonksiyon 4.4: Ortanca hesaplama fonksiyonu
 calc.median <- function(x, indices){
   if(is.vector(x)) x <- as.data.frame(x)
   if(missing(indices)) indices <- 1:nrow(x)
   xmed <- median(x[indices,1])
   return(xmed)
 }

# Fonksiyon 4.5: Budanmýþ ortalama hesaplama fonksiyonu
 calc.trmean <- function(x, indices, gama=0.1){
   if(is.vector(x)) x <- as.data.frame(x)
   if(missing(indices)) indices <- 1:nrow(x)
   xtrmean <- mean(x[indices,1], trim=gama)
   return(xtrmean)
 }

# Fonksiyon 4.6: Açýklayýcý istatistikler fonksiyonu 
 calc.stats <- function(x, indices, gama=0.1){
   if(is.vector(x)) x <- as.data.frame(x)
   if(missing(indices)) indices <- 1:nrow(x)
   x <- x[indices,]
   xmean <- mean(x)               #Aritmetik ortalama
   xmed <- median(x)              #Ortanca
   xtrmean <- mean(x, trim=gama)  #Budanmýþ ortalama
   xmin <- min(x)                 #Minimum
   xmax <- max(x)                 #Maksimum
   xiqr <- IQR(x)                 #Çeyreklikler arasý geniþlik
   xvar <- var(x)                 #Varyans
   xsd <- sd(x)                   #Std. sapma
   xse <- xsd/sqrt(length(x))     #Std. Hata
   xcv <- xsd / xmean * 100       #Varyasyon katsayýsý (%)
   sonuc <- c(xmean, xmed, xtrmean, xmin, xmax, xiqr, xvar, 
      xsd, xse, xcv)
   return(sonuc)
 }

# Fonksiyon 4.7: Bootstrap sonuçlarýný özetleme fonksiyonu
 bs.sonuc <- function(bootres, index=1){
   est <- bootres$t0[index]  # Ýstatistik kestirimi
   sbst <- sort(bootres$t[,index]) #Sýralý bootstrap kestirimleri 
   if(!is.null(bootres$tvals)) 
     tvals <- bootres$tvals[,index] else tvals <- NULL
   bsest <- mean(sbst) #Bootstrap kestirimleri ortalamasý
   se <- sd(sbst)  # Standart hata kestirimi
   bias <- bsest-est # Yanlýlýk kestirimi
   mse <- mean((sbst-est)^2)  # Hata kareler ortalamasý (MSE)
   R <- length(sbst)
   return(list(est=est, bsest=bsest, bias=bias, se=se, mse=mse,  
     sbst=sbst, tvals=tvals, R=R))
 }

# Fonksiyon 4.8: Normal GA (normal CI)
# Baðýmlýlýk – Fonksiyon: 4.7
#
 ga.norm <- function(bootres, index=1, conf.level=0.95){
   alpha <- 1-conf.level
   bsonuc <- bs.sonuc(bootres, index)
   Z <- qnorm(1-alpha/2)
   c(bsonuc$est-Z*bsonuc$se, bsonuc$est+Z*bsonuc$se)
 }

# Fonksiyon 4.9: Basit GA (basic CI) 
# Baðýmlýlýk – Fonksiyon 4.7; Paket: boot
#
 ga.temel <- function(bootres, index=1, conf.level=0.95){
   if(class(bootres)=="boot"){
      gares <- boot.ci(bootres, conf=conf.level, type="basic",
        index=index)
      ga <- gares$basic[4:5]
   }else{  
   alpha <- 1-conf.level
   bsonuc <- bs.sonuc(bootres, index)
   ga <- c(2*bsonuc$est-quantile(bsonuc$sbst, 
       prob=c(1-alpha/2, alpha/2), names=FALSE))
   }
   return(ga)
 }

# Fonksiyon 4.10: Yüzdelik GA (percentile CI)
# Baðýmlýlýk – Fonksiyon 4.7
#
 ga.yuzde <- function(bootres, index=1, conf.level=0.95){
   if(class(bootres)=="boot"){
      gares <- boot.ci(bootres, conf=conf.level, type="perc",
        index=index)
      ga <- gares$perc[4:5]
   }else{
      alpha <- 1-conf.level
      bsonuc <- bs.sonuc(bootres, index)
      yuzdelik <- c(alpha/2, 1-alpha/2)
      ga <- c(quantile(bsonuc$sbst, prob=yuzdelik, names=FALSE))
   }
   return(ga)
 }

# Fonksiyon 4.11: Yanlýlýk düzeltmeli GA (BC: Bias-corrected CI)
# Baðýmlýlýk – Fonksiyon 4.7
#
 ga.yd <- function(bootres, index=1, conf.level=0.95){
   alpha <- 1-conf.level
   bsonuc <- bs.sonuc(bootres, index)
  # Yanlýlýk düz. fak.
   Z0 <- qnorm((sum(bsonuc$sbst <= bsonuc$est)/2)/bsonuc$R) 
  # Standart normal sýnýrlarý sapta
   Z <- qnorm(c(alpha/2, 1-alpha/2)) 
  # Yanlýlýk düzeltmeli yuzdelikleri hesapla
   yuzdelik <- pnorm(Z-2*Z0) 
  # Yanlýlýk düzeltmeli güven aralýðýný hesapla
   return(quantile(bsonuc$sbst, prob=yuzdelik, names=FALSE)) 
 }

# Fonksiyon 4.12: Hýzlandýrýlmýþ-Yanlýlýk düzeltmeli GA (BCa)
# Baðýmlýlýk – Fonksiyon 4.7; Paket: boot
#
 ga.ydh <- function(bootres, index=1, conf.level=0.95){
   alpha <- 1-conf.level
   if(class(bootres)=="boot"){
      gares <- boot.ci(bootres, conf=conf.level, type="bca",
         index=index)
      ga <- gares$bca[4:5]
   } else{ 
      bsonuc <- bs.sonuc(bootres, index)
     # Std. norm sapma
      bias <- qnorm((sum(bsonuc$sbst >= bsonuc$est)/2)/bsonuc$R) 
     # Hýzlandýrma (a) sabitini hesapla
      n <- nrow(bootres$data)
      v <- n-1 
      nt <- n*bsonuc$sbst
      pseudoval <- c()
      for(i in 1:n)
        pseudoval[i] <- nt[i]-v*bootres$statistic(
          bootres$data[-i,])[index]
      jackest <- mean(pseudoval)-pseudoval
      a <- sum(jackest^3)/(6*sum(jackest^2))^(3/2)
     # Standart normal limitleri hesaplama
      Z <- qnorm(c(alpha/2, 1-alpha/2)) 
     # Yanlýlýðý düzeltme ve oranlara dönüþtürme 
      yuzdelik <- pnorm((Z-bias)/(1-a*(Z-bias))-bias)
     # BCa yüzdelik limitleri hesaplama 
      ga <- quantile(bsonuc$sbst, prob=yuzdelik, names=FALSE) 
   }
   return(ga)
 }

# Fonksiyon 4.13: Studentlaþtýrýlmýþ (Bootstrap-t) GA hesaplama 
# Baðýmlýlýk – Fonksiyon 4.7; Paket: boot
 ga.stud <- function(bootres, index=1, sigteta=NULL,
   conf.level=0.95){
   if(class(bootres)=="boot"){
      if(!is.null(sigteta)){
        gares <- boot.ci(bootres, conf=conf.level, type="stud", 
        var.t0=sigteta^2, index=index)
      }else{
        gares <- boot.ci(bootres, conf=conf.level, type="stud", 
        var.t=bootres$t[,index+1])
      }
      ga <- gares$stud[4:5]
   }else{ 
      alpha <- 1-conf.level
      bsonuc <- bs.sonuc(bootres, index)
      if(!is.null(bsonuc$tvals)){
         stvals <- sort(bsonuc$tvals)
         alt <- floor(length(stvals) * (1-(alpha/2)))
         ust <- floor(length(stvals) * alpha/2)
         ga <- c(bsonuc$est-stvals[alt]*bsonuc$se,
        bsonuc$est-stvals[ust]*bsonuc$se)
      }else if(!is.null(sigteta)){
         se <- sigteta/sqrt(bsonuc$R)
         ga <- c(bsonuc$est-qt(df=bsonuc$R-1, c(1-alpha/2,
           alpha/2))*se)
      }else{
        cat("Studentlaþtýrýlmýþ GA için çifte bootstrap nesnesi
          veya sigteta argümanýnda standart sapma kestirimi
          gerekir\n")
        ga <- c(NA, NA)
      }
    }
    return(ga)
 }

# Fonksiyon 4.14: Bootstrap güven aralýklarý fonksiyonu  
# Baðýmlýlýk – Fonksiyon: 4.8, 4.9, 4.10, 4.11, 4.12, 4.13
#
 ga.bs <- function(bootres, index=1, conf.level=0.95){
   norm.ga <- ga.norm(bootres,index, conf.level)
   basit.ga <- ga.temel(bootres, index, conf.level)
   yuzde.ga <- ga.yuzde(bootres, index, conf.level)
   yd.ga <- ga.yd(bootres, index, conf.level)
   ydh.ga <- ga.ydh(bootres, index, conf.level)
   stud.ga <- ga.stud(bootres, index, conf.level)
   ga <- cbind(norm.ga, basit.ga, yuzde.ga, yd.ga, 
      ydh.ga, stud.ga)
   gamat <- matrix(ga, nrow=dim(ga)[1], ncol=dim(ga)[2])
   colnames(gamat) <- c("gs.norm", "ga.temel", 
     "ga.yuzde", "ga.ydh", "ga.ydh", "ga.stud")
   rownames(gamat) <- c("GS.alt", "GS.üst")
   (gamat)
 }

# Fonksiyon 4.15: Bootstrap kestirimleri daðýlýþ grafiði
# Baðýmlýlýk – Fonksiyon: 4.7, 4.8, 4.9, 4.10, 4.11, 4.12, 4.13
 bs.dist <- function(bootres, index=1, gayontem=NULL,
   sigteta=NULL, conf.level=0.95){
   bsonuc <- bs.sonuc(bootres, index)
   cnum <- 2*(IQR(bsonuc$sbst)/length(bsonuc$sbst)^(1/3)) #Histogram sýnýf sayýsý
   ymax <- max(density(bsonuc$sbst)$y)
   hist(bsonuc$sbst, prob=TRUE, w=cnum, border=NA, col="gray", 
     ylim=c(0, 1.2*ymax), ylab="Yoðunluk", xlab="t*",
     main=paste("Bootstrap=",formalArgs(bs.dist)[1], ", index=",
     index), sub=paste0("GA Yöntemi: ", gayontem))
   lines(density(bsonuc$sbst), lwd=2, col=1)
   segments(bsonuc$est,0, bsonuc$est, 1.1*ymax, col=2)
   segments(bsonuc$bsest,0, bsonuc$bsest,1.05*ymax, col=4)
   text(bsonuc$est, 1.15*ymax, "t0")
   text(bsonuc$bsest, 1.1*ymax, "t*")
   if(!is.null(gayontem)){
    if(gayontem =="yuzde")
      ga <- ga.yuzde(bootres, index, conf.level)
    else if(gayontem =="temel")
      ga <- ga.temel(bootres, index, conf.level)
    else if(gayontem =="norm")
      ga <- ga.norm(bootres, index, conf.level)
    else if(gayontem =="yd")
      ga <- ga.yd(bootres, index, conf.level)
    else if(gayontem =="ydh")
      ga <- ga.ydh(bootres, index, conf.level)
    else if(gayontem =="stud")
      ga <- ga.stud(bootres, index=index, 
         sigteta=sigteta, conf.level)
    segments(ga[1],0, ga[1], 1.1*ymax, col=3, lwd=2)
    segments(ga[2],0, ga[2], 1.1*ymax, col=3, lwd=2)
    text(ga[1], 1.15*ymax, srt=0, "GS.alt")
    text(ga[2], 1.15*ymax, srt=0, "GS.üst")
   }
 }

# Fonksiyon 4.16: Bootstrap diyagnostik grafikleri
# Baðýmlýlýk – Fonksiyon: 4.7; Paket: vioplot
 bs.plot <- function(bootres, index=1, gayontem="yuzde",
   conf.level=0.95){
   opar <- par(mfrow=c(3,2))
   if(!require(vioplot, quietly=TRUE)){
     install.packages("vioplot"); 
     require(vioplot, quietly=TRUE)}
   bsonuc <- bs.sonuc(bootres, index)
 # Grafik 1: Bootstrap kestirimleri (t*) kutu-býyýk grafiði
   boxplot(bsonuc$sbst, col=3, bg=3, pch=21,
    ylab="t*", main="t* Kutu-býyýk Grafiði")
 # Grafik 2: Bootstrap kestirimleri (t*) keman grafiði
   vioplot::vioplot(bsonuc$sbst, col="gray90", 
     rectCol=3, colMed=1, xaxt="n", ylab="t*",
     main="t* Keman Grafiði", names=c(" "))
 # Grafik 3: Bootstrap kestirimleri (t*) QQ grafiði
   qqnorm(bsonuc$sbst, col=3, bg="gray90", pch=21, cex=1.1, 
     main="t* Norm QQ Grafiði", xlab="Teorik Kantiller",
      ylab="Örneklem Kantilleri")
   qqline(bsonuc$sbst, col=1, lwd=2) 
 # Grafik 4: Bootstrap kestirimleri (t*) ECDF grafiði
   plot(ecdf(bsonuc$sbst), lwd=1, xlab="t*", ylab="Fn(t*)", 
     main="t* Ampirik CDF")
 # Grafik 5: Bootstrap kestirimleri (t*) simetri grafiði
   plot(bsonuc$sbst, type="l", lwd=1, col=1, 
    xlab="Sýra#", ylab="t*", 
    main="t* kestirimleri simetri grafiði")
   abline(h= bsonuc$est, col=2, lty=2)
   abline(h= bsonuc$bsest, col=4, lty=3)
 # Grafik 6: Bootstrap kestirimleri (t) R'ye göre deðiþim grafiði
   plot(bootres$t[,index], type="l", lwd=1, col=3, 
    xlab="R", ylab="t*", main="Tekrarlara göre t* kestirimleri")
   abline(h= bsonuc$est, col=2, lty=2)
   abline(h= bsonuc$bsest, col=4, lty=3)
   par(opar)
 }

# Fonksiyon 4.17a: Çifte bootstrap için istatistik (ortalama) 
# hesaplama fonksiyonu
# Baðýmlýlýk – Fonksiyon: 4.18, Paket: boot
#
 calc.dbmean <- function(x, indices){
   if(is.vector(x)) x <- as.data.frame(x)
   if(missing(indices)) indices <- 1:nrow(x)
   xmean <- mean(x[indices,])
 #Ýç bootstrap çaðrýsý
   icboot <- boot::boot(x, statistic=calc.var, R=300) 
   xvar <- var(icboot$t)
   return(c(xmean, xvar))
 }

# Fonksiyon 4.17b: Çifte bootstrap için varyans hesaplama
# Baðýmlýlýk – Paket: boot
#
 calc.var <- function(x, indices){
   if(is.vector(x)) x <- as.data.frame(x)
   if(missing(indices)) indices <- 1:nrow(x)
   xvar <- var(x[indices,]) 
   return(xvar)
 }

# Fonksiyon 4.18a: Güven aralýklarý matrisi fonksiyonu
# Baðýmlýlýk – Paket: boot
# 
 ci.matrix <- function(cires){
   cm <- matrix(NA, nrow=5, ncol=2)
   rownames(cm) <- c("Yüzdelik", "Temel", "Normal", 
     "BCa", "Stud")
   colnames(cm) <- c("GS.Alt","GS.Üst")
   if(class(cires)=="bootci"){
     cm[1,] <- c(cires$percent[4],cires$percent[5])
     cm[2,] <- c(cires$basic[4],cires$basic[5])
     cm[3,] <- c(cires$normal[2], cires$normal[3])
     cm[4,] <- c(cires$bca[4], cires$bca[5])
     if(!is.null(cires$stud[4])) 
       cm[5,] <- c(cires$stud[4], cires$stud[5])
     else 
       cm <- cm[-5,]
   }
   return(list(cm=cm, t0=cires$t0))
 }

# Fonksiyon 4.18b: Güven aralýklarý grafiði fonksiyonu
 ci.plot <- function(cimat){
   rpal <- c(1, 2, 3, 4, 6)
   cm <- cimat$cm
   plot(NA,NA, xlim=c(min(cm)-0.2*min(cm), max(cm)+0.2*max(cm)),
     ylim=c(1, nrow(cm)+0.5), xlab="Aralýk", ylab="", 
     yaxt="n", bty="n", main="Güven Aralýklarý")
   axis(1, labels=FALSE)
   axis(2, at=1:nrow(cm), labels=rownames(cm), 
     tck=1, lwd=0, lty=2)
   for(i in 1:nrow(cm)){
     segments(cm[i,1], i, cm[i,2], i,  lwd=2, col=rpal[i])
     text(cimat$t0, nrow(cm)+0.3, expression(hat(theta)))
     points(cm[i,1], i, col=rpal[i], cex=0.8, pch="[")
     points(cm[i,2], i, col=rpal[i], cex=0.8, pch="]")
   }
   segments(cimat$t0, 0, cimat$t0, nrow(cm)+0.2,  
    lwd=2, col="gray", lty=2)
 }

# Fonksiyon 4.18c: Farklý örneklem büyüklüklerinde GA deðiþimi
# Baðýmlýlýk – Fonksiyon: 4.18a; Paket: boot
#
 cichange <- function(dset, ssize, statistic, 
    iade=TRUE, conf.level=0.95, R=2000){
  if(is.vector(dset)) dset <- as.data.frame(dset)
  n <- dim(dset)[1]
  if(missing(ssize)) ssize <- c(5, 15, 30, 50)
  xcimat <- data.frame()
  for(i in ssize){
   sdset <- dset[sample(1:n, i, replace=iade),]
   bootres <- boot(data=sdset, statistic, R=R)
   ci.bootres <- boot.ci(bootres, conf=conf.level)
   rn <- nrow(ci.matrix(ci.bootres)$cm)
   xcimat <- rbind(xcimat, cbind(rep(i,rn),
     rep(unname(ci.bootres$t0),rn),
     unname(ci.matrix(ci.bootres)$cm)))
  }
  xcimat <- cbind(rep(1:rn,length(ssize)),xcimat)
  colnames(xcimat) <- c("GA.Ynt.","n", "t0", "GS.Alt", "GS.Üst")
  (xcimat)
 }


# Fonksiyon 4.18d: GA deðiþim grafiði fonksiyonu
 cichange.plot <- function(x, idx=0, rpal=c(1,2,3,4,6)){
   n <- max(x[,2])
   ganames <- c("Yüzdelik","Temel", "Normal", "BCa","Stud")
   plot(NA, NA, xlim=c(0,n+5), ylim=c(min(x[,4]), max(x[,5])), 
   xlab="", ylab="Güven Aralýðý", main="Güven Aralýðý Deðiþimi")
   if(idx>=1 & idx<=5) x <- x[which(x[,1]==idx),]
   j <- 0; t <- ceiling(n/100)
   for(i in 1:nrow(x)){
     segments(x[i,2]+j, x[i,4], x[i,2]+j, x[i,5],
     col=rpal[x[i,1]], lwd=2)
     points(x[i,2]+1.5*t, x[i,3], col=1, pch=19, cex=2)
     j <- j+t; if(j>3*t) j <- 0
   }
   lines(x[,2]+1.5*t, x[,3], col="gray", lty=3, lwd=3)
   legend("top", legend=ganames, bty="n", horiz=TRUE, col=rpal,
    lwd=2, lty=1)
 }

# Fonksiyon 4.18e: Farklý R tekrarlarýna göre GA deðiþimleri
# Baðýmlýlýk –Fonksiyon: 4.18a; Paket: boot
#
 cichangeR <- function(dset, ssize, statistic, 
    iade=TRUE, conf.level=0.95, R=2000){
  if(is.vector(dset)) dset <- as.data.frame(dset)
  n <- dim(dset)[1]
  if(missing(ssize)) ssize <- 30
  sdset <- dset[sample(1:n, ssize, replace=iade),]
  xcimat <- data.frame()
  for(i in R){
   bootres <- boot(data=sdset, statistic, R=i)
   ci.bootres <- boot.ci(bootres, conf=conf.level)
   rn <- nrow(ci.matrix(ci.bootres)$cm)
   xcimat <- rbind(xcimat, cbind(rep(i,rn),
     rep(unname(ci.bootres$t0),rn),
     unname(ci.matrix(ci.bootres)$cm)))
  }
  xcimat <- cbind(rep(1:rn,length(ssize)),xcimat)
  colnames(xcimat) <- c("GA.Ynt.","n", "t0", "GS.Alt", "GS.Ust")
  (xcimat)
 }

# Fonksiyon 4.19: Farklý örneklem büyüklükleri için bootstrap 
# performans analizi
# Baðýmlýlýk – Fonksiyon: 4.3a, 4.20
#
 bs.perform1 <- function(x, teta, statistic,  
   gayontem=ga.yuzde, m=1000, n=c(30), dbs=FALSE, 
     R=2000, conf.level=0.95){
   if(is.vector(x)) x <- as.data.frame(x)
  # Güven aralýklarý 
   gvnara <- array(0, dim=c(m, 2, length(teta))) 
  # Baþarým ölçüleri dizisi
   perfarr <- array(0, dim=c(length(n), 6, length(teta))) 
   dimnames(perfarr) <- list(n,c("n", "KO(%)",
     "AKO(%)","ÜKO(%)","OG","BD"), names(teta))
   for(i in 1:length(n)){
     for(j in 1:m){
       sidx <- sample(1:nrow(x), size=n[i], replace=TRUE)
       bsres <- dbs.nonpar1(x[sidx,], statistic=statistic,
         double=dbs, R=R)
       for(k in 1:length(teta)){
         gvnara[j,,k] <- gayontem(bsres, index=k,
           conf.level=conf.level)
         cat("Ýlerleme durumu: n=",n[i], "Örneklem=",j, "\n")
       }
     }
     for(k in 1:length(teta)){
       perfarr[i,,k] <- c(n[i], kapsama(ga=gvnara[,,k],
         teta=teta[k], conf.level))
     }
   }
   (perfarr)
 }


# Fonksiyon 4.20: Baþarým ölçüleri hesaplama fonksiyonu
 kapsama <- function(ga, teta, conf.level=0.95){
  nrga <- nrow(ga)
 # Kapsama oraný
  ko <- sum(ga[,1] <= teta & ga[,2] >= teta)/nrga * 100 
  ako <- sum(ga[,2] < teta)/nrga * 100  # teta<GA.alt oraný
  uko <- sum(ga[,1] > teta)/nrga * 100  # teta>GA.üst oraný
  og <- mean(ga[,2]-ga[,1])  # Ortalama GA geniþliði
 # Nominal GD sýnýrlarýna göre baþarý durumu
 # (1=baþarýlý, 0=baþarýsýz)
  bsg <- 0 
  nomgd <- c(conf.level-0.2*(1-conf.level),
    conf.level+0.2*(1-conf.level))*100
  if(ko >= conf.level*100) bsg <- 1
  else if(ko >= nomgd[1] & ko <= nomgd[2]) bsg <- 1
  c(round(c(ko, ako, uko),2), og, bsg)
 }

# Fonksiyon 4.21: Konum istatistikleri hesaplama fonksiyonu 
 calc.location <- function(x, indices, gama=c(0.1,0.2), q=0.5){
   if(is.vector(x)) x <- as.data.frame(x)
   if(missing(indices)) indices <- 1:nrow(x)
   x <- x[indices,]
   xmean <- mean(x)               # Aritmetik ortalama
   xmed <- median(x)              # Ortanca
   xtrmean1 <- mean(x, trim=gama[1])  # Bud. ort. (gama=0.1)
   xtrmean2 <- mean(x, trim=gama[2])  # Bud. ort. (gama=0.2)
   sonuc <- c(xmean, xmed, xtrmean1, xtrmean2)
   names(sonuc) <- c("Ortalama", "Ortanca", 
     "BdOrt-0.1", "BdOrt-0.2")
   return(sonuc)
 }

# Fonksiyon 4.22a: Kapsama oraný grafiði 
 ko.plot <- function(perfarr, index=NULL){
   colpal <- c("red", "blue", "violet", "green", "orange",
     "magenta", "tomato")
   if(is.null(index)) index <- 1: (dim(perfarr)[3])
   plot(perfarr[,1,index[1]], perfarr[,2,index[1]], 
     type="l", lwd=2, col=colpal[index[1]],
     ylab="KO(%)",ylim=c(50,100), xlab=dimnames(perfarr)[[2]][1],
     main="Kapsama Oraný")
   if(length(index) > 1){
      for(i in index[-1]){
        lines(perfarr[,1,index[1]], perfarr[,2,i], 
          lwd=2, col=colpal[i])
      }
   }
   legendnames <- dimnames(perfarr)[[3]]
   legend("bottomright", lty=1, col=colpal[index], 
      legend=legendnames[index], horiz=FALSE)
   abline(h=94, col="gray", lty=2)
   abline(h=96, col="gray", lty=2)
 }

# Fonksiyon 4.22b: Altta kalma oraný grafiði 
 ako.plot <- function(perfarr, index=NULL){
   colpal <- c("red", "blue", "violet", "green", "orange",
     "magenta", "tomato")
   if(is.null(index)) index <- 1: (dim(perfarr)[3])
   plot(perfarr[,1,index[1]], perfarr[,3,index[1]], 
     type="l", lwd=2, col=colpal[index[1]],
     ylab="%AKO", ylim=c(0,10), xlab=dimnames(perfarr)[[2]][1],
     main="Altta Kalma Oraný")
   if(length(index) > 1){
      for(i in index[-1]){
        lines(perfarr[,1,index[1]], perfarr[,3,i], 
          lwd=2, col=colpal[i])
      }
   }
   legendnames <- dimnames(perfarr)[[3]]
   legend("topright", lty=1, col=colpal[index], 
      legend=legendnames[index], horiz=FALSE)
   abline(h=2,col="gray")
   abline(h=3,col="gray")
 }

# Fonksiyon 4.22c: Üstte kalma oraný grafiði 
 uko.plot <- function(perfarr, index=NULL){
   colpal <- c("red", "blue", "violet", "green", "orange",
     "magenta", "tomato")
   if(is.null(index)) index <- 1: (dim(perfarr)[3])
   plot(perfarr[,1,index[1]], perfarr[,4,index[1]],
     type="l", lwd=2, col=colpal[index[1]],
     ylab="%ÜKO", ylim=c(0,10), xlab=dimnames(perfarr)[[2]][1],
     main="Üstte Kalma Oraný")
   if(length(index) > 1){
      for(i in index[-1]){
        lines(perfarr[,1,index[1]], perfarr[,4,i], 
          lwd=2, col=colpal[i])
      }
   }
   legendnames <- dimnames(perfarr)[[3]]
   legend("topright", lty=1, col=colpal[index], 
      legend=legendnames[index], horiz=FALSE)
   abline(h=2,col="gray", lty=2)
   abline(h=3,col="gray", lty=2)
 }

# Fonksiyon 4.22d: Ortalama geniþlik grafiði 
 og.plot <- function(perfarr, index=NULL){
   colpal <- c("red", "blue", "violet", "green", "orange",
     "magenta", "tomato")
   if(is.null(index)) index <- 1: (dim(perfarr)[3])
   ogmax <- max(perfarr[,5,index])
   plot(perfarr[,1,index[1]], perfarr[,5,index[1]],
     type="l", lwd=2, col=colpal[index[1]],
     ylab="OG", ylim=c(0,ogmax), xlab=dimnames(perfarr)[[2]][1], 
     main="Ortalama Geniþlik")
   if(length(index) > 1){
      for(i in index[-1]){
        lines(perfarr[,1,index[1]], perfarr[,5,i], 
         lwd=2, col=colpal[i])
      }
   }
   legendnames <- dimnames(perfarr)[[3]]
   legend("topright", lty=1, col=colpal[index], 
     legend=legendnames[index], horiz=FALSE)
 }

# Fonksiyon 4.22e: Baþarý durumu karþýlaþtýrma grafiði 
 bd.plot <- function(perfarr, index=NULL){
   colpal <- c("red", "blue", "violet", "green", "orange",
     "magenta", "tomato")
   if(is.null(index)) index <- 1: (dim(perfarr)[3])
   bd <- c(); j<-1
   for(i in index){
     bd[j] <- round(sum(perfarr[,6,i]) / dim(perfarr)[1] *100,2)
     j <- j+1
   }
   barnames <- dimnames(perfarr)[[3]]
   barplot(bd, col=colpal, ylab="%Baþarý", xlab="Kestirici", 
     names.arg=barnames[index], 
     main="Kestiricilerin Bootstrap Baþarýsý") 
 }

# Fonksiyon 4.23: Farklý bootstrap tekrarlarý için bootstrap 
# baþarým analizi
# Baðýmlýlýk – Fonksiyon: 4.3a, 4.20
#
 bs.perform2 <- function(x, teta, statistic,  
   gayontem=ga.norm, m=1000, n=30, R=c(2000), dbs=dbs,
     conf.level=0.95){
   if(is.vector(x)) x <- as.data.frame(x)
 # Güven aralýklarý 
   gvnara <- array(0, dim=c(m, 2, length(teta))) 
 # Kapsama oraný
   perfarr <- array(0, dim=c(length(R), 6, length(teta))) 
   dimnames(perfarr) <- list(R,c("R", "KO(%)",
     "AKO(%)","ÜKO(%)","OG","BD"),
   names(teta))
   for(i in 1:length(R)){
     for(j in 1:m){
       sidx <- sample(1:nrow(x), size=n, replace=TRUE)
       bsres <- dbs.nonpar1(x[sidx,], statistic=statistic,
         double=dbs, R=R[i])
       for(k in 1:length(teta)){
         gvnara[j,,k] <- gayontem(bsres, index=k,
           conf.level=conf.level)
         cat("Ýlerleme durumu: R=",R[i], "Örneklem=",j, "\n")
       }
     }
     for(k in 1:length(teta)){
       perfarr[i,,k] <- c(R[i], kapsama(ga=gvnara[,,k],
         teta=teta[k], conf.level))
     }
   }
   (perfarr)
 }

# Fonksiyon 4.24: Farklý güven aralýðý yöntemleri için bootstrap
# baþarým analizi
# Baðýmlýlýk – Fonksiyon: 4.3a, 4.20
#
 bs.perform3 <- function(x, teta, statistic,  
   gayontem=ga.yuzde, m=1000, n=30, dbs=dbs, 
   R=2000, conf.level=0.95){
   if(is.vector(x)) x <- as.data.frame(x)
 # Güven aralýklarý 
   gvnara <- array(0, dim=c(m, 2, length(teta))) 
   perfarr <- array(0, dim=c(length(gayontem), 6, length(teta))) + # Kapsama oraný
   dimnames(perfarr) <-list(names(gayontem),
     c("GA","KO(%)", "AKO(%)","ÜKO(%)","OG","BD"), names(teta))
    for(i in 1:length(gayontem)){
      gacalc <- gayontem[[i]]
      for(j in 1:m){
        sidx <- sample(1:nrow(x), size=n, replace=TRUE)
        bsres <- dbs.nonpar1(x[sidx,], statistic=statistic,
          double=dbs, R=R)
        for(k in 1:length(teta)){
          gvnara[j,,k] <- gacalc(bsres, index=k,
            conf.level=conf.level)
          cat("Ýlerleme durumu: GA=",names(gayontem)[i],
            "Örneklem=",j, "\n")
        }
      }
      for(k in 1:length(teta)){
        perfarr[i,,k] <- c(i, kapsama(ga=gvnara[,,k],
          teta=teta[k], conf.level))
      }
    }
  (perfarr)
 }

# Fonksiyon 4.25: g&h daðýlýþýndan rastlantýsal örneklem çekme
 rgh <- function(n=30, A=0, B=1, g=0, h=0){
   if(sign(B)!=1 | sign(h)==-1) 
     stop("B pozitif; h sýfýr veya pozitif olmalý!")
   x <- rnorm(n)
   if (g!=0){
     ghval <- (exp(g*x)-1)*exp(h*x^2/2)/g
   }else{
     ghval <- x*exp(h*x^2/2)
   }
   ghval <- A+B*ghval
   (ghval)
 }

# Fonksiyon 4.26a: Huber'in M-kestiricisi fonksiyonu 1
 calc.huber <- function(x, indices, c=1.28, 
   niter=25, conv=1e-04){
   if(is.vector(x)) x <- as.data.frame(x)
   if(missing(indices)) indices <- 1:nrow(x)
   h <- numeric(niter)
   h[1] <- median(x[indices,1])
   for(i in 1:niter){
     if(sign(mad(x[indices,1]))!= 0){
        A <- (x[indices,1]-h[i]) / mad(x[indices,1])
        nA <- sum(sapply(A,function(x){max(-c, min(c,x))}))
        B <- (A >= c | A <= -c)
        nB <- length(B[!B])
        h[i+1]<- h[i] + ((mad(x[indices,1])*nA) / nB)
     }
     if((abs(h[i+1])-abs(h[i]))<=conv) break
   }
   return(h)
  }

# Fonksiyon 4.26b: Huber'in Tek Adým M-kestirim fonksiyonu 2
# Baðýmlýlýk – Paket: asbio
#
 calc.hosme <- function(x, indices, c=1.28){
   if(!require(asbio, quietly = TRUE)) {
     install.packages("asbio",
      repo="https://cloud.r-project.org/")
     require(asbio, quietly = TRUE)
   }
   if(is.vector(x)) x <- as.data.frame(x)
   if(missing(indices)) indices <- 1:nrow(x)
   (huber.one.step(x[indices,1], c))
 }

# Fonksiyon 4.27: Harrell-Davis kestirimi (HD) fonksiyonu 1  
 calc.hde <- function(x, indices, q=0.5){
   if(is.vector(x)) x <- as.data.frame(x)
   n <- nrow(x)
   if(missing(indices)) indices <- 1:n
    b1 <- q*(n+1)
    b2 <- (1-q)*(n+1)
    w <- pbeta((1:n)/n, b1, b2)- pbeta(0: (n-1)/n, b1, b2)
    hde <- sum(w*sort(x[indices,1]))
    (hde)
 }

# Fonksiyon 4.28: Harrell-Davis kestirimi (HD) fonksiyonu 2  
# Baðýmlýlýk – Paket: Hmisc
#
 calc.hde <- function(x, indices, q=0.5){
   if(!require(Hmisc, quietly = TRUE)) {
     install.packages("Hmisc",
      repo="https://cloud.r-project.org/")
     require(Hmisc, quietly = TRUE)
   }
   if(is.vector(x)) x <- as.data.frame(x)
   if(missing(indices)) indices <- 1:nrow(x)
   (unname(hdquantile(x[indices,1], probs=q)))
 }

# Fonksiyon 4.29: MLE hesaplama fonksiyonu 
# Baðýmlýlýk – Paket: fitdistrplus 
#
 calc.mle <- function(x, xdist, indices){
   if(!require(fitdistrplus, quietly=TRUE)) {
     install.packages("fitdistrplus",
       repo="https://cloud.r-project.org/")
     require(fitdistrplus, quietly=TRUE)
   }
   if(is.vector(x)) x <- as.data.frame(x)
   if(missing(indices)) indices <- 1:nrow(x)
   mle <- fitdistr(x[indices,], xdist)$estimate
   (mle)
 }

# Fonksiyon 4.30a: Normal daðýlýþ benzetim fonksiyonu
 gen.norm <- function(x, mle){
  (rnorm(length(x), mean=mle[1], sd=mle[2]))
 }

# Fonksiyon 4.30b: Log-Normal daðýlýþ benzetim fonksiyonu
 gen.lognorm <- function(x, mle){
  (rlnorm(length(x), meanlog=mle[1], sdlog=mle[2]))
 }

# Fonksiyon 4.30c: Weibull daðýlýþý benzetim fonksiyonu
 gen.weibull <- function(x, mle){
   (rweibull(n=length(x), shape=mle[1], scale=mle[2]))
 }

# Fonksiyon 4.30d: Gama daðýlýþý benzetim fonksiyonu
 gen.gamma <- function(x, mle){
   (rweibull(n=length(x), shape=mle[1], rate=mle[2]))
 }

# Fonksiyon 4.30e:  Üssel daðýlýþ simülasyon fonksiyonu
 gen.expo <- function(x, mle){
   (rexp(n=length(x), rate=mle[1]))

# Fonksiyon 4.31: Parametrik bootstrap fonksiyonu 
 bs.par <- function(dset, statistic, ran.gen, mle, R=2000){
   if(!is.data.frame(dset) | !is.matrix(dset)) 
     dset <- as.data.frame(dset)
   n <- nrow(dset)  # Örneklem büyüklüðü
 # Orijinal örneklem kestirimleri
   tetahat <- statistic(dset[1:n,]) 
   nstats <- length(tetahat) # Ýstatistik sayýsý
 # Bootstrap kestirimleri
   tetahatstar <- matrix(0, nrow=R, ncol=nstats)  
   for (i in 1:R){
     tetahatstar[i,] <- statistic(ran.gen(dset[1:n,], mle=mle))
   }
   return(list(t0=tetahat, t=tetahatstar, statistic=statistic,
     ran.gen=ran.gen, mle=mle, R=R, call=match.call()))
 }

# Fonksiyon 4.32 : Yarý parametrik bootstrap fonksiyonu 1
 dbs.sempar <- function(dset, statistic, double=FALSE, 
   R=2000, R2=NULL){
   n <- length(dset)  # Örneklem büyüklüðü
 # Orijinal örneklem kestirimleri
   tetahat <- statistic(dset, indices=1:n) 
   nstats <- length(tetahat) # Ýstatistik sayýsý
 # Bootstrap kestirimleri matrisi
   tetastar <- matrix(NA, nrow=R, ncol=nstats)  
 # Bootstrap t deðerleri matrisi
   tvals <- matrix(NA, nrow=R, ncol=nstats)  
   for (i in 1:R){
     drset <- dset+rnorm(n)  # Veriye gürültü ekle
     indices <- sample(1:n, replace=TRUE)
     tetastar[i,] <- statistic(drset,indices)
     if(double){
     ## Ýç bootstrap baþý
        if(is.null(R2) | !is.numeric(R2)) R2 <- floor(R*0.1)
        dset2 <- drset[indices]
     #Ýç bootstrap kestirimleri 
        tetastar2 <- matrix(NA, nrow=R2, ncol=nstats)  
        for(j in 1:R2){
          indices <- sample(1:n, replace=TRUE)
          tetastar2[j,] <- statistic(dset2,indices)
        }
        tvals[i,] <- (tetastar[i,]-tetahat)/
           apply(tetastar2, 2, sd)
     } else{tvals <- NULL}
     ## Ýç bootstrap sonu
   }
   return(list(t0=tetahat, t=tetastar, tvals=tvals, data=dset,
     statistic=statistic, double=double, R=R, R2=R2,
     call=match.call()))
 }


#### BÖLÜM 5 ##########################################################

# Fonksiyon 5.1: Jackknife fonksiyonu 
> jackknife <- function(dset, indices, statistic){
+   if(!is.data.frame(dset) | !is.matrix(dset)) 
+     dset <- as.data.frame(dset)
+   n <- nrow(dset)  # Örneklem büyüklüðü
+   if(missing(indices)) indices <- 1:n
+   teta <- statistic(dset[indices, ])
+   nstats <- length(teta)
+   tetajack <- matrix(0, nrow=n, ncol=nstats)
+   colnames(tetajack) <- names(teta)
+   for(i in 1:n){
+     looidx <- indices[-i]
+     tetajack[i,] <- statistic(dset, indices=looidx)
+   }
+   tetabar <- apply(tetajack, 2, mean)
+   jack.bias <- (n-1) * (tetabar-teta)
+   jack.se <- sqrt((n-1) * colMeans(t(t(tetajack)-tetabar)^2))
+   return(list(teta=teta, tetabar=tetabar, 
+     jack.bias=jack.bias, jack.se=jack.se,  
+     jack.values=tetajack, data=dset, call=match.call()))
+ }

# Fonksiyon 5.2: Jackknife kestirimleri daðýlýþ grafiði
# Baðýmlýlýk – Fonksiyon: 5.3
#
jack.dist <- function(jackres, index=1, conf.level=0.95){
    jv <- (as.data.frame(jackres$jack.values))[,index]
    ga <- jack.ga(jackres, index, conf.level)
    hist(jv, prob=TRUE, col="gray90", 
      w=2*(IQR(jv)/length(jv)^(1/3)), 
      ylab="Yoðunluk",  ylim=c(0, 1.1*max(density(jv)$y)),
      xlab="Jackknife kestirimi", xlim=c(ga[1]-0.01*ga[1],
        ga[2]+0.01*ga[2]),
      main= paste("Jackknife örnekleme daðýlýþý – ",index))
    lines(density(jv), col=4)
    segments(mean(jv),0, mean(jv), max(density(jv)$y)+0.40,
      col=2, lwd=2)
    text(mean(jv), max(density(jv)$y)+0.50, srt=0, "t*")
    segments(ga[1],0, ga[1], max(density(jv)$y)+0.40, 
      col=3, lwd=2)
    segments(ga[2],0, ga[2], max(density(jv)$y)+0.40, 
      col=3, lwd=2)
    text(ga[1], max(density(jv)$y)+0.50, srt=0, "GS.alt")
    text(ga[2], max(density(jv)$y)+0.55, srt=0, "GS.üst")
}

# Fonksiyon 5.3: Jackknife kestirimleri için güven aralýðý
jack.ga <- function(jackres, index=1, conf.level=0.95){
    alpha <- 1-conf.level
    jv <- (as.data.frame(jackres$jack.values))[,index]
    jse <- jackres$jack.se[index]
    df <- length(jv)-1
    tetabar <- mean(jv) 
    jack.ga <- tetabar-qt(c(1-alpha/2, alpha/2), df)*jse
    (jack.ga)
}

# Fonksiyon 5.4: Açýklayýcý istatistikler kestirim fonksiyonu 2
calc.stats2 <- function(x, indices){
   if(!is.data.frame(x)) x <- as.data.frame(x)
   if(missing(indices)) indices <- 1:nrow(x)
   x <- x[indices,]
   xmean <- mean(x)  #Aritmetik ortalama
   xmed <- median(x) #Ortanca
   xQ1 <- unname(quantile(x, prob=c(0.25))) #Q1
   xQ3 <- unname(quantile(x, prob=c(0.75))) #Q3
   xIQR <- IQR(x)    #Çeyreklikler arasý geniþlik (Q3-Q1)
   xvar <- var(x)    #Varyans
   sonuc <- c(xmean, xmed, xQ1, xQ3, xIQR, xvar)
   return(sonuc)
}

#### BÖLÜM 6 #############################################################

# Fonksiyon 6.1a: Geniþ-Uzun veri dönüþümü
wide2long <- function(x, y=NULL){
   if(is.null(y)){
     y <- x[,1]
     x <- x[,2]
   }
   nx <- length(x)
   ny <- length(y)
   grup <- rep(c("G1", "G2"), c(nx, ny))
   uvset <- data.frame(gozlem=c(x,y), grup)
   (uvset)
}

# Fonksiyon 6.1b: Uzun-Geniþ veri dönüþümü
long2wide <- function(x){
   grup <- unique(x[,2])
   xtemp <- x[which(x[,2]==grup[1]), 1]
   ytemp <- x[which(x[,2]==grup[2]), 1]
   (data.frame(V1=xtemp, V2=ytemp))
}

# Fonksiyon 6.1c: Uzun-Vektörel veri dönüþümü fonksiyonu
long2vect <- function(x){
   grup <- unique(x[,2])
   xtemp <- x[which(x[,2]==grup[1]), 1]
   ytemp <- x[which(x[,2]==grup[2]), 1]
   list(V1=xtemp, V2=ytemp)
}

# Fonksiyon 6.2: Ýki deðiþken için karþýlaþtýrma grafikleri
# Baðýmlýlýk – Paket: vioplot
#
bivar.plot <- function(x, y=NULL){
   opar <- par(mfrow=c(3,2))
   if(is.null(y)) {
     if(is.factor(x[,2])){
       grup <- unique(x[,2])
       y <- x[which(x[,2]==grup[2]), 1]
       x <- x[which(x[,2]==grup[1]), 1]
     }else if(is.numeric(x[,2])){
       y <- x[,2]
       x <- x[,1]
     }else{
       stop("Ýki deðiþken gerekli")
     }
   }
   minx <- min(density(x)$x, density(y)$x)
   maxx <- max(density(x)$x, density(y)$x)
   miny <- min(density(x)$y, density(y)$y)
   maxy <- max(density(x)$y, density(y)$y)
   minecdfx <- min(min(x), min(y))
   maxecdfx <- max(max(x), max(y))
   plot(density(x), lwd=1, col=2, 
     xlim=c(minx,maxx), ylim=c(miny,maxy),
     xlab="", ylab="yoðunluk", main="Daðýlýþ grafiði")
   lines(density(y), lwd=1, col=4)
   legend("topright", col=c(2,4), legend=c("x","y"), lty=c(1,1))
   plot(ecdf(x), lwd=1, col=2, xlab="",
     xlim=c(minecdfx, maxecdfx), main="ECDF Grafiði")
   lines(ecdf(y), lwd=1, col=4)
   legend("topleft", col=c(2,4), legend=c("x","y"), 
     lty=c(1,1), bty="n")
   boxplot(x,y, cex=1, col=c(2,4), names=c("x","y"), pch=21,
     medcol="white", medlwd=1, main="Kutu-býyýk Grafiði")
   if(!require(vioplot, quietly=TRUE)){
      install.packages("vioplot")
       require(vioplot, quietly=TRUE)
   }
   vioplot(x, y, names=c("x","y"), col="gray90", border=4, 
   rectCol=c(2,4), colMed=7, pchMed=19, horizontal=FALSE)
   title("Keman grafiði")
   qqnorm(x, pch=21, bg="gray90", xlab="teorik kantiller", 
     ylab="örneklem kantilleri", main="QQ Norm(x)")
   qqline(x, col=2)
   qqnorm(y, pch=21,bg="gray90", xlab="teorik kantiller", 
     ylab="örneklem kantilleri", main="QQ Norm(y)")
   qqline(y, col=4)
   par(opar)
}

# Fonksiyon 6.3: Ýki örneklemin karþýlaþtýrýlmasý için permütasyon # testi 
# Baðýmlýlýk – Fonksiyon 6.1a
#
permute.dif <- function(x, y=NULL, statistic, 
   alternative="two.sided", R=10000){
 # Veri uzun formatta deðilse uzun biçime dönüþtür
   if(!is.null(y)){
     x <- wide2long(x,y)
   }
   n <- tapply(x[,1], x[,2], length)
   grplev <- unique(x[,2])
 # Orijinal örneklemde farklarý hesapla
   tetahat <- statistic(x)
   nstats <- length(tetahat)
   tetahatstar <- matrix(NA, nrow=R, ncol=nstats)
   pv <- numeric(nstats)
   for(i in 1:R) {
     indices <- sample(1: (n[1]+n[2]), replace=FALSE)
     resampx <- x[indices,]
     resampx[1:n[1],2] <- grplev[1]
     resampx[-c(1:n[1]),2] <- grplev[2]
     tetahatstar[i,] <- statistic(resampx)
  }
  for(i in 1:nstats){  
    if(alternative == "less"){
      pv[i] <- sum(tetahatstar[,i] < tetahat[i])/R
    }else if(alternative == "greater"){
      pv[i] <- sum(tetahatstar[,i] > tetahat[i])/R
    }else{
      pv[i] <- sum(abs(tetahatstar[,i]) >= abs(tetahat[i]))/R
    } 
  }
  return(list(t0=tetahat, t=tetahatstar, p.value=pv,
    alternative=alternative, R=R, call=match.call()))
}

# Fonksiyon 6.4: Uzun veride ortalamalarýn farkýný hesaplama 
calc.meandif <- function(x, indices){
   if(missing(indices)) indices <- 1:nrow(x)
   delta <- diff(rev(tapply(x[indices,1], x[indices,2], mean)))
   #names(delta) <- "d.Ort"
   return(delta)
}

# Fonksiyon 6.5: Uzun veride ortancalarýn farkýný hesaplama 
calc.meddif <- function(x, indices){
   if(missing(indices)) indices <- 1:nrow(x)
   delta <- diff(rev(tapply(x[indices,1], x[indices,2], median)))
   #names(delta) <- "d.Ortc"
   return(delta)
}

# Fonksiyon 6.6. HD kantil kestirimi farklarýný hesaplama 
# Baðýmlýlýk – Fonksiyon: 3.28b
#
calc.hddif <- function(x, indices){
   if(missing(indices)) 
     indices <- 1:nrow(x)
   delta <- diff(rev(tapply(x[indices,1], x[indices,2], calc.hde)))
   #names(delta) <- "d.HD"
   return(delta)
}

# Fonksiyon 6.7: Uzun veride çoklu istatistiklerin farkýný 
# hesaplama fonksiyonu
calc.multdif <- function(x, indices){
   if(missing(indices)) indices <- 1:nrow(x)
   dMean <- diff(rev(tapply(x[indices,1], x[indices,2], mean)))
   dMed <- diff(rev(tapply(x[indices,1], x[indices,2], median)))
   dIQR <- diff(rev(tapply(x[indices,1], x[indices,2], IQR)))
   dVar <- diff(rev(tapply(x[indices,1], x[indices,2], var)))
   delta <- c(dMean, dMed, dIQR, dVar)
   #names(delta) <- c("d.Ort", "d.Ortc", "d.IQR", "d.Var")
   return(delta)
}

# Fonksiyon 6.8: Eþleþtirilmiþ örneklemlerde ortalamalar farkýný
# hesaplama fonksiyonu
calc.meandif2 <- function(x, indices){
   n <- nrow(x)
   if(missing(indices)) indices <- 1:n
   x <- x[indices,]
   delta <- mean(x[c(1: (n/2)),1]-x[-c(1: (n/2)),1])
   #names(delta) <- "d.Ort"
   return(delta)
}

# Fonksiyon 6.9: Eþleþtirilmiþ örneklemlerde ortancalar farkýnýn
# hesaplanmasý
calc.meddif2 <- function(x, indices){
   n <- nrow(x)
   if(missing(indices)) indices <- 1:n
   x <- x[indices,]
   delta <- median(x[c(1: (n/2)),1]-x[-c(1: (n/2)),1])
   #names(delta) <- "d.Ortc"
   return(delta)
}

# Fonksiyon 6.10: Eþleþtirilmiþ örneklemlerde istatistikler arasý
# farklarý hesaplama
calc.multdif2 <- function(x, indices){
   n <- nrow(x)
   if(missing(indices)) indices <- 1:n
   x <- x[indices,]
   difx <- x[c(1: (n/2)),1]-x[-c(1: (n/2)),1]
   dMean <- mean(difx)
   dMed <- median(difx)
   dIQR <- IQR(difx)
   dVar <- var(difx)
   dCv <- sd(difx)/mean(difx)
   dCv2 <- IQR(difx)/median(difx)
   delta <- c(dMean, dMed, dIQR, dVar, dCv, dCv2)
   return(delta)
}

#### BÖLÜM 7 #############################################################

# Fonksiyon 7.1a: Korelasyon ve regresyon analizi için veri
# benzetimi
# Baðýmlýlýk – Fonksiyon: 7.1b
#
gen.funcdata <- function(x, fx, a=0, b=1, nmu=0, nsigma=0){
   noise <- rnorm(length(x), nmu, nsigma)
   if(missing(fx)) fx <- fLinear
   y <- fx(x, b, a, noise)
   (data.frame(y=y, x=x))
}

# Fonksiyon 7.1b: Çeþitli iliþkisel fonksiyonlar
fLinear <- function(x,b,a, noise) (b*x + a + noise)
fQuad <- function(x,b,a, noise) (b*x^2 + a + noise)
fCubic <- function(x,b,a, noise) (b*x^3 + a + noise)
fPoly <- function(x,b,a, noise) (b*x+b*x^2+b*x^3 + a + noise)
fCos <- function(x,b,a, noise) (b*cos(x*pi) + a + noise)
fSin <- function(x,b,a, noise) (b*sin(x*pi) + a + noise)
fExp <- function(x,b,a, noise) (exp(-x) + a + noise)
fVolcano <- function(x,b,a, noise) (log(abs(x)) + a + noise)
fSemCirc <- function(x,b,a, noise) (b*sqrt(max(x^2)-(x^2)) + a + noise)
fNcVar <- function(x,b,a, noise) (b*x + a + noise) # x ile deðiþen varyanslar

# Fonksiyon 7.2: Ýki deðiþken arasýndaki iliþkiler grafikler
cor.plot <- function(x, y=NULL){
   opar <- par(mfrow=c(3,2))
   if(is.null(y)) {
     if(is.factor(x[,2])){
       grup <- unique(x[,2])
       y <- x[which(x[,2]==grup[2]), 1]
       x <- x[which(x[,2]==grup[1]), 1]
     }else if(is.numeric(x[,2])){
       y <- x[,1]
       x <- x[,2]
     }else{
       stop("Ýki deðiþken gerekli")
     }
   }
   minx <- min(density(x)$x, density(y)$x)
   maxx <- max(density(x)$x, density(y)$x)
   miny <- min(density(x)$y, density(y)$y)
   maxy <- max(density(x)$y, density(y)$y)
   minecdfx <- min(min(x), min(y))
   maxecdfx <- max(max(x), max(y))
   plot(density(x), lwd=1, col=2, xlim=c(minx,maxx),
     ylim=c(miny,maxy),xlab="", ylab="yoðunluk", 
     main="Daðýlýþ grafiði")
   lines(density(y), lwd=1, lty=2, col=4)
   legend("topright", col=c(2,4), legend=c("x","y"), lty=c(1,2),
     bty="n", horiz=TRUE)
   if(!require(vioplot)) 
     {install.packages("vioplot"); require(vioplot)}
   vioplot(x, y, names=c("x","y"), col="gray90", border=1, 
     rectCol=c(2,4), colMed=7, pchMed=19, horizontal=FALSE)
   title("Keman grafiði")
   qqnorm(y, pch=21,bg="gray90", 
     xlab="Teorik kantiller", ylab="Örneklem kantilleri",
     main="Q-Q Norm(y)")
   qqline(y, col=4)
   qqnorm(x, pch=21, bg="gray90", xlab="Teorik kantiller",
     ylab="Örneklem kantilleri", main="Q-Q Norm(x)")
   qqline(x, col=2)
   plot(x,y, cex=1, col=1, bg="gray90", xlab="x",ylab="y",
     pch=21, main="Serpilme Grafiði")
   abline(lm(y~x), col=3)
   lines(lowess(x, y), col=2, lwd=1)
   plot(x, type="l", lwd=1, col=2, yaxt='n', xlab="Gözlem",
     ylab="x", main="Birlikte Deðiþim Grafiði")
   axis(2, pretty(c(min(x), 1.1*max(x))), col=1)
   par(new = TRUE)
   plot(y, type="l", lwd=1, lty=2, col=4, xaxt='n', yaxt='n',
     xlab="", ylab="") 
   axis(4, pretty(c(min(y), 1.1*max(y))), col=1)
   mtext("y", side=4, cex=0.8, line=0)
   legend("bottomright", col=c(2,4), legend=c("x","y"),
     lty=c(1,2), bty="n", horiz=TRUE)
   par(opar)
}

# Fonksiyon 7.3: Çembersel veri benzetimi
gen.circle <- function(n, noise=0){
    t <- runif (n, 0, 2*pi)
    x <- append(cos(t) , cos(t))
    t <- t+runif(n, 0, noise*pi) 
    y <- append(sin(t), -sin(t))
   (data.frame(y=y, x=x))
}

# Fonksiyon 7.4: Sýnýf etiketli spiral veri üretimi (Ma, 2016)
gen.cspiral <- function(n, k, seed=NULL){
   if(!is.null(seed)) set.seed(seed)
   X <- data.frame() # Veri matrisi 
   y <- data.frame() # sýnýf etiketleri
   for (j in (1:k)){
     r <- seq(0.05,1,length.out = n) # yarýçap
     t <- seq((j-1)*4.7,j*4.7, length.out = n) + rnorm(n, 
       sd = 0.3) # teta
     xtemp <- data.frame(x =r*sin(t) , y = r*cos(t)) 
     ytemp <- data.frame(matrix(j, n, 1))
     X <- rbind(X, xtemp)
     y <- rbind(y, ytemp)
   }
   dset <- cbind(X,y)
   colnames(dset) <- c(colnames(X), 'sinif')
   (dset)
}

# Fonksiyon 7.5: Spiral veri oluþturma fonksiyonu
gen.spiral <- function(n, aci, ad){
   t <- (1:n)*aci
   x <- sin(t)
   y <- cos(t)
   df <- data.frame(x=t*x, y=t*y)  # Veri kümesi
   return(list(df=df, ad=ad))
}

# Fonksiyon 7.6: Ýki deðiþkenli korelasyonlu veri benzetimi
gen.cordata <- function(n, r=0, mu=0, sd=1, a=0, seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    xt <- rnorm(n, mu, sd)
    x <- rnorm(n, mu, sd)
    y <- scale(x) * r + scale(residuals(lm(xt~x))) * 
      sqrt(1-r*r)+a
    (data.frame(y=y, x=x))
}

# Fonksiyon 7.7: Çok deðiþkenli korelasyonlu veri benzetimi
# Baðýmlýlýk – Paket: MASS
#
gen.mvcordata <- function(n, mu, covmat, seed=NULL){
    if(!require(MASS, quietly=TRUE)){
      install.packages("MASS",
        repo="https://cloud.r-project.org/")
      require(MASS, quietly=TRUE)
    }
    if(!is.null(seed)) set.seed(seed)
    df <- mvrnorm(n, mu=mu, Sigma=covmat, empirical=TRUE)
    # Korelasyonu kesin olarak üretmede empirical TRUE olmalýdýr
    (as.data.frame(df))
}

# Fonksiyon 7.8a: HHG paketi ile veri benzetimi
# Baðýmlýlýk – Paket: HHG
#
gen.hhgdata <- function(n, type="W"){
 # HHG paketini çalýþma alanýna yükle, kurulu deðilse kur
   if(!require(HHG, quietly=TRUE)){
     install.packages("HHG", repo="https://cloud.r-project.org/")
     require(HHG, quietly=TRUE)
   }
   X <- hhg.example.datagen(n, type)
   (data.frame(y=X[2,], x=X[1,]))
}

# Fonksiyon 7.8b: mlbench paketi ile veri benzetimi fonksiyonu
# Baðýmlýlýk – Paket: mlbench
#
gen.benchdata <- function(n, rs1=1, rs2=1, rs3=1){
# mlbench paketini çalýþma alanýna yükle, kurulu deðilse kur
   if(!require(mlbench, quietly=TRUE)){
      install.packages("mlbench",
        repo="https://cloud.r-project.org/")
     require(mlbench, quietly=TRUE)
   }
   df <- mlbench.cassini(n, relsize=c(rs1,rs2,rs3))
   return(list(m=df$x,k=df$classes))
}

# Fonksiyon 7.9a: Pearson korelasyonu hesaplama
cor.pearson <- function(x, y=NULL, indices){
   if(is.null(y)) {
     if(is.factor(x[,2])){
       grup <- unique(x[,2])
       y <- x[which(x[,2]==grup[2]), 1]
       x <- x[which(x[,2]==grup[1]), 1]
     }else if(is.numeric(x[,2])){
       y <- x[,2]
       x <- x[,1]
     }else{
       stop("Ýki deðiþken gerekli")
     }
   }
   n <- length(x)
   if(missing(indices)) indices <- 1:n
   y <- y[indices]
   x <- x[indices]
   sx <- sum(x)
   sy <- sum(y)
   sxy <- sum(x*y)
   sxx <- sum(x^2)
   syy <- sum(y^2)
   corxy <- (n*sxy-sx*sy)/(sqrt(n*sxx-sx^2)*sqrt(n*syy-sy^2))
   return(round(corxy,4))
}

# Fonksiyon 7.9b: Pearson korelasyonu hesaplama fonksiyonu 2
calc.pearson <- function(dset, indices){
   if(missing(indices)) indices <- 1:nrow(dset)
    dset <- dset[indices,]
    y <- dset[,1];  x <- dset[,2]
   (cor(x,y, method="pearson"))
}

# Fonksiyon 7.10a: Spearman korelasyonu hesaplama fonksiyonu
cor.spearman <- function(x, y=NULL, indices){
   if(is.null(y)) {
     if(is.factor(x[,2])){
       grup <- unique(x[,2])
       y <- x[which(x[,2]==grup[2]), 1]
       x <- x[which(x[,2]==grup[1]), 1]
     }else if(is.numeric(x[,2])){
       y <- x[,2]
       x <- x[,1]
     }else{
       stop("Ýki deðiþken gerekli")
     }
   }
   n <- length(x)
   if(missing(indices)) indices <- 1:n
   y <- y[indices]
   x <- x[indices]
   d <- rank(x)-rank(y)
   corxy <- 1-(6*sum(d^2)/(n^3-n))
   return(round(corxy, 4))
}

# Fonksiyon 7.10b: Spearman korelasyonu hesaplama fonksiyonu 2
calc.spearman <- function(dset, indices){
   if(missing(indices)) indices <- 1:nrow(dset)
    dset <- dset[indices,]
    y <- dset[,1];  x <- dset[,2]
   (cor(x,y, method="spearman"))
}

# Fonksiyon 7.11a: Kendall korelasyonu hesaplama fonksiyonu
cor.kendall <- function(x, y=NULL, indices){
   if(is.null(y)) {
     if(is.factor(x[,2])){
       grup <- unique(x[,2])
       y <- x[which(x[,2]==grup[2]), 1]
       x <- x[which(x[,2]==grup[1]), 1]
     }else if(is.numeric(x[,2])){
       y <- x[,2]
       x <- x[,1]
     }else{
       stop("Ýki deðiþken gerekli")
     }
   }
   n <- length(x)
   if(missing(indices)) indices <- 1:n
   y <- y[indices]
   x <- x[indices]
   nc <- nd <- 0
   for(j in 1:length(indices)){
     nc <- nc + sum(((x[j] < x) & (y[j] < y)) | ((x[j] > x)
        & (y[j] > y)))
     nd <- nd + sum(((x[j] < x) & (y[j] > y)) | ((x[j] > x) 
       & (y[j] < y)))
   }
   nc <- nc/2; nd <- nd/2
   corxy <- (nc-nd)/(nc+nd)
   # corxy <- (nc-nd)/(n*(n-1)/2) # Alternatif formül
   (corxy)
}

# Fonksiyon 7.11b: Kendall-tau korelasyonu hesaplama fonksiyonu 2
calc.kendall <- function(dset, indices){
   if(missing(indices)) indices <- 1:nrow(dset)
    dset <- dset[indices,]
    y <- dset[,1];  x <- dset[,2]
   (cor(x,y, method="kendall"))
}

# Fonksiyon 7.12: Uzaklýklar matrisi oluþturma fonksiyonu
distMat <- function(x, dm="euclidean", pwr=2){
   x <- as.matrix(dist(x, method=dm, p=pwr))
   xsutort <- colMeans(x)
   xsatort <- rowMeans(x)
   xgenort <- mean(x)
   xstar <- sweep (x, 2, xsutort, FUN="-")
   xstar <- sweep(xstar, 1, xsatort, FUN="-")
   xstar <- xstar + xgenort
   (xstar)
}

# Fonksiyon 7.13 : Kovaryans fonksiyonu
dCov <- function(x, y) (sqrt(mean(x*y)))

# Fonksiyon 7.14: Uzaklýk korelasyonu hesaplama fonksiyonu 1
# Baðýmlýlýk – Fonksiyon: 7.12, 7.13
#
cor.dCor <- function(dset, indices, dm="euclidean", pwr=2){
   if(missing(indices)) indices <- 1:nrow(dset)
   dset <- dset[indices,]
   y <- dset[,1];  x <- dset[,2]
   xstar <- distMat(x, dm, pwr)
   ystar <- distMat(y, dm, pwr)
   dCovxy <- dCov(xstar, ystar)
   dCovxx <- dCov(xstar, xstar)
   dCovyy <- dCov(ystar, ystar)
   if(dCovxx*dCovyy==0) 
     dCorxy <- 0
   else 
     dCorxy <- dCovxy / sqrt(dCovxx*dCovyy)
   (dCorxy)
} 

# Fonksiyon 7.15: Uzaklýk korelasyonu hesaplama fonksiyonu 2
# Baðýmlýlýk – Paket: energy
#
calc.dCor <- function(dset, indices){
   if(missing(indices)) indices <- 1:nrow(dset)
   dset <- dset[indices,]
   y <- dset[,1];  x <- dset[,2]
   if(!require(energy, quietly=TRUE)) {
     install.packages("energy",
       repo="https://cloud.r-project.org/")
     require(energy, quietly=TRUE)
   }
   (energy::dcor(x,y))
}

# Fonksiyon 7.16a: Hoeffding D hesaplama fonksiyonu
cor.hoefd <- function(dset, indices){
   if(missing(indices)) indices <- 1:nrow(dset)
   dset <- dset[indices,]
   n <- nrow(dset)
   y <- dset[,1];  x <- dset[,2]
   R <- rank(x) ;  S <- rank(y)
   Q <- matrix(0, nrow=n, ncol=1)
   for(i in 1:n){
  # Qi, hem x hem de y bakýmýndan j. gözlem çiftinden küçük
  # olanlarýn sayýsýna 1 eklenerek baþlatýlýr.
     Q[i] <- 1 + sum(R < R[i] & S < S[i])
  # Hem x hem de y de baðlý gözlem çiftleri Qi'ye 
  # 1/4 katký saðlar
     Q[i] <- Q[i] + 0.25 * (sum(R == R[i] & S == S[i])-1)
  # Sadece x'de baðlý veya y'de baðlý gözlem çiftleri Qi'ye
  # 1/2 katký saðlar. Ancak deðiþkenlerden baðlý olmayanýn
  # j. gözlem çiftinin i.'den küçük olmasý gerekir.
     Q[i] <- Q[i] + 0.5 * sum(R == R[i] & S < S[i])  
     Q[i] <- Q[i] + 0.5 * sum(R < R[i] & S == S[i])  
   }
   D1 <- sum((Q-1)*(Q-2))
   D2 <- sum((R-1) *(R-2) *(S-1) *(S-2))
   D3 <- sum((R-2) *(S-2) *(Q-1))
   D <- 30*((n-2)*(n-3)*D1+D2-2*(n-2)*D3) / 
     (n*(n-1)*(n-2)*(n-3)*(n-4))
   (D)
}

# Fonksiyon 7.16b: Hmisc ile Hoeffding D hesaplama fonksiyonu 
calc.hoefd <- function(dset, indices){
   if(missing(indices)) indices <- 1:nrow(dset)
   dset <- dset[indices,]
   y <- dset[,1];  x <- dset[,2]
   if(!require(Hmisc, quietly=TRUE)) {
     install.packages("Hmisc",
       repo="https://cloud.r-project.org/")
     require(Hmisc, quietly=TRUE)
   }
   (hoeffd(x, y)$D[2,1])
}

# Fonksiyon 7.17: Shannon entropisi hesaplama fonksiyonu
calc.entropy <- function(x, base=2){
   px <- prop.table(x)
   logpx <- log(px, base)
   logpx[which(logpx==-Inf)] <- 0
   Hx <- -1*sum(px*logpx) 
   (Hx)
}

# Fonksiyon 7.18: Koþullu entropi hesaplama fonksiyonu
calc.condentropy <- function(x, base=2){
   pxy <- prop.table(x) # Ortak olasýlýklar
   px <- rowSums(pxy) # x marjinal olasýlýklar
   py <- colSums(pxy) # y marjinal olasýlýklar
   Hxy <- calc.entropy(pxy, base)
   Hx <- calc.entropy(px, base)
   HCxy <- Hxy – Hx
   return(HCxy)
}

# Fonksiyon 7.19: Kullback-Leibler ýraksamasý hesaplama fonksiyonu
# Baðýmlýlýk – Fonksiyon: 7.17
#
calc.DKL <- function(x, base=2){
   pxy <- prop.table(x)
   px <- rowSums(x)
   py <- colSums(x)
   DKL <- calc.entropy(pxy, base)-calc.entropy(px, base)
   return(DKL)
}

# Fonksiyon 7.20: Çapraz entropi hesaplama fonksiyonu
calc.crossentropy <- function(x, base=2){
   pxy <- unname(prop.table(x))
   px <- rowSums(pxy)
   py <- colSums(pxy)
   logpy <- log(py, base)
   logpy[which(logpy==-Inf)] <- 0
   #Hxy <- 0
   #for(i in 1:length(x)){
   #  Hxy <- Hxy + px[i]*logpy[i]
   #}
   HCxy <- sum(px*logpy)
   return(-1*HCxy)
}

# Fonksiyon 7.21a: Karþýlýklý bilgi hesaplama fonksiyonu
# Baðýmlýlýk – Fonksiyon: 7.17
#
calc.MI <- function(x, base=2){
   pxy <- prop.table(x) # Pi frekanslarý
   px <- rowSums(pxy) #Marjinal x
   py <- colSums(pxy) #Marjinal y
   Hx <- calc.entropy(px, base)
   Hy <- calc.entropy(py, base)
   Hxy <- calc.entropy(pxy, base)
   I <- Hx+Hy-Hxy
   return(I)
}

# Fonksiyon 7.21b: Koþullu entropi hesaplama fonksiyonu 2
calc.MI2 <- function(x, base=2){
   pxy <- prop.table(x) # Pi frekanslarý
   px <- rowSums(pxy) #Marjinal x
   py <- colSums(pxy) #Marjinal y
   pxpy <- px %o% py # Qi frekanslarý 
   logpxpy <- log(pxy/pxpy, base)
   logpxpy[which(logpxpy==-Inf|logpxpy==Inf|logpxpy=="NaN")]<-0
   I <- sum(pxy*logpxpy)
   return(I)
}

# Fonksiyon 7.22a: MIC hesaplama fonksiyonu
# Baðýmlýlýk – Fonksiyon: 7.21a
cor.MIC <- function(x, y, alpha=0.6, base=2){
   n <- length(x)
   bmax <- max(ceiling(n^alpha),4)
   Istars <- c()
   for(bx in 2:bmax) {
     for(by in 2:bmax){
       if(bx*by > bmax){
         next
       }
       xdisc <- cut(x, breaks=bx, labels=1:bx)
       ydisc <- cut(y, breaks=by, labels=1:by)
       txy <- table(xdisc, ydisc)
       Istar <- calc.MI(txy, base)
       normIstar <- Istar / log(min(bx,by), base)
       Istars <- append(Istars, normIstar)
     }
   }
   MIC <- max(Istars)
   (MIC) 
}

# Fonksiyon 7.22b: MIC hesaplama fonksiyonu 1 
# Baðýmlýlýk – Fonksiyon: 7.21a
#
calc.mic <- function(dset, indices, base=2){
   if(missing(indices)) indices <- 1:nrow(dset)
   dset <- dset[indices,]
   y <- dset[,1];  x <- dset[,2]
   (cor.MIC(x, y, base))
}

# Fonksiyon 7.22c: MIC hesaplama fonksiyonu 2 
# Baðýmlýlýk – Paket: minerva
#
calc.MIC <- function(dset, indices){
   if(missing(indices)) indices <- 1:nrow(dset)
   dset <- dset[indices,]
   y <- dset[,1];  x <- dset[,2]
   if(!require(minerva, quietly=TRUE)) {
     install.packages("minerva", 
       repo="https://cloud.r-project.org/")
     require(minerva, quietly=TRUE)
   }
   (minerva::mine(x, y, n.cores=4)$MIC)
}

# Fonksiyon 7.23a: RDC hesaplama fonksiyonu 1
cor.rdc <- function(dset, indices, k=20, s=1/6, func=sin){
   if(missing(indices)) indices <- 1:nrow(dset)
   dset <- dset[indices,]
   y <- dset[,1];  x <- dset[,2]
   avgrvar <- function(var) (rank(var)/length(var))
   x <- cbind(apply(as.matrix(x), 2, avgrvar),1)
   y <- cbind(apply(as.matrix(y), 2, avgrvar),1)
   x <- s / ncol(x)* x %*% matrix(rnorm(ncol(x) * k), ncol(x))
   y <- s / ncol(y)* y %*% matrix(rnorm(ncol(y) * k), ncol(y))
   ccor <- cancor(cbind(func(x), 1), cbind(func(y), 1))$cor
   rdc <- ccor[1]
   nccor <- length(ccor)
   chisq <- (((2 * nccor + 3) / 2)-nrow(x)) * log(prod(1-ccor^2))
   p.value <- pchisq(chisq, nccor^2, lower.tail=FALSE)
   rdc <- list(rdc=rdc, p.value=p.value)
   names(rdc) <- c("rdc", "p.value")
   (rdc)
}

# Fonksiyon 7.23b: RDC hesaplama fonksiyonu 2
# Baðýmlýlýk – Paket: AlterCorr
#
calc.rdc <- function(dset, indices, k=20, s=1/6, f=sin){
   if(missing(indices)) indices <- 1:nrow(dset)
   dset <- dset[indices,]
   y <- dset[,1]; x <- dset[,2]
  # AlterCorr paketinin yüklenmesi, kurulu deðilse kurulmasý
   if(!require("AlterCorr", quietly=TRUE)){
     if(!require(remotes))
       install.packages("remotes")
     remotes::install_github("AnaBPazos/AlterCorr")
     require(AlterCorr, quietly=TRUE) 
   }
   rdc <- unname(AlterCorr(x, y, type="RDC"))[[1]]
   (rdc)
}

# Fonksiyon 7.24: nlcor ile doðrusal olmayan korelasyon hesaplama
# Baðýmlýlýk – Paket: nclor
#
calc.nlcor <- function(dset, indices){
   if(missing(indices)) indices <- 1:nrow(dset)
   dset <- dset[indices,]
   y <- dset[,1]; x <- dset[,2]
   if(!require(nlcor, quietly=TRUE)){
     if(!require(devtools)){
       install.packages("devtools")
       library(devtools)
     }
     devtools::install_github("ProcessMiner/nlcor")
     require(nlcor, quietly=TRUE) 
   }
   nlcor <- unname(nlcor::nlcor(x,y))[[1]]
   (nlcor)
}

# Fonksiyon 7.25:acepack ile doðrusal olmayan korelasyon hesaplama
# Baðýmlýlýk – Paket: acepack
#
calc.acecor <- function(dset, indices){
   if(missing(indices)) indices <- 1:nrow(dset)
   dset <- dset[indices,]
   y <- dset[,1]; x <- dset[,2]
   if(!require(acepack, quietly=TRUE)) {
     install.packages("acepack",
       repo="https://cloud.r-project.org/")
     require(acepack, quietly=TRUE)
   }
   argmax <- acepack::ace(x, y)
   acecor <- cor(argmax$tx, argmax$ty)[1]
   (acecor)
}

# Fonksiyon 7.26: NNS ile doðrusal olmayan korelasyon hesaplama
# Baðýmlýlýk – Paket: NNS
#
calc.nnscor <- function(dset, indices){
   if(missing(indices)) indices <- 1:nrow(dset)
   dset <- dset[indices,]
   y <- dset[,1]; x <- dset[,2]
   if(!require(NNS, quietly=TRUE)) {
      install.packages("NNS",
        repo="https://cloud.r-project.org/")
      require(NNS, quietly=TRUE)
   }
   nnscor <- NNS::NNS.cor(x, y)
   (nnscor)
}

# Fonksiyon 7.27: Farklý veri kümeleri için korelasyonlar tablosu 
# Baðýmlýlýk – Fonksiyon: 7.9b, 7.10b, 7.11b, 7.16b, 7.22b, 7.23b,
#   7.24, 7.25, 7.26; Paket: ggplot2
#
cor.table <- function(dsets, plot=FALSE){
   if(!require(ggplot2, quietly=TRUE)) {
     install.packages("ggplot2",
       repo="https://cloud.r-project.org/")
     require(ggplot2, quietly=TRUE)
   }
   ccor <- matrix(NA, nrow=length(dsets), ncol=10)
   colnames(ccor) <- c("pear", "sper", "kend", "dcor", "hoefd",
    "mic", "rdc", "ace", "nlcor", "nnscor")
   rownames(ccor) <- names(dsets)
   i=1
   for(dset in dsets){
      ccor[i,1] <- calc.pearson(dset)
      ccor[i,2] <- calc.spearman(dset)
      ccor[i,3] <- calc.kendall(dset)
      ccor[i,4] <- calc.dCor(dset)
      ccor[i,5] <- calc.hoefd(dset)
      ccor[i,6] <- calc.MIC(dset)
      ccor[i,7] <- calc.rdc(dset)
      ccor[i,8] <- calc.acecor(dset)
      ccor[i,9] <- calc.nlcor(dset)
      ccor[i,10] <- calc.nnscor(dset)
      if(plot){  
        # Çizilen grafikleri çalýþma klasörüne kaydedilir.
        graphics.off() # Açýk grafikleri kapat
        jpeg(paste("ds_",i,".jpeg"))
          print(ggplot2::ggplot(dset, aes(x, y)) + 
            geom_point(color = "orange") +
            geom_smooth(method = "loess", se=TRUE))
        dev.off()
      }
      i <- i+1
   }
   (ccor)
}

# Fonksiyon 7.28: Korelasyonlar arasý iliþkiler grafiði
# Baðýmlýlýk – Fonksiyon: 7.9b, 7.22b
#
cor.compare <- function(dset, ns=30, cor1=calc.pearson,
   cor2=calc.MIC, R=100){
   opar <- par(mfrow=c(2,1))
   cormat <- matrix(0, nrow=R,ncol=2)
   cormat <- data.frame(cor1=cormat[,1], cor2=cormat[,2])
   for(i in 1:R){
     indices <- sample(1:nrow(dset), size=n, replace=TRUE)
     cormat[i,] <- c(cor1(dset,indices), cor2(dset,indices))
   }
   plot(cormat, pch=21, col=1, bg="gray", cex=1, 
     xlab="Cor1", ylab="Cor2", main="Serpilme grafiði")
   abline(v=mean(cormat[,1]), col=2)
   abline(h=mean(cormat[,2]), col=4)
   abline(v=median(cormat[,1]), col=2, lty=2)
   abline(h=median(cormat[,2]), col=4, lty=2)
   abline(lm(cor2~cor1, data=cormat), col="orange")
   plot(cormat[,1], pch=21, col=1, bg="gray90", cex=0.2,
     ylim=c(min(cormat),max(cormat)), xlab="Örneklemler",
     ylab="Korelasyon", main="Örneklemlere göre deðiþim")
   points(cormat[,2], pch=21, col=1, bg="gray90", cex=0.2)
   lines(cormat[,1], pch=21, col=2)
   lines(cormat[,2], pch=21, col=4)
   par(opar)
   (cormat)
}

# Fonksiyon 7.29: Korelasyon matrisi fonksiyonu
# Baðýmlýlýk- Fonksiyon:7.9b, 7.10b, 7.11b, 7.15, 7.16b,
#   7.22a,7.32b, Paket: corrplot
#
cor.matrix <- function(x, method="pearson", plot=TRUE){
   if(!(is.matrix(x)|is.data.frame(x)))  
      stop("x sayýsal bir matris veya veri çerçevesi olmalý")
   if(!require(corrplot, quietly=TRUE)) {
     install.packages("corrplot", 
      repo="https://cloud.r-project.org/"); 
     require(corrplot, quietly=TRUE) }
   nv <- ncol(x)
   cormat <- diag(nv)
   rownames(cormat) <- colnames(cormat) <- colnames(x)
   if(method=="pearson") statistic <- calc.pearson
   else if(method=="spearman") statistic <-  calc.spearman
   else if(method=="kendall") statistic <-  calc.kendall
   else if(method=="dcor") statistic <- calc.dCor
   else if(method=="hoefd") statistic <- calc.hoefd
   else if(method=="mic") statistic <- calc.MIC
   else if(method=="rdc") statistic <- calc.rdc
   else stop("method argümaný pearson, spearman, 
       kendall, dcor, hoefd, mic veya rdc olmalý")
   for(i in 1: (nv-1))
      for(j in i:nv){
       dset <- data.frame(x[,i], x[,j])
        cormat[i,j] <- cormat[j,i] <- statistic(dset)
    }
    if(plot){
       corrplot.mixed(cormat,  lower="number", upper="ellipse",
         lower.col="black", number.cex=0.8)
       mtext(paste("Yöntem:", method), line=2, side=2, col=2)
    }
    (round(cormat, 4))
}

# Fonksiyon 7.30: Çok sayýda korelasyon için hesaplama fonksiyonu
# Baðýmlýlýk – Fonksiyon: 7.9b, 7.10b, 7.11b, 7.16b, 7.22b, 7.23b,
#   7.24, 7.25, 7.26
#
calc.multcor <- function(dset, indices){
   if(missing(indices)) indices <- 1:nrow(dset)
   corpear <- calc.pearson(dset, indices)
   corspear <- calc.spearman(dset, indices)
   corkend <- calc.kendall(dset, indices)
   cordCor <- calc.dCor(dset, indices)
   corhoefd <- calc.hoefd(dset, indices)
   corMIC <- calc.MIC(dset, indices)
   corrdc <- calc.rdc(dset, indices)
   corace <- calc.acecor(dset, indices)
   cornlcor <- calc.nlcor(dset, indices)
   cornnscor <- calc.nnscor(dset, indices)
   return(c(corpear, corspear, corkend, cordCor, corhoefd, 
     corMIC, corrdc, corace, cornlcor, cornnscor))
}

# Fonksiyon 7.31: Farklý örneklem büyüklüðünde
# korelasyon yöntemlerini karþýlaþtýrma
# Baðýmlýlýk – Paket: boot, MASS
#
cor.perf <- function(statistic, n=c(30), r=c(0.2, 0.8), m=10,
   R=2000, dplot=FALSE){
   k <- length(n)
   l <- length(r)
   cormean <- numeric(m)
   corstdev <- numeric(m)
   corperci <- matrix(NA, nrow=m, ncol=2)
   corbias <- numeric(m)
   cormeans <- matrix(0, nrow=k, ncol=l)
   corbiases <- matrix(0, nrow=k, ncol=l)
   corstdevs <- matrix(0, nrow=k, ncol=l)
   corpercis1 <- matrix(0, nrow=k, ncol=l)
   corpercis2 <- matrix(0, nrow=k, ncol=l)
   rownames(cormeans) <- paste0("n=",n)
   colnames(cormeans) <- paste0("r=", r)
   rownames(corbiases) <- rownames(corstdevs) <-
    rownames(cormeans)
   rownames(corpercis1) <- rownames(corpercis1) <-
    rownames(cormeans)
   colnames(corbiases) <- colnames(corstdevs) <-
    colnames(cormeans)
   colnames(corpercis1) <- colnames(corpercis1) <-
    colnames(cormeans)
   for(i in 1:k){
     for(j in 1:l){
       for(z in 1:m){
          #anakitle korelasyonlarý
          covmat <- matrix(c(1, r[j], r[j], 1), ncol=2) 
          vset <- MASS::mvrnorm(n[i], 
          mu=c(0,0), Sigma=covmat, empirical=TRUE)
          plot(vset, pch=21, col=3, bg="gray90", 
          xlab="V.1", ylab="V.2",
          main=paste0("n=", n[i], " & r=", r[j]),
          sub=paste0("Örneklem#: ",z))
          corboot <- boot::boot(vset, statistic=statistic, R=R)
          cormean[z] <- mean(corboot$t)
          corstdev[z] <- sd(corboot$t)
          corperci[z,] <- quantile(corboot$t, 
            probs=c(0.025, 0.975))
          corbias[z] <- cormean[z]-r[j]
        }
        cormeans[i,j] <- round(mean(cormean, na.rm=TRUE),2)
        corbiases[i,j] <- round(mean(corbias, na.rm=TRUE),2)
        corstdevs[i,j] <- round(mean(corstdev, na.rm=TRUE),2)
        corpercis1[i,j] <- round(mean(corperci[,1], 
         na.rm=TRUE),2)
        corpercis2[i,j] <- round(mean(corperci[,2], 
         na.rm=TRUE),2)
     }
   }
   return(list(cormeans=cormeans, corbiases=corbiases, 
     corstdevs=corstdevs, corpercis1=corpercis1,
     corpercis2=corpercis2))
}

# Fonksiyon 7.32: Çok sayýda korelasyon için hesaplama fonksiyonu
# Baðýmlýlýk – Fonksiyon: 7.9b, 7.10b, 7.11b, 7.16b, 7.22b, 
#   7.23b, 7.24
#
calc.multcor2 <- function(dset, indices){
   if(missing(indices)) indices <- 1:nrow(dset)
   corpear <- calc.pearson(dset, indices)
   corspear <- calc.spearman(dset, indices)
   corkend <- calc.kendall(dset, indices)
   cordCor <- calc.dCor(dset, indices)
   corMIC <- calc.MIC(dset, indices)
   corrdc <- calc.rdc(dset, indices)
   corace <- calc.acecor(dset, indices)
   return(c(corpear, corspear, corkend, cordCor, corMIC,
     corrdc, corace))
}

# Fonksion 7.33: Korelasyon için permütasyon testi fonksiyonu 1
permute.cor <- function(dset, statistic,
   alternative="two.sided", R=3000){
   if(!is.data.frame(dset) | !is.matrix(dset)) 
     dset <- as.data.frame(dset)
   n <- nrow(dset)  # Örneklem büyüklüðü
   teta0 <- statistic(dset, indices=1:n) # Orijinal kestirim
   tetas <- c()  # Permütasyon kestirimleri vektörü
   perdset <- matrix(nrow=n, ncol=2)
   perdset[,1] <- dset[,1]
   for (i in 1:R){
     perdset[,2] <- sample(dset[,2], n, replace=FALSE)
     teta <- statistic(perdset, indices=1:n)
     tetas <- c(tetas, teta)
   }
   if(alternative == "less"){
     pv <- sum(tetas < teta0)/R
   }else if(alternative == "greater"){
     pv <- sum(tetas > teta0)/R
   }else{
     pv <- sum(abs(tetas) >= abs(teta0))/R
   } 
   return(list(t0=teta0, t=tetas, p.value=pv,  
     n=n, R=R, call=match.call()))
}

# Fonksiyon 7.34: Permütasyon kestirimleri daðýlýþ grafiði
permute.plot <- function(perres){
    hist(perres$t, col="gray", breaks=50, prob=TRUE, 
      xlab = "r*", ylab="Yoðunluk", main= "Örnekleme Daðýlýþý") 
    lines(density(perres$t), lwd=2, col=4)
    abline(v=0, col=2)
}

# Fonksiyon: 7.35: Korelasyon için permütasyon testi fonksiyonu 2
permute.cor2 <- function(dset, statistic,
   alternative="two.sided", R=1000){
   if(!is.data.frame(dset) | !is.matrix(dset)) 
     dset <- as.data.frame(dset)
   n <- nrow(dset)  # Örneklem büyüklüðü
   teta0 <- statistic(dset, indices=1:n) # Orijinal istatistik
   tetas <- replicate (R, expr=statistic(cbind(dset[,1],
     sample(dset[,2])), indices=1:n))
   if(alternative == "less"){
     pv <- sum(tetas < teta0)/R
   }else if(alternative == "greater"){
     pv <- sum(tetas > teta0)/R
   }else{
     pv <- sum(abs(tetas) >= abs(teta0))/R
   } 
   return(list(t0=teta0, t=tetas, p.value=pv,  n=n, R=R, 
   call=match.call()))
}

# Fonksiyon: 7.36: Korelasyon için permütasyon testi fonksiyonu 3
permute.cor3 <- function(dset, statistic,
    alternative="two.sided", R=10000){
   if(!is.data.frame(dset) | !is.matrix(dset)) 
     dset <- as.data.frame(dset)
   n <- nrow(dset)  # Örneklem büyüklüðü
   teta0 <- statistic(dset, indices=1:n) # Orijinal istatistik
    tetas <- sapply(1:R,FUN=function(i) statistic(cbind(dset[,1],
     sample(dset[,2])), indices=1:n))
   if(alternative == "less"){
     pv <- sum(tetas < teta0)/R
   }else if(alternative == "greater"){
     pv <- sum(tetas > teta0)/R
   }else{
     pv <- sum(abs(tetas) >= abs(teta0))/R
   } 
   return(list(t0=teta0, t=tetas, p.value=pv,  n=n, R=R, 
   call=match.call()))
}

# Fonksiyon 7.37: Jackknife ile korelasyon kestirimleri
jack.cortest <- function (dset, statistic,
    alternative="two.sided"){
   if(!is.data.frame(dset) | !is.matrix(dset)) 
    dset <- as.data.frame(dset)
    n <- nrow(dset)  # Örneklem büyüklüðü
    teta0 <- statistic(dset, indices=1:n) 
    tetas <- rep(0, n)
    for(i in 1:n) {
      #Jackknife kestirimleri
      tetas[i] <- statistic(dset, indices=(1:n)[-i]) 
    }
    pseudos <- n*teta0-(n-1)* tetas
    pseudos.mean <- mean(pseudos)
    pseudos.se <- sqrt(sum((pseudos-pseudos.mean)^2)/(n*(n-1)))
    pseudos.bias <- (n-1)*(teta0-mean(tetas)) 
    if(alternative == "less"){
      pv <- sum(tetas > 0)/n
    }else if(alternative == "greater"){
      pv <- sum(tetas < 0)/n
    }else{
      pv <- sum(abs(tetas == 0))/n
    } 
    return(list(
      t0 = teta0,       
      t = tetas, 
      pseudos.mean = pseudos.mean,     
      pseudos.bias = pseudos.bias,
      pseudos.se = pseudos.se,            
      pseudos = pseudos,
      data = dset,
      p.value = pv,
      call = match.call()))
}

#### BÖLÜM 8 ###################################################################

# Fonksiyon 8.1: Doðrusal regresyon katsayýlarý hesaplama 
calc.linreg <- function(x, formula, indices){
   if(missing(indices)) indices <- 1:nrow(x)
   model <- lm(formula, data=x[indices,])
   return(coef(model))
}

# Fonksiyon 8.2. R2 (R-kare) hesaplama fonksiyonu
calc.R2 <- function(x, formula, indices){
   if(missing(indices)) indices <- 1:nrow(x)
   model <- lm(formula, data=x[indices,])
   summodel <- summary(model)
   return(summodel$r.squared)
}

# Fonksiyon 8.3: Çoklu doðrusal regresyon katsayýlarýný hesaplama 
#
> calc.multreg <- function(x, formula, indices){
   if(missing(indices)) indices <- 1:nrow(x)
    model <- lm(formula, data=x[indices,])
   (coef(model))
}

# Fonksiyon 8.4. Basit Doðrusal regresyon için 
# çeþitli istatistikler
calc.mslinreg <- function(x, formula, indices){
   if(missing(indices)) indices <- 1:nrow(x)
   model <- lm(formula, data=x[indices,])
   rmse <- sqrt(mean(model$residuals)^2)
   summodel <- summary(model)
   return(e(summodel$coefficients[,1], summodel$r.squared,
     summodel$adj.r.squared, summodel$fstatistic[1],
     summodel$sigma, rmse))
}

# Fonksiyon 8.5: Kantil regresyonu katsayýlarý fonksiyonu
calc.kantreg <- function(x, formula, indices, tau=0.5){
if(!require(quantreg)) 
    {install.packages("quantreg"); library(quantreg)}
   if(missing(indices)) indices <- 1:nrow(x)
   kantregmodel <- rq(formula, data=x[indices,], tau=tau)
   return(coef(kantregmodel))  
}

# Fonksiyon 8.6: LOOCV ile çekirdek regresyonu tahmin hatalarý
loocv <- function(x, y, h=1.0){
   n <- length(x)
   tahmin.hata <- rep(NA, n)
   for(i in 1:n){
     # Test verisi
     x.test <- x[i] ; y.test <- y[i]
     # Eðitim verisi
     x.egitim <- x[-i] ; y.egitim <- y[-i]
     # Çekirdek regresyonu uygula
     y.tahmin <- ksmooth(x=x.egitim, y=y.egitim, 
        kernel="normal", bandwidth=h, x.points = x.test)
     # Tahmin hatasýný hesapla
     tahmin.hata[i] <- (y.test-y.tahmin$y)^2
   }
   #Tahmin hata ortalamasý
   return(mean(tahmin.hata, na.rm=TRUE))
}

# Fonksiyon 8.7: Lojistik regresyon katsayýlarýný hesaplama 
calc.logreg <- function(x, formula, indices){
   if(missing(indices)) indices <- 1:nrow(x)
   logmodel <- glm(formula, data=x[indices,], family="binomial")  
   return(coef(logmodel))  
}

# Fonksiyon 8.8: Chapman-Richard büyüme modeli fonksiyonu
chapricgm <- function(t, A, b, c) A*(1-exp(-b*t))^c

# Fonksiyon 8.9: Büyüme modeli katsayýlarý hesaplama fonksiyonu
calc.nlsreg <- function(x, formula, indices){
   if(missing(indices)) indices <- 1:nrow(x)
   nlsmodel <- nls(formula, data=x[indices,])
   return(coef(nlsmodel))  
}

# Fonksiyon 8.10: GAM modeli katsayýlarý hesaplama fonksiyonu
# Baðýmlýlýk – Paket: mgcv
calc.gamcoef <- function(x, formula, indices, 
   family=gaussian(), link=NULL){
   if(!require(mgcv)) 
      {install.packages("mgcv"); library(mgcv)}
   if(missing(indices)) indices <- 1:nrow(x)
   gammodel <- gam(formula, data=x[indices,], family=family)
   return(coef(gammodel))  
}

# Fonksiyon 8.11: Eklem (hinge) fonksiyonu grafiði fonksiyonu
hingeplot <- function(x, knot, d=1){
    x <- sort(x)
    t <- x[knot]
    h1 <- h2 <- rep(0,length(x)) 
    idx <- which(x>t)
    h1[idx] <- (x[idx]-t)^d
    h2[-idx] <- (t-x[-idx])^d
    plot(x, h1, col=2, type="l", lty=1, lwd=2,  
      ylab="h(x)")
    par(new=TRUE)
    plot(x, h2, col=4, type="l",lty=2, lwd=2, ylab="", yaxt="n")
    legend("top", legend=c("(x-t)+","(t-x)+"), col=c(2,4),
      lty=c(1,2), lwd=c(2,2), bty="n", horiz=TRUE)
}


#### BÖLÜM 9 ########################################################

# Fonksiyon 9.1a: Eðitim/Doðrulama/Test verisi bölümleme
split.swor <- function(dset, tp=.8){
   n <- nrow(dset)
   indices <- sample(1:n, size=n, replace=F) # Veriyi karýþtýr
   dset <- dset[indices,]
   ntr <- ceiling(n*tp)
   nvl <- ceiling(n*(1-tp)/2)
   train.set <- dset[1:ntr,]
   val.set <- dset[(ntr+1) : (ntr+nvl+1),]
   test.set <- dset[(ntr+nvl+2):n,]
   sets <- list(train=train.set, valid=val.set, test=test.set)
   return(sets)
}

# Fonksiyon 9.1b: Tabakalý bölümleme
split.strata <- function(dset, strata, vs=TRUE, alg=1, tp=.7){
   n <- nrow(dset)
   indices <- sample(1:n, size=n, replace=F) # Veriyi karýþtýr
   dset <- dset[indices,]
 # rminer ile bölümleme
   if(alg==1){
     if(!require(rminer)) 
       {install.packages("rminer"); require(rminer)}
     H <- holdout(y=dset[,names(dset)[strata]], ratio=tp,
       internalsplit=vs, mode="stratified")
     tr <- dset[H$tr,]
     ts <- dset[H$ts,]
     tv <- NULL
     if(!is.null(H$val)){
        tr <- dset[H$itr,]
        tv <- dset[H$val,]
     }
   }
   else{
 # rsample ile bölümleme
     if(!require(rsample)) 
       {install.packages("rsample"); require(rsample)}
     H <-initial_split(dset,strata=names(dset)[strata],prop=tp)
     tr <- training(H)
     ts <- testing(H)
     tv <- NULL
   }
   (list(train=tr, valid=tv, test=ts))
}

# Fonksiyon 9.2: Bölünmüþ Örneklemle Çapraz Doðrulama
split.cv <- function(dset, formula, tp=0.8, seed=NULL){
    if(tp<0.5 | tp>0.8) stop("tp oraný 0.5-0.9 olmalý")
    n <- nrow(dset) ; r <- ncol(dset); k <- r-1
    if(!is.null(seed)) set.seed(seed)  
    trainidx <- sample(1:n, size=ceiling(n*tp), replace=FALSE)
    trainset <- dset[trainidx,]
    valset <- dset[-trainidx,]
    model <- lm(formula, data=trainset) # Modeli eðit
    iv <- toString(formula(model)[[2]]) # Baðýmlý deðiþkeni bul
    valpreds <- predict(model, valset)  # Tahminleri hesapla
    errors <- valset[,iv]-valpreds      # Hatalarý hesapla
    ess <- sum(errors^2) # Hata kareler toplamý
    rss  <- sum((valpreds-mean(dset[,iv]))^2) # Reg.Kare.Top.
    tss <- ess + rss  # Genel kareler toplamý
    #tss <- sum((valset[,iv]-mean(dset[,iv]))^2) # Altrntv.GKT
    R2 <- 1-ess/tss 
    adjR2 <- 1-((n-1)/(n-k-1))*(1-R2)
    RMSE <- sqrt(mean(errors^2))
    MAE <- mean(abs(errors))
    PER <- round(RMSE/mean(valset[,iv])*100, 2)
    return(data.frame(R2=R2,adjR2=adjR2, 
      RMSE=RMSE, MAE=MAE, PER=PER))
}
# Fonksiyon 9.3: Monte Carlo Çapraz Doðrulamasý
montecarlo.cv <- function(dset, formula, tp=0.8, 
   R=1000, seed=NULL){
   if(!is.null(seed)) set.seed(seed)  # Yeniden üretilebilirlik 
   if(tp < 0.5 | tp > 0.9) stop("tp oraný 0.5-0.9 olmalý")
   n <- nrow(dset)
   lmmod <- lm(formula, data=dset)
   iv <- toString(formula(lmmod)[[2]]) 
   rmse <- mae <- mape <- numeric(R )
   for(i in 1:R) {
      ntr <- floor(n*tp)
      trainidx <- sample(1:n, ntr, replace=FALSE)
      trainmodel <- lm(formula, data=dset[trainidx,])
      ypreds <- predict(trainmodel, newdata=dset[-trainidx,])
      mae[i] <- mean(abs(dset[-trainidx,iv]-ypreds))
      rmse[i] <- sqrt(mean((dset[-trainidx,iv]-ypreds)^2))
      mape[i] <- mean(abs(dset[-trainidx,iv]-ypreds)/
        dset[-trainidx,iv])*100
    }
    MAE <- mean(mae, na.rm=TRUE)
    MAPE <- mean(mape, na.rm=TRUE)
    RMSE <- mean(rmse, na.rm=TRUE)
    return(data.frame(RMSE=RMSE, MAE=MAE, MAPE=MAPE))
}

# Fonksiyon 9.4a: Birini atarak çapraz doðrulama (LOOCV)
#
loo.cv <- function(dset, formula){
    n <- nrow(dset); r <- ncol(dset); k <- r-1
    hatmat <- lm.influence(lm(formula, data=dset))$hat
    errors <- numeric(n); hatvals <- numeric(n)
    for(i in 1:n){
      trainset <- dset[-i,]
      valset <- dset[i,]
      model <- lm(formula, data=trainset) # Modeli eðit
      iv <- toString(formula(model)[[2]]) # Baðýmlý deðiþken
      valpreds <- predict(model, valset)  # Tahminleri hesapla 
      errors[i] <- (valset[,iv]-valpreds) # Tahmin hatalarý
      hatvals[i] <- hatmat[i]		  # Þapka deðeri
    }
    ess <- sum(errors^2) # Hata kareler toplamý
    tss <- sum((dset[,iv]-mean(dset[,iv]))^2)  # Genel K.T.
    R2 <- 1-ess/tss
    adjR2 <- 1-((n-1)/(n-k-1))*(1-R2)
    RMSE <- sqrt(mean(errors^2))
    MAE <- mean(abs(errors))
    MAPE = mean(abs(errors)/dset[,iv])*100
    PER <- round(RMSE/mean(dset[,iv])*100, 2)
    PRESS <- sum((errors/(1-hatvals))^2)
    CV <- sqrt(PRESS/length(errors))
    resdf <- data.frame(R2=R2, adjR2=adjR2, RMSE=RMSE, MAE=MAE,
      MAPE=MAPE, PER=PER, PRESS=PRESS, CV=CV)
    return(resdf)
}

# Fonksiyon 9.4b: caret paketi ile LOOCV
# Baðýmlýlýk – Paket: caret
#
loo.cv2 <- function(data, formula, method="lm", seed=NULL){
   if(!require(caret, quietly=TRUE)){
     install.packages("caret",
      repo="https://cloud.r-project.org/")
     require(caret, quietly=TRUE)
   }
   if(!is.null(seed)) set.seed(seed)
   tkontrol <- trainControl(method="LOOCV") #Eðitim ayarlarý
   model <- train(formula=formula, data=data, 
     method=method, trControl=tkontrol)
   modelres <- unname(model$results) # Performans sonuçlarý
   return(data.frame(R2=modelres[3], RMSE=modelres[2],
     MAE=modelres[4]))
}

# Fonksiyon 9.5a: Kat oluþturma fonksiyonu 1
gen.folds <- function(indices, k){
    n <- length(indices)
    if(k<2 | k>n) stop(paste0("k, [2,",n, "] aralýðýnda olmalý"))
    foldsize <- rep(n%/%k,k)
    nrem <- n-k*n%/%k
    if(nrem!=0){
      rfold <- sample(1:k, nrem, replace=FALSE)
      foldsize[c(rfold)] <- foldsize[c(rfold)]+1
    }
    foldindices <- list()
    foldindices[[1]] <- indices[1: foldsize[1]]
    for(i in 2:k)
      foldindices[[i]] <- indices[(cumsum(foldsize)[i-1]+1):
        (cumsum(foldsize)[i])]
    return(list(foldsize=foldsize, foldindices=foldindices))
}

# Fonksiyon 9.5b: Kat oluþturma fonksiyonu 2
gen.folds2 <- function(indices, k){
    n <- length(indices)
    if(k<2 | k>n) stop(paste0("k, [2,",n, "] aralýðýnda olmalý"))
    folds <- cut(1:n, breaks=k, labels=FALSE)
    foldsize <- table(folds)
    foldindices <- list()
    foldindices[[1]] <- indices[1: foldsize[1]]
    for(i in 2:k)
      foldindices[[i]] <- indices[(cumsum(foldsize)[i-1]+1):
        (cumsum(foldsize)[i])]
    return(list(foldsize=foldsize, foldindices=foldindices))
}

# Fonksiyon 9.6a: K-katlý çapraz doðrulama 
# Baðýmlýlýk – Fonksiyon: 9.5a
#
kfold.cv <- function(dset, formula, k=5, seed=NULL){
    n <- nrow(dset); np <- ncol(dset)-1
    if(!is.null(seed)) set.seed(seed)  
    indices <- sample(1:n,size=n,replace=FALSE) #Veriyi karýþtýr
    folds <- gen.folds(indices, k)
    foldsize <- folds$foldsize
    foldindices <- folds$foldindices
    R2 <- adjR2 <- RMSE <- MAE <- PER <- numeric(k)
    for(i in 1:k){
      valset <- dset[foldindices[[i]],]
      trainset <- dset[-foldindices[[i]],]
      model <- lm(formula, data=trainset) # Modeli tanýmla
 # Baðýmlý deðiþkeni sapta
      iv <- toString(formula(model)[[2]]) 
      valpreds <- predict(model, newdata=valset) 
 # Test tahminleri hesapla
      errors <- valset[,iv]-valpreds
 # Hata kareler toplamý 
      ess <- sum(errors^2)  
 # Genel kareler toplamý
      tss <- sum((valset[,iv]-mean(trainset[,iv]))^2)  
 # Regresyon kareler toplamý
      rss  <- sum((valpreds-mean(dset[,iv]))^2) 
 # Performans ölçütlerini hesapla 
      R2[i] <- 1-ess/(ess+rss)    
      adjR2[i] <- 1-((n-1)/(n-np-1))*(1-R2[i])
      RMSE[i] <- sqrt(mean(errors^2))
      MAE[i] <- mean(abs(errors))
      PER[i] <- round(RMSE[i]/mean(valset[,iv])*100, 2)
    }
    kR2 <- mean(R2)
    kadjR2 <- mean(adjR2)
    kRMSE <- mean(RMSE)
    kMAE <- mean(MAE)
    kPER <- mean(PER)
    sonuc <- list(foldsize=foldsize, foldindices=foldindices,
      results=data.frame(R2=kR2, adjR2=kadjR2, 
      RMSE=kRMSE, MAE=kMAE, PER=kPER))
    return(sonuc)
}

# Fonksiyon 9.6b: caret paketi ile k-katlý çapraz doðrulama 
# Baðýmlýlýk – Paket: caret
#
kfold.cv2 <- function(data, formula, method="lm",
   k=5, seed=NULL){
   if(!require(caret, quietly=TRUE)){
     install.packages("caret",
     repo="https://cloud.r-project.org/")
     require(caret, quietly=TRUE)
   }
   if(!is.null(seed)) set.seed(seed)
   # Eðitim ayarlarý
   tkontrol <- trainControl(method="CV", number=k) 
   # Modelin eðitimi ve testi
   model <- train(formula, data, method, trControl=tkontrol)
   modelres <- unname(model$results) # Baþarým ölçütleri
   return(data.frame(R2=modelres[3], 
      RMSE=modelres[2], MAE=modelres[4]))
}

# Fonksiyon 9.7a: Tekrarlanan k-katlý çapraz doðrulama 
# Baðýmlýlýk – Fonksiyon: 9.6a
#
repkfold.cv <- function(data, formula, k=5, rep=5, seed=NULL){
  n <- nrow(data)
  foldindices  <- list()
  results <- data.frame()
  for(i in 1:rep){
    if(!is.null(seed)) set.seed(seed)  
    indices <- sample(1:n, size=n, replace=FALSE)
    data <- data[indices,]
    kfres<- kfold.cv(data, formula, k=k, seed=seed)
    results <- rbind(results, kfres$results)
    foldsize <- kfres$foldsize
    foldindices[[i]] <- kfres$foldindices
  }
  R2 <- mean(results[,1])
  adjR2 <- mean(results[,2])
  RMSE <- mean(results[,3])
  MAE <- mean(results[,4])
  PER <- mean(results[,5])
  reslist <- list(foldsize=foldsize, foldindices=foldindices,
    results=data.frame(R2=R2, adjR2=adjR2, 
    RMSE=RMSE, MAE=MAE, PER=PER))
  (reslist)
}

# Fonksiyon 9.7b: caret paketi ile tekrarlanan k-CV 
# Baðýmlýlýk – Paket: caret
#
repkfold.cv2 <- function(data, formula, method="lm", 
   k=5, rep=5, seed=NULL){
   if(!require(caret, quietly=TRUE)){
     install.packages("caret",
       repo="https://cloud.r-project.org/")
     require(caret, quietly=TRUE)
   }
   if(!is.null(seed)) set.seed(seed)
   tkontrol <- trainControl(method="repeatedcv", 
     number=k, repeats=rep) 
   model <- train(formula, data, method, trControl=tkontrol)
   modelres <- unname(model$results) # Performans sonuçlarý
   return(data.frame(R2=modelres[3], RMSE=modelres[2],
     MAE=modelres[4]))
}

# Fonksiyon 9.8: caret ile bootstrap 
# Baðýmlýlýk – Paket: caret
#
boot.cv <- function(dset, formula, method="lm", 
   R=1000, seed=NULL){
   if(!require(caret, quietly=TRUE)){
     install.packages("caret",
       repo="https://cloud.r-project.org/")
     require(caret, quietly=TRUE)
   }
   if(!is.null(seed)) set.seed(seed)
   tkontrol <- trainControl(method="boot", number=R) 
   model <- train(formula=formula, data=dset, method=method,
     trControl=tkontrol)
   modelres <- unname(model$results) # Performans sonuçlarý
   return(data.frame(R2=modelres[3], RMSE=modelres[2],
     MAE=modelres[4]))
}

# Fonksiyon 9.9: Bootstrap örnekleminde benzersiz gözlem oraný
test632 <- function(x, R=10000, plot=TRUE){
    nux <- function(x) length(unique(x))
    nx <- length(x)
    rsampx <- matrix(sample(x, size=nx*R, replace=TRUE),
      ncol=nx, nrow=R)
    uniquex <- apply(rsampx, 1, nux)/nx
    if(plot) hist(uniquex, col="gray", main="Bootstrap 
       örnekleminde eþi olmayan gözlem sayýsý")
    (mean(uniquex))
}


# Fonksiyon 9.10: bootstrap .632 ile CV 
boot.632 <- function(dset, formula, R=1000) {
   n <- nrow(dset)
   lmmod <- lm(formula, data=dset)
   iv <- toString(formula(lmmod)[[2]]) 
   errors <- lmmod$residuals
   MSE <- mean(errors^2)
   MAE <- mean(abs(errors))
   APE <- mean(abs(errors)/dset[,iv])
   y <- lmmod$model[,iv]
   x <- subset(lmmod$model, colnames(lmmod$model)==c(iv))
   bs.sqrerr <- bs.abserr <- bs.pererr <- numeric(R )
   for (i in 1:R){
     resmp <- sample(1:n, n, replace=T)
     residx <- sort(unique(resmp))
     bs.model <- lm(formula, data=dset[resmp,])
     ypreds <- predict(bs.model,newdata=dset[-residx,])
     bs.err <- y[-residx]-ypreds
     bs.sqrerr[i] <- mean(bs.err^2)
     bs.abserr[i] <- mean(abs(bs.err))
     bs.pererr[i] <- mean(abs(bs.err)/y[-residx])
   }
   MSEopt <- (mean(bs.sqrerr)-MSE)* 0.632
   MAEopt <- (mean(bs.abserr)-MAE)* 0.632
   APEopt <- (mean(bs.pererr)-APE)* 0.632
   RMSE <- sqrt(MSE + MSEopt)
   MAE <- MAE + MAEopt
   MAPE <- (APE + APEopt)*100
   return(data.frame(RMSE=RMSE, MAE=MAE, MAPE=MAPE))  
}

#### BÖLÜM 10 ##################################################

# Fonksiyon 10.1: Boþ model kalýntýlarýný hesaplama
calc.anova <- function(dset, formula){
    model <- lm(formula, data=dset)
    vat <- anova(model)
    F <- vat[toString(formula(model)[[3]]), "F value"]
 # Boþ model tanýmla 
    bosmodel <- lm(as.formula(paste0(formula(model)[[2]],"~1")),
      data=dset)
 # Boþ model ile tahminler 
    ytahmin <- fitted(bosmodel)   
 # Boþ model kalýntýlar 
    bosmodel.kalinti <- residuals(bosmodel)             
 # Düzeltilmiþ kalýntýlar
    cv <- bosmodel.kalinti / sqrt(1-hatvalues(bosmodel)) 
   return(list(formula=formula, ytahmin=ytahmin, 
      cv=cv, model.F=F))
}

# Fonksiyon 10.2: F deðerini hesaplama fonksiyonu
model.anovaF <- function(data, indices, formula, ytahmin, cv){
    if(missing(indices)) indices <- 1:nrow(data)
    model <- lm(formula, data=data)
    ystar <- ytahmin + cv[indices]
    formula2 <- as.formula(paste0("ystar~",formula(model)[[3]]))
    vat <- anova(lm(formula2, data=data))
    return(vat[toString(formula(model)[[3]]), "F value"])
}

# Fonksiyon 10.3: ANOVA Bootstrap F kestirimleri grafiði
fcdf.plot <- function(bootres){
    plot(bootres$t, ecdf(bootres$t)(bootres$t), 
      col=2, bg="gray90", cex=1.1,pch=21, 
      xlab="f,f*", ylab="Pr(F <= f)",
      main="Teorik (F) ve Ampirik (F*) CDF Grafiði")
    nlev <- length(unique(bsAnova1$data[,2]))
    curve(pf(x, nlev-1, sum(nrow(bootres$data))-nlev),
      lwd=2, col=4, add=TRUE)
    legend("bottomright", pch=c(21, NA), lty=c(NA, 1), 
      lwd=c(NA, 2), col=c(2,4), bty="n", 
      legend=c("Ampirik F CDF","Teorik F CDF"))
}

# Fonksiyon 10.4: ANOVA önemlilik olasýlýk deðeri hesaplama
anova.pval <- function(bootres, model.F){
    p.value  <- (sum(bootres$t >= model.F) + 1) / 
      (length(bootres$t) + 1)
   return(p.value)
}

# Fonksiyon 10.5: ANOVA Wild bootstrap F hesaplama fonksiyonu
wild.anovaF <- function(data, indices, formula, ytahmin, cv, uopt=1){
   #if(missing(indices)) indices <- 1:nrow(data)
   n <- length(indices)
   model <- lm(formula, data=data)
   if(uopt==1){
      Ur <- sample(c(-1, 1), size=n, replace=TRUE, 
         prob=c(0.5, 0.5))
      ystar <- ytahmin + (cv*ruv)[indices] 
   }else{
      Uf <- sample(c(-(sqrt(5)-1)/2, (sqrt(5) + 1)/2),
         size=n, replace=TRUE, prob=c((sqrt(5)+ 1 / (2*sqrt(5))),
        (sqrt(5)-1)/(2*sqrt(5))))
      ystar <- ytahmin + (cv*rfv)[indices] 
    }
    formula2 <- as.formula(paste0("ystar~",formula(model)[[3]]))
    vat <- anova(lm(formula2, data=data))
    return(vat[toString(formula(model)[[3]]), "F value"])
}

# Fonksiyon 10.6: Ýkili karþýlaþtýrma için permütasyon testi
# Baðýmlýlýk – Fonksiyon: 6.3
#
pairwise.permute <- function(x, statistic,
   alternative="two.sided", R=3000){
 # Uzun format kontrolü
   if(!is.factor(x[,2])){
     stop("Veri uzun formatta olmalý")
   }
   grplev <- unique(x[,2])
   pv <- diag(x=1, ncol=length(grplev), nrow=length(grplev))
   colnames(pv) <- rownames(pv) <- grplev
   for(i in 1: (length(grplev)-1)){
     for(j in (i+1):length(grplev)){
        gpair <- droplevels(x[which(x[,2]==
          grplev[i] | x[,2]==grplev[j]),])
        perres <- permute.dif(x=gpair, statistic=statistic, 
          alternative=alternative, R=R)
        pv[i,j] <- perres$p.value
     }
   }
   pv[lower.tri(pv)] <- t(pv)[lower.tri(pv)]
   return(list(p.value=pv,
     alternative=alternative, R=R, call=match.call()))
}
