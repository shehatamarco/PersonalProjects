
rm(list=ls())

load("/Users/marcoshehata/Desktop/UniPD/Magistrale/AnnoII/AnalisiFinanza/Homework/Quantile_Regression/mr1m.RDATA")

RV_funct <- function(matr) as.vector(apply(matr^2, 2, sum))

BPV_funct <- function(matr, n = 78){
  abs.r <- abs(matr)
  abs.rr <- abs.r[2:n,] * abs.r[1:(n-1),]
  mu1 <- sqrt(2/pi)
  BPV <- colSums(abs.rr) * (n/(n-1)) * ((mu1)^(-2))
  return(as.vector(BPV))
}

r <- as.vector(apply(mr1m, 2, sum))
RV <- RV_funct(mr1m)
BPV <- BPV_funct(mr1m)

library(zoo)

n <- length(r)

dati <- data.frame(
  r = r[2:n],
  RV = RV[1:(n-1)],
  BPV = BPV[1:(n-1)],
  RVW = rollmean(RV, 5, align = "right", fill = NA)[1:(n-1)],
  RVM = rollmean(RV, 22, align = "right", fill = NA)[1:(n-1)],
  BPVW = rollmean(BPV, 5, align = "right", fill = NA)[1:(n-1)],
  BPVM = rollmean(BPV, 22, align = "right", fill = NA)[1:(n-1)]
)

dati <- dati[complete.cases(dati),]

head(dati)

# HARQ QR: BAckCasting

alpha <- 0.01 # ordine del VAR 
n.info <- 300 # numero di unità usate per stima

library(quantreg)

# BackCasting

v.har.qr <- NULL

for(i in 1:(nrow(dati) - n.info)){
  
  tmp.used <- dati[i:(n.info+i-1),]
  
  QR.mod <- rq(
    r ~ RV + RVW + RVM, 
    data = tmp.used,
    tau = alpha
  )
  
  VaR.est <- predict(QR.mod, newdata = dati[(n.info+i),])
  
  v.har.qr <- c(v.har.qr, VaR.est)
  
  cat("Done:", i, "\n")
  
}

# HARQ: REQ

v.har.req <- NULL

library(extremefit)

VaR.REQ.H <- function(z.est, q.est, a = 0.01, q = round(0.1*length(z.est)) ){
  
  options(warn = -1)
  
  r <- q/length(z.est)
  xi <- extremefit::hill(z.est)$hill[q]
  zk <- extremefit::hill(z.est)$xsort[q]
  res <- q.est * zk * (r/a)^(xi)
  
  options(warn = 0)
  
  return(res)
}

for(i in 1:(nrow(dati) - n.info)){
  
  tmp.used <- dati[i:(n.info+i-1),]
  
  QR.mod <- rq(
    r ~ RV + RVW + RVM, 
    data = tmp.used,
    tau = alpha
  )
  
  q.t <- predict(QR.mod, newdata = tmp.used)
  
  q.new <- predict(QR.mod, newdata = dati[(n.info+i),]) # previsione quantile 1 passo in avanti
  
  z.t <- tmp.used$r / q.t # calcolo i quantili depurati
  
  VaR.est <- VaR.REQ.H(z.t, q.new, a = alpha) # stimo il var
  
  v.har.req <- c(v.har.req, as.numeric(VaR.est)) # lo salvo
  
  cat("Done:", i, "\n")
  
}

load("/Users/marcoshehata/Desktop/UniPD/Magistrale/AnnoII/AnalisiFinanza/Homework/Quantile_Regression/VaRQRGEVREQ.RDATA")

head(tab.res)

nrow(tab.res) - length(v.har.qr)

NA.s <- rep(NA, 21)

v.har.qr <- c(NA.s, v.har.qr)
v.har.req <- c(NA.s, v.har.req)

tab.res$HAR.QR <- v.har.qr
tab.res$HAR.REQ <- v.har.req

head(tab.res)

tab.res <- tab.res[complete.cases(tab.res),]

tab.res

# Analisi Grafica

par(mfrow = c(1,2))

plot(tab.res$r, type = "l", ylim = c(-0.3, 0.15), xlab = "t", ylab = "rendimenti", main = "HAR QR vs QR")
lines(tab.res$QR, col = 2, type = "l")
lines(tab.res$HAR.QR, col = 3, type = "l")

plot(tab.res$r, type = "l", ylim = c(-0.3, 0.15), xlab = "t", ylab = "rendimenti", main = "HAR REQ vs REQ")
lines(tab.res$REQ, col = 2, type = "l")
lines(tab.res$HAR.REQ, col = 3, type = "l")

alpha*nrow(tab.res)
sum(tab.res$r < tab.res$GEV)
sum(tab.res$r < tab.res$QR)
sum(tab.res$r < tab.res$REQ)
sum(tab.res$r < tab.res$HAR.QR)
sum(tab.res$r < tab.res$HAR.REQ)

# C'è una cosa che non mi torna: proviamo a diminuire q

v.har.req <- NULL

for(i in 1:(nrow(dati) - n.info)){
  
  tmp.used <- dati[i:(n.info+i-1),]
  
  QR.mod <- rq(
    r ~ RV + RVW + RVM, 
    data = tmp.used,
    tau = alpha
  )
  
  q.t <- predict(QR.mod, newdata = tmp.used)
  
  q.new <- predict(QR.mod, newdata = dati[(n.info+i),]) # previsione quantile 1 passo in avanti
  
  z.t <- tmp.used$r / q.t # calcolo i quantili depurati
  
  VaR.est <- VaR.REQ.H(z.t, q.new, a = alpha, q = round(0.05*length(z.t)) ) # stimo il var
  
  v.har.req <- c(v.har.req, as.numeric(VaR.est)) # lo salvo
  
  cat("Done:", i, "\n")
  
}

par(mfrow=c(1,2))

plot(tab.res$r, type = "l", ylim = c(-0.3, 0.15), xlab = "t", ylab = "rendimenti", main = "REQ: k = 30")
lines(tab.res$REQ, col = 2, type = "l")
lines(tab.res$HAR.REQ, col = 3, type = "l")

plot(tab.res$r, type = "l", ylim = c(-0.3, 0.15), xlab = "t", ylab = "rendimenti", main = "REQ: K = 15")
lines(tab.res$REQ, col = 2, type = "l")
lines(v.har.req, col = 3, type = "l")

sum(tab.res$r < v.har.req)
alpha*nrow(tab.res)
sum(tab.res$r < tab.res$GEV)
sum(tab.res$r < tab.res$QR)
sum(tab.res$r < tab.res$REQ)
sum(tab.res$r < tab.res$HAR.QR)
sum(tab.res$r < tab.res$HAR.REQ)

par(mfrow=c(1,1))
plot(tab.res$r, type = "l", ylim = c(-0.3, 0.15), xlab = "t", ylab = "rendimenti", main = "Simple QR vs HARQ REQ")
lines(tab.res$QR, col = 2, type = "l")
lines(v.har.req, col = 3, type = "l")

tab.res$HAR.REQ <- v.har.req

head(tab.res)

alpha*nrow(tab.res)
sum(tab.res$r < tab.res$GEV)
sum(tab.res$r < tab.res$QR)
sum(tab.res$r < tab.res$REQ)
sum(tab.res$r < tab.res$HAR.QR)
sum(tab.res$r < tab.res$HAR.REQ)

losses.GEV <- tab.res$GEV
losses.GEV[tab.res$r < tab.res$GEV] <- (1 + (tab.res$r - tab.res$GEV)^2)[tab.res$r < tab.res$GEV]
losses.GEV[tab.res$r >= tab.res$GEV] <- ( (tab.res$r < 0) * (abs(tab.res$GEV) - abs(tab.res$r)) )[tab.res$r >= tab.res$GEV]

losses.QR <- tab.res$QR
losses.QR[tab.res$r < tab.res$QR] <- (1 + (tab.res$r - tab.res$QR)^2)[tab.res$r < tab.res$QR]
losses.QR[tab.res$r >= tab.res$QR] <- ( (tab.res$r < 0) * (abs(tab.res$QR) - abs(tab.res$r)) )[tab.res$r >= tab.res$QR]

losses.REQ <- tab.res$REQ
losses.REQ[tab.res$r < tab.res$REQ] <- (1 + (tab.res$r - tab.res$REQ)^2)[tab.res$r < tab.res$REQ]
losses.REQ[tab.res$r >= tab.res$REQ] <- ( (tab.res$r < 0) * (abs(tab.res$REQ) - abs(tab.res$r)) )[tab.res$r >= tab.res$REQ]

losses.HAR.QR <- tab.res$HAR.QR
losses.HAR.QR[tab.res$r < tab.res$HAR.QR] <- (1 + (tab.res$r - tab.res$HAR.QR)^2)[tab.res$r < tab.res$HAR.QR]
losses.HAR.QR[tab.res$r >= tab.res$HAR.QR] <- ( (tab.res$r < 0) * (abs(tab.res$HAR.QR) - abs(tab.res$r)) )[tab.res$r >= tab.res$HAR.QR]

losses.HAR.REQ <- tab.res$HAR.REQ
losses.HAR.REQ[tab.res$r < tab.res$HAR.REQ] <- (1 + (tab.res$r - tab.res$HAR.REQ)^2)[tab.res$r < tab.res$HAR.REQ]
losses.HAR.REQ[tab.res$r >= tab.res$HAR.REQ] <- ( (tab.res$r < 0) * (abs(tab.res$HAR.REQ) - abs(tab.res$r)) )[tab.res$r >= tab.res$HAR.REQ]

library(forecast)

res.caporin <- c(
  dm.test(losses.HAR.REQ, losses.GEV, alternative = "less", power = 1)$p.value,
  dm.test(losses.HAR.REQ, losses.QR, alternative = "less", power = 1)$p.value,
  dm.test(losses.HAR.REQ, losses.REQ, alternative = "less", power = 1)$p.value,
  dm.test(losses.HAR.REQ, losses.HAR.QR, alternative = "less", power = 1)$p.value
)

# Solo Lopez

losses.GEV <- (tab.res$r < tab.res$GEV) * (1 + (tab.res$r - tab.res$GEV)^2)
losses.QR <- (tab.res$r < tab.res$QR) * (1 + (tab.res$r - tab.res$QR)^2)
losses.REQ <- (tab.res$r < tab.res$REQ) * (1 + (tab.res$r - tab.res$REQ)^2)
losses.HAR.QR <- (tab.res$r < tab.res$HAR.QR) * (1 + (tab.res$r - tab.res$HAR.QR)^2)
losses.HAR.REQ <- (tab.res$r < tab.res$HAR.REQ) * (1 + (tab.res$r - tab.res$HAR.REQ)^2)

res.lopez <- c(
  dm.test(losses.HAR.REQ, losses.GEV, alternative = "less", power = 1)$p.value,
  dm.test(losses.HAR.REQ, losses.QR, alternative = "less", power = 1)$p.value,
  dm.test(losses.HAR.REQ, losses.REQ, alternative = "less", power = 1)$p.value,
  dm.test(losses.HAR.REQ, losses.HAR.QR, alternative = "less", power = 1)$p.value
)

res <- cbind(res.lopez, res.caporin)
res

rownames(res) <- c("GEV", "QR", "REQ", "HAR.QR")

res

library(xtable)

xtable(res)

### MCS procedure

losses.GEV <- tab.res$GEV
losses.GEV[tab.res$r < tab.res$GEV] <- (1 + (tab.res$r - tab.res$GEV)^2)[tab.res$r < tab.res$GEV]
losses.GEV[tab.res$r >= tab.res$GEV] <- ( (tab.res$r < 0) * (abs(tab.res$GEV) - abs(tab.res$r)) )[tab.res$r >= tab.res$GEV]

losses.QR <- tab.res$QR
losses.QR[tab.res$r < tab.res$QR] <- (1 + (tab.res$r - tab.res$QR)^2)[tab.res$r < tab.res$QR]
losses.QR[tab.res$r >= tab.res$QR] <- ( (tab.res$r < 0) * (abs(tab.res$QR) - abs(tab.res$r)) )[tab.res$r >= tab.res$QR]

losses.REQ <- tab.res$REQ
losses.REQ[tab.res$r < tab.res$REQ] <- (1 + (tab.res$r - tab.res$REQ)^2)[tab.res$r < tab.res$REQ]
losses.REQ[tab.res$r >= tab.res$REQ] <- ( (tab.res$r < 0) * (abs(tab.res$REQ) - abs(tab.res$r)) )[tab.res$r >= tab.res$REQ]

losses.HAR.QR <- tab.res$HAR.QR
losses.HAR.QR[tab.res$r < tab.res$HAR.QR] <- (1 + (tab.res$r - tab.res$HAR.QR)^2)[tab.res$r < tab.res$HAR.QR]
losses.HAR.QR[tab.res$r >= tab.res$HAR.QR] <- ( (tab.res$r < 0) * (abs(tab.res$HAR.QR) - abs(tab.res$r)) )[tab.res$r >= tab.res$HAR.QR]

losses.HAR.REQ <- tab.res$HAR.REQ
losses.HAR.REQ[tab.res$r < tab.res$HAR.REQ] <- (1 + (tab.res$r - tab.res$HAR.REQ)^2)[tab.res$r < tab.res$HAR.REQ]
losses.HAR.REQ[tab.res$r >= tab.res$HAR.REQ] <- ( (tab.res$r < 0) * (abs(tab.res$HAR.REQ) - abs(tab.res$r)) )[tab.res$r >= tab.res$HAR.REQ]

library(MCS)
matr.losses <- cbind(losses.GEV, losses.QR, losses.REQ, losses.HAR.QR, losses.HAR.REQ)
colnames(matr.losses) <- c("GEV", "QR", "REQ", "HAR.QR", "HAR.REQ")
MCSprocedure(Loss = matr.losses, alpha = 0.3, statistic = "TR")

# Solo Lopez

losses.GEV <- (tab.res$r < tab.res$GEV) * (1 + (tab.res$r - tab.res$GEV)^2)
losses.QR <- (tab.res$r < tab.res$QR) * (1 + (tab.res$r - tab.res$QR)^2)
losses.REQ <- (tab.res$r < tab.res$REQ) * (1 + (tab.res$r - tab.res$REQ)^2)
losses.HAR.QR <- (tab.res$r < tab.res$HAR.QR) * (1 + (tab.res$r - tab.res$HAR.QR)^2)
losses.HAR.REQ <- (tab.res$r < tab.res$HAR.REQ) * (1 + (tab.res$r - tab.res$HAR.REQ)^2)

matr.losses <- cbind(losses.GEV, losses.QR, losses.REQ, losses.HAR.QR, losses.HAR.REQ)
colnames(matr.losses) <- c("GEV", "QR", "REQ", "HAR.QR", "HAR.REQ")
MCSprocedure(Loss = matr.losses, alpha = 0.3, statistic = "TR")

colnames(mr1m)[1]
colnames(mr1m)[ncol(mr1m)]







