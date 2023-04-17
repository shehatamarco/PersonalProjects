
rm(list=ls())

load("/Users/marcoshehata/Desktop/UniPD/Magistrale/AnnoII/AnalisiFinanza/Homework/Quantile_Regression/mr1m.RDATA")

RV_funct <- function(matr) as.vector(apply(matr^2, 2, sum)) # RV
RV <- RV_funct(mr1m)
r <- as.vector(apply(mr1m, 2, sum)) # rendimenti
l <- -r # perdite
alpha <- 0.01 # ordine del VAR 
length(r) # 878
n.info <- 300 # numero di unità usate per stima

# EVT: 

# modello completo

library(evir)

tot.mod <- gev(l, block = 22)
par(mfrow=c(1,1))
plot(tot.mod, main = "QQPLOT GEV")

# BackCasting

VaR.GEV <- function(mod, alpha = 0.01){
  
  mu <- mod$par.ests[3]
  sigma <- mod$par.ests[2]
  xi <- mod$par.ests[1]
  m <- mod$block
  
  res <- mu - sigma/xi * (1 - (-m*log(1-alpha))^(-xi) )
  
  return(-res)
  
} 

# calcolo VAR

VaR.GEV.s <- NULL

for(i in 1:(length(r) - n.info)){
  
  l.used <- l[i:(n.info+i-1)]
  
  EVT.mod <- gev(l.used, block = 22)
  
  VaR.GEV.s <- c(VaR.GEV.s, VaR.GEV(EVT.mod, alpha = alpha))
  
  cat("Done:", i, "\n")
  
}

sum(is.na(VaR.GEV.s))

length(VaR.GEV.s) # OK

# Calcolo Regressione Quantile

library(quantreg)

# stima Totale

r <- r[2:878]
RV <- RV[1:877]

QR.mod <- rq(
  r ~ RV, 
  tau = alpha
)

summary(QR.mod)

# BackCasting

tmp <- data.frame(
  r = r[2:length(r)],
  RV = RV[1:(length(r) - 1)]
) 

VaR.QR.s <- NULL

for(i in 1:(nrow(tmp) - n.info)){
  
  tmp.used <- tmp[i:(n.info+i-1),]
  
  QR.mod <- rq(
    r ~ RV, 
    data = tmp.used,
    tau = alpha
  )
  
  VaR.est <- predict(QR.mod, newdata = tmp[(n.info+i),])
  
  VaR.QR.s <- c(VaR.QR.s, VaR.est)
  
  cat("Done:", i, "\n")
  
}

VaR.QR.s <- c(NA, VaR.QR.s) # il primo non si può stimare con 300 unità

# vediamo come esce fin ora

plot(r[(n.info+2):(length(r))], type = "l", ylim = c(-0.2, max(r)+0.01), xlab = "t", ylab = "rendimenti")
lines(VaR.GEV.s[-1], type = "l", col = 3) # GEV in verde
lines(VaR.QR.s[-1], type = "l", col = 2) # QR in rosso

# eccedenze
length((n.info+2):(length(r)))*alpha # 5.77, ovvero 6 circa

sum(r[(n.info+2):(length(r))] < VaR.GEV.s[-1]) # entrambi sforano di 9 
sum(r[(n.info+2):(length(r))] < VaR.QR.s[-1])
# Vediamo che succede con REQ

# REQ: xi stimato via HILL

VaR.REQ.H.s <- NULL

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

for(i in 1:(nrow(tmp) - n.info)){
  
  tmp.used <- tmp[i:(n.info+i-1),]
  
  QR.mod <- rq(
    r ~ RV, 
    data = tmp.used,
    tau = alpha
  )
  
  q.t <- predict(QR.mod, newdata = tmp.used)
  
  q.new <- predict(QR.mod, newdata = tmp[(n.info+i),]) # previsione quantile 1 passo in avanti
  
  z.t <- tmp.used$r / q.t # calcolo i quantili depurati
  
  VaR.est <- VaR.REQ.H(z.t, q.new, a = alpha) # stimo il var
  
  VaR.REQ.H.s <- c(VaR.REQ.H.s, as.numeric(VaR.est)) # lo salvo
  
  cat("Done:", i, "\n")
  
} # ci mette molto di più ovviamente, ma non ci preoccupa questo

VaR.REQ.H.s <- c(NA, VaR.REQ.H.s)

# Numero di sforamenti?

sum(r[(n.info+2):(length(r))] < VaR.REQ.H.s[-1]) # arriva a 6

# plot delle varie cosine:

plot(r[(n.info+2):(length(r))], type = "l", ylim = c(-0.2, max(r)+0.01), xlab = "t", ylab = "rendimenti",
     xlim = c(0, 520))
lines(VaR.GEV.s[-1], type = "l", col = 3) # GEV in verde
lines(VaR.QR.s[-1], type = "l", col = 2) # QR in rosso
lines(VaR.REQ.H.s[-1], type = "l", col = 4) # QR in rosso

# coincidono posizioni eccedenze?

tab.res <- data.frame(
  data = colnames(mr1m)[(n.info+2):(length(r))],
  r = r[(n.info+2):(length(r))],
  GEV = VaR.GEV.s[-1],
  QR = VaR.QR.s[-1],
  REQ = VaR.REQ.H.s[-1]
)

head(tab.res)

tab.res$data[tab.res$r < tab.res$GEV]
tab.res$data[tab.res$r < tab.res$QR]
tab.res$data[tab.res$r < tab.res$REQ]

# TEST DM:

tab.res$data <- NULL

losses.GEV <- (tab.res$r < tab.res$GEV) * (1 + (tab.res$r - tab.res$GEV)^2)
losses.QR <- (tab.res$r < tab.res$QR) * (1 + (tab.res$r - tab.res$QR)^2)
losses.REQ <- (tab.res$r < tab.res$REQ) * (1 + (tab.res$r - tab.res$REQ)^2)

library(forecast)
dm.test(losses.REQ, losses.QR, alternative = "greater", power = 1) 
dm.test(losses.REQ, losses.GEV, alternative = "greater", power = 1) 

# REQ su tutti i dati, definisco GEV  e vedo come si adatta sugli Z_t

QR.mod <- rq(
  r[2:(length(r))] ~ RV[1:(length(r)-1)], 
  tau = alpha
)

z.t <- r[2:(length(r))] / fitted(QR.mod)

gev.mod <- gev(z.t, 22)

plot(gev.mod, main = "QQPLOT GEV su Z")

save(tab.res, file = "/Users/marcoshehata/Desktop/UniPD/Magistrale/AnnoII/AnalisiFinanza/Homework/Quantile_Regression/VaRQRGEVREQ.RDATA")





