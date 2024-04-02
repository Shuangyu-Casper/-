# 階段沖堤的預測模型
PRED <- function(x){
  VR <- c(0, 0, 0)
  VRadj <- c(rep(0,length(VE)))
  P <- c(rep(1,length(VE)))
  ER <- c(rep(0,length(VE)))
  state <- c(rep("已洗脫",length(VE)))
  RM <- log10(1/x-1)
  
  for (i in 1:length(VE)){
    VRadj[i] <- P[i] * Ktr * (1/(x[i]*RFadjF)-1) * VM
    
    if(VRadj[i]!=0){
      VR[i] <- VM + VRadj[i]}else{
        VR[i] <-0
      }
    
    if(VE[i]/VR[i] < 1){
      ER[i] <- VE[i]/VR[i]
      if(i+1<=length(VE)) P[i+1] <- P[i]-ER[i]*P[i]
    }else{
      ER[i] <- 1
      if(i+1<=length(VE)) P[i+1] <- 0
    }
    
    if(ER[i] == 1){
      state[i]<-"已沖提出"
    }else{
      state[i]<-"未沖提出"
    }
  }
  data.frame(x, x*RFadjF, RM, VR, VRadj, state)
}

VEPred <- function(x){
  e=0
  VEPred <- 0 
  for(i in 1:length(VE)){
    if(x$state[i]=="已沖提出"){
      e <- i
      break
    }
  }
  if(e>1){
    for(i in 1:(e-1)){
      if(is.finite(VE[i])){
        VEPred <- VEPred + VE[i]}
    }
  }
  VEPred <- VEPred + x$VR[e]
  return(VEPred)
}

Vb <- function(x, N){
  4*x/sqrt(N)
}

# 請使用者輸入檔案路徑
filepath <- readline("請輸入檔案：")
DATA <- read.csv(filepath)


# 設定死體積和理論板數
VM <- DATA[1,6]
N <- DATA[1,9]

# 指定校正RF值時使用的倍率和Ktr值
Ktr <- DATA[1,7]
RFadjF <- DATA[1,8]

# 輸入各沖堤體積和對應的RF值
VE <- DATA[,3]
RFA <- DATA[,4]
RFB <- DATA[,5]

# 預測A化合物滯留體積
compoundA <- PRED(RFA)
VEPredA <- VEPred(compoundA)
VbA <- Vb(VEPredA, N)

# 預測B化合物滯留體積
compoundB <- PRED(RFB)
VEPredB <- VEPred(compoundB)
VbB <- Vb(VEPredB, N)


result <- data.frame(VE, compoundA, compoundB)
colnames(result) <- c("VE", "RF(A)", "RF'(A)","RM(A)", "VR(A)", "VR'(A)", "狀態(A)",
                      "RF(B)", "RF'(B)","RM(B)", "VR(B)", "VR'(B)", "狀態(B)")


Rs <- abs(2*(VEPredB-VEPredA)/(VbA+VbB))

summary <- rbind(VEPredA, VbA, VEPredB, VbB, Rs)
sumaary <- data.frame(summary)


# 繪製示意圖
# 請使用者輸入圖檔儲存路徑
filepath <- readline("請輸入儲存結果的檔案路徑：")
setwd(filepath)
png("chromotagram.png", width = 1200, height = 900, res=300)
par(mar=c(4, 0.5, 0.5, 0.5))
x <- seq(0, VEPredB+VbB+10, by = 0.333)
yA <- dnorm(x, mean = VEPredA, sd = VbA/4)
yB <- dnorm(x, mean = VEPredB, sd = VbB/4)
plot(x, yA, type="l", xlab="沖堤體積 (mL)", ylab="", yaxt="n")
par(new=TRUE)
plot(x, yB, type="l", xlab="沖堤體積 (mL)", ylab="", yaxt="n")
dev.off()

# 匯出CSV檔案
write.csv(result, file="result.csv", row.names = T,  fileEncoding = "Big5")
write.csv(summary, file="summary.csv", row.names = T,  fileEncoding = "Big5")