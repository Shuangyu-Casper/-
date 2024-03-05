# 文獻中的死體積預測公式
VMpred <- function (x){
  1.8*x+0.3
}

# 階段沖堤的預測模型
PRED <- function(x){
  VR <- c(0, 0, 0)
  VRadj <- c(rep(0,length(VE)))
  P <- c(rep(1,length(VE)))
  ER <- c(rep(0,length(VE)))
  state <- c(rep("已洗脫",length(VE)))
  
  for (i in 1:length(VE)){
    VRadj[i] <- P[i] * Ktr * (1/x[i]-1) * VM
    
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
      state[i]<-"已洗脫"
    }else{
      state[i]<-"未洗脫"
    }
  }
  data.frame(x, x*RFadjF, VR, VRadj, state)
}

VRTotal <- function(x){
  e=0
  VRtotal <- 0 
  for(i in 1:length(VE)){
    if(x$state[i]=="已洗脫"){
      e <- i
      break
    }
  }
  if(e>1){
    for(i in 1:(e-1)){
      VRtotal <- VRtotal + VE[i]
    }
  }
  VRtotal <- VRtotal + x$VR[e]
  return(VRtotal)
}

Vb <- function(x, N){
  4*x/sqrt(N)
}

# 請使用者輸入檔案路徑
filepath <- readline("請輸入檔案：")
DATA <- read.csv(filepath)


# 設定矽膠重、死體積和理論板數
W <- 40
VM <- VMpred(10)
N <- 90

# 指定校正RF值時使用的倍率和Ktr值
Ktr <- 1
RFadjF <- 1.5

# 輸入各沖堤體積和對應的RF值
VE <- c(120, 120, 120)
RFA <- c(0.2, 0.3, 0.4)
RFB <- c(0.1, 0.1, 0.2)

# 計算校正後RF值
RFAadj <- RFA*RFadjF
RFBadj <- RFB*RFadjF

# 預測A化合物滯留體積
compoundA <- PRED(RFA)
VRTotalA <- VRTotal(compoundA)
VbA <- Vb(VRTotalA, N)

# 預測B化合物滯留體積
compoundB <- PRED(RFB)
VRTotalB <- VRTotal(compoundB)
VbB <- Vb(VRTotalB, N)


result <- data.frame(VE, compoundA, compoundB)
colnames(result) <- c("VE", "RF(A)", "RF'(A)", "VR(A)", "VR'(A)", "狀態(A)",
                      "RF(B)", "RF'(B)", "VR(B)", "VR'(B)", "狀態(B)")


Rs <- abs(2*(VRTotalB-VRTotalA)/(VbA+VbB))

summary <- rbind(VR(A), Vb(A), VR(B), Vb(B), Rs)
sumaary <- data.frame(summary)

# 繪製示意圖
# 請使用者輸入圖檔儲存路徑
filepath <- readline("請輸入儲存結果的檔案路徑：")
setwd(filepath)
png("chromotagram.png", width = 1200, height = 900, res=300)
par(mar=c(4, 0.5, 0.5, 0.5))
x <- seq(0, VRTotalB+VbB+10, by = 0.333)
yA <- dnorm(x, mean = VRTotalA, sd = VbA/4)
yB <- dnorm(x, mean = VRTotalB, sd = VbB/4)
plot(x, yA, type="l", xlab="沖堤體積 (mL)", ylab="", yaxt="n")
par(new=TRUE)
plot(x, yB, type="l", xlab="沖堤體積 (mL)", ylab="", yaxt="n")
dev.off()

# 匯出CSV檔案
write.csv(result, file="result.csv",row.names = T,  fileEncoding = "Big5")
write.csv(summary, file="summary.csv",row.names = T,  fileEncoding = "UTF-8")
