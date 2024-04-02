PRED <- function(x){
  k <- Ktr * (1/x-1)
  RM <- log10(1/x-1)
  VRadj <- Ktr * (1/(x*RFadjF)-1) * VM
  VR　<- VRadj + VM
  Vb <- 4*VR/sqrt(N)
  ER <- c(rep(0,2))
  state <- c(rep("已洗脫",2))
  
  for(i in 1:2){
    if(VE/VR[i] < 1){
      ER[i] <- VE/VR[i]
    }else{
      ER[i] <- 1
    }
    
    if(ER[i] == 1){
      state[i]<-"已沖提出"
    }else{
      state[i]<-"未沖提出"
    }
  }
  
  data.frame(x, x*RFadjF, RM, k, VR, VRadj, Vb, state)
}


# 請使用者輸入檔案路徑
filepath <- readline("請輸入檔案：")
DATA <- read.csv(filepath)


# 設定沖提體積、死體積和理論板數
VE <- DATA[1,2]
VM <- DATA[1,5]
N <- DATA[1,8]

# 指定校正RF值時使用的倍率和Ktr值
Ktr <- DATA[1,6]
RFadjF <- DATA[1,7]

# 輸入各化合物和對應的RF值
RF <- c(DATA[1,3], DATA[1,4])

# 預測化合物滯留因子和體積等參數
compound <- PRED(RF)


Rs <- abs(2*(compound[2,4]-compound[1,4])/(compound[2,6]+compound[1,6]))
result <- data.frame(c(VE, VE), compound, c(Rs,Rs))
colnames(result) <- c("VE", "RF", "RF'", "RM", "k", "VR", "VR adj", "Vb", "狀態","Rs")


# 繪製示意圖
# 請使用者輸入圖檔儲存路徑
filepath <- readline("請輸入儲存結果的檔案路徑：")
setwd(filepath)
png("chromotagram.png", width = 1200, height = 900, res=300)
par(mar=c(4, 0.5, 0.5, 0.5))
x <- seq(0, compound[2,4]+compound[2,6]+10, by = 0.333)
yA <- dnorm(x, mean = compound[1,4], sd = VbA/4)
yB <- dnorm(x, mean = compound[2,4], sd = VbB/4)
plot(x, yA, type="l", xlab="沖堤體積 (mL)", ylab="", yaxt="n")
par(new=TRUE)
plot(x, yB, type="l", xlab="沖堤體積 (mL)", ylab="", yaxt="n")
dev.off()

# 匯出CSV檔案
write.csv(result, file="result.csv", row.names = T,  fileEncoding = "Big5")
