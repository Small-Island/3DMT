#######　分散分析ｓＡ(１要因参加者内)　#######

## js-STARの入力: sA
levA <- 3
n    <- c(15,15,15)

data <- c(
 1, 2, 1,
 1, 2, 1,
 1, 2, 1,
 2, 3, 4,
 3, 2, 2,
 1, 1, 2,
 1, 2, 1,
 2, 3, 2,
 2, 3, 2,
 1, 3, 1,
 2, 3, 2,
 1, 2, 1,
 3, 3, 2,
 2, 2, 2,
 1, 2, 3
)
## js-STARの入力おわり

## スタック形式のデータセット
s  <- gl(mean(n), levA)
A  <- gl(levA, 1, n[1]*levA)
sA <- data.frame(s, A, data)
N  <- mean(n)

## 基本統計量の計算
hk  <- tapply(sA$data, sA$A, mean)
sd  <- tapply(sA$data, sA$A, sd)
mi  <- tapply(sA$data, sA$A, min)
ma  <- tapply(sA$data, sA$A, max)

tx1 <- matrix(round(c(
  n,
  hk,
  sd,
  mi,
  ma
  ), 4), nc=5 )
colnames(tx1) <- c(
  "   Ｎ",
  "     Mean",
  "       SD",
  "   Min",
  "   Max" )
rownames(tx1) <- c(
  paste("A.", 1:levA, sep="")
  )


## 作図：プロット＆ラインズ
suij <- 1:levA           # 水準番号
xjik <- c(1-.3, levA+.3) # ｘ軸の拡幅
bai=2                    # ｙ軸の拡幅倍率
yjik <- c(
  min(hk)-max(sd)*bai,
  max(hk)+max(sd)*bai)

plot(suij, hk,                 # 描点指定
  bty="l", xli=xjik, yli=yjik, # 枠の設定
  tck=0.03,                    # 目盛り突起
  las=1,                       # 目盛り正立
  xlab=c("水準番号"),
  ylab=c(""),
  xaxt="n",
  main=c("標準偏差= 不偏分散の平方根")
  )
axis(side=1,
     at=1:levA,
     tck=0.03,
     labels=c("A1","A2","A3")
     )
arrows(
  suij,hk, suij,hk+sd, # ＳＤ上向き
  ang=90, le=0.15 )    # 矢羽の角度と長さ
arrows(
  suij,hk, suij,hk-sd, # ＳＤ下向き
  ang=90, le=0.15 )

lines (suij,hk, lwd=2,col=1) # 線幅、線色
points(suij,hk, pch=21,      # 描点21番
       cex=3, bg=1 )         # ４倍、黒色


## 分散分析
kx2 <- c()
kx2 <- anova(lm(data~s+A, sA))
tx2 <- c()
tx2 <- matrix(round( c(
  kx2$S,kx2$D,kx2$M,kx2$F,kx2$P
  ), 4), nr=3)
colnames(tx2) <- c(
  "       SS",
  "  df",
  "       MS",
  "      Ｆ",
  "     ｐ")
rownames(tx2) <- c(
  "　ｓ　",
  "要因Ａ",
  "ｓ×Ａ")

kx3 <- c()
kx3 <- pairwise.t.test( # 多重比較
  sA$data, sA$A,
  p.ad="BH", pair=1)    # 対応あり
tx3 <- c()
tx3 <- round(kx3$p.v, 4)
dimnames(tx3) <- list(
  c(paste("A.", 2:levA,   sep="")),
  c(paste("　　 A.", 2:levA-1, sep="")) )

## パワーアナリシス
df1 <- kx2$D[2]
df2 <- kx2$D[3]
Fa  <- kx2$F[2]
es  <- sqrt(Fa*df1/df2) # ES.f


# 検出力
epsi <- 1  # 球面性ε=1

dmat <- matrix(sA$data,nr=N,by=1)
fugo <- sign(cor(dmat))
vsum <- sum(cor(dmat)^2*fugo)-levA
corr <- sqrt(  # 平均相関
  abs(vsum)/(levA^2-levA)
  )*sign(vsum)

corr <- abs(sum(cor(dmat))-levA)/(levA^2-levA) # 単純平均

NCP <-  # 非心度
  es^2*levA*N*(1/(1-corr))*epsi
crF <- qf(0.05,
  df1=df1,df2=df2,low=0)
pow <- 1-pf(crF,
  df1=df1,df2=df2,ncp=NCP)

ncp0 <- es^2*levA*N
pow0 <- 1-pf(
  crF,df1=df1,df2=df2,ncp=ncp0)

ncpSPSS <- es^2*df2
powSPSS <- 1-pf(
  crF,df1=df1,df2=df2,ncp=ncpSPSS)

# α
crF  <- qf(0.20,df1=df1,df2=df2,ncp=NCP)
alph <- pf(crF,df1=df1,df2=df2,low=0)

# Ｎのシミュレーション
DFv <- NA
kos <- 0   # 出力回数
if (pow<0.80)
  for(i in 1:100000){
    if (kos>1) i=99999
    else {
      DF2  <- df2+(i*0.01)
      DF2  <- min(DF2, 1000)
      ncp <- es^2*levA*(DF2/df1+1)*(1/(1-corr))*epsi
      crF  <- qf(0.95,df1=df1,df2=DF2)
      beta <- pf(crF, df1=df1,df2=DF2,ncp=ncp)

      if (beta<0.20) kos=kos+1
      if (kos==1)    DFv <- c(DF2, DFv) }
  }
if (pow>0.80)
  for(i in 1:100000){
    if (kos>1) i=99999
    else {
      DF2  <- df2-(i*0.01)
      DF2  <- max(DF2, 2)
      ncp <- es^2*levA*(DF2/df1+1)*(1/(1-corr))*epsi
      crF  <- qf(0.95,df1=df1,df2=DF2)
      beta <- pf(crF, df1=df1,df2=DF2,ncp=ncp)

      if (beta>0.20) kos=kos+1
      if (kos==1)    DFv <- c(DF2, DFv) }
  }

tx4 <- c()
tx4 <- matrix(round(c(
  Fa, es, pow, pow0, powSPSS
  ), 4), nr=1 )
colnames(tx4) <- c(
  "Ｆ値",
  " 効果量ｆ",
  " 検出力1",
  " 検出力0",
  " 検出力2" )
rownames(tx4) <- c(
  "要因Ａ" )

tz8 <- c()
tz8 <- matrix(round(c(
  es, 0.80, alph, df1, df2,
  es, 0.80, 0.05, df1, DFv[1]
  ), 4), nr=2, by=1 )
colnames(tz8) <- c(
  "効果量ｆ",
  " 検出力",
  "     α",
  "df1",
  "    df2")
rownames(tz8) <- c(
  " α の計算",
  " Ｎ の計算")

tz9 <- c()
tz9 <- matrix(round(c(
   ceiling(DFv[1]/df1+1),
   N,
   corr
   ), 4), nr=1)
rownames(tz9) <- c("Ｎ=df2/df1+1")
colnames(tz9) <- c(
  "次回総数",
  "  今回のＮ",
  "  水準間相関" )

# Mauchly.test
kx5 <- c()
kx5 <- mauchly.test(lm(dmat~1), X=~1)
tx5 <- c()
tx5 <- matrix(round(c(
  kx5$stat,
  kx5$p.val,
  1/(levA-1)
  ), 4), nr=1 )
rownames(tx5) <- c("要因Ａ")
colnames(tx5) <- c(
  "Mauchly's W",
  "    ｐ値",
  "  下限ε" )

kx6 <- c()
kx6 <- anova(lm(dmat~1), X=~1,
  test="Spherical" )
tx6 <- c()
tx6 <- matrix(round(c(
  kx6$F[1],
  kx6$Pr[1],
  kx6$G[1],
  kx6$H[1]
  ), 4), nr=1)
colnames(tx6) <- c(
  "   Ｆ値",
  "   ｐ値",
  "  G-G調整値",
  " H-F調整値" )
rownames(tx6) <- c(
  "要因Ａ" )

options(digits=5)
options(scipen=5)
###########################
#  分散分析 sＡ-design：  #
#　　１要因 参加者内計画  #
###########################
tx1 # 基本統計量（SD=不偏分散の平方根）

tx2 # 分散分析表

tx4 # 効果量ｆと検出力(1-β)
    # Part.η2(偏イータ２乗)= f^2/(f^2+1)
    # 検出力1 は非心度推定= f^2*DataSize, r=r
    # 検出力0 は非心度推定= f^2*DataSize, r=0
    # 検出力2 は非心度推定= f^2*DFsxa

tx3 # 多重比較（参加者内ｔ検定）
# 数値は調整後のｐ値（両側確率）
# ｐ値の調整は Benjamini & Hochberg(1995) による

tz8 # パワーアナリシス
    # ■ NA が出力されたら計算不能
tz9 # Ｎ(total sample size)の計算

tx5 # 球面性の検定（水準=２なら不要）
    # p<α なら球面性不成立！
tx6 # 球面性不成立のときの自由度調整Ｆ検定

# ○オプション：［↑］⇒行頭の♯を消す⇒［Enter］
# windows();boxplot(data~A,d=sA,las=1) # 箱ひげ図
# kx6[0,] # 自由度調整係数ε(epsilon)

# _/_/_/ Powered by js-STAR _/_/_/
