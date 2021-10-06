library(glmnet)
data.Y<- read.table('C:/Users/xuewu/Downloads/OneDrive_1_10-6-2021/GEUVADIS_normalized_expression_chr20', header=TRUE)
data.X<- read.table('C:/Users/xuewu/Downloads/OneDrive_1_10-6-2021/GEUVADIS_chr20_processed.traw', header=TRUE)
data.Y<- as.data.frame(data.Y)
data.X<- as.data.frame(data.X)
####### data cleaning based on genetic variants +/- 500,000 basepairs ########
gene_idx<- 1 
lower_pos<- data.Y$start[gene_idx] - 500000
higher_pos<- data.Y$end[gene_idx] + 500000
data.X_sub<- data.X[data.X$POS <= higher_pos&data.X$POS>=lower_pos,]
dim(data.X_sub)
data.X_sub<- data.X_sub[,-c(1, 3:6)]
dim(data.X_sub)
data.Y_sub<- data.Y[gene_idx,-(2:4)]

rownames(data.Y_sub)<- data.Y_sub$gene_id
y<- apply(data.Y_sub[,-1], 2, as.numeric)
y_scale<- scale(y)

x<- data.X_sub
x<- t(x)
colnames(x)<- x[1,]
x<- x[-1,]
x<- apply(x, 2, as.numeric)
x_scale<- apply(x, 2, scale)

N<- length(y)



soft_threshold<- function(rho, lambda){
  return(ifelse(rho>lambda, rho-lambda, ifelse(rho< -lambda, rho+lambda, 0)))
}

lambda<- 0.1
counter<- 0 
t<- proc.time()
b<- 1
##### Firstly, checking the scenario without intercept
intercept=FALSE
if(intercept==TRUE){
  x_scale<- cbind(1, x_scale)
  
  #x_scale.sub<- cbind(1, x_scale.sub)
  p<- ncol(x_scale_int)
}else{p<- ncol(x_scale)}
#x_scale.sub<- cbind(1, x_scale.sub)
p<- ncol(x_scale)
beta_init<- rep(1, p)

repeat{
  beta<- beta_init
  for (i in 1:p){
    xy<- t(y_scale) %*% x_scale[,i] 
    #xb<- x_scale[,-i] %*% beta[-i]
    xxb<- t(x_scale[,i]) %*% (x_scale[,-i] %*% beta[-i])
    #  xx<- x_scale[,i] %*% x_scale[,i]
    #ope_pos<- ((xy-xxb)-N*lambda)/xx
    #ope_neg<- ((xy-xxb)+N*lambda)/xx
    rho<- (xy-xxb)/N
    if(intercept==TRUE){
      if(i==1){
        beta[i]<- rho
      }else{beta[i]<- soft_threshold(rho, lambda)}
    }
    #ope<- ifelse((xy-xxb)>(N*lambda), ((xy-xxb)-N*lambda)/xx,
    #                     ifelse((xy-xxb)< -(N*lambda), ((xy-xxb)+N*lambda)/xx, 0))
    # beta[i]<- beta[i]-alpha*ope
    if(intercept==FALSE){
      beta[i]<- soft_threshold(rho, lambda)
    }
  }
  b<-b+1
  #if(max(abs(beta-beta_init))<0.1){
  # break
  #}
  if(b>1000){break}
  beta_init<- beta
  counter <- counter + 1 ####### FIXME
  cat("Iteration:",counter,"time elapsed:",(proc.time()-t)[3],"\n") ####### FIXME
}
############# Without intercept, we selected 14 variables with lambda=0.1 and the corresponding estimated coefficients as following. 
colnames(x_scale)[which(beta!=0)]
#[1] "20_49064671_G_A_b37"     "20_49078164_C_CAG_b37"   "20_49101662_A_G_b37"     "20_49374183_G_GA_b37"    "20_49424564_C_T_b37"    
#[6] "20_49647371_C_T_b37"     "20_49693755_T_C_b37"     "20_49710044_G_A_b37"     "20_49714264_G_A_b37"     "20_49715109_C_T_b37"    
#[11] "20_49936124_G_A_b37"     "20_49938899_G_A_b37"     "20_49974518_G_GTGGA_b37" "20_50064221_G_A_b37"
beta[beta!=0]
#[1]  0.0312533426  0.0056028048  0.0027816559 -0.0167844356 -0.0065485858 -0.0009462111  0.0132520156 -0.0253218738 -0.0120617099
#[10] -0.0095599004 -0.0082330625  0.0431645887 -0.0476844724  0.0244480452

################# Comparing the my results with glmnet results##################
g<- glmnet(x_scale, y_scale, intercept=FALSE, standardize=FALSE, lambda=0.1)
coef(g)@Dimnames[[1]][which(coef(g)!=0)]
# [1] "20_49064671_G_A_b37"     "20_49078164_C_CAG_b37"   "20_49101662_A_G_b37"     "20_49374183_G_GA_b37"    "20_49424564_C_T_b37"    
#[6] "20_49647371_C_T_b37"     "20_49693755_T_C_b37"     "20_49710044_G_A_b37"     "20_49714264_G_A_b37"     "20_49715109_C_T_b37"    
#[11] "20_49936124_G_A_b37"     "20_49938899_G_A_b37"     "20_49974518_G_GTGGA_b37" "20_50064221_G_A_b37" 
coef(g)[which(coef(g)!=0)]
#[1]  0.0313443015  0.0056045137  0.0027784027 -0.0168020935 -0.0065517145 -0.0009215823  0.0132782017 -0.0254068816 -0.0120501940
#[10] -0.0096074391 -0.0082409674  0.0433045782 -0.0478473965  0.0245129233

#### As can be seen, the selected SNPs from my codes are identical to the results from glmnet function, and the estimated coefficients
#### are approximately equal.


##### Secondly, checking the scenario with intercept
intercept=TRUE
if(intercept==TRUE){
  x_scale<- cbind(1, x_scale)
  
  #x_scale.sub<- cbind(1, x_scale.sub)
  p<- ncol(x_scale_int)
}else{p<- ncol(x_scale)}
#x_scale.sub<- cbind(1, x_scale.sub)
p<- ncol(x_scale)
beta_init<- rep(1, p)

repeat{
  beta<- beta_init
  for (i in 1:p){
    xy<- t(y_scale) %*% x_scale[,i] 
    #xb<- x_scale[,-i] %*% beta[-i]
    xxb<- t(x_scale[,i]) %*% (x_scale[,-i] %*% beta[-i])
    #  xx<- x_scale[,i] %*% x_scale[,i]
    #ope_pos<- ((xy-xxb)-N*lambda)/xx
    #ope_neg<- ((xy-xxb)+N*lambda)/xx
    rho<- (xy-xxb)/N
    if(intercept==TRUE){
      if(i==1){
        beta[i]<- rho
      }else{beta[i]<- soft_threshold(rho, lambda)}
    }
    #ope<- ifelse((xy-xxb)>(N*lambda), ((xy-xxb)-N*lambda)/xx,
    #                     ifelse((xy-xxb)< -(N*lambda), ((xy-xxb)+N*lambda)/xx, 0))
    # beta[i]<- beta[i]-alpha*ope
    if(intercept==FALSE){
      beta[i]<- soft_threshold(rho, lambda)
    }
  }
  b<-b+1
  #if(max(abs(beta-beta_init))<0.1){
  # break
  #}
  if(b>1000){break}
  beta_init<- beta
  counter <- counter + 1 ####### FIXME
  cat("Iteration:",counter,"time elapsed:",(proc.time()-t)[3],"\n") ####### FIXME
}
############# Without intercept, we selected 15 variables with lambda=0.1 and the corresponding estimated coefficients as following. 
colnames(x_scale)[which(beta!=0)]
#[1] ""                        "20_49064671_G_A_b37"     "20_49078164_C_CAG_b37"  
#[4] "20_49101662_A_G_b37"     "20_49374183_G_GA_b37"    "20_49424564_C_T_b37"    
#[7] "20_49647371_C_T_b37"     "20_49693755_T_C_b37"     "20_49710044_G_A_b37"    
#[10] "20_49714264_G_A_b37"     "20_49715109_C_T_b37"     "20_49936124_G_A_b37"    
#[13] "20_49938899_G_A_b37"     "20_49974518_G_GTGGA_b37" "20_50064221_G_A_b37"    
beta[which(beta!=0)]
#[1]  5.967158e-16  3.125334e-02  5.602805e-03  2.781656e-03 -1.678444e-02 -6.548586e-03 -9.462111e-04  1.325202e-02 -2.532187e-02
#[10] -1.206171e-02 -9.559900e-03 -8.233062e-03  4.316459e-02 -4.768447e-02  2.444805e-02

################# Comparing the my results with glmnet results##################
g<- glmnet(x_scale, y_scale, intercept=TRUE, standardize=FALSE, lambda=0.1)
coef(g)@Dimnames[[1]][which(coef(g)!=0)]
#[1] "(Intercept)"             "20_49064671_G_A_b37"     "20_49078164_C_CAG_b37"   "20_49101662_A_G_b37"     "20_49374183_G_GA_b37"   
#[6] "20_49424564_C_T_b37"     "20_49647371_C_T_b37"     "20_49693755_T_C_b37"     "20_49710044_G_A_b37"     "20_49714264_G_A_b37"    
#[11] "20_49715109_C_T_b37"     "20_49936124_G_A_b37"     "20_49938899_G_A_b37"     "20_49974518_G_GTGGA_b37" "20_50064221_G_A_b37"
coef(g)[which(coef(g)!=0)]
#[1]  5.715633e-16  3.134430e-02  5.604514e-03  2.778403e-03 -1.680209e-02 -6.551715e-03 -9.215823e-04  1.327820e-02 -2.540688e-02
#[10] -1.205019e-02 -9.607439e-03 -8.240967e-03  4.330458e-02 -4.784740e-02  2.451292e-02

#### As can be seen, the selected SNPs from my codes are identical to the results from glmnet function, and the estimated coefficients
#### are approximately equal.