
library(lpSolve)
library(ggplot2)
source("disparate.R")
#####################################################
# Obtaining Fairness using Optimal Transport Theory #
#####################################################
#
#Simulations: Geometric and Random repair
#########################################
n0<-600
n1<-400
n<-n0+n1


X0<-cbind(rnorm(n0,3,1),rnorm(n0,3,1),rnorm(n0,2,0.5),rnorm(n0,2.5,0.5),rnorm(n0,3.5,1))
X1<-cbind(rnorm(n1,4,1),rnorm(n1,4,1),rnorm(n1,3,0.5),rnorm(n1,3.5,0.5),rnorm(n1,4.5,1))


#Simulamos el valor de la Y

beta0<-c(1,-1,-0.5,1,-1,1)
beta1<-c(1,-0.4,1,-1,1,-0.5)
Y0<-1*(exp(beta0[1]+as.matrix(X0)%*%beta0[2:6])/(1+exp(beta0[1]+as.matrix(X0)%*%beta0[2:6]))>0.5)
Y1<-1*(exp(beta1[1]+as.matrix(X1)%*%beta1[2:6])/(1+exp(beta1[1]+as.matrix(X1)%*%beta1[2:6]))>0.5)
X<-rbind(X0,X1)
Y<-c(Y0,Y1)
S<-c(rep(0,n0), rep(1,n1))
#Y<-sample(c(0,1),n,replace=T, prob=c(0.5,0.5))
data_original<-as.data.frame(cbind(X,S,Y))


# #REGLA QUE ENTRENAMOS PARA SIMULAR EL RESULTADO DE LA CLASIFICACION DADA
# #learning sample
# n_sim<-floor(0.3*n)
# datasim<-data_original[sample(c(1:n), size=n_sim, replace=F),c(1:5,6)]
# logit_sim <- glm(S ~., family='binomial', data=datasim)
# #summary(logit)
# 
# #add a column with the prediction of the logit model for the subsample
# Y.pred_0<-1*(predict(logit_sim, newdata = data_original[S==0,c(1:5,6)], type="response")>0.6)
# Y.pred_1<-1*(predict(logit_sim, newdata = data_original[S==1,c(1:5,6)], type="response")>0.4)
# data_original$Y.pred<-rep(2,n)
# data_original[S==0,]$Y.pred<-Y.pred_0
# data_original[S==1,]$Y.pred<-Y.pred_1

#REGLA QUE ENTRENAMOS PARA PREDECIR
#learning sample
n_learn<-floor(0.3*n)
index_learn<-sample(c(1:n), size=n_learn, replace=F)
dataLearn<-data_original[index_learn,c(1:5,7)]
logit <- glm(Y ~., family='binomial', data=dataLearn)
#summary(logit)

#add a column with the prediction of the logit model for the subsample
data_study<-data_original[-index_learn,]
data_study$Y.logit<-1*(predict(logit, newdata = data_study[,c(1:5,7)], type="response")>0.5)


#DI with respect to the Y.pred (logit dependiendo de S)
#DI0.pred<-disparate(data_original,S = 6,Y = 7,0.05)
#DI with respect to the Ylogit
DI0.logit<-disparate(data_study,S = 6,Y = 8,0.05)
#confusion matrix Ytrue(column $income) vs. Ylogit
conf_matrix<-table(data_study$Y,data_study$Y.logit)
#estimated probability or error with the original data and the logit model vs the true Y
error0.logit<-(conf_matrix[1,2]+conf_matrix[2,1])/sum(conf_matrix)

##############
#Random forest
##############
# library(ranger)
# # #####################################
# # #we use the whole dataset -index_learn
# # 
#  rf.fit=ranger(as.factor(Y.pred)~., data=dataLearn) 
#  pred.rf=predict(rf.fit,data_study[,c(1:5,7)])
#  data_study$Y.rf<-as.numeric(pred.rf$predictions)-1
#  DI0.rf<-disparate(data_study,6,9,0.05)
#  conf_mat.rf<-table(data_study$Y,data_study$Y.rf)
#  error0.rf<-(conf_mat.rf[1,2]+conf_mat.rf[2,1])/sum(conf_mat.rf)

#The repair is only done for the testing sample
X0<-data_study[data_study$S==0,c(1:5)]
X1<-data_study[data_study$S==1,c(1:5)]
n0_test<-dim(X0)[1]
n1_test<-dim(X1)[1]
pi0<-n0_test/(n0_test+n1_test)
pi1<-1-pi0
##########
#transport
costes<-matrix(0,ncol=n1_test,nrow=n0_test)
for (i in 1:n0_test){for (j in 1:n1_test) costes[i,j]<-sum((X0[i,]-X1[j,])^2)}
t.opt<-lp.transport(costes,"min",row.rhs=rep(1/n0_test,n0_test),row.signs=rep("=",n0_test),col.rhs=rep(1/n1_test,n1_test),col.signs=rep('=',n1_test),integers=NULL)
##########

###########################################

#1) TOTAL REPAIR as described in Figure 4 left

###########################################

X0_repair<-as.data.frame(pi0*X0+pi1*n0_test*(t.opt$solution%*%as.matrix(X1)))
X0_repair$protected<-rep(0,n0_test)
X1_repair<-as.data.frame(pi0*n1_test*(t(t.opt$solution)%*%as.matrix(X0))+pi1*X1)
X1_repair$protected<-rep(1,n1_test)
data_repair<-rbind(X0_repair, X1_repair)
#names(data_repair)[1]<-"visible"
#predict Y with the logit model
#totalrepair_prediction<-predict(logit, newdata = data_repair, type="response")
data_repair$Y.logit<-1*(predict(logit, newdata = data_repair, type="response")>0.5)
#predict Y with the rf model
#data_repair$Y.rf<-as.numeric(predict(rf.fit, data_repair)$predictions)-1

#sort repaired data as in the original subsample data frame
#data_repair<-data_repair[rownames(data_original),]

#DI with the totally repaired data
#we observe that this procedure does not achieve full impredictability, DI!=1
#logit DI
DI_TR.logit<-disparate(data_repair,S = 6, Y = 7,0.05)
#RF DI
#DI_TR.rf<-disparate(data_repair,2,3,0.05)

#perdida en la prediccion logit, con respecto a Y verdadero
conf_matrix_TRepair<-table(data_study$Y, data_repair$Y.logit)
#conf_matrix_TRepair<-rbind(c(sum(data_original$Y==0),sum(data_original$Y==1)),c(sum(data_repair$Y.logit==0),sum(data_repair$Y.logit==1)))
error_TRepair.logit<-(conf_matrix_TRepair[1,2]+conf_matrix_TRepair[2,1])/sum(conf_matrix_TRepair)
error_TRepair.logit
#perdida en la prediccion RF, con respecto a Y verdadero
#conf_matrix_TRepair<-rbind(c(sum(data_original$Y==0),sum(data_original$Y==1)),c(sum(data_repair$Y.rf==0),sum(data_repair$Y.rf==1)))
#conf_matrix_TRepair<-table(data_original$Y, data_repair$Y.rf)
#error_TRepair.rf<-(conf_matrix_TRepair[1,2]+conf_matrix_TRepair[2,1])/sum(conf_matrix_TRepair)
#error_TRepair.rf
#
##########

#####################################

#2)TOTAL REPAIR : FULL IMPREDICTABILITY (Figure 4 right)

#####################################
id_X0<-c(1:n0_test)
id_X1<-c((n0_test+1):(n0_test+n1_test))
#matrix with the non zero positions of the transport solution
index_positive<-which(t.opt$solution!=0, arr.ind = T)
X_repair<-as.data.frame(pi0*X0[index_positive[,1],]+pi1*X1[index_positive[,2],])
X_repair$index_X0<-id_X0[index_positive[,1]]
X_repair$index_X1<-id_X1[index_positive[,2]]
X_repair$weight<-t.opt$solution[index_positive]
#names(X_repair)[1]<-"visible"
X_repair_ext<-rbind(X_repair,X_repair)
X_repair_ext$protected<-c(rep(0,nrow(index_positive)),rep(1,nrow(index_positive)))

#predict Y with the logit model
X_repair_ext$Y.logit<-1*(predict(logit, newdata = X_repair_ext, type="response")>0.5)
#predict Y with the rf model
#X_repair_ext$Y.rf<-as.numeric(predict(rf.fit, X_repair_ext)$predictions)-1


#DI logit with the totally repaired data with full impredictability
DI_FI.logit<-disparate(X_repair_ext,9,10,0.05)

#DI rf with the totally repaired data with full impredictability
#DI_FI.rf<-disparate(X_repair_ext,5,6,0.05)


#extended array with the true values of Y
Y_ext<-c(data_study[X_repair$index_X0,]$Y,data_study[X_repair$index_X1,]$Y)
#loss logit
conf_matrix_FIRepair<-table(Y_ext,X_repair_ext$Y.logit)
#conf_matrix_FIRepair<-rbind(c(sum(Y_ext==0),sum(Y_ext==1)),c(sum(X_repair_ext$Y.logit==0),sum(X_repair_ext$Y.logit==1)))
error_FIRepair.logit<-(conf_matrix_FIRepair[1,2]+conf_matrix_FIRepair[2,1])/sum(conf_matrix_FIRepair)
error_FIRepair.logit
#loss RF
#???conf_matrix_FIRepair<-table(Y_ext, X_repair_ext$Y.rf)
#conf_matrix_FIRepair<-rbind(c(sum(Y_ext==0),sum(Y_ext==1)),c(sum(X_repair_ext$Y.rf==0),sum(X_repair_ext$Y.rf==1)))
#error_FIRepair.rf<-(conf_matrix_FIRepair[1,2]+conf_matrix_FIRepair[2,1])/sum(conf_matrix_FIRepair)
#error_FIRepair.rf


##########
#####################################

  
#####################################

#4)PARTIAL REPAIR, Feldmann approach

#####################################
#tau<-disparate(data_original,S = 3,Y = 4,0.05)[2]

lambda<-seq(0,1,by=0.025)
i<-1
DI.logit<-vector()
error_PRepair.logit<-vector()
#DI.rf<-vector()
#error_PRepair.rf<-vector()
for(k in 1:length(lambda)){
  # if(sum(as.character(lambda_ev)=="1")==1) break
  
  X0_P_repair<-as.data.frame((1-lambda[k]*pi1)*X0+lambda[k]*pi1*(n0_test*t.opt$solution%*%as.matrix(X1)))
  X1_P_repair<-as.data.frame(lambda[k]*pi0*(n1_test*t(t.opt$solution)%*%as.matrix(X0))+(1-lambda[k]*pi0)*X1)
# names(X0_P_repair)<-"visible"
# names(X1_P_repair)<-"visible"
   data_P_repair<-rbind(X0_P_repair, X1_P_repair)
  data_P_repair$protected<-c(rep(0,n0_test), rep(1,n1_test))
  #data_P_repair$Y.prediction<-predict(logit, newdata = data_P_repair, type="response")
  data_P_repair$Y.logit<-1*(predict(logit, newdata = data_P_repair, type="response")>0.5)
  
  #data_P_repair$Y.rf<-as.numeric(predict(rf.fit, data_P_repair)$predictions)-1
  DI_PR.logit<-disparate(data_P_repair,6,7,0.05)
  #DI_PR.rf<-disparate(data_P_repair,2,3,0.05)
  
  #sort repaired data as in the original subsample data frame
  #data_P_repair<-data_P_repair[rownames(data_original),]
  
  #perdida en la prediccion logit
  conf_matrix_PRepair<-table(data_study$Y, data_P_repair$Y.logit)
  # conf_matrix_PRepair<-rbind(c(sum(data_original$Y==0),sum(data_original$Y==1)),c(sum(data_P_repair$Y.logit==0),sum(data_P_repair$Y.logit==1)))
   error_PRepair.logit[k]<-(conf_matrix_PRepair[1,2]+conf_matrix_PRepair[2,1])/sum(conf_matrix_PRepair)
   DI.logit<-rbind(DI.logit,DI_PR.logit)
  # 
  #perdida en la prediccion RF
  #conf_matrix_PRepair<-table(data_original$Y, data_P_repair$Y.rf)
 # conf_matrix_PRepair<-rbind(c(sum(data_original$Y==0),sum(data_original$Y==1)),c(sum(data_P_repair$Y.rf==0),sum(data_P_repair$Y.rf==1)))
  #error_PRepair.rf[k]<-(conf_matrix_PRepair[1,2]+conf_matrix_PRepair[2,1])/sum(conf_matrix_PRepair)
  #DI.rf<-rbind(DI.rf,DI_PR.rf)
  
 # tau<-DI_PR.logit[2]
#  tau<-DI_PR.rf[2]
}

# 
rownames(DI.logit)<-as.character(c(1:nrow(DI.logit)))
DI.logit<-as.data.frame(DI.logit)
names(DI.logit)<-c("linf","DI_hat","lsup", "BER")
DI.logit$lambda<-lambda
DI.logit$error<-error_PRepair.logit[1:nrow(DI.logit)]

ggplot(DI.logit, aes(x=DI_hat, y=BER))+geom_point(size=2)+geom_line(color="red", lwd=1)+scale_x_continuous()+scale_y_continuous()

p1<-ggplot(DI.logit, aes(x=lambda, y=DI_hat))+geom_point(size=2)+geom_line(color="red", lwd=1)+geom_ribbon(aes(ymin=linf,ymax=lsup),alpha=0.3, fill="blue")
p1<-p1+scale_x_continuous(name ="Amount of repair", breaks=lambda)+scale_y_continuous(name ="Disparate Impact", breaks=seq(0,1.5,by=0.1))+ggtitle("95% Confidence Interval for DI of Logistic Regression")
p1+geom_hline(yintercept = c(0.8,1,DI_TR.logit[2]), lty=2)

p2<-ggplot(DI.logit, aes(x=lambda, y=error))+geom_point(size=2)+geom_line(color="green", lwd=1)+scale_x_continuous(name ="Amount of repair", breaks=lambda)+scale_y_continuous()
p2

p3<-ggplot(DI.logit, aes(x=lambda, y=BER))+geom_line(color="green")+scale_x_continuous(name ="Amount of repair", breaks=lambda)
p3+geom_hline(yintercept = 0.5, lty=2)


###RANDOM FOREST

# rownames(DI.rf)<-as.character(c(1:nrow(DI.rf)))
# DI.rf<-as.data.frame(DI.rf)
# names(DI.rf)<-c("linf","DI_hat","lsup", "BER")
# DI.rf$lambda<-lambda
# DI.rf$error<-error_PRepair.rf[1:nrow(DI.rf)]
# 
# ggplot(DI.rf, aes(x=DI_hat, y=BER))+geom_point(size=2)+geom_line(color="red", lwd=1)+scale_x_continuous()+scale_y_continuous()
# 
# p1<-ggplot(DI.rf, aes(x=lambda, y=DI_hat))+geom_point(size=2)+geom_line(color="red", lwd=1)+geom_ribbon(aes(ymin=linf,ymax=lsup),alpha=0.3, fill="blue")
# p1<-p1+scale_x_continuous(name ="Amount of repair", breaks=lambda)+scale_y_continuous(name ="Disparate Impact", breaks=seq(0,1.5,by=0.1))+ggtitle("95% Confidence Interval for DI of Random Forest")
# p1+geom_hline(yintercept = c(0.8,1,DI_TR.rf[2]), lty=2)
# 
# p2<-ggplot(DI.rf, aes(x=lambda, y=error))+geom_point(size=2)+geom_line(color="green", lwd=1)+scale_x_continuous(name ="Amount of repair", breaks=lambda)+scale_y_continuous()
# 
# p3<-ggplot(DI.rf, aes(x=lambda, y=BER))+geom_point(size=2)+geom_line(color="red", lwd=1)+scale_x_continuous(name ="Amount of repair", breaks=lambda)
# p3+geom_hline(yintercept = 0.5, lty=2)

###############################################

#6) RANDOM REPAIR FULL IMPREDICTABILITY

#####################################

library(Rlab)

##logit

DI_RR_L.logit<-list()
for(k in 1:100){
  lambda<-seq(0,1,by=0.025)
  i<-1
  DI_RR.logit<-vector()
  error_RR.logit<-vector()
  for(r in 1:length(lambda)){
    B0<-rbern(n0_test,lambda[r])
    B1<-rbern(n1_test,lambda[r])
    # x0[B0==0,] and x1[B1==0,] are the instances that will not be repaired
    X0_repair<-as.data.frame(pi0*X0[index_positive[index_positive[,1]%in%which(B0==1),1],]+pi1*X1[index_positive[index_positive[,1]%in%which(B0==1),2],])
   # names(X0_repair)<-"visible"
    X0_repair$names<-as.character(id_X0[index_positive[index_positive[,1]%in%which(B0==1),1]])
    X1_repair<-as.data.frame(pi0*X0[index_positive[index_positive[,2]%in%which(B1==1),1],]+pi1*X1[index_positive[index_positive[,2]%in%which(B1==1),2],])
    #names(X1_repair)<-"visible"
    X1_repair$names<-as.character(id_X1[index_positive[index_positive[,2]%in%which(B1==1),2]])
    X0_no_repair<-as.data.frame(X0[B0==0,])
    #names(X0_no_repair)<-"visible"
    X0_no_repair$names<-as.character(id_X0[B0==0])
    X1_no_repair<-as.data.frame(X1[B1==0,])
    #names(X1_no_repair)<-"visible"
    X1_no_repair$names<-as.character(id_X1[B1==0])
    X_repair<-rbind(X0_no_repair,X0_repair,X1_no_repair,X1_repair)
    X_repair$protected<-c(rep(0,nrow(X0_no_repair) + nrow(X0_repair)), rep(1,nrow(X1_no_repair)+ nrow(X1_repair)))
    
    X_repair$Y.logit<-1*(predict(logit, newdata = X_repair, type="response")>0.5)
    
    DI_RR<-disparate(X_repair,7,8,0.05)
    
    #perdida logit
    Y_ext<-data_study[as.numeric(X_repair$names),]$Y
    conf_matrix_RRepair<-table(Y_ext, X_repair$Y.logit)
    #conf_matrix_RRepair<-rbind(c(sum(Y_ext==0),sum(Y_ext==1)),c(sum(X_repair$Y.rf==0),sum(X_repair$Y.rf==1)))
    error_RR.logit[i]<-(conf_matrix_RRepair[1,2]+conf_matrix_RRepair[2,1])/sum(conf_matrix_RRepair)
    
    
    #error_RRepair<-(conf_matrix_RRepair[1,2]+conf_matrix_RRepair[2,1])/sum(conf_matrix_RRepair)
    DI_RR.logit<-rbind(DI_RR.logit,DI_RR)
    
    i<-i+1
    
  }
  
  rownames(DI_RR.logit)<-as.character(c(1:nrow(DI_RR.logit)))
  DI_RR.logit<-as.data.frame(DI_RR.logit)
  names(DI_RR.logit)<-c("linf","DI_hat","lsup", "BER")
  DI_RR.logit$lambda<-lambda
  DI_RR.logit$error<-error_RR.logit
  
  ###
  # ggplot(DI_RR.rf, aes(x=DI_hat, y=BER))+geom_point(size=2)+geom_line(color="red", lwd=1)+scale_x_continuous()+scale_y_continuous()
  # 
  # p1<-ggplot(DI_RR.rf, aes(x=lambda, y=DI_hat))+geom_point(size=2)+geom_line(color="red", lwd=1)+geom_ribbon(aes(ymin=linf,ymax=lsup),alpha=0.3, fill="blue")
  # p1<-p1+scale_x_continuous(name ="Amount of repair", breaks=lambda)+scale_y_continuous(name ="Disparate Impact", breaks=seq(0,1.5,by=0.1))+ggtitle("95% Confidence Interval for DI of Random Forest")
  # p1+geom_hline(yintercept = c(0.8,1,DI_TR.rf[2]), lty=2)
  # 
  # p2<-ggplot(DI_RR.rf, aes(x=lambda, y=error))+geom_point(size=2)+geom_line(color="green", lwd=1)+scale_x_continuous(name ="Amount of repair", breaks=lambda)+scale_y_continuous()
  # 
  # p3<-ggplot(DI_RR.rf, aes(x=lambda, y=BER))+geom_line(color="green")+scale_x_continuous(name ="Amount of repair", breaks=lambda)
  # p3+geom_hline(yintercept = 0.5, lty=2)
  ###
  
  DI_RR_L.logit[[k]]<-DI_RR.logit
}

LINF<-lapply(DI_RR_L.logit, "[[", 1)[[1]]
DIHAT<-lapply(DI_RR_L.logit, "[[", 2)[[1]]
LSUP<-lapply(DI_RR_L.logit, "[[", 3)[[1]]
BER<-lapply(DI_RR_L.logit, "[[", 4)[[1]]
ERROR<-lapply(DI_RR_L.logit, "[[", 6)[[1]]
for(i in 2:100){
  LINF<-cbind(LINF,lapply(DI_RR_L.logit, "[[", 1)[[i]])
  DIHAT<-cbind(DIHAT,lapply(DI_RR_L.logit, "[[", 2)[[i]])
  LSUP<-cbind(LSUP,lapply(DI_RR_L.logit, "[[", 3)[[i]])
  BER<-cbind(BER,lapply(DI_RR_L.logit, "[[", 4)[[i]])
  ERROR<-cbind(ERROR,lapply(DI_RR_L.logit, "[[", 6)[[i]])
}
mLINF<-apply(LINF,MARGIN = 1, FUN = mean)
mDIHAT<-apply(DIHAT,MARGIN = 1, FUN = mean)
mLSUP<-apply(LSUP,MARGIN = 1, FUN = mean)
mBER<-apply(BER,MARGIN = 1, FUN = mean)
mERROR<-apply(ERROR,MARGIN = 1, FUN = mean)
mDI_RR.logit<-as.data.frame(cbind(mLINF,mDIHAT,mLSUP,mBER,mERROR,seq(0,1,by=0.025)))
names(mDI_RR.logit)[6]<-"lambda"

ggplot(mDI_RR.logit, aes(x=mDIHAT, y=mBER))+geom_point(size=2)+geom_line(color="red", lwd=1)+scale_x_continuous()+scale_y_continuous()

p1<-ggplot(mDI_RR.logit, aes(x=lambda, y=mDIHAT))+geom_point(size=1)+geom_line(color="red", lwd=1)+geom_ribbon(aes(ymin=mLINF,ymax=mLSUP),alpha=0.3, fill="blue")
p1<-p1+scale_x_continuous(name ="Amount of repair", breaks=seq(0,1,by=0.025))+scale_y_continuous(name ="Disparate Impact")+ggtitle("95% Confidence Interval for DI of logit")
p1+geom_hline(yintercept = c(0.8,1), lty=2)
p2<-ggplot(mDI_RR.logit, aes(x=lambda, y=mERROR))+geom_point(size=1)+geom_line(color="red", lwd=1)+scale_x_continuous(name ="Amount of repair", breaks=lambda)+scale_y_continuous(name ="error")+ggtitle("Classification Error with logit")
p2+geom_hline(yintercept = c(error_FIRepair.logit,error0.logit), lty=2, color=c("black"))
p3<-ggplot(mDI_RR.logit, aes(x=lambda, y=mBER))+geom_point(size=1)+geom_line(color="red", lwd=1)+scale_x_continuous(name ="Amount of repair", breaks=lambda)
p3+geom_hline(yintercept = 0.5, lty=2)

#plotting together
DI.logit$repair<-"GR"
mDI_RR.logit$repair<-"RR"
mDI_RR.logit<-mDI_RR.logit[,c(1,2,3,4,6,5,7)]
names(mDI_RR.logit)<-names(DI.logit)
TT<-rbind(DI.logit,mDI_RR.logit)

ggplot(TT, aes(x=DI_hat, y=BER, color=repair))+geom_point(size=2)+geom_line(color="red", lwd=1)+scale_x_continuous()+scale_y_continuous()

p1<-ggplot(TT, aes(x=lambda, y=DI_hat,color=repair))+geom_point(size=2)+geom_line(aes(color=repair), lwd=1)+geom_ribbon(aes(ymin=linf,ymax=lsup, fill=repair),alpha=0.3)
p1<-p1+scale_x_continuous(name ="Amount of repair", breaks=lambda)+scale_y_continuous(name ="Disparate Impact", breaks=seq(0,1.5,by=0.1))+ggtitle("95% Confidence Interval for DI of Logistic Regression")
p1+geom_hline(yintercept = c(0.8,1,DI_TR.logit[2]), lty=2)

p2<-ggplot(TT, aes(x=lambda, y=error, color=repair))+geom_point(size=2)+geom_line(aes(color=repair), lwd=1)+scale_x_continuous(name ="Amount of repair", breaks=lambda)+scale_y_continuous()
p2+geom_hline(yintercept = c(error0.logit,error_FIRepair.logit,error_TRepair.logit), lty=2)

p3<-ggplot(TT, aes(x=lambda, y=BER, color=repair))+geom_point(size=2)+geom_line(aes(color=repair), lwd=1)+scale_x_continuous(name ="Amount of repair", breaks=lambda)
p3+geom_hline(yintercept = 0.5, lty=2)




# ##ranfom forests
# 
# DI_RR_L.rf<-list()
# for(k in 1:100){
#   lambda<-seq(0,1,by=0.025)
#   i<-1
#   DI_RR.rf<-vector()
#   error_RR.rf<-vector()
#   for(r in 1:length(lambda)){
#     B0<-rbern(n0,lambda[r])
#     B1<-rbern(n1,lambda[r])
#     # x0[B0==0,] and x1[B1==0,] are the instances that will not be repaired
#     X0_repair<-as.data.frame(pi0*X0[index_positive[index_positive[,1]%in%which(B0==1),1]]+pi1*X1[index_positive[index_positive[,1]%in%which(B0==1),2]])
#     names(X0_repair)<-"visible"
#     X0_repair$names<-as.character(id_X0[index_positive[index_positive[,1]%in%which(B0==1),1]])
#     X1_repair<-as.data.frame(pi0*X0[index_positive[index_positive[,2]%in%which(B1==1),1]]+pi1*X1[index_positive[index_positive[,2]%in%which(B1==1),2]])
#     names(X1_repair)<-"visible"
#     X1_repair$names<-as.character(id_X1[index_positive[index_positive[,2]%in%which(B1==1),2]])
#     X0_no_repair<-as.data.frame(X0[B0==0])
#     names(X0_no_repair)<-"visible"
#     X0_no_repair$names<-as.character(id_X0[B0==0])
#     X1_no_repair<-as.data.frame(X1[B1==0])
#     names(X1_no_repair)<-"visible"
#     X1_no_repair$names<-as.character(id_X1[B1==0])
#     X_repair<-rbind(X0_no_repair,X0_repair,X1_no_repair,X1_repair)
#     X_repair$protected<-c(rep(0,nrow(X0_no_repair) + nrow(X0_repair)), rep(1,nrow(X1_no_repair)+ nrow(X1_repair)))
#     
#     X_repair$Y.rf<-as.numeric(predict(rf.fit, X_repair)$predictions)-1
#     
#   
#     DI_RR<-disparate(X_repair,3,4,0.05)
#     
#     #perdida rf
#     Y_ext<-data_original[X_repair$names,]$Y
#     conf_matrix_RRepair<-table(Y_ext, X_repair$Y.rf)
#     #conf_matrix_RRepair<-rbind(c(sum(Y_ext==0),sum(Y_ext==1)),c(sum(X_repair$Y.rf==0),sum(X_repair$Y.rf==1)))
#     error_RR.rf[i]<-(conf_matrix_RRepair[1,2]+conf_matrix_RRepair[2,1])/sum(conf_matrix_RRepair)
#     
#     
#     #error_RRepair<-(conf_matrix_RRepair[1,2]+conf_matrix_RRepair[2,1])/sum(conf_matrix_RRepair)
#     DI_RR.rf<-rbind(DI_RR.rf,DI_RR)
#     
#     i<-i+1
# 
#   }
#   
#   rownames(DI_RR.rf)<-as.character(c(1:nrow(DI_RR.rf)))
#   DI_RR.rf<-as.data.frame(DI_RR.rf)
#   names(DI_RR.rf)<-c("linf","DI_hat","lsup", "BER")
#   DI_RR.rf$lambda<-lambda
#   DI_RR.rf$error<-error_RR.rf
#   
#   ###
#   # ggplot(DI_RR.rf, aes(x=DI_hat, y=BER))+geom_point(size=2)+geom_line(color="red", lwd=1)+scale_x_continuous()+scale_y_continuous()
#   # 
#   # p1<-ggplot(DI_RR.rf, aes(x=lambda, y=DI_hat))+geom_point(size=2)+geom_line(color="red", lwd=1)+geom_ribbon(aes(ymin=linf,ymax=lsup),alpha=0.3, fill="blue")
#   # p1<-p1+scale_x_continuous(name ="Amount of repair", breaks=lambda)+scale_y_continuous(name ="Disparate Impact", breaks=seq(0,1.5,by=0.1))+ggtitle("95% Confidence Interval for DI of Random Forest")
#   # p1+geom_hline(yintercept = c(0.8,1,DI_TR.rf[2]), lty=2)
#   # 
#   # p2<-ggplot(DI_RR.rf, aes(x=lambda, y=error))+geom_point(size=2)+geom_line(color="green", lwd=1)+scale_x_continuous(name ="Amount of repair", breaks=lambda)+scale_y_continuous()
#   # 
#   # p3<-ggplot(DI_RR.rf, aes(x=lambda, y=BER))+geom_line(color="green")+scale_x_continuous(name ="Amount of repair", breaks=lambda)
#   # p3+geom_hline(yintercept = 0.5, lty=2)
#   ###
#   
#   DI_RR_L.rf[[k]]<-DI_RR.rf
# }
# 
# LINF<-lapply(DI_RR_L.rf, "[[", 1)[[1]]
# DIHAT<-lapply(DI_RR_L.rf, "[[", 2)[[1]]
# LSUP<-lapply(DI_RR_L.rf, "[[", 3)[[1]]
# BER<-lapply(DI_RR_L.rf, "[[", 4)[[1]]
# ERROR<-lapply(DI_RR_L.rf, "[[", 6)[[1]]
# for(i in 2:100){
#   LINF<-cbind(LINF,lapply(DI_RR_L.rf, "[[", 1)[[i]])
#   DIHAT<-cbind(DIHAT,lapply(DI_RR_L.rf, "[[", 2)[[i]])
#   LSUP<-cbind(LSUP,lapply(DI_RR_L.rf, "[[", 3)[[i]])
#   BER<-cbind(BER,lapply(DI_RR_L.rf, "[[", 4)[[i]])
#   ERROR<-cbind(ERROR,lapply(DI_RR_L.rf, "[[", 6)[[i]])
# }
# mLINF<-apply(LINF,MARGIN = 1, FUN = mean)
# mDIHAT<-apply(DIHAT,MARGIN = 1, FUN = mean)
# mLSUP<-apply(LSUP,MARGIN = 1, FUN = mean)
# mBER<-apply(BER,MARGIN = 1, FUN = mean)
# mERROR<-apply(ERROR,MARGIN = 1, FUN = mean)
# mDI_RR.rf<-as.data.frame(cbind(mLINF,mDIHAT,mLSUP,mBER,mERROR,seq(0,1,by=0.025)))
# names(mDI_RR.rf)[6]<-"lambda"
# 
# ggplot(mDI_RR.rf, aes(x=mDIHAT, y=mBER))+geom_point(size=2)+geom_line(color="red", lwd=1)+scale_x_continuous()+scale_y_continuous()
# 
# p1<-ggplot(mDI_RR.rf, aes(x=lambda, y=mDIHAT))+geom_point(size=1)+geom_line(color="red", lwd=1)+geom_ribbon(aes(ymin=mLINF,ymax=mLSUP),alpha=0.3, fill="blue")
# p1<-p1+scale_x_continuous(name ="Amount of repair", breaks=seq(0,1,by=0.025))+scale_y_continuous(name ="Disparate Impact")+ggtitle("95% Confidence Interval for DI of Random Forests")
# p1+geom_hline(yintercept = c(0.8,1), lty=2)
# p2<-ggplot(mDI_RR.rf, aes(x=lambda, y=mERROR))+geom_point(size=1)+geom_line(color="red", lwd=1)+scale_x_continuous(name ="Amount of repair", breaks=lambda)+scale_y_continuous(name ="error")+ggtitle("Classification Error with Random Forests")
# p2+geom_hline(yintercept = c(error_FIRepair.rf,error0.rf), lty=2, color=c("black"))
# p3<-ggplot(mDI_RR.rf, aes(x=lambda, y=mBER))+geom_point(size=1)+geom_line(color="red", lwd=1)+scale_x_continuous(name ="Amount of repair", breaks=lambda)
# p3+geom_hline(yintercept = 0.5, lty=2)
# 
# 
