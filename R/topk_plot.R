library(ggplot2)

n=10

se_u = rep(0,15)
se_o = rep(0,15)
sp_u = rep(0,15)
sp_o = rep(0,15)

for (k in seq(1,15)){
  ses=rep(0,n)
  sps=rep(0,n)
  for (i in seq(1,n)){
    a=read.table(paste('topk/te',k,'_',i,'.csv', sep = ''), sep = ",")
    ses[i]=a[1,1]/sum(a[,1])
    sps[i]=a[2,2]/sum(a[,2])
  }
  se_u[k]=mean(ses)
  se_o[k]=sd(ses)
  sp_u[k]=mean(sps)
  sp_o[k]=sd(sps)
}
# Default line plot
df<-data.frame("Mean_sp"=sp_u,"Mean_se"=se_u,"Std_sp"=sp_o,"Std_se"=se_o,"num_genes"=1:n)
ggplot(df, aes(x=num_genes, y=Mean_sp, title("Specificity"))) + 
  geom_errorbar(aes(ymin=Mean_sp-Std_sp, ymax=Mean_sp+Std_sp), width=.1) +
  geom_line() +
  geom_point() +
  ggtitle("Specificity")
ggplot(df, aes(x=num_genes, y=Mean_se)) + 
  geom_errorbar(aes(ymin=Mean_se-Std_se, ymax=Mean_se+Std_se), width=.1) +
  geom_line() +
  geom_point() +
  ggtitle("Sensitivity")

df<-data.frame("Acc"=c(se_u,sp_u),"num_genes"=1:n,"err"=c(se_o,sp_o),group=c(rep("Sensitivity",n),rep("Specificity",n)))
ggplot(df, aes(x=num_genes, y=Acc, color=group)) + 
  geom_errorbar(aes(ymin=Acc-err, ymax=Acc+err), width=.1) +
  geom_line() +
  geom_point()