p1=ggplot(comb,aes(x=-Estimate_ds,y=Estimate_sn,color=CellType))+geom_point()+xlab("Effect size mixed model (glmmTMB)")+ylab("Effect size DAESC_bb")+geom_abline(slope=1,linetype="dotted")
p2=ggplot(comb,aes(x=-Estimate_ds,y=Estimate_pb,color=CellType))+geom_point()+xlab("Effect size pseudobulk")+ylab("Effect size DAESC_bb")+geom_abline(slope=1,linetype="dotted")
p3=ggplot(comb,aes(x=Estimate_sn,y=Estimate_pb,color=CellType))+geom_point()+xlab("Effect size pseudobulk")+ylab("Effect size mixed model (glmmTMB)")+geom_abline(slope=1,linetype="dotted")
p1 | p2 | p3 | plot_layout(guides = "collect") 


p1=ggplot(comb,aes(x=-log(pval_sn,10),y=-log(pval_ds,10),color=CellType))+geom_point()+xlab("-log10 p-values mixed model (glmmTMB)")+ylab("-log10 p-values DAESC_bb")+geom_abline(slope=1,linetype="dotted")
p2=ggplot(comb,aes(x=-log(pval_pb,10),y=-log(pval_ds,10),color=CellType))+geom_point()+xlab("-log10 p-values pseudobulk")+ylab("-log10 p-values DAESC_bb")+geom_abline(slope=1,linetype="dotted")
p3=ggplot(comb,aes(x=-log(pval_pb,10),y=-log(pval_sn,10),color=CellType))+geom_point()+xlab("-log10 p-values pseudobulk")+ylab("-log10 p-values mixed model (glmmTMB)")+geom_abline(slope=1,linetype="dotted")
p1 | p2 | p3 | plot_layout(guides = "collect")

ggplot(dat,aes(x=method,y=time,color=CellType))+geom_point()+scale_y_log10()+theme(axis.text.x = element_text(angle = 45,vjust=.5))+xlab("")+ylab("Runtime (seconds)")
