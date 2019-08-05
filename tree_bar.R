library(ape)
library(phylobase)
library(reshape2)
library(ggplot2)



sig<- read.delim("all_patients_DP_contribution.txt", sep="\t", stringsAsFactors = F)
sig$sampleID<- rownames(sig)
tree2<- read.delim("summary_fro_timing.txt", sep="\t", header=T, stringsAsFactors = F) 

melter <- function(stri) {

s=sig[grep(stri,sig$sampleID),]
s$SBS1=s$Signature.Subs.01*s$mutations
s$SBS2=s$Signature.Subs.02*s$mutations
s$SBS5=s$Signature.Subs.05*s$mutations
s$SBS8=s$Signature.Subs.08*s$mutations
s$SBS9=s$Signature.Subs.09*s$mutations
s$SBS13=s$Signature.Subs.13*s$mutations
s$SBS18=s$Signature.Subs.18*s$mutations
s$SBS.MM1=s$MM1*s$mutations
s=s[,c("sampleID","SBS1","SBS2","SBS5","SBS8","SBS9","SBS13","SBS18","SBS.MM1")]

m=melt(s)
head(m)
pdf(sprintf("~/Desktop/%s_bargg.pdf",stri), useDingbats=FALSE)
p = ggplot(m,aes(sampleID,value,fill=variable))+geom_bar(stat="identity")+theme_minimal()+scale_fill_manual(values = col)
print(p)
#ggsave(p, file=sprintf("~/Desktop/%s_bargg.pdf",stri))
dev.off()

return(m)
}



m = melter("PD26400")
stri="PD26400"


a <- "(a:1,b:1);"
t <- read.tree(text="(((PD26400.2:504,B:1)PD26400.3:135,PD26400.na12:12),(PD26400.5:510,PD26400.na:30)PD26400.4:991):10;")
a = "((t3:0.2,t4:0.2):0.6,(t2:0.7,(t5:0.2,t6:0.2):0.5):0.1):0.2;"



a = "((PD26400.5:510, PD26400.na10:22):991, ((PD26400.1:235,P:2)PD26400.2:504, PD26400.12:4):135)PD26400.6:5560;"
#a = "(dog:20, (elephant:30, horse:60):20):50;"
tr <- read.tree(text = a)

#pdf(sprintf("~/Desktop/%s_tree.pdf",stri), useDingbats=FALSE)
plot(tr,show.node.label = T,, root.edge=T, edge.width = 3)
#axisPhylo()
axisPhylo(1,root.time=0,  backward=FALSE)

#dev.off()



a = "(((PD26400.5:510, PD26400.na10:22):991, ((PD26400.1:235,P:2)PD26400.2:504, PD26400.12:4):135)PD26400.6:5560)PPP:3;"
tr <- read.tree(text = a)

pdf("~/Desktop/PD26400.pdf", useDingbats=FALSE)
plot(tr,show.node.label = T,, root.edge=T, edge.width = 3)
#axisPhylo()
axisPhylo(1,root.time=0,  backward=FALSE)

dev.off()

############################# stri="PD26403" ############################
stri="PD26403"
m = melter(stri)
sig[grep(stri,sig$sampleID),]
tree2[grep(stri,tree2$patient),c("tree_structure")]


a = "((PD26403.3:282,(PD26403.2:1248,PD26403.10:1)PD26403.4:173,PD26403.1:3105)PD26403.7:5678)PPP:3;"
tr <- read.tree(text = a)

pdf(sprintf("~/Desktop/%s_tree.pdf",stri), useDingbats=FALSE,width=10,height=5)
plot(tr,show.node.label = T, root.edge=T, edge.width = 3)
#axisPhylo()
axisPhylo(1,root.time=0,  backward=FALSE)
dev.off()

#########################################################################

stri = "PD26414"
a=read.table("~/Downloads/PD26414_DP_redo.txt", sep="\t",header=T)
#a$code_sc = "PD26414"
mm=melt(a[,c(1:9)])
head(a)
head(mm)

pdf(sprintf("~/Desktop/%s_bargg.pdf",stri), useDingbats=FALSE)
p = ggplot(mm,aes(code,value,fill=variable))+geom_bar(stat="identity")+theme_minimal()+scale_fill_manual(values = col)
print(p)
dev.off()

a = "((PD26414.2:9230,PD26414.10:1)PD26414.1:643)PP1:1;"
tr <- read.tree(text = a)

pdf(sprintf("~/Desktop/%s_tree.pdf",stri), useDingbats=FALSE,width=10,height=5)
plot(tr,show.node.label = T, root.edge=T, edge.width = 3)
axisPhylo(1,root.time=0,  backward=FALSE)
dev.off()


############################# PD26411 ############################
stri="PD26411"
m = melter(stri)
sig[grep(stri,sig$sampleID),]
tree2[grep(stri,tree2$patient),c("tree_structure")]


#a = "((PD26403.3:282,(PD26403.2:1248,PD26403.10:1)PD26403.4:173,PD26403.1:3105)PD26403.7:5678)PPP:3;"
#a = "((PD26403.3:282,(PD26403.2:1248,PD26403.10:1)PD26403.4:173,PD26403.1:3105)PD26403.7:5678)PPP:3;"

a = "((PD26411.10:1,(PD26411.3:470,PD26411.11:1)PD26411.4:3530)PD26411.5:4410)PPP:1;"
tr <- read.tree(text = a)

#pdf(sprintf("~/Desktop/%s_tree.pdf",stri), useDingbats=FALSE,width=12,height=4)
plot(tr,show.node.label = T, root.edge=T, edge.width = 3, las=2)
axisPhylo(1,root.time=0,  backward=FALSE)
#dev.off()

#########################################################################



############################# PD26412 ############################
stri="PD26412"
m = melter(stri)
sig[grep(stri,sig$sampleID),]
tree2[grep(stri,tree2$patient),c("tree_structure")]


#a = "((PD26403.3:282,(PD26403.2:1248,PD26403.10:1)PD26403.4:173,PD26403.1:3105)PD26403.7:5678)PPP:3;"
#a = "((PD26403.3:282,(PD26403.2:1248,PD26403.10:1)PD26403.4:173,PD26403.1:3105)PD26403.7:5678)PPP:3;"

a = "((PD26411.10:1,(PD26411.3:470,PD26411.11:1)PD26411.4:3530)PD26411.5:4410)PPP:1;"
a = "(((PD26412.5:2187,PD26412.10:1)PD26412.7:2058,(PD26412.2:225,PD26412.1:154)PD26412.4:1820)PD26412.6:3221)PP1:1;"
tr <- read.tree(text = a)

pdf(sprintf("~/Desktop/%s_tree.pdf",stri), useDingbats=FALSE,width=10,height=5)
plot(tr,show.node.label = T, root.edge=T, edge.width = 3, las=2)
axisPhylo(1,root.time=0,  backward=FALSE)
dev.off()

#########################################################################



############################# PD26423 ############################
stri="PD26423"
m = melter(stri)
sig[grep(stri,sig$sampleID),]
tree2[grep(stri,tree2$patient),c("tree_structure")]


a = "(((PD26423.6:1128,PD26423.10:1)PD26423.7:2903,PD26423.3:347,PD26423.5:420,PD26423.2:212)PD26423.8:3567)PP1:1;"

tr <- read.tree(text = a)

pdf(sprintf("~/Desktop/%s_tree.pdf",stri), useDingbats=FALSE,width=10,height=5)
plot(tr,show.node.label = T, root.edge=T, edge.width = 3, las=2)
axisPhylo(1,root.time=0,  backward=FALSE)
dev.off()

#########################################################################



############################# PD26420 ############################
stri="PD26420"
m = melter(stri)
sig[grep(stri,sig$sampleID),]
tree2[grep(stri,tree2$patient),c("tree_structure")]



a = "((PD26420.2:183,PD26420.10:1)PD26420.5:3141)PP1:1;"
tr <- read.tree(text = a)

pdf(sprintf("~/Desktop/%s_tree.pdf",stri), useDingbats=FALSE,width=10,height=5)
plot(tr,show.node.label = T, root.edge=T, edge.width = 3, las=2)
axisPhylo(1,root.time=0,  backward=FALSE)
dev.off()

#########################################################################






