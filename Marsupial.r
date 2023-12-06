library(readxl)
library(phyloseq)
library(ape)
library(plyr)
library(vegan)
library(rbiom)
library(tidyverse)
library(ggpubr)
library(rstatix)
install.packages("FSA")
library(FSA)
options(max.print=999999)


otu_matMar <-read_excel("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/PhyseqMarsup.xlsx", sheet="OTU")
tax_matMar <- read_excel("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/PhyseqMarsup.xlsx", sheet="taxon")
MetaMar <-read_excel("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/PhyseqMarsup.xlsx", sheet="Samples")
otu_matMar <- otu_matMar %>%
  tibble::column_to_rownames("#OTU ID")
tax_matMar <- tax_matMar %>%
  tibble::column_to_rownames("#OTU ID")
MetaMar <- MetaMar  %>%
  tibble::column_to_rownames("Sample")
otu_matMar <- as.matrix(otu_matMar)
tax_matMar <- as.matrix(tax_matMar)
OTUMar = otu_table(otu_matMar, taxa_are_rows = TRUE)
TAXMar = tax_table(tax_matMar)
samplesMar = sample_data(MetaMar)
TreeMar <-ape::read.tree(file="~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Genera.nwk")
PhyseqMar <-phyloseq(OTUMar, TAXMar, samplesMar, TreeMar)
PhyseqMar
##Fig. 3A, and FigS1A
AlphaMar <-estimate_richness(PhyseqMar, measures=c("Shannon", "Simpson", "InvSimpson"))
AlphaMarComb <-cbind(samplesMar, AlphaMar)
sample_data(PhyseqMar) <- sample_data(PhyseqMar) %>%
  reorder_levels(Animal, order = c("Eastern_Grey_Kangaroo", "Red_Kangaroo", "Red-necked_Wallaby", "Koala", "Common_Wombat", "Southern_Hairy_Nosed_Wombat"))

p2= plot_richness(PhyseqMar, x="Family", measures=c("Shannon", "Simpson", "InvSimpson"))
p2 + geom_boxplot(data =p2$data, alpha = 0.1)+
  theme(axis.text.x = element_text(angle=90, hjust=1))

p3= plot_richness(PhyseqMar, x="GutType", measures=c("Shannon", "Simpson", "InvSimpson"))
p3 + geom_boxplot(data =p3$data, alpha = 0.1)+
  theme(axis.text.x = element_text(angle=90, hjust=1))

p4= plot_richness(PhyseqMar, x="Habitat", measures=c("Shannon", "Simpson", "InvSimpson"))
p4 + geom_boxplot(data =p4$data, alpha = 0.1)+
  theme(axis.text.x = element_text(angle=90, hjust=1))

p5= plot_richness(PhyseqMar, x="Nutrition", measures=c("Shannon", "Simpson", "InvSimpson"))
p5 + geom_boxplot(data =p5$data, alpha = 0.1)+
  theme(axis.text.x = element_text(angle=90, hjust=1))

p6= plot_richness(PhyseqMar, x="Animal", measures=c("Shannon", "Simpson", "InvSimpson"))
p6 + geom_boxplot(data =p6$data, alpha = 0.1)+
  theme(axis.text.x = element_text(angle=90, hjust=1))

otu_matMix <-read_excel("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/PhyseqMix.xlsx", sheet="OTU")
tax_matMix <- read_excel("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/PhyseqMix.xlsx", sheet="taxon")
MetaMix <-read_excel("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/PhyseqMix.xlsx", sheet="Samples")
otu_matMix <- otu_matMix %>%
  tibble::column_to_rownames("#OTU ID")
tax_matMix <- tax_matMix %>%
  tibble::column_to_rownames("#OTU ID")
MetaMix <- MetaMix  %>%
  tibble::column_to_rownames("Sample")
otu_matMix <- as.matrix(otu_matMix)
tax_matMix <- as.matrix(tax_matMix)
OTUMix = otu_table(otu_matMix, taxa_are_rows = FALSE)
TAXMix = tax_table(tax_matMix)
samplesMix = sample_data(MetaMix)
TreeMix <-ape::read.tree(file="~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Genera_mafft_aligned.nwk")
PhyseqMix <-phyloseq(OTUMix, TAXMix, samplesMix, TreeMix)
PhyseqMix


## running KruskalWallis on the alpha diversity results grouped by host factors
## A high F-value and low P-value mean that the variance between groups is significant. 
## The Dunn gives the pairwise p-values for all possible pairs of groups if significant

### on Shannon-Animal
kruskal.test(AlphaMar$Shannon ~ sample_data(PhyseqMar)$Animal, data=AlphaMar)
### on Shannon-Family
kruskal.test(AlphaMar$Shannon ~ sample_data(PhyseqMar)$Family, data=AlphaMar)
### on Shannon-GutType
kruskal.test(AlphaMar$Shannon ~ sample_data(PhyseqMar)$GutType, data=AlphaMar)
### on Shannon-Habitat
kruskal.test(AlphaMar$Shannon ~ sample_data(PhyseqMar)$Habitat, data=AlphaMar)
### on Shannon-Nutr
kruskal.test(AlphaMar$Shannon ~ sample_data(PhyseqMar)$Nutrition, data=AlphaMar)
### on Simpson-Animal
kruskal.test(AlphaMar$Simpson ~ sample_data(PhyseqMar)$Animal, data=AlphaMar)
### on Simpson-Family
kruskal.test(AlphaMar$Simpson ~ sample_data(PhyseqMar)$Family, data=AlphaMar)
### on Simpson-GutType
kruskal.test(AlphaMar$Simpson ~ sample_data(PhyseqMar)$GutType, data=AlphaMar)
### on Simpson-Habitat
kruskal.test(AlphaMar$Simpson ~ sample_data(PhyseqMar)$Habitat, data=AlphaMar)
### on Simpson-Nutr
kruskal.test(AlphaMar$Simpson ~ sample_data(PhyseqMar)$Nutrition, data=AlphaMar)
### on InvSimpson-Animal
kruskal.test(AlphaMar$InvSimpson ~ sample_data(PhyseqMar)$Animal, data=AlphaMar)
### on InvSimpson-Family
kruskal.test(AlphaMar$InvSimpson ~ sample_data(PhyseqMar)$Family, data=AlphaMar)
### on InvSimpson-GutType
kruskal.test(AlphaMar$InvSimpson ~ sample_data(PhyseqMar)$GutType, data=AlphaMar)
### on InvSimpson-Habitat
kruskal.test(AlphaMar$InvSimpson ~ sample_data(PhyseqMar)$Habitat, data=AlphaMar)
### on InvSimpson-Nutr
kruskal.test(AlphaMar$InvSimpson ~ sample_data(PhyseqMar)$Nutrition, data=AlphaMar)

kruskal.test(AlphaMix$Shannon ~ AlphaMix$GutType_Class, data=AlphaMix)
dunnTest(AlphaMix$Shannon ~ AlphaMix$GutType_Class, data=AlphaMix, method="holm")
kruskal.test(AlphaMix$Simpson ~ AlphaMix$GutType_Class, data=AlphaMix)
dunnTest(AlphaMix$Simpson ~ AlphaMix$GutType_Class, data=AlphaMix, method="holm")
kruskal.test(AlphaMix$InvSimpson ~ AlphaMix$GutType_Class, data=AlphaMix)
dunnTest(AlphaMix$InvSimpson ~ AlphaMix$GutType_Class, data=AlphaMix, method="holm")
sample_data(PhyseqMix) <- sample_data(PhyseqMix) %>%
  reorder_levels(Animal, order = c("Eastern_Grey_Kangaroo", "Red_Kangaroo", "Red-necked_Wallaby", "Koala", "Common_Wombat", "Southern_Hairy_Nosed_Wombat", "Cow", "Goat", "Sheep", "Elephant", "Rhino", "Zebra", "Horse"))

p7= plot_richness(PhyseqMix, x="Animal", measures=c("Shannon", "Simpson", "InvSimpson"))
p7 + geom_boxplot(data =p1$data, alpha = 0.1)+
  theme(axis.text.x = element_text(angle=90, hjust=1))




AlphaMix %>%
  group_by(Class)

pwcMix <- AlphaMix %>% 
  group_by(Class) %>%
  wilcox_test(Shannon ~ GutType, p.adjust.method = "bonferroni") 
P=ggboxplot(AlphaMix, x = "Class", y = "Shannon", color="GutType", bxp.errorbar=TRUE)
pwcMix <- pwcMix %>% add_xy_position(x = "Class")
P+stat_pvalue_manual(pwcMix)

pwcMix2 <- AlphaMix %>% 
  group_by(Class) %>%
  wilcox_test(Shannon ~ GutType, p.adjust.method = "bonferroni") 
P2=ggboxplot(AlphaMix, x = "Class", y = "Shannon", color="GutType", bxp.errorbar=TRUE)
pwcMix2 <- pwcMix2 %>% add_xy_position(x = "Class")
P2+stat_pvalue_manual(pwcMix2)

pwcMix3 <- AlphaMix %>% 
  group_by(Class) %>%
  wilcox_test(Simpson ~ GutType, p.adjust.method = "bonferroni") 
P3=ggboxplot(AlphaMix, x = "Class", y = "Simpson", color="GutType", bxp.errorbar=TRUE)
pwcMix3 <- pwcMix3 %>% add_xy_position(x = "Class")
P3+stat_pvalue_manual(pwcMix3)

##FIG5
PCoAMar_UnifracW <-ordinate(PhyseqMar, method="PCoA", distance="unifrac", weighted=TRUE)
PCoAMar3 <-plot_ordination(PhyseqMar, PCoAMar_UnifracW, "samples", color="Family", shape="GutType")
PCoAMar3=PCoAMar3+geom_point(size=2)+scale_shape_manual(values=seq(0,15))
PCoAMar3
write.table(PCoAMar_UnifracW$vectors, file="~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/PCoA_UnifracW")
PCoAMar4 <-plot_ordination(PhyseqMar, PCoAMar_UnifracW, "samples", color="Animal", shape="GutType")
PCoAMar4=PCoAMar4+geom_point(size=2)+scale_shape_manual(values=seq(0,15))
PCoAMar4

PCoAMar5 <-plot_ordination(PhyseqMar, PCoAMar_UnifracW, "samples", color="Habitat", shape="GutType")
PCoAMar5=PCoAMar5+geom_point(size=2)+scale_shape_manual(values=seq(0,15))
PCoAMar5

PCoAMar6 <-plot_ordination(PhyseqMar, PCoAMar_UnifracW, "samples", color="Nutrition", shape="GutType")
PCoAMar6=PCoAMar6+geom_point(size=2)+scale_shape_manual(values=seq(0,15))
PCoAMar6

Mar_UnifracW <-phyloseq::distance(PhyseqMar, method="wunifrac")
Samples <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Samples.txt", header=TRUE)
Adonis2 <-adonis(Mar_UnifracW ~ Animal, AlphaMar)
Adonis2$aov.tab
Adonis3 <-adonis(Mar_UnifracW ~ Family, AlphaMar)
Adonis3$aov.tab
Adonis4 <-adonis(Mar_UnifracW ~ Habitat, Samples)
Adonis4$aov.tab
Adonis5 <-adonis(Mar_UnifracW ~ Nutrition, AlphaMar)
Adonis5$aov.tab
Adonis6 <-adonis(Mar_UnifracW ~ GutType, AlphaMar)
Adonis6$aov.tab



PCoAMix_UnifracW <-ordinate(PhyseqMix, method="PCoA", distance="unifrac", weighted=TRUE)
PCoAMix3 <-plot_ordination(PhyseqMix, PCoAMix_UnifracW, "samples", color="Class", shape="GutType")
PCoAMix3=PCoAMix3+geom_point(size=2)+scale_shape_manual(values=seq(0,15))+stat_ellipse(show.legend=TRUE)
PCoAMix3



FactorsMixed <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/FactorsMixed.txt", header=TRUE)
Mix_UnifracW <-phyloseq::distance(PhyseqMix, method="wunifrac")
Adonis7 <-adonis(Mix_UnifracW ~ GutType, FactorsMixed) 
Adonis7$aov.tab
Adonis8 <-adonis(Mix_UnifracW ~ Class, FactorsMixed) 
Adonis8$aov.tab
Adonis9 <-adonis(Mix_UnifracW ~ GutType*Class, FactorsMixed) 
Adonis9$aov.tab
Adonis10 <-adonis(Mix_UnifracW ~ GutType+Class, FactorsMixed) 
Adonis10$aov.tab


##Fig4

library(NST)

DM <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/comm.txt", header=TRUE, row.names=1)
DM2 <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/comm2.txt", header=TRUE, row.names=1)
DM3 <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/comm3.txt", header=TRUE, row.names=1)
Family <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Family.txt", row.names=1, header=TRUE)
Animal <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Animal.txt", row.names=1, header=TRUE)
Gut <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Gut.txt", row.names=1, header=TRUE)
Habitat <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Habitat.txt", row.names=1, header=TRUE)
Nutr <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Nutr.txt", row.names=1, header=TRUE)


#Jaccard
tnstoutFamily_Jac=NST::tNST(comm=DM, dist.method="jaccard", group=Family,
                        abundance.weighted=FALSE, nworker=32,
                        null.model="PF", output.rand = TRUE,
                        SES = TRUE, RC = TRUE)
write.table(tnstoutFamily_Jac$index.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindex.Fam.jaccard.txt")
write.table(tnstoutFamily_Jac$index.pair, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpair.Fam.jaccard.txt")
write.table(tnstoutFamily_Jac$index.pair.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpairgrp.Fam.jaccard.txt")

tnstoutAnimal_Jac=NST::tNST(comm=DM2, dist.method="jaccard", group=Animal,
                             abundance.weighted=FALSE, nworker=32,
                             null.model="PF", output.rand = TRUE,
                             SES = TRUE, RC = TRUE)
write.table(tnstoutAnimal_Jac$index.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindex.Animal.jaccard.txt")
write.table(tnstoutAnimal_Jac$index.pair, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpair.Animal.jaccard.txt")
write.table(tnstoutAnimal_Jac$index.pair.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpairgrp.Animal.jaccard.txt")

tnst.btFam.Jac=NST::nst.boot(nst.result=tnstoutFamily_Jac, group=Family, nworker=32, out.detail=TRUE)
write.table(tnst.btFam.Jac$summary,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Fam.Jac.txt")
write.table(tnst.btFam.Jac$detail$NST.boot,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Fam.Jac.detail.txt")
write.table(tnst.btFam.Jac$compare,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Fam_Comp.Jac.txt")

tnst.btAnimal.Jac=NST::nst.boot(nst.result=tnstoutAnimal_Jac, group=Animal,nworker=32, out.detail=TRUE)
write.table(tnst.btAnimal.Jac$summary,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Animal.Jac.txt")
write.table(tnst.btAnimal.Jac$detail$NST.boot,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Animal.Jac.detail.txt")
write.table(tnst.btAnimal.Jac$compare,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Animal_Comp.Jac.txt")

#Bray
tnstoutFamily_Bray=NST::tNST(comm=DM, dist.method="bray", group=Family,
                             abundance.weighted=FALSE, nworker=32,
                             null.model="PF", output.rand = TRUE,
                             SES = TRUE, RC = TRUE)
write.table(tnstoutFamily_Bray$index.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindex.Fam.bray.txt")
write.table(tnstoutFamily_Bray$index.pair, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpair.Fam.bray.txt")
write.table(tnstoutFamily_Bray$index.pair.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpairgrp.Fam.bray.txt")

tnstoutAnimal_Bray=NST::tNST(comm=DM2, dist.method="bray", group=Animal,
                             abundance.weighted=FALSE, nworker=32,
                             null.model="PF", output.rand = TRUE,
                             SES = TRUE, RC = TRUE)
write.table(tnstoutAnimal_Bray$index.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindex.Animal.bray.txt")
write.table(tnstoutAnimal_Bray$index.pair, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpair.Animal.bray.txt")
write.table(tnstoutAnimal_Bray$index.pair.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpairgrp.Animal.bray.txt")

tnst.btFam.Bray=NST::nst.boot(nst.result=tnstoutFamily_Bray, group=Family, nworker=32, out.detail=TRUE)
write.table(tnst.btFam.Bray$summary,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Fam.Bray.txt")
write.table(tnst.btFam.Bray$detail$NST.boot,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Fam.Bray.detail.txt")
write.table(tnst.btFam.Bray$compare,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Fam_Comp.Bray.txt")

tnst.btAnimal.Bray=NST::nst.boot(nst.result=tnstoutAnimal_Bray, group=Animal,nworker=32, out.detail=TRUE)
write.table(tnst.btAnimal.Bray$summary,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Animal.Bray.txt")
write.table(tnst.btAnimal.Bray$detail$NST.boot,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Animaly.Bray.detail.txt")
write.table(tnst.btAnimal.Bray$compare,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Animal_Comp.Bray.txt")


#Jaccard
tnstoutGut_Jac=NST::tNST(comm=DM3, dist.method="jaccard", group=Gut,
                         abundance.weighted=FALSE, nworker=32,
                         null.model="PF", output.rand = TRUE,
                         SES = TRUE, RC = TRUE)
write.table(tnstoutGut_Jac$index.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindex.Gut.jaccard.txt")
write.table(tnstoutGut_Jac$index.pair, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpair.Gut.jaccard.txt")
write.table(tnstoutGut_Jac$index.pair.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpairgrp.Gut.jaccard.txt")

tnstoutHabitat_Jac=NST::tNST(comm=DM3, dist.method="jaccard", group=Habitat,
                             abundance.weighted=FALSE, nworker=32,
                             null.model="PF", output.rand = TRUE,
                             SES = TRUE, RC = TRUE)
write.table(tnstoutHabitat_Jac$index.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindex.Habitat.jaccard.txt")
write.table(tnstoutHabitat_Jac$index.pair, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpair.Habitat.jaccard.txt")
write.table(tnstoutHabitat_Jac$index.pair.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpairgrp.Habitat.jaccard.txt")

tnst.btGut.Jac=NST::nst.boot(nst.result=tnstoutGut_Jac, group=Gut, nworker=32, out.detail=TRUE)
write.table(tnst.btGut.Jac$summary,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Gut.Jac.txt")
write.table(tnst.btGut.Jac$detail$NST.boot,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Gut.Jac.detail.txt")
write.table(tnst.btGut.Jac$compare,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Gut_Comp.Jac.txt")

tnst.btHabitat.Jac=NST::nst.boot(nst.result=tnstoutHabitat_Jac, group=Habitat,nworker=32, out.detail=TRUE)
write.table(tnst.btHabitat.Jac$summary,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Habitat.Jac.txt")
write.table(tnst.btHabitat.Jac$detail$NST.boot,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Habitaty.Jac.detail.txt")
write.table(tnst.btHabitat.Jac$compare,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Habitat_Comp.Jac.txt")
#Bray
tnstoutGut_Bray=NST::tNST(comm=DM3, dist.method="bray", group=Gut,
                          abundance.weighted=FALSE, nworker=32,
                          null.model="PF", output.rand = TRUE,
                          SES = TRUE, RC = TRUE)
write.table(tnstoutGut_Bray$index.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindex.Gut.bray.txt")
write.table(tnstoutGut_Bray$index.pair, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpair.Gut.bray.txt")
write.table(tnstoutGut_Bray$index.pair.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpairgrp.Gut.bray.txt")

tnstoutHabitat_Bray=NST::tNST(comm=DM3, dist.method="bray", group=Habitat,
                              abundance.weighted=FALSE, nworker=32,
                              null.model="PF", output.rand = TRUE,
                              SES = TRUE, RC = TRUE)
write.table(tnstoutHabitat_Bray$index.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindex.Habitat.bray.txt")
write.table(tnstoutHabitat_Bray$index.pair, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpair.Habitat.bray.txt")
write.table(tnstoutHabitat_Bray$index.pair.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpairgrp.Habitat.bray.txt")

tnst.btGut.Bray=NST::nst.boot(nst.result=tnstoutGut_Bray, group=Gut, nworker=32, out.detail=TRUE)
write.table(tnst.btGut.Bray$summary,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Gut.Bray.txt")
write.table(tnst.btGut.Bray$detail$NST.boot,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Gut.Bray.detail.txt")
write.table(tnst.btGut.Bray$compare,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Gut_Comp.Bray.txt")

tnst.btHabitat.Bray=NST::nst.boot(nst.result=tnstoutHabitat_Bray, group=Habitat,nworker=32, out.detail=TRUE)
write.table(tnst.btHabitat.Bray$summary,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Habitat.Bray.txt")
write.table(tnst.btHabitat.Bray$detail$NST.boot,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Habitaty.Bray.detail.txt")
write.table(tnst.btHabitat.Bray$compare,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Habitat_Comp.Bray.txt")
#Jaccard
tnstoutNutr_Jac=NST::tNST(comm=DM3, dist.method="jaccard", group=Nutr,
                          abundance.weighted=FALSE, nworker=32,
                          null.model="PF", output.rand = TRUE,
                          SES = TRUE, RC = TRUE)
write.table(tnstoutNutr_Jac$index.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindex.Nutr.jaccard.txt")
write.table(tnstoutNutr_Jac$index.pair, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpair.Nutr.jaccard.txt")
write.table(tnstoutNutr_Jac$index.pair.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpairgrp.Nutr.jaccard.txt")
tnst.btNutr.Jac=NST::nst.boot(nst.result=tnstoutNutr_Jac, group=Nutr, nworker=32, out.detail=TRUE)
write.table(tnst.btNutr.Jac$summary,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Nutr.Jac.txt")
write.table(tnst.btNutr.Jac$detail$NST.boot,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Nutr.Jac.detail.txt")
write.table(tnst.btNutr.Jac$compare,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Nutr_Comp.Jac.txt")

#Bray
tnstoutNutr_Bray=NST::tNST(comm=DM3, dist.method="bray", group=Nutr,
                           abundance.weighted=FALSE, nworker=32,
                           null.model="PF", output.rand = TRUE,
                           SES = TRUE, RC = TRUE)
write.table(tnstoutNutr_Bray$index.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindex.Nutr.bray.txt")
write.table(tnstoutNutr_Bray$index.pair, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpair.Nutr.bray.txt")
write.table(tnstoutNutr_Bray$index.pair.grp, file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstindexpairgrp.Nutr.bray.txt")
tnst.btNutr.Bray=NST::nst.boot(nst.result=tnstoutNutr_Bray, group=Nutr, nworker=32, out.detail=TRUE)
write.table(tnst.btNutr.Bray$summary,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Nutr.Bray.txt")
write.table(tnst.btNutr.Bray$detail$NST.boot,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Nutr.Bray.detail.txt")
write.table(tnst.btNutr.Bray$compare,file = "~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/tnstboot.Nutr_Comp.Bray.txt")


##NST_Boxplots
Fam_NST <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Family_NST.txt", header=TRUE)
Anim_NST <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Animal_NST.txt", header=TRUE)
Gut_NST <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Gut_NST.txt", header=TRUE)
Habitat_NST <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Habitat_NST.txt", header=TRUE)
Nutr_NST <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Nutr_NST.txt", header=TRUE)

ggboxplot(Fam_NST, x = "Family", y = "NST_Bray", bxp.errorbar=TRUE)
ggboxplot(Fam_NST, x = "Family", y = "NST_Jac", bxp.errorbar=TRUE)

Anim_NST <- Anim_NST %>%
  reorder_levels(Animal, order = c("Eastern_Grey_Kangaroo", "Red_Kangaroo", "Red.necked_Wallaby", "Koala", "Common_Wombat", "Southern_Hairy_Nosed_Wombat"))

ggboxplot(Anim_NST, x = "Animal", y = "NST_Bray", bxp.errorbar=TRUE)
ggboxplot(Anim_NST, x = "Animal", y = "NST_Jac", bxp.errorbar=TRUE)

ggboxplot(Gut_NST, x = "Gut", y = "NST_Bray", bxp.errorbar=TRUE)
ggboxplot(Gut_NST, x = "Gut", y = "NST_Jac", bxp.errorbar=TRUE)

ggboxplot(Nutr_NST, x = "Nutrition", y = "NST_Bray", bxp.errorbar=TRUE)
ggboxplot(Nutr_NST, x = "Nutrition", y = "NST_Jac", bxp.errorbar=TRUE)

ggboxplot(Habitat_NST, x = "Habitat", y = "NST_Bray", bxp.errorbar=TRUE)
ggboxplot(Habitat_NST, x = "Habitat", y = "NST_Jac", bxp.errorbar=TRUE)


pwcNSTBray1 <- Fam_NST %>% 
  wilcox_test(NST_Bray ~ Family, p.adjust.method = "bonferroni") 
P10=ggboxplot(Fam_NST, x = "Family", y = "NST_Bray", bxp.errorbar=TRUE)
pwcNSTBray1 <- pwcNSTBray1 %>% add_xy_position(x = "Family")
P10+stat_pvalue_manual(pwcNSTBray1)

pwcNSTJac1 <- Fam_NST %>% 
  wilcox_test(NST_Jac ~ Family, p.adjust.method = "bonferroni") 
P11=ggboxplot(Fam_NST, x = "Family", y = "NST_Jac", bxp.errorbar=TRUE)
pwcNSTJac1 <- pwcNSTJac1 %>% add_xy_position(x = "Family")
P11+stat_pvalue_manual(pwcNSTJac1)

pwcNSTBray2 <- Anim_NST %>% 
  wilcox_test(NST_Bray ~ Animal, p.adjust.method = "bonferroni") 
P12=ggboxplot(Anim_NST, x = "Animal", y = "NST_Bray", bxp.errorbar=TRUE)
pwcNSTBray2 <- pwcNSTBray2 %>% add_xy_position(x = "Animal")
P12+stat_pvalue_manual(pwcNSTBray2)

pwcNSTJac2 <- Anim_NST %>% 
  wilcox_test(NST_Jac ~ Animal, p.adjust.method = "bonferroni") 
P13=ggboxplot(Anim_NST, x = "Animal", y = "NST_Jac", bxp.errorbar=TRUE)
pwcNSTJac2 <- pwcNSTJac2 %>% add_xy_position(x = "Animal")
P13+stat_pvalue_manual(pwcNSTJac2)


pwcNSTBray3 <- Gut_NST %>% 
  wilcox_test(NST_Bray ~ Gut, p.adjust.method = "bonferroni") 
P14=ggboxplot(Gut_NST, x = "Gut", y = "NST_Bray", bxp.errorbar=TRUE)
pwcNSTBray3 <- pwcNSTBray3 %>% add_xy_position(x = "Gut")
P14+stat_pvalue_manual(pwcNSTBray3)

pwcNSTJac3 <- Gut_NST %>% 
  wilcox_test(NST_Jac ~ Gut, p.adjust.method = "bonferroni") 
P15=ggboxplot(Gut_NST, x = "Gut", y = "NST_Jac", bxp.errorbar=TRUE)
pwcNSTJac3 <- pwcNSTJac3 %>% add_xy_position(x = "Gut")
P15+stat_pvalue_manual(pwcNSTJac3)


pwcNSTBray4 <- Habitat_NST %>% 
  wilcox_test(NST_Bray ~ Habitat, p.adjust.method = "bonferroni") 
P16=ggboxplot(Habitat_NST, x = "Habitat", y = "NST_Bray", bxp.errorbar=TRUE)
pwcNSTBray4 <- pwcNSTBray4 %>% add_xy_position(x = "Habitat")
P16+stat_pvalue_manual(pwcNSTBray4)

pwcNSTJac4 <- Habitat_NST %>% 
  wilcox_test(NST_Jac ~ Habitat, p.adjust.method = "bonferroni") 
P17=ggboxplot(Habitat_NST, x = "Habitat", y = "NST_Jac", bxp.errorbar=TRUE)
pwcNSTJac4 <- pwcNSTJac4 %>% add_xy_position(x = "Habitat")
P17+stat_pvalue_manual(pwcNSTJac4)

pwcNSTBray5 <- Nutr_NST %>% 
  wilcox_test(NST_Bray ~ Nutrition, p.adjust.method = "bonferroni") 
P18=ggboxplot(Nutr_NST, x = "Nutrition", y = "NST_Bray", bxp.errorbar=TRUE)
pwcNSTBray5 <- pwcNSTBray5 %>% add_xy_position(x = "Nutrition")
P18+stat_pvalue_manual(pwcNSTBray5)

pwcNSTJac5 <- Nutr_NST %>% 
  wilcox_test(NST_Jac ~ Nutrition, p.adjust.method = "bonferroni") 
P19=ggboxplot(Nutr_NST, x = "Nutrition", y = "NST_Jac", bxp.errorbar=TRUE)
pwcNSTJac5 <- pwcNSTJac5 %>% add_xy_position(x = "Nutrition")
P19+stat_pvalue_manual(pwcNSTJac5)


##Fig 6
library(ggpubr)
qPCR <- qPCR %>%
  reorder_levels(Family, order = c("Macropodidae", "Phascolarctidae", "Vombatidae"))
qPCR3 <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/Paper/qPCR3.txt", header=TRUE)
qPCR1 <-read.table("~/Desktop/GradStudents/Casey-Adrienne/Feces_Nonmammal/MarsupialsPaper/qPCR1.txt", header=TRUE)
qPCR1 <- qPCR1 %>%
  reorder_levels(Animal, order = c("Red_Kangaroo", "Eastern_Grey_Kangaroo", "Red-legged_Pademelon", "Red-necked_Wallaby", "Koala", "Southern_Hairy_Nosed_Wombat", "Common_Wombat"))


qPCR3 <- qPCR3 %>%
  reorder_levels(Animal, order = c("Red_Kangaroo", "Eastern_Grey_Kangaroo", "Red-legged_Pademelon", "Red-necked_Wallaby", "Koala", "Southern_Hairy_Nosed_Wombat", "Common_Wombat", "Cow", "Goat", "Sheep", "Horse"))
kruskal.test(Fungi ~ Animal, data=qPCR1)
dunnTest(Fungi ~ Animal, data=qPCR1, method="holm")
kruskal.test(Fungi ~ Family, data=qPCR1)
dunnTest(Fungi ~ Family, data=qPCR1, method="holm")
wilcox_test(Fungi ~ GutType, data=qPCR1)
kruskal.test(Fungi ~ GutType, data=qPCR1)
P1=ggboxplot(qPCR1, x = "Family", y = "Fungi", color="Family", bxp.errorbar=TRUE)
P1+yscale("log10", .format = TRUE)
P2=ggboxplot(qPCR1, x = "Animal", y = "Fungi", color="Animal", bxp.errorbar=TRUE)
P2+yscale("log10", .format = TRUE)
P3=ggboxplot(qPCR1, x = "GutType", y = "Fungi", color="GutType", bxp.errorbar=TRUE)
P3+yscale("log10", .format = TRUE)

P4=ggboxplot(qPCR3, x = "Animal", y = "Fungi", color="Class", bxp.errorbar=TRUE)
P4+yscale("log10", .format = TRUE)


P5=ggboxplot(qPCR3, x = "Class", y = "Fungi", color="Class", bxp.errorbar=TRUE)
P5+yscale("log10", .format = TRUE)


qPCR3 %>%
  group_by(Class) %>%
  get_summary_stats(Fungi, type = "mean_sd")
pwcQPCR3 <- qPCR3 %>%
  wilcox_test(Fungi ~ Class, p.adjust.method = "bonferroni")







