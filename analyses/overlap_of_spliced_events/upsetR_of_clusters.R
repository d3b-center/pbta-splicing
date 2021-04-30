################################################################################
# upsetR_of_clusters.R
# script that takes in "results/CC_memberships_expr" and CC_memberships_psi data 
#files and computes cluster membership/groupings overlap
# written by Ammar Naqvi
#
# usage: Rscript upsetR_of_clusters.R
################################################################################

library(UpSetR)
library(ggplot2)
library(gplots)

## upsetR plots for recurrent LSVs
#cc_expr_1 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_expr.Cl1.v2.txt", sep = "\t")
#cc_expr_2 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_expr.Cl2.v2.txt", sep = "\t")
#cc_expr_3 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_expr.Cl3.v2.txt", sep = "\t")
#cc_expr_4 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_expr.Cl4.v2.txt", sep = "\t")

cc_expr_1 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_tpm.Cl1.txt", sep = "\t")
cc_expr_2 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_tpm.Cl2.txt", sep = "\t")
cc_expr_3 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_tpm.Cl3.txt", sep = "\t")
cc_expr_4 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_tpm.Cl4.txt", sep = "\t")
cc_expr_5 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_tpm.Cl5.txt", sep = "\t")
cc_expr_6 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_tpm.Cl6.txt", sep = "\t")

cc_psi_1 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_psi.Cl1.txt", sep = "\t")
cc_psi_2 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_psi.Cl2.txt", sep = "\t")
cc_psi_3 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_psi.Cl3.txt", sep = "\t")
cc_psi_4 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_psi.Cl4.txt", sep = "\t")
cc_psi_5 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_psi.Cl5.txt", sep = "\t")
cc_psi_6 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/CC_memberships_psi.Cl6.txt", sep = "\t")


listInput <- list("expr_1" =cc_expr_1$V1, 
                  "expr_2" =cc_expr_2$V1,
                  "expr_3" =cc_expr_3$V1,
                  "expr_4" =cc_expr_4$V1,
                  "expr_5" =cc_expr_5$V1,
                  "expr_6" =cc_expr_6$V1,
                  "psi_1" =cc_psi_1$V1, 
                  "psi_2" =cc_psi_2$V1,
                  "psi_3" =cc_psi_3$V1,
                  "psi_4" =cc_psi_4$V1,
                  "psi_5" =cc_psi_5$V1,
                  "psi_6" =cc_psi_6$V1)
                  
k = upset(fromList(listInput), 
          mainbar.y.label = "Cluster Intersections", sets.x.label = "Clusters", order.by = "freq", line.size = 1.5, nsets = 12,   
          queries = list(list(query = intersects, params = list("psi_4"), color = "red", active = T), 
                         list(query = intersects, params = list("psi_5"), color = "red", active = T), 
                         list(query = intersects, params = list("psi_6"), color = "red", active = T),
                         list(query = intersects, params = list("psi_6","expr_6"), color = "red", active = T),
                         list(query = intersects, params = list("psi_4","expr_4"), color = "red", active = T),
                         list(query = intersects, params = list("psi_4","expr_2"), color = "red", active = T),
                         list(query = intersects, params = list("psi_4","expr_1"), color = "red", active = T),
                         list(query = intersects, params = list("psi_5","expr_4"), color = "red", active = T)
                         ))

          #mb.ratio = c(0.5,0.50), text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.4),point.size = 1.5, line.size = 1.5, nsets = 12,
        

# ##                            params = list("psi_4", "expr_1"), color = "red", active = T), list(query = intersects,
# 
# params = list("psi_4", "expr_2"), color = "red", active = T), list(query = intersects,
#                                                                    params = list("psi_4", "expr_4"), color = "red", active = T), list(query = intersects,
#                                                                                                                                       params = list("psi_5", "expr_1"), color = "red", active = T), list(query = intersects,                                                                   
#                                                                                                                                                                                                          params = list("psi_6", "expr_4"), color = "red", active = T), list(query = intersects,

list_1 = read.table("/Users/naqvia/Desktop/polyA_vs_stranded/BS_68KX6A42_BS_D7XRFE0R.polyA_vs_stranded.genes.txt", sep = "\t")
list_2 = read.table("/Users/naqvia/Desktop/polyA_vs_stranded/BS_7WM3MNZ0_BS_KABQQA0T.polyA_vs_stranded.genes.txt", sep = "\t")
list_3 = read.table("/Users/naqvia/Desktop/polyA_vs_stranded/BS_BYCX6VK1_BS_SB12W1XT.polyA_vs_stranded.genes.txt", sep = "\t")
list_4 = read.table("/Users/naqvia/Desktop/polyA_vs_stranded/BS_HWGWYCY7_BS_HE0WJRW6.polyA_vs_stranded.genes.txt", sep = "\t")
list_5 = read.table("/Users/naqvia/Desktop/polyA_vs_stranded/BS_QKT3TJVK_BS_8QB4S4VA.polyA_vs_stranded.genes.txt", sep = "\t")
list_6 = read.table("/Users/naqvia/Desktop/polyA_vs_stranded/BS_W4H1D4Y6_BS_FN07P04C.polyA_vs_stranded.genes.txt", sep = "\t")
list_7 = read.table("/Users/naqvia/Desktop/polyA_vs_stranded/BS_X0XXN9BK_BS_SHJA4MR0.polyA_vs_stranded.genes.txt", sep = "\t")

listInput <- list("list_1" =list_1$V1, 
                  "list_2" =list_2$V1,
                  "list_3" =list_3$V1,
                  "list_4" =list_4$V1,
                  "list_5" =list_5$V1, 
                  "list_6" =list_6$V1,
                  "list_7" =list_7$V1)
                       
upset(fromList(listInput), 
      mainbar.y.label = "Intersections", sets.x.label = "Events", order.by = "freq",,point.size = 1.5, line.size = 1.5, nsets = 7)


perc_hist_as.Cran.tsv
perc_hist_as.Epend.tsv
perc_hist_as.Gang.tsv
perc_hist_as.HGAT.tsv
perc_hist_as.LGAT.tsv
perc_hist_as.medul.tsv



list_1 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/perc_hist_as.ATRT.20thr.tsv", sep = "\t")
list_2 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/perc_hist_as.Cran.20thr.tsv", sep = "\t")
list_3 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/perc_hist_as.Epend.20thr.tsv", sep = "\t")
list_4 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/perc_hist_as.HGAT.20thr.tsv", sep = "\t")
list_5 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/perc_hist_as.LGAT.thr20.tsv", sep = "\t")
list_6 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/perc_hist_as.medul.thr20.tsv", sep = "\t")
list_7 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/perc_hist_as.Gang.thr20.tsv",sep = "\t")

listInput <- list("ATRT" =list_1$V2, 
                  "Cran" =list_2$V2,
                  "Epend" =list_3$V2,
                  "HGAT" =list_4$V2,
                  "LGAT" =list_5$V2, 
                  "Med" =list_6$V2,
                  "Gang" =list_7$V2)



upset(fromList(listInput), 
      mainbar.y.label = "", sets.x.label = "Clusters", order.by = "freq",
      mb.ratio = c(0.5,0.50), text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.4),point.size = 2, line.size = 1.5, nsets = 7)

#DMG specific
list_1 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/perc_hist_as.DMG_h3k28.thr20.tsv", sep = "\t")
list_2 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/perc_hist_as.DMG_h3k28_tp53loss.thr20.tsv", sep = "\t")
list_3 = read.table("/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/perc_hist_as.DMG_h3k28_tp53gain.thr20.tsv", sep = "\t")

listInput <- list("H3K28" =list_1$V2, 
                  "H3K28, TP53 Loss" =list_2$V2,
                  "H3K28, TP53 Activated" =list_3$V2)

upset(fromList(listInput), 
      mainbar.y.label = "", sets.x.label = "Subtype", order.by = "freq",
      #mb.ratio = c(0.5,0.50), text.scale = c(1.3, 1.3, 2, 1.3, 2, 2),
      point.size = 3, line.size =1, nsets = 3)

