library(tidyverse)
library(data.table)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)

#setwd("~/Dropbox (HMS)/__research-projects/Stenotrophomonas_maltophilia/__manuscript/revision_1/resistance-MCA")

# Classic palette BuPu, with 4 colors
 coul = brewer.pal(4, "PuOr") 

# I can add more tones to this palette :
 coul = colorRampPalette(coul)(25)

 library(viridis)
  A <- viridis_pal(option = "A")(25)
  B <- viridis_pal(option = "B")(25)
  C <- viridis_pal(option = "C")(25)
  D <- viridis_pal(option = "D")(25)
  E <- viridis_pal(option = "E")(25)


###MULTIPLE CORRESPONDENCE ANALYSIS OF THE SELECTED GENES###

data_mca <- read_rds("data_mca_cleaned.rds")

data_mca <- data_mca[,c(1:3,5:13,15,16,18,19)]
data_mca <- data_mca[,c(1:3,5:10,13,15,16)]
data_mca <- rename(data_mca, Lineage = group)

data_mca <- data_mca %>% 
  mutate_if(sapply(data_mca, is.factor), as.character)%>% 
  na.omit()

data_mca <- data_mca %>% 
  mutate_if(sapply(data_mca, is.numeric), as.factor)%>% 
  na.omit()

str(data_mca)

res_mca <- MCA(data_mca, ncp=10, quali.sup =c(1,2,3),  graph = FALSE)
res_mca_2 <- MCA(data_mca, ncp=10, quali.sup =c(1,2),  graph = FALSE)
eig_val <- get_eigenvalue(res_mca)
eig_val_2 <- get_eigenvalue(res_mca_2)

fviz_screeplot(res_mca, addlabels = TRUE, ylim = c(0, 30), title = "MCA variance by dimension")

var <- get_mca_var(res_mca)

fviz_mca_var(res_mca, choice = "var", col.var ="firebrick", 
             labelsize = 5,
             col.quanti.sup ="dodgerblue4", 
             repel = TRUE, # Avoid text overlapping (slow)
             ggtheme = theme_minimal())

fviz_mca_var(res_mca, 
             repel = TRUE, # Avoid text overlapping (slow)
             ggtheme = theme_minimal())

fviz_mca_var(res_mca, col.var ="firebrick", 
             col.quali.sup ="dodgerblue4",
             repel = TRUE,
             ggtheme = theme_minimal())

summary(res_mca)
dimdesc(res_mca)

plotellipses(res_mca, keepvar="quali", palette = A)

colors_stupid <- rep(c("firebrick","dodgerblue4"),17)
fviz_ellipses(res_mca,4:12, keepvar="quali", palette = colors_stupid,geom ="point",
              ggtheme=theme_gray())


fviz_mca_var(res_mca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
             ggtheme = theme_minimal(),select.var = list(contrib =30))

# Contributions of rows to dimension 1
fviz_contrib(res_mca, choice = "var", axes = 1, top = 15)
# Contributions of rows to dimension 2
fviz_contrib(res_mca, choice = "var", axes = 2, top = 15)

fviz_mca_var(res_mca, col.var = "contrib",invisible="quali.sup",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # avoid text overlapping (slow)
             ggtheme = theme_minimal())

ind <- get_mca_ind(res_mca)

fviz_mca_biplot(res_mca, 
                label="var",
                labelsize = 5,
                habillage = "Lineage",
                palette = C,
                col.var ="firebrick", 
                col.quali.sup ="dodgerblue4",
                select.var = list(cos2 = 10),
                shape.ind =19,
                geom.ind = "text", 
                addEllipses = TRUE, 
                ellipse.type = "confidence",
                ellipse.alpha=0.25,
                ellipse.level=0.99,
                invisible = "quali.sup",
                xlim = c(-1,3),
                ylim = c(-1,3))


# Cos2 of individuals
fviz_cos2(res_mca_groups, choice = "ind", axes = 1:2, top = 20)
# Contribution of individuals to the dimensions
fviz_contrib(res_mca_groups, choice = "ind", axes = 1:2, top = 20)

cl <- colors()
fviz_mca_ind(res_mca,
             axes = c(1,5),
             label = "none", # hide individual labels
             habillage = "origin_detailed", # color by groups 
             #palette = cl[c(20,40)],
             #palette = cl[c(5,10,30,35,40,45,50,55,60,65,70,75,80,85,90,
             #95,100,105,110,115,120,125,130,135)], #for groups
             palette = cl[c(5,10,30,35,40,45,50,55,60)],
             addEllipses = TRUE, ellipse.type = "confidence",
             ggtheme = theme_minimal()) 

fviz_mca_biplot(res_mca,
                habillage = "origin_detailed",
                label = "var",
                palette = cl[c(5,10,30,35,40,45,50,55,60)],
                col.var ="firebrick", 
                col.quali.sup ="dodgerblue4",
                select.var = list(cos2 = 10),
                shape.ind =19,
                geom.ind = "text", 
                addEllipses = TRUE, 
                ellipse.type = "confidence",
                ellipse.alpha=0.25,
                ellipse.level=0.99
)


res_desc <- dimdesc(res_mca, axes = c(1,2))


data_mca_origin <- data_mca


res_mca_origin_filtered <- MCA(data_mca_origin, ncp=10, quali.sup =c(1:3),  graph = FALSE)

fviz_mca_ind(res_mca,
             axes = c(1,2),
             label = c("var","quanti.sup"), # hide individual labels
             habillage = "Lineage", # color by groups 
             #palette = cl[c(20,40,60,80)],
             palette = cl[c(5,10,30,35,40,45,50,55,60,65,70,75,80,85,90,
                            95,100,105,110,115,120,125,130,135)], #for groups
             #palette = cl[c(5,10,30,35,40,45,50,55,60)],
             addEllipses = TRUE, ellipse.type = "confidence",ellipse.alpha=0.3,
             ellipse.level=0.995,
             ggtheme = theme_minimal()) 



eig.val <- res_mca$eig
barplot(eig.val[, 2], 
        names.arg = 1:nrow(eig.val), 
        main = "Variances Explained by Dimensions (%)",
        xlab = "Principal Dimensions",
        ylab = "Percentage of variances",
        col ="steelblue")
######
###just for environment plots###
######

F <- viridis_pal(option = "E")(3)
G <- viridis_pal(option = "E")(5)

data_mca_first <- data_mca

data_mca_first$origin_detailed[data_mca_first$origin_detailed == "respiratory"] <- "human-respiratory"
data_mca_first$origin_detailed[data_mca_first$origin_detailed == "sputum"] <- "human-respiratory"
data_mca_first$origin_detailed[data_mca_first$origin_detailed == "unkown"] <- "unknown"

data_mca_first <- data_mca_first %>% 
  filter(origin_detailed != "unknown") %>% 
  rename(Isolation_source =origin_detailed)

res_mca_first <- MCA(data_mca_first, ncp=10,quali.sup =c(1,2,3),  graph = FALSE)



fviz_mca_biplot(res_mca_first,
                habillage = "Isolation_source" ,
                labelsize=5,
                label = "var",
                palette = cl[c(417,142,121,491,115)],
                col.var ="firebrick", 
                col.quali.sup ="dodgerblue4",
                select.var = list(cos2 = 10),
                shape.ind =19,
                geom.ind = "text", 
                addEllipses = TRUE, 
                invisible = "quali.sup",
                ellipse.type = "confidence",
                ellipse.alpha=0.50,
                ellipse.level=0.99
)

###Fisher and proportion tests###
data_mca_second <- data_mca
data_mca_second$origin_detailed[data_mca$origin_detailed == "respiratory"] <- "human-respiratory"
data_mca_second$origin_detailed[data_mca$origin_detailed == "sputum"] <- "human-respiratory"
data_mca_second$origin_detailed[data_mca$origin_detailed == "unkown"] <- "unknown"

data_mca_second <- data_mca_second %>% 
  filter(origin_detailed != "unknown")

###Fisher tests for origin vs groups###

fisher_table <- table(data_mca_second$origin_detailed,data_mca_second$Lineage)


f_p_test <- function(x) {
  
  df <- data.frame(matrix(nrow = length(row.names(x)), 
                          ncol = length(colnames(x))))
  
  rownames(df) <- rownames(x)
  colnames(df) <- colnames(x)
  
  for (i in 1:length(rownames(x))){
    for (j in 1:length(colnames(x))){
      s <- x[i,j]
      m <- matrix(c(s, sum(x[i,]) - s,
                    sum(x[,j]) - s, s + sum(x) -
                      sum(x[i,]) - sum(x[,j])),
                  nrow = 2, byrow = T)
      if (any(m<5)){
        if (sum(x[i,j])/sum(x[i,]) > 
            sum(x[,j])/sum(x)) {
          df[i,j] <- fisher.test(m,alternative = "greater")$p.value  
        }else{
          df[i,j] <- fisher.test(m,alternative = "less")$p.value
        }
      }else{
        if (sum(x[i,j])/sum(x[i,]) > 
            sum(x[,j])/sum(x)) {
          df[i,j] <- prop.test(m,alternative = "greater")$p.value
        }else{
          df[i,j] <- prop.test(m,alternative = "less")$p.value
        }
      }
    }
  }
  
  assign(paste(deparse(substitute(x)),"p_values",sep="_"),
         df,envir = .GlobalEnv)
}

f_p_test(fisher_table)

x <- unlist(fisher_table_p_values)
p.adjust(x,"BH")


###Fisher tests for gene vs groups###
data_spread <- data_mca_second %>% 
  select(c(1,4:12))
# gather(Gene,Presence,c(2:7))

fisher_table <- table(data_mca_second$Lineage,data_mca_second$qacE)
f_p_test(fisher_table)
x <- unlist(fisher_table_p_values)
a <- as.data.frame(p.adjust(x[0:24],"BH"))

fisher_table <- table(data_mca_second$Lineage,data_mca_second$aac)
f_p_test(fisher_table)
x <- unlist(fisher_table_p_values)
a <- as.data.frame(p.adjust(x[0:24],"BH"))

fisher_table <- table(data_mca_second$Lineage,data_mca_second$aph)
f_p_test(fisher_table)
x <- unlist(fisher_table_p_values)
a <- as.data.frame(p.adjust(x[0:24],"BH"))

fisher_table <- table(data_mca_second$Lineage,data_mca_second$blaL1)
f_p_test(fisher_table)
x <- unlist(fisher_table_p_values)
a <- as.data.frame(p.adjust(x[0:24],"BH"))

fisher_table <- table(data_mca_second$Lineage,data_mca_second$blaL2)
f_p_test(fisher_table)
x <- unlist(fisher_table_p_values)
a <- as.data.frame(p.adjust(x[0:24],"BH"))

fisher_table <- table(data_mca_second$Lineage,data_mca_second$sul1)
f_p_test(fisher_table)
x <- unlist(fisher_table_p_values)
a <- as.data.frame(p.adjust(x[0:24],"BH"))

fisher_table <- table(data_mca_second$Lineage,data_mca_second$smoR)
f_p_test(fisher_table)
x <- unlist(fisher_table_p_values)
a <- as.data.frame(p.adjust(x[0:24],"BH"))

fisher_table <- table(data_mca_second$Lineage,data_mca_second$pilU)
f_p_test(fisher_table)
x <- unlist(fisher_table_p_values)
a <- as.data.frame(p.adjust(x[0:24],"BH"))

fisher_table <- table(data_mca_second$Lineage,data_mca_second$katA)
f_p_test(fisher_table)
x <- unlist(fisher_table_p_values)
a <- as.data.frame(p.adjust(x[0:24],"BH"))
