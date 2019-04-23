params <-
structure(list(child_path = "", setup_path = "../resources/"), .Names = c("child_path", 
"setup_path"))

## ----setup, include = FALSE, cahe = FALSE--------------------------------
source(paste0(params$setup_path, "knitr_setup.R"))
with_sol <- TRUE ## in order to control the output
with_course <- TRUE

## ----Rstudio, child = '../Rmd_files/Rstudio/presR.Rmd', cache = TRUE-----



## ----Rstudio, child = '../Rmd_files/Rstudio/presRStudio.Rmd', cache = TRUE----

## ----calcul_consol1------------------------------------------------------
1+1

## ----calcul_consol2------------------------------------------------------
c<-154*36
c

## ----install_package1, eval = FALSE--------------------------------------
## install.packages('ggplot2')

## ----install_package_2, eval = FALSE-------------------------------------
## if ( ! require('FactoMineR')) install.packages('FactoMineR)

## ----library_1-----------------------------------------------------------
library('FactoMineR')

## ----help_package_1,eval=FALSE-------------------------------------------
## help(dplyr) # lance l'aide associée à la commande dplyr
## help.start()  # lande l'aide HTML

## ----vignette_package,eval=FALSE-----------------------------------------
## browseVignettes("dplyr")


## ----Rstudio, child = '../Rmd_files/Rstudio/RapportAutomatise.Rmd', cache = TRUE----



## ----Rstudio, child = '../Rmd_files/Rstudio/LesObjetsR.Rmd', cache = TRUE----

## ----objets--------------------------------------------------------------
X<-c(1:5,10,12) # création d'un vecteur intitulé X
X # affichage de X
is.vector(X) # X est-il un vecteur?
class(X) # quelle est la classe de l'objet?

## ----opLogiq, eval= FALSE------------------------------------------------
## ?Comparison
## # et
## ?base::Logic


## ----Manip, child = '../Rmd_files/ManipulationDonnees/Importation.Rmd', cache = TRUE----

## ----import_1------------------------------------------------------------
mon_fichier <- "../../Datasets/SamaresEq.txt"
SamaresEq_base <- read.table(file = mon_fichier, sep = " ", header = TRUE, dec = '.')
head(SamaresEq_base, n=3)

## ----import_2------------------------------------------------------------
mon_fichier <- "../../Datasets/SamaresEq.txt"
SamaresEq_readr <- read_delim(file = mon_fichier, delim = " ")
SamaresEq_readr

## ----import_3------------------------------------------------------------
monFichier <- 'https://husson.github.io/img/decathlon.csv'
decathlon <- read.table(file = monFichier, header = TRUE, sep = ';')

## ----import_exo_2--------------------------------------------------------
library(readr)
vins_readr <- read_delim("http://factominer.free.fr/livre/vins.csv",
                                delim = ";",
                                escape_double = FALSE, locale = locale(encoding = "latin1"),
                                trim_ws = TRUE)

## ---- echo = FALSE-------------------------------------------------------
SamaresEq_readr <- read_delim(file = 'https://marieetienne.github.io/datasets/SamaresEq.txt', delim = " ")

SamaresEq_base <-  read.table(file = 'https://marieetienne.github.io/datasets/SamaresEq.txt', sep = " ")


## ----affich1, eval = FALSE-----------------------------------------------
## SamaresEq_readr <- read_delim(file = 'https://marieetienne.github.io/datasets/SamaresEq.txt', delim = " ")
## SamaresEq_readr
## 
## SamaresEq_base

## ----affich2-------------------------------------------------------------
head(SamaresEq_base, n = 2)
tail(SamaresEq_base, n = 3)

## ----affich3-------------------------------------------------------------
SamaresEq_readr %>% print(n  = 2)

## ----sql1, eval=FALSE----------------------------------------------------
## library(RODBC)
## # Liste les tables de la base de données connectée
## sqlTables(connect_base, tableType = "TABLE")
## # Liste les champs de la table  DonneesTotales
## sqlColumns(connect_base, sqtable = "DonneesTotales")

## ----sql2, eval=FALSE----------------------------------------------------
## # execute une requete SQL
## OtoYFT <- sqlQuery( channel = connect_base,
##                     query =
## "
## SELECT * FROM DonneesTotales
## WHERE (DonneesTotales.ProblemeSp = 'Ok' AND DonneesTotales.REC_Sp='Y'
## AND DonneesTotales.Otolithe = 'OT')
## ")
## # Liste les champs de la table  DonneesTotales
## sqlColumns(connect_base, sqtable = "DonneesTotales")


## ----Manip, child = '../Rmd_files/ManipulationDonnees/Manip_standard.Rmd', cache = TRUE----

## ---- echo = FALSE-------------------------------------------------------
list.files()
mon_fichier <- "../../Datasets/SamaresEq.txt"
SamaresEq_base <- read.table(file = mon_fichier, sep = " ", header = TRUE, dec = '.')

## ------------------------------------------------------------------------
dim(SamaresEq_base)

## ------------------------------------------------------------------------
colnames(SamaresEq_base)

## ------------------------------------------------------------------------
head(SamaresEq_base$Poids, n=5)

## ------------------------------------------------------------------------
SamaresEq_base[2, ]

## ------------------------------------------------------------------------
SamaresEq_base[c(2, 3, 7), ]


## ----Manip, child = '../Rmd_files/ManipulationDonnees/Manip_tidy.Rmd', cache = TRUE----

## ---- echo = FALSE-------------------------------------------------------
list.files()
mon_fichier <- "../../Datasets/SamaresEq.txt"
SamaresEq_base <- read.table(file = mon_fichier, sep = " ", header = TRUE, dec = '.')
SamaresEq_readr <- read_delim(file = mon_fichier, delim = " ")

## ---- eval = FALSE-------------------------------------------------------
## install.packages("tidyverse")

## ---- message = FALSE----------------------------------------------------
library(tidyverse)

## ------------------------------------------------------------------------
SamaresEq_base %>% filter( Surface > 0.75) -> Grand_samares
class(Grand_samares)
Grand_samares

## ----filter1-------------------------------------------------------------
SamaresEq_base %>% as_tibble %>% filter( Surface > 0.75) -> Grand_samares
class(Grand_samares)
Grand_samares

## ----filter2-------------------------------------------------------------
SamaresEq_readr %>% as_tibble %>% filter( Surface > 0.75) -> Grand_samares
class(Grand_samares)
Grand_samares

## ----  child= 'exo_filter_sol.Rmd', eval = with_sol----------------------

## ---- echo = FALSE-------------------------------------------------------
list.files()
mon_fichier <- "../../Datasets/SamaresEq.txt"
SamaresEq_base <- read.table(file = mon_fichier, sep = " ", header = TRUE, dec = '.')
SamaresEq_readr <- read_delim(file = mon_fichier, delim = " ")

## ----filter_ex_1, eval = FALSE-------------------------------------------
## SamaresEq_readr %>% filter(Site == 'Gornies')

## ----filter_ex_2, eval = FALSE-------------------------------------------
## SamaresEq_readr %>% filter(Largeur > 0.45 & Longueur < 2.32)

## ----filter_ex_3---------------------------------------------------------
SamaresEq_readr %>% filter(  ! NomSite %in% c('Gornies', 'StEtienne') )


## ----select1, eval = FALSE-----------------------------------------------
## SamaresEq_base %>% select(NomSite)

## ----select2, eval = FALSE-----------------------------------------------
## SamaresEq_base %>% select(-Site, -Arbre)

## ----select_ex1, eval = FALSE--------------------------------------------
## SamaresEq_base %>% filter( NomSite = 'Gornies') %>%  select(-Site, -Arbre) -> SamaresEq_base_gornies
## head(SamaresEq_base)

## ----mutate1-------------------------------------------------------------
SamaresEq_readr %>% 
  mutate(dispersion = Surface / Poids, 
         log_disp = log( dispersion )) -> SamaresEq_disp
SamaresEq_disp %>% select(-Site, -Arbre, -Distance, -CircArbre)

## ----mutate_ex1, eval = FALSE--------------------------------------------
## SamaresEq_readr %>% mutate(larg_x_long = Largeur * Longueur,
##                            diff_surf = larg_x_long - Surface)  %>%
##   select(-Site, -Arbre, -Distance, -CircArbre)

## ----mutate_ex2----------------------------------------------------------
SamaresEq_readr %>% mutate(larg_x_long = Largeur * Longueur,
                           diff_surf = larg_x_long - Surface)  %>%
  select(-Site, -Arbre, -Distance, -CircArbre)  %>%
  print(n = 3)


## ----summarise1----------------------------------------------------------
SamaresEq_readr %>% 
  summarise( longueur_m  = mean(Longueur, na.rm = TRUE)) 


## ----summarise2----------------------------------------------------------
SamaresEq_readr %>% 
  summarise_at( vars(Largeur, Longueur), funs(n(), median)) 


## ----summarise3----------------------------------------------------------
SamaresEq_readr %>% 
  summarise_if(is.numeric, mean, na.rm=TRUE) 

## ----summarise_ex1-------------------------------------------------------
SamaresEq_disp %>% 
  summarise_at( vars(Surface, dispersion),  funs( sd, mean), na.rm = TRUE) 

## ----group_by_1----------------------------------------------------------
SamaresEq_disp %>% group_by( NomSite) %>%
  summarise( Surface_m = mean (Surface)) %>%
  print(n = 3)

## ----group_by_2----------------------------------------------------------
  SamaresEq_disp %>% group_by( NomSite, Arbre ) %>%
  summarise( n_obs = n())  %>%   print(n = 3)

## ----group_by_ex1--------------------------------------------------------
  SamaresEq_disp %>% group_by( NomSite, Arbre ) %>%
  summarise( n_obs = n(), poids_m = mean(Poids))  %>%   print(n = 3)

## ----group_by_ex2--------------------------------------------------------
  SamaresEq_disp %>% group_by( NomSite) %>% summarise( n_Arbre = n_distinct(Arbre))

## ------------------------------------------------------------------------
write_csv(SamaresEq_disp, path = "../../Datasets/SamaresEq_disp.csv")


## ----Visu, child = '../Rmd_files/Visualisation/simple_plot.Rmd', cache = TRUE----

## ---- echo = FALSE, message = FALSE--------------------------------------
SamaresEq <- readr::read_delim(file = "../../Datasets/SamaresEq.txt", delim = " ")

invisible(SamaresEq_disp <- readr::read_delim(file = "../../Datasets/SamaresEq_disp.csv", , delim = ","))

## ----sp_scatter_1, out.height = '40%', out.width = '60%'-----------------
head(SamaresEq_disp)
plot(dispersion~Surface, data = SamaresEq_disp)

## ----simple_boxplot1, out.height = '40%', out.width = '60%', eval = FALSE----
## plot(dispersion~NomSite, data=SamaresEq_disp)

## ----simple_boxplot1_bis, out.height = '40%', out.width = '60%'----------
SamaresEq_disp %>% mutate(NomSite = as.factor(NomSite)) -> SamaresEq_disp
plot(dispersion~NomSite, data=SamaresEq_disp)


## ----Visu, child = '../Rmd_files/Visualisation/ggplot.Rmd', cache = TRUE----

## ---- echo = FALSE, message = FALSE--------------------------------------
SamaresEq <- readr::read_delim(file = "../../Datasets/SamaresEq.txt", delim = " ")

invisible(SamaresEq_disp <- readr::read_delim(file = "../../Datasets/SamaresEq_disp.csv", , delim = ","))

## ----scatter_1, out.height = '40%', out.width = '60%'--------------------
library(ggplot2)
ggplot(data = SamaresEq_disp) + aes( x = Surface, y = dispersion) + geom_point()

## ----boxplot1, out.height = '40%', out.width = '60%'---------------------
ggplot(data = SamaresEq_disp) + 
  aes( x = NomSite, y = dispersion) + 
  geom_boxplot()

## ----scatter_2,  out.height = '40%', out.width = '60%'-------------------
ggplot(data = SamaresEq_disp) + 
  aes( x = Surface, y = dispersion, col = NomSite) +
  geom_point( shape = 'a') 

## ----scatter_3,  out.height = '40%', out.width = '60%'-------------------
ggplot(data = SamaresEq_disp) + 
  aes( x = Surface, y = dispersion, col = NomSite) +
  geom_point( shape = 'a') +
  geom_smooth()

## ----scatter_4,  out.height = '40%', out.width = '60%'-------------------
ggplot(data = SamaresEq_disp) + 
  aes( x = Surface, y = dispersion, col = NomSite) + 
  geom_smooth(se = FALSE)

## ----scatter_5,  out.height = '40%', out.width = '60%'-------------------
ggplot(data = SamaresEq_disp) + 
  aes( x = Surface, y = dispersion, col = NomSite) + 
  geom_point( alpha= 0.3) 

## ----scatter_6,  out.height = '40%', out.width = '60%'-------------------
ggplot(data = SamaresEq_disp) +
  aes( x = Surface, y = dispersion, col = NomSite) +
  geom_point( alpha= 0.3) +
  geom_smooth( se = FALSE)

## ----scatter_7,  out.height = '40%', out.width = '60%'-------------------
ggplot(data = SamaresEq_disp) +
  aes( x = Surface, y = dispersion, col = NomSite) +
  geom_point( alpha= 0.3) +
  geom_smooth(method = 'lm', se = FALSE)

## ----color_1,  out.height = '40%', out.width = '60%'---------------------
ggplot(data = SamaresEq_disp) +
  aes( x = Surface, y = dispersion, col = NomSite) +
  geom_point( alpha= 0.3) +
  geom_smooth(method = 'lm', se = FALSE) +
  scale_color_viridis_d()

## ----grey_1,  out.height = '40%', out.width = '60%'----------------------
ggplot(data = SamaresEq_disp) +
  aes( x = Surface, y = dispersion, col = NomSite) +
  geom_point( alpha= 0.3) +
  geom_smooth(method = 'lm', se = FALSE) +
  scale_color_grey()

## ----color_2,  out.height = '40%', out.width = '60%'---------------------
palette <- c('#FF2F24', '#FF5B24', '#FF8324', '#FFBb00', '#FF9100', '#B0A64F', '#4FB06C' )

ggplot(data = SamaresEq_disp) +
  aes( x = Surface, y = dispersion, col = NomSite) +
  geom_point( alpha= 0.3) +
  geom_smooth(method = 'lm', se = FALSE) +
  scale_color_manual(values = palette)

## ----facet_1,  out.height = '50%'----------------------------------------
ggplot(data = SamaresEq_disp) +
  aes( x = Surface, y = dispersion, col = NomSite) +
  facet_wrap(~NomSite) + theme( legend.position = 'none' ) +
  geom_point( alpha= 0.3) +
  geom_smooth(method = 'lm', se = FALSE) +
  scale_color_manual(values = palette)

## ----axe_1,  out.height = '40%', out.width = '60%'-----------------------
ggplot(data = SamaresEq_disp) +
  aes( x = Surface, y = dispersion, col = NomSite) +
  geom_point( shape = 'a')  +
  labs(y = 'Dispersion du samare', x = 'Surface du samare')

## ----axe_2,  out.height = '40%', out.width = '60%'-----------------------
ggplot(data = SamaresEq_disp) +
  aes( x = Surface, y = dispersion, col = NomSite) +
  geom_point( shape = 'a')  +
  labs(y = expression ("Dispersion in"~g.cm^{-2}), x =  expression ("Surface in"~cm^{2}))

## ----exo_ggplot0---------------------------------------------------------
biomass <- readr::read_csv(file = '../../Datasets/Biomass_diversity.csv')


## ----exo_ggplot0_bis-----------------------------------------------------
biomass
class(biomass$YEAR)

## ----exo_ggplot0_ter-----------------------------------------------------
biomass %>%  mutate(Year_fact = as.factor(YEAR)) -> biomass

## ----exo_ggplot1, out.width = '60%'--------------------------------------
p1 <- biomass %>% ggplot() + aes(y = HARV_YIELD, x = COUNTRY ) + geom_boxplot() +
  xlab('') + theme( text = element_text(size=8))
p2 <- biomass %>% ggplot() +  aes(y = HARV_YIELD, x = COUNTRY) +
  geom_boxplot(aes(fill = Year_fact))  + theme( text = element_text(size=8))
library(ggpubr)
ggarrange(p1, p2, nrow= 2,  common.legend = TRUE, legend = 'right')

## ----exo_ggplot1_alt, out.width = '60%'----------------------------------
 biomass %>% ggplot() +  aes(y = HARV_YIELD, x = COUNTRY) +
  facet_wrap(~Year_fact)  +
  geom_boxplot(aes(fill = Year_fact))  + theme( text = element_text(size=8))

## ----exo_ggplot2, out.width = '50%'--------------------------------------
biomass %>% ggplot() + aes(y = HARV_YIELD, x = H ) + geom_point()

## ----exo_ggplot3, out.width = '50%'--------------------------------------
biomass %>% ggplot() +
  aes(y = HARV_YIELD, x = H ) +
  geom_point(aes(col = COUNTRY, shape = Year_fact))

## ----exo_ggplot4, out.width = '50%'--------------------------------------
biomass %>% ggplot() +   aes(y = HARV_YIELD, x = H) +
  geom_point(  aes(col = COUNTRY, shape = Year_fact ), alpha = 0.5 ) +
  geom_smooth(method="lm", se= F, size = 0.5, aes(col = COUNTRY, group = COUNTRY)) +
  geom_smooth(method = 'lm',size = 1, linetype = 'dashed', colour = 'black', se = F)

## ----exo_ggplot4_bis, out.width = '50%'----------------------------------
biomass %>% ggplot() +   aes(y = HARV_YIELD, x = H) +
  geom_point(  aes(col = COUNTRY, shape = Year_fact ), alpha = 0.5 ) + facet_wrap(~COUNTRY) +
  geom_smooth(method="lm", se= F, size = 0.5, aes(col = COUNTRY, group = COUNTRY)) -> p

p  +   geom_smooth(method = 'lm',size = 1, linetype = 'dashed', colour = 'black', se = F)

## ----exo_ggplot4_ter, out.width = '50%'----------------------------------
reg_coef <- coef(lm(HARV_YIELD ~ H, data = biomass))

p + geom_abline( intercept = reg_coef[1], slope = reg_coef[2], linetype = 'dashed', colour = 'black')

## ----exo_ggplot5, out.width = '50%'--------------------------------------
biomass %>% ggplot() +   aes(y = HARV_YIELD, x = H) +
    geom_point(  aes(col = COUNTRY, shape = Year_fact ), alpha = 0.5 ) +
  geom_smooth(method="lm", se= F, size = 0.5, aes(col = COUNTRY, group = COUNTRY)) +
  geom_smooth(method = 'lm',size = 1, linetype = 'dashed', colour = 'black', se = F)  +
  labs( x = 'Indice de diversité de Shannon', y = 'Rendement')


## ----test, child = '../Rmd_files/test/test.Rmd'--------------------------

## ----fig.height=2--------------------------------------------------------
poulpe <- read.table("https://r-stat-sc-donnees.github.io/poulpe.csv", header=TRUE, sep=";")
summary(poulpe)

## ---- out.height="80px", fig.height=3,fig.width=4------------------------
library(ggplot2)
poulpe %>% ggplot() + aes(x=Sexe,y=Poids) + geom_boxplot(fill=c("pink","lightblue"))

## ---- eval=FALSE, fig.height=3-------------------------------------------
## library(plotly)
## poulpe %>% ggplot() + aes(x=Sexe,y=Poids) + geom_boxplot(fill=c("pink","lightblue"))
## ggplotly()

## ---- eval=FALSE---------------------------------------------------------
## boxplot(Poids ~ Sexe, col=c("pink","lightblue"), data=poulpe)

## ---- eval=FALSE---------------------------------------------------------
## library(rAmCharts)
## amBoxplot(Poids ~ Sexe, col=c("pink","lightblue"), data=poulpe)

## ------------------------------------------------------------------------
by(poulpe$Poids, poulpe$Sexe, shapiro.test)

## ------------------------------------------------------------------------
var.test(Poids ~ Sexe, conf.level=.95, data=poulpe)

## ------------------------------------------------------------------------
res <- t.test(Poids~Sexe, alternative="two.sided", conf.level=.95, 
              var.equal=FALSE, data=poulpe)
res


## ----reg, child = '../Rmd_files/regression/regression.Rmd'---------------

## ----fig.height=2--------------------------------------------------------
ozone <- read.table("https://r-stat-sc-donnees.github.io/ozone.txt",header=TRUE)
library(tidyverse)
ozone.m <- ozone %>% select(1:11)
ozone.m %>% select(1:4) %>% summary()

## ----fig.height=2, eval=FALSE--------------------------------------------
## ozone <- read.table("https://r-stat-sc-donnees.github.io/ozone.txt",header=TRUE)
## ozone.m <- ozone[,1:11]
## summary(ozone.m[,1:4])

## ---- out.height="150px"-------------------------------------------------
library(GGally)
ozone.m %>% select(1:3) %>% ggpairs() 

## ---- eval=FALSE---------------------------------------------------------
## pairs(ozone.m[,1:3])

## ---- eval=FALSE---------------------------------------------------------
## reg.mul <- lm(maxO3~., data=ozone.m)
## summary(reg.mul)

## ---- echo=FALSE---------------------------------------------------------
reg.mul <- lm(maxO3~., data=ozone.m)
out <- capture.output(summary(reg.mul))
cat(out[-c(1,4:7)],sep="\n")

## ---- eval=FALSE---------------------------------------------------------
## library(FactoMineR)
## select <- RegBest(ozone.m$maxO3, ozone.m[,2:11])
## select$summary ; select$best

## ---- echo=FALSE---------------------------------------------------------
library(FactoMineR)
select <- RegBest(ozone.m$maxO3, ozone.m[,2:11]) 
select$summary
out <- capture.output(select$best)
cat(out[-c(1:7)],sep="\n")

## ------------------------------------------------------------------------
reg.fin <- lm(maxO3~T12+Ne9+Vx9+maxO3v, data=ozone.m)
summary(reg.fin)

## ------------------------------------------------------------------------
library(ggfortify)
autoplot(reg.fin)

## ---- fig.height=3,fig.width=4,out.height="120px"------------------------
residutib <- tibble(jour = 1:112, residu = rstudent(reg.fin))
residutib %>% ggplot() + aes(x=jour, y=residu) + geom_point() + 
  labs(x="Jour", y="Résidu", title = "Graphe des résidus studentisés") +
  geom_abline(slope=0, intercept=c(-2,0,2), linetype=c(2,1,2)) +
  geom_rect(aes(xmin=0, xmax=113, ymin=-2, ymax=2), alpha=0.002,fill="green") +
  geom_point(data = residutib %>% filter(abs(residu)>2), cex=2, col="red")

## ---- eval=FALSE---------------------------------------------------------
## plot(residu,pch=15,cex=.5,ylab="Résidus",main="Graphe des résidus studentisés",ylim=c(-3,3))
## abline(h=c(-2,0,2),lty=c(2,1,2))

## ------------------------------------------------------------------------
xnew <- matrix(c(19,8,2.05,70),nrow=1)
colnames(xnew) <- c("T12","Ne9","Vx9","maxO3v")
xnew <- as.data.frame(xnew)
predict(reg.fin,xnew,interval="pred")


## ----anova, child = '../Rmd_files/anova/anova.Rmd'-----------------------

## ----anova1, fig.height=2------------------------------------------------
ozone <- read.table("https://r-stat-sc-donnees.github.io/ozone.txt",header=TRUE)
summary(ozone[,c("maxO3","vent","pluie")])

## ---- fig.height=3, fig.width=10-----------------------------------------
library(ggplot2)
ozone %>% ggplot() + aes(y=maxO3, x=vent) + geom_point(aes(col=pluie, shape=vent))

## ---- fig.height=3,fig.width=10------------------------------------------
ozone %>% ggplot() + aes(pluie, maxO3) + geom_boxplot(aes(fill=vent)) +
  scale_fill_manual(values=c("lightblue","orange","green","grey"))

## ----anova_visu1,, eval=FALSE--------------------------------------------
## boxplot(maxO3~vent, data = ozone)

## ----anova_visu2, echo=FALSE,fig.height=3,fig.width=12-------------------
par(mar=c(2.5,2.5,0.5,0.5))
boxplot(maxO3~vent, data = ozone)

## ----anova_visu3, eval=FALSE---------------------------------------------
## boxplot(maxO3~vent*pluie, data = ozone, col=c(rep("Lightblue",4),rep("orange",4)))

## ----anova_visu4, echo=FALSE,fig.height=3,fig.width=12-------------------
par(mar=c(2.5,2.5,0.1,0.5))
boxplot(maxO3~vent*pluie, data = ozone, col=c(rep("lightblue",4),rep("orange",4)))

## ----anova_visu_gg3, fig.height=3, fig.width=6, out.height="120px"-------
ozone %>% ggplot() + aes(x = vent, y = maxO3, group = pluie) +
  geom_point(aes(color = pluie, shape=vent)) + 
  stat_summary(fun.y = mean, geom = "point", size=3, shape=15,aes(color = pluie)) +
  stat_summary(fun.y = mean, geom = "line", aes(color = pluie))

## ----anova_visu_gg4, fig.height=3, fig.width=6, out.height="70px"--------
ozone %>%  ggplot() + aes(x = pluie, y = maxO3, group = vent, color = vent, shape=pluie) +
  geom_point(alpha=0.5) + stat_summary(fun.y = mean, geom = "point", size=3, shape=15) +
  stat_summary(fun.y = mean, geom = "line")

## ---- eval=FALSE---------------------------------------------------------
## with(ozone,interaction.plot(vent,pluie,maxO3,col=1:nlevels(pluie)))
## with(ozone,interaction.plot(pluie,vent,maxO3,col=1:nlevels(vent)))

## ---- echo=FALSE, fig.width=12,fig.height=5------------------------------
par(mfrow=c(1,2))
with(ozone,interaction.plot(vent,pluie,maxO3,col=1:nlevels(pluie)))
with(ozone,interaction.plot(pluie,vent,maxO3,col=1:nlevels(vent)))

## ----anova_mod1_hypo-----------------------------------------------------
library(ggfortify)
mod.interaction <- lm(maxO3 ~ vent + pluie + vent:pluie, data=ozone)
autoplot(mod.interaction)

## ----anova_mod1_hypo_log-------------------------------------------------
library(ggfortify)
ozone %>% mutate(log_maxO3 = log(maxO3)) -> ozone
mod.interaction <- lm(log_maxO3 ~ vent + pluie + vent:pluie, data=ozone)
autoplot(mod.interaction)

## ----anova_mod1_modcomp--------------------------------------------------
mod.interaction <- lm(log_maxO3 ~ vent + pluie + vent:pluie, data=ozone)
mod.0 <- lm(log_maxO3 ~ 1, data=ozone)
anova(mod.0, mod.interaction)

## ----anova_mod1----------------------------------------------------------
anova(mod.interaction)
Anova(mod.interaction)

## ----anova_mod2----------------------------------------------------------
modele_12 <- lm(log_maxO3 ~ vent + pluie, data = ozone)
anova(modele_12)
Anova(modele_12)

## ------------------------------------------------------------------------
ozone %>%   ggplot() +  
  geom_point( mapping = aes(x=vent, y=log_maxO3))+  
  geom_boxplot( mapping = aes(x=vent, y=log_maxO3), alpha=0.3, fill='gray') 

## ------------------------------------------------------------------------
ozone %>%  ggplot() +  facet_wrap(~pluie)+
  geom_point( mapping = aes(x=vent, y=log_maxO3, col = pluie)) +   
  geom_boxplot( mapping = aes(x=vent, y=log_maxO3, fill = pluie), alpha=0.3) 

## ------------------------------------------------------------------------
ozone %>% group_by(vent, pluie) %>% summarise(n_obs = n()) 

## ----anova_mod3----------------------------------------------------------
summary(mod.interaction)

## ----anova_mod4----------------------------------------------------------
summary(modele_12)

## ------------------------------------------------------------------------
library('emmeans')
emmeans(modele_12,pairwise~pluie,adjust="hochberg")
emmeans(modele_12,pairwise~vent,adjust="hochberg")


## ----acp, child = '../Rmd_files/acp/acp.Rmd'-----------------------------

## ----config_acp, echo=FALSE----------------------------------------------
knitr::opts_chunk$set(cache = TRUE)

## ----fig.height=2--------------------------------------------------------
decath <- read.table("https://r-stat-sc-donnees.github.io/decathlon.csv",
                     sep=";",dec=".", header=TRUE, row.names=1, check.names=FALSE)

## ---- eval=FALSE---------------------------------------------------------
## library(FactoMineR)
## res.pca <- PCA(decath, quanti.sup=11:12, quali.sup=13)

## ---- eval=FALSE---------------------------------------------------------
## library(Factoshiny)
## res <- PCAshiny(decath)

## ---- echo=FALSE, fig.height=8,fig.width=8, out.width=c('49%', '49%'), fig.show='hold'----
library(FactoMineR)
res.pca <- PCA(decath, quanti.sup=11:12, quali.sup=13, graph=FALSE)
plot(res.pca,new.plot=FALSE)
plot(res.pca,choix="var",new.plot=FALSE)

## ---- fig.height=8,fig.width=8, out.width=c('49%', '49%'), fig.show='hold'----
plot(res.pca,habillage=13, cex=0.9, title="Graphe des individus")
plot(res.pca,choix="var", title="Graphe des variables")

## ---- eval=FALSE---------------------------------------------------------
## summary(res.pca, ncp=2)

## ---- echo=FALSE---------------------------------------------------------
out <- capture.output(summary(res.pca, ncp=2))
cat(out[2:27],sep="\n")

## ---- eval=FALSE---------------------------------------------------------
## summary(res.pca, ncp=2)

## ---- echo=FALSE---------------------------------------------------------
cat(out[29:50],sep="\n")

## ------------------------------------------------------------------------
dimdesc(res.pca, axes=1:2)


## ----paleoclimato, child = '../Rmd_files/exercice/paleoclimato.Rmd'------

## ----bioimport, result=FALSE---------------------------------------------
ss700 <- read.table("https://husson.github.io/img/ss700.csv", header=TRUE, 
                    sep=";", row.names=1)

## ----acp_bioclimato, fig.height=8,fig.width=8, out.width=c('48%', '48%'), fig.show='hold'----
library(FactoMineR)
res.pca <- PCA(ss700,quanti.sup=32:40, quali.sup=41,graph=FALSE)
plot(res.pca,hab=41,label="quali",cex=0.8)
plot(res.pca,choix="var",cex=0.8)

## ------------------------------------------------------------------------
dimdesc(res.pca)

## ----reg_bioclimato------------------------------------------------------
library(FactoMineR)
mod <- RegBest(ss700[,"tann"],ss700[,1:31])
mod$summary

## ----reg_bioclimato_best, eval=FALSE-------------------------------------
## mod$best

## ---- echo=FALSE---------------------------------------------------------
out <- capture.output(mod$best)
cat(out[-c(1:8)],sep="\n")

## ------------------------------------------------------------------------
library(FactoMineR)
mod <- AovSum(tann ~ biome,data=ss700)
mod

## ---- fig.width=12,fig.height=6, eval=FALSE------------------------------
## library(leaflet)
## pal <- colorNumeric(palette=c(low="blue",high="red"),domain=ss700["tann"])
## m <- leaflet() %>% addTiles() %>%
##   addCircles(ss700[,"long"],ss700[,"lati"], color=pal(ss700[,"tann"]),
##              fillOpacity=1, opacity=1)
## m


