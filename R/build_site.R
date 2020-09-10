
# Lister les packages utilises dans le site -------------------------------

library(stringr)

dta <-   readLines(con = 'used_package.txt') 

dta <- lapply(dta,  function(s_){
  prov  <- str_extract(string = s_, pattern = 'library\\([:graph:]+\\)') 
  prov <- str_remove(string = prov, pattern = "library\\(")
  prov <- str_remove(string = prov, pattern = "\\)")
  prov<- str_replace_all(string = prov, pattern = "[:punct:]+", replacement = '')
  prov <- str_subset( string = prov, pattern = "NA", negate = TRUE)
}
) 



# Installer les packages nécessaires au site ------------------------------

lapply( do.call('c', dta), function(s_){
  install.packages(s_)
})



lapply(
  list.files('.', recursive = TRUE, pattern = '.Rmd', 
                                  full.names = TRUE), rmarkdown::render)
rmarkdown::render_site(encoding = 'UTF-8')
