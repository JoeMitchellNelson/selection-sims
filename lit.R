lit <- read.csv("~/selection-sims/litsearch.csv")

str_detect(lit$Abstract,"representative.{1,50}sample|representative.{1,50}survey|representative of.{1,50}population") %>% sum()
str_detect(lit$Abstract,"selection") %>% sum()
str_detect(lit$Abstract,"non.{0,1}response bias") %>% sum()

a <- lit[which(str_detect(lit$Abstract,"selection")),c(9,22,31)]
View(a)
             