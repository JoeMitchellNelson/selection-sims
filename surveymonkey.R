require(pacman)
p_load(tidyverse,RODBC,readxl,data.table,openxlsx,stringdist,Rfast,foreign,grid,gridExtra,Rcpp,ggVennDiagram,gt)


devtools::install_github("tntp/surveymonkey")
#docs at: https://github.com/tntp/surveymonkey


library(surveymonkey)


# copy oauth key into r profile:
# usethis::edit_r_profile()
#options(sm_oauth_token = "BjVjlV9MiVBgfe1XpS2xPXr.FgBNQSlYu5mmQzr4i13EBhDS2qC-I30eIplY-WEf4JPycIjtTkQsalOTl91YR1552XMxb0tGwqyjTUHRP6kS9Lgt9hOZ769MsLnJWzlP")
#usethis::edit_r_profile()

# running this should return the key
# getOption("sm_oauth_token")


respath <- "S:/Offices/Salem (500 Summer St)/OHIT/EnvironmentalScan/CCOs/2021HITsurvey/Results/"
overlappath <- "S:/Offices/Salem (500 Summer St)/OHIT/EnvironmentalScan/CCOs/CCO_HIT_Reporting/HIT_DataFileV1a/Overlap"


# useful mainly for getting survey IDs
surveys <- browse_surveys(200)
surveys <- as.data.frame(surveys) %>% arrange(desc(nickname))
survey_ids <- surveys$id[which(surveys$title=="2021 Oregon Health IT Survey")]


# function to fix column names


fixcolnames <- function (x) {
  x <- trimws(x)
  x <- trimws(x)
  
  x <- x %>% str_remove_all("<strong>|<.strong>|<br>|<em>|<.em>| - NA$")
  x <- x %>% str_remove_all("<span style=\"text-decoration: underline;\">|</span>")
  x <- x %>% str_remove_all("</a>")
  x <- x %>% str_remove_all("<a href.{0,200}>")
  
  x <- ifelse(str_detect(x,"^Thank you"),"Thank you...additional Qs okay?",x)
  x <- ifelse(str_detect(x,"^Providers play"),"Health equity comments?",x)
  
  x
}
ohitfolks <- c(joem@uoregon.edu,Martamak@comcast.net)


linkcollector <- c(313886999,312776611)


resdflink <- NULL
resdf <- NULL


for (i in 1:length(survey_ids)) {
  a_survey_obj <- fetch_survey_obj(survey_ids[i])
  
  df <- parse_survey(a_survey_obj)
  
  if (survey_ids[i] %in% linkcollector & ncol(df)>1) {
    resdflink <- rbind(resdflink,df)
  } else {
    
    if (i == 1) {
      colnames <- names(df)
    }
    
    if (i > 1 & ncol(df)>1) {
      names(df) <- colnames
    }
    
    if (ncol(df)>1) {
      names(df) <- fixcolnames(names(df))
      resdf <- rbind(resdf,df)
    }
  }
}




names(resdflink) <- fixcolnames(names(resdflink))
resdflink$recipient_id <- NA
resdflink$custom_value2 <- NA
resdflink$email <- NA
names(resdflink)[which(names(resdflink)=="OrgID")] <- "custom_value"
names(resdflink)[which(names(resdflink)=="What type of HIT support would your organization like from the CCO/DCO(s) it contracts with or from OHA?")] <-
  "What type of HIT support would your organization like from the CCO(s) it contracts with or from OHA?"


resdflink <- resdflink %>% select(names(resdf))


results <- rbind(resdf,resdflink)
#rm(resdf)
results <- results %>% dplyr::filter(!email %in% ohitfolks)


results$minutes_spent <- ((results$date_modified - results$date_created)/60) %>% as.numeric()


results <- results %>% dplyr::select(custom_value2,custom_value,`Organization name`,email,response_status,minutes_spent,everything())
names(results)[1:2] <- c("CCO","OHA Mastername(s)")
names(results) <- names(results) %>% str_remove_all(\\t)


#results_unsplit <- results
#results <- results_unsplit


#results <- results_unsplit %>% separate_rows(`OHA Mastername(s)`,sep="; ")


results$`OHA Mastername(s)` <- results$`OHA Mastername(s)` %>% str_replace_all("__"," & ") %>% str_replace_all("_"," ")


results$CCO <- ifelse(results$survey_id==313886999,"CapitolDental",results$CCO)
results$CCO <- ifelse(results$survey_id==312776611,"AdvancedHealth",results$CCO)


######### get provider type from overlap table ############


overlap <- read.xlsx(paste0(overlappath,"/overlap_ASSIGNED_20210922.xlsx")) %>%
  dplyr::filter(!is.na(OrgName)) %>% dplyr::select(1:4)


results <- results %>% mutate(id=row_number())


lil <- results %>% dplyr::select(`OHA Mastername(s)`,id) %>%
  separate_rows(`OHA Mastername(s)`,sep="; ")


lil <- left_join(lil,overlap,by=c("OHA Mastername(s)"="OrgName"))
lil <- lil %>% group_by(id) %>% summarise(PH=max(PH),BH=max(BH),OH=max(OH))


results <- left_join(results,lil,by="id")


results$reported_type <-
  paste0(results$`Which provider types within your organization use this EHR (select all that apply)? - Physical`,"/",
         results$`Which provider types within your organization use this EHR (select all that apply)? - Behavioral`,"/",
         results$`Which provider types within your organization use this EHR (select all that apply)? - Oral`) %>%
  str_replace_all("NA/NA/NA","Did not answer") %>%
  str_remove_all("/NA") %>% str_remove_all("NA/")


results$reported_type %>% as.factor() %>% summary()


results$typecode <- paste0(results$PH,results$BH,results$OH)


codekey <- data.frame(typecode=c("100","010","001","110","101","011","111"),
                      determined_type=c("Physical","Behavioral","Oral",
                                        "Physical/Behavioral","Physical/Oral","Behavioral/Oral",
                                        "Physical/Behavioral/Oral"))


results <- left_join(results,codekey)


names(results)[which(str_detect(names(results),"About how many"))] <- c("num_sites","num_providers","num_patients")




#########################################
######## summary table ##################
############ (excel cover sheet) ########
#########################################


# cover <- results %>% mutate(ehr_sat=6-as.numeric(`Overall, how satisfied is your organization with its EHR?`)) %>%
#   group_by(CCO) %>%
#   summarise(`Number of responses`=n(),
#             `Physical Health`=sum(PH),
#             `Behavioral Health`=sum(BH),
#             `Oral Health`=sum(OH),
#             `Orgs represented`=n()+sum(str_count(`OHA Mastername(s)`,";")),
#             `% use EHR` = round(100*sum(str_detect(`Does your organization currently use an electronic health record (EHR)?`,"^Yes"),na.rm=T)/n()),
#             `Median time`=paste0(round(median(minutes_spent))," minutes"),
#             `Answered additional Qs`=sum(`Thank you...additional Qs okay?`=="Yes",na.rm=T),
#             `Satisfied with EHR (5-pt scale)`=round(mean(ehr_sat,na.rm=T),2))
# coverall <- results %>% mutate(ehr_sat=6-as.numeric(`Overall, how satisfied is your organization with its EHR?`)) %>%
#   summarise(`Number of responses`=n(),
#             `Physical Health`=sum(PH),
#             `Behavioral Health`=sum(BH),
#             `Oral Health`=sum(OH),
#             `Orgs represented`=n()+sum(str_count(`OHA Mastername(s)`,";")),
#             `% use EHR` = round(100*sum(str_detect(`Does your organization currently use an electronic health record (EHR)?`,"^Yes"),na.rm=T)/n()),
#             `Median time`=paste0(round(median(minutes_spent))," minutes"),
#             `Answered additional Qs`=sum(`Thank you...additional Qs okay?`=="Yes",na.rm=T),
#             `Satisfied with EHR (5-pt scale)`=round(mean(ehr_sat,na.rm=T),2))
# coverall$CCO <- "All/Total"
# coverall <- coverall %>% dplyr::select(CCO,everything())
# cover <- rbind(cover,coverall)


cover <- results %>% mutate(ehr_sat=6-as.numeric(`Overall, how satisfied is your organization with its EHR?`)) %>%
  group_by(determined_type) %>%
  summarise(`Number of responses`=n(),
            `Complete responses` = sum(response_status=="completed"),
            `Partial responses` = sum(response_status=="partial"),
            `% use EHR` = round(100*sum(str_detect(`Does your organization currently use an electronic health record (EHR)?`,"^Yes"),na.rm=T)/n()),
            `Median time`=paste0(round(median(minutes_spent))," minutes"),
            `Answered additional Qs`=sum(`Thank you...additional Qs okay?`=="Yes",na.rm=T),
            `Satisfied with EHR (5-pt scale)`=round(mean(ehr_sat,na.rm=T),2))
coverall <- results %>% mutate(ehr_sat=6-as.numeric(`Overall, how satisfied is your organization with its EHR?`)) %>%
  summarise(`Number of responses`=n(),
            `Complete responses` = sum(response_status=="completed"),
            `Partial responses` = sum(response_status=="partial"),
            `% use EHR` = round(100*sum(str_detect(`Does your organization currently use an electronic health record (EHR)?`,"^Yes"),na.rm=T)/n()),
            `Median time`=paste0(round(median(minutes_spent))," minutes"),
            `Answered additional Qs`=sum(`Thank you...additional Qs okay?`=="Yes",na.rm=T),
            `Satisfied with EHR (5-pt scale)`=round(mean(ehr_sat,na.rm=T),2))
coverall$determined_type <- "All"
coverall <- coverall %>% dplyr::select(determined_type,everything())
cover <- rbind(cover,coverall)
names(cover)[1] <- "Provider type"


#########################################
###### plot some summary stats ##########
#########################################


#### Do you use EHR? by service type (bars) #######


ehrPH <- results %>% filter(PH==1 & !is.na(`Does your organization currently use an electronic health record (EHR)?`))
nPH <- nrow(ehrPH)
ehrPH <- count(ehrPH,`Does your organization currently use an electronic health record (EHR)?`,.drop=F)
ehrPH$total <- nPH
ehrPH$type <- "Physical"


ehrBH <- results %>% filter(BH==1 & !is.na(`Does your organization currently use an electronic health record (EHR)?`))
nBH <- nrow(ehrBH)
ehrBH <- count(ehrBH,`Does your organization currently use an electronic health record (EHR)?`,.drop=F)
ehrBH$total <- nBH
ehrBH$type <- "Behavioral"




ehrOH <- results %>% filter(OH==1 & !is.na(`Does your organization currently use an electronic health record (EHR)?`))
nOH <- nrow(ehrOH)
ehrOH <- count(ehrOH,`Does your organization currently use an electronic health record (EHR)?`,.drop=F)
ehrOH$total <- nOH
ehrOH$type <- "Oral"


ehruse <- rbind(ehrPH,ehrBH,ehrOH)


ehruse$Answer <- rep(c('Yes, fully','Yes, partially','Yes, changing vendors',
                       'No, selection stage','No, info gathering stage','No, but planning to',
                       'No, and no plans'),3)






ehruse <- ehruse %>% filter(!is.na(Answer)) %>% select(Answer,type, n,total)
ehruse$proportion <- ehruse$n/ehruse$total


multiples_data <- ehruse %>%
  dplyr::mutate(type_f = factor(type, levels = c("Physical", "Behavioral","Oral")),
                proportion = round(100*proportion),
                label_position = dplyr::if_else(proportion < 15, proportion+1, proportion-17),
                label_color = "white")


colors <- data.frame(
  facet_values = c("Physical", "Behavioral","Oral"),
  colors = c("#005595", "#ec8902","#A01c3f")
)


## Merge colors together with the original data, readying us for the graph.
input_data <- merge(multiples_data,
                    colors,
                    all.x = TRUE,
                    by.x = "type",
                    by.y = "facet_values")


input_data$label_color <- ifelse(input_data$proportion < 15,input_data$colors,input_data$label_color)


axis_order <- multiples_data %>%
  dplyr::filter(type == "Physical") %>%
  dplyr::select(Answer)


useplot <- ggplot(
  data = input_data,
  aes(
    x = Answer,
    y = proportion
  )
) +
  geom_bar(
    aes(
      fill = colors
    ),
    stat = "identity",
    position = "identity",
    width = 3/4,
    show.legend = F
  ) +
  scale_fill_manual(
    values = c(
      "#005595" = "#005595",
      "#ec8902" = "#ec8902",
      "#A01c3f" = "#A01c3f"
    )
  ) +
  facet_wrap(~type_f,strip.position = "bottom") +
  coord_flip() +
  scale_x_discrete(
    limits = rev(axis_order$Answer)
  ) +
  ## Stretch the bar graphs to align completely with strip text
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill="white"),
    # strip.text = element_text(
    #   face = "bold",
    #   hjust = 0,
    #   size = 12
    # ),
    strip.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(
      face = "bold",
      size = 13,
      color = "#646464",
      margin = margin(
        r = 8
      )
    ),
    axis.text.x = element_blank(),
    panel.spacing = unit(3, "lines"),
    axis.line.x = element_blank(),
    axis.line.y = element_line()
  ) +
  labs(x=NULL,y=NULL,title="") +
  geom_text(
    aes(
      group = type_f,
      label = paste0(proportion, "%"),
      color = label_color,
      y = label_position,
      hjust = 0
    ),
    family="sans",
    fontface="bold",
    size=4.2,
    position = position_nudge(
      x = 0,
      y = 1.5
    )
  ) +
  scale_color_manual(
    values = c(
      "white" = "white",
      "#005595" = "#005595",
      "#ec8902" = "#ec8902",
      "#A01c3f" = "#A01c3f"
    )
  )


# t1 <- textGrob(expression("Most " *
#                             phantom(bold("physical")) * ", " *
#                             phantom(bold("behavioral")) * ", " *
#                             "and " * phantom(bold("oral")) *
#                             " health providers use EHR"),
#                x = unit(0.5, "lines"), y = unit(-0.5, "lines"), hjust = 0, vjust = 0.3,  gp = gpar(col = "#646464",fontsize = 14))
#
# t2 <- textGrob(expression(phantom("Most ") *
#                             bold("physical") * phantom(", ") *
#                             phantom(bold("behavioral")) * phantom(", ") *
#                             phantom("and ") * phantom(bold("oral")) *
#                             phantom(" health providers use EHR")),
#                x = unit(0.5, "lines"), y = unit(-0.5, "lines"), hjust = 0, vjust = 0.3,  gp = gpar(col = "#005595",fontsize = 14))
#
# t3 <- textGrob(expression(phantom("Most ") *
#                             phantom(bold("physical")) * phantom(", ") *
#                             bold("behavioral") * phantom(", ") *
#                             phantom("and ") * phantom(bold("oral")) *
#                             phantom(" health providers use EHR")),
#                x = unit(0.5, "lines"), y = unit(-0.5, "lines"), hjust = 0, vjust = 0.3,  gp = gpar(col = "#ec8902",fontsize = 14))
#
# t4 <- textGrob(expression(phantom("Most ") *
#                             phantom(bold("physical")) * phantom(", ") *
#                             phantom(bold("behavioral")) * phantom(", ") *
#                             phantom("and ") * bold("oral") *
#                             phantom(" health providers use EHR")),
#                x = unit(0.5, "lines"), y = unit(-0.5, "lines"), hjust = 0, vjust = 0.3,  gp = gpar(col = "#A01c3f",fontsize = 14))
#
# t <- grobTree(t1,t2,t3,t4)
#
# b1 <- textGrob(paste0("N Physical = ", nPH,", N Behavioral = ", nBH,", N Oral = ",nOH),
#                x = unit(0.5, "lines"), y = unit(0.5, "lines"), hjust = -0, vjust = -0,  gp = gpar(col = "#646464",fontsize = 11))
#
#
# gg <- arrangeGrob(useplot,top=t,bottom=b1)
# useplot2 <- grid.arrange(gg)


useplot


ggsave(paste0(respath,"Images/EHRbars.png"),useplot,height=3,width = 8.8,units="in")


#####################################################
############# num sites bars ########################
#####################################################


sitesPH <- results %>% filter(PH==1 & !is.na(num_sites))
sitesPH <- count(sitesPH,num_sites,.drop=F)
nPH <- sum(sitesPH$n)
sitesPH <- sitesPH %>% mutate(type="Physical",total=nPH)




sitesBH <- results %>% filter(BH==1 & !is.na(num_sites))
sitesBH <- count(sitesBH,num_sites,.drop=F)
nBH <- sum(sitesBH$n)
sitesBH <- sitesBH %>% mutate(type="Behavioral",total=nBH)


sitesOH <- results %>% filter(OH==1 & !is.na(num_sites))
sitesOH <- count(sitesOH,num_sites,.drop=F)
nOH <- sum(sitesOH$n)
sitesOH <- sitesOH %>% mutate(type="Oral",total=nOH)


sites <- rbind(sitesPH,sitesBH,sitesOH)


names(sites)[1] <- "Answer"


sites <- sites %>% filter(!is.na(Answer)) %>% select(Answer,type, n,total)
sites$proportion <- sites$n/sites$total


multiples_data <- sites %>%
  dplyr::mutate(type_f = factor(type, levels = c("Physical", "Behavioral","Oral")),
                proportion = round(100*proportion),
                label_position = dplyr::if_else(proportion < 15, proportion+1, proportion-19),
                label_color = "white")


colors <- data.frame(
  facet_values = c("Physical", "Behavioral","Oral"),
  colors = c("#005595", "#ec8902","#A01c3f")
)


## Merge colors together with the original data, readying us for the graph.
input_data <- merge(multiples_data,
                    colors,
                    all.x = TRUE,
                    by.x = "type",
                    by.y = "facet_values")


input_data$label_color <- ifelse(input_data$proportion < 15,input_data$colors,input_data$label_color)


axis_order <- multiples_data %>%
  dplyr::filter(type == "Physical") %>%
  dplyr::select(Answer)


sitesplot <- ggplot(
  data = input_data,
  aes(
    x = Answer,
    y = proportion
  )
) +
  geom_bar(
    aes(
      fill = colors
    ),
    stat = "identity",
    position = "identity",
    width = 3/4,
    show.legend = F
  ) +
  scale_fill_manual(
    values = c(
      "#005595" = "#005595",
      "#ec8902" = "#ec8902",
      "#A01c3f" = "#A01c3f"
    )
  ) +
  facet_wrap(~type_f,strip.position = "bottom") +
  coord_flip() +
  scale_x_discrete(
    limits = rev(axis_order$Answer)
  ) +
  ## Stretch the bar graphs to align completely with strip text
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0,90)
  ) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    # strip.text = element_text(
    #   face = "bold",
    #   hjust = 0,
    #   size = 12
    # ),
    strip.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(
      face = "bold",
      size = 13,
      color = "#646464",
      margin = margin(
        r = 10
      )
    ),
    axis.text.x = element_blank(),
    panel.spacing = unit(3, "lines"),
    axis.line.x = element_blank(),
    axis.line.y = element_line()
  ) +
  labs(x=NULL,y=NULL,title="") +
  geom_text(
    aes(
      group = type_f,
      label = paste0(proportion, "%"),
      color = label_color,
      y = label_position,
      hjust = 0
    ),
    family="sans",
    size=4.2,
    fontface="bold",
    position = position_nudge(
      x = 0,
      y = 1.5
    )
  ) +
  scale_color_manual(
    values = c(
      "white" = "white",
      "#005595" = "#005595",
      "#ec8902" = "#ec8902",
      "#A01c3f" = "#A01c3f"
    )
  )


sitesplot


#ggsave(paste0(respath,"Images/sitesbars.pdf"),sitesplot,height=5,width = 8.8,units="in")


ggsave(paste0(respath,"Images/sitesbars.png"),sitesplot,height=3,width = 7.8,units="in")


#write.xlsx(results,"~/results_11_3_formarta.xlsx")


##### Response completeness (donut) #########


# Create test data.
data <- results %>% group_by(response_status) %>% summarise(count=n())
data$category <- str_to_title(data$response_status) %>% str_replace("Completed","Complete")
# Compute percentages
data$fraction <- data$count / sum(data$count)


# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)


# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))


# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2


# Compute a good label
data$label <- paste0("   ",data$category, ": ", data$count,"     \n ")


# Make the plot
completenessplot <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text( x=1.3, aes(y=labelPosition, label=label, color=category), size=3.5) + # x here controls label position (inner / outer)
  scale_fill_manual(values=c("Complete"="#005595","Partial"= "#ec8902")) +
  scale_color_manual(values=c("Complete"="#005595","Partial"= "#ec8902")) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")


t <- textGrob(label="Responses",
              x = unit(0.5, "lines"), y = unit(0.0, "lines"), hjust = 0, vjust = 0.3,  gp = gpar(col = "#646464",fontsize = 14))


completenessplot <- grid.arrange(arrangeGrob(completenessplot,top=t))


ggsave(paste0(respath,"Images/completeness.png"),completenessplot,height=3,width=3,units="in")








#### which EHR benefits have you experienced? (bars) #####


#names(results)[which(str_detect(names(results), "Which benefits related to EHR use") & !str_detect(names(results),"Comments\\?"))]


# cols = which(str_detect(names(results), "Which benefits related to EHR use"))
# results[,cols] = apply(results[,cols], 2, function(x) as.character(x))


a <- results %>% dplyr::filter(`Thank you...additional Qs okay?`=="Yes") %>%
  select(CCO,starts_with("Which provider types within"),starts_with("Which benefits related to EHR use")) %>% as.data.frame()
a <- a %>% dplyr::select(-c(5:9))
a[,2:ncol(a)] <- apply(a[,2:ncol(a)], 2, function(x) as.character(x))


a[,2:ncol(a)] <- a[,2:ncol(a)] %>% replace(!is.na(.),1)
a[,2:ncol(a)] <- a[,2:ncol(a)] %>% replace(is.na(.),0)
names(a) <- str_remove_all(names(a),"Which benefits related to EHR use has your organization experienced\\? \\(Select all that apply\\) - ") %>%
  str_remove_all(" \\(please specify\\)") %>% trimws()
names(a) <- str_remove_all(names(a),"Which provider types within your organization use this EHR \\(select all that apply\\)\\? - ") %>%
  str_remove_all(" \\(please specify\\)") %>% trimws()
ehrbenefits <- a


nPH2 <- sum(as.numeric(a$Physical))
nBH2 <- sum(as.numeric(a$Behavioral))
nOH2 <- sum(as.numeric(a$Oral))


beneph <- ehrbenefits %>% dplyr::filter(Physical==1) %>%
  dplyr::select(-Physical,-Behavioral,-Oral,-CCO) %>%
  t() %>% as.data.frame() %>%  rownames_to_column(var="Benefit") %>%
  mutate(type="Physical") %>% dplyr::select(Benefit,type,everything())
beneph$proportion <- apply(beneph[,-(1:2)],1,function (x) mean(as.numeric(x)))
beneph <- beneph %>% dplyr::select(-starts_with("V"))


benebh <- ehrbenefits %>% dplyr::filter(Behavioral==1) %>%
  dplyr::select(-Physical,-Behavioral,-Oral,-CCO)  %>%
  t() %>% as.data.frame() %>%  rownames_to_column(var="Benefit") %>%
  mutate(type="Behavioral") %>% dplyr::select(Benefit,type,everything())
benebh$proportion <- apply(benebh[,-(1:2)],1,function (x) mean(as.numeric(x)))
benebh <- benebh %>% dplyr::select(-starts_with("V"))


beneoh <- ehrbenefits %>% dplyr::filter(Oral==1) %>%
  dplyr::select(-Physical,-Behavioral,-Oral,-CCO)  %>%
  t() %>% as.data.frame() %>%  rownames_to_column(var="Benefit") %>%
  mutate(type="Oral") %>% dplyr::select(Benefit,type,everything())
beneoh$proportion <- apply(as.data.frame(beneoh[,-(1:2)]),1,function (x) mean(as.numeric(x)))
beneoh <- beneoh %>% dplyr::select(-starts_with("V"))
if (is.nan(beneoh$proportion[1])) {
  beneoh$proportion <- runif(16,min=0,max=1)
}


ehrbenefits <- rbind(beneph,benebh,beneoh)


multiples_data <- ehrbenefits %>%
  dplyr::mutate(type_f = factor(type, levels = c("Physical", "Behavioral","Oral")),
                proportion = round(100*proportion),
                label_position = dplyr::if_else(proportion < 15, proportion, 0),
                label_color = "white")


colors <- data.frame(
  facet_values = c("Physical", "Behavioral","Oral"),
  colors = c("#005595", "#ec8902","#A01c3f")
)


## Merge colors together with the original data, readying us for the graph.
input_data <- merge(multiples_data,
                    colors,
                    all.x = TRUE,
                    by.x = "type",
                    by.y = "facet_values")


input_data$label_color <- ifelse(input_data$proportion < 15,input_data$colors,input_data$label_color)


axis_order <- multiples_data %>%
  dplyr::filter(type == "Physical") %>%
  dplyr::arrange(proportion) %>%
  dplyr::select(Benefit)


beneplot <- ggplot(
  data = input_data,
  aes(
    x = Benefit,
    y = proportion
  )
) +
  geom_bar(
    aes(
      fill = colors
    ),
    stat = "identity",
    position = "identity",
    width = 3/4,
    show.legend = F
  ) +
  scale_fill_manual(
    values = c(
      "#005595" = "#005595",
      "#ec8902" = "#ec8902",
      "#A01c3f" = "#A01c3f"
    )
  ) +
  facet_wrap(~type_f) +
  coord_flip() +
  scale_x_discrete(
    limits = axis_order$Benefit
  ) +
  ## Stretch the bar graphs to align completely with strip text
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    # strip.text = element_text(
    #   face = "bold",
    #   hjust = 0,
    #   size = 12
    # ),
    strip.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(
      face = "bold",
      size = 13,
      color = "#646464",
      margin = margin(
        r = 10
      )
    ),
    axis.text.x = element_blank(),
    panel.spacing = unit(3, "lines"),
    axis.line.x = element_blank(),
    axis.line.y = element_line()
  ) +
  labs(x=NULL,y=NULL,title="") +
  geom_text(
    aes(
      group = type_f,
      label = paste0(proportion, "%"),
      color = label_color,
      y = label_position,
      hjust = 0
    ),
    family="sans",
    size=4.2,
    fontface="bold",
    position = position_nudge(
      x = 0,
      y = 2
    )
  ) +
  scale_color_manual(
    values = c(
      "white" = "white",
      "#005595" = "#005595",
      "#ec8902" = "#ec8902",
      "#A01c3f" = "#A01c3f"
    )
  )


beneplot


#Reported benefits of EHRs for Physical, Behavioral, and Oral providers


# t1 <- textGrob(expression("Reported benefits of EHRs for " *
#                             phantom(bold("physical")) * ", " *
#                           phantom(bold("behavioral")) * ", " *
#                             "and " * phantom(bold("oral")) *
#                             " health providers"),
#                x = unit(0.5, "lines"), y = unit(-0.5, "lines"), hjust = 0, vjust = 0.3,  gp = gpar(col = "#646464",fontsize = 14))
#
# t2 <- textGrob(expression(phantom("Reported benefits of EHRs for ") *
#                             bold("physical") * phantom(", ") *
#                             phantom(bold("behavioral")) * phantom(", ") *
#                             phantom("and ") * phantom(bold("oral")) *
#                             phantom(" health providers")),
#                x = unit(0.5, "lines"), y = unit(-0.5, "lines"), hjust = 0, vjust = 0.3,  gp = gpar(col = "#005595",fontsize = 14))
#
# t3 <- textGrob(expression(phantom("Reported benefits of EHRs for ") *
#                             phantom(bold("physical")) * phantom(", ") *
#                             bold("behavioral") * phantom(", ") *
#                             phantom("and ") * phantom(bold("oral")) *
#                             phantom(" health providers")),
#                x = unit(0.5, "lines"), y = unit(-0.5, "lines"), hjust = 0, vjust = 0.3,  gp = gpar(col = "#ec8902",fontsize = 14))
#
# t4 <- textGrob(expression(phantom("Reported benefits of EHRs for ") *
#                             phantom(bold("physical")) * phantom(", ") *
#                             phantom(bold("behavioral")) * phantom(", ") *
#                             phantom("and ") * bold("oral") *
#                             phantom(" health providers")),
#                x = unit(0.5, "lines"), y = unit(-0.5, "lines"), hjust = 0, vjust = 0.3,  gp = gpar(col = "#A01c3f",fontsize = 14))
#
# t <- grobTree(t1,t2,t3,t4)
#
# b1 <- textGrob(paste0("N Physical = ", nPH2,", N Behavioral = ", nBH2,", N Oral = ",nOH2),
#                x = unit(0.5, "lines"), y = unit(0.5, "lines"), hjust = -0, vjust = -0,  gp = gpar(col = "#646464",fontsize = 11))
#
# gg <- arrangeGrob(beneplot,top=t,bottom=b1)
# beneplot2 <- grid.arrange(gg)


beneplot


ggsave(paste0(respath,"Images/benefits.png"),beneplot,height=5,width=12,units="in")


##################################################
######## response nested boxes ###################
##################################################


overlap <- read.xlsx(paste0(respath,"Overlap_responses_2021-11-02.xlsx"))


cap <- read.xlsx("S:/Offices/Salem (500 Summer St)/OHIT/EnvironmentalScan/CCOs/CCO_HIT_Reporting/HIT_DataFileV1a/Email Lists/Capitol Dental URLs.xlsx")
ah <- read.xlsx(paste0(overlappath,"/advancedhealth_unformatted.xlsx"))
doneso <- list.files(writepath,pattern="[[:digit:]]",full.names=T)
doneso <- doneso[which(!str_detect(doneso,"US_for_SM|US2_for_SM"))]


sent <- NULL
for (i in 1:length(doneso)) {
  temp <- read.csv(doneso[i])
  sent <- rbind(sent,temp)
}


sent$Email <- sent$Email %>% tolower() %>% trimws(whitespace = "[\\h\\v]")


nbig <- length(unique(overlap$Email)) + sum(is.na(overlap$Email))
nemails <- length(unique(sent$Email)) + sum(str_detect(cap$Link,"https"),na.rm=T) + sum(str_detect(ah$link,"https"),na.rm=T)
ncomplete <- sum(results$response_status=="completed")
nextra <- sum(results$`Thank you...additional Qs okay?`=="Yes",na.rm=T)


define_x <- function(
  survey_responses_count,
  parent_responses_count,
  affirmative_responses_count
) {
  parent_max <- parent_responses_count / survey_responses_count
  proportion_of_parent <- affirmative_responses_count / parent_responses_count
  x_max <- parent_max * proportion_of_parent
  
  x_output <- x_max * 0.5
  x_output
}


define_y <- function(
  tile_number,
  tile_margin = 0.15
) {
  y_max <- 1 - ((tile_number-1)*tile_margin)
  
  y_output <- y_max * 0.5
  y_output
}




plot <- ggplot() +
  geom_rect(
    aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1),
    fill = "#005595"
  ) + theme_void()








input_df <- data.frame(
  TEXT_STRING = c(
    " with email contacts",
    " submitted survey",
    " answered optional Q's"
  ),
  AFFIRM_RESPONSES = c(nemails,ncomplete,nextra),
  PARENT_RESPONSES = c(nbig,nemails,ncomplete),
  SURVEY_RESPONSES = nrow(overlap),
  stringsAsFactors = FALSE
)


input_df$TILE_NUMBER <- seq(nrow(input_df))


input_df <- input_df %>%
  dplyr::mutate(
    X_POSITION = define_x(SURVEY_RESPONSES, PARENT_RESPONSES, AFFIRM_RESPONSES),
    X_WIDTH = X_POSITION*2,
    Y_POSITION = define_y(TILE_NUMBER),
    Y_HEIGHT = Y_POSITION*2,
    PERCENT = paste0(round((AFFIRM_RESPONSES/PARENT_RESPONSES)*100), "%"),
    ORG_COUNT = paste0("(",AFFIRM_RESPONSES, " out of ", PARENT_RESPONSES,"", ")")
  )




input_df$PERCENT <- ifelse(
  as.numeric(str_remove(input_df$PERCENT,"%")) < 10,
  paste0("  ",input_df$PERCENT),
  input_df$PERCENT
)


input_df$TEXT_STRING <- paste0(input_df$TEXT_STRING," ",input_df$ORG_COUNT)


input_df$Y_HEIGHT <- sqrt(input_df$AFFIRM_RESPONSES/input_df$SURVEY_RESPONSES)


input_df$X_WIDTH <- sqrt(input_df$AFFIRM_RESPONSES/input_df$SURVEY_RESPONSES)
input_df$X_POSITION <- input_df$X_WIDTH/2
input_df$Y_POSITION <- input_df$Y_HEIGHT/2




for (n in 1:nrow(input_df)) {
  
  input_row <- input_df[n,]
  
  plot <- plot +
    geom_tile(
      data = input_row,
      aes(x = X_POSITION, y = Y_POSITION),
      width = input_row$X_WIDTH, height = input_row$Y_HEIGHT,
      fill = "#EEE0FF", alpha = 0.12*input_row$TILE_NUMBER
    )
}
plot


plot <- plot +
  annotate(
    "text", x = 0.02, y = input_df$Y_HEIGHT - 0.025, label = input_df$PERCENT,
    size = 7.5*1.15, fontface = "bold", color = "white", hjust = 0
  ) +
  annotate(
    "text", x = 0.0784, y = input_df$Y_HEIGHT - 0.03, label = input_df$TEXT_STRING,
    size = 5*1.15, fontface = "bold", color = "white", hjust = 0
  ) +
  annotate(
    "text", x = 0.02, y = .965, label = nbig,
    size = 8*1.15, fontface = "bold", color = "white", hjust = 0
  ) +
  annotate(
    "text", x = 0.089, y = .9615, label = " healthcare organizations",
    size = 6*1.15, fontface = "bold", color = "white", hjust = 0
  )
plot


ggsave(paste0(respath,"Images/nested_response.png"),plot,height=8,width=12.45,units="in")


results$day <- results$date_modified %>% as_date()
results$weekday <- results$day %>% weekdays()
results$weekday <- ifelse(results$weekday=="Tuesday","Tuesday","Not Tuesday")
ggplot(results) +
  geom_histogram(aes(x=day,fill=weekday),bins=24)


##################################################
#### response venn diagram provider type #########
##################################################






resps <- results %>% dplyr::filter(response_status=="completed")
resps$type <- paste0(resps$PH,resps$BH,resps$OH)
summary(factor(resps$type))


resps$OrgName <- resps$`OHA Mastername(s)`


x <- list(
  "Physical" = resps$OrgName[which(resps$PH==1)],
  "          Behavioral" = resps$OrgName[which(resps$BH==1)], 
  "Oral   " = resps$OrgName[which(resps$OH==1)]
)


venn <- x %>% Venn() %>% process_data()


venn_region(venn)






vcolors <- c("Physical"="#005595",
             "Behavioral"= "#ec8902",
             "Oral" = "#A01c3f" ,
             "Physical..Behavioral"="#6a6c53",
             "Physical..Oral"="#50386a",
             "Behavioral..Oral"="#ca581d",
             "Physical..Behavioral..Oral"="#6f5753")


vcolors <- unname(vcolors)


vennd <- ggplot() +
  geom_sf(aes(), data = venn_region(venn),fill=vcolors,color=NA) +
  scale_fill_manual(values=vcolors) +
  geom_sf_text(aes(label = name), data = venn_setlabel(venn),color="#646464",fontface="bold",size=7) +
  geom_sf_label(aes(label=count), fontface = "bold",
                data = venn_region(venn),fill="white", size=7,
                alpha=0,color="white",label.size=0) +
  theme_void()

vennd


ggsave(paste0(respath,"Images/response_venn.png"),vennd,height=5,width=5,units="in")
#ggsave(paste0(respath,"Images/response_venn.pdf"),vennd,height=5,width=5,units="in")




###############################################
############## EHR satisfaction ###############
###############################################










satPH <- results %>% filter(PH==1 & !is.na(`Overall, how satisfied is your organization with its EHR?`))
satPH <- count(satPH,`Overall, how satisfied is your organization with its EHR?`,.drop=F)
nPH <- sum(satPH$n)
satPH <- satPH %>% mutate(type="Physical",total=nPH)




satBH <- results %>% filter(BH==1 & !is.na(`Overall, how satisfied is your organization with its EHR?`))
satBH <- count(satBH,`Overall, how satisfied is your organization with its EHR?`,.drop=F)
nBH <- sum(satBH$n)
satBH <- satBH %>% mutate(type="Behavioral",total=nBH)


satOH <- results %>% filter(OH==1 & !is.na(`Overall, how satisfied is your organization with its EHR?`))
satOH <- count(satOH,`Overall, how satisfied is your organization with its EHR?`,.drop=F)
nOH <- sum(satOH$n)
satOH <- satOH %>% mutate(type="Oral",total=nOH)


sat <- rbind(satPH,satBH,satOH)


names(sat)[1] <- "Answer"


sat <- sat %>% filter(!is.na(Answer)) %>% select(Answer,type, n,total)
sat$proportion <- sat$n/sat$total


multiples_data <- sat %>%
  dplyr::mutate(type_f = factor(type, levels = c("Physical", "Behavioral","Oral")),
                proportion = round(100*proportion),
                label_position = dplyr::if_else(proportion < 17, proportion+1, proportion-21),
                label_color = "white")


colors <- data.frame(
  facet_values = c("Physical", "Behavioral","Oral"),
  colors = c("#005595", "#ec8902","#A01c3f")
)


## Merge colors together with the original data, readying us for the graph.
input_data <- merge(multiples_data,
                    colors,
                    all.x = TRUE,
                    by.x = "type",
                    by.y = "facet_values")


input_data$label_color <- ifelse(input_data$proportion < 17,input_data$colors,input_data$label_color)


axis_order <- multiples_data %>%
  dplyr::filter(type == "Physical") %>%
  dplyr::select(Answer)


satplot <- ggplot(
  data = input_data,
  aes(
    x = Answer,
    y = proportion
  )
) +
  geom_bar(
    aes(
      fill = colors
    ),
    stat = "identity",
    position = "identity",
    width = 3/4,
    show.legend = F
  ) +
  scale_fill_manual(
    values = c(
      "#005595" = "#005595",
      "#ec8902" = "#ec8902",
      "#A01c3f" = "#A01c3f"
    )
  ) +
  facet_wrap(~type_f,strip.position = "bottom") +
  coord_flip() +
  scale_x_discrete(
    limits = rev(axis_order$Answer)
  ) +
  ## Stretch the bar graphs to align completely with strip text
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    # strip.text = element_text(
    #   face = "bold",
    #   hjust = 0,
    #   size = 12
    # ),
    strip.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(
      face = "bold",
      size = 13,
      color = "#646464",
      margin = margin(
        r = 10
      )
    ),
    axis.text.x = element_blank(),
    panel.spacing = unit(3, "lines"),
    axis.line.x = element_blank(),
    axis.line.y = element_line()
  ) +
  labs(x=NULL,y=NULL,title="") +
  geom_text(
    aes(
      group = type_f,
      label = paste0(proportion, "%"),
      color = label_color,
      y = label_position,
      hjust = 0
    ),
    family="sans",
    fontface="bold",
    size=4.2,
    position = position_nudge(
      x = 0,
      y = 2
    )
  ) +
  scale_color_manual(
    values = c(
      "white" = "white",
      "#005595" = "#005595",
      "#ec8902" = "#ec8902",
      "#A01c3f" = "#A01c3f"
    )
  )


satplot


#ggsave(paste0(respath,"Images/satbars.pdf"),satplot,height=5,width = 9.8,units="in")


ggsave(paste0(respath,"Images/satbars.png"),satplot,height=3,width = 8.8,units="in")




###############################################
######### response times ######################
###############################################


timePH <- results %>% filter(PH==1) %>% group_by(`Thank you...additional Qs okay?`) %>%
  summarise(medtime=median(minutes_spent,na.rm=T),n=n())
timePH$type <- "Physical"
timeOH <- results %>% filter(OH==1) %>% group_by(`Thank you...additional Qs okay?`) %>%
  summarise(medtime=median(minutes_spent,na.rm=T),n=n())
timeOH$type <- "Oral"
timeBH <- results %>% filter(BH==1) %>% group_by(`Thank you...additional Qs okay?`) %>%
  summarise(medtime=median(minutes_spent,na.rm=T),n=n())
timeBH$type <- "Behavioral"


rtime <- rbind(timePH,timeBH,timeOH)
rtime <- rtime %>% filter(!is.na(`Thank you...additional Qs okay?`))


rtime$medminutes <- floor(rtime$medtime)
rtime$medsecounds <- round((rtime$medtime - rtime$medminutes)*60)


rtime$timelabel <- round(rtime$medtime)
rtime$`Thank you...additional Qs okay?` <- as.character(rtime$`Thank you...additional Qs okay?`)




plot <- ggplot(data = rtime, aes(x = medtime, y = type, color = `Thank you...additional Qs okay?`)) +
  geom_point(size = 12, shape = 19) +
  theme_bw() +
  theme(
    legend.position = "none"
  ) +
  labs(x = "", y = "") +
  # scale_y_discrete(limits = type) +
  scale_x_continuous(limits = c(0,25), breaks = seq(0,25,25)) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_color_manual(values = c(Yes = "#536D60", No = "#7DA830")) +
  annotate("text", x = rtime$medtime, y = rtime$type, label = rtime$timelabel, size = 4, fontface = "bold", color = "white") +
  
  
  
  
  
  annotate("text", x = min(rtime$medtime[rtime$type=="Physical"])-3.2,
           y = "Physical", label = "Physical", size = 6, fontface = "bold", color = "#005595") +
  annotate("text", x = min(rtime$medtime[rtime$type=="Behavioral"])-3.6,
           y = "Behavioral", label = "Behavioral", size = 6, fontface = "bold", color = "#ec8902") +
  annotate("text", x = min(rtime$medtime[rtime$type=="Oral"])-2,
           y = "Oral", label = "Oral", size = 6, fontface = "bold", color = "#A01c3f") +
  
  theme(
    panel.border = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#EEEEEE", size = 0.5),
    axis.line.x = element_line(color = "#999999", size = 0.5),
    axis.line.y = element_blank(),
    axis.text = element_text(face="bold",size=12)
  )


plot


t1 <- textGrob(expression("Median completion time " *
                            phantom(bold("without")) * " and " *
                            phantom(bold("with")) *
                            " optional questions"),
               x = unit(0.5, "lines"), y = unit(-0.5, "lines"), hjust = 0, vjust = 0.3,  gp = gpar(col = "#646464",fontsize = 17))


t2 <- textGrob(expression(phantom("Median completion time ") *
                            bold("without") * phantom(" and ") *
                            phantom(bold("with")) *
                            phantom(" optional questions")),
               x = unit(0.5, "lines"), y = unit(-0.5, "lines"), hjust = 0, vjust = 0.3,  gp = gpar(col = "#7DA830",fontsize = 17))


t3 <- textGrob(expression(phantom("Median completion time ") *
                            phantom(bold("without")) * phantom(" and ") *
                            bold("with") *
                            phantom(" optional questions")),
               x = unit(0.5, "lines"), y = unit(-0.5, "lines"), hjust = 0, vjust = 0.3,  gp = gpar(col = "#536D60",fontsize = 17))






t <- grobTree(t1,t2,t3)


# b1 <- textGrob(paste0("N Physical = ", nPH2,", N Behavioral = ", nBH2,", N Oral = ",nOH2),
#                x = unit(0.5, "lines"), y = unit(0.5, "lines"), hjust = -0, vjust = -0,  gp = gpar(col = "#646464",fontsize = 11))


gg <- arrangeGrob(plot,top=t)
plot2 <- grid.arrange(gg)

plot


ggsave(paste0(respath,"Images/response_time.png"),plot2,height=5,width=7,units="in")




{
  
  rtime$`Median time` <- paste0(rtime$medminutes,"m ",rtime$medsecounds,"s")
  
  rtime <- rtime %>% select(type,`Median time`,`Thank you...additional Qs okay?`) %>%
    pivot_wider(names_from="type",values_from = "Median time")
  
  
  rtime <- rtime %>% arrange(`Thank you...additional Qs okay?`)
  
  rtime$`Thank you...additional Qs okay?` <- c("Completed first section","Completed both sections")
  rtime %>% gt() %>%
    tab_header("Median completion time") %>%
    cols_label(
      `Thank you...additional Qs okay?` = "",
      Physical = html("<center><strong><font color=#005595>Physical</font></strong></center>"),
      Behavioral = html("<center><strong><font color=#ec8902>Behavioral</font></strong></center>"),
      Oral = html("<center><strong><font color=#A01c3f>Oral</font></strong></center>")
    ) %>%
    opt_table_font(font="Arial") %>%
    tab_style(
      style = list(
        cell_fill(color = "#e6f4ff"),
        cell_text(weight = "bold")
      ),
      locations=cells_body(
        columns=Physical
      )) %>%
    tab_style(
      style = list(
        cell_fill(color = "#fff4e6"),
        cell_text(weight = "bold")
      ),
      locations=cells_body(
        columns=Behavioral
      )) %>%
    tab_style(
      style = list(
        cell_fill(color = "#fbe9ee"),
        cell_text(weight = "bold")
      ),
      locations=cells_body(
        columns=Oral
      )) %>%
    cols_align(align="center",columns=everything())
  
  
}


##################################################
############# HIE tool access ####################
##################################################
{
  hie <- results %>%  select(72:85,PH,OH,BH)
  hie$missings <- apply(hie,1,function (x) {sum(is.na(x))})
  hie <- hie %>% filter(missings<14)
  names(hie) <- names(hie) %>%
    str_remove_all("Does your organization have access to the following health information exchange \\(HIE\\) tools\\? - ") %>%
    str_remove_all(\\(.{0,1000}\\)$) %>% trimws()
  
  
  hiePH <- hie %>% filter(PH==1) %>% select(-Clara,-missings,-PH,-BH,-OH,-`Other HIE Tool`)
  nPH <- nrow(hiePH)
  PHnone <- sum(apply(hiePH,1,function (x) {sum(x!="Yes" | is.na(x))})==12)
  
  hiePH <- t(hiePH) %>% as.data.frame()
  hiePH$n <- apply(hiePH,1,function (x) {sum(x=="Yes",na.rm=T)})
  hiePH$proportion <- hiePH$n/nPH
  hiePH <- hiePH %>% select(n,proportion) %>% rownames_to_column(var="Answer") %>%
    mutate(type="Physical")
  
  
  
  
  hiePH <- rbind(hiePH,c("None",PHnone,PHnone/nPH,"Physical"))
  hiePH <- rbind(hiePH,c("At least one",nPH-PHnone,(nPH-PHnone)/nPH,"Physical"))
  
  
  
  
  hieOH <- hie %>% filter(OH==1) %>% select(-Clara,-missings,-PH,-BH,-OH,-`Other HIE Tool`)
  nOH <- nrow(hieOH)
  OHnone <- sum(apply(hieOH,1,function (x) {sum(x!="Yes" | is.na(x))})==12)
  
  
  hieOH <- t(hieOH) %>% as.data.frame()
  hieOH$n <- apply(hieOH,1,function (x) {sum(x=="Yes",na.rm=T)})
  hieOH$proportion <- hieOH$n/nOH
  hieOH <- hieOH %>% select(n,proportion) %>% rownames_to_column(var="Answer") %>%
    mutate(type="Oral")
  
  
  hieOH <- rbind(hieOH,c("None",OHnone,OHnone/nOH,"Oral"))
  hieOH <- rbind(hieOH,c("At least one",nOH-OHnone,(nOH-OHnone)/nOH,"Oral"))
  
  
  hieBH <- hie %>% filter(BH==1) %>% select(-Clara,-missings,-PH,-BH,-OH,-`Other HIE Tool`)
  nBH <- nrow(hieBH)
  BHnone <- sum(apply(hieBH,1,function (x) {sum(x!="Yes" | is.na(x))})==12)
  
  
  hieBH <- t(hieBH) %>% as.data.frame()
  hieBH$n <- apply(hieBH,1,function (x) {sum(x=="Yes",na.rm=T)})
  hieBH$proportion <- hieBH$n/nBH
  hieBH <- hieBH %>% select(n,proportion) %>% rownames_to_column(var="Answer") %>%
    mutate(type="Behavioral")
  
  
  hieBH <- rbind(hieBH,c("None",BHnone,BHnone/nBH,"Behavioral"))
  hieBH <- rbind(hieBH,c("At least one",nBH-BHnone,(nBH-BHnone)/nBH,"Behavioral"))
  
  
  hie <- rbind(hiePH,hieOH,hieBH)
  hie$proportion <- as.numeric(hie$proportion)
  
  
  hie <- hie %>% dplyr::filter(Answer!="None")
  
  
  multiples_data <- hie %>%
    dplyr::mutate(type_f = factor(type, levels = c("Physical", "Behavioral","Oral")),
                  proportion = round(100*proportion),
                  label_position = dplyr::if_else(proportion <= 13, proportion, proportion-14),
                  label_color = "white")
  
  
  #multiples_data$label_position <- ifelse(multiples_data$proportion>6 & multiples_data$proportion<10,multiples_data$label_position-6,multiples_data$label_position)
  #multiples_data$label_position <- ifelse(multiples_data$proportion==0,multiples_data$proportion,multiples_data$label_position)
  
  
  axis_order <- multiples_data %>%
    dplyr::filter(type == "Physical") %>%
    dplyr::arrange(proportion) %>%
    dplyr::select(Answer)
  
  
  
  
  colors <- data.frame(
    facet_values = c("Physical", "Behavioral","Oral"),
    colors = c("#005595", "#ec8902","#A01c3f")
  )
  
  
  ## Merge colors together with the original data, readying us for the graph.
  input_data <- merge(multiples_data,
                      colors,
                      all.x = TRUE,
                      by.x = "type",
                      by.y = "facet_values")
  
  
  input_data$label_color <- ifelse(input_data$proportion <= 13,input_data$colors,input_data$label_color)
  
  
  axis_order <- multiples_data %>%
    dplyr::filter(type == "Physical") %>%
    arrange(-proportion) %>%
    dplyr::select(Answer)
  
  
  
  
  hieplot <- ggplot(
    data = input_data,
    aes(
      x = Answer,
      y = proportion
    )
  ) +
    geom_bar(
      aes(
        fill = colors
      ),
      stat = "identity",
      position = "identity",
      width = 3/4,
      show.legend = F
    ) +
    scale_fill_manual(
      values = c(
        "#005595" = "#005595",
        "#ec8902" = "#ec8902",
        "#A01c3f" = "#A01c3f"
      )
    ) +
    facet_wrap(~type_f,strip.position = "bottom") +
    coord_flip() +
    scale_x_discrete(
      limits = rev(axis_order$Answer)
    ) +
    ## Stretch the bar graphs to align completely with strip text
    scale_y_continuous(
      expand = c(0, 0),
      limits=c(0,65)
    ) +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      # strip.text = element_text(
      #   face = "bold",
      #   hjust = 0,
      #   size = 12
      # ),
      strip.text = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.y = element_text(
        face = "bold",
        size = 13,
        color = "#646464",
        margin = margin(
          r = 10
        )
      ),
      axis.text.x = element_blank(),
      panel.spacing = unit(3, "lines"),
      axis.line.x = element_blank(),
      axis.line.y = element_line()
    ) +
    labs(x=NULL,y=NULL,title="") +
    geom_text(
      aes(
        group = type_f,
        label = paste0(proportion, "%"),
        color = label_color,
        y = label_position,
        hjust = 0
      ),
      family="sans",
      fontface="bold",
      size=4.2,
      position = position_nudge(
        x = 0,
        y = 2
      )
    ) +
    scale_color_manual(
      values = c(
        "white" = "white",
        "#005595" = "#005595",
        "#ec8902" = "#ec8902",
        "#A01c3f" = "#A01c3f"
      )
    )
}
hieplot


#ggsave(paste0(respath,"Images/satbars.pdf"),satplot,height=5,width = 9.8,units="in")


ggsave(paste0(respath,"Images/HIEbars.png"),hieplot,height=3.5,width = 8.8,units="in")


###############################################
########### patient engagement ###############
###############################################
{
  patient <- results %>%  select(117:123,PH,OH,BH)
  patient$missings <- apply(patient,1,function (x) {sum(is.na(x))})
  patient <- patient %>% filter(missings<7)
  names(patient) <- names(patient) %>%
    str_remove_all("Does your EHR system allow patients to. - ") %>%
    str_remove_all("\\?") %>%
    str_remove_all(\\(.{0,1000}\\)$) %>% trimws()
  
  
  names(patient)[1:5] <- c("Access a patient portal",
                           "Download medical record",
                           "Send records to 3rd party",
                           "Upload device/app info",
                           "Send secure messages")
  
  
  patientPH <- patient %>% filter(PH==1) %>% select(-missings,-PH,-BH,-OH)
  nPH <- nrow(patientPH)
  patientPH <- t(patientPH) %>% as.data.frame()
  patientPH$n <- apply(patientPH,1,function (x) {sum(x=="Yes",na.rm=T)})
  patientPH$proportion <- patientPH$n/nPH
  patientPH <- patientPH %>% select(n,proportion) %>% rownames_to_column(var="Answer") %>%
    mutate(type="Physical")
  
  
  patientOH <- patient %>% filter(OH==1) %>% select(-missings,-PH,-BH,-OH)
  nOH <- nrow(patientOH)
  patientOH <- t(patientOH) %>% as.data.frame()
  patientOH$n <- apply(patientOH,1,function (x) {sum(x=="Yes",na.rm=T)})
  patientOH$proportion <- patientOH$n/nOH
  patientOH <- patientOH %>% select(n,proportion) %>% rownames_to_column(var="Answer") %>%
    mutate(type="Oral")
  
  
  patientBH <- patient %>% filter(BH==1) %>% select(-missings,-PH,-BH,-OH)
  nBH <- nrow(patientBH)
  patientBH <- t(patientBH) %>% as.data.frame()
  patientBH$n <- apply(patientBH,1,function (x) {sum(x=="Yes",na.rm=T)})
  patientBH$proportion <- patientBH$n/nBH
  patientBH <- patientBH %>% select(n,proportion) %>% rownames_to_column(var="Answer") %>%
    mutate(type="Behavioral")
  
  
  patient <- rbind(patientPH,patientOH,patientBH)
  
  
  multiples_data <- patient %>%
    dplyr::mutate(type_f = factor(type, levels = c("Physical", "Behavioral","Oral")),
                  proportion = round(100*proportion),
                  label_position = dplyr::if_else(proportion < 10, proportion-1, proportion-18),
                  label_color = "white")
  
  
  multiples_data$label_position <- ifelse(multiples_data$proportion>6 & multiples_data$proportion<10,multiples_data$label_position-5,multiples_data$label_position)
  
  
  axis_order <- multiples_data %>%
    dplyr::filter(type == "Physical") %>%
    dplyr::arrange(proportion) %>%
    dplyr::select(Answer)
  
  
  
  
  colors <- data.frame(
    facet_values = c("Physical", "Behavioral","Oral"),
    colors = c("#005595", "#ec8902","#A01c3f")
  )
  
  
  ## Merge colors together with the original data, readying us for the graph.
  input_data <- merge(multiples_data,
                      colors,
                      all.x = TRUE,
                      by.x = "type",
                      by.y = "facet_values")
  
  
  input_data$label_color <- ifelse(input_data$proportion <= 6,input_data$colors,input_data$label_color)
  
  
  axis_order <- multiples_data %>%
    dplyr::filter(type == "Physical") %>%
    arrange(-proportion) %>%
    dplyr::select(Answer)
  
  
  
  
  patientplot <- ggplot(
    data = input_data,
    aes(
      x = Answer,
      y = proportion
    )
  ) +
    geom_bar(
      aes(
        fill = colors
      ),
      stat = "identity",
      position = "identity",
      width = 3/4,
      show.legend = F
    ) +
    scale_fill_manual(
      values = c(
        "#005595" = "#005595",
        "#ec8902" = "#ec8902",
        "#A01c3f" = "#A01c3f"
      )
    ) +
    facet_wrap(~type_f,strip.position = "bottom") +
    coord_flip() +
    scale_x_discrete(
      limits = rev(axis_order$Answer)
    ) +
    ## Stretch the bar graphs to align completely with strip text
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      # strip.text = element_text(
      #   face = "bold",
      #   hjust = 0,
      #   size = 12
      # ),
      strip.text = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.y = element_text(
        face = "bold",
        size = 13,
        color = "#646464",
        margin = margin(
          r = 10
        )
      ),
      axis.text.x = element_blank(),
      panel.spacing = unit(3, "lines"),
      axis.line.x = element_blank(),
      axis.line.y = element_line()
    ) +
    labs(x=NULL,y=NULL,title="") +
    geom_text(
      aes(
        group = type_f,
        label = paste0(proportion, "%"),
        color = label_color,
        y = label_position,
        hjust = 0
      ),
      family="sans",
      fontface="bold",
      size=4.2,
      position = position_nudge(
        x = 0,
        y = 2
      )
    ) +
    scale_color_manual(
      values = c(
        "white" = "white",
        "#005595" = "#005595",
        "#ec8902" = "#ec8902",
        "#A01c3f" = "#A01c3f"
      )
    )
}
patientplot


#ggsave(paste0(respath,"Images/satbars.pdf"),satplot,height=5,width = 9.8,units="in")


ggsave(paste0(respath,"Images/patientbars.png"),patientplot,height=2.7,width = 10,units="in")


################################################
######### make a pretty workbook ###############
################################################


names(results)[which(str_length(names(results))>255)] <- gsub(\\s*\\([^\\)]+\\),"",names(results)[which(str_length(names(results))>255)])


ccos <- unique(results$CCO) %>% sort()
wb <- createWorkbook()
addWorksheet(wb,"Overview")


for (i in 1:length(ccos)) {
  addWorksheet(wb,ccos[i])
}


writeDataTable(wb,1,cover,startRow=1,startCol=1,tableStyle="TableStyleMedium1",withFilter = F)


for (i in 1:length(ccos)) {
  writeDataTable(wb,sheet=i+1,results[which(results$CCO==ccos[i]),],startRow=1,startCol=1,
                 tableStyle="TableStyleMedium1",withFilter=F)
}


style1 <- createStyle(wrapText = TRUE,valign="center",halign="center")
style2 <- createStyle(wrapText = F,valign="center",halign="center")


addStyle(wb, sheet = 1, style1, rows = 1, cols = 1:ncol(cover), gridExpand = TRUE)
addStyle(wb, sheet = 1, style2, rows = 2:(nrow(cover)+1), cols = 1:ncol(cover), gridExpand = TRUE)


setRowHeights(wb, sheet=1, rows=1, heights=35)
setColWidths(wb,sheet=1,cols=1:ncol(cover),widths=15)


insertImage(wb,sheet=1,file=paste0(respath,"Images/EHRbars.png"),width=8.8,height=3,units="in",startRow = nrow(cover)+3,startCol = 1)
insertImage(wb,sheet=1,file=paste0(respath,"Images/completeness.png"),width=3,height=3,units="in",startRow = nrow(cover)+3,startCol = 7)
insertImage(wb,sheet=1,file=paste0(respath,"Images/benefits.png"),width=12*.8,height=5*.8,units="in",startRow = 1,startCol = 11)


openXL(wb)




results$`Health equity comments?` %>% unique() %>% cat(sep="\n\n")
results$`Is there anything else you'd like to share?` %>% unique() %>% cat(sep="\n\n")