#BST260 final project - 2020FALL
#Date:12/13/2020
#Group 9
# COVID-19: US County-Level Analyses & Predictions
# Zhaoxun Hou
# Zechen Liu
# Jiayue Guan
# Haoyuan Li

### Note: in my environment, it takes about 3-4 minutes to run this rhiny app, 
### feel free to grab a cup of coffee :)

#setwd("D:/OneDrive - Harvard University/Desktop/files/material/BST260/BST260final project")

#packages used in Tab1: 
library(tidyverse)
library(lubridate)
library(shinythemes)
library(shiny)
library(RColorBrewer)
library(maps)
library(usmap)
library(zoo) #rolling mean
library(sp)   #spatial nearest neighbour
library(rgeos)#spatial nearest neighbour

#packages used in Tab2: GLM
library(caret)
library(MASS)

#packages used in Tab3: Instantaneous Rt
library(EpiDynamics)
library(EpiEstim) 
library(R0)
library(distcrete)
library(epitrix)
library(projections)
library(magrittr)
library(incidence)
library(splitstackshape)
###MASS and dplyr both have select
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

######################
# Covid-19-dashboard #
######################
us_counties <- read_csv("us-counties.csv")
us_counties[us_counties$fips=="35013",]$county=="dona ana" #Do??a Ana

allcounty <- map_data("county") %>%
    mutate(subregion=str_remove_all(subregion,pattern="\\s"))
allcounty[allcounty$region=="south dakota"&allcounty$subregion=="shannon",]$subregion<-"Oglala Lakota"


us_counties<-us_counties %>%
    mutate(region=tolower(state),
           subregion=tolower(county),
           subregion=str_remove_all(subregion,pattern="\\s"),# counties named "st. xxx"
           subregion=str_remove_all(subregion,pattern="\\.")) %>% # counties named "st. xxx"
    group_by(county,state) %>%
    arrange(date) %>%
    mutate(newcases=cases-lag(cases),
           newcases=ifelse(newcases<0,0,newcases),
           newcases7dayavg=rollmean(newcases, k = 7, fill = NA),
           newdeaths=deaths-lag(deaths),
           newdeaths=ifelse(newdeaths<0,0,newdeaths),
           newdeaths7dayavg=rollmean(newdeaths, k = 7, fill = NA)) %>%
    ungroup()
#######################################
#computing 6 spatial nearest neighbour#
#######################################
time_series_covid19_confirmed_US <- read_csv("time_series_covid19_confirmed_US.csv")
time_series_covid19_confirmed_US[which(time_series_covid19_confirmed_US$FIPS=="36061"),]$Admin2 <-"New York City"

    
spatialdata <- time_series_covid19_confirmed_US %>%
    mutate(region=tolower(Province_State),
           subregion=tolower(Admin2),
           subregion=str_remove_all(subregion,pattern="\\."),#remove . from place named "st. xxx"
           subregion=str_remove_all(subregion,pattern="\\s"),
           subregion=str_remove_all(subregion,pattern="\\'")) %>%
    select(region,subregion,Lat,Long_)
sp.spatialdata<-spatialdata
coordinates(sp.spatialdata) <- ~Long_+Lat
d <- gDistance(sp.spatialdata, byid=T)
#min0.d is 1st nearest county: itself
#minx.d is xth nearest county
min0.d <- apply(d, 1, function(x) order(x, decreasing=F)[1])
min1.d <- apply(d, 1, function(x) order(x, decreasing=F)[2])
min2.d <- apply(d, 1, function(x) order(x, decreasing=F)[3])
min3.d <- apply(d, 1, function(x) order(x, decreasing=F)[4])
min4.d <- apply(d, 1, function(x) order(x, decreasing=F)[5])
min5.d <- apply(d, 1, function(x) order(x, decreasing=F)[6])
min6.d <- apply(d, 1, function(x) order(x, decreasing=F)[7])

#
newdata <- cbind(spatialdata,
                 spatialdata[min0.d,c("region","subregion")],
                 spatialdata[min1.d,c("region","subregion")],
                 spatialdata[min2.d,c("region","subregion")],
                 spatialdata[min3.d,c("region","subregion")],
                 spatialdata[min4.d,c("region","subregion")],
                 spatialdata[min5.d,c("region","subregion")])
colnames(newdata) <- c(colnames(spatialdata),
                       "neighbor0nregion", "neighbor0nsubregion",
                       "neighbor1nregion", "neighbor1nsubregion",
                       "neighbor2nregion", "neighbor2nsubregion",
                       "neighbor3nregion", "neighbor3nsubregion",
                       "neighbor4nregion", "neighbor4nsubregion",
                       "neighbor5nregion", "neighbor5nsubregion")
# convering the data format from wide data to long data,
# then to wide data with "n","region","subregion"
spatialneighbordata<-newdata %>%
    gather(key, value,"neighbor0nregion":"neighbor5nsubregion") %>%
    tidyr::extract(key, c("n", "regionsubregion"), "(\\d)([a-z]+)")%>%
    spread("regionsubregion", value)

#joining captialized names of state and county into "spatialneighbordata"
spatialneighbordata<-us_counties %>%
    mutate(placeforselect=str_c(county,state,sep=", ")) %>%
    select(placeforselect,county,state,region,subregion) %>%
    distinct(placeforselect,county,state,region,subregion) %>%
    left_join(spatialneighbordata,by=c("region","subregion"))

#################################################



#########################
#GLM                    #
#########################

#datasets used
apple <- read.csv("apple.csv")

apple[which(apple$countyFIPS=="36061"),]$Region<-"new york city" # "new york city and new york county"

infection <- read.csv("us-counties.csv") %>% 
    mutate(placeforselect=paste(county,state, sep=", ")) 


########################

###########################
# Instantaneous Rt        #
###########################
counties <- read.csv("us-counties.csv")
counties$date<-ymd(counties$date)
counties<-counties%>%
    mutate(key=paste(county,state, sep=", "))%>%
    mutate(date=ymd(counties$date))%>%
    group_by(key)%>%
    arrange(date)%>%
    mutate(new_case=cases-lag(cases,default = 0))%>%
    mutate(new_death=deaths-lag(deaths,default = 0))%>%
    mutate(average=rollmean(new_case, k = 7, fill = NA))
mu <- 5.49
sigma <- 2.54
cv <- sigma / mu
params <- gamma_mucv2shapescale(mu, cv)
si <- distcrete("gamma", shape = params$shape,
                scale = params$scale,
                interval = 1L, w = 0.5)
si<-c(0,si$d(1:20))




#############
# rnn model #
#############

#read model results
#rnn model used pre-runned results, which have been writen into "RNN_result.csv"
rnn<-read_csv("rnn_results.csv")
# colnames(rnn) <- c("fips","2020-11-19","2020-11-20",
#                    "2020-11-21","2020-11-22","2020-11-23","2020-11-24",
#                    "2020-11-25","2020-11-26","2020-11-27","2020-11-28",
#                    "2020-11-29","2020-11-30","2020-12-01","2020-12-02")
rnn <- rnn %>% 
    select(-name) %>% 
    gather("date","pred_newcases",-fips)
rnn_fips<-rnn %>% select(fips) %>% distinct()
rnn <- infection %>%
    left_join(rnn,c("fips","date")) %>% 
    right_join(rnn_fips,"fips") %>% 
    mutate(placeforselect = paste(county,state, sep=", ")) %>%
    select(-c(county,state,fips,deaths)) %>% 
    group_by(placeforselect) %>% 
    arrange(date) %>% 
    mutate(date=ymd(date),
           newcases=cases-lag(cases),
           newcases=ifelse(newcases<0,0,newcases),
           newcases7dayavg=rollmean(newcases, k = 7, fill = NA))
rm(rnn_fips)
########
###ui###
########
ui <- fluidPage(theme = shinytheme("yeti"),
                
                tabsetPanel(
                    #################
                    #First tab
                    tabPanel(
                        # Application title
                        titlePanel("US COVID-19 Dashboard"),

                        sidebarLayout(
                            sidebarPanel(
                                # radioButtons("tab1_choice",
                                #              "How would you like to select the counites:",
                                #              choices = c("5 spatial neasest neighbor" = "choice1",
                                #                          "Manually select up to 6 counties"="choice2")),

                                selectizeInput("tab1_county","Select 1 county to compare with 5 spatially nearest neighbors",
                                               choices=c(distinct(spatialneighbordata,placeforselect)),
                                               selected=c("Middlesex, Massachusetts")),

                                radioButtons("tab1_var","Select which variable you want to display:",
                                             choices = c("Cases" = "cases",
                                                         "New cases" = "newcases",
                                                         "New cases (7 day average)"="newcases7dayavg",
                                                         "Deaths" = "deaths",
                                                         "New deaths" = "newdeaths",
                                                         "New deaths (7 day average)"="newdeaths7dayavg")),

                                dateInput("tab1_date",
                                          "Select the date:",
                                          value="2020-11-29",
                                          min="2020-01-24",
                                          max="2020-11-29"),

                                textOutput("tab1_tabletitle"),
                                tableOutput("tab1_table")

                            ), #sidebarPanel END

                            # Show a plot of the generated distribution
                            mainPanel(
                                plotOutput("tab1_lineplot"),
                                plotOutput("tab1_map")
                            )
                        )


                    ),# End of First tab
                    ##########



                    # Second tab
                    tabPanel(
                        # Application title
                        titlePanel("GLM"),
                        sidebarLayout(
                            sidebarPanel(
                                selectizeInput("tab2_county",label = "Select 1 county",
                                               choices=c(distinct(infection,placeforselect)),
                                               selected=c("Middlesex, Massachusetts")),
                                p("We used linear regression and negative binomial distribution to try to predict new cases and seven-day roll mean of new cases per day, with predictors: "),br(),
                                p("1.Direction searches trends in Apple Map to characterize social distancing pattern;) "),br(),
                                p("2.New cases or seven-day roll mean of new cases per day for the past five days.) "),br(),
                                p("Here, you can type in a specific county and obtain a plot of true cases (as red dots) and predicted cases (as orange lines) against dates (all from a specific test set instead of the whole dataset). Moreover, the MSE of our predicted values against the true values is given with our plots."), br(),
                                p("Notes: "),br(),
                                p("1.Dates for different counties may differ in length because of our resource data."), br(),
                                p("2.Negative values of cases per day were set to zero before regression."),br(),
                                p("3.Some counties may return no outcome because they have too few cases across the whole year. (e.g., Walla Walla, Washington)."),
                            ),#End of sidebarPanel


                            mainPanel(
                                 plotOutput("tab2_plot1"),
                                 plotOutput("tab2_plot2")
                            )#End of mainPanel

                        )#End of sidebarLayout

                    ),# End of second tab
                    
                    
                    # Third tab
                    tabPanel(
                        titlePanel("Instantaneous Rt"),
                        sidebarLayout(
                            sidebarPanel(
                                selectInput("tab3_county", label = "Select 1 county",
                                            choices = as.list(unique(counties$key))),
                                p("The plot shows the estimated instantaneous reproduction number (Rt) of an epidemic in every county by the incidence time series and the serial interval distribution. In other word, it shows that how many people one cases can infected on average at a given time." )),
                            mainPanel(
                                plotOutput("tab3_plot1")
                            )
                        )
                    ),
                    
                    # Fourth tab
                    tabPanel(
                        # Application title
                        titlePanel("RNN"),
                        sidebarLayout(
                            sidebarPanel(
                                selectizeInput("tab4_county",label = "Select 1 county",
                                               choices=c(distinct(rnn,placeforselect)),
                                               selected=c("Middlesex, Massachusetts")),
                                p("Single-layer LSTM network was used to predict the new cases in next day give the new cases in previous days. The model is trained by minimizing the mean square error of predicted new cases and actual new cases. After that, trained LSTM was used to predict the new cases from November 16th to November 29th given the cases in previous days."),
                                p("Noted that only counties in testing sets are included here, you may not find your interested counties.")
                            ),#End of sidebarPanel


                            mainPanel(
                                 plotOutput("tab4_plot1"),
                                 plotOutput("tab4_plot2")
                            )#End of mainPanel

                        )#End of sidebarLayout

                    )# End of fourth tab
                    
                    
                    
                    
                    
                ) #END of tabsetPanel
) #END of fluidPage







############
###server###
############
server <- function(input, output) {
##############
    #"tab1_tabletitle" is the interactive title of the following "tab1_table"
    output$tab1_tabletitle<-renderText({
        paste("Information of each selected counties in",input$tab1_date)
    })#END of render text: tabletitle


    #"tab1_table" shows cumulative/new cases/deathes of selected counties on selected date
    output$tab1_table<-renderTable({
        spatialneighbordata %>%
            #filter the selected county
            filter(placeforselect==input$tab1_county) %>%
            #joining dataset to show 6 nearest counties of selected county
            left_join(us_counties,by=c("nregion"="region","nsubregion"="subregion")) %>%
            #delect useless column
            select(placeforselect, date:newdeaths7dayavg) %>%
            #mutating place: convert county name into "County, State" formate for different color
            mutate(place=str_c(county.y,state.y,sep=", ")) %>%
            filter(date==input$tab1_date) %>%
            select(place,
                   cases,newcases,newcases7dayavg,
                   deaths,newdeaths,newdeaths7dayavg) %>%
            rename("County"="place",
                   "Cases"="cases",
                   "New cases"="newcases",
                   "Avg new cases"="newcases7dayavg",
                   "Deaths"="deaths",
                   "New deaths"="newdeaths",
                   "Avg New deaths"="newdeaths7dayavg")
    })#END of render table



    #reactive line plot the year trend of selected counties and 5 nearest counties
    tab1_lineplotselect<-reactive({
        spatialneighbordata %>%
            #filter the selected county
            filter(placeforselect==input$tab1_county) %>%
            #joining dataset to show 6 nearest counties of selected county
            left_join(us_counties,by=c("nregion"="region","nsubregion"="subregion")) %>%
            #delect useless column
            select(placeforselect, date:newdeaths7dayavg) %>%
            #mutating place: convert county name into "County, State" formate for different color
            mutate(place=str_c(county.y,state.y,sep=", ")) %>%
            ggplot(aes_string(x="date",y=input$tab1_var,color="place"))+
            geom_line(alpha=0.9,size=0.5)+
            geom_vline(xintercept=input$tab1_date,color="red",size=1)+
            scale_x_date(limits=c(),date_breaks = "1 month", labels = scales::date_format("%Y-%m"))+
            theme(axis.text.x=element_text(angle=45, hjust=1))+
            theme_bw()+
            theme(legend.position = "top")+
            scale_color_manual(name=NULL,values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                                                  "#A65628","#000000"))
    })

    output$tab1_lineplot <-renderPlot({

        if(input$tab1_var == "cases") {
            tab1_lineplotselect() +
                labs(x="",y="Cases")
        }

        else if(input$tab1_var == "newcases") {
            tab1_lineplotselect() +
                labs(x="",y="New cases")
        }

        else if(input$tab1_var == "newcases7dayavg") {
            tab1_lineplotselect() +
                labs(x="",y="New cases (7 day average)")
        }

        else if(input$tab1_var == "deaths") {
            tab1_lineplotselect() +
                labs(x="",y="Deaths")
        }

        else if(input$tab1_var == "newdeaths") {
            tab1_lineplotselect() +
                labs(x="",y="New deaths")
        }
        else if(input$tab1_var == "newdeaths7dayavg") {
            tab1_lineplotselect() +
                labs(x="",y="New deaths (7 day average)")
        }

    }) #END of output$tab1_lineplot <-renderPlot


    #reactive map for the data selected
    tab1_mapselect <-  reactive({
        us_counties %>%
            filter(date==input$tab1_date) %>%
            full_join(allcounty,by=c("region","subregion")) %>%
            ggplot(aes(x = long, y = lat, group = group))+
            geom_polygon(aes_string(fill=input$tab1_var),color = "gray")+
            theme(panel.grid.major = element_blank(),
                  panel.background = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank())+
            coord_fixed(1.3)
    })

    
    output$tab1_map<-renderPlot({

        if(input$tab1_var == "cases") {
            tab1_mapselect() +
                scale_fill_gradientn(name="Cases",
                                     colors = brewer.pal(9,"Reds"),
                                     trans="log10",
                                     na.value="white"
                                     )
        }

        else if(input$tab1_var == "newcases") {
            tab1_mapselect() +
                scale_fill_gradientn(name="New Cases",
                                     colors = brewer.pal(9,"Reds"),
                                     trans="log10",
                                     na.value="white",
                                     limits = c(1,15000)
                                     )
        }

        else if(input$tab1_var == "newcases7dayavg") {
            tab1_mapselect() +
                scale_fill_gradientn(name="New cases (7 day average)",
                                     colors = brewer.pal(9,"Reds"),
                                     trans="log10",
                                     na.value="white",
                                     limits = c(1,1500)
                )
        }

        else if(input$tab1_var == "deaths") {
            tab1_mapselect() +
                scale_fill_gradientn(name="Deaths",
                                     colors = brewer.pal(9,"Reds"),
                                     trans="log10",
                                     na.value="white",
                                     limits = c(1,15000)
                                     )
        }

        else if(input$tab1_var == "newdeaths") {
                tab1_mapselect() +
                scale_fill_gradientn(name="New Deaths",
                                     colors = brewer.pal(9,"Reds"),
                                     trans="log10",
                                     na.value="white",
                                     limits = c(1,1500)
                                     )
        }

        else if(input$tab1_var == "newdeaths7dayavg") {
            tab1_mapselect() +
                scale_fill_gradientn(name="New deaths (7 day average)",
                                     colors = brewer.pal(9,"Reds"),
                                     trans="log10",
                                     na.value="white",
                                     limits = c(1,1500)
                                     )
        }
    })#End of output$map<-renderPlot
#######################
#    END of TAB 1     #
#######################


###########################
# Function for TAB 2: GLM #
###########################
    #tab2_glmdat is the dataset for GLM
    tab2_glmdat<-reactive({
        #select the data for selected county
        county_code_choose <- unique(infection$fips[which(infection$placeforselect==input$tab2_county)])
        #data wrangling for glm
        infection_choose <- infection %>% filter(fips == county_code_choose) %>% dplyr::select(-deaths)
        apple_choose <- apple %>% filter(countyFIPS == county_code_choose) %>% dplyr::select(-c(X,State.Province,Region,countyFIPS)) %>% spread(Sector,Percent.Change)
        #calculate 7day avg of new cases as predictors
        dat <- right_join(apple_choose,infection_choose,by = c("Date" = "date")) %>%
            mutate(new_cases = cases-lag(cases),
                   new_cases=ifelse(new_cases<0,0,new_cases),
                   one = lag(new_cases,1),
                   two = lag(new_cases,2),
                   three = lag(new_cases,3),
                   four = lag(new_cases,4),
                   five = lag(new_cases,5))
        dat$sdm <- rollmean(dat$new_cases, k = 7, fill = NA)
        dat$sdm_1 <- rollmean(dat$one, k = 7, fill = NA)
        dat$sdm_2 <- rollmean(dat$two, k = 7, fill = NA)
        dat$sdm_3 <- rollmean(dat$three, k = 7, fill = NA)
        dat$sdm_4 <- rollmean(dat$four, k = 7, fill = NA)
        dat$sdm_5 <- rollmean(dat$five, k = 7, fill = NA)
        dat$Date <- as.Date(dat$Date)
        dat <- dat%>% filter(is.na(dat$new_cases)==F)
    })

        #tab2_glm1 make line plot using linear models for new cases
        tab2_glm1<-reactive({

            #data wrangling after a county is chosen
            dat<-tab2_glmdat()

            #glm model
            set.seed(520)
            train_index <- createDataPartition(dat$new_cases, p = 0.7, list = F)
            training <- dat[train_index,]
            testing <- dat[-train_index,]
            #linear model for new cases
            plain_lm <- lm(new_cases ~ driving + sdm_1 + sdm_2 + sdm_3 + sdm_4 + sdm_5, data = training)

            #predictions using the linear models
            pred_plain <- predict(plain_lm, testing)
            mse_plain<-sqrt(mean(((testing$new_cases-(pred_plain)))^2, na.rm = T))
                ggplot(testing)+
                    geom_point(aes(x=Date,y=new_cases),color = "red")+
                    scale_x_date(date_breaks = "1 month",labels = scales::date_format("%Y-%m"))+
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
                    geom_line(aes(x=Date, y=pred_plain),color="darkgoldenrod2")+
                    xlab(paste("MSE: ",format(mse_plain,digits = 3)))+
                    ylab("Daily New Cases")+
                    ggtitle("Linear Model for Daily New Case")+
                    theme_bw()+
                    scale_colour_manual(values = c("red", "darkgoldenrod2"),
                                        guide = guide_legend(override.aes = list(
                                            linetype = c("blank", "solid"),
                                            shape = c(16, NA))))
    })

            #ploting linear models for new cases
            output$tab2_plot1<-renderPlot({
                tab2_glm1()
            })

            #tab2_glm1 make line plot using linear models for seven-day roll mean of new cases
            tab2_glm2<-reactive({
                dat<-tab2_glmdat()
                set.seed(520)
                train_index <- createDataPartition(dat$new_cases, p = 0.7, list = F)
                training <- dat[train_index,]
                testing <- dat[-train_index,]
                #linear model for seven-day roll mean of new cases
                mean_lm <- lm(sdm~ driving + sdm_1 + sdm_2 + sdm_3 + sdm_4 + sdm_5, data = training)
                pred_mean <- predict(mean_lm, testing)

                #linear model for seven-day roll mean of new cases
                mean_lm <- lm(sdm~ driving + sdm_1 + sdm_2 + sdm_3 + sdm_4 + sdm_5, data = training)
                #predictions using the two models
                pred_mean <- predict(mean_lm, testing)
                mse_mean<-sqrt(mean(((testing$sdm-(pred_mean)))^2, na.rm = T))

                    ggplot(testing)+geom_point(aes(x=Date,y=sdm),color = "red")+
                        scale_x_date(date_breaks = "1 month",labels = scales::date_format("%Y-%m"))+
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
                        geom_line(aes(x=Date, y=pred_mean),color="darkgoldenrod2")+
                        xlab(paste("MSE: ",format(mse_mean,digits = 3)))+
                        ylab("Seven Day Roll Mean of Daily New Cases")+
                        ggtitle("Linear Model for Seven Day Roll Mean of Daily New Cases")+
                        theme_bw()+
                        scale_colour_manual(values = c("red", "darkgoldenrod2"),
                                            guide = guide_legend(override.aes = list(
                                                linetype = c("blank", "solid"),
                                                shape = c(16, NA))))
            })

            #ploting linear models for seven-day roll mean of new cases
            output$tab2_plot2<-renderPlot({
                tab2_glm2()
            })

#######################
#    END of TAB 2     #
#######################

    

########################################
# Function for TAB 3: Instantaneous Rt #
########################################

    subdata<-reactive({
        counties%>%
            filter(key==input$tab3_county)%>%
            drop_na(average)%>%
            filter(cumsum(new_case)>1)%>%
            filter(new_case>=0)
    })

    sub_plot<-reactive({
        subdata=subdata()
        Civ.R=EstimateR(subdata$average, T.Start=15:(nrow(subdata())-15), T.End=30:nrow(subdata()), method="NonParametricSI",SI.Distr=si, plot=TRUE,leg.pos=xy.coords(1,3))
        as.data.frame(cbind(subdata$date,Civ.R[["R"]][["Mean(R)"]],Civ.R[["R"]][["Quantile.0.975(R)"]],Civ.R[["R"]][["Quantile.0.025(R)"]]))%>%setNames(.,c("date","mean","high","low"))%>%mutate(date=as.Date(date, origin = "1970-01-01"))
    })

    i<-reactive({
        sub_pred<-subdata()%>%
            filter(new_case>0)%>%
            expandRows("new_case",drop = FALSE)

        incidence(sub_pred$date)
    })

    pred_1<-reactive({
        i()[1:(nrow(subdata()-14-7))]%>%
            project(R=list(sub_plot()$mean[(nrow(subdata())-14-7):(nrow(subdata())-14)],sub_plot$mean[(nrow(subdata())-7):nrow(subdata())]), si = si, n_days = 30, n_sim = 50,time_change = 21)
    })

    output$tab3_plot1=renderPlot({

        ggplot(sub_plot())+
            geom_line(aes(x=sub_plot()$date,y=sub_plot()$mean),col="red")+
            geom_ribbon(aes(x=sub_plot()$date,ymin=sub_plot()$low,ymax=sub_plot()$high),alpha=0.3,fill="coral1")+
            scale_x_date(date_breaks = "2 month")+
            labs(y="Rt",x="",title = paste0("Estimate Rt of ",input$tab3_county))+
            geom_hline(yintercept = 1,col="black", linetype="dashed")+
            scale_x_date(limits=c(),date_breaks = "1 month", labels = scales::date_format("%Y-%m"))+
            theme_bw()
    })
    
    

#########################
#Function for TBA4: rnn #
#########################
    output$tab4_plot1=renderPlot({
        rnn %>%
            #filter the selected county
            filter(placeforselect==input$tab4_county) %>%
            ggplot()+
            geom_line(aes(x=date,y=newcases7dayavg),color="black")+
            geom_line(aes(x=date,y=pred_newcases),color="red")+
            scale_x_date(limits=c(),date_breaks = "1 month", labels = scales::date_format("%Y-%m"))+
            theme(axis.text.x=element_text(angle=45, hjust=1))+
            theme_bw()+
            labs(x="",
                y="7-day-average New Cases",
                title="7-day-average New Cases vs RNN model prediction about last 14 days")
    })
    
    output$tab4_plot2=renderPlot({
        rnn %>%
            #filter the selected county
            filter(placeforselect==input$tab4_county,
                   date>"2020-11-01") %>%
            ggplot()+
            geom_line(aes(x=date,y=newcases7dayavg),color="black")+
            geom_line(aes(x=date,y=pred_newcases),color="red")+
            scale_x_date(limits=c(),date_breaks = "1 day", labels = scales::date_format("%m-%d"))+
            theme(axis.text.x=element_text(angle=45, hjust=1))+
            theme_bw()+expand_limits(y=0)+
            labs(x="",y="7-day-average New Cases",title="Zoomed in")
    })

} #End of server

# Run the application 
shinyApp(ui = ui, server = server)