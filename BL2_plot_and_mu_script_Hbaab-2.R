#Script for transforming BioLector output data


## Initial setup and loading of data for analysis

rm(list = ls(all.names = TRUE))

graphics.off()
#installing/loading the latest installr package:
install.packages("installr")
library(installr)

#Required packages are loaded into the current session. 
library(tidyverse)
library(reshape)
library(doBy)
library(htmf)
library(MASS)
library(gridExtra)
library(Biobase)
library(BiocGenerics)
library(parallel)

#If one or more pacage is not installed on your computer, uncomment and run appropriate line
#install.packages("tidyverse")
#install.packages("reshape")
#install.packages("doBy")
#install.packages("MASS")
#install.packages("gridExtra")
#install.packages("Biobase")
#install.packages("BiocGenerics")
#install.packages("parallel")

#Prior to first time use you will need to install Bioconductor using:
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install()

#And htmf custom package made by Chris Workman
install.packages("htmf_0.5.1.tar.gz", repos = NULL, type="source")

#After installation, load "htmf" in the library section

#Show files in working directory
setwd("")
dir()

#csvFiles is the output file from the BioLector in csv format
csvFiles = "filename.csv"

#desfiles is a table with information on each well.
#Desfile should be .txt format and tab separated. Column 1 should be well numbers A01, A02 etc
desFiles = "design_file.txt"


#Outfiles are the new files created by the script, that contains the processed data. 
#One file is created for each filter in the BioLector; e.g. Biomass and RFP
outfile1 = "Plate.txt"
#outfile2 = "Outfile_RFP.txt"

#Defining values for baseline and outlier correction 
baselineThr <- 5
outliersThr <- 7

read.table(file=desFiles)

pr2 = read_biolector2(file=csvFiles, designfile=desFiles, measureType="AmpRef_1",  blVersion=3)
cData(pr2) #write pr1 or pr2 from here,depending on the one used above

pInt <- plate_intervals(pr2)
print(pInt)

#The following determines the time window for the output data
#Defining the end time-point (in hours)
tMax = 100
#Trim data to this time window (in hours)
prs <- sliceByTimeIntervals(pr2, intervals = list(c(0, tMax)))
plate_intervals(prs)

plot_biolector(x=prs, spar=0.4) # plot each (raw) well
plot_biolector(x=prs, spar=0.4, log="y")

#Correct baseline and remove outliers 
prsn = plate_normalize(x=prs, filtersets="Biomass", spar=0.5,
                       valign=TRUE, targetTime=2, targetAmp=0.5,
                       baseline=FALSE, baselineThr=baselineThr,
                       outliers=TRUE, outliersThr=outliersThr)
range(measure(prsn$Biomass))

plot_biolector(x=prsn, spar=0.4, log="y", ylim=c(0.5, 150))
names(pr2)
#Construct output file

odat = plateRun2df(prsn, filtersets = "Biomass")
dim(odat)
biomass = odat[odat$type=="Biomass",]
dim(biomass) 

head(biomass)
write.table(file=outfile1, biomass, quote=FALSE, sep="\t", row.names=FALSE, dec=",")
#write.table(file=outfile2, biomass[,c(2:5,7:12)], quote=FALSE, sep="\t", row.names=FALSE, dec=",")
print(outfile1)


#Plotting and µ determination
growth_data<-biomass[,c(8,2,9,4,5,10)] %>% mutate(plate=1)

# head(growth_data)

growth_data$Condition<-as.factor(growth_data$Condition)

## Data pruning: Removal of read data from empty wells and removal of data
## points from plate that has been removed. If any additional points needs to be
## removed, it can be done by adding logical statements to this function.

growth_data <- growth_data %>% filter(Strain != 'None')

log_dat <- growth_data$measure

log_dat[log_dat>0] <- log(log_dat[log_dat>0])
log_dat[log_dat<=0] <- 0

growth_data <- mutate(growth_data, ln_Read = log_dat)

growth_avg_data <- growth_data %>% group_by(Strain, Condition, time) %>%
  summarise(Avg_measure = mean(measure), sd_measure = sd(measure), CI95_read = 1.96*sd(measure))

## Fitting of smoothing spline for detection of growth phases and calculation 
## of growth rates

growth_data_noblanks <- growth_data %>% filter(Strain!= 'Blank')

strain_vec <- unique(growth_data_noblanks$Strain_name)
interval_vec <- unique(growth_data_noblanks$interval)
condition_vec <- unique(growth_data_noblanks$Condition)

#Explanation: high: higher confidence, but it might hit a non-growth related phase 
#look at data in BioLection and decide on a suitable timeframe. E.g. 3.5 h
#in example here that is 73 cycles
window_size <- 73 ## Must be an odd number

#if time has been limited to include 1 phase only, the below will have to be set at 1
max_num_phases <- 1 ## Number of growth phases that will be detected for each strain (1 or higher)

growth_rate_tib <- tibble()

growth_phase_data_tib <- tibble()


n <- 0

for (i in strain_vec) {
  
  for (j in interval_vec) {
    
    for (k in condition_vec) {
      
      n <- n+1
      
      tmp_data <- filter(growth_data_noblanks, Strain == i, interval == j, Condition == k)
      
      if (dim(tmp_data)[1]>0) {
        
        # Fitting of the spline to the data - Growth data is log transformed to 
        # get the correct estimation of the highest growth rate.
        
        fit_data <- smooth.spline(x = tmp_data$time, y = log(10+tmp_data$measure),df = 10)
        
        predict_data <- predict(fit_data, deriv = 0)
        predict_deriv_data <- predict(fit_data, deriv = 1)
        
        peak_ind <- which(splus2R::peaks(predict_deriv_data$y, span = window_size))
        
        peak_mat <- matrix(data = c(predict_deriv_data$x[peak_ind],predict_deriv_data$y[peak_ind]), ncol = 2)
        
        peak_top <- peak_mat[order(peak_mat[,2], decreasing = TRUE)[1:max_num_phases],1]
        
        for (l in 1:length(peak_top)) {
          
          time_vec <- unique(tmp_data$time) %>% sort()
          
          growth_ind <- which(time_vec == peak_top[l])
          
          growth_time_low <- time_vec[growth_ind-floor(window_size/2)]
          
          growth_time_hi <- time_vec[growth_ind+floor(window_size/2)]
          
          gr_phase_data <- filter(tmp_data, time >= growth_time_low, time <= growth_time_hi) %>%
            mutate(Phase = l)
          
          growth_phase_data_tib <- rbind(growth_phase_data_tib,gr_phase_data)
          
          summary_growth_phase <- summary(lm(ln_Read~time,data = gr_phase_data))
          
          tmp_tib <- tibble(Strain = i, 
                            interval = j, 
                            pH = k,
                            Phase_size_order = l,
                            Phase_mid_point = peak_top[l],
                            growth_rate = summary_growth_phase$coefficients[2,1], 
                            p_val = summary_growth_phase$coefficients[2,4], 
                            Rsqr = summary_growth_phase$adj.r.squared)
          
          growth_rate_tib <- rbind(growth_rate_tib,tmp_tib)
          
        }
        
      }
      
    }
    
  }
  
}

## Data visualisation

# To save plots as files, change print_bin to 1

print_bin <- 1
print_h <- 2000
print_w <- 4000

p1<- ggplot(data = growth_data, mapping = aes(x = time, y = measure)) + 
  geom_point(mapping = aes(color = Condition)) + 
  facet_wrap(~Strain)

#if (print_bin==1) {

p2 <- p1 + scale_color_manual(values=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933")) +
  ggtitle("") + xlab("Time (h)") + ylab("Biomass") +
  theme_bw()

ggsave(p2, filename = "plot.pdf",height = 5, width = 20)

#}
