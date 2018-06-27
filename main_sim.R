#main file
# variable names use a suffix of "_mult" for the multiplicative model and
# "_add" for the additive model

library("SimDesign")
library("pbapply") # required by SimDesign
library("plyr")
library("ggplot2")
library("psych")
library("boot")

setwd("/home/mark/sims")

# Define levels of each factor
sample_sizes <- c(10, 30, 50, 100, 1000)
num_preds <- c(2, 3, 4, 5, 6)
mult_colls <- c(3, 6, 9)
res_dists <- c('norm','cont_norm', 'high_kurt', 'log_norm')

# test factor levels
#sample_sizes <- c(1000)
#num_preds <- c(3)
#mult_colls <- c(3)
#res_dists <- c('norm')

# create Design matrix for all permutations of conditions
Design <- expand.grid(
  sample_size = sample_sizes,
  num_pred = num_preds,
  mult_coll = mult_colls,
  res_dist = res_dists)

Replications <- 1000

#Run the simulation

# define simulation functions: "Generate", "Analyse", "Summarise", plus helper functions
source('functions_v5.R')

# define population covariance matrices using num_pred and condition value
# defines a list called fixed_objects to pass onto runSimulation
# also define simulation constants: alpha, sigma_sq, beta, IV_sds
source('Pop_Matrices_Adkins.R')

results <- runSimulation(design=Design, replications=Replications,
                         parallel=TRUE, generate=Generate, analyse=Analyse, summarise=Summarise,
                         fixed_objects=fixed_objects, packages = c('psych','fungible','car','plyr','MASS','boot'),
                         edit='none', save=TRUE, save_results=TRUE, progress = TRUE,
                         filename='results_09_22',save_generate_data = FALSE,save_seeds =FALSE)

# create a copy of results and add a column to track row number for graphing later
temp <- results
temp <- cbind(temp,index=1:nrow(temp))

#################################

# forgot to calculate efficiency of both bootstrap methods
temp$bs.bca.efficiency1 <- temp$bs.bca.upper1 - temp$bs.bca.lower1
temp$bs.bca.efficiency2 <- temp$bs.bca.upper2 - temp$bs.bca.lower2
temp$bs.bca.efficiency3 <- temp$bs.bca.upper3 - temp$bs.bca.lower3
temp$bs.bca.efficiency4 <- temp$bs.bca.upper4 - temp$bs.bca.lower4
temp$bs.bca.efficiency5 <- temp$bs.bca.upper5 - temp$bs.bca.lower5
temp$bs.bca.efficiency6 <- temp$bs.bca.upper6 - temp$bs.bca.lower6
temp$bs.bca.efficiency7 <- temp$bs.bca.upper7 - temp$bs.bca.lower7

temp$bs.perc.efficiency1 <- temp$bs.perc.upper1 - temp$bs.perc.lower1
temp$bs.perc.efficiency2 <- temp$bs.perc.upper2 - temp$bs.perc.lower2
temp$bs.perc.efficiency3 <- temp$bs.perc.upper3 - temp$bs.perc.lower3
temp$bs.perc.efficiency4 <- temp$bs.perc.upper4 - temp$bs.perc.lower4
temp$bs.perc.efficiency5 <- temp$bs.perc.upper5 - temp$bs.perc.lower5
temp$bs.perc.efficiency6 <- temp$bs.perc.upper6 - temp$bs.perc.lower6
temp$bs.perc.efficiency7 <- temp$bs.perc.upper7 - temp$bs.perc.lower7

#################################

# turn results into long format by CI_type,
# the first estimate is done outside the for loop to initialize the df "final_results"

# new column to create 
keycol <- "CI_type"

# master list of every column containing lower proportions, upper proportions, coverage, and efficiency
gathercols_lower <- c("wald.Lower Prob","bs.perc.Lower Prob", "bs.bca.Lower Prob")
gathercols_upper <- c("wald.Upper Prob","bs.perc.Upper Prob", "bs.bca.Upper Prob")
gathercols_coverage <- c("wald.Coverage","bs.perc.Coverage", "bs.bca.Coverage")
gathercols_efficiency <- c("wald.efficiency","bs.perc.efficiency", "bs.bca.efficiency")

# create final_results with lowerINT as the first new column
temp_lower <- paste(gathercols_lower,as.character(1),sep="")
final_results <- tidyr::gather_(data = temp, key=keycol, value = "lowerINT",temp_lower)

# create and extract upperINT
temp_upper <- paste(gathercols_upper,as.character(1),sep="")
upperINT <- tidyr::gather_(data = temp, key=keycol, value = "upperINT",temp_upper)$upperINT

# create and extract coverageINT
temp_coverage <- paste(gathercols_coverage,as.character(1),sep="")
coverageINT <- tidyr::gather_(data = temp, key=keycol, value = "coverageINT",temp_coverage)$coverageINT

# create and extract efficiencyINT
temp_efficiency <- paste(gathercols_efficiency,as.character(1),sep="")
efficiencyINT <- tidyr::gather_(data = temp, key=keycol, value = "efficiencyINT",temp_efficiency)$efficiencyINT

# combine master results(which has lower 1) and upper 1 and coverage
final_results <- cbind(final_results,upperINT,coverageINT,efficiencyINT)

# remove the old columns that were gathered together to remove redundant info and slim down the df
final_results <- final_results[, ! names(final_results) %in% c(temp_upper,temp_coverage,temp_efficiency), drop = F]

# do this same process for the the other remaining estimates
for(num in 2:I(max(as.numeric(as.character(temp$num_pred)))+1)){
  
  temp_lower <- paste(gathercols_lower,as.character(num),sep="")
  lower <- tidyr::gather_(data = temp, key=keycol, value = "lower",temp_lower)$lower
  
  temp_upper <- paste(gathercols_upper,as.character(num),sep="")
  upper <- tidyr::gather_(data = temp, key=keycol, value = "upper",temp_upper)$upper
  
  temp_coverage <- paste(gathercols_coverage,as.character(num),sep="")
  coverage <- tidyr::gather_(data = temp, key=keycol, value = "coverage",temp_coverage)$coverage
  
  (temp_efficiency <- paste(gathercols_efficiency,as.character(num),sep=""))
  efficiency <- tidyr::gather_(data = temp, key=keycol, value = "efficiency",temp_efficiency)$efficiency
  
  final_results <- cbind(final_results,lower,upper,coverage,efficiency)
  
  # rename columns using the appropriate estimate number after making estimate 1 into intercept
  names(final_results)[ncol(final_results)-3] <- paste("lowerX",num-1,sep="")
  names(final_results)[ncol(final_results)-2] <- paste("upperX",num-1,sep="")
  names(final_results)[ncol(final_results)-1] <- paste("coverageX",num-1,sep="")
  names(final_results)[ncol(final_results)] <- paste("efficiencyX",num-1,sep="")
  final_results <- final_results[, ! names(final_results) %in% c(temp_lower,temp_upper,temp_coverage,temp_efficiency), drop = F]
  names(final_results)
}

# convert vars into factors for plotting 
final_results$CI_type <- as.factor(final_results$CI_type)
final_results$sample_size <- as.factor(final_results$sample_size)
final_results$num_pred <- as.factor(final_results$num_pred)
final_results$mult_coll <- as.factor(final_results$mult_coll)

# rename levels for better labels in plotting
levels(final_results$CI_type) <- c("bs.bca.Lower Prob1"="BCa","bs.perc.Lower Prob1"="Perc","wald.Lower Prob1"="Wald")
levels(final_results$res_dist) <- c("norm"="Normal Distribution","comt_norm"="Contaminated-Normal Distribution","high_kurt"="Highly Kurtotic","log_norm"="Log-normal")

#################################

# get average coverage of all CIs within each model
coverage_vars <- c("coverageX1","coverageX2","coverageX3","coverageX4","coverageX5","coverageX6")
average_coverage <- adply(.data = final_results[,c("num_pred",coverage_vars)],
                          .margins = 1, 
                          .fun = function(x){"average_coverage"=sum(x[coverage_vars],na.rm=TRUE)/as.numeric(as.character(x$num_pred))})$V1
final_results <- cbind(final_results,average_coverage)

# get average efficiency of all CIs within each model
efficiency_vars <- c("efficiencyX1","efficiencyX2","efficiencyX3","efficiencyX4","efficiencyX5","efficiencyX6")
average_efficiency <- adply(.data = final_results[,c("num_pred",efficiency_vars)],
                          .margins = 1, 
                          .fun = function(x){"average_efficiency"=sum(x[efficiency_vars],na.rm=TRUE)/as.numeric(as.character(x$num_pred))})$V1
final_results <- cbind(final_results,average_efficiency)

# get average lower prob of all CIs within each model
lower_vars <- c("lowerX1","lowerX2","lowerX3","lowerX4","lowerX5","lowerX6")
average_lower <- adply(.data = final_results[,c("num_pred",lower_vars)],
                            .margins = 1, 
                            .fun = function(x){"average_lower"=sum(x[lower_vars],na.rm=TRUE)/as.numeric(as.character(x$num_pred))})$V1
final_results <- cbind(final_results,average_lower)

# get average upper prob of all CIs within each model
upper_vars <- c("upperX1","upperX2","upperX3","upperX4","upperX5","upperX6")
average_upper <- adply(.data = final_results[,c("num_pred",upper_vars)],
                       .margins = 1, 
                       .fun = function(x){"average_upper"=sum(x[upper_vars],na.rm=TRUE)/as.numeric(as.character(x$num_pred))})$V1
final_results <- cbind(final_results,average_upper)

#################################

# the first for columns are the design matrix
# last chunk grabs three columns (lower, upper, and coverage) for each estimate and "CI_type"
(num_vars_per_estimate <- (max(num_preds)+1)*4) # include (6 Xs + 1 intercept) * 4 variables each = 28
(num_vars_averages <- 4) # for CI_type
(total_vars <- length(colnames(final_results)))

var_nums_to_retain <- c(1:4,(total_vars - num_vars_averages - num_vars_per_estimate):total_vars)
trimmed_results <- final_results[,var_nums_to_retain]
names(trimmed_results)

# create subset which only includes the following levels
dists <- c("Normal Distribution", "Contaminated-Normal Distribution",
            "Highly Kurtotic","Log-normal")
preds <- c(2, 3, 4, 5, 6)
sams <- c(10,30, 50, 100,1000)
mult <- c(3, 6, 9)

graph_data <- subset(trimmed_results, res_dist %in% dists & 
                       num_pred %in% preds & 
                       sample_size %in% sams &
                       mult_coll %in% mult)

# add a row number variable
graph_data <- cbind(graph_data,1:nrow(graph_data))
colnames(graph_data)[length(colnames(graph_data))] <- "index_full"

# confidence band
CI_level <- 1 - fixed_objects[["alpha"]]
size <- 2 * sqrt(CI_level * fixed_objects[["alpha"]] / Replications)
num_conditions <- length(row.names(graph_data))

# dichtomize data using the empirical rule as a guide
final_results$within_CI_band <- as.numeric(final_results$average_coverage >= (CI_level - size) &
                               final_results$average_coverage <= (CI_level + size))
sum(final_results$within_CI_band)

# define range of indexs to plot the CI band
CI_band_xmin <- 0
CI_band_xmax <- length(row.names(graph_data))
wald_xmin <- min(which(graph_data$CI_type == "Wald"))
wald_xmax <- max(which(graph_data$CI_type == "Wald"))
perc_xmin <- min(which(graph_data$CI_type == "Perc"))
perc_xmax <- max(which(graph_data$CI_type == "Perc"))
bca_xmin <- min(which(graph_data$CI_type == "BCa"))
bca_xmax <- max(which(graph_data$CI_type == "BCa"))

ggplot2::ggplot() +
  geom_line(data=graph_data,aes(y=average_coverage,x=index_full,col=res_dist,group=1),size = 1) +
  geom_point(data=graph_data,aes(y=average_coverage,x=index_full,shape=CI_type),size = 1.5) + 
  scale_x_continuous(limits = c(CI_band_xmin,CI_band_xmax),breaks = c(wald_xmax,perc_xmax,bca_xmax)) +
  geom_hline(yintercept = 0.95, color = "red") +
  geom_rect(aes(xmin=CI_band_xmin,xmax=CI_band_xmax,ymin=CI_level - size, ymax = CI_level + size),colour = "blue",alpha = .2,fill="blue") +
  geom_rect(aes(xmin=wald_xmin,xmax=wald_xmax,ymin=0.80,ymax=1.00),colour = "grey",alpha = 0.1) +
  geom_rect(aes(xmin=perc_xmin,xmax=perc_xmax,ymin=0.80,ymax=1.00),colour = "grey",alpha = 0.2) +
  geom_rect(aes(xmin=bca_xmin,xmax=bca_xmax,ymin=0.80,ymax=1.00),colour = "grey",alpha = 0.3) +
  labs(x="Design Condition",y="Average Coverage") +
  ggtitle("Average Confidence Interval Coverage by Design Condition") + 
  scale_color_discrete(name="Residual Distribution Type") +
  theme(plot.background = element_rect(size = 1, color = "white", fill = "white"),
        legend.position = "bottom",
        legend.margin = margin(10, 10, 10, 10),
        legend.background = element_rect(color = 'black', fill = 'white', size = 1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=13),
        plot.title = element_text(size = 15,face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size=10,hjust=0.5),
        axis.text.x = element_text(size=10,hjust=0.5),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.grid.minor.x = element_line(),
        panel.grid.major.x = element_line(colour="gray50", size=0.5),
        strip.background = element_rect(fill = "white",colour="black",size=1.5))


final_results_sub <- final_results

# the subset function is used to remove cases such as ss == 10
final_results_sub <- subset(final_results, sample_size != 10)

# This is the variable name to summarize data for
var_interest <- "average_coverage"

marginals_sample_size <- get_summary(var_interest,c("sample_size","CI_type"))
marginals_sample_size <- cbind(marginals_sample_size,get_summary("average_efficiency",c("sample_size","CI_type"))[,4:5])
marginals_sample_size <- cbind(marginals_sample_size,get_summary("average_lower",c("sample_size","CI_type"))[,4:5])
marginals_sample_size <- cbind(marginals_sample_size,get_summary("average_upper",c("sample_size","CI_type"))[,4:5])

summary_sample_size <- get_summary(var_interest,c("sample_size","CI_type"))
colnames(marginals_sample_size) <- c("Sample Size", "Confidence Interval type","N",rep(c("Mean","SD"),3))

summary_sample_size_eff <- get_summary("average_efficiency",c("sample_size","CI_type"))
colnames(summary_sample_size_eff) <- c("Sample Size", "Confidence Interval type","N","Mean","SD")

summary_res_dist <- get_summary(var_interest,c("res_dist","CI_type"))
colnames(summary_res_dist) <- c("Residual Distribution Type", "Confidence Interval type","N","Mean","SD")

summary_num_pred <- get_summary(var_interest,c("num_pred","CI_type"))
colnames(summary_num_pred) <- c("Number of Indepenent Variables", "Confidence Interval type","N","Mean","SD")

summary_mult_coll <- get_summary(var_interest,c("mult_coll","CI_type"))
colnames(summary_mult_coll) <- c("Multi-collinearity", "Confidence Interval type","N","Mean","SD")

summary_all <- get_summary(var_interest,c("res_dist","sample_size","mult_coll","num_pred","CI_type"))
colnames(summary_all) <- c("Residual Distribution Type",
                           "Sample Size",
                           "Multicollinearity",
                           "Number of Indepenent Variables",
                           "Confidence Interval type","N","Mean","SD")


# breakdown by design condition and find best coverage and which approach was used
new_list <- ddply(final_results_sub,.variables = c("res_dist","sample_size","mult_coll","num_pred"),
      summarise,
      best_coverage = max(average_coverage),
      # find best method, but in case of ties concatenate the methods together
      best_method = paste(CI_type[which(best_coverage == average_coverage)],collapse ="/"),
      worst_coverage = min(average_coverage),
      worst_method = paste(CI_type[which(worst_coverage == average_coverage)],collapse ="/")
      )


freqs <- count(new_list,vars=c("best_method","worst_method"))
freqs$perc <- paste(c(perc = round(100 * freqs$freq / sum(freqs$freq),2)),"%",sep="")
names(freqs) <- c("Best Method", "Worst Method","Freq","Percent")

(freqs <- freqs[order(-freqs$Freq),])
# examine patterns of average coverage below a certain threshold
vnames <- list(set_varnames = c(sample_size="Sample Size", 
                                num_pred="Number of IVs",
                                mult_coll="Multicollinearity",
                                res_dist= "Distribution",
                                within_CI_band = "Within Confidence Band"))
lnames <- list(sample_size = sample_sizes,
               num_pred = num_preds,
               mult_coll = mult_colls,
               res_dist = res_dists,
               within_CI_band = c("Out", "In"))

#vcd::mosaic(within_CI_band~ sample_size + num_pred + mult_coll + res_dist, 
#            data = final_results[,c("sample_size","num_pred","mult_coll","res_dist","average_coverage","within_CI_band")],
#            shade=TRUE, labeling_args=vnames, set_labels=lnames)

count(final_results,c("within_CI_band","CI_type","res_dist"))

mytable <- xtabs(~within_CI_band+CI_type+res_dist,data=final_results_sub)
table_df <- ftable((mytable))
str(table_df)
table_df[["ftable"]]

# create a subset of final_results with average coverage within the CI band
final_results_good_coverage <- subset(final_results_sub, within_CI_band == 1)
# breakdown by design condition and find best efficiency and which approach was used
efficiency_summary <- ddply(final_results_good_coverage,.variables = c("res_dist","sample_size","mult_coll","num_pred"),
                            summarise,
                            lowest_efficiency = min(average_efficiency),
                            # find best method, but in case of ties concatenate the methods together
                            best_method = paste(CI_type[which(lowest_efficiency == average_efficiency)],collapse ="/"),
                            highest_efficiency = max(average_efficiency),
                            worst_method = paste(CI_type[which(highest_efficiency == average_efficiency)],collapse ="/")
)


freqs_eff <- count(efficiency_summary,vars=c("best_method","worst_method"))
freqs_eff$perc <- paste(c(perc = round(100 * freqs_eff$freq / sum(freqs_eff$freq),2)),"%",sep="")
names(freqs_eff) <- c("Best Method", "Worst Method","Freq","Percent")

efficiency_summary$highest_efficiency
eff_perc_conv <- subset(efficiency_summary,efficiency_summary$best_method == "Perc" & efficiency_summary$worst_method == "Wald")
count(eff_perc_conv,vars="res_dist")
count(eff_perc_conv,vars="sample_size")
count(eff_perc_conv,vars="mult_coll")
count(eff_perc_conv,vars="num_pred")

eff_conv_conv <- subset(efficiency_summary,efficiency_summary$best_method == "Wald" & efficiency_summary$worst_method == "Wald")
count(eff_conv_conv,vars="res_dist")
count(eff_conv_conv,vars="sample_size")
count(eff_conv_conv,vars="mult_coll")
count(eff_conv_conv,vars="num_pred")


xtabs(~within_CI_band+CI_type+res_dist,data=final_results_sub)
eff_perc_conv <- subset(efficiency_summary,efficiency_summary$best_method == "Perc")
eff_conv_conv <- subset(efficiency_summary,efficiency_summary$best_method == "Wald")

(perc_conv_xtab <- ftable(xtabs(~res_dist+sample_size+mult_coll+num_pred,eff_perc_conv)))
sum(perc_conv_xtab[,])
(conv_conv_xtab <- ftable(xtabs(~res_dist+sample_size+mult_coll+num_pred,eff_conv_conv)))
sum(conv_conv_xtab[,])
conv_conv_xtab[conv_conv_xtab[,] == 1] <- 2
temp_cont_table <- perc_conv_xtab + conv_conv_xtab
temp_cont_table[temp_cont_table[,] == 0] <- "-"
temp_cont_table[temp_cont_table[,] == 1] <- "X"
temp_cont_table[temp_cont_table[,] == 2] <- "O"


vcd::mosaic(vcd::structable(conv_conv_xtab))
