#----------------------------------------------------------------------------
##############################################################################
#@@@@@@@@@@@@@@@ Sweetrain @@@@@@@@@@@@@@@#
# author: Evan Heintz

# description: reads in Sweetrain data, expands points, and fits logistic 
# regression models to the data.

# requirements: file path to 'sweetrain_processed.csv'.

# outputs: expanded Sweetrain dataframe, logistic regression models, and
# robust standard errors.
##############################################################################
#----------------------------------------------------------------------------

#------ load libraries ------#
packages <- c("stats","tidyverse","tiff","ggforce","gridExtra","viridis",
              "ggthemes","sandwich","clubSandwich","fdrtool")
lapply(packages, library, character.only = TRUE)
theme_set(theme_minimal())


if (file.exists('/Users/evanheintz')) {
  sweetrain <- read.csv("/Users/evanheintz/Dropbox/Evan_Thesis/math470-atoh7-modeling-thesis/Data/sweetrain_processed.csv")
} else {
  # sweetrain <- enter file path for sweetrain
}

sweetrain <- sweetrain %>% mutate("intensity" = 1)
# ggplot(sweetrain, aes(x=x,y=y,color=t))+geom_point()+scale_color_viridis_c()+labs(title="Sweetrain")

current <- list()
for(i in 1:max(sweetrain$t)){
  #convert to data frames
  current[[i]] <- as.data.frame(sweetrain %>% filter(t == i))
  #rename columns
  colnames(current[[i]]) <- c("t", "x", "y", "intensity")
  #convert to numeric values
  current[[i]]$x <- as.numeric(current[[i]]$x)
  current[[i]]$y <- as.numeric(current[[i]]$y)
}

current_new <- list()
for(i in 1:(length(current)-1)){
  moment <- current[[i]]
  for(j in (i+1):length(current)){
    moment <- rbind(moment, (current[[j]] %>% 
                               mutate(intensity = replace(intensity, 0:length(current[[j]]$intensity), 0)) %>%
                               mutate(t = replace(t, 0:length(current[[j]]$intensity), i))))
  }
  current_new[[i]] <- as.data.frame(moment)
}
current_new[[length(current)]] <- current[[length(current)]]


# e1 <- ggplot(sample_n(sweetrain, 100000), aes(x=x,y=y,color=t)) + geom_point() + scale_color_viridis_c() + 
#  labs(x="x", y="y", title="Sweetrain Sample")
# e1
# ggplot(current[[length(current)]], aes(x=y,y=x,color=t)) + geom_point() 


#create dataframe with results
#takes around 4-5 minutes
st <- Sys.time()
data_new <- data.frame()
for (i in 1:length(current_new)){
  data_new <- rbind(data_new,current_new[[i]])
}
print(Sys.time()-st)

#function for finding distance to nearest "on" point
euclidean_dist <- function(x,y,x0,y0){
  sqrt((x-x0)^2 + (y-y0)^2)
}

timesteps = length(current)
#for each off point, minimize the distance to all "on" points
dist <- function(data,timesteps){
  data_final <- data.frame()
  for (i in (min(data_new$t)+1):timesteps) {
    on_t <- data %>% filter(t == (i-1) & intensity == 1)
    data_t <- data %>% filter(t == i)
    
    distance <- numeric(length(data_t$t))
    for (j in 1:length(data_t$t)) {
      if(length(on_t$t) == 0) {distance[j] <- 0}
      else {distance[j] <- min(euclidean_dist(x = data_t$x[j], y = data_t$y[j], 
                                              x0 = on_t$x, y0 = on_t$y), na.rm=TRUE)} 
    }
    data_t <- data_t %>% 
      mutate(distance = distance)
    
    data_final <- rbind(data_final,data_t)
  }
  return(data_final)
}

#-----functions-----#
#skip to line 191 to compute hull_center
# radial
radial <- function(x,y){
  sqrt((x - hull_center[1])^2 + (y - hull_center[2])^2)
}
# rotational
angle <- function(x,y){
  atan2((y-hull_center[2]), -1*(x-hull_center[1])) - pi/4
}

#apply function to data!
#takes around 17 minutes to run :/ but yay!
st <- Sys.time()
data_final_new <- dist(data=data_new,timesteps=timesteps)
print(Sys.time()-st)


# visualize
ggplot(sample_n(data_final_new %>% filter(intensity==1),50000), aes(x=x,y=y,color=t)) + geom_point()

### if x and y values are swapped, then run below
# x <- data_final_new$x
# y <- data_final_new$y
# data_final_new$x <- y
# data_final_new$y <- x

set.seed(470)
split <- sample(c(TRUE, FALSE), nrow(data_final_new), replace=TRUE, prob=c(0.7,0.3))
data_train_new <- data_final_new[split,]
data_test_new <- data_final_new[!split,]


st <- Sys.time()
logistic_new <- glm(intensity ~ x + y + t + distance + angle(x=x,y=y) + radial(x=x,y=y), 
                    data=data_train_new, family="binomial")
print(Sys.time()-st)

st <- Sys.time()
logistic_interaction_new <- glm(intensity ~ x + y + t + distance + angle(x=x,y=y):t + radial(x=x,y=y):t, 
                                data=data_train_new, family="binomial")
print(Sys.time()-st)

# final logistic model with all components is 'logistic_all_new'
st <- Sys.time()
logistic_all_new <- glm(intensity ~ x + y + t + distance + angle(x=x,y=y) + 
                          radial(x=x,y=y) + angle(x=x,y=y):t + radial(x=x,y=y):t, 
                        data=data_train_new, family="binomial")
print(Sys.time()-st)


#summaries
summary(logistic_new)
summary(logistic_interaction_new)
summary(logistic_all_new)


#predictions
set.seed(470)

data_test_preds <- data_test_new


preds <- predict(logistic_new, data_test_preds, type="response")
data_preds <- data_test_preds %>% select(-intensity) %>% mutate(preds = as.numeric(ifelse(preds < mean(preds), 0, 1))) #quantile(preds,probs=0.25)

p1 <- ggplot(data=sample_n((data_preds %>%  filter(preds == 1)), 50000), aes(x = x, y = y, color = t)) + geom_point(alpha=0.1) + 
  scale_color_viridis_c() + labs(title="Regular")

p2 <- ggplot(sample_n(sweetrain, 50000), aes(x = x, y = y, color = t)) + geom_point(alpha=0.1) + 
  scale_color_viridis_c() + labs(title="Sample") 


preds <- predict(logistic_interaction_new, data_test_preds, type="response")
data_preds_interaction <- data_test_preds %>% select(-intensity) %>% mutate(preds = as.numeric(ifelse(preds < mean(preds), 0, 1)))

p3 <- ggplot(data=sample_n((data_preds_interaction %>%  filter(preds == 1)), 50000), aes(x = x, y = y, color = t)) + geom_point(alpha=0.1) + 
  scale_color_viridis_c() + labs(title="Interaction")


# select which quantile (or mean) to use for predictions based on error:
preds <- predict(logistic_all_new, data_test_preds, type="response")

error <- c()
for(i in 1:20){
  if(i %in% 1:19){
    data_preds_all <- data_test_preds %>% select(-intensity) %>% mutate(preds = as.numeric(ifelse(preds < quantile(preds, i/20), 0, 1)))
    error[i] <- mean(sum(data_preds_all$preds - data_test_new$intensity)^2)
  } else{ # where i=20 will use the mean as a cutoff
    data_preds_all <- data_test_preds %>% select(-intensity) %>% mutate(preds = as.numeric(ifelse(preds < mean(preds), 0, 1)))
    error[i] <- mean(sum(data_preds_all$preds - data_test_new$intensity)^2)
  }
}
plot(error)

quantile <- which.min(error)/20
data_preds_all <- data_test_preds %>% select(-intensity) %>% mutate(preds = as.numeric(ifelse(preds < quantile(preds, 0.91), 0, 1))) # quantile(preds,quantile)
min <- min((data_preds_all %>%  filter(preds == 1))$t)
# find smallest error that includes all time points
# 0.91 in the quantile() function gives the smallest error while still including all timepoints
  
p4 <- ggplot(data=sample_n(data_preds_all %>% filter(preds==1), 10000), 
             aes(x = x, y = y, color = t)) + geom_point(alpha=0.8) + 
  scale_color_viridis_c() # + labs(title="All") 
p4

og <- ggplot(sample_n(sweetrain %>% filter(t>min & t<15), 500), aes(x = x, y = y, color = t)) + geom_point(alpha=0.1) + 
  scale_color_viridis_c() + labs(title="Sample") 
og

grid.arrange(p4, og)

grid.arrange(p1,p3,p4,p2)



#cluster robust SE
crse <- sqrt(diag(vcovCL(logistic_all_new, cluster=data_train_new$t, type="HC0")))

coefs <- numeric(9)
for(i in 1:9){
  coefs[i] <- as.numeric(logistic_all_new$coefficients[i])
}

z <- coefs/crse
pvalues <- 2*pnorm(-1*abs(z))

#new p values


#%%%%%%%%%%% Expansion %%%%%%%%%%%#
#----------------------------------

#convex hull

filtered<-current[[length(current)]]

#calculate the minimal convex polygon that contains them
hull<-chull(filtered$x,filtered$y)
hull_df<-data.frame(x=filtered$x[hull],y=filtered$y[hull])

#centers and radii
hull_center <- c((min(filtered$x[hull])+max(filtered$x[hull]))/2, 
                 (min(filtered$y[hull])+max(filtered$y[hull]))/2)
xradius <- (max(hull_df$x) - min(hull_df$x))/2 
yradius <- (max(hull_df$y) - min(hull_df$y))/2 

plot(hull_df, type='l')
points(hull_df)
points(hull_center[1], hull_center[2], col='red')


current_expanded_new <- list()
for(i in 1:max(sweetrain$t)){
  #convert to data frames
  current_expanded_new[[i]] <- as.data.frame(sweetrain %>% filter(t == i))
  #rename columns
  colnames(current_expanded_new[[i]]) <- c("t", "x", "y", "intensity")
  #convert to numeric values
  x <- as.numeric(current_expanded_new[[i]]$x)
  y <- as.numeric(current_expanded_new[[i]]$y)
  current_expanded_new[[i]]$x <- y
  current_expanded_new[[i]]$y <- x
}

ggplot(current_expanded_new[[timesteps]],aes(x=x,y=y,color=t))+geom_point()





#%%%%%%%%%% approach 2, equal ratios
#%%---------------------------------------

st<- Sys.time()
major_radius <- c()
minor_radius <- c()
ratio <- xradius / yradius
for(i in 1:timesteps){
  #minor radius, since minimum y always seems to be right value to use in videos
  minor_radius[i] <- hull_center[2] - min(current_expanded_new[[i]]$y)
    #if(length(current_expanded_new[[i]]$y) == 0) 0 else min(current_expanded_new[[i]]$y)
  
  #assuming radius ratios are the same throughout
  major_radius[i] <- minor_radius[i] * ratio
}
print(Sys.time()-st)


#plots
plot(minor_radius)
plot(major_radius)


# expansion for approach 2
st<- Sys.time()
expand_x <- list()
expand_y <- list()
for(i in 1:timesteps){  
  #get centered at origin
  current_expanded_new[[i]]$x <- current_expanded_new[[i]]$x - hull_center[1]
  current_expanded_new[[i]]$y <- current_expanded_new[[i]]$y - hull_center[2]
  
  expand_x[[i]] <- current_expanded_new[[i]]$x * xradius / major_radius[i]
  expand_y[[i]] <- current_expanded_new[[i]]$y * yradius / minor_radius[i] 
  
  #recenter
  current_expanded_new[[i]]$x <- expand_x[[i]] + hull_center[1]
  current_expanded_new[[i]]$y <- expand_y[[i]] + hull_center[2] 
  
}
print(Sys.time()-st)

#%%%%% create data

current_new <- list()
for(i in 1:(length(current_expanded_new)-1)){
  moment <- current_expanded_new[[i]]
  for(j in (i+1):length(current)){
    moment <- rbind(moment, (current_expanded_new[[j]] %>% 
                               mutate(intensity = replace(intensity, 0:length(current_expanded_new[[j]]$intensity), 0)) %>%
                               mutate(t = replace(t, 0:length(current_expanded_new[[j]]$intensity), i))))
  }
  current_new[[i]] <- as.data.frame(moment)
}
current_new[[length(current_expanded_new)]] <- current_expanded_new[[length(current_expanded_new)]]

#create dataframe with results
#takes around 4-5 minutes
st <- Sys.time()
data_new<-data.frame()
for (i in 1:length(current_new)){
  data_new<-rbind(data_new,current_new[[i]])
}
print(Sys.time()-st)


###### visualize
ggplot(sample_n(data_new %>% filter(intensity==1), 100000), aes(x=x,y=y,color=t)) + geom_point(alpha=0.5) + scale_color_viridis_c()






#%%%%%%%%%% approach 1, calculating max_index
#%%---------------------------------------

# re-run lines 256-267 for the correct 'current_expanded_new' list!!

max_index <- c() 
major_radius_1 <- c()
for(i in 4:timesteps){
  max_index <- which(sqrt((current_expanded_new[[i]]$x - hull_center[1])^2 + (current_expanded_new[[i]]$y - hull_center[2])^2) == 
                       max(sqrt((current_expanded_new[[i]]$x - hull_center[1])^2 + (current_expanded_new[[i]]$y - hull_center[2])^2)))
  
  x <- current_expanded_new[[i]]$x[max_index]
  y <- current_expanded_new[[i]]$y[max_index]
  
  major_radius_1[i] <- (x - hull_center[1]) / (sqrt(1 - ((y - hull_center[2]) / minor_radius[i])^2))
}
plot(major_radius_1)

# expansion using approach 1
st<- Sys.time()
expand_x <- list()
expand_y <- list()
for(i in 1:timesteps){  
  #get centered at origin
  current_expanded_new[[i]]$x <- current_expanded_new[[i]]$x - hull_center[1]
  current_expanded_new[[i]]$y <- current_expanded_new[[i]]$y - hull_center[2]
  
  expand_x[[i]] <- current_expanded_new[[i]]$x * xradius / major_radius_1[i]
  expand_y[[i]] <- current_expanded_new[[i]]$y * yradius / minor_radius[i] 
  
  #recenter
  current_expanded_new[[i]]$x <- expand_x[[i]] + hull_center[1]
  current_expanded_new[[i]]$y <- expand_y[[i]] + hull_center[2] 
  
}
print(Sys.time()-st)


#%%%%% create data

current_new <- list()
for(i in 1:(length(current_expanded_new)-1)){
  moment <- current_expanded_new[[i]]
  for(j in (i+1):length(current)){
    moment <- rbind(moment, (current_expanded_new[[j]] %>% 
                               mutate(intensity = replace(intensity, 0:length(current_expanded_new[[j]]$intensity), 0)) %>%
                               mutate(t = replace(t, 0:length(current_expanded_new[[j]]$intensity), i))))
  }
  current_new[[i]] <- as.data.frame(moment)
}
current_new[[length(current_expanded_new)]] <- current_expanded_new[[length(current_expanded_new)]]

#create dataframe with results
#takes around 4-5 minutes
st <- Sys.time()
data_new<-data.frame()
for (i in 1:length(current_new)){
  data_new<-rbind(data_new,current_new[[i]])
}
print(Sys.time()-st)


###### visualize
ggplot(sample_n(data_new %>% filter(intensity==1 & t > 36), 100000), aes(x=x,y=y,color=t)) + geom_point(alpha=0.5) + scale_color_viridis_c()
ggplot(sample_n(sweetrain, 100000), aes(x=x,y=y,color=t)) + geom_point(alpha=0.5) + scale_color_viridis_c()








# final approach, using equal ratios assumption for smoothness
#%%---------------------------------------

# look at graph, choose reasonable time point where fluctuation stops
# using major_radius from approach 2

# for minor radius at each timepoint
xmin <- c()
min <- c()
t <- c()
for(i in 1:timesteps){
  data <- sweetrain %>% filter(t==i)
  min[i] <- min(data$y)
  xmin[i] <- data$x[which(data$y == min(data$y))]
  t[i] <- i
}

ggplot() +
  geom_point(data=data_final_new %>% filter(t==timesteps | t == timesteps-1 | t == timesteps-2), aes(x=x, y=y), color="lightgrey") +
  geom_point(data=data.frame(x=xmin, y=min, t=t), aes(x=x, y=y, color=t),size=3,alpha=0.8) +
  scale_color_paletteer_c("viridis::cividis") +
  labs(color = "t")



diff <- c()
for(i in 4:(timesteps-1)){
  diff[i] <- major_radius[i+1]-major_radius[i]
}
plot(diff)
which(diff==max(diff, na.rm=TRUE))

which(abs(minor_radius[4:timesteps] - mean(minor_radius[4:timesteps])) == min(abs(minor_radius[4:timesteps] - mean(minor_radius[4:timesteps]))))
which(abs(major_radius[4:timesteps] - mean(major_radius[4:timesteps])) == min(abs(major_radius[4:timesteps] - mean(major_radius[4:timesteps]))))

#assign that radius to earlier timepoints
for(j in 1:34){
  major_radius[j] <- major_radius[34]
  minor_radius[j] <- minor_radius[34]
}
#plots
plot(minor_radius)
plot(major_radius)



#expansion using final approach
st<- Sys.time()
expand_x <- list()
expand_y <- list()
for(i in 1:timesteps){  
  #get centered at origin
  current_expanded_new[[i]]$x <- current_expanded_new[[i]]$x - hull_center[1]
  current_expanded_new[[i]]$y <- current_expanded_new[[i]]$y - hull_center[2]
  
  expand_x[[i]] <- current_expanded_new[[i]]$x * xradius / major_radius[i]
  expand_y[[i]] <- current_expanded_new[[i]]$y * yradius / minor_radius[i] 
  
  #recenter
  current_expanded_new[[i]]$x <- expand_x[[i]] + hull_center[1]
  current_expanded_new[[i]]$y <- expand_y[[i]] + hull_center[2] 
  
}
print(Sys.time()-st)



#logistic regression on expanded data

current_new <- list()
for(i in 1:(length(current_expanded_new)-1)){
  moment <- current_expanded_new[[i]]
  for(j in (i+1):length(current)){
    moment <- rbind(moment, (current_expanded_new[[j]] %>% 
                               mutate(intensity = replace(intensity, 0:length(current_expanded_new[[j]]$intensity), 0)) %>%
                               mutate(t = replace(t, 0:length(current_expanded_new[[j]]$intensity), i))))
  }
  current_new[[i]] <- as.data.frame(moment)
}
current_new[[length(current_expanded_new)]] <- current_expanded_new[[length(current_expanded_new)]]

#create dataframe with results
#takes around 4-5 minutes
st <- Sys.time()
data_new<-data.frame()
for (i in 1:length(current_new)){
  data_new<-rbind(data_new,current_new[[i]])
}
print(Sys.time()-st)


ggplot(sample_n(data_new %>% filter(intensity==1 & t < 35 | t == timesteps), 70000), aes(x=x,y=y,color=t)) + geom_point(alpha=0.5) + scale_color_viridis_c() 

#apply distance function to data!
#takes around 17 minutes to run :/ but yay!
st <- Sys.time()
data_expanded_final_new <- dist(data=data_new,timesteps=timesteps)
print(Sys.time()-st)

ggplot(sample_n(data_expanded_final_new %>% filter(intensity==1),100000), aes(x=x,y=y,color=t)) + geom_point()+
  scale_color_viridis_c()
data_expanded_final_new <- data_expanded_final_new %>% filter(t < 74)


x <- data_expanded_final_new$x
y <- data_expanded_final_new$y
data_expanded_final_new$x <- y
data_expanded_final_new$y <- x


set.seed(470)
split <- sample(c(TRUE, FALSE), nrow(data_expanded_final_new), replace=TRUE, prob=c(0.7,0.3))
data_expanded_train_new <- data_expanded_final_new[split,]
data_expanded_test_new <- data_expanded_final_new[!split,]

st <- Sys.time()
logistic_expanded_new <- glm(intensity ~ x + y + t + distance + angle(x=x,y=y) + 
                                   radial(x=x,y=y), 
                                 data=data_expanded_train_new, family="binomial")
print(Sys.time()-st)

st <- Sys.time()
logistic_expanded_interaction_new <- glm(intensity ~ x + y + t + distance + angle(x=x,y=y):t + radial(x=x,y=y):t, 
                                 data=data_expanded_train_new, family="binomial")
print(Sys.time()-st)

st <- Sys.time()
logistic_expanded_all_new <- glm(intensity ~ x + y + t + distance + angle(x=x,y=y) + 
                                   radial(x=x,y=y) + angle(x=x,y=y):t + radial(x=x,y=y):t, 
                                 data=data_expanded_train_new, family="binomial")
print(Sys.time()-st)




#summaries
summary(logistic_expanded_new)
summary(logistic_expanded_interaction_new)
summary(logistic_expanded_all_new)


#predictions
set.seed(470)

data_test_preds <- data_expanded_test_new
preds <- predict(logistic_expanded_new, data_test_preds, type="response")
data_preds <- data_test_preds %>% mutate(preds = as.numeric(ifelse(preds < mean(preds), 0, 1)))

p1 <- ggplot(data=sample_n((data_preds %>%  filter(preds == 1)), 50000), aes(x = x, y = y, color = t)) + geom_point(alpha=0.3) + 
  scale_color_viridis_c() + labs(title="Regular") + theme_minimal()

p2 <- ggplot(sample_n(data_expanded_final_new, 50000), aes(x = x, y = y, color = t)) + geom_point(alpha=0.3) + 
  scale_color_viridis_c() + labs(title="Sample") + theme_minimal()


preds <- predict(logistic_expanded_interaction_new, data_test_preds, type="response")
data_preds_interaction <- data_test_preds %>% mutate(preds = as.numeric(ifelse(preds < mean(preds), 0, 1)))

p3 <- ggplot(data=sample_n((data_preds_interaction %>%  filter(preds == 1)), 50000), aes(x = x, y = y, color = t)) + geom_point(alpha=0.3) + 
  scale_color_viridis_c() + labs(title="Interaction") + theme_minimal()


# final logistic logistic model used for expanded data predictions
# select which quantile (or mean) to use for predictions based on error:
preds <- predict(logistic_expanded_all_new, data_test_preds, type="response")

error <- c()
for(i in 1:20){
  if(i %in% 1:19){
    data_preds_all <- data_test_preds %>% select(-intensity) %>% mutate(preds = as.numeric(ifelse(preds < quantile(preds, i/20), 0, 1)))
    error[i] <- sqrt(sum(data_preds_all$preds - data_expanded_test_new$intensity)^2)
  } else{ # where i=20 will use the mean as a cutoff
    data_preds_all <- data_test_preds %>% select(-intensity) %>% mutate(preds = as.numeric(ifelse(preds < mean(preds), 0, 1)))
    error[i] <- sqrt(sum(data_preds_all$preds - data_expanded_test_new$intensity)^2)
  }
}
plot(error)

data_preds_all <- data_test_preds %>% select(-intensity) %>% mutate(preds = as.numeric(ifelse(preds < mean(preds), 0, 1))) # quantile(preds,which.min(error)/20)
min <- min((data_preds_all %>%  filter(preds == 1))$t)
# find smallest error that includes all time points

p4 <- ggplot(data=sample_n((data_preds_all %>%  filter(preds == 1)), 10000), aes(x = x, y = y, color = t)) + geom_point(alpha=0.8) + 
  scale_color_viridis_c() #+ labs(title="All") 
p4

grid.arrange(p1,p3,p4,p2)

#cluster robust SE
crse <- sqrt(diag(vcovCL(logistic_expanded_all_new, cluster=data_expanded_train_new$t, type="HC0")))

coefs <- numeric(9)
for(i in 1:9){
  coefs[i] <- as.numeric(logistic_expanded_all_new$coefficients[i])
}

z <- coefs/crse
pvalues <- 2*pnorm(-1*abs(z))

#new p values



###%%%%%%% logistic regression with ventral patch included
# include ventral patch indicator in data
ventral_indicator <- ifelse((((data_expanded_final_new$x - 50)^2 / 25^2) + ((data_expanded_final_new$y - 150)^2 / 50^2)) < 1 & data_expanded_final_new$t < 10, 1, 0)
time_indicator <- ifelse(data_expanded_final_new$t < 10, 1, 0)

data_expanded_final_new$ventral <- ventral_indicator
data_expanded_final_new$time_ind <- time_indicator

#plots
ggplot(data_expanded_final_new %>% filter(ventral == 1), aes(x=x,y=y,color=t)) +
  geom_point() + scale_color_viridis_c()
ggplot((data_expanded_final_new %>% filter(time_ind == 1)), aes(x=x,y=y,color=t)) +
  geom_point() + scale_color_viridis_c()

set.seed(470)
split <- sample(c(TRUE, FALSE), nrow(data_expanded_final_new), replace=TRUE, prob=c(0.7,0.3))
data_expanded_train_new <- data_expanded_final_new[split,]
data_expanded_test_new <- data_expanded_final_new[!split,]


# fit model
st <- Sys.time()
updated_log_reg <- glm(intensity ~ x + y + t + distance + angle(x=x,y=y) + 
                         radial(x=x,y=y) + angle(x=x,y=y):t + radial(x=x,y=y):t +
                         ventral*(x + y + t + distance + angle(x=x,y=y) + 
                                    radial(x=x,y=y) + angle(x=x,y=y):t + radial(x=x,y=y):t), 
                       data=data_expanded_new, family="binomial")
print(Sys.time()-st)


p2 <- ggplot(sample_n(data_expanded_final_new %>% filter(intensity==1),10000), aes(x = x, y = y, color = t)) + geom_point(alpha=0.3) + 
  scale_color_viridis_c() + labs(title="Sample")

data_test_preds <- data_expanded_test_new
preds <- predict(updated_log_reg, data_test_preds, type="response")
data_preds_all <- data_test_preds %>% mutate(preds = as.numeric(ifelse(preds < mean(preds), 0, 1)))

p4 <- ggplot(data=sample_n((data_preds_all %>%  filter(preds == 1)), length(moonbeam$t)), aes(x = x, y = y, color = t)) + geom_point(alpha=0.3) + 
  scale_color_viridis_c() + labs(title="New")

grid.arrange(p4,p2)

summary(updated_log_reg)
summary(logistic_all_new)


summary(preds)





#UPDATED LOGISITC REGRESSIONS


### logistic regression predictions at each time point expanded

ventral_indicator <- ifelse((((data_final_new$x - 50)^2 / 25^2) + ((data_final_new$y - 150)^2 / 50^2)) < 1 & data_final_new$t < 10, 1, 0)
time_indicator <- ifelse(data_final_new$t < 10, 1, 0)

data_final_new$ventral <- ventral_indicator
data_final_new$time_ind <- time_indicator

updated_log_reg_ventral <- glm(intensity ~ x + y + t + angle(x=x,y=y) + 
                         radial(x=x,y=y) + angle(x=x,y=y):t + radial(x=x,y=y):t +
                         ventral*(x + y + t + angle(x=x,y=y) + 
                                    radial(x=x,y=y) + angle(x=x,y=y):t + radial(x=x,y=y):t), 
                       data=data_expanded_train_new, family="binomial")

summary(updated_log_reg_ventral)



unique_coords_v <- unique(data.frame(x = data_final_new$x, y= data_final_new$y,
                                     ventral = data_final_new$ventral))
test <- unique_coords_v
unique_coords_v <- test

data <- data.frame()
for(i in 2:timesteps){
  
  time_data <- unique_coords_v %>% mutate(t = i)
  
  expanded_time_data <- time_data
  expanded_time_data$x <- expanded_time_data$x - hull_center[1]
  expanded_time_data$y <- expanded_time_data$y - hull_center[2]
  
  expand_x <- expanded_time_data$x * xradius / major_radius[i]
  expand_y <- expanded_time_data$y * yradius / minor_radius[i]
  
  expanded_time_data$x <- expand_x + hull_center[1]
  expanded_time_data$y <- expand_y + hull_center[2]
  
  #predictions
  preds <- predict(updated_log_reg_ventral, expanded_time_data, type="response")
  
  #finding optimal cutoff value
  error <- rep(NA,100)
  for(j in 1:100){
    preds_new <- as.numeric(ifelse(preds < j/100, 0, 1))
    #vector indicating which cells/observations are predicted to be "on"
    if(sum(preds_new) > 0){
      points <- data.frame(x = time_data$x[which(preds_new == 1)], 
                           y = time_data$y[which(preds_new == 1)], 
                           t = time_data$t[which(preds_new == 1)])
      
      error[j] <- mean((points$t - i)^2)
      #creating some level of error in predictions, somewhat resembling a mean squared error
    }
  }
  best_cutoff <- which.min(error)/100
  
  preds_new <- as.numeric(ifelse(preds < best_cutoff, 0, 1))
  #predictions using best cutoff value
  
  data_preds <- data.frame(x = time_data$x[which(preds_new == 1)], y = time_data$y[which(preds_new == 1)])
  
  if(dim(data_preds)[1]!=0){
    data_preds <- data_preds %>% mutate(t = i)
    data <- rbind(data, data_preds)
    unique_coords_v <- data.frame(x = time_data$x[which(preds_new == 0)], 
                                  y = time_data$y[which(preds_new == 0)],
                                  ventral = time_data$ventral[which(preds_new == 0)])
  }
}

ggplot(sample_n(data,5000), aes(x=x,y=y,color=t))+geom_point()+scale_color_viridis_c(limits=c(0,timesteps))



#%%% distance component added
# doesn't seem to work well

# adding distance into this, for each point that is predicted to be on, calculate the minimum distance for 
# the remaining points

dist_v <- function(data, data_preds){
    distance <- vector("numeric", length(data$x))
    for (j in 1:length(data$x)) {
      distance[j] <- min(euclidean_dist(x = data$x[j], y = data$y[j],
                                        x0 = data_preds$x, y0 = data_preds$y))
    }

  return(distance)
}

#selecting distance, but will be updated through iterations
#only need distances for first timepoint
unique_coords_v <- unique(data.frame(x = data_final_new$x, y= data_final_new$y,
                                     ventral = data_final_new$ventral,
                                     distance = data_final_new$distance))

updated_log_reg_distance <- glm(intensity ~ x + y + t + distance + angle(x=x,y=y) + 
                                 radial(x=x,y=y) + angle(x=x,y=y):t + radial(x=x,y=y):t +
                                 ventral*(x + y + t + distance + angle(x=x,y=y) + 
                                            radial(x=x,y=y) + angle(x=x,y=y):t + radial(x=x,y=y):t), 
                               data=data_expanded_train_new, family="binomial")

st <- Sys.time()
data <- data.frame()
for(i in 4:timesteps){
  
  time_data <- unique_coords_v %>% mutate(t = i)
  
  expanded_time_data <- time_data
  expanded_time_data$x <- expanded_time_data$x - hull_center[1]
  expanded_time_data$y <- expanded_time_data$y - hull_center[2]
  
  expand_x <- expanded_time_data$x * xradius / major_radius[i]
  expand_y <- expanded_time_data$y * yradius / minor_radius[i]
  
  expanded_time_data$x <- expand_x + hull_center[1]
  expanded_time_data$y <- expand_y + hull_center[2]
  
  #predictions
  preds <- predict(updated_log_reg_distance, expanded_time_data, type="response")
  
  #finding optimal cutoff value
  error <- rep(NA,100)
  for(j in 1:100){
    preds_new <- as.numeric(ifelse(preds < j/100, 0, 1))
    #vector indicating which cells/observations are predicted to be "on"
    if(sum(preds_new) > 0){
      points <- data.frame(x = time_data$x[which(preds_new == 1)], 
                           y = time_data$y[which(preds_new == 1)], 
                           t = time_data$t[which(preds_new == 1)])
      
      error[j] <- mean((points$t - i)^2)
      #creating some level of error in predictions, somewhat resembling a mean squared error
    }
  }
  best_cutoff <- which.min(error)/100
  
  preds_new <- as.numeric(ifelse(preds < best_cutoff, 0, 1))
  #predictions using best cutoff value
  
  data_preds <- data.frame(x = time_data$x[which(preds_new == 1)],
                           y = time_data$y[which(preds_new == 1)])
  
  if(dim(data_preds)[1]!=0){
    data_preds <- data_preds %>% mutate(t = i)
    data <- rbind(data, data_preds)
    if(i < timesteps){
      indicator <- which(unique_coords_v$x %in% data_preds$x & unique_coords_v$y %in% data_preds$y)
      x <- unique_coords_v$x[-1*indicator]
      y <- unique_coords_v$y[-1*indicator]
      ventral <- unique_coords_v$ventral[-1*indicator]
      distance <- dist_v(unique_coords_v, data_preds)[-1*indicator]
      unique_coords_v <- data.frame(x=x, y=y, ventral=ventral, distance=distance)
    }
  } else{
    unique_coords_v <- time_data
  }
}
print(Sys.time()-st)

ggplot(data, aes(x=x,y=y,color=t))+geom_point()+scale_color_viridis_c()


# survival analysis
# kaplan meier estimator of survival curve
library(survival)
library(ggfortify)

kaplan_meier <- survfit(Surv(t, intensity) ~ ventral, data=sample_n(data_final_new %>% filter(intensity==1),10000))
# filtering because 1 indicates *when* the subject (cell) expressed atoh7, it's off before that time
autoplot(kaplan_meier)
summary(kaplan_meier)

summary(data_final_new$intensity)
