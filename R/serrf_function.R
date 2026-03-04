library(shiny)
library(data.table)
library(parallel)
library(ranger)
library(openxlsx)
library(pacman)
library(bootstrap)
library(shinyjs)
library(shiny)


serrfR = function(train = e[,p$sampleType == 'qc'],
target = e[,p$sampleType == 'sample'],
num = 10,
batch. = factor(c(batch[p$sampleType=='qc'],batch[p$sampleType=='sample'])),
time. = c(time[p$sampleType=='qc'],time[p$sampleType=='sample']),
sampleType. = c(p$sampleType[p$sampleType=='qc'],p$sampleType[p$sampleType=='sample']),minus=minus,cl){


all = cbind(train, target)
normalized = rep(0, ncol(all))

for(j in 1:nrow(all)){
for(b in 1:length(unique(batch.))){
current_batch = levels(batch.)[b]
# fixing zeros.
all[j,batch.%in%current_batch][all[j,batch.%in%current_batch] == 0] = rnorm(length(all[j,batch.%in%current_batch][all[j,batch.%in%current_batch] == 0]),mean = min(all[j,batch.%in%current_batch][!is.na(all[j,batch.%in%current_batch])])+1,sd = 0.1*(min(all[j,batch.%in%current_batch][!is.na(all[j,batch.%in%current_batch])])+.1))
# fixing NAs.
all[j,batch.%in%current_batch][is.na(all[j,batch.%in%current_batch])] = rnorm(length(all[j,batch.%in%current_batch][is.na(all[j,batch.%in%current_batch])]),mean = 0.5*min(all[j,batch.%in%current_batch][!is.na(all[j,batch.%in%current_batch])])+1,sd = 0.1*(min(all[j,batch.%in%current_batch][!is.na(all[j,batch.%in%current_batch])])+.1))
}
}

corrs_train = list()
corrs_target = list()
for(b in 1:length(unique(batch.))){

current_batch = levels(batch.)[b]

train_scale = t(apply(train[,batch.[sampleType.=='qc']%in%current_batch],1,scale))
if(is.null(dim(target[,batch.[!sampleType.=='qc']%in%current_batch]))){ # !!!
target_scale = scale(target[,batch.[!sampleType.=='qc']%in%current_batch])#!!!
}else{
# target_scale = scale(target[,batch.[!sampleType.=='qc']%in%current_batch])#!!!
target_scale = t(apply(target[,batch.[!sampleType.=='qc']%in%current_batch],1,scale))
}

# all_scale = cbind(train_scale, target_scale)

# e_current_batch = all_scale
corrs_train[[current_batch]] = cor(t(train_scale), method = "spearman")
corrs_target[[current_batch]] = cor(t(target_scale), method = "spearman")
# corrs[[current_batch]][is.na(corrs[[current_batch]])] = 0
}



pred = matrix(nrow = nrow(all), ncol = length(sampleType.))

# cs = c()
withProgress(message = "Normalization in Progress.", value = 0, { # HERE

for(j in 1:nrow(all)){
incProgress(1/nrow(all), detail = paste("Working on compound", j,"/", nrow(all))) # HERE
# print(j) # HERE
# if(j == 2){
#   stop("fine")
# }
# j = j+1
# print(j)
normalized  = rep(0, ncol(all))
qc_train_value = list()
qc_predict_value = list()
sample_value = list()
sample_predict_value = list()

for(b in 1:length(levels(batch.))){
current_batch = levels(batch.)[b]
current_time = time.[batch. %in% current_batch]
e_current_batch = all[,batch.%in%current_batch]
corr_train = corrs_train[[current_batch]]
corr_target = corrs_target[[current_batch]]


corr_train_order = order(abs(corr_train[,j]),decreasing = TRUE) #!!! try using corr_target_order[1:num]
corr_target_order = order(abs(corr_target[,j]),decreasing = TRUE)

sel_var = c()
l = num
# print(intersect(corr_train_order[1:l], corr_target_order[1:l]))
while(length(sel_var)<(num)){
sel_var = intersect(corr_train_order[1:l], corr_target_order[1:l])
sel_var = sel_var[!sel_var == j]
l = l+1
}



train.index_current_batch = sampleType.[batch.%in%current_batch]
# train_data_y = scale(e_current_batch[j, train.index_current_batch=='qc'],scale=F) #!!! trying to use different scale. Scale to the test_y

# remove_outlier(e_current_batch[j, train.index_current_batch=='qc'])

factor = sd(e_current_batch[j, train.index_current_batch=='qc'])/sd(e_current_batch[j, !train.index_current_batch=='qc'])
if(factor==0 | is.nan(factor) | factor<1 | is.na(factor)){#!!!
train_data_y = scale(e_current_batch[j, train.index_current_batch=='qc'],scale=F) 
}else{
# print(j)
# print("!!")
if(sum(train.index_current_batch=='qc')*2>=sum(!train.index_current_batch=='qc')){
train_data_y = (e_current_batch[j, train.index_current_batch=='qc'] - mean(e_current_batch[j, train.index_current_batch=='qc']))/factor ### need to be careful with outlier!
}else{
train_data_y = scale(e_current_batch[j, train.index_current_batch=='qc'],scale=F)
}
}


train_data_x = apply(e_current_batch[sel_var, train.index_current_batch=='qc'],1,scale)

if(is.null(dim(e_current_batch[sel_var, !train.index_current_batch=='qc']))){
test_data_x = t(scale(e_current_batch[sel_var, !train.index_current_batch=='qc']))
}else{
test_data_x = apply(e_current_batch[sel_var, !train.index_current_batch=='qc'],1,scale)
}

train_NA_index  = apply(train_data_x,2,function(x){
sum(is.na(x))>0
})

train_data_x = train_data_x[,!train_NA_index]
test_data_x = test_data_x[,!train_NA_index]

if(!"matrix" %in% class(test_data_x)){ # !!!
test_data_x = t(test_data_x)
}

good_column = apply(train_data_x,2,function(x){sum(is.na(x))==0}) & apply(test_data_x,2,function(x){sum(is.na(x))==0})
train_data_x = train_data_x[,good_column]
test_data_x = test_data_x[,good_column]
if(!"matrix" %in% class(test_data_x)){ # !!!
test_data_x = t(test_data_x)
}
train_data = data.frame(y = train_data_y,train_data_x )

if(ncol(train_data)==1){# some samples have all QC constent.
norm = e_current_batch[j,]
normalized[batch.%in%current_batch] = norm
}else{
colnames(train_data) = c("y", paste0("V",1:(ncol(train_data)-1)))

set.seed(1)
model = ranger(y~., data = train_data)

test_data = data.frame(test_data_x)
colnames(test_data) = colnames(train_data)[-1]

norm = e_current_batch[j,]



# plot(e_current_batch[j, train.index_current_batch=='qc'], ylim = c(min(c(e_current_batch[j, train.index_current_batch=='qc'],e_current_batch[j, !train.index_current_batch=='qc'])),max(c(e_current_batch[j, train.index_current_batch=='qc'],e_current_batch[j, !train.index_current_batch=='qc']))))
# points(predict(model, data = train_data)$prediction+mean(e_current_batch[j,train.index_current_batch=='qc'],na.rm=TRUE), col = 'red')
# plot(e_current_batch[j,!train.index_current_batch=='qc'], ylim = c(min(c(e_current_batch[j, train.index_current_batch=='qc'],e_current_batch[j, !train.index_current_batch=='qc'])),max(c(e_current_batch[j, train.index_current_batch=='qc'],e_current_batch[j, !train.index_current_batch=='qc']))))
# points(predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE) - mean(predict(model,data = test_data)$predictions), col = 'red')
# 
# plot(y =  e_current_batch[j,], x = current_time, col = factor(train.index_current_batch), ylim = range(c(e_current_batch[j,], norm)))
# plot(y = norm, x = current_time, col = factor(train.index_current_batch), ylim = range(c(e_current_batch[j,], norm)))





if(minus){


norm[train.index_current_batch=='qc'] = e_current_batch[j, train.index_current_batch=='qc']-((predict(model, data = train_data)$prediction+mean(e_current_batch[j,train.index_current_batch=='qc'],na.rm=TRUE))-mean(all[j,sampleType.=='qc'],na.rm=TRUE))

norm[!train.index_current_batch=='qc'] = e_current_batch[j,!train.index_current_batch=='qc']-((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE))-(median(all[j,!sampleType.=='qc'],na.rm = TRUE)))


}else{
norm[train.index_current_batch=='qc'] = e_current_batch[j, train.index_current_batch=='qc']/((predict(model, data = train_data)$prediction+mean(e_current_batch[j,train.index_current_batch=='qc'],na.rm=TRUE))/mean(all[j,sampleType.=='qc'],na.rm=TRUE))

# norm[!train.index_current_batch=='qc'] = e_current_batch[j,!train.index_current_batch=='qc']/((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE))/(median(all[j,!sampleType.=='qc'],na.rm = TRUE)))


norm[!train.index_current_batch=='qc'] = e_current_batch[j,!train.index_current_batch=='qc']/((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE)- mean(predict(model,data = test_data)$predictions))/(median(all[j,!sampleType.=='qc'],na.rm = TRUE)))




# norm[!train.index_current_batch=='qc'] = e_current_batch[j,!train.index_current_batch=='qc']/((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE))/(median(e_current_batch[j,!train.index_current_batch=='qc'],na.rm = TRUE)))

}

# norm[!train.index_current_batch=='qc'] =(e_current_batch[j,!train.index_current_batch=='qc'])/((predict(model, data = test_data)$prediction + mean(e_current_batch[j,!train.index_current_batch=='qc'],na.rm=TRUE))/mean(e_current_batch[j,!train.index_current_batch=='qc'],na.rm=TRUE))


norm[!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']<0]=e_current_batch[j,!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']<0] # fix negative value


# plot(p$time[batch.%in%b][!train.index_current_batch=='qc'], (e_current_batch[j,!train.index_current_batch=='qc'])/((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, train.index_current_batch=='qc'],na.rm=TRUE))/(median(e_current_batch[j,!train.index_current_batch=='qc'],na.rm = TRUE))))


norm[train.index_current_batch=='qc'] = norm[train.index_current_batch=='qc']/(median(norm[train.index_current_batch=='qc'],na.rm=TRUE)/median(all[j,sampleType.=='qc'],na.rm=TRUE)) #!!! putting all to the same batch level.


norm[!train.index_current_batch=='qc'] = norm[!train.index_current_batch=='qc']/(median(norm[!train.index_current_batch=='qc'],na.rm=TRUE)/median(all[j,!sampleType.=='qc'],na.rm=TRUE))


norm[!is.finite(norm)] = rnorm(length(norm[!is.finite(norm)]),sd = sd(norm[is.finite(norm)],na.rm=TRUE)*0.01)# fix infinite




out = boxplot.stats(norm, coef = 3)$out
attempt = ((e_current_batch[j,!train.index_current_batch=='qc'])-((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE))-(median(all[j,!sampleType.=='qc'],na.rm = TRUE))))[norm[!train.index_current_batch=='qc']%in%out];
if(length(out)>0 & length(attempt)>0){
if(mean(out)>mean(norm)){
if(mean(attempt)<mean(out)){
norm[!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']%in%out] =  attempt# !!! this may not help deal with outlier effect..
}
}else{
if(mean(attempt)>mean(out)){
norm[!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']%in%out] =  attempt# !!! this may not help deal with outlier effect..
}
}
}




norm[!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']<0]=e_current_batch[j,!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']<0]


normalized[batch.%in%current_batch] = norm


}


}



c = (median(normalized[sampleType.=="sample"])+(median(all[j,sampleType.=="qc"])-median(all[j,!sampleType.=="qc"]))/sd(all[j,!sampleType.=="qc"]) * sd(normalized[sampleType.=="sample"]))/median(normalized[!sampleType.=="sample"])
normalized[sampleType.=="qc"] = normalized[sampleType.=="qc"] * ifelse(c>0,c,1)
# cs[j] = c

pred[j,] = normalized
# print(j)
}#!!!!


}) # HERE










normed = pred

normed_target = normed[,!sampleType.=='qc']


for(i in 1:nrow(normed_target)){ # fix NA
normed_target[i,is.na(normed_target[i,])] = rnorm(sum(is.na(normed_target[i,])), mean = min(normed_target[i,!is.na(normed_target[i,])], na.rm = TRUE), sd = sd(normed_target[i,!is.na(normed_target[i,])])*0.1)
}
for(i in 1:nrow(normed_target)){ # fix negative value
normed_target[i,normed_target[i,]<0] = runif(1) * min(normed_target[i,normed_target[i,]>0], na.rm = TRUE)
}


normed_train = normed[,sampleType.=='qc']


for(i in 1:nrow(normed_train)){ # fix NA
normed_train[i,is.na(normed_train[i,])] = rnorm(sum(is.na(normed_train[i,])), mean = min(normed_train[i,!is.na(normed_train[i,])], na.rm = TRUE), sd = sd(normed_train[i,!is.na(normed_train[i,])])*0.1)
}
for(i in 1:nrow(normed_train)){ # fix negative value
normed_train[i,normed_train[i,]<0] = runif(1) * min(normed_train[i,normed_train[i,]>0], na.rm = TRUE)
}
return(list(normed_train=normed_train,normed_target=normed_target))
}

# serrf_normalized = e
# train = e[,p$sampleType == 'qc']
# target = e[,p$sampleType == 'sample']
# batch. = factor(c(batch[p$sampleType=='qc'],batch[p$sampleType=='sample']))
# time. = c(time[p$sampleType=='qc'],time[p$sampleType=='sample'])
# sampleType. = c(p$sampleType[p$sampleType=='qc'],p$sampleType[p$sampleType=='sample'])

serrf_normalized = e

serrf_normalized_modeled = serrfR(train = e_qc, target = e_sample, num = 10,batch. = factor(c(p_qc$batch, p_sample$batch)),time. = c(p_qc$time, p_sample$time),sampleType. = c(p_qc$sampleType, p_sample$sampleType),minus,cl)