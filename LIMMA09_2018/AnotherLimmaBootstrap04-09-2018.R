library(reshape2)
library(data.table)
library(dplyr)
library(limma)
library(splines)
library(parallel)
library(combinat)

## We are primarily interest in three questions:
##    1) Does gene g respond to TGFb? - Compare TGFb vs Control groups
##    2) Is the gene g response different in adult/senescent compared to neonatal cell lines? - Compare Control groups for cell lines 
##    3) Is the gene g response to TGFb different in adult/senescent compared to neonatal cell lines? - Compare treated groups for cell lines


## normalize for the reference genes. 
calc_dct2 = function(df){
  df = data.frame(t(df))
  g=df$PPIA
  return(2^(-(df-g)))
}


do_limma_q1 = function(dct, comparison = 'adult'){
  ## this function uses limma to answer the question: 
  ##    - does gene g respond to TGFb in any of the cell lines? 
  
  df = dct[which(dct$treatment != 'Baseline'),]
  
  ## create new column name format treat
  new_names = paste(df$cell_line, df$treatment, df$time, sep='_')

  for (i in c('time', 'cell_line', 'replicate', 'treatment')){
    df[,i] = NULL
  }
  df = t(df)
  colnames(df) = new_names

    ## create pd matrix
  pd = data.frame(matrix(unlist(strsplit(new_names, split='_')), nrow=length(new_names), byrow=T))
  colnames(pd) =c('cell_line', 'treatment', 'time')
  pd$time = as.numeric(as.character(pd$time))
  pd$cell_treat = paste(pd$cell_line, pd$treatment, sep='_')
  pd$cell_treat = as.factor(pd$cell_treat)
 
  X = ns(pd$time, df=4)

  design = model.matrix(~0 + pd$cell_treat*X)
  print(head(design))
  print(dim(design))

  fit = lmFit(df, design)

  fit = eBayes(fit)
  top = topTable(fit, number = 100, coef=23:90)
  
  fname = file.path(saved_objects_path, 'TopTableForQuestion1.csv')
  write.csv(top, file=fname)
  return (top)
  
}


do_limma = function(dct, n_cell_line, other_cell_line, treatment='Control'){

  
  dct = dct[which(dct$treatment != 'Baseline'),]
  
  df = dct[which(dct$cell_line == n_cell_line | dct$cell_line == other_cell_line),]
  
  if (treatment == 'Control'){
    df = df[df$treatment == 'Control',]
  }
  else if (treatment == 'TGFb') {
    df = df[df$treatment == 'TGFb', ]
  }

  new_names = paste(df$cell_line, df$treatment, df$time, sep='_')
  
  for (i in c('time', 'cell_line', 'replicate', 'treatment')){
    df[,i] = NULL
  }
  df = t(df)
  colnames(df) = new_names
  
  # ## create pd matrix
  pd = data.frame(matrix(unlist(strsplit(new_names, split='_')), nrow=length(new_names), byrow=T))
  colnames(pd) =c('cell_line', 'treatment', 'time')
  pd$time = as.numeric(as.character(pd$time))

  X = ns(pd$time, df=4)
  
  design = model.matrix(~0 + pd$cell_line*X)

  fit = lmFit(df, design)

  fit = eBayes(fit)
  top = topTable(fit, number = 100, coef=7:10)

  # fname = file.path(saved_objects_path, 'TopTableForQuestion1.csv')
  # write.csv(top, file=fname)
  return (top)
  
}

limma_combinations_between_groups = function(dct, treatment='Control'){
  adult = list()
  sen = list()
  for (i in c('A', 'B', 'C')){
    for (j in c('D', 'E', 'F') ){

      current = paste(i, j, sep='_')
      sen[[current]] = do_limma(dct, n_cell_line = i, other_cell_line = j, treatment = treatment)
    }
  }
  for (i in c('A', 'B', 'C')){
    for (j in c('G', 'H', 'I') ){

      current = paste(i, j, sep='_')
      adult[[current]] = do_limma(dct, n_cell_line = i, other_cell_line = j, treatment = treatment)
    }
  }
  q = list(adult=adult, sen=sen)
  return(q)
}

limma_combinations_within_groups = function(dct, cell_line='neonatal'){
  if (cell_line == 'neonatal') {
    l = list(list('A', 'B'), list('A', 'C'), list('B', 'C'))
  } else if (cell_line == 'senescent'){
    l = list(list('D', 'E'), list('D', 'F'), list('E', 'F'))
  } else if (cell_line == 'adult'){
    l = list(list('G', 'H'), list('G', 'I'), list('H', 'I'))
  }
  
  tgf = list()
  control = list()
  
  for (i in l){
    current = paste(i[[1]], i[[2]], sep='_')
    tgf[[current]] = do_limma(dct, n_cell_line = i[[1]], other_cell_line = i[[2]], treatment = 'TGFb')
    
    control[[current]] = do_limma(dct, n_cell_line = i[[1]], other_cell_line = i[[2]], treatment = 'Control')
      # write.csv(tgf[[current]])
    }
  return(list(TGFb = tgf, Control = control))
}

## Take one element (second level) of the list produced by do_limma
## if adjusted p val is < x, 1 else 0. 
check_pval = function(limma_output, pval=0.05){
  vec = vector()
  for (i in rownames(limma_output)){
    if (is.na(limma_output[[i, 'adj.P.Val']]) ){
      vec = c(vec, 0)
    }
    else if (limma_output[[i, 'adj.P.Val']] < pval) {
      vec = c(vec, 1)
    }
    else {
      vec = c(vec, 0)
    }
  }
  names(vec) = as.factor(rownames(limma_output))
  
  df = data.frame(vec)
  df[, 'gene'] = rownames(df)
  df = df[order(df$gene),]
  df$gene = NULL
  return(df)
}



calculate_percentages_between_groups = function(l, pval=0.05) {
  ad = vector()
  sen = vector()
  for (i in names(l$adult)){
      ad[[i]] = check_pval(l$adult[[i]], pval)
  }
  
  for (i in names(l$sen)){
    sen[[i]] = check_pval(l$sen[[i]], pval)
  }
  ad_df = do.call(cbind, ad)
  colnames(ad_df) = names(ad)

  sen_df = do.call(cbind, sen)
  colnames(sen_df) = names(sen)
  
  sen_perc = rowSums(sen_df) / 9 * 100
  ad_perc = rowSums(ad_df) / 9 * 100
  return (list(adult=ad_df, sen=sen_df, sen_perc=sen_perc, ad_perc=ad_perc))
}


calculate_percentages_within_groups = function(l, pval=0.05) {
  tgfb = list()
  ctrl = list()
  for (i in names(l$TGFb)){
    tgfb[[i]] = check_pval(l$TGFb[[i]], pval)
  }
  
  for (i in names(l$Control)){
    ctrl[[i]] = check_pval(l$Control[[i]], pval)
  }
  
  tgfb_df = do.call(cbind, tgfb)
  colnames(tgfb_df) = names(tgfb)
  print(tgfb_df)
  ctrl_df = do.call(cbind, ctrl)
  colnames(ctrl_df) = names(ctrl)

  tgfb_perc = rowSums(tgfb_df) / 3 * 100
  ctrl_perc = rowSums(ctrl_df) / 3 * 100
  return (list(tgfb_perc = tgfb_perc, ctrl_perc = ctrl_perc))
}


##############################
##  main section

## set some printing options
options(scipen = 999)

## sort out directories. Change the below directory variable to appropriate location on your computer
directory = '/home/b3053674/Documents/LargeStudy/LIMMA09-2018'
saved_objects_path = file.path(directory, 'SavedObjects')
raw_data_file = file.path(saved_objects_path, 'FullDataFrameRawCT_16_003_2018.csv')

## read desired pval from settings file
pval_path = file.path(directory, 'pval')
PVAL = read.csv(pval_path, header = F)$V1

pval_dir = file.path(saved_objects_path, paste0('pval_less_than_', gsub('\\.', '_', as.character(PVAL))  ) )

pval_dir
dir.create(pval_dir, showWarnings = FALSE)

# PVAL = 0.0001

## read data
df = read.csv(raw_data_file)

## rearrange the data
df = dcast(df, formula = Assay ~ Sample, value.var='Ct')
rownames(df) = df$Assay
df$Assay = NULL
dim(df)


dct = calc_dct2(df)

## sanity check
sanity = df[,'Baseline_0_A_1']
names(sanity) = rownames(df)
2 ^ (-(sanity['ACTA2'] - sanity['PPIA'])) #== 0.2702207 

dct['Baseline_0_A_1','ACTA2'] #== 0.2702207


## add experimental factors  
exp_factors = data.frame(matrix(unlist(strsplit(as.character(rownames(dct)), '_')), nrow = dim(dct)[1], byrow=T))
colnames(exp_factors) = c('treatment', 'time', 'cell_line', 'replicate')
exp_factors$time = as.numeric(as.character(exp_factors$time))
exp_factors$replicate = as.numeric(as.character(exp_factors$replicate))
##concatonate factors into df
dct = cbind.data.frame(dct, exp_factors)
##convert back to data.frame
dct = data.frame(dct)
dct
############
## between code

## do limma for q2 and q2. The only difference is whether the treatment or control data is used.
## limma_combinations does limma for all combinations of adult/senescent with neonatal cell lines. 

between_control = limma_combinations_between_groups(dct, treatment='Control')

between_tgfb = limma_combinations_between_groups(dct, treatment='TGFb')


## concationate frames and write to file
# do.call(cbind, calculate_percentages_between_groups(between_control))

# do.call(cbind, calculate_percentages_between_groups(between_tgfb))

write.csv(
  do.call(cbind, calculate_percentages_between_groups(between_control, pval = PVAL)), 
  file.path(pval_dir, 'between_control_statistics.csv')
)

write.csv(
  do.call(cbind, calculate_percentages_between_groups(between_tgfb, pval = PVAL)), 
  file.path(pval_dir, 'between_tgfb_statistics.csv')
)




###############
## 'within' code


neo = limma_combinations_within_groups(dct, cell_line='neonatal')

adult = limma_combinations_within_groups(dct, cell_line = 'adult')

sen = limma_combinations_within_groups(dct, cell_line = 'senescent')


neo
## concationate frames and write to file
do.call(cbind, calculate_percentages_within_groups(neo, pval = PVAL))

do.call(cbind, calculate_percentages_within_groups(adult, pval = PVAL))

write.csv(
  do.call(cbind, calculate_percentages_within_groups(neo, pval = PVAL)), 
  file.path(pval_dir, 'within_neonatal.csv')
)

write.csv(
  do.call(cbind, calculate_percentages_within_groups(adult, pval = PVAL)), 
  file.path(pval_dir, 'within_adult.csv')
)

write.csv(
  do.call(cbind, calculate_percentages_within_groups(sen, pval = PVAL)), 
  file.path(pval_dir, 'within_senescent.csv')
)












