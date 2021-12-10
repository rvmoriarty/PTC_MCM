options(stringsAsFactors = FALSE)

if("vegan" %in% installed.packages() == F){install.packages("vegan")}
if("stringr" %in% installed.packages() == F){install.packages("stringr")}
if("readxl" %in% installed.packages() == F){install.packages("readxl")}

library(vegan)
library(stringr)
library(readxl)

to.df = function(longlist){
  df = data.frame()
  for(item in 1:length(longlist)){
    df = rbind(df, longlist[[item]])
  }
  return(df)
}
read_excel_allsheets = function(filename, tibble = FALSE){
  sheets =readxl::excel_sheets(filename)
  x = lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x = lapply(x, as.data.frame)
  names(x) = sheets
  x
}

setwd("/Users/ryanmoriarty/Documents/BarcodeAnalysisTool/analysis_outputs")

cyno_n803 = c("cy1035", "cy1036", "cy1039", "cy1043")
cyno_cont = c("cy1037", "cy1040", "cy1044", "cy1045")

animals = c(cyno_cont, cyno_n803)

analysis_outputs = list.files(pattern = "Analysis.xlsx")

analysis_output_xlsx = list()
analysis_output_xlsx_names = c()
for(i in 1:length(analysis_outputs)){
  analysis_output_xlsx[[i]] = read_excel_allsheets(analysis_outputs[i])
  analysis_output_xlsx_names[i] = str_split(analysis_outputs[i], " ")[[1]][2]
}
analysis_output_xlsx_names = unlist(analysis_output_xlsx_names)
names(analysis_output_xlsx) = analysis_output_xlsx_names

analysis_output_samples = list()
for (a in 1:length(analysis_output_xlsx)) {
  analysis_output_samples[[a]] = analysis_output_xlsx[[a]][grep("samples", names(analysis_output_xlsx[[a]]))]
}
names(analysis_output_samples) = names(analysis_output_xlsx)

sample_sdis = list()
sample_sdi_dfs = data.frame()
d = 1 
for(a in 1:length(analysis_output_samples)){
  animal_dfs = analysis_output_samples[[a]]
  for(b in 1:length(animal_dfs)){
    data = animal_dfs[[b]]
    count_cols = grep("Counts",data)
    for(c in 1:length(count_cols)){
      sample_sdi_dfs[d, "animal"] = gsub(" samples","", names(animal_dfs)[b])
      sample_sdi_dfs[d, "sample_date"] = colnames(data)[count_cols[c]]
      sample_sdi_dfs[d, "number_sequences"] = as.numeric(data[2, count_cols[c]-1])
      sample_sdi_dfs[d, "number_barcodes_found"] = as.numeric(data[2, count_cols[c]])
      sample_sdi_dfs[d, "simpsdi"] = diversity(as.numeric(data[5:data[2, count_cols[c]], count_cols[c]]), "simpson")
      sample_sdi_dfs[d, "shandi"] = diversity(as.numeric(data[5:data[2, count_cols[c]], count_cols[c]]), "shannon")
      d = d + 1
    }
  }
}

write.csv(number_barcodes, "number_of_barcodes.csv")
write.csv(sample_sdi_dfs, "sample_sdi_dfs.csv")

barcodes_found = data.frame()
sample_barcodes_found = list()
for(b in 1:length(barcode_csvs_list)){
  data = barcode_csvs_list[[b]][grep("SIVmac239M", barcode_csvs_list[[b]][,1]), c(1:4)]
  data = cbind(new_names[b], data)
  data$animal = unlist(str_split(barcode_csvs_names[b], "_"))[1]
  data$dpi = as.numeric(gsub("day", "", unlist(str_split(barcode_csvs_names[b], "_"))[2]))
  data$rep = unlist(str_split(barcode_csvs_names[b], "_"))[3]
  colnames(data) = c("sample", "barcode", "sequence", "count", "proportion", "animal", "dpi", "rep")
  barcodes_found = rbind(barcodes_found, data)
  sample_barcodes_found[[b]] = data
}
names(sample_barcodes_found) = new_names

by_animal = list()
for(a in 1:length(animals)){
  by_animal[[a]] = barcodes_found[grep(animals[a], barcodes_found$sample), ]
}
names(by_animal) = animals

samples = sort(unique(barcodes_found$sample))
sample_list = list()
for(i in 1:length(samples)){
  sample_list[[i]] = barcodes_found[which(barcodes_found$sample == samples[i]), ]
}
names(sample_list) = samples

for(i in 1:length(sample_list)){
  for(j in 1:nrow(sample_list[[i]])){
    sample_list[[i]][j, "proportion"] = as.numeric(sample_list[[i]][j, "count"])/sum(as.numeric(sample_list[[i]][,"count"]))
  }
}

sample_sdis = data.frame()
for(i in 1:length(sample_list)){
  sample_sdis[i, "sample"] = samples[i]
  sample_sdis[i, "sdi"] = diversity(as.numeric(sample_list[[i]][,"count"]), index = "simpson")
  sample_sdis[i, "total_sequences"] = number_barcodes[i, "number_sequences"]
  sample_sdis[i, "animal"] = unlist(str_split(samples[i], "_"))[1]
  sample_sdis[i, "dpi"] = as.numeric(gsub("day", "", unlist(str_split(samples[i], "_"))[2]))
  sample_sdis[i, "rep"] = unlist(str_split(samples[i], "_"))[3]
}
write.csv(sample_sdis, "sdi_to_date.csv")

barcode_1pc_cutoff = list()
for(s in 1:length(sample_list)){
  sample_list[[s]]$proportion = as.numeric(sample_list[[s]]$proportion)
  data = sample_list[[s]][which(sample_list[[s]]$proportion > 0.01), ]
  otherrow = data.frame(names(sample_list)[s], "other", "other", sum(as.numeric(sample_list[[s]]$count)) - sum(as.numeric(data$count)), 1-sum(data$proportion), "animal", "dpi", "rep")
  colnames(otherrow) = colnames(data)
  barcode_1pc_cutoff[[s]] = rbind(data, otherrow)
}
names(barcode_1pc_cutoff) = names(sample_list)

barcode_3pc_cutoff = list()
for(s in 1:length(sample_list)){
  sample_list[[s]]$proportion = as.numeric(sample_list[[s]]$proportion)
  data = sample_list[[s]][which(sample_list[[s]]$proportion > 0.03), ]
  otherrow = data.frame(names(sample_list)[s], "other", "other", sum(as.numeric(sample_list[[s]]$count)) - sum(as.numeric(data$count)), 1-sum(data$proportion), "animal", "dpi", "rep")
  colnames(otherrow) = colnames(data)
  barcode_3pc_cutoff[[s]] = rbind(data, otherrow)
}
names(barcode_3pc_cutoff) = names(sample_list)

barcode_3pc_cutoff_df = to.df(barcode_3pc_cutoff)
write.csv(barcode_3pc_cutoff_df, "barcode_3pc_cutoff.csv")

by_animal_3pc_cutoff = list()
for(a in 1:length(animals)){
  by_animal_3pc_cutoff[[a]] = barcode_3pc_cutoff_df[grep(animals[a], barcode_3pc_cutoff_df$sample),]
}
names(by_animal_3pc_cutoff) = animals

number_barcodes_over_3 = data.frame()
for(a in 1:length(barcode_3pc_cutoff)){
  number_barcodes_over_3[a, "sample"] = unique(barcode_3pc_cutoff[[a]]$sample)
  number_barcodes_over_3[a, "number_over_3"] = nrow(barcode_3pc_cutoff[[a]])-1
  number_barcodes_over_3[a, "other_percent"] = barcode_3pc_cutoff[[a]][which(barcode_3pc_cutoff[[a]]$barcode == "other"), "proportion"]
  number_barcodes_over_3[a, "total_count"] = sum(as.integer(barcode_3pc_cutoff[[a]]$count))
}

barcode_1pc_cutoff_df = to.df(barcode_1pc_cutoff)

by_animal_1pc_cutoff = list()
for(a in 1:length(animals)){
  by_animal_1pc_cutoff[[a]] = barcode_1pc_cutoff_df[grep(animals[a], barcode_1pc_cutoff_df$sample),]
  by_animal_1pc_cutoff[[a]]$rep = str_extract(by_animal_1pc_cutoff[[a]]$sample, "rep1|rep2")
  by_animal_1pc_cutoff[[a]]$dpi = str_sub(by_animal_1pc_cutoff[[a]]$sample, start = 7, end = -10)
}
names(by_animal_1pc_cutoff) = animals



reorganized_byanimal_3pc_cutoff = list()
for(a in 1:length(by_animal_3pc_cutoff)){
  animaldf = by_animal_3pc_cutoff[[a]]
  dpi = unique(as.numeric(str_sub(animaldf$sample, start = 7, end = -10)))
  unique_barcodes = unique(animaldf$barcode)
  unique_barcodes = sort(unique_barcodes, decreasing = TRUE)
  newdf = data.frame()
  for(d in 1:length(dpi)){
    for(b in 1:length(unique_barcodes)){
      barcode_dpi_comb = animaldf[which(animaldf$barcode == unique_barcodes[b] & animaldf$dpi == dpi[d]), ]
      newdf[d, "animal"] = unique(str_sub(animaldf$sample, end = 5))
      newdf[d, "dpi"] = dpi[d]
      
      if("rep1" %in% barcode_dpi_comb$rep){
        newdf[d, paste0(unique_barcodes[b], "_rep1_proportion")] = barcode_dpi_comb[which(barcode_dpi_comb$rep == "rep1"), "proportion"]
        newdf[d, paste0(unique_barcodes[b], "_rep1_count")] = barcode_dpi_comb[which(barcode_dpi_comb$rep == "rep1"), "count"]
      }else{
        newdf[d, paste0(unique_barcodes[b], "_rep1_proportion")] = NA
        newdf[d, paste0(unique_barcodes[b], "_rep1_count")] = NA
      }
      if("rep2" %in% barcode_dpi_comb$rep){
        newdf[d, paste0(unique_barcodes[b], "_rep2_proportion")] = barcode_dpi_comb[which(barcode_dpi_comb$rep == "rep2"), "proportion"]
        newdf[d, paste0(unique_barcodes[b], "_rep2_count")] = barcode_dpi_comb[which(barcode_dpi_comb$rep == "rep2"), "count"]
      }else{
        newdf[d, paste0(unique_barcodes[b], "_rep2_proportion")] = NA
        newdf[d, paste0(unique_barcodes[b], "_rep2_count")] = NA
      }
    }
  }
  reorganized_byanimal_3pc_cutoff[[a]] = newdf
}
names(reorganized_byanimal_3pc_cutoff) = names(by_animal_3pc_cutoff)

by_barcode_by_animal_proportions_3pc = list()
for(a in 1:length(reorganized_byday_byanimal_3pc_cutoff)){
  by_barcode_by_animal_proportions_3pc[[a]] = reorganized_byanimal_3pc_cutoff[[a]][, c(1:3, grep("proportion", colnames(reorganized_byanimal_3pc_cutoff[[a]])))]
  write.csv(by_barcode_by_animal_proportions_3pc[[a]], paste0(names(reorganized_byanimal_3pc_cutoff)[a], "_bybarcode_reorganized_proportions.csv"))
}

reorganized_byday_byanimal_3pc_cutoff = list()
for(a in 1:length(by_animal_3pc_cutoff)){
  animaldf = by_animal_3pc_cutoff[[a]]
  dpi = unique(as.numeric(str_sub(animaldf$sample, start = 7, end = -10)))
  unique_barcodes = unique(animaldf$barcode)
  newdf = data.frame()
  for(d in 1:length(dpi)){
    for(b in 1:length(unique_barcodes)){
      barcode_dpi_comb = animaldf[which(animaldf$barcode == unique_barcodes[b] & animaldf$dpi == dpi[d]), ]
      newdf[b, "animal"] = unique(str_sub(animaldf$sample, end = 5))
      newdf[b, "barcode"] = unique_barcodes[b]
      
      if("rep1" %in% barcode_dpi_comb$rep){
        newdf[b, paste0(dpi[d], "_rep1_proportion")] = barcode_dpi_comb[which(barcode_dpi_comb$rep == "rep1"), "proportion"]
        newdf[b, paste0(dpi[d], "_rep1_count")] = barcode_dpi_comb[which(barcode_dpi_comb$rep == "rep1"), "count"]
      }else{
        newdf[b, paste0(dpi[d], "_rep1_proportion")] = NA
        newdf[b, paste0(dpi[d], "_rep1_count")] = NA
      }
      if("rep2" %in% barcode_dpi_comb$rep){
        newdf[b, paste0(dpi[d], "_rep2_proportion")] = barcode_dpi_comb[which(barcode_dpi_comb$rep == "rep2"), "proportion"]
        newdf[b, paste0(dpi[d], "_rep2_count")] = barcode_dpi_comb[which(barcode_dpi_comb$rep == "rep2"), "count"]
      }else{
        newdf[b, paste0(dpi[d], "_rep2_proportion")] = NA
        newdf[b, paste0(dpi[d], "_rep2_count")] = NA
      }
    }
  }
  reorganized_byday_byanimal_3pc_cutoff[[a]] = newdf
}
names(reorganized_byday_byanimal_3pc_cutoff) = names(by_animal_3pc_cutoff)

by_day_by_animal_proportions_3pc = list()
for(a in 1:length(reorganized_byday_byanimal_3pc_cutoff)){
  by_day_by_animal_proportions_3pc[[a]] = reorganized_byday_byanimal_3pc_cutoff[[a]][, c(1:3, grep("proportion", colnames(reorganized_byday_byanimal_3pc_cutoff[[a]])))]
  write.csv(by_day_by_animal_proportions_3pc[[a]], paste0(names(reorganized_byday_byanimal_3pc_cutoff)[a], "_reorganized_proportions.csv"))
}


barcode_5pc_cutoff = list()
for(s in 1:length(sample_list)){
  sample_list[[s]]$proportion = as.numeric(sample_list[[s]]$proportion)
  data = sample_list[[s]][which(sample_list[[s]]$proportion > 0.05), ]
  otherrow = data.frame(names(sample_list)[s], "other", "other", sum(as.numeric(sample_list[[s]]$count)) - sum(as.numeric(data$count)), 1-sum(data$proportion))
  colnames(otherrow) = colnames(data)
  barcode_3pc_cutoff[[s]] = rbind(data, otherrow)
}
names(barcode_3pc_cutoff) = names(sample_list)
