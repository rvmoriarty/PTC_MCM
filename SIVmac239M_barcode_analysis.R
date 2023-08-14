### Import Excel sheets and format data for Viral Evolution Study - ###
### Which lineages are present over the study, and at what proportions/concentrations?  ###
### Written in R version 4.2.1 ###
### Written by Ryan V Moriarty, August 2023 ### 
options(stringsAsFactors = FALSE)

#### Packages to install #### 
if("vegan" %in% installed.packages() == F){install.packages("vegan")}
if("stringr" %in% installed.packages() == F){install.packages("stringr")}
if("stringdist" %in% installed.packages() == F){install.packages("stringdist")}
if("readxl" %in% installed.packages() == F){install.packages("readxl")}
if("dplyr" %in% installed.packages() == F){install.packages("dplyr")}
if("plyr" %in% installed.packages() == F){install.packages("plyr")}
if("DescTools" %in% installed.packages() == F){install.pacakges("DescTools")}

library(vegan)
library(stringr)
library(stringdist)
library(readxl)
library(plyr)
library(dplyr)
library(DescTools)

#### Define functions used to analyze data #### 
# Read all pages of the excel sheet made by BAF.DistributedScript.R and SIVBarcodeFinder_loop.R #
read_excel_allsheets = function(filename, tibble = FALSE){
  sheets =readxl::excel_sheets(filename)
  x = lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x = lapply(x, as.data.frame)
  names(x) = sheets
  x
}

# Rename barcodes that are not part of the known list to something more informative 
fix_unid_barcode_names = function(unid_seq_table){
  for(row in 1:nrow(unid_seq_table)){
    closest = which(stringdist(unid_seq_table[row, "Sequence"], known_barcodes$Sequence, method = "hamming") == min(stringdist(unid_seq_table[row, "Sequence"], known_barcodes$Sequence, method = "hamming")))
    if(length(closest) != 1){
      unid_seq_table[row, "Closest_match"] = paste0("Multiple close matches: ", str_flatten(known_barcodes[closest, "Name"], ", "), ". Discard.")
      unid_seq_table[row, "Num_mismatches"] = str_flatten(unique(stringdist(unid_seq_table[row, "Sequence"], known_barcodes[closest, "Sequence"])), ", ")
      unid_seq_table[row, "NewName"] = "Multiple close matches - discard."
    }else{
      unid_seq_table[row, "Closest_match"] = known_barcodes[closest, "Name"]
      unid_seq_table[row, "Num_mismatches"] = stringdist(unid_seq_table[row, "Sequence"], known_barcodes[closest, "Sequence"])
      
      mismatches = sort(which(str_split(known_barcodes[closest, "Sequence"], "")[[1]] != str_split(unid_seq_table[row, "Sequence"], "")[[1]]))
      if(length(mismatches[c(TRUE,diff(mismatches) != 1)]) > 1){
        unid_seq_table[row, "NewName"] = paste0("Discarded due to non-consecutive mismatches")
      }else{
        start_pos = mismatches[1]
        stop_pos = mismatches[length(mismatches)]
        
        if(length(mismatches) == 1 & length(closest) == 1){
          unid_seq_table[row, "NewName"] = paste0(known_barcodes[closest, "Name"], ".", start_pos, str_split(known_barcodes[closest, "Sequence"], "")[[1]][start_pos], ">", str_split(unid_seq_table[row, "Sequence"], "")[[1]][start_pos])
        }else{
          unid_seq_table[row, "NewName"] = paste0(known_barcodes[closest, "Name"], ".", start_pos, "_", stop_pos, "delins", str_flatten(str_split(unid_seq_table[row, "Sequence"], "")[[1]][start_pos:stop_pos]))
        }
      }
    }
  }
  return(unid_seq_table)
}

#### Identify working directory and declare sets and animals ####
setwd("~/Documents/BarcodeAnalysisTool/PTC_239M_lineages")
sets = paste0("set", c(seq(1:16))) 

## Identify the study ranges 
ranges = data.frame(preART = c(0, 15), postART = c(268, 442), postRC = c(443, 497), postDep = c(497, 560))

## Define animals and the infection date 
animals = c("cy1035", "cy1036", "cy1039", "cy1043", "cy1040", "cy1044", "cy1045")
blippers = c("cy1035", "cy1040", "cy1043", "cy1045")
infx_date = as.Date("2020-01-27") 

## Definie your known barcodes 
known_barcodes = read.delim("~/Documents/BarcodeAnalysisTool/KnownBarcodes.tsv", header = F)
colnames(known_barcodes) = c("Name", "Sequence")

## Import a csv of your viral loads - 
## columns should be labeled "animal", "sampletype", "date" (as YYYY-MM-DD format), "vl", and "days post SIV"
viralloads = read.csv("~/Documents/BarcodeAnalysisTool/viralloads2.csv")
viralloads$sampletype = gsub("LN|\\s+", "", viralloads$sampletype)

#### Pull in the Analysis Sheets #### 
analysis_outputs = list()
for(s in 1:length(sets)){
  analysis_file = list.files(path = paste0(PATH_TO_ANALYSIS_FILE, sets[s]), pattern = "Analysis.xlsx")
  analysis_outputs[[s]] = read_excel_allsheets(paste0(PATH_TO_ANALYSIS_FILE, sets[s], "/", analysis_file))
}
names(analysis_outputs) = sets

#### Separate by "Samples" or "Matrix", keep "runinfo" sheets #### 
analysis_output_samples = list()
analysis_output_matrix = list()
analysis_runinfo = list()
for(s in 1:length(analysis_outputs)){
  analysis_output_xlsx = analysis_outputs[[s]]
  analysis_output_samples[[s]] = analysis_output_xlsx[grep("samples", names(analysis_output_xlsx))]
  analysis_output_matrix[[s]] = analysis_output_xlsx[grep("matrix", names(analysis_output_xlsx))]
  analysis_runinfo[[s]] = analysis_output_xlsx[[grep("runinfo", names(analysis_output_xlsx))]]
}
names(analysis_output_samples) = names(analysis_output_matrix) = names(analysis_runinfo) = names(analysis_outputs)

analysis_outputs_sep = list(analysis_output_samples, analysis_output_matrix)
names(analysis_outputs_sep) = c("samples", "matrix")

#### Make a full list of the run info sheets that I can export and make notes on ####
## Sometimes you need to manually go through the data and make sure things make sense ##  
comb_runinfo = data.frame()
for(a in 1:length(analysis_runinfo)){
  df = analysis_runinfo[[a]]
  df = df[which(!is.na(df$Animal)), ]
  df$`Run Number` = names(analysis_runinfo)[a]
  df = df[, c("Run Number", "Animal", "Sample", "Date", "Barcodes", "# of 5' Indexing Sequences", "Total Sequences Extracted per Primer", "# of Sequences Matching Known Barcode", "# of Known Barcodes", "Sequencing Input")]
  comb_runinfo = rbind(comb_runinfo, df)
}

#### Filter out samples with less than 5000 reads and less than 10 input templates (initially set to 200 templates, but adjusted for some of the little blips  #### 
analysis_output_samples_filtered = list()
for(a in 1:length(analysis_outputs)){
  setinfo = analysis_outputs[[a]][[grep("runinfo", names(analysis_outputs[[a]]))]]
  sample_dfs = analysis_outputs[[a]][grep("samples", names(analysis_outputs[[a]]))]
  
  ## Adjust the sample info sheet so it's easier to match to the matrix ##  
  ## Also make new "Samples" and "Matrix" Sheets that only have samples that have sufficient sequencing depth 
  ok_dfs = list()
  for(b in 1:length(sample_dfs)){
    num_samples = ncol(sample_dfs[[b]])/8
    sample_ids = data.frame()
    goodcols = c()
    for(c in 1:num_samples){ ## This might need to be adjusted depending on how your columns in the analysis_outputs "Samples" sheets are named 
      sample_ids[c, "an_id"] = gsub("\\...\\d{1,2}", "", colnames(sample_dfs[[b]])[c*8-7])
      sample_ids[c, "sampletype"] = gsub("\\...\\d{1,2}", "", colnames(sample_dfs[[b]])[c*8-6])
      sample_ids[c, "date_id"] = gsub("\\...\\d{1,2}", "", colnames(sample_dfs[[b]])[c*8-5])
      sample_ids[c, "num_seqs"] = sample_dfs[[b]][2, c*8-6]
      sample_ids[c, "input"] = as.numeric(sample_dfs[[b]][2, c*8-7])
      
      setinfo[which(setinfo$Animal == sample_ids[c, "an_id"] & as.character(setinfo$Date) == sample_ids[c, "date_id"] & setinfo$`# of Sequences Matching Known Barcode` == sample_ids[c, "num_seqs"]), "Sample"] = sample_ids[c, "sampletype"]
      
      if(as.numeric(sample_ids[c, "num_seqs"]) >= 5000 & sample_ids[c, "input"] >= 10){
        goodcols = c(goodcols, (c*8-7):(c*8))
      }
    }
    ok_dfs[[b]] = sample_dfs[[b]][, goodcols]
  }
  names(ok_dfs) = names(sample_dfs)
  
  analysis_output_samples_filtered[[a]] = ok_dfs
}
names(analysis_output_samples_filtered) = names(analysis_outputs)
analysis_outputs_sep = list(analysis_output_samples_filtered)

#### Remove a layer of lists so it's easier to loop through #### 
analysis_outputs_sep_flat = list()
for(a in 1:length(analysis_outputs_sep)){
  analysis_outputs_sep_flat[[a]] = unlist(analysis_outputs_sep[[a]], recursive = F)
}
names(analysis_outputs_sep_flat) = c("samples")

#### Group each sample by animal ####
animal_data = list()
for(a in 1:length(animals)){
  animal_data[[a]] = analysis_outputs_sep_flat$samples[grep(animals[a], names(analysis_outputs_sep_flat$samples))]
}
names(animal_data) = animals

#### Identify and combine replicates #### 
analysis_output_replicates = list()
combined_sample_dfs = list()
for(a in 1:length(animal_data)){
  animal_id = names(animal_data)[a]
  animal_samples = animal_data[[a]]
  animal_samples = animal_samples[which(lengths(animal_samples) > 0)]
  
  ## Separate all the samples into their own "list", named by what the matrix would be named ##
  separated_samples = list()
  c = 1 
  samplenames = c()
  for(d in 1:length(animal_samples)){
    num_samples = ncol(animal_samples[[d]])/8
    for(e in 1:num_samples){
      separated_samples[[c]] = animal_samples[[d]][, (e*8-7):(e*8)]
      samplenames[c] = gsub("\\...\\d{1,4}", "", str_flatten(colnames(animal_samples[[d]][, c((e*8-7),(e*8-6), (e*8-5))]), " "))
      c = c + 1 
    }
  }
  names(separated_samples) = samplenames
  
  ## find unique Sampletype/Date sets ## 
  samplenames = gsub("Plasma \\d{1}\\s", "Plasma ", samplenames)
  all_samples = data.frame(str_extract(samplenames, "cy\\d{4}"), unlist(lapply(str_split(samplenames, " "), "[[", 2)), str_extract(samplenames, "\\d{4}-\\d{2}-\\d{2}"))
  colnames(all_samples) = c("animal", "sampletype", "date")
  all_samples = unique(all_samples)
  
  ## Match up replicates and make combined data frames for each sample  ## 
  sampletypes = unique(all_samples$sampletype)
  alldates = sort(unique(all_samples$date))
  combined_dfs = list()
  z = 1 
  all_bcs_comb = c()
  min_input_templates = Inf
  for(d in 1:length(sampletypes)){
    for(e in 1:length(alldates)){
      dfs = separated_samples[which(str_detect(names(separated_samples), sampletypes[d]) & str_detect(names(separated_samples), alldates[e]))]
      if(length(dfs) > 0){
        ## Adjust the data frames so we only have the barcodes that passed filters ##
        ok_bcs = c()
        adjusted_dfs = list()
        for(f in 1:length(dfs)){
          newdf = dfs[[f]][5:(which(is.na(dfs[[f]][,1]))[2]-1), ]
          input_templates = as.numeric(dfs[[f]][2, 1])
          if(input_templates < min_input_templates){
            min_input_templates = input_templates
          }
          colnames(newdf) = dfs[[f]][4, ]
          
          # Fix any unidentified barcodes # 
          newdf_unid = newdf[grep("Unique", newdf$Barcode), c("Sequence", "Barcode")]
          if(nrow(newdf_unid) > 0){
            newdf_unid = fix_unid_barcode_names(newdf_unid)
            for(g in 1:nrow(newdf_unid)){
              newdf[which(newdf$Sequence == newdf_unid[g, "Sequence"]), "Barcode"] = newdf_unid[g, "NewName"]
            }
          }
          
          # Remove barcodes that have multiple matches or non-consecutive mismatches 
          newdf = newdf[grep("SIV", newdf$Barcode), ]
          newdf$Proportion = as.numeric(newdf$Proportion)/sum(as.numeric(newdf$Proportion))
          ok_bcs = unique(c(ok_bcs, newdf$Barcode))
          adjusted_dfs[[f]] = newdf
        }
        names(adjusted_dfs) = names(dfs)
        
        comb_data = data.frame()
        for(f in 1:length(adjusted_dfs)){
          if(any(duplicated(names(adjusted_dfs)))){
            names(adjusted_dfs)[which(duplicated(names(adjusted_dfs)))] = paste(names(adjusted_dfs)[which(duplicated(names(adjusted_dfs)))], "2")
          }
          for(g in 1:length(ok_bcs)){
            check.exists = adjusted_dfs[[f]][which(adjusted_dfs[[f]]$Barcode == ok_bcs[g]), ]
            if(nrow(check.exists) == 1){
              comb_data[g, paste(names(adjusted_dfs[f]), "counts")] = as.numeric(check.exists$Counts)
              comb_data[g, paste(names(adjusted_dfs[f]), "prop")] = as.numeric(check.exists$Proportion)
            }else if(nrow(check.exists) > 1){
              print("wtf")
            }else{
              comb_data[g, paste(names(adjusted_dfs[f]), "counts")] = NA
              comb_data[g, paste(names(adjusted_dfs[f]), "prop")] = NA 
            }
          }
        }
        rownames(comb_data) = ok_bcs
        
        if(length(grep("count", colnames(comb_data))) > 1){
          comb_data$count_sum = rowSums(comb_data[, grep("count", colnames(comb_data))], na.rm = T)
        }else{
          comb_data$count_sum = comb_data[, grep("count", colnames(comb_data))]
        }
        comb_data["comb_sum", ] = colSums(comb_data, na.rm = T)
        comb_data$comb_prop = comb_data$count_sum/comb_data["comb_sum", "count_sum"]
        
        combined_dfs[[z]] = comb_data
        names(combined_dfs)[z] = paste(animal_id, sampletypes[d], alldates[e], "comb")
        z = z + 1 
        
        all_bcs_comb = unique(c(all_bcs_comb, ok_bcs))
      }
    }
  }
  combined_sample_dfs[[a]] = combined_dfs 
  
  ## Actually combine the replicate data ##  
  animal_matrix_repcomb = data.frame()
  for(d in 1:length(combined_dfs)){
    for(e in 1:length(all_bcs_comb)){
      check.exists = combined_dfs[[d]][which(rownames(combined_dfs[[d]]) == all_bcs_comb[e]), ]
      if(nrow(check.exists) == 1){
        animal_matrix_repcomb[all_bcs_comb[e], names(combined_dfs[d])] = check.exists[1, "comb_prop"]
      }else{
        animal_matrix_repcomb[all_bcs_comb[e], names(combined_dfs[d])] = NA
      }
    }
  }
  
  # Adjust the proportions such that any sample post-ART has a filter based on minimum input
  animal_matrix_repcomb_inputfilter = animal_matrix_repcomb
  min_prop = max(1/min_input_templates, 1/200) # just in case there were fewer than 200 input templates, as was the instance for some blips 
  for(d in 1:ncol(animal_matrix_repcomb_inputfilter)){
    dpi = as.numeric(as.Date(str_extract(colnames(animal_matrix_repcomb_inputfilter)[d], "\\d{4}-\\d{2}-\\d{2}")) - infx_date)
    if(dpi > 14){
      animal_matrix_repcomb_inputfilter[which(animal_matrix_repcomb_inputfilter[,d] < min_prop), d] = NA
    }
  }
  analysis_output_replicates[[a]] = list(animal_matrix_repcomb, animal_matrix_repcomb_inputfilter)
  names(analysis_output_replicates[[a]]) = c("repcomb", "inputfilter")
}
names(combined_sample_dfs) = names(animal_data)
names(analysis_output_replicates) = names(animal_data)

#### Find all barcodes in all animals #### 
total_bcs_found = c()
for(a in 1:length(analysis_output_replicates)){
  total_bcs_found = unique(c(total_bcs_found, rownames(analysis_output_replicates[[a]][[1]])))
}
total_bcs_found_nounique = total_bcs_found[grep(">|_", total_bcs_found, invert = T)]

#### Change frequencies to copies in the analysis_output_replicates #### 
analysis_output_replicates_copies = list()
for(a in 1:length(analysis_output_replicates)){
  midlist = list()
  for(z in 1:length(analysis_output_replicates[[a]])){
    df = analysis_output_replicates[[a]][[z]]
    df2 = data.frame(row.names = rownames(df))
    for(b in 1:ncol(df)){
      coldate = str_extract(colnames(df)[b], "\\d{4}-\\d{2}-\\d{2}")
      colsample = str_trim(gsub("cy\\d{4}|\\d{4}-\\d{2}-\\d{2}|comb", "", colnames(df)[b]))
      vl = viralloads[which(viralloads$animal == names(analysis_output_replicates[a]) & viralloads$sampletype == colsample & viralloads$date == coldate), "vl"]
      if(length(vl) == 0){
        stop("Viral load not found!")
      }
      df2[,b] = as.numeric(df[,b])*vl
    }
    colnames(df2) = colnames(df)
    midlist[[z]] = df2
  }
  names(midlist) = c("repcomb", "inputfilter")
  analysis_output_replicates_copies[[a]] = midlist
}
names(analysis_output_replicates_copies) = names(analysis_output_replicates)

# Which barcodes were present during peak viremia?
bcs_at_peak_vl_11dpi = list()
for(a in 1:length(animals)){
  if(animals[a] == "cy1035"){
    date = "2020-02-10"
  }else{
    date = "2020-02-07"
  }
  bcs_at_peak_vl_11dpi[[a]] = rownames(analysis_output_replicates[[animals[a]]][[1]][!is.na(analysis_output_replicates[[animals[a]]][[1]][,grep(date, colnames(analysis_output_replicates[[animals[a]]][[1]]))]), ])
}
names(bcs_at_peak_vl_11dpi) = animals 

bcs_at_peak_postdep = list()
for(a in 1:length(animals)){
  vls_post_dep = viralloads[which(viralloads$days.post.SIV > 497 & viralloads$sampletype == "Plasma"), ]
  postdep_peak = max(vls_post_dep[which(vls_post_dep$animal == animals[a]), "vl"], na.rm = T)
  
  midlist = list()
  for(b in 1:length(analysis_output_replicates[[animals[a]]])){
    midlist[[b]] = rownames(analysis_output_replicates[[animals[a]]][[b]][!is.na(analysis_output_replicates[[animals[a]]][[b]][,grep(vls_post_dep[which(vls_post_dep$vl == postdep_peak & vls_post_dep$animal == animals[a]), "date"], colnames(analysis_output_replicates[[animals[a]]][[b]]))]), ])
  }
  names(midlist) = c("repcomb", "inputfilter")
  bcs_at_peak_postdep[[a]] = midlist
}
names(bcs_at_peak_postdep) = animals

#### Only include barcodes that were present during peak viremia, either acute or post-depletion #### 
aor_bcs_at_peaks = list()
for(a in 1:length(analysis_output_replicates)){
  midlist = list()
  for(c in 1:length(analysis_output_replicates[[a]])){
    df = analysis_output_replicates[[a]][[c]]
    comb_bcs_at_peak = unique(c(bcs_at_peak_vl_11dpi[[names(analysis_output_replicates)[a]]], bcs_at_peak_postdep[[names(analysis_output_replicates)[[a]]]][[c]]))
    newdf = df[which(rownames(df) %in% comb_bcs_at_peak), ]
    midlist[[c]] = newdf[which(rownames(newdf) %in% comb_bcs_at_peak), ]
  }
  names(midlist) = c("repcomb", "inputfilter")
  aor_bcs_at_peaks[[a]] = midlist
}
names(aor_bcs_at_peaks) = names(analysis_output_replicates)

## Convert the proportions to copies/mL of plasma ## 
aor_bcs_at_peaks_copies = list()
for(a in 1:length(analysis_output_replicates_copies)){
  midlist = list()
  for(b in 1:length(analysis_output_replicates_copies[[a]])){
    df = analysis_output_replicates_copies[[a]][[b]]
    comb_bcs_at_peak = unique(c(bcs_at_peak_vl_11dpi[[names(analysis_output_replicates_copies)[a]]], bcs_at_peak_postdep[[names(analysis_output_replicates_copies)[a]]][[b]]))
    midlist[[b]] = df[which(rownames(df) %in% comb_bcs_at_peak), ]
  }
  names(midlist) = c("repcomb", "inputfilter")
  aor_bcs_at_peaks_copies[[a]] = midlist
}
names(aor_bcs_at_peaks_copies) = names(analysis_output_replicates_copies)

## Convert the YYYY-MM-DD to days post infection - proportion of population ##
aor_bcs_at_peaks_dpi = list()
for(a in 1:length(aor_bcs_at_peaks)){
  midlist = list()
  for(b in 1:length(aor_bcs_at_peaks[[a]])){
    df = aor_bcs_at_peaks[[a]][[b]]
    dpi = as.character(as.Date(str_extract(colnames(df)[grep("Plasma", colnames(df))], "\\d{4}-\\d{2}-\\d{2}")) - infx_date)
    tissues = str_trim(gsub("cy\\d{4}|comb|\\d{4}-\\d{2}-\\d{2}", "", colnames(df)[grep("Plasma", colnames(df), invert = T)]))
    colnames(df) = c(dpi, tissues)
    midlist[[b]] = df[, c(as.character(sort(as.numeric(colnames(df)[grep("\\d", colnames(df))]))), tissues)]  
  }
  names(midlist) = c("repcomb", "inputfilter")
  aor_bcs_at_peaks_dpi[[a]] = midlist
}
names(aor_bcs_at_peaks_dpi) = names(aor_bcs_at_peaks)

## Convert the YYYY-MM-DD to days post infection ##
aor_bcs_at_peaks_copies_dpi = list()
for(a in 1:length(aor_bcs_at_peaks_copies)){
  midlist = list()
  for(b in 1:length(aor_bcs_at_peaks_copies[[a]])){
    df = aor_bcs_at_peaks_copies[[a]][[b]]
    
    dpi = as.character(as.Date(str_extract(colnames(df)[grep("Plasma", colnames(df))], "\\d{4}-\\d{2}-\\d{2}")) - infx_date)
    tissues = str_trim(gsub("cy\\d{4}|comb|\\d{4}-\\d{2}-\\d{2}", "", colnames(df)[grep("Plasma", colnames(df), invert = T)]))
    colnames(df) = c(dpi, tissues)
    
    midlist[[b]] = df[, c(as.character(sort(as.numeric(colnames(df)[grep("\\d", colnames(df))]))), tissues)]  
    midlist[[b]]["total", ] = colSums(df, na.rm = T)
    
    for(c in 1:ncol(midlist[[b]])){
      midlist[[b]][which(midlist[[b]][, c] < 1), c] = NA
    }
  }
  names(midlist) = c("repcomb", "inputfilter")
  aor_bcs_at_peaks_copies_dpi[[a]] = midlist
}
names(aor_bcs_at_peaks_copies_dpi) = names(aor_bcs_at_peaks_copies)

# Summarize the number of barcodes over a given proportion # 
prop_summaries = data.frame()
for(a in 1:length(aor_bcs_at_peaks_dpi)){
  df = aor_bcs_at_peaks_dpi[[a]]$inputfilter[, grep("\\d", colnames(aor_bcs_at_peaks_dpi[[a]]$inputfilter))]
  summarydf = data.frame()
  for(b in 1:ncol(df)){
    summarydf[b, "dpi"] = as.numeric(colnames(df)[b])
    for(c in 1:length(summary(df[,b]))){
      summarydf[b, names(summary(df[,b]))[c]] = summary(df[,b])[c]
    }
    summarydf[b, "num_nonNA"] = length(which(!is.na(df[,b])))
    summarydf[b, "num_over0.01"] = length(which(df[,b] >= 0.01))
    summarydf[b, "num_over0.02"] = length(which(df[,b] >= 0.02))
    summarydf[b, "num_over0.03"] = length(which(df[,b] >= 0.03))
    summarydf[b, "num_over0.04"] = length(which(df[,b] >= 0.04))
    summarydf[b, "num_over0.05"] = length(which(df[,b] >= 0.05))
  }
  summarydf$animal = names(aor_bcs_at_peaks_dpi)[a]
  if(names(aor_bcs_at_peaks_dpi)[a] %in% blippers){
    summarydf$cohort = "viremic"
  }else{
    summarydf$cohort = "aviremic"
  }
  prop_summaries = rbind(prop_summaries, summarydf)
}

#### Make list of barcodes that are over given thresholds #### 
## If the barcode is ever above the given threshold ##
aor_bcs_at_peaks_thresholds = list()
thresholds = c(0.01, 0.03, 0.05)
for(a in 1:length(aor_bcs_at_peaks_dpi)){
  midlist = list()
  for(z in 1:length(aor_bcs_at_peaks_dpi[[a]])){
    df = aor_bcs_at_peaks_dpi[[a]][[z]]
    
    tlist = list()
    for(b in thresholds){
      tdf = data.frame()
      keep_rows = c()
      for(c in 1:nrow(df)){
        if(any(df[c,] >= b, na.rm = T)){
          keep_rows = c(keep_rows, c)
        }
      }
      tdf = df[keep_rows, ]
      tdf[paste0("other_", b), ] = 1 - colSums(tdf, na.rm = T)
      tlist[[as.character(b)]] = tdf
    }
    midlist[[z]] = tlist
  }
  names(midlist) = c("repcomb", "inputfilter")
  aor_bcs_at_peaks_thresholds[[a]] = midlist
}
names(aor_bcs_at_peaks_thresholds) = names(aor_bcs_at_peaks_dpi)

## Convert proportions of barcodes over threshold to copies/mL ## 
aor_bcs_at_peaks_copies_thresholds = list()
for(a in 1:length(aor_bcs_at_peaks_thresholds)){
  midlist = list()
  for(z in 1:length(aor_bcs_at_peaks_thresholds[[a]])){
    copy_tlist = list()
    for(b in 1:length(aor_bcs_at_peaks_thresholds[[a]][[z]])){
      df = aor_bcs_at_peaks_thresholds[[a]][[z]][[b]][, grep("\\d", colnames(aor_bcs_at_peaks_thresholds[[a]][[z]][[b]]))]
      all_vl = viralloads[which(viralloads$animal == names(aor_bcs_at_peaks_thresholds[a]) & viralloads$sampletype == "Plasma"), c("vl", "days.post.SIV")]
      
      for(d in 1:nrow(all_vl)){
        df["total", as.character(all_vl[d, "days.post.SIV"])] = all_vl[d, "vl"]
      }
      
      for(c in 1:ncol(df)){
        df[grep("total", rownames(df), invert = T), c] = df[grep("total", rownames(df), invert = T),c]*df["total", c]
        df[which(df[, c] < 1), c] = NA
      }
      df = df[, c(as.character(all_vl$days.post.SIV), colnames(df)[grep("\\d", colnames(df), invert = T)])]
      
      #write.csv(df, paste0(names(aor_bcs_at_peaks_thresholds[a]), "_", names(aor_bcs_at_peaks_thresholds[[a]][z]), "_", names(aor_bcs_at_peaks_thresholds[[a]][[z]][b]), "_bcs_at_peaks.csv"))
      copy_tlist[[b]] = df
    }
    names(copy_tlist) = names(aor_bcs_at_peaks_thresholds[[a]][[z]])
    midlist[[z]] = copy_tlist
  }
  names(midlist) = names(aor_bcs_at_peaks_thresholds[[a]])
  aor_bcs_at_peaks_copies_thresholds[[a]] = midlist
}
names(aor_bcs_at_peaks_copies_thresholds) = names(aor_bcs_at_peaks_thresholds)

## Only when the barcode is at or above the given threshold -- NOT REALLY USED ##
aor_bcs_at_peaks_thresholds_hard = list()
for(a in 1:length(aor_bcs_at_peaks_dpi)){
  midlist = list()
  for(z in 1:length(aor_bcs_at_peaks_dpi[[a]])){
    df = aor_bcs_at_peaks_dpi[[a]][[z]]
    
    tlist = list()
    for(b in thresholds){
      tdf = df
      for(c in 1:nrow(tdf)){
        tdf[c, which(tdf[c, ]< b)] = NA 
      }
      tdf[paste0("other_", b), ] = 1 - colSums(tdf, na.rm = T)
      tlist[[as.character(b)]] = tdf
    }
    midlist[[z]] = tlist
  }
  names(midlist) = c("repcomb", "inputfilter")
  aor_bcs_at_peaks_thresholds_hard[[a]] = midlist
}
names(aor_bcs_at_peaks_thresholds_hard) = names(aor_bcs_at_peaks_dpi)

#### MANUALLY INSPECT DATA AND MAKE SURE THINGS MAKE BIOLOGICAL SENSE #### 
## Some things had to be manually adjusted after looking at the data ##
aor_bcs_at_peaks_copies_mod = list()
for(a in animals){
  df = read.csv(list.files(path = "manually_adjusted_lineages", pattern = a, full.names = T), header = T)
  rownames(df) = df$Days.Post.SIVmac239M
  tdf = t(df)

  aor_bcs_at_peaks_copies_mod[[a]] = tdf[grep("Days", rownames(tdf), invert = T), ]
}
names(aor_bcs_at_peaks_copies_mod) = animals

#### Identify the number of rebounding lineages - total lineages present pre-ART to total lineages post-ART ####
# Use the manually adjusted values from above with the input filter (non-thresholded) values to get all bcs present
rebound_data = list()
for(a in 1:length(aor_bcs_at_peaks_copies_mod)){
  moddf = as.data.frame(aor_bcs_at_peaks_copies_mod[[a]][grep("SIV|total", rownames(aor_bcs_at_peaks_copies_mod[[a]])), ])
  moddf$bc = rownames(moddf)
  minor_df = aor_bcs_at_peaks_copies_dpi[[names(aor_bcs_at_peaks_copies_mod[a])]][["inputfilter"]][, grep("\\d", colnames(aor_bcs_at_peaks_copies_dpi[[names(aor_bcs_at_peaks_copies_mod[a])]][["inputfilter"]]))]
  minor_df$bc = rownames(minor_df)
  
  rebound_data[[a]] = rbind.fill(moddf, minor_df[which(rownames(minor_df) %in% rownames(moddf) == F), ])
  rownames(rebound_data[[a]]) = rebound_data[[a]]$bc
}
names(rebound_data) = names(aor_bcs_at_peaks_copies_mod)

# For cy1043, remove the individual 1246 variants and just keep them pooled as one "Variants" row 
rebound_data[["cy1043"]] = rebound_data[["cy1043"]][grep("1246.\\d{1,2}[ACTG]", rownames(rebound_data[["cy1043"]]), invert = T), ]

num_barcodes_overtime = data.frame()
for(a in 1:length(aor_bcs_at_peaks_dpi)){
  df = aor_bcs_at_peaks_dpi[[a]][["inputfilter"]]
  for(b in 1:ncol(df)){
    num_barcodes_overtime[as.character(colnames(df)[b]), names(aor_bcs_at_peaks_dpi)[a]] = length(which(!is.na(df[, b])))
  }
}

#### Identify the barcodes present (in log10 copies/mL) at different intervals, their cumulative VL, and AUC #### 
total_copy_vls = list()
for(a in 1:length(rebound_data)){
  df = rebound_data[[a]][, grep("\\d", colnames(rebound_data[[a]]))]
  if(names(rebound_data)[a] == "cy1035"){
    peak_preART = "14" # due to sample availability, the peak pre-ART sample for cy1035 is 14dpi 
  }else{
    peak_preART = "11"
  }
  
  ## Identify the barcodes present at each time point 
  ## Day 442 is day of SIVmac239 rechallenge 
  ## Day 497 is day of CD8a depletion 
  acute_bcs = rownames(df[which(rowSums(is.na(df[,which(as.numeric(colnames(df)) <= 14)])) != length(which(as.numeric(colnames(df)) <= 14))), which(as.numeric(colnames(df)) <= 14)])
  peak_preART_bcs = rownames(df[which(!is.na(df[, peak_preART])), ])
  post_dep_bcs = rownames(df[which(rowSums(is.na(df[,which(as.numeric(colnames(df)) >= 497)])) != length(which(as.numeric(colnames(df)) >= 497))), which(as.numeric(colnames(df)) >= 497)])
  
  ## add up the total VL (in log10 copies/mL) for each barcode during a given interval 
  total_vl_acute = log10(rowSums(df[, which(as.numeric(colnames(df)) <= 14)], na.rm = T))
  peak_preART_vl = log10(df[, peak_preART])
  total_vl_postART = log10(rowSums(as.data.frame(df[, which(as.numeric(colnames(df)) > 268 & as.numeric(colnames(df)) < 442)]), na.rm = T))
  total_vl_postRC = log10(rowSums(as.data.frame(df[, which(as.numeric(colnames(df)) > 442 & as.numeric(colnames(df)) < 497)]), na.rm = T))
  total_vl_postdep = log10(rowSums(df[, which(as.numeric(colnames(df)) >= 497)], na.rm = T))
  total_vl_predep = log10(10^total_vl_acute + 10^total_vl_postART + 10^total_vl_postRC)
  
  ## Find the AUC for each barcode during each interval 
  auc_df = data.frame()
  for(b in 1:nrow(df)){
    acute_cols = colnames(df[which(as.numeric(colnames(df)) <= 14)])
    postART_cols = colnames(df)[which(as.numeric(colnames(df)) > 268 & as.numeric(colnames(df)) < 442)]
    postRC_cols = colnames(df)[which(as.numeric(colnames(df)) > 442 & as.numeric(colnames(df)) < 497)]
    postdep_cols = colnames(df)[which(as.numeric(colnames(df)) >= 497)]
    auc_df[rownames(df)[b], "acute_auc"] = log10(AUC(x = as.numeric(acute_cols), y = as.numeric(df[b, acute_cols]), na.rm = T))
    auc_df[rownames(df)[b], "postART_auc"] = log10(AUC(x = as.numeric(postART_cols), y = as.numeric(df[b, postART_cols]), na.rm = T))
    auc_df[rownames(df)[b], "postRC_auc"] = log10(AUC(x = as.numeric(postRC_cols), y = as.numeric(df[b, postRC_cols]), na.rm = T))
    auc_df[rownames(df)[b], "postdep_auc"] = log10(AUC(x = as.numeric(postdep_cols), y = as.numeric(df[b, postdep_cols]), na.rm = T))
    auc_df[rownames(df)[b], "predep_auc"] = log10(AUC(x = as.numeric(c(acute_cols, postART_cols, postRC_cols)), y = as.numeric(df[b, c(acute_cols, postART_cols, postRC_cols)]), na.rm = T))
  }
  
  total_copy_vls[[a]] = as.data.frame(cbind(total_vl_acute, peak_preART_vl, total_vl_postART, total_vl_postRC, total_vl_postdep, total_vl_predep, auc_df, names(rebound_data)[a]))
}
names(total_copy_vls) = names(rebound_data)

#### Bin each barcode by pre-ART peak viral load, then determine if it was present post-depletion or not #### 
preART_bc_ranks = list()
vl_bin_list = list()
rebound_bybin = list()
logvl_bin_list = list()
for(a in 1:length(rebound_data)){
  df = rebound_data[[a]][,grep("\\d", colnames(rebound_data[[a]]))]
  
  # Convert to log10 copies 
  for(c in 1:ncol(df)){
    df[,c] = log10(df[,c])
  }
  
  if(names(rebound_data)[a] == "cy1035"){
    peak = "14"
  }else{
    peak = "11"
  }
  
  ordered = order(df[grep("SIV", rownames(df)), ][, peak], decreasing = T)
  
  preART_bc_ranks[[a]] = data.frame("rank" = as.numeric(seq(1, length(ordered))), "BC" = rownames(df)[ordered], "log10VL" = as.numeric(df[ordered, peak]))
  
  # Turn the log10 copies into 0.5log10 bins and group by this value
  bindf = ZeroIfNA(df[grep("SIV", rownames(df)), ])
  maxvl = round_any(max(bindf[, peak], na.rm = T), 0.5, ceiling)
  vl_bins = c(0, seq(3, maxvl, 0.5))
  by_vl_bin = list()
  for(c in 2:length(vl_bins)-1){
    preART_bcs_in_bin = rownames(bindf)[which(bindf[, peak] <= vl_bins[c+1] & bindf[, peak] > vl_bins[c])]
    by_vl_bin[[paste0("vl_", vl_bins[c], "-", vl_bins[c+1])]] = bindf[preART_bcs_in_bin,]
  }
  
  maxlogvl = ceiling(max(bindf[, peak]))
  logvl_bins = c(0, seq(3, maxlogvl))
  by_logvl_bin = list()
  for(c in 2:length(logvl_bins)-1){
    preART_bcs_in_bin = rownames(bindf)[which(bindf[, peak] <= logvl_bins[c+1] & bindf[, peak] > logvl_bins[c])]
    by_logvl_bin[[paste0("vl_", logvl_bins[c], "-", logvl_bins[c+1])]] = bindf[preART_bcs_in_bin,]
  }
  
  # Calculate how many barcodes are in each pre-ART bin and how many are detected in each post-ART interval 
  rebound_bybin_intervals = data.frame()
  for(c in 1:length(by_vl_bin)){
    rebound_bybin_intervals[names(by_vl_bin)[c], "num_peak_bcs"] = length(rownames(by_vl_bin[[c]]))
    rebound_bybin_intervals[names(by_vl_bin)[c], "num_detected_postART"] = length(which(rowSums(by_vl_bin[[c]][, which(as.numeric(colnames(by_vl_bin[[c]])) >= 268 & as.numeric(colnames(by_vl_bin[[c]])) < 442)]) > 0))
    rebound_bybin_intervals[names(by_vl_bin)[c], "num_detected_postRC"] = length(which(rowSums(by_vl_bin[[c]][, which(as.numeric(colnames(by_vl_bin[[c]])) > 442 & as.numeric(colnames(by_vl_bin[[c]])) < 497)]) > 0))
    rebound_bybin_intervals[names(by_vl_bin)[c], "num_detected_postDep"] = length(which(rowSums(by_vl_bin[[c]][, which(as.numeric(colnames(by_vl_bin[[c]])) >= 497)]) > 0))
    rebound_bybin_intervals[names(by_vl_bin)[c], "num_detected_postAny"] =  length(which(rowSums(by_vl_bin[[c]][, which(as.numeric(colnames(by_vl_bin[[c]])) >= 268 & as.numeric(colnames(by_vl_bin[[c]])) < 442 | as.numeric(colnames(by_vl_bin[[c]])) > 442)]) > 0))
    
    rebound_bybin_intervals[names(by_vl_bin)[c], "prop_detected_postART"] = rebound_bybin_intervals[names(by_vl_bin)[c], "num_detected_postART"]/rebound_bybin_intervals[names(by_vl_bin)[c], "num_peak_bcs"]
    rebound_bybin_intervals[names(by_vl_bin)[c], "prop_detected_postRC"] = rebound_bybin_intervals[names(by_vl_bin)[c], "num_detected_postRC"]/rebound_bybin_intervals[names(by_vl_bin)[c], "num_peak_bcs"]
    rebound_bybin_intervals[names(by_vl_bin)[c], "prop_detected_postDep"] = rebound_bybin_intervals[names(by_vl_bin)[c], "num_detected_postDep"]/rebound_bybin_intervals[names(by_vl_bin)[c], "num_peak_bcs"]
    rebound_bybin_intervals[names(by_vl_bin)[c], "prop_detected_postAny"] = rebound_bybin_intervals[names(by_vl_bin)[c], "num_detected_postAny"]/rebound_bybin_intervals[names(by_vl_bin)[c], "num_peak_bcs"]
  }
  rebound_bybin[[a]] = rebound_bybin_intervals
  
  for(c in 1:length(by_vl_bin)){
    by_vl_bin[[c]] = rbind(by_vl_bin[[c]], df["total", ])
  }
  vl_bin_list[[a]] = by_vl_bin
  
  for(c in 1:length(by_logvl_bin)){
    by_logvl_bin[[c]] = rbind(by_logvl_bin[[c]], df["total", ])
  }
  logvl_bin_list[[a]] = by_logvl_bin
}
names(vl_bin_list) = names(rebound_bybin) = names(logvl_bin_list) = names(rebound_data) 

#### Simulate the number of a that would be expected to rebound based on the pre-ART viral load #### 
## Simulations based on Taina's response: 
# number of barcodes pre-ART:  N
# proportion of each barcode: p = c(p1, p2, ... pN)
# pre-ART vL category: VL.cat = c(c1, c2, ... cN)
# number of rebounding lineages: M 

simdata = list()
#pdf("Simdata.pdf")
rebound_stats = list()
chisq.data = data.frame()
for(a in 1:length(rebound_data)){
  df = ZeroIfNA(rebound_data[[a]][,grep("\\d", colnames(rebound_data[[a]]))])
  if(names(rebound_data)[a] == "cy1035"){
    peak = "14"
  }else{
    peak = "11"
  }
  
  ## Calculate the frequencies of each barcode 
  freqdf = data.frame(row.names = rownames(df))
  for(b in 1:ncol(df)){
    copynum_other_bcs = df["total", b] - round(sum(df[grep("SIV", rownames(df)), b]))
    if(copynum_other_bcs < 0){
      # re-scale so that the total copies sum up to the VL
      df[, b] = (df[,b]/round(sum(df[grep("SIV", rownames(df)), b])))*df["total", b]
      copynum_other_bcs = 0
    }
    
    freqdf[grep("SIV", rownames(freqdf)),b] = df[grep("SIV", rownames(df)), b]/df["total", b]
    freqdf["total", b] = df["total", b]
    freqdf["other", b] = copynum_other_bcs/df["total", b]
  }
  colnames(freqdf) = colnames(df)
  
  ## Identify which lineages belong in which bin 
  vlbin_df = vl_bin_list[[names(rebound_data[a])]]
  binlineages = list()
  for(b in 1:length(vlbin_df)){
    binlineages[[b]] = rownames(vlbin_df[[b]])[grep("SIV", rownames(vlbin_df[[b]]))]
    binlineages[[b]] = binlineages[[b]][grep("SIVmac239M.1921|Variants", binlineages[[b]], invert = T)]
  }
  
  ## Find the total copies present in each interval 
  interval_df = total_copy_vls[[names(rebound_data[a])]][grep("SIV",rownames(total_copy_vls[[names(rebound_data[a])]])), grep("vl", colnames(total_copy_vls[[names(rebound_data[a])]]))]
  for(b in 1:ncol(interval_df)){
    interval_df[,b] = as.numeric(interval_df[,b])
    interval_df[which(interval_df[,b] < 0),b] = 0
  }
  
  ## list which barcodes are present during peak pre-ART viremia and which appear post-depletion 
  peak_preART_bcs = rownames(interval_df)[which(interval_df$peak_preART_vl > 0)]
  preART_freq = freqdf[peak_preART_bcs, peak]
  post_dep_bcs = rownames(interval_df)[which(interval_df$total_vl_postdep > 0)]
  post_dep_bcs = post_dep_bcs[grep("SIVmac239M.1921|Variants", post_dep_bcs, invert = T)] ## Make sure to exclude the rechallenge virus
  
  ## Determine how many barcodes from each bin were present post-depletion 
  num_rebound_lineages = unlist(lapply(binlineages, function(x){length(which(post_dep_bcs %in% x))}))
  
  ## Simulate reactivations based on pre-ART viral load 
  ## Pull observed # clonotypes with probability based on pre-ART frequency 
  num.react.sim = c()
  for(sim in 1:1000){
    sim.rebounder = sample(peak_preART_bcs, sum(num_rebound_lineages), prob = preART_freq, replace = F)
    # how many of sampled barcodes in each pre-ART category?
    num.react.sim = rbind(num.react.sim,unlist(sapply(binlineages,function(x){length(which(x%in%sim.rebounder))})))
  }
  
  ## Get mean and quantile data for number of sampled barcodes 
  ## Add 90% credibility intervals and means of rebounder proportion to barplots 
  meansim = colMeans(num.react.sim)
  lowsim = c()
  highsim = c()
  for(i in 1:ncol(num.react.sim)){
    lowsim[i] = sort(num.react.sim[,i])[50]
    highsim[i] = sort(num.react.sim[,i])[950]
  }
  median.sim = apply(num.react.sim,2,median)
  
  # Chi-squared test of fit: is number of reactivated clonotypes in
  # each pre-ART VL category consistent with simulated "expected" proportion?
  mysig = signif(chisq.test(num_rebound_lineages,median.sim,simulate.p.value = TRUE)$p.value,2)
  
  ## Use the predepletion totals instead of just peak pre-ART 
  predep_bcs = rownames(interval_df)[which(interval_df$total_vl_predep > 0)]
  predep_freq = (10^interval_df[predep_bcs, "total_vl_predep"])/sum(10^interval_df[predep_bcs, "total_vl_predep"])
    
  num.react.sim.predep = c()
  for(sim in 1:1000){
    sim.rebounder = sample(predep_bcs, sum(num_rebound_lineages), prob = predep_freq, replace = F)
    # how many of sampled barcodes in each pre-ART category?
    num.react.sim.predep = rbind(num.react.sim,unlist(sapply(binlineages,function(x){length(which(x%in%sim.rebounder))})))
  }
  
  ## Get mean and quantile data for number of sampled barcodes 
  ## Add 90% credibility intervals and means of rebounder proportion to barplots 
  meansim.predep = colMeans(num.react.sim.predep)
  lowsim.pd = c()
  highsim.pd = c()
  for(i in 1:ncol(num.react.sim.predep)){
    lowsim.pd[i] = sort(num.react.sim.predep[,i])[50]
    highsim.pd[i] = sort(num.react.sim.predep[,i])[950]
  }
  median.sim.pd = apply(num.react.sim.predep,2,median)
  chisq.test(num_rebound_lineages, median.sim, simulate.p.value = T)
  
  rd = cbind("num_lineages" = as.numeric(num.tot), "num_rebound_lineages" = num_rebound_lineages, "prop_rebounding" = num_rebound_lineages/num.tot,
             "mean_sim.num.rebound" = meansim, "median.sim" = median.sim, "low.sim" = lowsim, "high.sim" = highsim, 
             "median.sim.prop" = median.sim/num.tot, "low.sim.prop" = lowsim/num.tot, "high.sim.prop" = highsim/num.tot,
             "mean_sim.num.rebound.pd" = meansim.predep, "median.sim.pd" = median.sim.pd, "low.sim.pd" = lowsim.pd, "high.sim.pd" = highsim.pd, 
             "median.sim.prop.pd" = median.sim.pd/num.tot, "low.sim.prop.pd" = lowsim.pd/num.tot, "high.sim.prop.pd" = highsim.pd/num.tot 
             )
  rd = cbind(rd, names(rebound_data)[a])
  rownames(rd) = names(vlbin_df)
  
  rebound_stats[[a]] = rd
  
  chisq.data[names(rebound_data)[a], "preART_xsq"] = signif(chisq.test(num_rebound_lineages,median.sim,simulate.p.value = TRUE)$statistic,3)
  chisq.data[names(rebound_data)[a], "preART_pval"] = signif(chisq.test(num_rebound_lineages,median.sim,simulate.p.value = TRUE)$p.value,3)
  chisq.data[names(rebound_data)[a], "preDep_xsq"] = signif(chisq.test(num_rebound_lineages,median.sim.pd,simulate.p.value = TRUE)$statistic,3)
  chisq.data[names(rebound_data)[a], "preDep_pval"] = signif(chisq.test(num_rebound_lineages,median.sim.pd,simulate.p.value = TRUE)$p.value,3)
  
  # If you want to change the log intervals 
  # logvlbin_df = logvl_bin_list[[names(rebound_data[a])]]
  # logbinlineages = list()
  # for(b in 1:length(logvlbin_df)){
  #   logbinlineages[[b]] = rownames(logvlbin_df[[b]])[grep("SIV", rownames(logvlbin_df[[b]]))]
  #   logbinlineages[[b]] = logbinlineages[[b]][grep("SIVmac239M.1921|Variants", logbinlineages[[b]], invert = T)]
  # }
  #num_rebound_lineages_log = unlist(lapply(logbinlineages, function(x){length(which(post_dep_bcs %in% x))}))
  
  #num.react.sim.log = c()
  #for(sim in 1:1000){
    #sim.rebounder = sample(peak_preART_bcs, sum(num_rebound_lineages_log), prob = preART_freq, replace = F)
    # how many of sampled barcodes in each pre-ART category?
    #num.react.sim.log = rbind(num.react.sim.log,unlist(sapply(logbinlineages,function(x){length(which(x%in%sim.rebounder))})))
  #}
  
  ## Get mean and quantile data for number of sampled barcodes 
  ## Add 90% credibility intervals and means of rebounder proportion to barplots 
  #meansim.log = colMeans(num.react.sim.log)
  #lowsim.log = c()
  # highsim.log = c()
  # for(i in 1:ncol(num.react.sim.log)){
  #   lowsim.log[i] = sort(num.react.sim.log[,i])[50]
  #   highsim.log[i] = sort(num.react.sim.log[,i])[950]
  # }
  # median.sim.log = apply(num.react.sim.log,2,median)
}
names(rebound_stats) = names(rebound_data)

rebound_stats_df = do.call(rbind, rebound_stats)
write.csv(rebound_stats_df, "Rebound_stats.csv")

#### Do a logistic regression for each animal - is the pre-ART viral load predictive of if the barcode will be detected post-depletion? ####
logistic_regressions = list()
for(a in 1:length(total_copy_vls)){
  df = ZeroIfNA(total_copy_vls[[a]][grep("total", rownames(total_copy_vls[[a]]), invert = T), grep("vl|auc", colnames(total_copy_vls[[a]]))])
  
  # Remove any "-Inf" values and make sure everything is numeric 
  for(b in 1:ncol(df)){
    df[which(df[,b] < 0), b] = 0
  }
  df = as.data.frame(df)
  for(b in 1:ncol(df)){
    df[,b] = as.numeric(df[,b])
  }
  
  # Add a column for each post- interval indicating if that lineage was present in the interval 
  postcols = grep("total_vl_post", colnames(df))
  for(b in postcols){
    df[which(df[,b] != 0 ), paste0(colnames(df)[b], "_binary")] = 1
    df[which(df[,b] == 0 ), paste0(colnames(df)[b], "_binary")] = 0
  }
  
  regs = list()
  # Baseline that RCR = VL: 
  reg = glm(total_vl_postdep_binary~peak_preART_vl, df, family = "binomial")
  regs = list("RCR = VL" = summary(reg))
  # For each interval, does adding this factor in help the model? 
  if(length(unique(df$total_vl_postART_binary)) > 1){
    reg2 = glm(total_vl_postdep_binary~peak_preART_vl+total_vl_postART_binary, df, family = "binomial") # is peak pre-ART + reactivating post-ART predictive?
    regs[["Inc_postART"]] = summary(reg2)
    
    reg2_scale = glm(total_vl_postdep_binary~peak_preART_vl+scale(total_vl_postART_binary), df, family = "binomial")
    regs[["Inc_Scaled_postART"]] = summary(reg2_scale)
    
    reg6 = glm(total_vl_postdep_binary~peak_preART_vl+postART_auc, df, family = "binomial") # does adding the post-ART auc improve the model? (the acutal value, not just a yes/no)
    regs[["Inc_postART_auc"]] = summary(reg6)
  }
  
  if(length(unique(df$total_vl_postRC_binary)) > 1){
    reg3 = glm(total_vl_postdep_binary~peak_preART_vl+total_vl_postRC_binary, df, family = "binomial") # is peak pre-ART + reactivating post-RC predictive?
    regs[["Inc_postRC"]] = summary(reg3)
    
    reg3_scale = glm(total_vl_postdep_binary~peak_preART_vl+scale(total_vl_postRC_binary), df, family = "binomial") 
    regs[["Inc_Scaled_postRC"]] = summary(reg3_scale)  
    
    reg7 = glm(total_vl_postdep_binary~peak_preART_vl+postRC_auc, df, family = "binomial") #does adding the post-RC AUC improve the model? 
    regs[["Inc_postRC_auc"]] = summary(reg7)
  }
  
  # Model as a function of viral replication using AUC values 
  reg4 = glm(total_vl_postdep_binary~acute_auc, df, family = "binomial") # is the viral replication (AUC) pre-ART predictive? 
  regs[["acute_auc"]] = summary(reg4)
  
  reg5 = glm(total_vl_postdep_binary~predep_auc, df, family = "binomial") # is the viral replication (AUC) pre-depletion predictive?
  regs[["predep_auc"]] = summary(reg5)
  
  reg_summaries = data.frame()
  for(b in 1:length(regs)){
    reg_summaries[names(regs)[b], "AIC"] = regs[[b]]$aic
    
    for(c in 1:nrow(regs[[b]]$coefficients)){
      reg_summaries[names(regs)[b], paste0(rownames(regs[[b]]$coefficients)[c], "_Estimate")] = regs[[b]]$coefficients[c,1]
      reg_summaries[names(regs)[b], paste0(rownames(regs[[b]]$coefficients)[c], "_pval")] = regs[[b]]$coefficients[c,4]
      reg_summaries[names(regs)[b], paste0(rownames(regs[[b]]$coefficients)[c], "_comb")] = paste0(signif(regs[[b]]$coefficients[c,1], 3), " (", signif(regs[[b]]$coefficients[c,4], 3), ")")
    }
  }
  reg_summaries$animal = names(total_copy_vls)[a]
  reg_summaries$model = names(regs)
  
  # Compare relative likelihoods of the models (lower AIC score - higher AIC score)
  # If dif < 0.05, evidence that the higher AIC model is significantly worse fitting than the lower AIC model
  aic_comps = combn(names(regs), 2)
  dif_df = data.frame()
  for(b in 1:ncol(aic_comps)){
    dif_df[str_flatten(aic_comps[,b], " vs "), paste0(aic_comps[1,b], "_AIC")] = reg_summaries[aic_comps[1,b], "AIC"]
    dif_df[str_flatten(aic_comps[,b], " vs "), paste0(aic_comps[2,b], "_AIC")] = reg_summaries[aic_comps[2,b], "AIC"]
    
    lower  = min(reg_summaries[aic_comps[1,b], "AIC"], reg_summaries[aic_comps[2,b], "AIC"])
    higher = max(reg_summaries[aic_comps[1,b], "AIC"], reg_summaries[aic_comps[2,b], "AIC"])
    dif_df[str_flatten(aic_comps[,b], " vs "), "dif"] = exp((lower - higher)/2)
  }
  logistic_regressions[[a]] = list(regs, reg_summaries, dif_df)
  names(logistic_regressions[[a]]) = c("models", "model_summary", "dif")
}
names(logistic_regressions) = names(total_copy_vls)

lr_summaries = lapply(logistic_regressions, "[[", 2)
lr_summaries_df = do.call(rbind.fill, lr_summaries)

dif = lapply(logistic_regressions, "[[", 3)
dif_dfs = do.call(rbind, dif)

#### Correlate pre-ART viral loads to the post-dep viral loads #### 
log10_lineage_cor = data.frame()
for(a in 1:length(total_copy_vls)){
  df = ZeroIfNA(total_copy_vls[[a]][grep("total", rownames(total_copy_vls[[a]]), invert = T), grep("vl|auc", colnames(total_copy_vls[[a]]))])
  for(b in 1:ncol(df)){
    df[which(df[,b] < 0), b] = 0
  }
  df = as.data.frame(df)
  for(b in 1:ncol(df)){
    df[,b] = as.numeric(df[,b])
  }
  
  df = df[which(df$peak_preART_vl != 0), ]
  log10_lineage_cor[names(total_copy_vls)[a], "peak_preART_vs_postdep_cor"] = cor.test(df$peak_preART_vl, df$total_vl_postdep)$estimate
  log10_lineage_cor[names(total_copy_vls)[a], "peak_preART_vs_postdep_pval"] = cor.test(df$peak_preART_vl, df$total_vl_postdep)$p.value
  
  log10_lineage_cor[names(total_copy_vls)[a], "acute_vs_postdep_cor"] = cor.test(df$total_vl_acute, df$total_vl_postdep)$estimate
  log10_lineage_cor[names(total_copy_vls)[a], "acute_vs_postdep_pval"] = cor.test(df$total_vl_acute, df$total_vl_postdep)$p.value 
  
  log10_lineage_cor[names(total_copy_vls)[a], "acute_AUC_vs_postdep_AUC_cor"] = cor.test(df$acute_auc, df$postdep_auc)$estimate
  log10_lineage_cor[names(total_copy_vls)[a], "acute_AUC_vs_postdep_AUC__pval"] = cor.test(df$acute_auc, df$postdep_auc)$p.value
  
  log10_lineage_cor[names(total_copy_vls)[a], "predep_AUC_vs_postdep_AUC_cor"] = cor.test(df$predep_auc, df$postdep_auc)$estimate
  log10_lineage_cor[names(total_copy_vls)[a], "predep_AUC_vs_postdep_AUC__pval"] = cor.test(df$predep_auc, df$postdep_auc)$p.value
}

