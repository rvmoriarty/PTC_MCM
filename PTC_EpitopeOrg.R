# Organize the epitope CSVs made by the PTC_Epitope_ID.ipynb files
# Written by Ryan V. Moriarty in R 4.2.3, August 2023

#### Pull in all packages that will be needed #### 
options(stringsAsFactors = F)
if("vegan" %in% installed.packages() == F){install.packages("vegan")}
if("stringr" %in% installed.packages() == F){install.packages("stringr")}
if("stringdist" %in% installed.packages() == F){install.packages("stringdist")}
if("Biostrings" %in% installed.packages() == F){install.packages("Biostrings")}

library(vegan)
library(stringr)
library(stringdist)
library(Biostrings)
library(DescTools)
library(plyr)

#### Define some constants and custom functions #### 
import_csvs = function(folderpath, suffix){
  files = list.files(path = folderpath, pattern = suffix, full.names = T)
  filelist = lapply(files, read.csv)
  names(filelist) = gsub(suffix, "", basename(files))
  
  return(filelist)
}
epitope_mutation_names = function(dataframe, epitope_name){
  ## Take the wild type epitope sequence and find where the other detected epitopes differ
  ## Once differences are found, replace the matching AAs with a period and leave the mismatch 
  repsep = FALSE
  if(any(str_detect(colnames(dataframe), "freq|count"))){
    repsep = TRUE
    suffix = str_extract(colnames(dataframe), "freq|count")
    colnames(dataframe) = gsub("\\s+freq|\\s+count|\\s+\\(WT\\)", "", colnames(dataframe)) 
  }
  
  wt_seq = str_split(epitope_seq_ref[which(epitope_seq_ref$EpitopeName == epitope_name), "EpitopeSeq"], "")[[1]]
  other_seqs = colnames(dataframe)[str_detect(colnames(dataframe), epitope_seq_ref[which(epitope_seq_ref$EpitopeName == epitope_name), "EpitopeSeq"], negate = T)]
  other_seqs = other_seqs[grep("other", other_seqs, invert = T)]
  other_seqs = unique(str_trim(gsub("rep\\d{1}|comb", "", other_seqs)))
  original_seqs = other_seqs
  other_seqs = str_split(other_seqs, "")
  if(length(other_seqs) > 0){
    modified_seqs = data.frame()
    for(seq in 1:length(other_seqs)){
      if(length(wt_seq) != length(other_seqs[[seq]])){
        print(paste(length(wt_seq), length(other_seqs[[seq]]), str_flatten(wt_seq), str_flatten(other_seqs[[seq]])))
        stop()
      }
      matches = which(wt_seq == other_seqs[[seq]])
      for(match in matches){
        other_seqs[[seq]][match] = "."
      }
      modified_seqs[seq, "original"] = original_seqs[seq]
      modified_seqs[seq, "mod"] = str_flatten(other_seqs[[seq]])
    }
    
    if(length(grep("\\*", modified_seqs$original)) > 0){
      modified_seqs$original = gsub("\\*", "\\\\*", modified_seqs$original)
      print(paste0("Stop codon found - ", epitope_name, ". Check."))
    }
    for(row in 1:nrow(modified_seqs)){
      if(repsep){
        colnames(dataframe) = str_replace(colnames(dataframe), modified_seqs[row, "original"], modified_seqs[row, "mod"])
      }else{
        colnames(dataframe) = str_replace(colnames(dataframe), modified_seqs[row, "original"], modified_seqs[row, "mod"])
      }
      
    }
    
    if(repsep){
      colnames(dataframe) = paste(colnames(dataframe), suffix)
    }
    
    colnames(dataframe) = gsub(str_flatten(wt_seq), paste(str_flatten(wt_seq), "(WT)"), colnames(dataframe))
  }
  
  return(dataframe)
}

epitope_seq_ref = read.csv("pre_and_post_epitope_seqs.csv")
viralloads = read.csv("viralloads2.csv")
viralloads$sampletype = gsub("LN|\\s+", "", viralloads$sampletype)

#lateART_vls = read.csv("LateART_vls.csv")

animal_pattern = "cy\\d{4}"
lateART = c("cy0719", "cy0721", "cy0723", "cy0724", "cy0726")
earlyVir = c("cy1035", "cy1040", "cy1043", "cy1045")
earlyAv = c("cy1036", "cy1039", "cy1044")
reps = c("rep1", "rep2")
animals = c(lateART, earlyVir, earlyAv)

#### Grab the .csv (merged + filtered and unfiltered) and .txt files you made ####
filtered_csvs = import_csvs(PATH, "merged.filtered.csv")
unfiltered_csvs  = import_csvs(PATH, "_WGS.csv")

txts = list.files(path = PATH, pattern = ".txt", full.names = T)
txt_list = lapply(txts, read.delim)
names(txt_list) = gsub(".barcode_counts.txt", "", basename(txts))

#### For each filtered csv (which is grouped by animal), separate into individual epitopes with dpi in rows and epitope sequences in columns ####
reformatted_csvs = list()
silent_mutations = list()
z = 1
for(a in 1:length(filtered_csvs_list)){
  animal = str_extract(names(filtered_csvs_list[a]), animal_pattern)
  df = filtered_csvs_list[[a]]
  
  epitopes_found = unique(df$epitope)
  dpi_found = unique(df$dpi)
  
  epitope_list = list()
  for(b in epitopes_found){
    edf = df[which(df$epitope == b), ]
    vars = unique(edf$barcode_sequence)
    
    tdf = data.frame()
    for(c in dpi_found){
      for(d in vars){
        for(e in reps){
          check.exists = edf[which(edf$dpi == c & edf$barcode_sequence == d & edf$rep == e), ]
          if(nrow(check.exists) == 1){
            tdf[c, paste(d, e, "count")] = check.exists$barcode_count
            tdf[c, paste(d, e, "freq")] = check.exists$freq
          }else if(nrow(check.exists) == 0){
            tdf[c, paste(d, e, "count")] = NA
            tdf[c, paste(d, e, "freq")] = NA
          }else{
            print(paste("Silent mutations found for:", animal, c, d , e))
            silent_mutations[[z]] = check.exists
            z = z + 1
            tdf[c, paste(d, e, "count")] = sum(check.exists$barcode_count)
            tdf[c, paste(d, e, "freq")] = sum(check.exists$freq)
            #stop("Check check.exists!")
          }
        }
      }
    }
    
    tdf_repcomb = data.frame(row.names = dpi_found)
    tdf_counts =  tdf[, grep("count", colnames(tdf))]
    for(d in vars){
      tdf_repcomb[, paste(d, "count")] = rowSums(tdf_counts[, grep(d, colnames(tdf_counts), fixed = T)], na.rm = T)
    }
    tdf_repcomb_freq = tdf_repcomb/rowSums(tdf_repcomb, na.rm = T)
    
    epitope_list[[b]] = list("reorg" = epitope_mutation_names(tdf, b), "repcomb_counts" = epitope_mutation_names(tdf_repcomb,b), "repcomb_freq" = epitope_mutation_names(tdf_repcomb_freq,b))
  }
  reformatted_csvs[[animal]] = epitope_list
}

#### Regroup by epitope then animal ####  
epitope_freqs = list()
epitopes = epitope_seq_ref$EpitopeName
for(a in epitopes){
  epitope_freqs[[a]] = lapply(reformatted_csvs, "[[", a)
}

#### We are focusing on Gag GW9, Nef RM9, and Rev SP10 ####
epitope_freq_subset = epitope_freqs[grep(str_flatten(c("Gag_GW9", "Nef_RM9", "Rev_SP10"), "|"), names(epitope_freqs))]

#### This is where the .txt files can be reformatted - this should just be redoing what we did in jupyter notebool #### 
reformatted_txts = list()
for(a in 1:length(txt_list)){
  df = txt_list[[a]]
  if(sum(df$barcode_count, na.rm = T) >= 1000){
    df$freq = df$barcode_count/sum(df$barcode_count)
    df$animal = str_extract(df$sample_name, animal_pattern)
    df$dpi = lapply(str_split(names(txt_list)[a], "-"), "[[", 2)[[1]]
    df$rep = str_extract(df$sample_name, "rep\\d{1}|comb")
    df$epitope = lapply(str_split(names(txt_list)[a], "-"), "[[", 4)[[1]]
    
    reformatted_txts[[a]] = df
  }else{
    reformatted_txts[[a]] = NULL
  }
}
names(reformatted_txts) = names(txt_list)
reformatted_txts = reformatted_txts[which(lengths(reformatted_txts) > 0)]
reformatted_txts_df = do.call(rbind, reformatted_txts)

low_txts = list()
for(a in 1:length(txt_list)){
  df = txt_list[[a]]
  if(sum(df$barcode_count, na.rm = T) < 1000 & nrow(df) > 0){
    df$freq = df$barcode_count/sum(df$barcode_count)
    df$animal = str_extract(df$sample_name, animal_pattern)
    df$dpi = lapply(str_split(names(txt_list)[a], "-"), "[[", 2)[[1]]
    df$rep = str_extract(df$sample_name, "rep\\d{1}|comb")
    df$epitope = lapply(str_split(names(txt_list)[a], "-"), "[[", 4)[[1]]
    
    low_txts[[a]] = df
  }else{
    low_txts[[a]] = NULL
  }
}
names(low_txts) = names(txt_list)
low_txts = low_txts[which(lengths(reformatted_txts) > 0)]
low_txts_df = do.call(rbind, low_txts)
low_txts_df = low_txts_df[which(low_txts_df$epitope != "ARF1_QL11"), ]

all_nt_seqs = unique(reformatted_txts_df$barcode_sequence)

reformatted_txts_df_withlow = rbind(reformatted_txts_df, low_txts_df)

#### For each animal and epitope: combine replicates, calculate frequencies, apply minimum frequency threshold, and identify "other" proportion #### 
epitope_repcomb = list()
check_fs = data.frame()
for(a in animals){
  adf = reformatted_txts_df[which(reformatted_txts_df$animal==a), ]
  
  epitopes = unique(adf$epitope)
  epitope_list = list()
  epitope_seq_adjust = data.frame()
  for(b in epitopes){
    all_epitopes = data.frame()
    edf = adf[which(adf$epitope == b), ]
    dates = unique(edf$dpi)
    allnt_seqs = unique(edf$barcode_sequence)
    dlist = list()
    for(c in dates){
      ddf = edf[which(edf$dpi == c), ]
      dnt_seqs = unique(ddf$barcode_sequence)
      
      comb_df = data.frame()
      for(d in dnt_seqs){
        for(e in reps){
          nt_df = ddf[which(ddf$barcode_sequence == d & ddf$rep == e), ]
          comb_df[d, paste0("nt_count_", e)] = sum(nt_df[, "barcode_count"])
          comb_df[d, paste0("nt_freq_", e)] = sum(nt_df[, "freq"])
        }
      }
      rm(nt_df)
      rm(ddf)
      
      comb_df$nt_count_sum = rowSums(comb_df[, grep("count", colnames(comb_df))], na.rm = T)
      comb_df$nt_freq_avg = comb_df$nt_count_sum/sum(comb_df$nt_count_sum)
      
      # Remove anything below a minimum frequency threshold - 
      # our minimum sequencing threshold was 1000 copies/mL 
      # (and we know not all of that gets reverse transcribed, etc), 
      # but we're going to start with a 1/1000 threshold first
      other_freqs = t(as.data.frame(colSums(comb_df[which(comb_df$nt_freq_avg < 1/1000), ], na.rm = T)))
      rownames(other_freqs) = "TGATGATGATGATGATGATGATGATGA" # we want it to translate, but still be able to identify that it's the "other"
      
      comb_df = comb_df[which(comb_df$nt_freq_avg >= 1/1000), ]
      
      # recalculate the frequencies 
      comb_df = rbind(comb_df, other_freqs)
      
      comb_df$dpi = c
      comb_df$aa_seq = unlist(lapply(lapply(lapply(rownames(comb_df), DNAString), translate), as.character))
      comb_df$nt_seq = gsub("\\d{1,3}dpi.","", rownames(comb_df))
      
      all_epitopes = rbind(all_epitopes, comb_df[, c("aa_seq", "nt_seq")])
       
      dlist[[c]] = comb_df
      
      #rm(comb_df)
    }
    
    epitope_list[[b]] = do.call(rbind, dlist)
    rm(edf)
    all_epitopes = all_epitopes[which(all_epitopes$aa_seq != "*********"),]
    
    epitope_rename = t(as.data.frame(unique(all_epitopes[, 1])))
    colnames(epitope_rename) = epitope_rename[1, ]
    
    seq_adjust = as.data.frame(t(epitope_mutation_names(epitope_rename, b)))
    seq_adjust$renamed = rownames(seq_adjust)
    colnames(seq_adjust)[1] = "fullseq"
    s = 1
    seq_adjust2 = data.frame()
    for(seq in 1:length(unique(seq_adjust$fullseq))){
      nt_seqs_found = unique(all_epitopes[which(all_epitopes$aa_seq == unique(seq_adjust$fullseq)[seq]), ])
      nt_seqs_found$renamed = seq_adjust[which(seq_adjust$fullseq == unique(seq_adjust$fullseq)[seq]), "renamed"]
      seq_adjust2 = rbind(seq_adjust2, nt_seqs_found)
    }
    
    seq_adjust2$epitope = b
    epitope_seq_adjust = unique(rbind(epitope_seq_adjust, seq_adjust2))
  }
  rm(adf)
  
  check_fs = rbind(check_fs, epitope_seq_adjust)
  # Reformat the epitope lists so that each epitope is a column and date is a row
  # epitope_list_reformat = list()
  # for(b in names(epitope_list)){
  #   alldates = unique(epitope_list[[b]]$dpi)
  #   allaas = unique(epitope_list[[b]]$aa_seq)
  #   
  #   rdf = data.frame()
  #   for(c in allaas){
  #     for(d in alldates){
  #       check.exists = epitope_list[[b]][which(epitope_list[[b]]$aa_seq == c & epitope_list[[b]]$dpi == d), "nt_freq_avg"]
  #       if(length(check.exists) == 1){
  #         rdf[c, d] = check.exists
  #       }else if(length(check.exists) == 0){
  #         rdf[c,d] = NA
  #       }else{
  #         rdf[c,d] = sum(check.exists)
  #         print(paste("silent mutations for:", c, d))
  #       }
  #     }
  #   }
  #   
  #   rdf$barcode_sequence = rownames(rdf)
  #   rdf$EpitopeSeq = find_epitope_wt(rdf)
  #   rdf2 = t(rdf[, grep("\\d", colnames(rdf))])
  #   colnames(rdf2)[which(colnames(rdf2) == "*********")] = "other"
  #   rdf3 = epitope_mutation_names(rdf2, b)
  #   
  #   # Reorganize so the "other" is first and wild type is second 
  #   wt_index = grep("\\.|other", colnames(rdf3), invert = T)
  #   other_index = which(colnames(rdf3) == "other")
  #   variants = grep("\\.", colnames(rdf3))
  #   
  #   epitope_list_reformat[[b]] = as.data.frame(rdf3[, c(other_index, wt_index, variants)])
  #   rm(rdf, rdf2, rdf3)
  # }

  epitope_repcomb[[a]] = epitope_list
}
check_fs2 = unique(check_fs) 

epitope_repcomb_withlow = list()
check_fs_withlow = data.frame()
for(a in animals){
  adf = reformatted_txts_df_withlow[which(reformatted_txts_df_withlow$animal==a), ]
  
  epitopes = unique(adf$epitope)
  epitope_list = list()
  epitope_seq_adjust = data.frame()
  for(b in epitopes){
    all_epitopes = data.frame()
    edf = adf[which(adf$epitope == b), ]
    dates = unique(edf$dpi)
    allnt_seqs = unique(edf$barcode_sequence)
    dlist = list()
    for(c in dates){
      ddf = edf[which(edf$dpi == c), ]
      dnt_seqs = unique(ddf$barcode_sequence)
      
      comb_df = data.frame()
      for(d in dnt_seqs){
        for(e in reps){
          nt_df = ddf[which(ddf$barcode_sequence == d & ddf$rep == e), ]
          comb_df[d, paste0("nt_count_", e)] = sum(nt_df[, "barcode_count"])
          comb_df[d, paste0("nt_freq_", e)] = sum(nt_df[, "freq"])
        }
      }
      rm(nt_df)
      rm(ddf)
      
      comb_df$nt_count_sum = rowSums(comb_df[, grep("count", colnames(comb_df))], na.rm = T)
      comb_df$nt_freq_avg = comb_df$nt_count_sum/sum(comb_df$nt_count_sum)
      
      # Remove anything below a minimum frequency threshold - 
      # our minimum sequencing threshold was 1000 copies/mL 
      # (and we know not all of that gets reverse transcribed, etc), 
      # but we're going to start with a 1/1000 threshold first
      other_freqs = t(as.data.frame(colSums(comb_df[which(comb_df$nt_freq_avg < 1/1000), ], na.rm = T)))
      rownames(other_freqs) = "TGATGATGATGATGATGATGATGATGA" # we want it to translate, but still be able to identify that it's the "other"
      
      comb_df = comb_df[which(comb_df$nt_freq_avg >= 1/1000), ]
      
      # recalculate the frequencies 
      comb_df = rbind(comb_df, other_freqs)
      
      comb_df$dpi = c
      comb_df$aa_seq = unlist(lapply(lapply(lapply(rownames(comb_df), DNAString), translate), as.character))
      comb_df$nt_seq = gsub("\\d{1,3}dpi.","", rownames(comb_df))
      
      all_epitopes = rbind(all_epitopes, comb_df[, c("aa_seq", "nt_seq")])
      
      dlist[[c]] = comb_df
      
      #rm(comb_df)
    }
    
    epitope_list[[b]] = do.call(rbind, dlist)
    rm(edf)
    all_epitopes = all_epitopes[which(all_epitopes$aa_seq != "*********"),]
    
    epitope_rename = t(as.data.frame(unique(all_epitopes[, 1])))
    colnames(epitope_rename) = epitope_rename[1, ]
    
    seq_adjust = as.data.frame(t(epitope_mutation_names(epitope_rename, b)))
    seq_adjust$renamed = rownames(seq_adjust)
    colnames(seq_adjust)[1] = "fullseq"
    s = 1
    seq_adjust2 = data.frame()
    for(seq in 1:length(unique(seq_adjust$fullseq))){
      nt_seqs_found = unique(all_epitopes[which(all_epitopes$aa_seq == unique(seq_adjust$fullseq)[seq]), ])
      nt_seqs_found$renamed = seq_adjust[which(seq_adjust$fullseq == unique(seq_adjust$fullseq)[seq]), "renamed"]
      seq_adjust2 = rbind(seq_adjust2, nt_seqs_found)
    }
    
    seq_adjust2$epitope = b
    epitope_seq_adjust = unique(rbind(epitope_seq_adjust, seq_adjust2))
  }
  rm(adf)
  
  check_fs_withlow = rbind(check_fs_withlow, epitope_seq_adjust)
  
  
  epitope_repcomb_withlow[[a]] = epitope_list
}
check_fs2_withlow = unique(check_fs_withlow) 

#### For epitopes with weird frameshifty looking stuff, let's figure out where those are and then go and manually fix them 
for(a in epitopes){
  check_fs2[which(check_fs2$epitope == a), "preEpitopeNt"] = gsub("U","T", epitope_seq_ref[which(epitope_seq_ref$EpitopeName == a), "preEpitopeSeq"])
  check_fs2[which(check_fs2$epitope == a), "EpitopeNt"] = gsub("U","T", epitope_seq_ref[which(epitope_seq_ref$EpitopeName == a), "epitopeNtSeq"])
  check_fs2[which(check_fs2$epitope == a), "postEpitopeNt"] = gsub("U","T", epitope_seq_ref[which(epitope_seq_ref$EpitopeName == a), "PostEpitopeSeq"])
  check_fs2[which(check_fs2$epitope ==a), "EpitopeAA"] = epitope_seq_ref[which(epitope_seq_ref$EpitopeName == a), "EpitopeSeq"]
}

epitope_fs_list = list()
for(a in unique(check_fs2$epitope)){
  ref_aa = unique(check_fs2[which(check_fs2$epitope == a), "EpitopeAA"])
  epitope_fs_list[[a]] = check_fs2[which(check_fs2$epitope == a), ]
  epitope_fs_list[[a]] = epitope_fs_list[[a]][which(epitope_fs_list[[a]]$aa_seq != ref_aa), ]
  
  check_seqs = epitope_fs_list[[a]]$renamed
  mismatches = lengths(str_match_all(check_seqs, "\\."))
  adjust_ntseq = check_seqs[which(mismatches < nchar(ref_aa)/2)]
  
  epitope_fs_list[[a]] = epitope_fs_list[[a]][which(epitope_fs_list[[a]]$renamed %in% adjust_ntseq), ]
}
write.csv(epitope_fs_list$Gag_GW9, "GagGW9_frameshift_vars.csv")

epitope_fs_list$Nef_RM9$nt_seq = (aligned(pairwiseAlignment(pattern = epitope_fs_list$Nef_RM9$nt_seq, subject = paste0("GTG",	"AGGCCAAAAGTTCCCCTAAGAACAATG",	"AGT"), type = "overlap")))
write.csv(epitope_fs_list$Nef_RM9, "NefRM9_frameshift_vars.csv")

epitope_fs_list$Rev_SP10$nt_seq = (aligned(pairwiseAlignment(pattern = epitope_fs_list$Rev_SP10$nt_seq, subject = paste0("TAT",	"TCATTTCCTGATCCGCCAACTGATACGCCT",	"CTT"), type = "overlap")))
write.csv(epitope_fs_list$Rev_SP10, "RevSP10_frameshift_vars.csv")

gw9_rm9_sp10_fs_adjust = read.csv("GW9_RM9_SP10_fs_adjusted.csv")
gw9_rm9_sp10_fs_adjust$adjusted_aa = unlist(lapply(lapply(lapply(gw9_rm9_sp10_fs_adjust$adjusted_nt, DNAString), translate), as.character))

#### Go through each epitope_repcomb amd adjust the GW9, RM9, and SP10 frameshifts ####
# Reformat the epitope lists so that each epitope is a column and date is a row
# epitope_list_reformat = list()
# for(b in names(epitope_list)){
#   alldates = unique(epitope_list[[b]]$dpi)
#   allaas = unique(epitope_list[[b]]$aa_seq)
#   
#   rdf = data.frame()
#   for(c in allaas){
#     for(d in alldates){
#       check.exists = epitope_list[[b]][which(epitope_list[[b]]$aa_seq == c & epitope_list[[b]]$dpi == d), "nt_freq_avg"]
#       if(length(check.exists) == 1){
#         rdf[c, d] = check.exists
#       }else if(length(check.exists) == 0){
#         rdf[c,d] = NA
#       }else{
#         rdf[c,d] = sum(check.exists)
#         print(paste("silent mutations for:", c, d))
#       }
#     }
#   }
#   
#   rdf$barcode_sequence = rownames(rdf)
#   rdf$EpitopeSeq = find_epitope_wt(rdf)
#   rdf2 = t(rdf[, grep("\\d", colnames(rdf))])
#   colnames(rdf2)[which(colnames(rdf2) == "*********")] = "other"
#   rdf3 = epitope_mutation_names(rdf2, b)
#   
#   # Reorganize so the "other" is first and wild type is second 
#   wt_index = grep("\\.|other", colnames(rdf3), invert = T)
#   other_index = which(colnames(rdf3) == "other")
#   variants = grep("\\.", colnames(rdf3))
#   
#   epitope_list_reformat[[b]] = as.data.frame(rdf3[, c(other_index, wt_index, variants)])
#   rm(rdf, rdf2, rdf3)
# }

#### Reformat epitope_repcomb #### 

epitope_repcomb_reformat = epitope_repcomb
for(a in 1:length(epitope_repcomb)){
  epitope_list_reformat = list()
  for(b in 1:length(epitope_repcomb[[a]])){
    if(names(epitope_repcomb[[a]])[b] %in% c("Gag_GW9", "Nef_RM9", "Rev_SP10")){
      adjustdf = epitope_repcomb[[a]][[b]]
      check_adjusted = gw9_rm9_sp10_fs_adjust[which(gw9_rm9_sp10_fs_adjust$epitope == names(epitope_repcomb[[a]])[b]), ]
      
      rownames(adjustdf) = NULL
      check_aas = unique(adjustdf[, c("nt_seq", "aa_seq")])
      for(c in 1:nrow(check_aas)){
        if(check_aas[c, "aa_seq"] %in% check_adjusted$fullseq){
          if(length(unique(check_adjusted[which(check_adjusted$fullseq == check_aas[c, "aa_seq"]), "adjusted_aa"])) > 1){
           
            ntseq_match = which(str_detect(check_adjusted[which(check_adjusted$fullseq == check_aas[c, "aa_seq"]),"aligned_nt"],  check_aas[c, "nt_seq"]) == T)
            if(length(ntseq_match) !=0){
              adjustdf[which(adjustdf$aa_seq == check_aas[c, "aa_seq"]), "aa_seq"] = unique(check_adjusted[which(check_adjusted$fullseq == check_aas[c, "aa_seq"]), "adjusted_aa"][ntseq_match])
            }else{
              print("Nt substring not found in aligned nt")
              
            }
            
          }else{
            adjustdf[which(adjustdf$aa_seq == check_aas[c, "aa_seq"]), "aa_seq"] = unique(check_adjusted[which(check_adjusted$fullseq == check_aas[c, "aa_seq"]), "adjusted_aa"])
          }
          
        }
      }
    }else{
      adjustdf =  epitope_repcomb[[a]][[b]]
    }
    
    alldates = unique(adjustdf$dpi)
    allaas = unique(adjustdf$aa_seq)
    
    rdf = data.frame()
    for(c in allaas){
      for(d in alldates){
        check.exists = adjustdf[which(adjustdf$aa_seq == c & adjustdf$dpi == d), "nt_freq_avg"]
        if(length(check.exists) == 1){
          rdf[c, d] = check.exists
        }else if(length(check.exists) == 0){
          rdf[c,d] = NA
        }else{
          rdf[c,d] = sum(check.exists)
          print(paste("silent mutations for:", c, d))
        }
      }
    }
    
    rdf$barcode_sequence = rownames(rdf)
    rdf$EpitopeSeq = find_epitope_wt(rdf)
    rdf2 = t(rdf[, grep("\\d", colnames(rdf))])
    colnames(rdf2)[which(colnames(rdf2) == "*********")] = "other"
    rdf3 = epitope_mutation_names(rdf2, names(epitope_repcomb[[a]])[b])
    
    
    # Reorganize so the "other" is first and wild type is second
    wt_index = grep("\\.|other", colnames(rdf3), invert = T)
    other_index = which(colnames(rdf3) == "other")
    variants = grep("\\.", colnames(rdf3))
    
    epitope_list_reformat[[b]] = as.data.frame(rdf3[, c(other_index, wt_index, variants)])
    rm(rdf, rdf2, rdf3)
  }
  names(epitope_list_reformat) = names(epitope_repcomb[[a]])
  epitope_repcomb_reformat[[a]] = epitope_list_reformat
}


epitope_repcomb_reformat_wl = epitope_repcomb_withlow
for(a in 1:length(epitope_repcomb_withlow)){
  epitope_list_reformat = list()
  for(b in 1:length(epitope_repcomb_withlow[[a]])){
    if(names(epitope_repcomb_withlow[[a]])[b] %in% c("Gag_GW9", "Nef_RM9", "Rev_SP10")){
      adjustdf = epitope_repcomb_withlow[[a]][[b]]
      check_adjusted = gw9_rm9_sp10_fs_adjust[which(gw9_rm9_sp10_fs_adjust$epitope == names(epitope_repcomb_withlow[[a]])[b]), ]
      
      rownames(adjustdf) = NULL
      check_aas = unique(adjustdf[, c("nt_seq", "aa_seq")])
      for(c in 1:nrow(check_aas)){
        if(check_aas[c, "aa_seq"] %in% check_adjusted$fullseq){
          if(length(unique(check_adjusted[which(check_adjusted$fullseq == check_aas[c, "aa_seq"]), "adjusted_aa"])) > 1){
            
            ntseq_match = which(str_detect(check_adjusted[which(check_adjusted$fullseq == check_aas[c, "aa_seq"]),"aligned_nt"],  check_aas[c, "nt_seq"]) == T)
            if(length(ntseq_match) !=0){
              adjustdf[which(adjustdf$aa_seq == check_aas[c, "aa_seq"]), "aa_seq"] = unique(check_adjusted[which(check_adjusted$fullseq == check_aas[c, "aa_seq"]), "adjusted_aa"][ntseq_match])
            }else{
              print("Nt substring not found in aligned nt")
              
            }
            
          }else{
            adjustdf[which(adjustdf$aa_seq == check_aas[c, "aa_seq"]), "aa_seq"] = unique(check_adjusted[which(check_adjusted$fullseq == check_aas[c, "aa_seq"]), "adjusted_aa"])
          }
          
        }
      }
    }else{
      adjustdf =  epitope_repcomb_withlow[[a]][[b]]
    }
    
    alldates = unique(adjustdf$dpi)
    allaas = unique(adjustdf$aa_seq)
    
    rdf = data.frame()
    for(c in allaas){
      for(d in alldates){
        check.exists = adjustdf[which(adjustdf$aa_seq == c & adjustdf$dpi == d), "nt_freq_avg"]
        if(length(check.exists) == 1){
          rdf[c, d] = check.exists
        }else if(length(check.exists) == 0){
          rdf[c,d] = NA
        }else{
          rdf[c,d] = sum(check.exists)
          print(paste("silent mutations for:", c, d))
        }
      }
    }
    
    rdf$barcode_sequence = rownames(rdf)
    rdf$EpitopeSeq = find_epitope_wt(rdf)
    rdf2 = t(rdf[, grep("\\d", colnames(rdf))])
    colnames(rdf2)[which(colnames(rdf2) == "*********")] = "other"
    rdf3 = epitope_mutation_names(rdf2, names(epitope_repcomb_withlow[[a]])[b])
    
    
    # Reorganize so the "other" is first and wild type is second
    wt_index = grep("\\.|other", colnames(rdf3), invert = T)
    other_index = which(colnames(rdf3) == "other")
    variants = grep("\\.", colnames(rdf3))
    
    epitope_list_reformat[[b]] = as.data.frame(rdf3[, c(other_index, wt_index, variants)])
    rm(rdf, rdf2, rdf3)
  }
  names(epitope_list_reformat) = names(epitope_repcomb_withlow[[a]])
  epitope_repcomb_reformat_wl[[a]] = epitope_list_reformat
}

# Reorganize to be ordered by epitope THEN animal 
by_epitope = list()
for(a in epitopes){
  for(b in 1:length(epitope_repcomb_reformat)){
    if(!is.null(epitope_repcomb_reformat[[b]][[a]])){
      if(nrow(epitope_repcomb_reformat[[b]][[a]])> 0){
        by_epitope[[a]][[b]] = epitope_repcomb_reformat[[b]][[a]]
        by_epitope[[a]][[b]]$dpi = rownames(by_epitope[[a]][[b]])
        by_epitope[[a]][[b]]$animal = names(epitope_repcomb_reformat)[b]
        names(by_epitope[[a]])[b] = names(epitope_repcomb_reformat)[b]
      }
    }
  }
}

# Make a large data frame for each epitope that has all the dates, variants, and animals 
by_epitope_dfs = list()
for(a in 1:length(by_epitope)){
  by_epitope_dfs[[a]] = do.call(rbind.fill, by_epitope[[a]])
}
names(by_epitope_dfs) = names(by_epitope)

by_epitope_wl = list()
for(a in epitopes){
  for(b in 1:length(epitope_repcomb_withlow)){
    if(!is.null(epitope_repcomb_withlow[[b]][[a]])){
      if(nrow(epitope_repcomb_withlow[[b]][[a]])> 0){
        by_epitope_wl[[a]][[b]] = epitope_repcomb_withlow[[b]][[a]]
        by_epitope_wl[[a]][[b]]$dpi = rownames(by_epitope_wl[[a]][[b]])
        by_epitope_wl[[a]][[b]]$animal = names(epitope_repcomb_withlow)[b]
        names(by_epitope_wl[[a]])[b] = names(epitope_repcomb_withlow)[b]
      }
    }
  }
}

by_epitope_wl_dfs = list()
for(a in 1:length(by_epitope_wl)){
  by_epitope_wl_dfs[[a]] = do.call(rbind.fill, by_epitope_wl[[a]])
}
names(by_epitope_wl_dfs) = names(by_epitope_wl)

#### Right now we only care about GW9, RM9, and SP10, so pull those out ####
eoi = c("Gag_GW9", "Nef_RM9", "Rev_SP10")
by_epitope_dfs_mini = by_epitope_dfs[grep(str_flatten(eoi, "|"), names(by_epitope_dfs))]

#### Convert proportions to log10 copies based on the sample viral load #### 
epitope_repcomb_freqfilter = list()
thresholds = c(0, 0.01, 0.03, 0.05)
for(a in 1:length(epitope_repcomb_reformat)){
  threshold_list = list()
  for(b in 1:length(epitope_repcomb_reformat[[a]])){
    df = epitope_repcomb_reformat[[a]][[b]][, grep("other", colnames(epitope_repcomb_reformat[[a]][[b]]), invert = T)]
    
    tlist = list()
    for(c in thresholds){
      tdf = df 
      if(!is.null(ncol(tdf)) & length(tdf) > 1){
        for(d in 1:ncol(tdf)){
          if(all(tdf[,d] < c, na.rm = T)){
            tdf[, d] = NA
          }
        }
        tdf[, "other"] = 1-rowSums(tdf, na.rm = T)
      }else{
        tdf = epitope_repcomb_reformat[[a]][[b]]
      }
      tdf = tdf[, which(colSums(tdf, na.rm = T) > 0)]
      
      tlist[[as.character(c)]] = tdf
    }
    threshold_list[[b]] = tlist
  }
  names(threshold_list) = names(epitope_repcomb_reformat[[a]])
  epitope_repcomb_freqfilter[[a]] = threshold_list
}
names(epitope_repcomb_freqfilter) = names(epitope_repcomb)

epitope_copies = list()
for(a in 1:length(epitope_repcomb_freqfilter)){
  
  if(names(epitope_repcomb_freqfilter)[a] %in% lateART){
    animal_vls = lateART_vls[, c("Days.post.SIV", names(epitope_repcomb_freqfilter)[a])]
  }else{
    animal_vls = viralloads[which(viralloads$animal == names(epitope_repcomb_freqfilter)[a] & viralloads$sampletype == "Plasma" ), ]
    colnames(animal_vls) = c("animal", "sampletype", "date", names(epitope_repcomb_freqfilter)[a], "Days.post.SIV")
  }
  
  copylists = list()
  for(b in 1:length(epitope_repcomb_freqfilter[[a]])){
    sublist = list()
    for(d in 1:length(epitope_repcomb_freqfilter[[a]][[b]])){
      df = epitope_repcomb_freqfilter[[a]][[b]][[d]]
      
      if(is.data.frame(df)){
        if(nrow(df) > 0){
          df$dpi = as.numeric(str_extract(rownames(df), "\\d{1,3}"))
          copydf = df
          for(c in 1:nrow(df)){
            copydf[c, grep("dpi", colnames(copydf), invert = T)] = df[c, grep("dpi", colnames(df), invert = T)]*animal_vls[which(animal_vls$Days.post.SIV == df[c, "dpi"]), names(epitope_repcomb_freqfilter)[a]]
            copydf[c, "vl"] = animal_vls[which(animal_vls$Days.post.SIV == df[c, "dpi"]), names(epitope_repcomb_freqfilter)[a]]
          }
          
          copydf = copydf[, which(colSums(df, na.rm = T) > 0)]
          sublist[[names(epitope_repcomb_freqfilter[[a]][[b]][d])]] = copydf
        }else{
          sublist[[names(epitope_repcomb_freqfilter[[a]][[b]][d])]] = NULL
        }
      }
      
    }
    copylists[[names(epitope_repcomb_freqfilter[[a]][b])]] = sublist 
  }
  epitope_copies[[a]] = copylists
}
names(epitope_copies) = names(epitope_repcomb_reformat)

epitope_copies_log10 = vector("list", length(epitope_copies))
for(a in 1:length(epitope_copies)){
  for(b in 1:length(epitope_copies[[a]])){
    sublist = list()
    if(length(epitope_copies[[a]][[b]]) > 0){
      for(d in 1:length(epitope_copies[[a]][[b]])){
        
        df = epitope_copies[[a]][[b]][[d]][, grep("dpi", colnames(epitope_copies[[a]][[b]][[d]]), invert = T)]
        
        if(!is.null(nrow(df))){
          log10df = log10(df)
          total_vl = log10(rowSums(df, na.rm = T))
          for(c in 1:nrow(log10df)){
            log10df[c, which(log10df[c, ] < 0)] = 0
          }
          log10df$dpi = as.numeric(str_extract(rownames(log10df), "\\d{1,3}"))
          log10df$total_vl = total_vl
        }else{
          log10df = log10(epitope_copies[[a]][[b]][[d]])
          log10df$dpi = epitope_copies[[a]][[b]][[d]]$dpi
          log10df$total_vl = log10df[, 1]
        }
        
        sublist[[d]] = log10df
      }
      names(sublist) = names(epitope_copies[[a]][[b]])
      epitope_copies_log10[[a]][[names(epitope_copies[[a]][b])]] = sublist
    }else{
      
    }
    
  }
}
names(epitope_copies_log10) = names(epitope_copies)

ec_oi_log10 = list()
for(a in 1:length(epitope_copies_log10)){
  ec_oi_log10[[a]] = epitope_copies_log10[[a]][grep(str_flatten(eoi, "|"), names(epitope_copies_log10[[a]]))]
}
names(ec_oi_log10) = names(epitope_copies_log10)

ec_oi_prop = list()
for(a in 1:length(epitope_repcomb_freqfilter)){
  ec_oi_prop[[a]] = epitope_repcomb_freqfilter[[a]][grep(str_flatten(eoi, "|"), names(epitope_repcomb_freqfilter[[a]]))]
}
names(ec_oi_prop) = names(epitope_repcomb_freqfilter)

#### Do all the above, but with the low txts included #### 
epitope_repcomb_freqfilter_wl = list()
for(a in 1:length(epitope_repcomb_reformat_wl)){
  threshold_list = list()
  for(b in 1:length(epitope_repcomb_reformat_wl[[a]])){
    df = epitope_repcomb_reformat_wl[[a]][[b]][, grep("other", colnames(epitope_repcomb_reformat_wl[[a]][[b]]), invert = T)]
    
    tlist = list()
    for(c in thresholds){
      tdf = df 
      if(!is.null(ncol(tdf)) & length(tdf) > 1){
        for(d in 1:ncol(tdf)){
          if(all(tdf[,d] < c, na.rm = T)){
            tdf[, d] = NA
          }
        }
        tdf[, "other"] = 1-rowSums(tdf, na.rm = T)
      }else{
        tdf = epitope_repcomb_reformat_wl[[a]][[b]]
      }
      tdf = tdf[, which(colSums(tdf, na.rm = T) > 0)]
      
      tlist[[as.character(c)]] = tdf
    }
    threshold_list[[b]] = tlist
  }
  names(threshold_list) = names(epitope_repcomb_reformat_wl[[a]])
  epitope_repcomb_freqfilter_wl[[a]] = threshold_list
}
names(epitope_repcomb_freqfilter_wl) = names(epitope_repcomb_withlow)

epitope_copies_wl = list()
for(a in 1:length(epitope_repcomb_freqfilter_wl)){
  
  if(names(epitope_repcomb_freqfilter_wl)[a] %in% lateART){
    animal_vls = lateART_vls[, c("Days.post.SIV", names(epitope_repcomb_freqfilter_wl)[a])]
  }else{
    animal_vls = viralloads[which(viralloads$animal == names(epitope_repcomb_freqfilter_wl)[a] & viralloads$sampletype == "Plasma" ), ]
    colnames(animal_vls) = c("animal", "sampletype", "date", names(epitope_repcomb_freqfilter_wl)[a], "Days.post.SIV")
  }
  
  copylists = list()
  for(b in 1:length(epitope_repcomb_freqfilter_wl[[a]])){
    sublist = list()
    for(d in 1:length(epitope_repcomb_freqfilter_wl[[a]][[b]])){
      df = epitope_repcomb_freqfilter_wl[[a]][[b]][[d]]
      
      if(is.data.frame(df)){
        if(nrow(df) > 0){
          df$dpi = as.numeric(str_extract(rownames(df), "\\d{1,3}"))
          copydf = df
          for(c in 1:nrow(df)){
            copydf[c, grep("dpi", colnames(copydf), invert = T)] = df[c, grep("dpi", colnames(df), invert = T)]*animal_vls[which(animal_vls$Days.post.SIV == df[c, "dpi"]), names(epitope_repcomb_freqfilter_wl)[a]]
            copydf[c, "vl"] = animal_vls[which(animal_vls$Days.post.SIV == df[c, "dpi"]), names(epitope_repcomb_freqfilter_wl)[a]]
          }
          
          copydf = copydf[, which(colSums(df, na.rm = T) > 0)]
          sublist[[names(epitope_repcomb_freqfilter_wl[[a]][[b]][d])]] = copydf
        }else{
          sublist[[names(epitope_repcomb_freqfilter_wl[[a]][[b]][d])]] = NULL
        }
      }
      
    }
    copylists[[names(epitope_repcomb_freqfilter_wl[[a]][b])]] = sublist 
  }
  epitope_copies_wl[[a]] = copylists
}
names(epitope_copies_wl) = names(epitope_repcomb_reformat_wl)

epitope_copies_log10_wl = vector("list", length(epitope_copies_wl))
for(a in 1:length(epitope_copies_wl)){
  for(b in 1:length(epitope_copies_wl[[a]])){
    sublist = list()
    if(length(epitope_copies_wl[[a]][[b]]) > 0){
      for(d in 1:length(epitope_copies_wl[[a]][[b]])){
        
        df = epitope_copies_wl[[a]][[b]][[d]][, grep("dpi", colnames(epitope_copies_wl[[a]][[b]][[d]]), invert = T)]
        
        if(!is.null(nrow(df))){
          log10df = log10(df)
          total_vl = log10(rowSums(df, na.rm = T))
          for(c in 1:nrow(log10df)){
            log10df[c, which(log10df[c, ] < 0)] = 0
          }
          log10df$dpi = as.numeric(str_extract(rownames(log10df), "\\d{1,3}"))
          log10df$total_vl = total_vl
        }else{
          log10df = log10(epitope_copies_wl[[a]][[b]][[d]])
          log10df$dpi = epitope_copies_wl[[a]][[b]][[d]]$dpi
          log10df$total_vl = log10df[, 1]
        }
        
        sublist[[d]] = log10df
      }
      names(sublist) = names(epitope_copies_wl[[a]][[b]])
      epitope_copies_log10_wl[[a]][[names(epitope_copies_wl[[a]][b])]] = sublist
    }else{
      
    }
    
  }
}
names(epitope_copies_log10_wl) = names(epitope_copies_wl)

ec_oi_log10_wl = list()
for(a in 1:length(epitope_copies_log10_wl)){
  ec_oi_log10_wl[[a]] = epitope_copies_log10_wl[[a]][grep(str_flatten(eoi, "|"), names(epitope_copies_log10_wl[[a]]))]
}
names(ec_oi_log10_wl) = names(epitope_copies_log10_wl)

ec_oi_prop_wl = list()
for(a in 1:length(epitope_repcomb_freqfilter_wl)){
  ec_oi_prop_wl[[a]] = epitope_repcomb_freqfilter_wl[[a]][grep(str_flatten(eoi, "|"), names(epitope_repcomb_freqfilter_wl[[a]]))]
}
names(ec_oi_prop_wl) = names(epitope_repcomb_freqfilter_wl)

#### Complete chisquare tests to compare the pre-ART and post-ART epitope sequences ####
prevspostART_times = data.frame("cy0719" = c(49, 659), "cy0721" = c(49, 693), 
                                "cy0723" = c(49, 715), "cy0724" = c(49, 686), 
                                "cy0726" = c(49, 686), 
                                "cy1036" = c(11, 507),"cy1039" = c(11, 504),
                                "cy1044" = c(14, 504), 
                                "cy1035" = c(14, 288), "cy1040" = c(14, 316), 
                                "cy1043" = c(11, 336), "cy1045" = c(11, 470))

chisq_data_wl = data.frame()
for(a in 1:length(ec_oi_log10_wl)){
  animal = names(ec_oi_log10_wl)[a]
  for(b in 1:length(ec_oi_log10_wl[[a]])){
    for(c in 1:length(ec_oi_log10_wl[[a]][[b]])){
      df = ec_oi_log10_wl[[a]][[b]][[c]]
      
      dates = prevspostART_times[, animal]
      pre_row = which(df$dpi == dates[1])
      post_row = which(df$dpi == dates[2])
      df = df[, which(colSums(df[c(pre_row, post_row), grep("dpi|vl", colnames(df),invert = T)], na.rm = T) > 0)]
      
      if(is.data.frame(df)){
        pre = ZeroIfNA(as.numeric(df[pre_row, grep("dpi|vl", colnames(df), invert = T)]))
        post = ZeroIfNA(as.numeric(df[post_row, grep("dpi|vl", colnames(df), invert = T)]))
        
        df_prop = ec_oi_prop_wl[[a]][[b]][[c]]
        df_prop$dpi = as.numeric(str_extract(rownames(df_prop), "\\d{1,3}"))
        df_prop = df_prop[, which(colSums(df_prop[c(pre_row, post_row),grep("dpi|vl", colnames(df_prop), invert = T)], na.rm = T) > 0)]
        
        pre_prop = ZeroIfNA(as.numeric(df_prop[pre_row, grep("dpi|vl", colnames(df_prop), invert = T)]))
        post_prop = ZeroIfNA(as.numeric(df_prop[post_row, grep("dpi|vl", colnames(df_prop), invert = T)]))
        
        if(length(pre) & length(post) != 0 & !all(is.nan(pre)) & !all(is.nan(post)) & length(pre_row) > 0 & length(post_row) > 0){
          chsq = chisq.test(post, p = 10^pre, simulate.p.value = T, rescale.p = T)
          chisq_data_wl[paste0(animal, "-", names(ec_oi_log10_wl[[a]])[b], "-", names(ec_oi_log10_wl[[a]][[b]])[c]), "xsq_copies"] = chsq$statistic
          chisq_data_wl[paste0(animal, "-", names(ec_oi_log10_wl[[a]])[b], "-", names(ec_oi_log10_wl[[a]][[b]])[c]), "pval_copies"] = chsq$p.value
          
          chsq_prop = chisq.test(post_prop, p = pre_prop, rescale.p = T)
          chisq_data_wl[paste0(animal, "-", names(ec_oi_log10_wl[[a]])[b], "-", names(ec_oi_log10_wl[[a]][[b]])[c]), "xsq_prop"] = chsq_prop$statistic
          chisq_data_wl[paste0(animal, "-", names(ec_oi_log10_wl[[a]])[b], "-", names(ec_oi_log10_wl[[a]][[b]])[c]), "pval_prop"] = chsq_prop$p.value
        }
      }
      
    }
  }
}
epitopes = epitope_seq_ref$EpitopeName
chisq_data_wl$animal = str_extract(rownames(chisq_data_wl), animal_pattern)
chisq_data_wl$epitope = str_extract(rownames(chisq_data_wl), str_flatten(epitopes, "|"))
chisq_data_wl$freqfilter = str_extract(rownames(chisq_data_wl), str_flatten(c(0.01, 0.03, 0.05), "|"))

chisq_data_eoi = chisq_data[grep("GW9|RM9|SP10", chisq_data$epitope), ]

#### There were missing values in the above, so let's try to fill those in now ####
cy0721_rm9_xsq = epitope_freq_recalc$cy0721$Nef_RM9[c(2,3), c(2:6)] 
cy0721_rm9_xsq_cp = epitope_copies_log10$cy0721$Nef_RM9[c(2,3), c(2:6)]
cy0721_rm9_xsq[1, 4] = 0.001 # Failed because one of the epitopes post-ART is not present pre-ART 

chisq_data_eoi["cy0721-Nef_RM9", "xsq_prop"] = chisq.test(as.numeric(cy0721_rm9_xsq[2, ]), p = as.numeric(cy0721_rm9_xsq[1, ]), rescale.p = T, simulate.p.value = T)$statistic
chisq_data_eoi["cy0721-Nef_RM9", "pval_prop"] = chisq.test(as.numeric(cy0721_rm9_xsq[2, ]), p = as.numeric(cy0721_rm9_xsq[1, ]), rescale.p = T, simulate.p.value = T)$p.value

chisq_data_eoi["cy0721-Nef_RM9", "xsq_copies"] = chisq.test(as.numeric(cy0721_rm9_xsq_cp[2, ]), p = as.numeric(cy0721_rm9_xsq[1, ]), rescale.p = T)$statistic
chisq_data_eoi["cy0721-Nef_RM9", "pval_copies"] = chisq.test(as.numeric(cy0721_rm9_xsq_cp[2, ]), p = as.numeric(cy0721_rm9_xsq[1, ]), rescale.p = T)$p.value

cy0726_gw9_xsq = epitope_freq_recalc$cy0726$Gag_GW9[c(2,3), ]
cy0726_gw9_xsq_cp = epitope_copies_log10$cy0726$Gag_GW9[c(2,3), c(1:5)]
cy0726_gw9_xsq[1, 4] = 0.001

chisq_data_eoi["cy0726-Gag_GW9", "xsq_prop"] = chisq.test(as.numeric(cy0726_gw9_xsq[2, ]), p = as.numeric(cy0726_gw9_xsq[1, ]), rescale.p = T, simulate.p.value = T)$statistic
chisq_data_eoi["cy0726-Gag_GW9", "pval_prop"] = chisq.test(as.numeric(cy0726_gw9_xsq[2, ]), p = as.numeric(cy0726_gw9_xsq[1, ]), rescale.p = T, simulate.p.value = T)$p.value

chisq_data_eoi["cy0726-Gag_GW9", "xsq_copies"] = chisq.test(as.numeric(cy0726_gw9_xsq_cp[2, ]), p = as.numeric(cy0726_gw9_xsq[1, ]), rescale.p = T)$statistic
chisq_data_eoi["cy0726-Gag_GW9", "pval_copies"] = chisq.test(as.numeric(cy0726_gw9_xsq_cp[2, ]), p = as.numeric(cy0726_gw9_xsq[1, ]), rescale.p = T)$p.value

cy0726_sp10_xsq = epitope_freq_recalc$cy0726$Rev_SP10[c(2,3), ]
cy0726_sp10_xsq_cp = epitope_copies_log10$cy0726$Rev_SP10[c(2,3), c(1:8)]
cy0726_sp10_xsq[1, c(6,7)] = 0.001 

chisq_data_eoi["cy0726-Rev_SP10", "xsq_prop"] = chisq.test(as.numeric(cy0726_sp10_xsq[2, ]), p = as.numeric(cy0726_sp10_xsq[1, ]), rescale.p = T, simulate.p.value = T)$statistic
chisq_data_eoi["cy0726-Rev_SP10", "pval_prop"] = chisq.test(as.numeric(cy0726_sp10_xsq[2, ]), p = as.numeric(cy0726_sp10_xsq[1, ]), rescale.p = T, simulate.p.value = T)$p.value

chisq_data_eoi["cy0726-Rev_SP10", "xsq_copies"] = chisq.test(as.numeric(cy0726_sp10_xsq_cp[2, ]), p = as.numeric(cy0726_sp10_xsq[1, ]), rescale.p = T)$statistic
chisq_data_eoi["cy0726-Rev_SP10", "pval_copies"] = chisq.test(as.numeric(cy0726_sp10_xsq_cp[2, ]), p = as.numeric(cy0726_sp10_xsq[1, ]), rescale.p = T)$p.value

cy1039_sp10_xsq = epitope_freq_recalc$cy1039$Rev_SP10[c(1,3), c(1,2,6)]
cy1039_sp10_xsq_cp = epitope_copies_log10$cy1039$Rev_SP10[c(1,3), c(1,2,6)]
cy1039_sp10_xsq[1, 2] = 0.001 

chisq_data_eoi["cy1039-Rev_SP10", "xsq_prop"] = chisq.test(as.numeric(cy1039_sp10_xsq[2, ]), p = as.numeric(cy1039_sp10_xsq[1, ]), rescale.p = T, simulate.p.value = T)$statistic
chisq_data_eoi["cy1039-Rev_SP10", "pval_prop"] = chisq.test(as.numeric(cy1039_sp10_xsq[2, ]), p = as.numeric(cy1039_sp10_xsq[1, ]), rescale.p = T, simulate.p.value = T)$p.value

chisq_data_eoi["cy1039-Rev_SP10", "xsq_copies"] = chisq.test(as.numeric(cy1039_sp10_xsq_cp[2, ]), p = as.numeric(cy1039_sp10_xsq[1, ]), rescale.p = T)$statistic
chisq_data_eoi["cy1039-Rev_SP10", "pval_copies"] = chisq.test(as.numeric(cy1039_sp10_xsq_cp[2, ]), p = as.numeric(cy1039_sp10_xsq[1, ]), rescale.p = T)$p.value

cy1043_rm9_xsq = epitope_freq_recalc$cy1043$Nef_RM9[c(1,2), c(1,2,9)]
cy1043_rm9_xsq_cp = epitope_copies_log10$cy1043$Nef_RM9[c(1,2), c(1,2,9)]
cy1043_rm9_xsq[1, 2] = 0.001 

chisq_data_eoi["cy1043-Nef_RM9", "xsq_prop"] = chisq.test(as.numeric(cy1043_rm9_xsq[2, ]), p = as.numeric(cy1043_rm9_xsq[1, ]), rescale.p = T, simulate.p.value = T)$statistic
chisq_data_eoi["cy1043-Nef_RM9", "pval_prop"] = chisq.test(as.numeric(cy1043_rm9_xsq[2, ]), p = as.numeric(cy1043_rm9_xsq[1, ]), rescale.p = T, simulate.p.value = T)$p.value

chisq_data_eoi["cy1043-Nef_RM9", "xsq_copies"] = chisq.test(as.numeric(cy1043_rm9_xsq_cp[2, ]), p = as.numeric(cy1043_rm9_xsq[1, ]), rescale.p = T)$statistic
chisq_data_eoi["cy1043-Nef_RM9", "pval_copies"] = chisq.test(as.numeric(cy1043_rm9_xsq_cp[2, ]), p = as.numeric(cy1043_rm9_xsq[1, ]), rescale.p = T)$p.value

cy1044_rm9_xsq = epitope_freq_recalc$cy1044$Nef_RM9[c(1,4), c(1,3,7)]
cy1044_rm9_xsq_cp = epitope_copies_log10$cy1044$Nef_RM9[c(1,4), c(1,3,7)]
cy1044_rm9_xsq[1, 2] = 0.001 

chisq_data_eoi["cy1044-Nef_RM9", "xsq_prop"] = chisq.test(as.numeric(cy1044_rm9_xsq[2, ]), p = as.numeric(cy1044_rm9_xsq[1, ]), rescale.p = T, simulate.p.value = T)$statistic
chisq_data_eoi["cy1044-Nef_RM9", "pval_prop"] = chisq.test(as.numeric(cy1044_rm9_xsq[2, ]), p = as.numeric(cy1044_rm9_xsq[1, ]), rescale.p = T, simulate.p.value = T)$p.value

chisq_data_eoi["cy1044-Nef_RM9", "xsq_copies"] = chisq.test(as.numeric(cy1044_rm9_xsq_cp[2, ]), p = as.numeric(cy1044_rm9_xsq[1, ]), rescale.p = T)$statistic
chisq_data_eoi["cy1044-Nef_RM9", "pval_copies"] = chisq.test(as.numeric(cy1044_rm9_xsq_cp[2, ]), p = as.numeric(cy1044_rm9_xsq[1, ]), rescale.p = T)$p.value

cy1045_rm9_xsq = epitope_freq_recalc$cy1045$Nef_RM9[c(1,4), c(1,3,4,6)]
cy1045_rm9_xsq_cp = epitope_copies_log10$cy1045$Nef_RM9[c(1,4), c(1,3,4,6)]
cy1045_rm9_xsq[1, c(2,3)] = 0.001 

chisq_data_eoi["cy1045-Nef_RM9", "xsq_prop"] = chisq.test(as.numeric(cy1045_rm9_xsq[2, ]), p = as.numeric(cy1045_rm9_xsq[1, ]), rescale.p = T, simulate.p.value = T)$statistic
chisq_data_eoi["cy1045-Nef_RM9", "pval_prop"] = chisq.test(as.numeric(cy1045_rm9_xsq[2, ]), p = as.numeric(cy1045_rm9_xsq[1, ]), rescale.p = T, simulate.p.value = T)$p.value

chisq_data_eoi["cy1045-Nef_RM9", "xsq_copies"] = chisq.test(as.numeric(cy1045_rm9_xsq_cp[2, ]), p = as.numeric(cy1045_rm9_xsq[1, ]), rescale.p = T)$statistic
chisq_data_eoi["cy1045-Nef_RM9", "pval_copies"] = chisq.test(as.numeric(cy1045_rm9_xsq_cp[2, ]), p = as.numeric(cy1045_rm9_xsq[1, ]), rescale.p = T)$p.value

cy1045_gw9_xsq = epitope_freq_recalc$cy1045$Gag_GW9[c(1,4), c(1,2,4)]
cy1045_gw9_xsq_cp = epitope_copies_log10$cy1045$Gag_GW9[c(1,4), c(1,2,4)]
cy1045_gw9_xsq[1, 2] = 0.001 

chisq_data_eoi["cy1045-Gag_GW9", "xsq_prop"] = chisq.test(as.numeric(cy1045_gw9_xsq[2, ]), p = as.numeric(cy1045_gw9_xsq[1, ]), rescale.p = T, simulate.p.value = T)$statistic
chisq_data_eoi["cy1045-Gag_GW9", "pval_prop"] = chisq.test(as.numeric(cy1045_gw9_xsq[2, ]), p = as.numeric(cy1045_gw9_xsq[1, ]), rescale.p = T, simulate.p.value = T)$p.value

chisq_data_eoi["cy1045-Gag_GW9", "xsq_copies"] = chisq.test(as.numeric(cy1045_gw9_xsq_cp[2, ]), p = as.numeric(cy1045_gw9_xsq[1, ]), rescale.p = T)$statistic
chisq_data_eoi["cy1045-Gag_GW9", "pval_copies"] = chisq.test(as.numeric(cy1045_gw9_xsq_cp[2, ]), p = as.numeric(cy1045_gw9_xsq[1, ]), rescale.p = T)$p.value


#### After looking at the chisq data, we are looking back 

#### Calculate bray-curtis dissimilarity for pre- and post-ART epitope sequences ####
braycurtis_list = list()
shannon_df = data.frame()
for(a in 1:length(epitope_freq_recalc)){
  nonempty = epitope_freq_recalc[[a]][which(lengths(epitope_freq_recalc[[a]]) > 0)]
  sublist = list()
  for(b in 1:length(nonempty)){
    edf = ZeroIfNA(nonempty[[b]])
    edf = edf[which(rowSums(edf) > 0), grep("total_counts", colnames(edf), invert = T)]
    edf = edf[grep("442", rownames(edf), invert = T), ]
    
    sublist[[b]] = vegdist(edf, method = "bray") 
    
    shannon_df = rbind.fill(shannon_df, data.frame(cbind(diversity(edf, "shannon"), dpi = rownames(edf), epitope = names(nonempty[b]), animal = names(epitope_freq_recalc[a]))))
  }
  bcdf$animal = names(by_animal_epitope_freqs_repcomb[a])
  
  names(sublist) = names(nonempty)
  braycurtis_list[[a]] = sublist
}
names(braycurtis_list) = names(epitope_freq_recalc)
shannon_df[which(shannon_df$animal %in% mcm8wk), "cohort"] = "mcm8wk"
shannon_df[which(shannon_df$animal %in% mcm2wk_vir), "cohort"] = "viremic"
shannon_df[which(shannon_df$animal %in% mcm2wk_avir), "cohort"] = "aviremic"
shannon_df = shannon_df[which(shannon_df$animal != "cy1037"), ]

shannon_byepitope = list()
for(a in 1:length(unique(shannon_df$epitope))){
  shannon_byepitope[[a]] = shannon_df[which(shannon_df$epitope == unique(shannon_df$epitope)[a]), ]
  shannon_byepitope[[a]]$V1 = as.numeric(shannon_byepitope[[a]]$V1)
  shannon_byepitope[[a]]$dpi = as.numeric(shannon_byepitope[[a]]$dpi)
}
names(shannon_byepitope) = unique(shannon_df$epitope)

shannon_byepitope_reorg = list()
for(a in 1:length(shannon_byepitope)){
  alldpi = sort(unique(shannon_byepitope[[a]]$dpi))
  org_animals = c(mcm2wk_avir, mcm2wk_vir, mcm8wk)
  df = data.frame()
  for(b in 1:length(org_animals)){
    for(c in 1:length(alldpi)){
      check.exists = shannon_byepitope[[a]][which(shannon_byepitope[[a]]$animal == org_animals[b] & shannon_byepitope[[a]]$dpi == alldpi[c]), "V1"]
      if(is.numeric(check.exists) & length(check.exists) == 1){
        df[paste0("d", alldpi[c]), org_animals[b]] = check.exists
      }else{
        df[paste0("d", alldpi[c]), org_animals[b]] = NA
      }
    }
  }
  shannon_byepitope_reorg[[a]] = df 
}
names(shannon_byepitope_reorg) = names(shannon_byepitope)

intervals_8wk = list(preART_8wk = c(0,56), postART_8wk = c(600, 800))
intervals_2wk = list(preART_2wk = c(0, 14), postART_2wk = c(260, 441), postRC = c(443, 496), postDep = c(497, 560))


braycurtis_2wk = braycurtis_list[which(names(braycurtis_list) %in% mcm2wk)]
braycurtis_8wk = braycurtis_list[which(names(braycurtis_list) %in% mcm8wk)]


bc_reorg_df = data.frame()
for(a in 1:length(braycurtis_2wk)){
  animaldf = data.frame()
  for(b in 1:length(braycurtis_2wk[[a]])){
    epitopedf = as.matrix(braycurtis_2wk[[a]][[b]])
    z = 1 
    tdf = data.frame()
    for(c in 1:ncol(epitopedf)){
      for(d in 1:nrow(epitopedf)){
        if(colnames(epitopedf)[c] != rownames(epitopedf)[d]){
          tdf[z, "t1"] = colnames(epitopedf)[c]
          tdf[z, "t2"] = rownames(epitopedf)[d]
          tdf[z, "braycurtis"] = epitopedf[d, c]
          z = z + 1 
        }
        
        
      }
    }
    tdf$epitope = names(braycurtis_2wk[[a]][b])
    animaldf = rbind(animaldf, tdf)
  }
  animaldf$animal = names(braycurtis_2wk[a])
  bc_reorg_df = rbind(bc_reorg_df, animaldf)
}

bc_reorg_df_mini = bc_reorg_df[grep("GW9|RM9|SP10", bc_reorg_df$epitope), ]

braycurtis_byepitope = list()
for(a in 1:length(epitopes)){
  epitope_rows = braycurtis_df[which(braycurtis_df$epitope == epitopes[a]), ]
}

