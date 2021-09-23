# RO analysis subtractiong the S0
library(tidyverse)
library(foreach)
library(doSNOW)
library(mixtools)

  #Run 150521
  s0 = read_delim(
    '/Volumes/Storage2/RT_DATA/Old_DATA/Repliseq_profiles_50kbadjacent_15Mreads/Run_150429-150521/UnsortedNTrep1.sorted.mappUni.delDupl.bam.bg',
    skip = 1,
    col_names = c('chr', 'start', 'end', 'S0'),
    delim = ' '
  )
  
  th_first_peak = function(x, cont) {
    cont = cont[x > 0]
    x = x[x > 0]
    x = x - cont
    
    library(mixtools)
    
    m = normalmixEM(x, k = 2, maxit = 2000)
    m = min(m$mu)
    
    return(m)
  }
  #load files
  dir = c(
    '/Volumes/Storage2/RT_DATA/RO_data_second_analysis/50kb_adj_15M/',
    '/Volumes/Storage2/RT_DATA/New_DATA/Run_200214/50kb_adj_15M/'
  )
  
  AphRO = foreach(
    i = 1:2,
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach')
  ) %do% {
    Files = list.files(dir[i])
    Files = Files[grep(x = Files, pattern = 'Rep5-AphiRO16h|aph-RO-16h')]
    
    foreach(File = Files,
            .combine = 'rbind',
            .packages = 'tidyverse') %do% {
              phase = read_delim(
                paste0(dir[i], File),
                col_names = c('chr', 'start', 'end', 'reads'),
                delim = ' ',
                skip = 1
              ) %>%
                mutate(
                  Condition = paste0('AphRO', i),
                  phase = str_extract(File, pattern = '[SG][1-4]'),
                  rep = i
                ) %>%
                inner_join(s0) %>%
                mutate(reads_s0 = reads - S0) %>%
                group_by(Condition, rep, phase) %>%
                mutate(reads = ifelse(reads > quantile(reads, .9999)[[1]], 0, reads),
                       th = th_first_peak(reads, S0))
              
              d = density(phase$reads_s0[phase$reads > 0])
              Th = unique(phase$th)
              
              TITLE = paste('AphRO', i, str_extract(File, pattern = '[SG][1-4]'))
              sub_TITLE = paste('Th = ', Th)
              plot(d, main = TITLE, sub = sub_TITLE)
              abline(v = Th, col = "purple")
              print(Th)
              
              input = NA
              while (input != '' | is.na(input)) {
                input = readline()
                Th = as.numeric(input)
                phase = phase %>%
                  mutate(th = ifelse(!is.na(Th), Th, th))
                Th = unique(phase$th)
                
                sub_TITLE = paste('Th = ', Th)
                plot(d, main = TITLE, sub = sub_TITLE)
                abline(v = Th, col = "purple")
                
                
                
              }
              phase %>% mutate(
                reads_s0 = reads_s0 - th,
                reads_s0 = ifelse(reads <= 0 |
                                    reads_s0 <= 0 , 0, reads_s0),
                reads = reads_s0
              ) %>%
                ungroup()
            }
  }
  
  dir = c(
    '/Volumes/Storage2/RT_DATA/RO_data_second_analysis/50kb_adj_15M/',
    '/Volumes/Storage2/RT_DATA/Old_DATA/Repliseq_profiles_50kbadjacent_15Mreads//Run_171107/',
    '/Volumes/Storage2/RT_DATA/Old_DATA/Repliseq_profiles_50kbadjacent_15Mreads/Run_15-381_run160308/',
    '/Volumes/Storage2/RT_DATA/Old_DATA/Repliseq_profiles_50kbadjacent_15Mreads/Run_15-235-1/'
  )
  
  
  Aph = foreach(
    i = 1:4,
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach')
  ) %do% {
    Files = list.files(dir[i])
    Files = Files[grepl(x = Files, pattern = 'Aph16hrep|Aph-600nM-16h|Aph16hDT4')]
    foreach(File = Files,
            .combine = 'rbind',
            .packages = 'tidyverse') %do% {
              phase = read_delim(
                paste0(dir[i], File),
                skip = 1,
                col_names = c('chr', 'start', 'end', 'reads'),
                delim = ' '
              ) %>%
                mutate(
                  Condition = paste0('Aph', i),
                  phase = str_extract(File, pattern = '[SG][1-4]'),
                  rep = i
                ) %>%
                inner_join(s0) %>%
                mutate(reads_s0 = reads - S0) %>%
                group_by(Condition, rep, phase) %>%
                mutate(reads = ifelse(reads > quantile(reads, .9999)[[1]], 0, reads),
                       th = th_first_peak(reads, S0))
              
              d = density(phase$reads_s0[phase$reads > 0])
              Th = unique(phase$th)
              
              TITLE = paste('Aph', i, str_extract(File, pattern = '[SG][1-4]'))
              sub_TITLE = paste('Th = ', Th)
              plot(d, main = TITLE, sub = sub_TITLE)
              abline(v = Th, col = "purple")
              print(Th)
              
              input = NA
              while (input != '' | is.na(input)) {
                input = readline()
                Th = as.numeric(input)
                phase = phase %>%
                  mutate(th = ifelse(!is.na(Th), Th, th))
                Th = unique(phase$th)
                
                sub_TITLE = paste('Th = ', Th)
                plot(d, main = TITLE, sub = sub_TITLE)
                abline(v = Th, col = "purple")
                
                
                
              }
              
              phase %>% mutate(
                reads_s0 = reads_s0 - th,
                reads_s0 = ifelse(reads <= 0 |
                                      reads_s0 <= 0 , 0, reads_s0),
                reads = reads_s0
              ) %>%
                ungroup()
              
            }
  }
  
  dir = c(
    '/Volumes/Storage2/RT_DATA/RO_data_second_analysis/50kb_adj_15M/',
    '/Volumes/Storage2/RT_DATA/Old_DATA/Repliseq_profiles_50kbadjacent_15Mreads/Run_15-381-2/',
    '/Volumes/Storage2/RT_DATA/New_DATA/Run_180719/50kb_adj_15M/',
    '/Volumes/Storage2/RT_DATA/Old_DATA/Repliseq_profiles_50kbadjacent_15Mreads/Run_150429-150521/'
  )
  
  NT = foreach(
    i = 1:4,
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach')
  ) %do% {
    Files = list.files(dir[i])
    Files = Files[grepl(x = Files, pattern = 'NT.*[GS]')]
    Files = Files[!grepl(x = Files, pattern = 'Aph2h')]
    
    foreach(File = Files,
            .combine = 'rbind',
            .packages = 'tidyverse') %do% {
              phase = read_delim(
                paste0(dir[i], File),
                skip = 1,
                col_names = c('chr', 'start', 'end', 'reads'),
                delim = ' '
              ) %>%
                mutate(
                  Condition = paste0('NT', i),
                  phase = str_extract(File, pattern = '[SG][1-4]'),
                  rep = i
                ) %>%
                inner_join(s0) %>%
                mutate(reads_s0 = reads - S0) %>%
                group_by(Condition, rep, phase) %>%
                mutate(reads = ifelse(reads > quantile(reads, .9999)[[1]], 0, reads),
                       th = th_first_peak(reads, S0))
              
              d = density(phase$reads_s0[phase$reads > 0])
              Th = unique(phase$th)
              
              TITLE = paste('NT', i, str_extract(File, pattern = '[SG][1-4]'))
              sub_TITLE = paste('Th = ', Th)
              plot(d, main = TITLE, sub = sub_TITLE)
              abline(v = Th, col = "purple")
              print(Th)
              
              input = NA
              while (input != '' | is.na(input)) {
                input = readline()
                Th = as.numeric(input)
                phase = phase %>%
                  mutate(th = ifelse(!is.na(Th), Th, th))
                Th = unique(phase$th)
                
                sub_TITLE = paste('Th = ', Th)
                plot(d, main = TITLE, sub = sub_TITLE)
                abline(v = Th, col = "purple")
                
                
                
              }
              
              
              phase %>% mutate(
                reads_s0 = reads_s0 - th,
                reads_s0 = ifelse(reads <= 0 |
                                      reads_s0 <= 0 , 0, reads_s0),
                reads = reads_s0
              ) %>%
                ungroup()
              
            }
  }
  
  dir = '/Volumes/Storage2/RT_DATA/New_DATA/Run_181210/50kb_adj_15M/'
  
  HU = foreach(
    i = 1,
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach')
  ) %do% {
    Files = list.files(dir[i])
    Files = Files[grepl(x = Files, pattern = '^HU')]
    
    foreach(File = Files,
            .combine = 'rbind',
            .packages = 'tidyverse') %do% {
              phase = read_delim(
                paste0(dir[i], File),
                skip = 1,
                col_names = c('chr', 'start', 'end', 'reads'),
                delim = ' '
              ) %>%
                mutate(
                  Condition = paste0('HU', i),
                  phase = str_extract(File, pattern = '[SG][1-4]'),
                  rep = i
                ) %>%
                inner_join(s0) %>%
                mutate(reads_s0 = reads - S0) %>%
                group_by(Condition, rep, phase) %>%
                mutate(reads = ifelse(reads > quantile(reads, .9999)[[1]], 0, reads),
                       th = th_first_peak(reads, S0))
              
              d = density(phase$reads_s0[phase$reads > 0])
              Th = unique(phase$th)
              
              TITLE = paste('HU', i, str_extract(File, pattern = '[SG][1-4]'))
              sub_TITLE = paste('Th = ', Th)
              plot(d, main = TITLE, sub = sub_TITLE)
              abline(v = Th, col = "purple")
              print(Th)
              
              input = NA
              while (input != '' | is.na(input)) {
                input = readline()
                Th = as.numeric(input)
                phase = phase %>%
                  mutate(th = ifelse(!is.na(Th), Th, th))
                Th = unique(phase$th)
                sub_TITLE = paste('Th = ', Th)
                plot(d, main = TITLE, sub = sub_TITLE)
                abline(v = Th, col = "purple")
                
                
              }
              
              
              phase %>% mutate(
                reads_s0 = reads_s0 - th,
                reads_s0 = ifelse(reads <= 0 |
                                      reads_s0 <= 0 , 0, reads_s0),
                reads = reads_s0
              ) %>%
                ungroup()
              
            }
  }
  
  ATRiHU = foreach(
      i = 1,
      .combine = 'rbind',
      .packages = c('tidyverse', 'foreach')
  ) %do% {
      Files = list.files(dir[i])
      Files = Files[grepl(x = Files, pattern = 'InhATR_HU')]
      
      foreach(File = Files,
              .combine = 'rbind',
              .packages = 'tidyverse') %do% {
                  phase = read_delim(
                      paste0(dir[i], File),
                      skip = 1,
                      col_names = c('chr', 'start', 'end', 'reads'),
                      delim = ' '
                  ) %>%
                      mutate(
                          Condition = paste0('ATRiHU', i),
                          phase = str_extract(File, pattern = '[SG][1-4]'),
                          rep = i
                      ) %>%
                      inner_join(s0) %>%
                      mutate(reads_s0 = reads - S0) %>%
                      group_by(Condition, rep, phase) %>%
                      mutate(reads = ifelse(reads > quantile(reads, .9999)[[1]], 0, reads),
                             th = th_first_peak(reads, S0))
                  
                  d = density(phase$reads_s0[phase$reads > 0])
                  Th = unique(phase$th)
                  
                  TITLE = paste('ATRiHU', i, str_extract(File, pattern = '[SG][1-4]'))
                  sub_TITLE = paste('Th = ', Th)
                  plot(d, main = TITLE, sub = sub_TITLE)
                  abline(v = Th, col = "purple")
                  print(Th)
                  
                  input = NA
                  while (input != '' | is.na(input)) {
                      input = readline()
                      Th = as.numeric(input)
                      phase = phase %>%
                          mutate(th = ifelse(!is.na(Th), Th, th))
                      Th = unique(phase$th)
                      sub_TITLE = paste('Th = ', Th)
                      plot(d, main = TITLE, sub = sub_TITLE)
                      abline(v = Th, col = "purple")
                      
                      
                  }
                  
                  
                  phase %>% mutate(
                      reads_s0 = reads_s0 - th,
                      reads_s0 = ifelse(reads <= 0 |
                                            reads_s0 <= 0 , 0, reads_s0),
                      reads = reads_s0
                  ) %>%
                      ungroup()
                  
              }
  }
  
  dir = '/Volumes/Storage2/RT_DATA/New_DATA/Run_190402/50kb_adj_15M/'
  
  
  ATRi = foreach(
      i = 1,
      .combine = 'rbind',
      .packages = c('tidyverse', 'foreach')
  ) %do% {
      Files = list.files(dir[i])
      Files = Files[grepl(x = Files, pattern = 'IA')]
      
      foreach(File = Files,
              .combine = 'rbind',
              .packages = 'tidyverse') %do% {
                  phase = read_delim(
                      paste0(dir[i], File),
                      skip = 1,
                      col_names = c('chr', 'start', 'end', 'reads'),
                      delim = ' '
                  ) %>%
                      mutate(
                          Condition = paste0('ATRi', i),
                          phase = str_extract(File, pattern = '[SG][1-4]'),
                          rep = i
                      ) %>%
                      inner_join(s0) %>%
                      mutate(reads_s0 = reads - S0) %>%
                      group_by(Condition, rep, phase) %>%
                      mutate(reads = ifelse(reads > quantile(reads, .9999)[[1]], 0, reads),
                             th = th_first_peak(reads, S0))
                  
                  d = density(phase$reads_s0[phase$reads > 0])
                  Th = unique(phase$th)
                  
                  TITLE = paste('ATRi', i, str_extract(File, pattern = '[SG][1-4]'))
                  sub_TITLE = paste('Th = ', Th)
                  plot(d, main = TITLE, sub = sub_TITLE)
                  abline(v = Th, col = "purple")
                  print(Th)
                  
                  input = NA
                  while (input != '' | is.na(input)) {
                      input = readline()
                      Th = as.numeric(input)
                      phase = phase %>%
                          mutate(th = ifelse(!is.na(Th), Th, th))
                      Th = unique(phase$th)
                      sub_TITLE = paste('Th = ', Th)
                      plot(d, main = TITLE, sub = sub_TITLE)
                      abline(v = Th, col = "purple")
                      
                      
                  }
                  
                  
                  phase %>% mutate(
                      reads_s0 = reads_s0 - th,
                      reads_s0 = ifelse(reads <= 0 |
                                            reads_s0 <= 0 , 0, reads_s0),
                      reads = reads_s0
                  ) %>%
                      ungroup()
                  
              }
  }
  all_phase = bind_rows(NT, Aph, AphRO, HU,ATRi,ATRiHU)

  
  all_phase %>%
    group_by(Condition, phase)%>%
    dplyr::select(chr, start, end, reads, Condition, phase, th, rep) %>%
    write_tsv('All_normalised_Tracks_50adj_ATRIadd.tsv')


