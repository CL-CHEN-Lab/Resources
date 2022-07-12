# RO analysis subtraction the S0
library(tidyverse)
library(foreach)
library(doSNOW)
library(mixtools)

dir.create('noise_remove')

#Run 150521
s0 = read_delim(
    '/Volumes/Storage/RT_DATA/Run_150521/50kb_adj_15M/UnsortedNTrep1_delDupl.bam.bedgraph',
    skip = 1,
    col_names = c('chr', 'start', 'end', 'S0'),
    delim = ' '
)

smoothing = function(x,
                     bins = 10, sigma = 15) {
    n = round((max(x) - min(x)) / 10 ^ -2)
    
    hist = hist(x, n, plot = F)
    
    z = -(sigma * bins):(sigma * bins)
    convolve_y <- (1 / sqrt(2 * pi * sigma ^ 2)) * exp(-(z ^ 2 / (2 * sigma ^
                                                                      2)))
    Gaussian = c(
        rep(NA, bins * sigma),
        convolve(hist$counts, convolve_y, type = 'filter'),
        rep(NA, bins * sigma)
    )
    
    return(list(x = hist$mids, y = Gaussian))
    
}

th_first_peak = function(x,
                         cont,
                         bins = 10,
                         sigma = 15) {
    cont = cont[x > 0]
    x = x[x > 0]
    x = x - cont
    
    n = round((max(x) - min(x)) / 10 ^ -2)
    
    x = hist(x, n, plot = F)
    
    z = -(sigma * bins):(sigma * bins)
    
    convolve_y <-
        (-z / (sigma ^ 3 * sqrt(2 * pi * sigma ^ 2))) * exp(-(z ^ 2 / (2 * sigma ^
                                                                           2)))
    
    x$first_der = c(
        rep(NA, bins * sigma),
        convolve(x$counts, convolve_y, type = 'filter'),
        rep(NA, bins * sigma)
    )
    x$lagged_first_der = lag(x$first_der)
    
    convolve_y <-
        ((sigma ^ 2 - z ^ 2) / sigma ^ 5) * (1 / sqrt(2 * pi)) * exp(-(z ^ 2 /
                                                                           (2 * sigma ^ 2)))
    x$second_der = c(
        rep(NA, bins * sigma),
        convolve(x$counts, convolve_y, type = 'filter'),
        rep(NA, bins * sigma)
    )
    x$lagged_second_der = lag(x$second_der)
    
    
    x = tibble(
        intensity = x$counts,
        th = x$mids,
        maxima = case_when(x$first_der > 0 &
                               x$lagged_first_der < 0 ~ 'Max',
                           T ~ 'Rem'),
        flex = case_when(
            x$second_der > 0 & x$lagged_second_der < 0 ~ 'start',
            x$second_der < 0 &
                x$lagged_second_der > 0 ~ 'stop',
            T ~ 'Rem'
        )
    )
    
    local_maxima = x %>% filter(maxima == 'Max', intensity > 0.01 * max(intensity)) %>% pull(th)
    
    flexes = x %>%
        filter(flex != 'Rem') %>%
        mutate(flex_located_before = lag(flex)) %>%
        filter(flex != flex_located_before |
                   is.na(flex_located_before)) %>%
        filter(!(row_number() == 1 &
                     flex == 'stop'),
               !(row_number() == n() & flex == 'start'))
    
    if (nrow(flexes) != 0) {
        flexes = flexes %>%
            mutate(peak_id = (rep(1:(n(
            ) / 2), each = 2))) %>%
            dplyr::select(th, peak_id, flex) %>%
            tidyr::spread(key = flex, value = th)
        
        results = sapply(local_maxima, function(y) {
            result = (flexes$start <= y & flexes$stop >= y)
            result = ifelse(sum(result) == 1, y, NA)
            return(result)
        })
        results = results[!is.na(results)]
        
        if (length(results) == 0) {
            return(NULL)
        } else{
            return(results[1])
            
        }
    } else{
        return(NULL)
    }
    
}



#load files

dir = c(
    '/Volumes/Storage/RT_DATA/Run_160308//50kb_adj_15M/',
    '/Volumes/Storage/RT_DATA/Run_200214/50kb_adj_15M/'
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
                x = phase$reads_s0[phase$reads > 0]
                d = smoothing(x)
                Th = unique(phase$th)
                
                TITLE = paste('AphRO', i, str_extract(File, pattern = '[SG][1-4]'))
                sub_TITLE = paste('Th = ', Th)
                plot(
                    d$x,
                    d$y,
                    main = paste(TITLE, File),
                    sub = sub_TITLE,
                    type = 'l',
                    xlab = 'Normalized reads',
                    ylab = 'counts per bin'
                )
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
                    jpeg(filename = file.path('noise_remove', paste0(TITLE, '.jpeg')))
                    
                    plot(
                        d$x,
                        d$y,
                        main = paste(TITLE, File),
                        sub = sub_TITLE,
                        type = 'l',
                        xlab = 'Normalized reads',
                        ylab = 'counts per bin'
                    )
                    abline(v = Th, col = "purple")
                    
                    dev.off()
                    
                    plot(
                        d$x,
                        d$y,
                        main = paste(TITLE, File),
                        sub = sub_TITLE,
                        type = 'l',
                        xlab = 'Normalized reads',
                        ylab = 'counts per bin'
                    )
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
    '/Volumes/Storage/RT_DATA/Run_141216/50kb_adj_15M/',
    '/Volumes/Storage/RT_DATA/Run_150724/50kb_adj_15M/'
)


Aph = foreach(
    i = 1:2,
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach')
) %do% {
    Files = list.files(dir[i])
    Files = Files[grepl(x = Files, pattern = 'Aph16hrep|Aph-600nM-16h|Aph16hDT4|Rep5-Aphi16h-DT40')]
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
                
                x = phase$reads_s0[phase$reads > 0]
                d = smoothing(x)
                Th = unique(phase$th)
                
                TITLE = paste('Aph', i, str_extract(File, pattern = '[SG][1-4]'))
                sub_TITLE = paste('Th = ', Th)
                plot(
                    d$x,
                    d$y,
                    main = paste(TITLE, File),
                    sub = sub_TITLE,
                    type = 'l',
                    xlab = 'Normalized reads',
                    ylab = 'counts per bin'
                )
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
                    jpeg(filename = file.path('noise_remove', paste0(TITLE, '.jpeg')))
                    
                    plot(
                        d$x,
                        d$y,
                        main = paste(TITLE, File),
                        sub = sub_TITLE,
                        type = 'l',
                        xlab = 'Normalized reads',
                        ylab = 'counts per bin'
                    )
                    abline(v = Th, col = "purple")
                    
                    dev.off()
                    
                    plot(
                        d$x,
                        d$y,
                        main = paste(TITLE, File),
                        sub = sub_TITLE,
                        type = 'l',
                        xlab = 'Normalized reads',
                        ylab = 'counts per bin'
                    )
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
    '/Volumes/Storage/RT_DATA/Run_170125/50kb_adj_15M/',
    '/Volumes/Storage/RT_DATA/Run_150521/50kb_adj_15M/',
    '/Volumes/Storage/RT_DATA/Run_180719/50kb_adj_15M/'
)
NT = foreach(
    i = 1:3,
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
                
                x = phase$reads_s0[phase$reads > 0]
                d = smoothing(x)
                Th = unique(phase$th)
                
                TITLE = paste('NT', i, str_extract(File, pattern = '[SG][1-4]'))
                sub_TITLE = paste('Th = ', Th)
                plot(
                    d$x,
                    d$y,
                    main = paste(TITLE, File),
                    sub = sub_TITLE,
                    type = 'l',
                    xlab = 'Normalized reads',
                    ylab = 'counts per bin'
                )
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
                    jpeg(filename = file.path('noise_remove', paste0(TITLE, '.jpeg')))
                    
                    plot(
                        d$x,
                        d$y,
                        main = paste(TITLE, File),
                        sub = sub_TITLE,
                        type = 'l',
                        xlab = 'Normalized reads',
                        ylab = 'counts per bin'
                    )
                    abline(v = Th, col = "purple")
                    
                    dev.off()
                    
                    plot(
                        d$x,
                        d$y,
                        main = paste(TITLE, File),
                        sub = sub_TITLE,
                        type = 'l',
                        xlab = 'Normalized reads',
                        ylab = 'counts per bin'
                    )
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

dir = '/Volumes/Storage/RT_DATA/Run_181210/50kb_adj_15M/'

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
                
                x = phase$reads_s0[phase$reads > 0]
                d = smoothing(x)
                Th = unique(phase$th)
                
                TITLE = paste('HU', i, str_extract(File, pattern = '[SG][1-4]'))
                sub_TITLE = paste('Th = ', Th)
                plot(
                    d$x,
                    d$y,
                    main = paste(TITLE, File),
                    sub = sub_TITLE,
                    type = 'l',
                    xlab = 'Normalized reads',
                    ylab = 'counts per bin'
                )
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
                    jpeg(filename = file.path('noise_remove', paste0(TITLE, '.jpeg')))
                    
                    plot(
                        d$x,
                        d$y,
                        main = paste(TITLE, File),
                        sub = sub_TITLE,
                        type = 'l',
                        xlab = 'Normalized reads',
                        ylab = 'counts per bin'
                    )
                    abline(v = Th, col = "purple")
                    
                    dev.off()
                    
                    plot(
                        d$x,
                        d$y,
                        main = paste(TITLE, File),
                        sub = sub_TITLE,
                        type = 'l',
                        xlab = 'Normalized reads',
                        ylab = 'counts per bin'
                    )
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
                
                x = phase$reads_s0[phase$reads > 0]
                d = smoothing(x)
                Th = unique(phase$th)
                
                TITLE = paste('ATRiHU', i, str_extract(File, pattern = '[SG][1-4]'))
                sub_TITLE = paste('Th = ', Th)
                plot(
                    d$x,
                    d$y,
                    main = paste(TITLE, File),
                    sub = sub_TITLE,
                    type = 'l',
                    xlab = 'Normalized reads',
                    ylab = 'counts per bin'
                )
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
                    jpeg(filename = file.path('noise_remove', paste0(TITLE, '.jpeg')))
                    
                    plot(
                        d$x,
                        d$y,
                        main = paste(TITLE, File),
                        sub = sub_TITLE,
                        type = 'l',
                        xlab = 'Normalized reads',
                        ylab = 'counts per bin'
                    )
                    abline(v = Th, col = "purple")
                    
                    dev.off()
                    
                    plot(
                        d$x,
                        d$y,
                        main = paste(TITLE, File),
                        sub = sub_TITLE,
                        type = 'l',
                        xlab = 'Normalized reads',
                        ylab = 'counts per bin'
                    )
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


all_phase = bind_rows(NT, Aph, AphRO, HU, ATRiHU)


all_phase %>%
    dplyr::select(chr, start, end, reads, Condition, phase, th, rep) %>%
    write_tsv('All_normalised_Tracks_50adj.tsv')
