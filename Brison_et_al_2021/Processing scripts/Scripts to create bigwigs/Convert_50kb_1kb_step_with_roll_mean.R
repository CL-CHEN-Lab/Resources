folder='~/Desktop/RO_paper_bw_tracks/'
sub_f='50Kb_1kb_step_smoothed_with_roll_mean/'

system(paste0('mkdir ',folder,' ',folder,sub_f))
#Normalize with simulation data
#produce files for further analysis
library(tidyverse)
library(zoo)
library(foreach)
library(doSNOW)

cl=makeCluster(8)
registerDoSNOW(cl)

#load FACS profile

facs=read_csv('FACS/AphRO.csv',col_names = F)

gates=tibble(
    phase=c('G1_all','G1/S1','S2','S3','S4','S5','S6/G2/M'),
    start=c(21.6,24.2519169380677,27.282955036018755,30.61406625807228,34.19265411124849,37.47504529876629,40.60739992423285),
    end=c(26.4,27.282955036018755,30.61406625807228,34.19265411124849,37.47504529876629,40.60739992423285,47.27059148295092),
    G1_peak = 24.15,
    G2_peak = 44.2)

gates_to_plot=gates%>%filter(phase!='G1_all')%>%mutate(m=(start+end)/2)%>%
    mutate(percent_start=100*(start-G1_peak)/(G2_peak-G1_peak),
           percent_end=100*(end-G1_peak)/(G2_peak-G1_peak))

#load simylation results
results=read_tsv('Results_periodic_simulationXY_1kb.tsv')

# calculate normalization factor
results=results%>%
    group_by(percent)%>%
    summarise(Forks=mean(Forks))

gates=gates%>%
    mutate(start=round(100*(((start-G1_peak)/(G2_peak-G1_peak))),2),
           end=round(100*(((end-G1_peak)/(G2_peak-G1_peak))),2))

results$phase=NA

for (i in 2:length(gates$phase)){
    
    results=results%>%
        mutate(
            phase=ifelse(percent > gates$start[i] & percent < gates$end[i],gates$phase[i],phase)
        )
}

results=results%>%drop_na()%>%
    group_by(phase)%>%
    summarise(Forks=mean(Forks))%>%
    ungroup()%>%
    mutate(Normalize_factor=90000000*Forks/sum(Forks))

to_norm = read_tsv('All_normalised_Tracks_50adj_normalised_with_sim_periodic_1kbXY.tsv') %>% group_by(Condition, phase, rep, Normalize_factor) %>% summarise(th = min(th))
#load 50kb smoothed tracks and average them
s0 = read_delim(
    '/Volumes/Storage/RT_DATA/Run_150521/1kb/UnsortedNTrep1_delDupl.bam.bedgraph',
    skip = 1,
    col_names = c('chr', 'start', 'end', 'S0'),
    delim = ' '
) %>% 
    group_by(chr) %>% 
    mutate(S0 = rollmean(S0, k = 50, fill = NA)) %>%
    drop_na() %>% ungroup() %>%
    mutate(S0 = (15 * 10 ^ 6) * S0 / sum(S0))


#load files
dir = c(
    '/Volumes/Storage/RT_DATA/Run_160308//1kb/',
    '/Volumes/Storage/RT_DATA/Run_200214/1kb/'
)


AphRO = foreach(
    i = 1:2,
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach','zoo')
) %do% {
    Files = list.files(dir[i])
    Files = Files[grep(x = Files, pattern = 'Rep5-AphiRO16h|aph-RO-16h')]
    Files = Files[grep(x = Files, pattern = '.bg$|.bedgraph$')]
    
    foreach(File = Files,
            .combine = 'rbind',
            .packages = c('tidyverse','zoo')) %dopar% {
                
                C = paste0('AphRO', i)
                P= str_extract(File, pattern = '[SG][1-4]')
              read_delim(
                    paste0(dir[i], File),
                    col_names = c('chr', 'start', 'end', 'reads'),
                    skip = 1,
                    delim = ' '
                    )%>% group_by(chr) %>% 
                  mutate(reads = rollmean(reads, k = 50, fill = NA)) %>%
                  drop_na() %>% ungroup() %>%
                  mutate(reads = (15 * 10 ^ 6) * reads / sum(reads)) %>%
                    mutate(
                        Condition = C,
                        phase = P,
                        rep = i,
                        phase=case_when(
                            phase=='G1'~'G1/S1',
                            phase=='S1'~'S2',
                            phase=='S2'~'S3',
                            phase=='S3'~'S4',
                            phase=='S4'~'S5',
                            phase=='G2'~'S6/G2/M'

                        )
                    ) %>%
                    inner_join(to_norm) %>%
                    inner_join(s0) %>%
                    group_by(Condition, rep, phase) %>%
                    mutate(
                        reads_s0 = reads - S0 - th,
                        reads = ifelse(
                            reads <= 0 |
                                reads_s0 <= 0 |
                                reads > quantile(reads, .9999)[[1]],
                            0,
                            reads_s0
                        ),
                        reads = Normalize_factor* reads / sum(reads)

                    ) %>%
                    ungroup() %>%
                    dplyr::select(-reads_s0)%>%
                   dplyr::select(chr,start,end,reads)%>%
                   arrange(chr,start)%>%
                   write_tsv(paste0(folder,sub_f,C,'_',str_replace_all(string = P,pattern = '/',replacement = '-'),'.bedGraph'),col_names = F)

              tibble( Condition=C,
                      Phase=P,
                      Directory= paste0(folder,sub_f,C,'_',str_replace_all(string = P,pattern = '/',replacement = '-'),'.bedGraph'))
            }
}


dir = c(
    '/Volumes/Storage/RT_DATA/Run_141216/1kb/',
    '/Volumes/Storage/RT_DATA/Run_150724/1kb/'
)

Aph = foreach(
    i = 1:2,
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach','zoo')
) %do% {
    Files = list.files(dir[i])
    Files = Files[grepl(x = Files, pattern = 'Aph16hrep|Aph-600nM-16h|Aph16hDT4')]
    Files = Files[grep(x = Files, pattern = '.bg$|.bedgraph$')]
    
   
    foreach(File = Files,
            .combine = 'rbind',
            .packages = c('tidyverse','zoo')) %dopar% {
                
                C = paste0('Aph', i)
                P= str_extract(File, pattern = '[SG][1-4]')
                
                read_delim(
                    paste0(dir[i], File),
                    col_names = c('chr', 'start', 'end', 'reads'),
                    skip = 1,
                    delim = ' '
                ) %>%
                    group_by(chr) %>%
                    mutate(reads = rollmean(reads, k = 50, fill = NA)) %>%
                    drop_na() %>% ungroup() %>%
                    mutate(reads = (15 * 10 ^ 6) * reads / sum(reads)) %>%
                    mutate(
                        Condition = C,
                        phase = P,
                        rep = i,
                        phase=case_when(
                            phase=='G1'~'G1/S1',
                            phase=='S1'~'S2',
                            phase=='S2'~'S3',
                            phase=='S3'~'S4',
                            phase=='S4'~'S5',
                            phase=='G2'~'S6/G2/M'

                        )
                    ) %>%
                    inner_join(to_norm) %>%
                    inner_join(s0) %>%
                    group_by(Condition, rep, phase) %>%
                    mutate(
                        reads_s0 = reads - S0 - th,
                        reads = ifelse(
                            reads <= 0 |
                                reads_s0 <= 0 |
                                reads > quantile(reads, .9999)[[1]],
                            0,
                            reads_s0
                        ),
                        reads = Normalize_factor* reads / sum(reads)

                    ) %>%
                    ungroup() %>%
                    dplyr::select(-reads_s0)%>%
                    dplyr::select(chr,start,end,reads)%>%
                    arrange(chr,start)%>%
                    write_tsv(paste0(folder,sub_f,C,'_',str_replace_all(string = P,pattern = '/',replacement = '-'),'.bedGraph'),col_names = F)
                
                tibble( Condition=C,
                        Phase=P,
                        Directory= paste0(folder,sub_f,C,'_',str_replace_all(string = P,pattern = '/',replacement = '-'),'.bedGraph'))
            }
}



dir = c(
    '/Volumes/Storage/RT_DATA/Run_170125/1kb/',
    '/Volumes/Storage/RT_DATA/Run_150521/1kb/',
    '/Volumes/Storage/RT_DATA/Run_180719/1kb/'
)

NT = foreach(
    i = 1:3,
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach','zoo')
) %do% {
    Files = list.files(dir[i])
    Files = Files[grepl(x = Files, pattern = 'NT.*[GS]')]
    Files = Files[!grepl(x = Files, pattern = 'Aph2h')]
    Files = Files[grep(x = Files, pattern = '.bg$|.bedgraph$')]
    
    
    foreach(File = Files,
            .combine = 'rbind',
            .packages = c('tidyverse','zoo')) %dopar% {
                
                C = paste0('NT', i)
                P= str_extract(File, pattern = '[SG][1-4]')
                
                read_delim(
                    paste0(dir[i], File),
                    col_names = c('chr', 'start', 'end', 'reads'),
                    skip = 1,
                    delim = ' '
                    ) %>%
                    group_by(chr) %>% 
                    mutate(reads = rollmean(reads, k = 50, fill = NA)) %>%
                    drop_na() %>% ungroup() %>%
                    mutate(reads = (15 * 10 ^ 6) * reads / sum(reads)) %>%
                    mutate(
                        Condition =C,
                        phase = P,
                        rep = i,
                        phase=case_when(
                            phase=='G1'~'G1/S1',
                            phase=='S1'~'S2',
                            phase=='S2'~'S3',
                            phase=='S3'~'S4',
                            phase=='S4'~'S5',
                            phase=='G2'~'S6/G2/M'

                        )
                    ) %>%
                    inner_join(to_norm) %>%
                    inner_join(s0) %>%
                    group_by(Condition, rep, phase) %>%
                    mutate(
                        reads_s0 = reads - S0 - th,
                        reads = ifelse(
                            reads <= 0 |
                                reads_s0 <= 0 |
                                reads > quantile(reads, .9999)[[1]],
                            0,
                            reads_s0
                        ),
                        reads = Normalize_factor* reads / sum(reads)
                    ) %>%
                    ungroup() %>%
                    dplyr::select(-reads_s0)%>%
                    dplyr::select(chr,start,end,reads)%>%
                    arrange(chr,start)%>%
                    write_tsv(paste0(folder,sub_f,C,'_',str_replace_all(string = P,pattern = '/',replacement = '-'),'.bedGraph'),col_names = F)
                
                tibble( Condition=C,
                        Phase=P,
                        Directory= paste0(folder,sub_f,C,'_',str_replace_all(string = P,pattern = '/',replacement = '-'),'.bedGraph'))
            }
}


dir = '/Volumes/Storage/RT_DATA/Run_181210/1kb/'

HU = foreach(
    i = 1,
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach','zoo')
) %dopar% {
    Files = list.files(dir[i])
    Files = Files[grepl(x = Files, pattern = '^HU')]
    Files = Files[grep(x = Files, pattern = '.bg$|.bedgraph$')]
    
    foreach(File = Files,
            .combine = 'rbind',
            .packages = c('tidyverse','zoo')) %dopar% {
                C = paste0('HU', i)
                P= str_extract(File, pattern = '[SG][1-4]')

                read_delim(
                    paste0(dir[i], File),
                    col_names = c('chr', 'start', 'end', 'reads'),
                    skip = 1,
                    delim = ' '
                    ) %>%
                    group_by(chr) %>%
                    mutate(reads = rollmean(reads, k = 50, fill = NA)) %>%
                    drop_na() %>% ungroup() %>%
                    mutate(reads = (15 * 10 ^ 6) * reads / sum(reads)) %>%
                    mutate(
                        Condition = C,
                        phase = P,
                        rep = i,
                        phase=case_when(
                            phase=='G1'~'G1/S1',
                            phase=='S1'~'S2',
                            phase=='S2'~'S3',
                            phase=='S3'~'S4',
                            phase=='S4'~'S5',
                            phase=='G2'~'S6/G2/M'

                        )
                    ) %>%
                    inner_join(to_norm) %>%
                    inner_join(s0) %>%
                    group_by(Condition, rep, phase) %>%
                    mutate(
                        reads_s0 = reads - S0 - th,
                        reads = ifelse(
                            reads <= 0 |
                                reads_s0 <= 0 |
                                reads > quantile(reads, .9999)[[1]],
                            0,
                            reads_s0
                        ),
                        reads = Normalize_factor* reads / sum(reads)
                    ) %>%
                    ungroup() %>%
                    dplyr::select(-reads_s0)%>%
                    dplyr::select(chr,start,end,reads)%>%
                    arrange(chr,start)%>%
                    write_tsv(paste0(folder,sub_f,C,'_',str_replace_all(string = P,pattern = '/',replacement = '-'),'.bedGraph'),col_names = F)
                
                tibble( Condition=C,
                        Phase=P,
                        Directory= paste0(folder,sub_f,C,'_',str_replace_all(string = P,pattern = '/',replacement = '-'),'.bedGraph'))
            }
}

ATRiHU = foreach(
    i = 1,
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach','zoo')
) %do% {
    Files = list.files(dir[i])
    Files = Files[grepl(x = Files, pattern = 'InhATR_HU')]
    Files = Files[grep(x = Files, pattern = '.bg$|.bedgraph$')]
    
    foreach(File = Files,
            .combine = 'rbind',
            .packages = c('tidyverse','zoo')) %dopar% {
                
                C = paste0('ATRiHU', i)
                P= str_extract(File, pattern = '[SG][1-4]')
                read_delim(
                    paste0(dir[i], File),
                    col_names = c('chr', 'start', 'end', 'reads'),
                    skip = 1,
                    delim = ' ' 
                    ) %>%
                    group_by(chr) %>% 
                    mutate(reads = rollmean(reads, k = 50, fill = NA)) %>%
                    drop_na() %>% ungroup() %>%
                    mutate(reads = (15 * 10 ^ 6) * reads / sum(reads)) %>%
                    mutate(
                        Condition = C,
                        phase =P,
                        rep = i,
                        phase=case_when(
                            phase=='G1'~'G1/S1',
                            phase=='S1'~'S2',
                            phase=='S2'~'S3',
                            phase=='S3'~'S4',
                            phase=='S4'~'S5',
                            phase=='G2'~'S6/G2/M'

                        )
                    ) %>%
                    inner_join(to_norm) %>%
                    inner_join(s0) %>%
                    group_by(Condition, rep, phase) %>%
                    mutate(
                        reads_s0 = reads - S0 - th,
                        reads = ifelse(
                            reads <= 0 |
                                reads_s0 <= 0 |
                                reads > quantile(reads, .9999)[[1]],
                            0,
                            reads_s0
                        ),
                        reads = Normalize_factor* reads / sum(reads)
                    ) %>%
                    ungroup() %>%
                    dplyr::select(-reads_s0)%>%
                    dplyr::select(chr,start,end,reads)%>%
                    arrange(chr,start)%>%
                    write_tsv(paste0(folder,sub_f,C,'_',str_replace_all(string = P,pattern = '/',replacement = '-'),'.bedGraph'),col_names = F)
                
                tibble( Condition=C,
                        Phase=P,
                        Directory= paste0(folder,sub_f,C,'_',str_replace_all(string = P,pattern = '/',replacement = '-'),'.bedGraph'))
            }
}


all = rbind(
    Aph,
    AphRO,
    NT,
    HU ,
    ATRiHU
)


all=all%>%
    mutate(group_cond=str_remove(string = Condition,pattern = '[1-9]$'))%>%
    group_by(group_cond)%>%
    mutate(n=length(unique(Condition)))


data=foreach(CG=unique(all$group_cond),.combine = 'full_join')%do%{
    sub_all=all%>%
        filter(group_cond==CG)
    Line_C=foreach(C=unique(sub_all$Condition),.combine = 'rbind')%do%{
        sub_sub_all=sub_all%>%
            filter(Condition==C)
    Line=foreach(D=unique(sub_sub_all$Directory),.combine = 'rbind')%dopar%{
        
        read_tsv(D,col_names = F)
        
    }

    Line%>%
        group_by(X1,X2,X3)%>%
        summarise(X4=sum(X4,na.rm = T))%>%
        ungroup()

    }
    
    Line_C%>%
        group_by(X1,X2,X3)%>%
        summarise(X4=mean(X4,na.rm = T))%>%
        ungroup()%>%
        `colnames<-`(c('chr','start','end',CG))
    
}

Columns=names(data)

Columns=Columns[!Columns%in%c('start','end','chr','NT')]

#URI
foreach(C=Columns)%do%{
    
    sub=data%>%select_('chr','start','end','NT',C)
    sub=sub[sub[,C] >0 & sub[,'NT'] > 0 & !is.na(sub[,C]) & !is.na(sub[,'NT']),]
    sub$delta=pull(sub[,C])-sub$NT
    sub$URI=(sub$delta-mean(sub$delta,na.rm = T))/sd(sub$delta,na.rm = T)
    sub%>%
         dplyr::select(chr,start,end,URI)%>%
         arrange(chr,start)%>%
            drop_na()%>%
         write_tsv(paste0(folder,sub_f,'URI_',C,'.bedGraph'),col_names = F) 
    
    C
}


all_average = foreach(CG = unique(all$group_cond),.combine = 'rbind',.packages = 'tidyverse') %do% {
    sub_all = all %>%
        filter(group_cond == CG)
    foreach(P = unique(sub_all$Phase),.combine = 'rbind',.packages = 'tidyverse') %dopar% {
        sub_phase_all = sub_all %>%
            filter(Phase == P)
        Line = foreach(D = sub_phase_all$Directory, .combine = 'rbind',.packages = 'tidyverse') %do%
        {
            read_tsv(D, col_names = F)

        }

        Line %>%
            group_by(X1,X2,X3)%>%
            summarise(X4 = mean(X4,na.rm = T)) %>%
            dplyr::select(X1, X2, X3, X4) %>%
            arrange(X1, X2) %>%
            drop_na()%>%
            write_tsv(
                paste0(
                    folder,sub_f,
                    'Average_profile_',
                    CG,
                    '_',
                    P,
                    '.bedGraph'
                ),
                col_names = F
            )
        
    
        tibble(Condition=CG,Phase=P,Directory= paste0(
            folder,sub_f,'Average_profile_',
            CG,
            '_',
            P,
            '.bedGraph'
        ))
    }
    

    
}


data = foreach(CG = unique(all_average$Condition),.packages = c('tidyverse','foreach')) %do% {
    sub_all = all_average %>%
        filter(Condition == CG)
    
    all_phases=foreach(P=sub_all$Phase,.combine = 'full_join',.packages = c('tidyverse','foreach'))%do%{
     
        read_tsv(sub_all$Directory[sub_all$Phase==P], col_names =c('chr','start','end',P))
           
    }
    
    chunks=trunc(seq(1,length(all_phases$chr)+1,length.out =  300))
    
    S50=foreach(i=1:299,.combine = 'rbind',.packages = c('tidyverse','foreach'))%dopar%{
        
        Lines=all_phases[chunks[i]:(chunks[i+1]-1),]
        x=Lines%>%
            select(G1,S1,S2,S3,S4,G2)%>%
            as.matrix()
        x[is.na(x)] <- 0

       x = t(apply(x, 1, FUN = function(x) rep(x,each=1000)))
        s50=apply(x,1,FUN = function(x) min(which(cumsum(x) >= sum(x)/2))/length(x))
        
        Lines[,c('chr','start','end')]%>%
            mutate(S50=s50)
        
    }
    
    S50 %>%
        arrange(chr, start) %>%
        write_tsv(
            paste0(
                folder,sub_f,'S50_',
                CG,
                '.bedGraph'
            ),
            col_names = F
        )
    
    
}

stopCluster(cl)


#files renaming 

system(paste0("rename 's/AphRO/ARO/g' ",folder,sub_f,"*"))
system(paste0("rename 's/ATRiHU/HU+ATRi/g' ",folder,sub_f,"*"))

cov=tibble(
    old=c('_G1.bedGraph','_G2.bedGraph','_S4.bedGraph','_S3.bedGraph','_S2.bedGraph','_S1.bedGraph'),
    new=c('_G1-S1.bedGraph','_S6-G2-M.bedGraph','_S5.bedGraph','_S4.bedGraph','_S3.bedGraph','_S2.bedGraph')
)

for (i in 1:length(cov$old)){
    
    
    system(  paste0('rename -s ',cov$old[i],' ',cov$new[i],' ',folder,sub_f,'*') )
    
    
}

#create bw
 system(
     paste0('for i in ', folder,sub_f,'*.bedGraph; do /Users/sgnan/anaconda3/bin/bedGraphToBigWig $i ~/H19_chr_size.txt ${i/%.bedGraph/.bigWig}; rm $i; done')
 )


 system(
     paste0('mkdir ',folder,sub_f,'S50 ',folder,sub_f,'URI ',folder,sub_f,'Average ',folder,sub_f,'Single')
 )
 system(
     paste0('mv ',folder,sub_f,'URI_* ',folder,sub_f,'URI' )
 )
 system(
     paste0('mv ',folder,sub_f,'S50_* ',folder,sub_f,'S50' )
 )
 system(
     paste0('mv ',folder,sub_f,'Average_profile_* ',folder,sub_f,'Average' )
 )
 system(
     paste0('mv ',folder,sub_f,'*.bigWig ',folder,sub_f,'Single' )
 )

