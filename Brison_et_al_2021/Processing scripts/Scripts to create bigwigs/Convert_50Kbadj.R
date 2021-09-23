library(tidyverse)
folder='~/Desktop/RO_paper_bw_tracks/'
sub_f='50Kb_adj/'

system(paste0('mkdir ',folder,' ',folder,sub_f))

all_phase = read_tsv('All_normalised_Tracks_50adj_normalised_with_sim_periodic_1kbXY_ATRiadd.tsv')
URI = read_tsv('URi_simulation_periodic_1kbXY_ATRIadd.tsv')
S50 = read_tsv('S50_simulation_periodic_1kbXY_ATRiadd.tsv')


for (condition in unique(all_phase$Condition)){
    
    for (Phase in unique(all_phase$phase)){
        
        
        all_phase%>%
            filter(phase==Phase,
                   Condition==condition)%>%
            dplyr::select(chr,start,end,reads)%>%
            arrange(chr,start)%>%
            write_tsv(paste0(folder,sub_f,condition,'_',str_replace_all(string = Phase,pattern = '/',replacement = '-'),'.bedGraph'),col_names = F)
        
    }
    
}
for (condition in unique(URI$Condition)){
        
    URI%>%
            filter( Condition==condition)%>%
            dplyr::select(chr,start,end,URI)%>%
            arrange(chr,start)%>%
            drop_na()%>%
            write_tsv(paste0(folder,sub_f,'URI_',condition,'.bedGraph'),col_names = F)
        
}
for (condition in unique(S50$Condition)){
    
    S50%>%
        filter( Condition==condition)%>%
        dplyr::select(chr,start,end,s50)%>%
        arrange(chr,start)%>%
        write_tsv(paste0(folder,sub_f,'S50_',condition,'.bedGraph'),col_names = F)
    
}
#average

tracks = all_phase %>% ungroup() %>%
    mutate(phase = factor(phase, levels = c(
        'G1/S1', 'S2', 'S3', 'S4', 'S5', 'S6/G2/M'
    ))) %>%
    mutate(Condition = str_remove(Condition, '[1-9]{1,2}$')) %>%
    group_by(Condition, chr, start, end, phase) %>%
    summarise(reads = mean(reads,na.rm = T))%>%
    group_by(Condition,chr,phase)%>%
    ungroup() %>%
    drop_na()

for (condition in unique(tracks$Condition)){
    
    for (Phase in unique(all_phase$phase)){
        
        
        tracks%>%
            filter(phase==Phase,
                   Condition==condition)%>%
            dplyr::select(chr,start,end,reads)%>%
            arrange(chr,start)%>%
            write_tsv(paste0(folder,sub_f,condition,'_Average_',str_replace_all(string = Phase,pattern = '/',replacement = '-'),'.bedGraph'),col_names = F)
        
    }
    
}

#files renaming 

system(paste0("rename 's/AphRO/ARO/g' ",folder,sub_f,"*"))
system(paste0("rename 's/ATRiHU/HU+ATRi/g' ",folder,sub_f,"*"))


#convert to bw
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
paste0('mv ',folder,sub_f,'*_Average_* ',folder,sub_f,'Average' )
)
system(
paste0('mv ',folder,sub_f,'*.bigWig ',folder,sub_f,'Single' )
)
