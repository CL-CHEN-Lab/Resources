#Normalize data with simulation
library(tidyverse)
library(zoo)
library(foreach)
library(corrplot)
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

# plot facs profile
facs%>%
    mutate(X1=round(X1,1))%>%
    group_by(X1)%>%
    summarise(X2=min(X2))%>%
    ungroup()%>%
    mutate(X2=rollmean(X2,2,fill = NA))%>%
    ggplot()+geom_line(aes(X1,X2))+
    geom_errorbarh(data = gates_to_plot,aes(xmin=start,xmax =end-0.15,y=270,color=phase , height = 30))+
    geom_text(data = gates_to_plot,aes(x=m,y=270,color=phase, label=phase ),vjust=-0.55)+
    scale_x_continuous(sec.axis = sec_axis(~100*(.-gates$G1_peak)/(gates$G2_peak-gates$G1_peak), breaks =c(0,25,50,75,100),labels = paste0(c(0,25,50,75,100),'%') ,name = '% of replicated DNA'))+xlab('fluorescent units')+ylab('cells counts')+theme_classic()+
    theme(legend.position = 'none',axis.line.x.top = element_blank())+
    geom_segment(x=23,xend=45,y=Inf,yend=Inf)

#load simulation data

# dir='SimulationXY/'
# files=paste0(dir,list.files(dir))
# 
# results=foreach(file=files,.combine = 'rbind',.packages = 'tidyverse')%do%{
#     x=read_delim(file,delim = ' ')
# 
#     if(max(x$percent)> 99.95){
#         x
#     }else{
#         tibble()
#     }
# 
# 
# }
# 
# results%>%write_tsv('Results_periodic_simulationXY_1kb.tsv')

results=read_tsv('Results_periodic_simulationXY_1kb.tsv')

MAX_ratio=results%>%group_by(cycles)%>%summarise(Initiations=mean(Initiations,na.rm = T),unreplicated_dna=mean(unreplicated_dna,na.rm = T))%>%
    mutate(Ratio_firing_unreplicated=Initiations/unreplicated_dna)%>%pull(Ratio_firing_unreplicated)%>%max(na.rm = T)


# firing has been smoothed
results %>%
    group_by(cycles)%>%
    summarise(time=mean(time),
              Initiations=mean(Initiations),
              unreplicated_dna=mean(unreplicated_dna),
              percent=mean(percent)/100)%>%
    mutate(Ratio_firing_unreplicated=Initiations/unreplicated_dna,
           Ratio_firing_unreplicated=rollmean(Ratio_firing_unreplicated,k = 15,fill = 0,align = 'left')
    )%>%
    ggplot() +
    geom_line(aes(
        time/60,
        Ratio_firing_unreplicated / MAX_ratio
    ), color = 'red') +
    geom_line(aes(time/60, percent), color = 'blue') +
    scale_y_continuous(
        name = 'Replicated DNA',
        labels = scales::percent,
        sec.axis = sec_axis(name = expression(paste(
            frac('Origins Firing', 'unreplicated DNA'),
            ' (',
            Mb ^ -1,
            min ^ -1,
            ')'
        )),trans = ~.*MAX_ratio))+
    xlab('time (h)')+theme_classic()+theme(axis.title.y.left  = element_text(color='blue'),
                                           axis.title.y.right  =  element_text(color='red'),
                                           text = element_text(family = 'Courier')
    )


results%>%
    group_by(cycles)%>%
    summarise(percent=mean(percent),
              Forks=mean(Forks))%>%
    ggplot()+
    geom_line(aes(percent/100,Forks/max(Forks)))+
    geom_errorbarh(data = gates_to_plot,aes(xmin=percent_start/100,xmax =(percent_end-0.6)/100,y=0.5,color=phase , height = 0.2))+
    geom_text(data = gates_to_plot,aes(x=(percent_start/100+percent_end/100)/2,y=0.5,color=phase, label=phase ),vjust=-0.55)+theme_classic()+
    theme(legend.position = 'none',axis.line.x.top = element_blank(),
          text = element_text(family = 'Courier'))+ylab('Active Forks')+xlab('Replicated DNA')+
    scale_y_continuous(labels = scales::percent)+
    scale_x_continuous(labels = scales::percent)


results%>%
    group_by(time)%>%
    summarise(total_limiting_factors=mean(total_limiting_factors),
              produced_limiting_factors=mean(produced_limiting_factors),
              Free_factor=mean(Free_factor))%>%
    ggplot()+
    geom_line(aes(time/60,Free_factor/max(total_limiting_factors),color='Free'))+
    geom_line(aes(time/60,total_limiting_factors/max(total_limiting_factors),color='Total'))+
    geom_line(aes(time/60,produced_limiting_factors/max(total_limiting_factors),color='Newly Produced'))+
    geom_line(aes(time/60,(total_limiting_factors-Free_factor)/max(total_limiting_factors),color='Bound to forks'))+
    ylab('Limiting factors')+ xlab('time (h)')+theme_classic()+theme(legend.title = element_blank(),legend.position = 'top',
                                                                     text = element_text(family = 'Courier'))+
    scale_y_continuous(labels = scales::percent)+
    guides(col = guide_legend(nrow = 2))

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

all_phase = read_tsv('All_normalised_Tracks_50adj.tsv') %>%
    mutate(
        phase=case_when(
            phase=='G1'~'G1/S1',
            phase=='S1'~'S2',
            phase=='S2'~'S3',
            phase=='S3'~'S4',
            phase=='S4'~'S5',
            phase=='G2'~'S6/G2/M'

        )
    )%>%
    inner_join(results)%>%
    group_by(Condition,phase)%>%
    mutate(
        reads =Normalize_factor* reads / sum(reads*50)
    )

all_phase%>%write_tsv('All_normalised_Tracks_50adj_normalised_with_sim_periodic_1kbXY.tsv')

x=all_phase%>%filter(reads>0)%>%dplyr::select(chr,start,end,reads,Condition,phase)%>%
    spread(Condition,reads,fill = NA)%>%
    drop_na()%>%
    ungroup()%>%
    dplyr::select(-chr,-start,-end,-phase)%>%
    cor()
r_names=rownames(x)
rownames(x)=case_when(r_names=='Aph1' ~ 'Aph rep 1',
                      r_names=='Aph2' ~ 'Aph rep 2',
                      r_names=='AphRO1' ~ 'ARO rep 1',
                      r_names=='AphRO2' ~ 'ARO rep 2',
                      r_names=='HU1' ~ 'HU',
                      r_names=='NT1' ~ 'NT rep 1',
                      r_names=='NT2' ~ 'NT rep 2',
                      r_names=='NT3' ~ 'NT rep 3',
                      r_names=='ATRiHU1' ~ 'Aph + HU'
                      )
c_names=colnames(x)
colnames(x)=case_when(r_names=='Aph1' ~ 'Aph rep 1',
                      r_names=='Aph2' ~ 'Aph rep 2',
                      r_names=='AphRO1' ~ 'Aph + RO rep 1',
                      r_names=='AphRO2' ~ 'Aph + RO rep 2',
                      r_names=='HU1' ~ 'HU',
                      r_names=='NT1' ~ 'NT rep 1',
                      r_names=='NT2' ~ 'NT rep 2',
                      r_names=='NT3' ~ 'NT rep 3',
                      r_names=='ATRiHU1' ~ 'Aph + HU'
)

corrplot::corrplot(x,method = 'color',addCoef.col = "red",is.corr = F)

#calculate URI
URI = all_phase %>%
    filter(reads > 0)%>%
    group_by(chr, start, end, Condition) %>%
    summarise(reads = sum(reads, na.rm = T)) %>%
    ungroup() %>%
    filter(reads > 0)%>%
    mutate(Condition = str_remove(Condition, '[0-9]$')) %>%
    group_by(chr, start, end, Condition) %>%
    summarise(reads = mean(reads, na.rm = T))%>%
    spread(Condition, reads,fill = NA) %>%
    ungroup() %>%
    mutate(
        delta_Aph_NT = Aph - NT,
        delta_AphRO_NT = AphRO - NT,
        delta_HU_NT = HU - NT,
        delta_ATRiHU_NT = ATRiHU - NT,
        URI_Aph = (delta_Aph_NT - mean(delta_Aph_NT, na.rm = T)) / sd(delta_Aph_NT, na.rm = T),
        URI_AphRO = (delta_AphRO_NT - mean(delta_AphRO_NT, na.rm = T)) /
            sd(delta_AphRO_NT, na.rm = T),
        URI_HU = (delta_HU_NT - mean(delta_HU_NT, na.rm = T)) / sd(delta_HU_NT, na.rm = T),
        URI_ATRiHU = (delta_ATRiHU_NT - mean(delta_ATRiHU_NT, na.rm = T)) / sd(delta_ATRiHU_NT, na.rm = T)
    ) %>%
    dplyr::select(chr, start, end, URI_Aph, URI_AphRO, URI_HU,URI_ATRiHU) %>%
    gather(Condition, URI, URI_Aph, URI_AphRO, URI_HU,URI_ATRiHU) %>%
    mutate(Condition = str_remove(Condition, 'URI_'))

URI%>%write_tsv('URi_simulation_periodic_1kbXY.tsv')

#calculate S50
S50 = all_phase %>%ungroup()%>%
    mutate(Condition = str_remove(Condition, '[0-9]$')) %>%
    group_by(chr, start, end, Condition,phase) %>%
    summarise(reads = mean(reads, na.rm = T))%>%
    dplyr::select(chr,start,end,Condition,phase,reads)%>%
    mutate(phase = factor(phase, levels = c('G1/S1', 'S2', 'S3', 'S4', 'S5', 'S6/G2/M'))) %>%
    spread(phase, reads, fill = 0)

calc_s50=function( G1, S1, S2, S3, S4, G2)    {
    list = rep(c(G1, S1, S2, S3, S4, G2), each = 1000)
    sum = sum(list)/2
    list =  cumsum(list)
    list = min(which(list >= sum))
    list = round(list / 6000, 3)
    return(list)
}

S50=foreach(COND=unique(S50$Condition),.combine = 'rbind')%do%{
    S50%>%
        filter(Condition==COND)%>%
        group_by(chr,start,end,Condition)%>%
        summarise(s50=calc_s50( `G1/S1`, S2, S3, S4, S5, `S6/G2/M`))
}

S50%>%write_tsv('S50_simulation_periodic_1kbXY.tsv')
