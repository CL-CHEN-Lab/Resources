output_dir='~/Desktop/Figure_paper_to_update/FigS3'
system(paste0('mkdir ',output_dir))
library(tidyverse)
library(zoo)
library(gridGraphics)
library(cowplot)

# general pic settings
tex_size = 7
label_text_size=11
font = 'Arial'
fontface = "plain"
line_size = 0.25

general_theme = theme_classic() + theme(
    plot.background = element_blank(),
    strip.text = element_text(
        size = tex_size,
        family = font,
        face = fontface,
        color = 'black'
    ),
    legend.position = 'top',
    panel.background = element_blank(),
    text = element_text(
        family = font,
        size = tex_size,
        face = fontface,
        color = 'black',
        margin = margin(b=0)
    ),
    strip.text.y = element_text(angle = 90),
    legend.text = element_text(
        family = font,
        size = tex_size,
        face = fontface,
        color = 'black'
    ),
    axis.text = element_text(
        family = font,
        size = tex_size,
        face = fontface,
        color = 'black'
    ),
    axis.title = element_text(
        family = font,
        size = tex_size,
        face = fontface,
        color = 'black'
    ),line = element_line(size = line_size),
    strip.background = element_blank(),
    panel.border = element_blank(),
    panel.spacing=unit(1,'mm'),
    legend.key.size = unit(8,'pt'),
    plot.margin=unit(c(0,0,0,0),"cm"),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(0,0,0,0), 
    legend.spacing = unit(0, 'cm'),
    legend.background = element_blank(),
    axis.title.y.right = element_text(
        family = font,
        size = tex_size,
        face = fontface,
        color = 'black',angle = 90
    ),
    legend.title = element_text(vjust = 1.25)
)


#load FACS profile
facs = read_csv('./FACS/AphRO.csv', col_names = F)

gates = tibble(
    phase = c('G1_all', 'G1/S1', 'S2', 'S3', 'S4', 'S5', 'S6/G2/M'),
    start = c(
        21.6,
        24.2519169380677,
        27.282955036018755,
        30.61406625807228,
        34.19265411124849,
        37.47504529876629,
        40.60739992423285
    ),
    end = c(
        26.4,
        27.282955036018755,
        30.61406625807228,
        34.19265411124849,
        37.47504529876629,
        40.60739992423285,
        47.27059148295092
    ),
    G1_peak = 24.15,
    G2_peak = 44.2
)

gates_to_plot = gates %>% filter(phase != 'G1_all') %>% mutate(m = (start +
                                                                        end) / 2) %>%
    mutate(
        percent_start = 100 * (start - G1_peak) / (G2_peak - G1_peak),
        percent_end = 100 * (end - G1_peak) / (G2_peak - G1_peak)
    )
# plot facs profile
plot_facs = facs %>%
    # round x
    mutate(X1 = round(X1, 1)) %>%
    group_by(X1) %>%
    # select min value of y for each x
    summarise(X2 = min(X2)) %>%
    ungroup() %>%
    # smooth the profile
    mutate(X2 = rollmean(X2, 3, fill = NA, align = 'center')) %>%
    #pkit
    ggplot() + geom_line(aes(X1, X2), size = line_size) +
    #add gates
    geom_errorbarh(
        data = gates_to_plot,
        aes(
            xmin = start,
            xmax = end - 0.15,
            y = 270 ,
            height = 30
        ),
        color = 'black',
        size = line_size
    ) +
    #add text
    geom_text(
        data = gates_to_plot,
        aes(
            x = m,
            y = 270,
            label = phase
        ),
        color = 'black',
        size = tex_size / 3,
        family = font,
        fontface = fontface,
        vjust = -0.55
    ) +
    # add second scale 
    scale_x_continuous(sec.axis = sec_axis(
        ~ 100 * (. - unique(gates$G1_peak)) / unique(gates$G2_peak - gates$G1_peak),
        breaks = c(0, 25, 50, 75, 100),
        labels = paste0(c(0, 25, 50, 75, 100), '%') ,
        name = 'Replicated DNA'
    )) + xlab('Fluorescent units') + ylab('cells counts') + general_theme +
    theme(legend.position = 'none', axis.line.x.top = element_blank()) +
    # dro second x ax
    geom_segment(x = 23,
                 xend = 45,
                 y = Inf,
                 yend = Inf)

#load results form simulation
results = read_tsv('Results_periodic_simulationXY_1kb.tsv')

# id highest point normalize 
MAX_ratio = results %>% group_by(cycles) %>% summarise(
    Initiations = mean(Initiations, na.rm = T),
    unreplicated_dna = mean(unreplicated_dna, na.rm = T)
) %>%
    mutate(Ratio_firing_unreplicated = Initiations / unreplicated_dna) %>%
    pull(Ratio_firing_unreplicated) %>% max(na.rm = T)


# firing has been smoothed
simulation_percent_firing = results %>%
    # per cicle calculate mean values of the parameters 
    group_by(cycles) %>%
    summarise(
        Initiations = mean(Initiations),
        unreplicated_dna = mean(unreplicated_dna),
        percent = mean(percent) / 100
    ) %>%
    # caclulate firing over unreplicated dna and smooth
    mutate(
        Ratio_firing_unreplicated = Initiations / unreplicated_dna,
        Ratio_firing_unreplicated = rollmean(
            Ratio_firing_unreplicated,
            k = 15,
            fill = 0,
            align = 'left'
        )
    ) %>%
    ggplot() +
    # plot cycle vs normalized ratio 
    geom_line(
        aes(cycles,
            Ratio_firing_unreplicated / MAX_ratio),
        color = 'red',
        size = line_size
    ) +
    # add percentage of rep DNA 
    geom_line(aes(cycles, percent), color = 'blue', size = line_size) +
    # add scales 
    scale_y_continuous(
        name = 'Replicated DNA',
        labels = scales::percent,
        sec.axis = sec_axis(name = expression(paste(
            frac('Origins Firing', 'unreplicated DNA'),
            ' (',
            Mb ^ -1,
            min ^ -1,
            ')'
        )), trans = ~ . * MAX_ratio*1000)
    ) +
    #add label
    xlab('Cycles') + general_theme + 
    #set colors
    theme(
        axis.title.y.left  = element_text(
            color = 'blue',
            family = font,
            size = tex_size,
            face = fontface
        ),
        axis.title.y.right  =  element_text(
            color = 'red',
            family = font,
            size = tex_size,
            face = fontface
        )
    )

# number of active forks
forks = results %>%
    # average percent and forks over cycle 
    group_by(cycles) %>%
    summarise(percent = mean(percent),
              Forks = mean(Forks)) %>%
    #plot percent vs normalise nforks 
    ggplot() +
    geom_line(aes(percent / 100, Forks / max(Forks)), size = line_size) +
    #add gates
    geom_errorbarh(
        data = gates_to_plot,
        aes(
            xmin = percent_start / 100,
            xmax = (percent_end - 0.6) / 100,
            y = 0.5,
           
            height = 0.05
        ),
        color = 'black' ,
        size = line_size
    ) +
    #add text 
    geom_text(
        data = gates_to_plot,
        aes(
            x = (percent_start / 100 + percent_end / 100) / 2,
            y = 0.5,
            label = phase
        ),color = 'black',
        vjust = -0.55,
        size = tex_size / 3,
        family = font,
        fontface = fontface
    ) + general_theme +
    theme(legend.position = 'none', axis.line.x.top = element_blank()) +
    #labels
    ylab('Active Forks density') + xlab('Replicated DNA') +
    scale_x_continuous(labels = scales::percent)

#load data

all_phase = read_tsv('All_normalised_Tracks_50adj_normalised_with_sim_periodic_1kbXY_ATRiadd.tsv')%>%
    filter(
        Condition %in% c("NT2","NT3","NT4","Aph1","Aph4","AphRO1","AphRO2","HU1","ATRiHU1")
    )

x = all_phase %>% 
    # remove missing values from each track
    filter(reads > 0) %>% 
    #select columns 
    dplyr::select(chr, start, end, reads, Condition, phase) %>%
    # exclude conditions without replcates
    filter(!Condition %in% c('HU1',"ATRiHU1","ATRi1"))%>%
    # spred and drop incomplete lines 
    spread(Condition, reads, fill = NA) %>%
    drop_na() %>%
    ungroup() %>%
    # keep only rads counts
    dplyr::select(-chr, -start, -end, -phase) %>%
    # calculate pearson
    cor()

#assign names
r_names = rownames(x)
rownames(x) = case_when(
    r_names == 'Aph1' ~ 'Aph 1',
    r_names == 'Aph4' ~ 'Aph 2',
    r_names == 'AphRO1' ~ 'ARO 1',
    r_names == 'AphRO2' ~ 'ARO 2',
    r_names== "ATRi1" ~ 'ATRi',
    r_names=="ATRiHU1" ~ 'ATRiHU',
    r_names == 'NT2' ~ 'NT 1',
    r_names == 'NT3' ~ 'NT 2',
    r_names == 'NT4' ~ 'NT 3'
)
c_names = colnames(x)
colnames(x) = case_when(
    c_names == 'Aph1' ~ 'Aph 1',
    c_names == 'Aph4' ~ 'Aph 2',
    c_names == 'AphRO1' ~ 'ARO 1',
    c_names == 'AphRO2' ~ 'ARO 2',
    c_names == 'NT2' ~ 'NT 1',
    c_names == 'NT3' ~ 'NT 2',
    c_names == 'NT4' ~ 'NT 3',
    c_names== "ATRi1" ~ 'ATRi',
    c_names=="ATRiHU1" ~ 'ATRiHU'
)

# plot heatmap
p = x %>% as_tibble() %>%
    # add row names 
    mutate(against = rownames(x)) %>%
    # reshape 
    gather(Sample, Value, -against) %>%
    # order by colnames levels
    mutate(
        Sample = factor(Sample, levels = colnames(x)),
        against = factor(against, levels = colnames(x))
    ) %>%
    #plot
    ggplot() +
    geom_tile(aes(against, Sample, fill = Value)) +
    # plot value
    geom_text(
        aes(against, Sample, label = round(Value, 2)),
        size = tex_size / 3,
        family = font,
        fontface = fontface,
        color = 'white'
    ) + 
    # style
    general_theme +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1,
            family = font,
            size = tex_size,
            face = fontface,
            color = 'black'
        ),
        axis.title = element_blank(),
        legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.position = 'right'
    ) +
    # color gradient 
    scale_fill_gradient(low = 'green', high = 'darkgreen', name = 'Pearson correlation', oob=squish,breaks=c(0.5,0.75,1)) +
    # impose square shape
    coord_equal()

#load data of conditions with replicas
all_phase = read_tsv('All_normalised_Tracks_50adj_normalised_with_sim_periodic_1kbXY_ATRiadd.tsv')%>%
    filter(
        Condition %in% c("NT2","NT3","NT4","Aph1","Aph4","AphRO1","AphRO2")
    )

#load URI of conditions with replicas
URI = read_tsv('URi_simulation_periodic_1kbXY_ATRIadd.tsv') %>%
    filter(!Condition %in% c('HU',"ATRiHU","ATRi"))
    

FHIT = tibble(
    chr = 'chr3',
    start = 59735035,
    end = 61237087,
    gene='FHIT'
)
WWOX = tibble(chr = 'chr16',
              start = 78133310,
              end = 79246567,
              gene='WWOX')

plot_fragile_site = function(Tracks,
                             URI,
                             position = tibble(chr = 'chr1', start = 0, end = 10000000),
                             line_size = 1.5,
                             style = general_theme,
                             legend_posistion = 'top',
                             font = 30,
                             tex_size = 'Courier',
                             fontface = 'bold',
                             flanking = 2 * 10 ^ 6,
                             pos_heatmap=0.5,
                             heatmap_limits=c(-3,3), 
                             n_fill_breaks=3,
                             more_or_less_label=c('N','L','M','B')) {
    
    more_or_less_label=more_or_less_label[1]
    n_fill_breaks=n_fill_breaks-1
    
    #invert start and end if needed
    position = position %>%
        mutate(s = ifelse(start > end, end, start),
               e = ifelse(start > end, start, end)) %>%
        dplyr::select(chr, 'start' = s, 'end' = 'e',gene)
    
    # select bins of interest over tracks ± flanking
    tracks = Tracks %>% ungroup() %>%
        filter(
            chr == position$chr,
            start >= position$start - flanking,
            end <= position$end + flanking
        ) %>%
        # convert phases into factors
        mutate(phase = factor(phase, levels = c(
            'G1/S1', 'S2', 'S3', 'S4', 'S5', 'S6/G2/M'
        ))) %>%
        # average conditions
        mutate(Condition_2 = str_remove(Condition, '[1-9]{1,2}$')) %>%
        group_by(Condition_2, chr, start, end, phase) %>%
        mutate(reads_mean = mean(reads)) %>%
        ungroup()%>%
        #Chnge names of ARO and ATRiHU
        mutate(Condition_2 = factor(
            case_when(Condition_2 == 'AphRO' ~ 'ARO',
                      Condition_2 == "ATRiHU" ~ 'HU+ATRi',
                      T ~ Condition_2),
            levels = c('NT', 'Aph', 'ARO')
        ),
        #calculae center of the bin
        cent = (start + end) / 2)
    
    #average track
    track_mean=tracks%>%dplyr::select(Condition,Condition_2,cent,chr,start,end,reads_mean,phase)%>%unique()
    # remove average from the other tracks 
    tracks=tracks%>%dplyr::select(-reads_mean)
   
     # calculate max to set ylim
    max = round(max(tracks %>% pull(reads)))
    
    #URI track 
    URi_toplot = URI %>%
        #select region of interest
        filter(chr == position$chr,
               start >= position$start -flanking,
               end <= position$end + flanking) %>%
        # create a NT track of NA and add it to the data
        mutate(NT = NA) %>%
        spread(Condition, URI) %>%
        gather(Condition, URI, -chr, -start, -end) %>%
        # add a name for the track
        mutate(phase = 'URI',
               rep = ' ') %>%
        #calculate the center of the track 
        mutate(cent = start + (end - start) / 2) %>%
        #rename and convert into factors
        mutate(Condition = factor(
            case_when(Condition == 'AphRO' ~ 'ARO',
                      T ~ Condition),
            levels = c('NT', 'Aph', 'ARO')
        ))
    # max uri to set limits
    URI_max=max(URi_toplot$URI,na.rm=T)
    
    
    p1 = 
        ggplot() +
        # prepare facet
        facet_grid(phase ~ Condition_2) +
        #plot single tracks
        geom_line(data=tracks,aes(x = cent, y = reads, group = Condition),
                  color = '#636466',
                  size = line_size) +
        #plot average tracks
        geom_line(data=track_mean,aes(x = cent, y = reads_mean, color = Condition_2), size =
                      line_size)+
        #genes positions
        geom_vline(xintercept = position$end,
                   color = 'purple', size =line_size) +
        geom_vline(xintercept = position$start,
                   color = 'purple', size =line_size) +
        #general style
        style +
        theme(
            strip.background.x = element_blank(),
            strip.text.x = element_blank(),
            
            legend.position = 'none',
            panel.background = element_blank(),
            axis.line.y = element_line(color = 'black'),
            axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                vjust = 1,
                family = font,
                size = tex_size,
                face = fontface,
                color = 'black'
            )
        ) +
        # plot x ax line
        geom_hline(yintercept = -Inf, size =line_size, color='black') +
        # format axis
        scale_x_continuous(
            labels =  function(x)
                paste(round(x / 1000000, 2), 'Mb')
        ) +
        # assign colors
        scale_color_manual(values = c(
            'Aph' = '#95C11F',
            'ARO' = '#F9B233',
            'HU' = '#CA9E67',
            'NT' = '#36A9E1',
            'HU+ATRi'='#7D4E24'
        )) + 
        # set xlab
        xlab('')+
        # fix coord
        coord_cartesian(xlim = c(position$start - flanking, position$end + flanking))+
        # add side line
        geom_vline(xintercept = c(-Inf,Inf), size =line_size)
    
    #limits
    segment=tibble(x=c(-Inf,Inf),y=c(-Inf,-Inf),xend=c(-Inf,Inf), yend=c(URI_max,URI_max),Condition=list(c('NT', 'Aph', 'ARO')))%>%
        unnest(cols = c(Condition))%>%
        mutate(Condition=factor(Condition,            
                                levels = c('NT', 'Aph', 'ARO')
        ))
    
    
    # color breaks for the heatmap 
    colo_breaks=seq(heatmap_limits[1], heatmap_limits[2], sum(abs(heatmap_limits)) /
                        n_fill_breaks)
    # labels for the colors of the heatmap ( if the signal is saturated over that value)
    labels_colors=tibble(color_breks=colo_breaks,more_or_less_label=more_or_less_label)%>%
        mutate(
            labels_colors=case_when(
                
                more_or_less_label %in% c('B','L') & color_breks==min(color_breks) ~ paste0('≤',colo_breaks),
                more_or_less_label %in% c('B','M') & color_breks==max(color_breks) ~ paste0('≥',colo_breaks),
                T ~ as.character(colo_breaks)
                
            ))%>%pull(labels_colors)
    
    # URI plot 
    p2 = URi_toplot %>%
        ggplot() +
        # draw box
        geom_segment(data=segment,aes(x=x,xend=xend, y=y,yend=yend), size =line_size,color='black')+
        # plot URI as line 
        geom_line(aes(x = cent, y = URI), color = 'black', size =
                      line_size) +
        #plot URI as heatmap
        ggplot2::geom_rect(aes(
            xmin = start,
            xmax = end,
            ymin = unique(segment$yend)+pos_heatmap,
            ymax = unique(segment$yend)+pos_heatmap+1,
            fill = URI
        )) +
        # color gradient
        scale_fill_gradient2(
            low = 'red',
            high = 'blue',
            mid = 'white',
            na.value = "grey", oob=squish,limits=heatmap_limits,breaks=colo_breaks,labels=labels_colors
        )+
        # y= -2 line 
        geom_hline(yintercept = -2, color = 'red', size =line_size) +
        # split conditiona and phase
        facet_grid(phase ~ Condition) +
        # add y label 
        ylab('URI') +
        # plot style
        style +
        theme(
            legend.position = 'none',
            axis.line.y = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x  = element_blank(),
            axis.text.x = element_blank(),
            axis.title = element_text(color = 'white')
        ) +
        # add lines for the axes
        geom_hline(yintercept = -Inf, size =line_size,col='black') +
        #add lines for the gene body
        geom_vline(xintercept = position$end,
                   color = 'purple', size =
                       line_size) +
        geom_vline(xintercept = position$start,
                   color = 'purple', size =
                       line_size) +
        #add line for the x axis
        geom_hline(yintercept = -Inf, size =line_size) + 
        scale_y_continuous(breaks = c(-2, 0, 2)) +
        #add limits
        coord_cartesian(xlim = c(position$start - flanking, position$end + flanking))
    

    # return plot and legend
    return(list(Plot=plot_grid(
        p2,
        p1,
        rel_heights = c (3, 15),
        ncol = 1,
        align = 'v',
        axis = "l"
    ),legend=get_legend(p2+theme(legend.position = 'top'))))
}

FRA3B = plot_fragile_site(
    all_phase,
    URI,
    FHIT,
    line_size = line_size,
    font = font,
    tex_size = tex_size,
    fontface = fontface,
    legend_posistion = 'none'
)
FRA16D = plot_fragile_site(
    all_phase,
    URI,
    WWOX,
    line_size = line_size,
    font = font,
    tex_size = tex_size,
    fontface = fontface,
    legend_posistion = 'none'
)


Fragiale_sites = plot_grid(
    FRA3B$legend,
    rel_heights =c(0.1,2),
    ncol=1,
    plot_grid(
    plot_grid(
    textGrob('FRA3B (FHIT)',gp = gpar(
        fontsize = tex_size,
        fontface = fontface,
        fontfamily = font
    )),
    FRA3B$Plot,
    ncol = 1,
    scale = 1,
    rel_heights =  c(0.05, 1)),
plot_grid(
    textGrob('FRA16D (WWOX)',gp = gpar(
        fontsize = tex_size,
        fontface = fontface,
        fontfamily = font
    )),
    FRA16D$Plot,
    ncol = 1,
    scale = 1,
    rel_heights =  c(0.05, 1)
)))

simulation_fig = plot_grid(
    plot_facs,
    simulation_percent_firing,
    forks,
    p ,
    labels = c("A", "B", "C", "D"),
    label_size = label_text_size,
    label_fontfamily = font,
    scale = 0.95
)

ggsave(
    plot = plot_grid(
        simulation_fig,
        Fragiale_sites,
        ncol = 1,
        rel_heights = c(1, 1),
        labels = c('', 'E'),
        scale = 0.99,
        label_size = label_text_size,
        label_fontfamily = font
    ),
    filename = paste0(output_dir,'/FigS3.pdf'),
    width = 174,
    height = 220,
    limitsize = F,
    units = 'mm',
    device = cairo_pdf
)

