output_dir='~/Desktop/Figure_paper_to_update/Fig3/Fig3C'
system(paste0('mkdir ',output_dir))

library(extrafont)
loadfonts(device = "postscript")

library(tidyverse)
library(gridGraphics)
library(cowplot)

#load data
all_phase = read_tsv('All_normalised_Tracks_50adj_normalised_with_sim_periodic_1kbXY_ATRiadd.tsv')%>%
    filter(
        Condition %in% c("NT2","NT3","NT4","Aph1","Aph4","AphRO1","AphRO2","HU1","ATRiHU1")
    )
URI_Fs = read_tsv('URi_simulation_periodic_1kbXY_ATRIadd.tsv')%>%
    filter(
        Condition %in% c("NT","Aph","AphRO","HU","ATRiHU")
    )

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
    axis.line = element_line(size = line_size),
    strip.background = element_blank(),
    panel.border = element_blank(),
    panel.spacing=unit(1,'mm'),
    legend.key.size = unit(8,'pt'),
    plot.margin=unit(c(0,0,0,0),"cm"),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(0,0,0,0), 
    legend.spacing = unit(0, 'cm'),
    legend.background = element_blank(),
    legend.title = element_text(vjust = 1.25)
)


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
                             font = 30,
                             tex_size = 'Courier',
                             fontface = 'bold',
                             line_size = 1.5,
                             style = general_theme,
                             legend_posistion = 'top',
                             flanking = 2 * 10 ^ 6,
                             pos_heatmap = 0.5,
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
        mutate(Condition = str_remove(Condition, '[1-9]{1,2}$')) %>%
        group_by(Condition, chr, start, end, phase) %>%
        summarise(reads = mean(reads))%>%
        ungroup() %>%
        #Chnge names of ARO and ATRiHU
        mutate(Condition = factor(
            case_when(Condition == 'AphRO' ~ 'ARO',
                      Condition == "ATRiHU" ~ 'HU+ATRi',
                      T ~ Condition),
            levels = c('NT', 'Aph', 'ARO', 'HU','HU+ATRi')
        ),
        #calculae center of the bin
        cent = (start + end) / 2)
    
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
                      Condition == "ATRiHU" ~ 'HU+ATRi',
                      T ~ Condition),
            levels = c('NT', 'Aph', 'ARO', 'HU','HU+ATRi')
        ))
    # max uri to set limits
    URI_max=max(URi_toplot$URI,na.rm=T)
    
    #plot tracks 
    p1 = tracks %>%
        ggplot() + geom_area(aes(x = cent, y = reads, fill = Condition),position="identity") +
        #lines to identify the gene body 
        geom_vline(xintercept = position$end,
                   color = 'purple', size =line_size) +
        geom_vline(xintercept = position$start,
                   color = 'purple', size =line_size) +
        #split base on condition and phase
        facet_grid(phase ~ Condition) +
        # general style
        style +
        theme(
            strip.background.x = element_blank(),
            strip.text.x = element_blank(),
            legend.position = 'none',
            panel.background = element_blank(),
            axis.line.y = element_blank(),
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
        # add axis line 
        geom_hline(yintercept = -Inf, size =line_size, color='black') +
        # convert x axis lables 
        scale_x_continuous(
            labels =  function(x)
                paste(round(x / 1000000, 2), 'Mb')
        ) +
        # assign colors
        scale_fill_manual(values = c(
            'Aph' = '#95C11F',
            'ARO' = '#F9B233',
            'HU' = '#CA9E67',
            'NT' = '#36A9E1',
            'HU+ATRi'='#7D4E24'
        )) + 
        #remouve xaxis label
        xlab('')+
        #assign limits
        coord_cartesian(xlim = c(position$start - flanking, position$end + flanking))+
        theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
        # add axis line 
        geom_vline(xintercept = c(-Inf,Inf), size =line_size)
    
    #limits
    segment=tibble(x=c(-Inf,Inf),y=c(-Inf,-Inf),xend=c(-Inf,Inf), yend=c(URI_max,URI_max),Condition=list(c('NT', 'Aph', 'ARO', 'HU','HU+ATRi')))%>%
        unnest(cols = c(Condition))%>%
        mutate(Condition=factor(Condition,            
                                levels = c('NT', 'Aph', 'ARO', 'HU','HU+ATRi')
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
        rel_heights = c( 4, 15),
        ncol = 1,
        align = 'v',
        axis = "l"
    ),legend=get_legend(p2+theme(legend.position = 'top'))))
}

FRA3B = plot_fragile_site(
    Tracks = all_phase,
    URI = URI_Fs,
    position = FHIT,
    tex_size = tex_size,
    font = font,
    fontface = fontface,
    line_size = line_size,
    style = general_theme,
    legend_posistion = 'none',more_or_less_label = 'B'
)
FRA16D = plot_fragile_site(
    Tracks = all_phase,
    URI = URI_Fs,
    position = WWOX,
    tex_size = tex_size,
    font = font,
    fontface = fontface,
    line_size = line_size,
    style = general_theme,
    legend_posistion = 'none',more_or_less_label = 'B'
)

#final assembly
Fragiale_sites = plot_grid(
    FRA3B$legend,
    
    textGrob('FRA3B (FHIT)',gp = gpar(
        fontsize = tex_size,
        fontface = fontface,
        fontfamily = font
    )),
    FRA3B$Plot,
    textGrob('FRA16D (WWOX)',gp = gpar(
        fontsize = tex_size,
        fontface = fontface,
        fontfamily = font
    )),
    FRA16D$Plot,
    ncol = 1,
    scale = 1,
    rel_heights =  c(0.1 ,0.03, 1,0.03,1)
)


ggsave(plot =Fragiale_sites ,filename = paste0(output_dir,'/Fig3C.pdf'),device = cairo_pdf,width = 17 ,height = 20,units = 'cm')
