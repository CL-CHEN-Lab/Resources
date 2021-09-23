output_dir='~/Desktop/Figure_paper_to_update/Fig3/Fig3D'
system(paste0('mkdir ',output_dir))

library(extrafont)
loadfonts(device = "postscript")

library(tidyverse)
library(GenomicRanges)

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
#load URI
URI = read_tsv('URi_simulation_periodic_1kbXY_ATRIadd.tsv')%>%
    filter(
        Condition %in% c("NT","Aph","AphRO","HU","ATRiHU")
    )
#load S50
S50 = read_tsv('S50_simulation_periodic_1kbXY_ATRiadd.tsv')%>%
    filter(
        Condition %in% c("NT","Aph","AphRO","HU","ATRiHU")
    )

#load SDR and t-SDR and merge the annotation
DR = read_delim(
    '/Volumes/Storage2/CFS_reference/SDR.bed',
    col_names = c('chr', 'start', 'end'),
    delim = '\t'
) %>%
    left_join(
        read_delim(
            '/Volumes/Storage2/CFS_reference/t-SDR.bed',
            col_names = c('chr', 'start', 'end'),
            delim = '\t'
        ) %>%
            mutate(DR = 't-SDR')
    ) %>%
    mutate(DR = ifelse(is.na(DR), 'SDR', DR))
#convert URI into grange and associate infos
URI = URI %>% inner_join(S50 %>% filter(Condition == 'NT') %>% dplyr::select(-Condition)) %>%
    mutate(s50 = ifelse(s50 < 0.5, 'Early', 'Late'))

# convert tracks into Gragnes
URI = URI %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
DR = DR %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)

#find overlaps 
hits = findOverlaps(query = DR, subject = URI)
#initialise colum
URI$DR = NA
# check width overlap
w = width(pintersect(URI[subjectHits(hits)],
                     DR[queryHits(hits)]))
# select overlpas bigger than 1bp
hits = hits[w > 1]

# assign DR type
URI$DR[subjectHits(hits)] = DR$DR[queryHits(hits)]

#convert to DF
URI_SDR = URI %>% as_tibble()


TSDR_SDR = URI_SDR %>%
    # Change DR == NA to NON SDR, change reformat names samples and convert to factor
    mutate(
        DR = ifelse(is.na(DR), 'NON SDR', DR),
        Condition = factor(
            case_when(Condition == 'AphRO' ~ 'ARO',
                      Condition == "ATRiHU" ~ 'HU+ATRi',
                      T ~ Condition),
            levels = c('NT', 'Aph', 'ARO', 'HU','HU+ATRi')
        )
    ) %>%
    # keep only late Rep
    filter(!(s50 == 'Early' & DR == 'NON SDR')) %>%
    # rename DR
    dplyr::rename(Region = DR) %>%
    drop_na() %>%
    #plot densities 
    ggplot(aes(x = URI)) + geom_density(aes(color = Region,linetype=Region), size = line_size) +
    coord_cartesian(xlim = c(-4, 4)) +
    facet_grid( ~ Condition) +
    general_theme +
    theme(
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.1, "cm"),
        legend.title = element_blank()
    )+
    # assign colors
    scale_color_manual(
        values = c( "t-SDR"='#E52620', "SDR"='#E52620' ,'NON SDR'='#2B4B9B' )
    )+
    #assign line type
    scale_linetype_manual(
        values = c( "t-SDR"='solid', "SDR"='dashed' ,'NON SDR'='solid' )
    )+theme(legend.position = 'top',
            legend.justification = 'center')+
    #add -2 line
    geom_vline(xintercept = -2,color='red', size = line_size)

ggsave(plot =TSDR_SDR ,filename = paste0(output_dir,'/Fig3D.pdf'),device = cairo_pdf,width = 8 ,height = 4,units = 'cm')
 