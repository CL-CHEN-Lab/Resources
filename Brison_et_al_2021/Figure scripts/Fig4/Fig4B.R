#########################################################
#########################################################

# RUN Fig4A.R first


#########################################################
#########################################################
output_dir='~/Desktop/Figure_paper_to_update/Fig4/Fig4B'
system(paste0('mkdir -p ',output_dir))

average_profile_small = function(matrix,
                                 color = 'black',
                                 name = 'signal',
                                 breaks = 'auto',
                                 labels = c('-100', 'TSS', 'TTS', '100'),
                                 limits = NA,
                                 font = 30,
                                 tex_size = 'Courier',
                                 fontface = 'bold',
                                 style = general_theme,
                                 line_size = 1.5) {
    
    #if the used does not provide info about the x axis breaks it looks  for the values into the matrix
    if (breaks == 'auto') {
        upstream = sum(grepl(x = colnames(matrix), pattern = 'u'))
        body = sum(grepl(x = colnames(matrix), pattern = 't'))
        dowstream = sum(grepl(x = colnames(matrix), pattern = 'd'))
        
        breaks = c(0, upstream + 1, upstream + body, upstream + body + dowstream) %>%
            unique()
    }
    
    #calculate mean of each colum in the matrix and assign position in a df
    x = tibble(value = colMeans(matrix),
               pos = colnames(matrix)) %>%
        # covert position in Factors and convert them into numeric
        mutate(pos = factor(pos, levels = pos),
               pos = as.numeric(pos)) %>%
        # assign a name to the plot
        mutate(name = name) %>%
        #plot
        ggplot(aes(pos, value)) + geom_line(color = color, size = line_size) +
        facet_grid( ~ name) +
        #set breaks in the x axis
        scale_x_continuous(breaks = breaks, labels = labels) +
        #add style 
        style +
        theme(
            legend.position = 'none',
            legend.title = element_blank(),
            axis.text.x = element_text(
                angle = 90,
                hjust = 1,
                vjust = 1,
                family = font,
                size = tex_size,
                face = fontface,
                color = 'black'
            )
        ) +
        #set ax labels to empty 
        xlab('') + ylab('') +
        # add dashed lines on breaks 
        geom_vline(xintercept = breaks[c(2, 3)],
                   linetype = 'dashed',
                   color = 'black', size =line_size)+
        # add rectable around plot 
        geom_rect(
            data = tibble(
                xmin = -Inf,
                xmax = Inf,
                ymin = -Inf,
                ymax = Inf
            ),
            aes(
                xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax
            ),
            fill = NA,
            color = 'black',
            size = line_size,
            inherit.aes =F
        )+
        # font settings
        theme(axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 1,
            family = font,
            size = tex_size,
            face = fontface,
            color = 'black'
        ),
        axis.line = element_blank())
    
    if (!any(is.na(limits))) {
        # impose limits if any
        x = x + coord_cartesian(ylim = limits)
    }
    return(x)
}

#assemble plot
Heat_maps_small = plot_grid(
        average_profile_small(
            exp_sdr,
            name = 'Gro-seq',
            tex_size = tex_size,
            font = font,
            fontface = fontface,
            line_size = line_size
        )+scale_y_continuous(breaks=c(0,1000,1500),labels = c(0,'1K','1.5K'))+ theme(strip.text = element_blank()),
    

        average_profile_small(
            rt_sdr,
            name = 'RT (NT)',
            limits = c(0.8, 0.4),
            tex_size = tex_size,
            font = font,
            fontface = fontface,
            line_size = line_size
        ) + theme(strip.text = element_blank() ),
        average_profile_small(
            mat_Aph_sdr,
            name = 'URI (Aph)',
            limits = c(-3, 3),
            tex_size = tex_size,
            font = font,
            fontface = fontface,
            line_size = line_size
        )+ theme(strip.text = element_blank()),

        average_profile_small(
            mat_AphRO_sdr,
            name = 'URI (ARO)',
            limits = c(-3, 3),
            tex_size = tex_size,
            font = font,
            fontface = fontface,
            line_size = line_size
        )+ theme(strip.text = element_blank()),
  
        average_profile_small(
            mat_HU_sdr,
            name = 'URI (HU)',
            limits = c(-3, 3),
            tex_size = tex_size,
            font = font,
            fontface = fontface,
            line_size = line_size
        )+ theme(strip.text = element_blank()),
   
     
        average_profile_small(
            mat_ATRiHU_sdr,
            name = 'URI (HU+ATRi)',
            limits = c(-3, 3),
            tex_size = tex_size,
            font = font,
            fontface = fontface,
            line_size = line_size
        )+ theme(strip.text = element_blank()),
      
    nrow = 1,
    rel_widths = c(1.1, 1.06, 1, 1, 1,1),scale = 0.9
)

#save
ggsave(
    plot =Heat_maps_small,
    filename = paste0(output_dir,'/Fig4B.pdf'),
    width = 17,
    height = 3,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)
 