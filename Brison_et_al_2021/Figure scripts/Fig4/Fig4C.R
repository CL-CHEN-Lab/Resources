#########################################################
#########################################################

# RUN Fig4A.R and Fig4B first


#########################################################
#########################################################

output_dir='~/Desktop/Figure_paper_to_update/Fig4/Fig4C'
system(paste0('mkdir -p ',output_dir))

average_profile = function(matrix,
                           name = 'signal',
                           group = 'genes',
                           breaks = 'auto',
                           labels = c('-100', 'TSS', 'TTS', '100'),
                           limits = NA,
                           tex_size = 30,
                           font = 'Courier',
                           fontface = 'bold',
                           style = general_theme,
                           line_size = 1.5,
                           legend_position = 'top',
                           Colors=NA) {
    
    #if the used does not provide info about the x axis breaks it looks  for the values into the matrix
    if (breaks == 'auto') {
        upstream = sum(grepl(x = colnames(matrix), pattern = 'u'))
        body = sum(grepl(x = colnames(matrix), pattern = 't'))
        dowstream = sum(grepl(x = colnames(matrix), pattern = 'd'))
        
        breaks = c(0, upstream + 1, upstream + body, upstream + body + dowstream) %>%
            unique()
    }
    
    # convert matrix into a df
    x = matrix %>% as_tibble() %>%
        #assign group and aname to each line
        mutate(group =group,
               name = name) %>%
        # reshape data
        gather(pos, value, -group, -name) %>%
        #calculate the mean value per position per group (the name is the same everywhere)
        group_by(group, name, pos) %>%
        summarise(value = mean(value, na.rm = T)) %>%
        ungroup() %>%
        #convert position to factor and than number
        mutate(pos = factor(pos, levels = colnames(matrix)),
               pos = as.numeric(pos)) %>%
        #arrange by group and position
        arrange(group, pos)

    #start plotting
    x = ggplot(x) +
        #line position vs value with color depending on group
        geom_line(aes(
            x = pos,
            y = value,
            color = group
        ), size = line_size) + 
        #facet grid given name 
        facet_grid( ~ name) +
        #set x axis breaks 
            scale_x_continuous(breaks = breaks, labels = labels) +
        # general style 
            style +
            theme(
                legend.position = legend_position,
                legend.title = element_blank(),
                axis.text.x = element_text(
                    angle = 90,
                    hjust = 0.5,
                    vjust = 0.5,
                    family = font,
                    size = tex_size,
                    face = fontface,
                    color = 'black'
                ),
                panel.spacing=unit(0,'mm')
            )+
        #set labels to empty
        xlab('') + ylab('') +
        # add dashed lines 
        geom_vline(xintercept = breaks[c(2, 3)],
                   linetype = 'dashed',
                   color = 'black', size =line_size)+
        #add black rectangles around the plot
        geom_rect(data = tibble(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill=NA,color='black', size =line_size)+
        theme( axis.line = element_blank())
    if (!any(is.na(limits))) {
        # if limits are give apply them
        x = x + coord_cartesian(ylim = limits)
    }
    if (!any(is.na(Colors))) {
        # if colors are given pass them to the plot
        x = x + scale_color_manual(values = Colors)
    }
    
    return(x)
}


#select genes that are bigger than 300 that are not SDR
region_300 = genes[genes$split_size == "≥300kb" & genes$split!='SDR']

# calculate matrixes 
rt_300 = normalizeToMatrix(
    target = region_300,
    signal =  NT,
    value_column = 's50',
    mean_mode = "w0",
    extend = 100000,
    target_ratio = 0.45
)
exp_300 = normalizeToMatrix(
    target = region_300,
    signal =  bins,
    value_column = 'counts',
    mean_mode = "w0",
    extend = 100000,
    target_ratio = 0.45
)
mat_Aph_300 = normalizeToMatrix(
    target = region_300,
    signal =  URI,
    value_column = 'Aph',
    extend = 100000,
    mean_mode = "w0",
    target_ratio = 0.45
)
mat_AphRO_300 = normalizeToMatrix(
    target = region_300,
    signal =  URI,
    value_column = 'AphRO',
    extend = 100000,
    mean_mode = "w0",
    target_ratio = 0.45
)
mat_HU_300 = normalizeToMatrix(
    target = region_300,
    signal =  URI,
    value_column = 'HU',
    extend = 100000,
    mean_mode = "w0",
    target_ratio = 0.45
)
mat_ATRiHU_300 = normalizeToMatrix(
    target = region_300,
    signal =  URI,
    value_column = 'ATRiHU',
    extend = 100000,
    mean_mode = "w0",
    target_ratio = 0.45
)

#colors
Colors=c(
    '1st'='#0000FF',
    '2nd'='#3AAA35',
    '3rd'='#FF00FF',
    '4th'='#951B81'
)
Heat_maps_300 = plot_grid(
    
    plot_grid(
      
            average_profile(
                exp_300,
                name = 'Gro-seq',
                group = region_300$split_expression,
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            )+scale_y_continuous(breaks=c(0,1000,2000),labels = c(0,'1K','2K'))+
                theme(strip.text  = element_blank()),
        
      
            average_profile(
                rt_300,
                name = 'RT (NT)',
                group = region_300$split_expression,
                limits = c(1, 0),
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            ) + scale_y_continuous(breaks=c(0,0.5,1))+
                theme(strip.text  = element_blank()),
      
            average_profile(
                mat_Aph_300,
                name = 'URI (Aph)',
                group = region_300$split_expression,
                limits = c(-3, 3),
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            )+
                theme(strip.text  = element_blank()),
  
            average_profile(
                mat_AphRO_300,
                name = 'URI (ARO)',
                group = region_300$split_expression,
                limits = c(-3, 3),
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            )+
                theme(strip.text  = element_blank()),

            average_profile(
                mat_HU_300,
                name = 'URI (HU)',
                group = region_300$split_expression,
                limits = c(-3, 3),
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            )+
                theme(strip.text  = element_blank()),
   
            average_profile(
                mat_ATRiHU_300,
                name = 'URI (HU+ATRi)',
                group = region_300$split_expression,
                limits = c(-3, 3),
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            )+
                theme(strip.text  = element_blank()),
        nrow = 1,
        rel_widths = c(1.03, 1.05, 1, 1, 1, 1, 1, 0.3)
    ),
    get_legend(
        average_profile(
            exp_300,
            name = 'Gro-seq',
            group = region_300$split_expression,
            tex_size = tex_size,
            font = font,
            fontface = fontface,
            Colors = Colors
        )
    ),
    ncol = 1,
    rel_heights = c( 2,0.1),scale = 0.9
)
#save
ggsave(
    plot =Heat_maps_300,
    filename = paste0(output_dir,'/Fig4C.pdf'),
    width = 17,
    height = 3.5,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)
