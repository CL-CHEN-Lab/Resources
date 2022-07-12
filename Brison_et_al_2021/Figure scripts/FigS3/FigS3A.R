#########################################################
#########################################################

# RUN Fig3A.R, Fig3B and Fig3C first


#########################################################
#########################################################
output_dir='~/Desktop/Figure_paper_to_update/FigS3/FigS3A'
system(paste0('mkdir -p ',output_dir))

heatmap = function(matrix,
                   line_order = NA,
                   name = 'heatmap',
                   limits = c(-Inf, Inf),
                   more_or_less_label=c('N','B','M','L'),
                   low_color = 'red',
                   high_color = 'blue',
                   mid_color = NA,
                   legend.position = 'bottom',
                   group = 'genes',
                   group_order = NA,
                   breaks = 'auto',
                   labels = c('-100', 'TSS', 'TTS', '100'),
                   tex_size = 7,
                   font = 'Courier',
                   fontface = 'bold',
                   style = general_theme,
                   n_fill_breaks = 5,
                   line_size=1.5) {
    
    # parameter used for the color scale labels
    more_or_less_label=more_or_less_label[1]
    
    #if the used does not provide info about the x axis breaks it looks  for the values into the matrix
    if (breaks == 'auto') {
        
        if(is.na(group_order)){
            group_order=unique(group)%>%sort()
        }
        
        upstream = sum(grepl(x = colnames(matrix), pattern = 'u'))
        body = sum(grepl(x = colnames(matrix), pattern = 't'))
        dowstream = sum(grepl(x = colnames(matrix), pattern = 'd'))
        
        breaks = c(0, upstream + 1, upstream + body, upstream + body + dowstream) %>%
            unique()
    }
    
    # convert matrix into df
    matrix = matrix %>% as_tibble() %>%
        # define order of plotting of the groups
        mutate(group = factor(group,levels = group_order)) %>%
        #assign line number
        mutate(line = 1:n())
    
    #if there is a line order it is used to reorder inside the groups
    if (!any(is.na(line_order))) {
        matrix = matrix  %>%
            mutate(line = factor(line, levels = line_order),
                   line = as.numeric(line)) %>%
            arrange(line) %>%
            group_by(group) %>%
            #new order
            mutate(line = 1:n())
    }
    
    
    #reshape for plotting
    matrix = matrix %>%
        gather(position, value, -line, -group) %>%
        mutate(
            position = factor(position, levels = colnames(matrix)),
            position = as.numeric(position),
            value = as.numeric(value)
        )
    
    
    # calculate min max to use to plot rectangles 
    max_min = matrix %>% group_by(group) %>% summarise(
        min_x = min(position),
        max_x = max(position),
        min_y = min(line),
        max_y = max(line)
    )%>%filter(!is.infinite(min_x),!is.infinite(max_x),!is.infinite(min_y),!is.infinite(max_y))
    
    #start plotting
    matrix = ggplot(matrix) +
        geom_tile(aes(
            y = line,
            x = position,
            #if the color value is higher than limits saturate it
            fill = ifelse(
                value < min(limits),
                min(limits),
                ifelse(value > max(limits), max(limits),
                       value)
            )
        ))
    #setting breaks 
    n_fill_breaks = n_fill_breaks - 1
    
    # if a mid color for the heatmap is provided 
    if (is.na(mid_color)) {
        
        # calculate breaks positions
        colo_breaks=seq(limits[1], limits[2], sum(abs(limits)) /
                            n_fill_breaks)
        #create breaks labels 
        labels_colors=tibble(color_breks=colo_breaks,more_or_less_label=more_or_less_label)%>%
            mutate(
                labels_colors=case_when(
                    
                    more_or_less_label %in% c('B','L') & color_breks==min(color_breks) ~ paste0('≤',colo_breaks),
                    more_or_less_label %in% c('B','M') & color_breks==max(color_breks) ~ paste0('≥',colo_breaks),
                    T ~ as.character(colo_breaks)
                    
                ))%>%pull(labels_colors)
        # assign color fill 
        matrix = matrix + scale_fill_gradient(
            low = low_color,
            high = high_color,
            limits = limits,
            name = name,
            breaks = colo_breaks,
            labels=labels_colors
        )
        # if no mid color is assign 
    } else{
        # calculate breaks positions
        
        colo_breaks=seq(limits[1], limits[2], sum(abs(limits)) /
                            n_fill_breaks)
        #create breaks labels 
        labels_colors=tibble(color_breks=colo_breaks,more_or_less_label=more_or_less_label)%>%
            mutate(
                labels_colors=case_when(
                    
                    more_or_less_label %in% c('B','L') & color_breks==min(color_breks) ~ paste0('≤',colo_breaks),
                    more_or_less_label %in% c('B','M') & color_breks==max(color_breks) ~ paste0('≥',colo_breaks),
                    T ~ as.character(colo_breaks)
                    
                ))%>%pull(labels_colors)
        # assign color fill 
        matrix = matrix + scale_fill_gradient2(
            low = low_color,
            high = high_color,
            mid = mid_color,
            midpoint = 0,
            limits = limits,
            name = name,
            breaks = colo_breaks,
            labels=labels_colors
        )
        
    }
    
    
    matrix = matrix +
        # add breaks to the plot
        scale_x_continuous(breaks = breaks, labels = labels) +
        #define style
        style +
        theme(
            axis.line.y = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            legend.position = legend.position,
            legend.title = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(
                hjust = 0.5,
                family = font,
                size = tex_size,
                face = fontface,
                color = 'black'
            )
        ) +
        #add rectangle around groups
        geom_rect(
            data = max_min,
            aes(
                xmin = min_x - 0.5,
                ymin = min_y - 0.5,
                ymax = max_y + 0.5,
                xmax = max_x + 0.5
            ),
            size = line_size,
            fill = NA,
            color = 'black'
        ) +
        #add title 
        ggtitle(name) +
        # divide groups
        facet_grid(group ~ ., scales = 'free', space = "free_y") +
        #add gene bodies dashed lines
        geom_segment(
            data = max_min %>%
                mutate(x_breaks = list(breaks[c(2, 3)])) %>%
                unnest(cols =x_breaks),
            aes(
                x = x_breaks ,
                xend = x_breaks,
                y = min_y - 0.5,
                yend = max_y + 0.5
            ),
            linetype = 'dashed',
            color = 'black',
            size =line_size
        )
    #return with final theme touches
    return(
        matrix + theme(
            axis.text.x = element_blank(),
            legend.key.width = unit(0.2, "cm"),
            legend.key.height = unit(0.1, "cm"),
            strip.background = element_blank(),
            strip.text.y = element_blank(),
            panel.spacing=unit(0,'mm')
        )
    )
    
}

#using the info of groups creates labels to put on the side of the heatmaps
side_banner = function(group,
                       group_order=NA,
                       tex_size = 7,
                       font = 'Courier',
                       fontface = 'bold',rotate=90) {
    
    if(is.na(group_order)){
        group_order=unique(group)%>%sort()
    }
    group = tibble(group = factor(group,levels = group_order)) %>%
        group_by(group) %>%
        summarise(n = n()) %>%
        ungroup() %>%
        arrange(group)
    
    vec_textGrob = Vectorize(FUN = textGrob,
                             vectorize.args = 'label',
                             SIMPLIFY = F)
    group_to_plot = vec_textGrob(
        label = group$group,
        rot = rotate,
        hjust = 0.5,
        vjust = 0.5,
        gp = gpar(
            fontsize = tex_size,
            fontface = fontface,
            fontfamily = font
        )
    )
    p = plot_grid(plotlist = group_to_plot,
                  ncol = 1,
                  rel_heights = group$n/sum(group$n))
    
    return(p)
}

#select genes that are bigger than 300 including SDRs
region_300 = genes[genes$split_size == "≥300kb" ]
# assign genes with SDR to a new group
region_300$split_expression[region_300$split=='SDR']='SDR'
# calculate matrices 
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

#lines order based on Groseq
line_order_300 = tibble(m = rowMeans(exp_300[, grepl(pattern = '^t', x = colnames(exp_300))])) %>%
    mutate(line = 1:n()) %>%
    arrange(m) %>% pull(line)


Gro_Heatmap_300=heatmap(
    matrix = exp_300,
    limits = c(0, 500),
    more_or_less_label = 'M',
    line_order = line_order_300,
    name = 'Gro-seq',
    high_color = 'black',
    low_color = 'white',
    group = region_300$split_expression,
    group_order = c( "SDR","1st","2nd", "3rd", "4th" ),
    tex_size = tex_size,
    font = font,
    fontface = fontface,
    n_fill_breaks = 2,
    line_size = line_size,
    legend.position = 'top'
)
RT_Heatmap_300= heatmap(
    matrix = rt_300,
    limits = c(0, 1),
    line_order = line_order_300,
    name = 'RT (NT)',
    high_color = 'red',
    low_color = 'green',
    group = region_300$split_expression,
    group_order = c( "SDR","1st","2nd", "3rd", "4th" ),
    tex_size = tex_size,
    font = font,
    fontface = fontface,
    n_fill_breaks = 3,
    line_size = line_size,
    legend.position = 'top')

Aph_Heatmap_300=heatmap(
    matrix = mat_Aph_300,
    limits = c(-3, 3),
    more_or_less_label = 'B',
    line_order = line_order_300,
    name = 'URI (Aph)',
    high_color = 'blue',
    low_color = 'red',
    mid_color = 'white',
    group = region_300$split_expression,
    group_order = c( "SDR","1st","2nd", "3rd", "4th" ),
    tex_size = tex_size,
    font = font,
    fontface = fontface,
    n_fill_breaks = 3,
    line_size = line_size,
    legend.position = 'top')
ARO_Heatmap_300=heatmap(
    matrix = mat_AphRO_300,
    limits = c(-3, 3),
    more_or_less_label = 'B',
    line_order = line_order_300,
    name = 'URI (ARO)',
    high_color = 'blue',
    low_color = 'red',
    mid_color = 'white',
    group = region_300$split_expression,
    group_order = c( "SDR","1st","2nd", "3rd", "4th" ),
    tex_size = tex_size,
    font = font,
    fontface = fontface,
    n_fill_breaks = 3,
    line_size = line_size,
    legend.position = 'top')
HU_Heatmap_300=heatmap(
    matrix = mat_HU_300,
    limits = c(-3, 3),
    more_or_less_label = 'B',
    line_order = line_order_300,
    name = 'URI (HU)',
    high_color = 'blue',
    low_color = 'red',
    mid_color = 'white',
    group = region_300$split_expression,
    group_order = c( "SDR","1st","2nd", "3rd", "4th" ),
    tex_size = tex_size,
    font = font,
    fontface = fontface,
    n_fill_breaks = 3,
    line_size = line_size,
    legend.position = 'top')


ATRIHU_Heatmap_300=heatmap(
    matrix = mat_ATRiHU_300,
    limits = c(-3, 3),
    more_or_less_label = 'B',
    line_order = line_order_300,
    name = 'URI (HU+ATRi)',
    high_color = 'blue',
    low_color = 'red',
    mid_color = 'white',
    group = region_300$split_expression,
    group_order = c( "SDR","1st","2nd", "3rd", "4th" ),
    tex_size = tex_size,
    font = font,
    fontface = fontface,
    n_fill_breaks = 3,
    line_size = line_size,
    legend.position = 'top')

Heat_maps_300 = plot_grid(
    Gro_Heatmap_300,
    RT_Heatmap_300,
    Aph_Heatmap_300,
    ARO_Heatmap_300,
    HU_Heatmap_300,
    ATRIHU_Heatmap_300,
            side_banner(
                group = region_300$split_expression,group_order = c( "SDR","1st","2nd", "3rd", "4th" ),
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                rotate=90
            ),
        
        nrow = 1,
        rel_widths = c(1.03, 1.05, 1, 1,1,1, 0.3),scale = 0.9
    )
#save
ggsave(
    plot =Gro_Heatmap_300,
    filename = paste0(output_dir,'/Gro_Heatmap_300_FigS3A.pdf'),
    width = 3,
    height = 17,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)
ggsave(
    plot =RT_Heatmap_300,
    filename = paste0(output_dir,'/RT_Heatmap_300_FigS3A.pdf'),
    width = 3,
    height = 17,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)
ggsave(
    plot =Aph_Heatmap_300,
    filename = paste0(output_dir,'/Aph_Heatmap_300_FigS3A.pdf'),
    width = 3,
    height = 17,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)
ggsave(
    plot =ARO_Heatmap_300,
    filename = paste0(output_dir,'/ARO_Heatmap_300_FigS3A.pdf'),
    width = 3,
    height = 17,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)
ggsave(
    plot =HU_Heatmap_300,
    filename = paste0(output_dir,'/HU_Heatmap_300_FigS3A.pdf'),
    width = 3,
    height = 17,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)

ggsave(
    plot =ATRIHU_Heatmap_300,
    filename = paste0(output_dir,'/ATRIHU_Heatmap_300_FigS3A.pdf'),
    width = 3,
    height = 17,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)
ggsave(
    plot =Heat_maps_300,
    filename = paste0(output_dir,'/FigS3A.pdf'),
    width = 17,
    height = 17,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)
