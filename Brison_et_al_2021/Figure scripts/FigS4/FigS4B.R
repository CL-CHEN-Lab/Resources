#########################################################
#########################################################

# RUN Fig4A.R, Fig4B, Fig4C and FigS4A first


#########################################################
#########################################################
output_dir='~/Desktop/Figure_paper_to_update/FigS4/FigS4B'
system(paste0('mkdir -p ',output_dir))

#250
region_250 = genes[genes$split_size == "≥250kb"]

rt_250 = normalizeToMatrix(
    target = region_250,
    signal =  NT,
    value_column = 's50',
    mean_mode = "w0",
    extend = 100000,
    target_ratio = 0.45
)
exp_250 = normalizeToMatrix(
    target = region_250,
    signal =  bins,
    value_column = 'counts',
    mean_mode = "w0",
    extend = 100000,
    target_ratio = 0.45
)
mat_Aph_250 = normalizeToMatrix(
    target = region_250,
    signal =  URI,
    value_column = 'Aph',
    extend = 100000,
    mean_mode = "w0",
    target_ratio = 0.45
)

#lines order based on Groseq
line_order_250 = tibble(m = rowMeans(exp_250[, grepl(pattern = '^t', x = colnames(exp_250))])) %>%
    mutate(line = 1:n()) %>%
    arrange(m) %>% pull(line)


Heat_maps_250 = plot_grid(
    
    plot_grid(
        plot_grid(
            
            heatmap(
                matrix = exp_250,
                limits = c(0, 500),
                more_or_less_label = 'M',
                line_order = line_order_250,
                name = 'Gro-seq',
                high_color = 'black',
                low_color = 'white',
                group = region_250$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                n_fill_breaks = 2,
                line_size = line_size,
                legend.position = 'top'),
            average_profile(
                exp_250,
                name = 'Gro-seq',
                group = region_250$split_expression,
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            )+scale_y_continuous(breaks=c(0,1000,2000),labels = c(0,'1K','2K'))+
                theme(strip.text  = element_blank()),
            ncol = 1,
            rel_heights = c(3,1),
            align = 'v',
            axis = "l"
        )  ,
        
        plot_grid(
            
            heatmap(
                matrix = rt_250,
                limits = c(0, 1),
                line_order = line_order_250,
                name = 'RT (NT)',
                high_color = 'red',
                low_color = 'green',
                group = region_250$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                n_fill_breaks = 3,
                line_size = line_size,
                legend.position = 'top'),
            average_profile(
                rt_250,
                name = 'RT (NT)',
                group = region_250$split_expression,
                limits = c(1, 0),
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            ) + scale_y_continuous(breaks=c(0,0.5,1))+
                theme(strip.text  = element_blank()),
            ncol = 1,
            rel_heights = c(3,1),
            align = 'v',
            axis = "l"
        ),
        
        
        plot_grid(
            
            heatmap(
                matrix = mat_Aph_250,
                limits = c(-3, 3),
                more_or_less_label = 'B',
                line_order = line_order_250,
                name = 'URI (Aph)',
                high_color = 'blue',
                low_color = 'red',
                mid_color = 'white',
                group = region_250$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                n_fill_breaks = 3,
                line_size = line_size,
                legend.position = 'top'),
            average_profile(
                mat_Aph_250,
                name = 'URI (Aph)',
                group = region_250$split_expression,
                limits = c(-3, 3),
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            )+
                theme(strip.text  = element_blank()),
            ncol = 1,
            rel_heights = c(3,1),
            align = 'v',
            axis = "l"
        ),
        
        plot_grid(
            textGrob(''),
            side_banner(
                group = region_250$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                rotate=90
            ),
            textGrob(''),
            ncol = 1,
            rel_heights = c(0.45, 2, 0.85)
        ),
        
        nrow = 1,
        rel_widths = c(1.03, 1.05, 1, 0.3)
    ),
    ncol = 1
)

#200
region_200 = genes[genes$split_size == "≥200kb"]

rt_200 = normalizeToMatrix(
    target = region_200,
    signal =  NT,
    value_column = 's50',
    mean_mode = "w0",
    extend = 100000,
    target_ratio = 0.45
)
exp_200 = normalizeToMatrix(
    target = region_200,
    signal =  bins,
    value_column = 'counts',
    mean_mode = "w0",
    extend = 100000,
    target_ratio = 0.45
)
mat_Aph_200 = normalizeToMatrix(
    target = region_200,
    signal =  URI,
    value_column = 'Aph',
    extend = 100000,
    mean_mode = "w0",
    target_ratio = 0.45
)


#lines order based on Groseq
line_order_200 = tibble(m = rowMeans(exp_200[, grepl(pattern = '^t', x = colnames(exp_200))])) %>%
    mutate(line = 1:n()) %>%
    arrange(m) %>% pull(line)


Heat_maps_200 = plot_grid(
    
    plot_grid(
        plot_grid(
            
            heatmap(
                matrix = exp_200,
                limits = c(0, 500),
                more_or_less_label = 'M',
                line_order = line_order_200,
                name = 'Gro-seq',
                high_color = 'black',
                low_color = 'white',
                group = region_200$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                n_fill_breaks = 2,
                line_size = line_size,
                legend.position = 'top'),
            average_profile(
                exp_200,
                name = 'Gro-seq',
                group = region_200$split_expression,
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            )+scale_y_continuous(breaks=c(0,1000,2000),labels = c(0,'1K','2K'))+
                theme(strip.text  = element_blank()),
            ncol = 1,
            rel_heights = c(3,1),
            align = 'v',
            axis = "l"
        )  ,
        
        plot_grid(
            
            heatmap(
                matrix = rt_200,
                limits = c(0, 1),
                line_order = line_order_200,
                name = 'RT (NT)',
                high_color = 'red',
                low_color = 'green',
                group = region_200$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                n_fill_breaks = 3,
                line_size = line_size,
                legend.position = 'top'),
            average_profile(
                rt_200,
                name = 'RT (NT)',
                group = region_200$split_expression,
                limits = c(1, 0),
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            ) + scale_y_continuous(breaks=c(0,0.5,1))+
                theme(strip.text  = element_blank()),
            ncol = 1,
            rel_heights = c( 3,1),
            align = 'v',
            axis = "l"
        ),
        
        plot_grid(
            
            heatmap(
                matrix = mat_Aph_200,
                limits = c(-3, 3),
                more_or_less_label = 'B',
                line_order = line_order_200,
                name = 'URI (Aph)',
                high_color = 'blue',
                low_color = 'red',
                mid_color = 'white',
                group = region_200$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                n_fill_breaks = 3,
                line_size = line_size,
                legend.position = 'top'),
            average_profile(
                mat_Aph_200,
                name = 'URI (Aph)',
                group = region_200$split_expression,
                limits = c(-3, 3),
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            )+
                theme(strip.text  = element_blank()),
            ncol = 1,
            rel_heights = c(3,1),
            align = 'v',
            axis = "l"
        ),
        
        plot_grid(
            textGrob(''),
            side_banner(
                group = region_200$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                rotate=90
            ),
            textGrob(''),
            ncol = 1,
            rel_heights = c(0.45, 2, 0.85)
        ),
        nrow = 1,
        rel_widths = c(1.03, 1.05, 1, 0.3)
    ), 

    ncol = 1
)


#150
region_150 = genes[genes$split_size == "≥150kb"]

rt_150 = normalizeToMatrix(
    target = region_150,
    signal =  NT,
    value_column = 's50',
    mean_mode = "w0",
    extend = 100000,
    target_ratio = 0.45
)
exp_150 = normalizeToMatrix(
    target = region_150,
    signal =  bins,
    value_column = 'counts',
    mean_mode = "w0",
    extend = 100000,
    target_ratio = 0.45
)
mat_Aph_150 = normalizeToMatrix(
    target = region_150,
    signal =  URI,
    value_column = 'Aph',
    extend = 100000,
    mean_mode = "w0",
    target_ratio = 0.45
)

#lines order based on Groseq
line_order_150 = tibble(m = rowMeans(exp_150[, grepl(pattern = '^t', x = colnames(exp_150))])) %>%
    mutate(line = 1:n()) %>%
    arrange(m) %>% pull(line)


Heat_maps_150 = plot_grid(
    
    plot_grid(
        plot_grid(
            
            heatmap(
                matrix = exp_150,
                limits = c(0, 500),
                more_or_less_label = 'M',
                line_order = line_order_150,
                name = 'Gro-seq',
                high_color = 'black',
                low_color = 'white',
                group = region_150$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                n_fill_breaks = 2,
                line_size = line_size,
                legend.position = 'top'),
            average_profile(
                exp_150,
                name = 'Gro-seq',
                group = region_150$split_expression,
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            )+scale_y_continuous(breaks=c(0,1000,2000),labels = c(0,'1K','2K'))+
                theme(strip.text  = element_blank()),
            ncol = 1,
            rel_heights = c(3,1),
            align = 'v',
            axis = "l"
        )  ,
        
        plot_grid(
            
            heatmap(
                matrix = rt_150,
                limits = c(0, 1),
                line_order = line_order_150,
                name = 'RT (NT)',
                high_color = 'red',
                low_color = 'green',
                group = region_150$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                n_fill_breaks = 3,
                line_size = line_size,
                legend.position = 'top'),
            average_profile(
                rt_150,
                name = 'RT (NT)',
                group = region_150$split_expression,
                limits = c(1, 0),
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            ) + scale_y_continuous(breaks=c(0,0.5,1))+
                theme(strip.text  = element_blank()),
            ncol = 1,
            rel_heights = c(3,1),
            align = 'v',
            axis = "l"
        ),
        
        plot_grid(
            
            heatmap(
                matrix = mat_Aph_150,
                limits = c(-3, 3),
                more_or_less_label = 'B',
                line_order = line_order_150,
                name = 'URI (Aph)',
                high_color = 'blue',
                low_color = 'red',
                mid_color = 'white',
                group = region_150$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                n_fill_breaks = 3,
                line_size = line_size,
                legend.position = 'top'),
            average_profile(
                mat_Aph_150,
                name = 'URI (Aph)',
                group = region_150$split_expression,
                limits = c(-3, 3),
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            )+
                theme(strip.text  = element_blank()),
            ncol = 1,
            rel_heights = c(3,1),
            align = 'v',
            axis = "l"
        ),
        
        plot_grid(
            textGrob(''),
            side_banner(
                group = region_150$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                rotate=90
            ),
            textGrob(''),
            ncol = 1,
            rel_heights = c(0.45, 2, 0.85)
        ),
        
        nrow = 1,
        rel_widths = c(1.03, 1.05, 1, 0.3)
    ),

    ncol = 1
)

#100
region_100 = genes[genes$split_size == "≥100kb"]

rt_100 = normalizeToMatrix(
    target = region_100,
    signal =  NT,
    value_column = 's50',
    mean_mode = "w0",
    extend = 100000,
    target_ratio = 0.45
)
exp_100 = normalizeToMatrix(
    target = region_100,
    signal =  bins,
    value_column = 'counts',
    mean_mode = "w0",
    extend = 100000,
    target_ratio = 0.45
)
mat_Aph_100 = normalizeToMatrix(
    target = region_100,
    signal =  URI,
    value_column = 'Aph',
    extend = 100000,
    mean_mode = "w0",
    target_ratio = 0.45
)


#lines order based on Groseq
line_order_100 = tibble(m = rowMeans(exp_100[, grepl(pattern = '^t', x = colnames(exp_100))])) %>%
    mutate(line = 1:n()) %>%
    arrange(m) %>% pull(line)


Heat_maps_100 = plot_grid(
    plot_grid(
        plot_grid(
            heatmap(
                matrix = exp_100,
                limits = c(0, 500),
                more_or_less_label = 'M',
                line_order = line_order_100,
                name = 'Gro-seq',
                high_color = 'black',
                low_color = 'white',
                group = region_100$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                n_fill_breaks = 2,
                line_size = line_size,
                legend.position = 'top'),
            average_profile(
                exp_100,
                name = 'Gro-seq',
                group = region_100$split_expression,
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            )+scale_y_continuous(breaks=c(0,1000,2000),labels = c(0,'1K','2K'))+
                theme(strip.text  = element_blank()),
            ncol = 1,
            rel_heights = c(3,1),
            align = 'v',
            axis = "l"
        )  ,
        
        plot_grid(
            
            heatmap(
                matrix = rt_100,
                limits = c(0, 1),
                line_order = line_order_100,
                name = 'RT (NT)',
                high_color = 'red',
                low_color = 'green',
                group = region_100$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                n_fill_breaks = 3,
                line_size = line_size,
                legend.position = 'top'),
            average_profile(
                rt_100,
                name = 'RT (NT)',
                group = region_100$split_expression,
                limits = c(1, 0),
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            ) + scale_y_continuous(breaks=c(0,0.5,1))+
                theme(strip.text  = element_blank()),
            ncol = 1,
            rel_heights = c(3,1),
            align = 'v',
            axis = "l"
        ),
        
        plot_grid(
            heatmap(
                matrix = mat_Aph_100,
                limits = c(-3, 3),
                more_or_less_label = 'B',
                line_order = line_order_100,
                name = 'URI (Aph)',
                high_color = 'blue',
                low_color = 'red',
                mid_color = 'white',
                group = region_100$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                n_fill_breaks = 3,
                line_size = line_size,
                legend.position = 'top'),
            average_profile(
                mat_Aph_100,
                name = 'URI (Aph)',
                group = region_100$split_expression,
                limits = c(-3, 3),
                legend_position = 'none',
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                line_size = line_size,
                Colors = Colors
            )+
                theme(strip.text  = element_blank()),
            ncol = 1,
            rel_heights = c(3,1),
            align = 'v',
            axis = "l"
        ),
        
        plot_grid(
            textGrob(''),
            side_banner(
                group = region_100$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                rotate=90
            ),
            textGrob(''),
            ncol = 1,
            rel_heights = c(0.45, 2, 0.85)
        ),
        
        nrow = 1,
        rel_widths = c(1.03, 1.05, 1, 0.3)
    ),

    ncol = 1
)

# assemble pic and save
ggsave(
    plot = plot_grid(
        textGrob(''),
        plot_grid(
            plot_grid(
                Heat_maps_250,
                textGrob('300Kb > genes ≥ 250kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=-0.50),
                
                rel_widths = c(1,0.1)),
            plot_grid(
                Heat_maps_200,
                textGrob('250Kb > genes ≥ 200kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            plot_grid(        
                Heat_maps_150,
                textGrob('200Kb > genes ≥ 150kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            plot_grid(
                Heat_maps_100,
                textGrob('150Kb > genes ≥ 100kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            ncol = 2,
            nrow = 2,
            labels = c("G", "H", "I", "J"),
            label_size = label_text_size,
            label_fontfamily = font
        ),
        get_legend(
            average_profile(
                exp_100,
                name = 'Gro-seq',
                group = region_100$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                Colors = Colors
            )
        ),rel_heights = c(0.03,1,0.03), ncol = 1),
    filename = paste0(output_dir,'/FigS4B.pdf'),
    width = 174,
    height = 210,
    limitsize = F,
    units = 'mm',
    device = cairo_pdf
)

#save single
ggsave(
    plot = plot_grid(
        textGrob(''),
        plot_grid(
            plot_grid(
                Heat_maps_250,
                textGrob('300Kb > genes ≥ 250kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=-0.50),
                
                rel_widths = c(1,0.1)),
            plot_grid(
                textGrob(''),
                textGrob('250Kb > genes ≥ 200kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            plot_grid(        
                textGrob(''),
                textGrob('200Kb > genes ≥ 150kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            plot_grid(
                textGrob(''),
                textGrob('150Kb > genes ≥ 100kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            ncol = 2,
            nrow = 2,
            labels = c("G", "H", "I", "J"),
            label_size = label_text_size,
            label_fontfamily = font
        ),
        get_legend(
            average_profile(
                exp_100,
                name = 'Gro-seq',
                group = region_100$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                Colors = Colors
            )
        ),rel_heights = c(0.03,1,0.03), ncol = 1),
    filename = paste0(output_dir,'/250_FigS4B.pdf'),
    width = 174,
    height = 210,
    limitsize = F,
    units = 'mm',
    device = cairo_pdf
)

ggsave(
    plot = plot_grid(
        textGrob(''),
        plot_grid(
            plot_grid(
                textGrob(''),
                textGrob('300Kb > genes ≥ 250kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=-0.50),
                
                rel_widths = c(1,0.1)),
            plot_grid(
                Heat_maps_200,
                textGrob('250Kb > genes ≥ 200kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            plot_grid(        
                textGrob(''),
                textGrob('200Kb > genes ≥ 150kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            plot_grid(
                textGrob(''),
                textGrob('150Kb > genes ≥ 100kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            ncol = 2,
            nrow = 2,
            labels = c("G", "H", "I", "J"),
            label_size = label_text_size,
            label_fontfamily = font
        ),
        get_legend(
            average_profile(
                exp_100,
                name = 'Gro-seq',
                group = region_100$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                Colors = Colors
            )
        ),rel_heights = c(0.03,1,0.03), ncol = 1),
    filename = paste0(output_dir,'/200_FigS4B.pdf'),
    width = 174,
    height = 210,
    limitsize = F,
    units = 'mm',
    device = cairo_pdf
)


ggsave(
    plot = plot_grid(
        textGrob(''),
        plot_grid(
            plot_grid(
                textGrob(''),
                textGrob('300Kb > genes ≥ 250kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=-0.50),
                
                rel_widths = c(1,0.1)),
            plot_grid(
                textGrob(''),
                textGrob('250Kb > genes ≥ 200kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            plot_grid(        
                Heat_maps_150,
                textGrob('200Kb > genes ≥ 150kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            plot_grid(
                textGrob(''),
                textGrob('150Kb > genes ≥ 100kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            ncol = 2,
            nrow = 2,
            labels = c("G", "H", "I", "J"),
            label_size = label_text_size,
            label_fontfamily = font
        ),
        get_legend(
            average_profile(
                exp_100,
                name = 'Gro-seq',
                group = region_100$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                Colors = Colors
            )
        ),rel_heights = c(0.03,1,0.03), ncol = 1),
    filename = paste0(output_dir,'/150_FigS4B.pdf'),
    width = 174,
    height = 210,
    limitsize = F,
    units = 'mm',
    device = cairo_pdf
)

ggsave(
    plot = plot_grid(
        textGrob(''),
        plot_grid(
            plot_grid(
                textGrob(''),
                textGrob('300Kb > genes ≥ 250kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=-0.50),
                
                rel_widths = c(1,0.1)),
            plot_grid(
                textGrob(''),
                textGrob('250Kb > genes ≥ 200kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            plot_grid(        
                textGrob(''),
                textGrob('200Kb > genes ≥ 150kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            plot_grid(
                Heat_maps_100,
                textGrob('150Kb > genes ≥ 100kb',gp = gpar(
                    fontsize = tex_size,
                    fontface = fontface,
                    fontfamily = font
                ),
                rot = 90,
                vjust=0),
                rel_widths = c(1,0.1)),
            ncol = 2,
            nrow = 2,
            labels = c("G", "H", "I", "J"),
            label_size = label_text_size,
            label_fontfamily = font
        ),
        get_legend(
            average_profile(
                exp_100,
                name = 'Gro-seq',
                group = region_100$split_expression,
                tex_size = tex_size,
                font = font,
                fontface = fontface,
                Colors = Colors
            )
        ),rel_heights = c(0.03,1,0.03), ncol = 1),
    filename = paste0(output_dir,'/100_FigS4B.pdf'),
    width = 174,
    height = 210,
    limitsize = F,
    units = 'mm',
    device = cairo_pdf
)

