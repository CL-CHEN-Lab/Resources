output_dir='~/Desktop/Figure_paper_to_update/Fig4/Fig4A'
system(paste0('mkdir -p ',output_dir))

library(extrafont)
loadfonts(device = "postscript")
library(cowplot)
library(tidyverse)
library(GenomicRanges)
library(EnrichedHeatmap)

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
# load genes
genes = read_tsv('./Gro-Seq/refseq_genes_h19') %>%
    # keep 1:22 and X
    filter(chrom %in% paste0('chr', c(1:22, 'X'))) %>%
    # convert to Granges
    makeGRangesFromDataFrame(
        seqnames.field = 'chrom',
        start.field = 'txStart',
        end.field = 'txEnd',
        strand.field = 'strand',
        keep.extra.columns = T
    )

# split isoforms into differen lists and reduce 
genes = split(genes, genes$name2) %>% GenomicRanges::reduce() %>% unlist()

# keep genes bigger than 100kb
genes = genes[width(genes) >= 100000]

#assign names
genes$name2 = genes@ranges@NAMES

# load groseq data
Gro = rbind(
    read_tsv(
        './Gro-Seq/GM12878_GROseq_plus.bedGraph',
        col_names = c('chr', 'start', 'end', 'counts')
    ) %>% mutate(strand = '+') ,
    read_tsv(
        './Gro-Seq/GM12878_GROseq_minus.bedGraph',
        col_names = c('chr', 'start', 'end', 'counts')
    ) %>% mutate(strand = '-',
                 counts = -counts)
) %>% filter(!chr %in% c('chrY', 'chrM')) %>%
    makeGRangesFromDataFrame(
        seqnames.field = 'chr',
        start.field = 'start',
        end.field = 'end',
        strand = 'strand',
        keep.extra.columns = T
    )

# from S50 recove NT rt and convert to GRanges
NT = S50 %>%
    #exclude 0
    filter(Condition == 'NT',
           !is.na(s50),
           s50 != 0) %>%
    makeGRangesFromDataFrame(
        seqnames.field = 'chr',
        start.field = 'start',
        end.field = 'end',
        keep.extra.columns = T
    )

library(BSgenome.Hsapiens.UCSC.hg19)

# breat a 10kb genome binning
bins = tileGenome(
    seqlengths = seqinfo(Hsapiens),
    tilewidth = 10000,
    cut.last.tile.in.chrom = T
)
# overlap between Groseq data and bins
hits = findOverlaps(query = Gro, subject = bins)
overlaps <-
    pintersect(Gro[queryHits(hits)], bins[subjectHits(hits)])

bins_m = bins[subjectHits(hits)] %>%
    as_tibble() %>%
    # percentage of count overlapping with the bin
    mutate(counts = Gro$counts[queryHits(hits)] * width(overlaps) / width(Gro[queryHits(hits)])) %>%
    group_by(seqnames, start, end, width) %>%
    # sum counts in bin
    summarise(counts = sum(counts))

#add counts to bins
bins = bins %>%
    as_tibble() %>%
    left_join(bins_m, by = c("seqnames", "start", "end", "width")) %>%
    mutate(counts = ifelse(is.na(counts), 0, counts)) %>% makeGRangesFromDataFrame(keep.extra.columns = T)

# overlap with genes
hits = findOverlaps(query = bins, subject = genes)
overlaps <-
    pintersect(bins[queryHits(hits)], genes[subjectHits(hits)])

# genes level of expression are a weighted median 
genes = genes[subjectHits(hits)] %>%
    as_tibble() %>%
    mutate(counts = bins$counts[queryHits(hits)],
           W = width(overlaps)) %>%
    group_by(seqnames, start, end, width, strand, name2) %>%
    summarise(counts = median(rep(counts, W))) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)

#assign RT
hits = findOverlaps(query = NT, subject = genes)
overlaps <-
    pintersect(NT[queryHits(hits)], genes[subjectHits(hits)])

# RT is the weighted median RT
genes = genes[subjectHits(hits)] %>%
    as_tibble() %>%
    mutate(RT = NT$s50[queryHits(hits)],
           W = width(overlaps)) %>%
    group_by(seqnames, start, end, strand, width, counts, name2) %>%
    summarise(RT = median(rep(RT, W))) %>%
    ungroup() 

# calculate quintiles once excluced conts ==0
q = quantile(genes$counts[genes$counts != 0], c(.33, .66))

#assign genes to differen expression groups and size 
genes = genes %>% dplyr::mutate(
    f = 0,
    s = q[[1]],
    t = q[[2]],
    # assign size group
    split_size=case_when(
        width >= 300000 ~ '≥300kb',
        width >= 250000 ~ '≥250kb',
        width >= 200000 ~ '≥200kb',
        width >= 150000 ~ '≥150kb',
        width >= 100000 ~ '≥100kb'
    ),
    #assign expression
    split_expression= case_when(
        counts <= f ~ '4th',
        counts > f & counts <= s ~ '3rd',
        counts > s & counts <= t  ~ '2nd',
        counts > t  ~ '1st'
    )
) %>%
    # filter out chrY and chrM if any and make Granges 
    dplyr::filter(!seqnames %in% c('chrY', 'chrM')) %>% makeGRangesFromDataFrame(keep.extra.columns = T)

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
    mutate(DR = ifelse(is.na(DR), 'SDR', DR))%>%
    makeGRangesFromDataFrame(keep.extra.columns = T)

# overlap with DR
hits = findOverlaps(query = genes, subject = DR)
overlaps <-
    pintersect(genes[queryHits(hits)], DR[subjectHits(hits)])
# check that the overlap of either DR or Genes is at least 25%
overlaps_1 = width(overlaps) / width(DR[subjectHits(hits)]) >= 0.25
overlaps_2 = width(overlaps) / width(genes[queryHits(hits)]) >= 0.25

#initialise a column
genes$split='NonSDR'
# select gens in which either feature overlaps for at least 25% and chagen NonSDR in SDR
genes[queryHits(hits)[overlaps_1|overlaps_2]]$split = 'SDR'

#Load and convert URI into grange and associate infos
URI = read_tsv('URi_simulation_periodic_1kbXY_ATRIadd.tsv')%>%
    filter(
        Condition %in% c("NT","Aph","AphRO","HU","ATRiHU")
    ) %>% inner_join(S50 %>% filter(Condition == 'NT') %>% dplyr::select(-Condition)) %>%
    mutate(s50 = ifelse(s50 < 0.5, 'Early', 'Late'))%>%
    spread(Condition, URI) %>% makeGRangesFromDataFrame(keep.extra.columns = T)

# select Genes that are bigger than 300kb and overlap with SDR
region_small = genes[genes$split == 'SDR' & genes$split_size =="≥300kb"]

# calculate matrixes for plotting of each URI signal as well as Gro-seq and RT (NT)
rt_sdr = normalizeToMatrix(
    target = region_small,
    signal =  NT,
    value_column = 's50',
    mean_mode = "w0",
    extend = 100000,
    target_ratio = 0.45
)
exp_sdr = normalizeToMatrix(
    target = region_small,
    signal =  bins,
    value_column = 'counts',
    mean_mode = "w0",
    extend = 100000,
    target_ratio = 0.45
)
mat_Aph_sdr = normalizeToMatrix(
    target = region_small,
    signal =  URI,
    value_column = 'Aph',
    extend = 100000,
    mean_mode = "w0",
    target_ratio = 0.45
)
mat_AphRO_sdr = normalizeToMatrix(
    target = region_small,
    signal =  URI,
    value_column = 'AphRO',
    extend = 100000,
    mean_mode = "w0",
    target_ratio = 0.45
)
mat_HU_sdr = normalizeToMatrix(
    target = region_small,
    signal =  URI,
    value_column = 'HU',
    extend = 100000,
    mean_mode = "w0",
    target_ratio = 0.45
)

mat_ATRiHU_sdr = normalizeToMatrix(
    target = region_small,
    signal =  URI,
    value_column = 'ATRiHU',
    extend = 100000,
    mean_mode = "w0",
    target_ratio = 0.45
)

#lines order based on Groseq
line_order_small = tibble(m = rowMeans(exp_sdr[, grepl(pattern = '^t', x = colnames(exp_sdr))])) %>%
    mutate(line = 1:n()) %>%
    arrange(m) %>% pull(line)

heatmap_small = function(matrix,
                         line_order = NA,
                         name = heatmap,
                         limits = c(-Inf, Inf),
                         low_color = 'red',
                         high_color = 'blue',
                         mid_color = NA,
                         legend.position = 'bottom',
                         breaks = 'auto',
                         n_fill_breaks = 5,
                         more_or_less_label=c('N','B','M','L'),
                         labels = c('-100', 'TSS', 'TTS', '100'),
                         font = 30,
                         tex_size = 'Courier',
                         fontface = 'bold',
                         style = general_theme,
                         gene_list = NA,
                         show_name = T,
                         line_size=1.5) {
    more_or_less_label=more_or_less_label[1]
    
    #if the used does not provide info about the x axis breaks it looks  for the values into the matrix
    if (breaks == 'auto') {
        upstream = sum(grepl(x = colnames(matrix), pattern = 'u'))
        body = sum(grepl(x = colnames(matrix), pattern = 't'))
        dowstream = sum(grepl(x = colnames(matrix), pattern = 'd'))
        
        breaks = c(0, upstream + 1, upstream + body, upstream + body + dowstream) %>%
            unique()
    }
    
    # convert matrix ito a tibble 
    matrix = matrix %>% as_tibble() %>%
        mutate(line = 1:n())
    
    # if the names of genes in each line are provided it adds a colum
    if (!any(is.na(gene_list))) {
        matrix = matrix  %>%
            mutate(genes = gene_list)
    }
    
    # if a specifc order is request, it convers line into a factor 
    # and then it assigns a new line index
    if (!any(is.na(line_order))) {
        matrix = matrix  %>%
            mutate(line = factor(line, levels = line_order),
                   line = as.numeric(line)) %>%
            arrange(line) %>%
            mutate(line = 1:n())
    }
    
    # gene list case
    if (!any(is.na(gene_list))) {
        #reformat data for plotting 
        matrix = matrix  %>%
            mutate(genes = factor(genes, levels = genes)) %>%
            gather(position, value, -line, -genes) %>%
            mutate(
                position = factor(position, levels = colnames(matrix)),
                position = as.numeric(position),
                value = as.numeric(value)
            )
        #non gene list case
    } else{
        #reformat data for plotting 
        matrix = matrix %>%
            gather(position, value, -line) %>%
            mutate(
                position = factor(position, levels = colnames(matrix)),
                position = as.numeric(position),
                value = as.numeric(value)
            )
    }
    
    # summarise limits
    max_min = matrix %>% summarise(
        min_x = min(position),
        max_x = max(position),
        min_y = min(line),
        max_y = max(line)
    )
    
    
    #gene list case 
    if (!any(is.na(gene_list))) {
        # plot 
        matrix = ggplot(matrix) +
            # tiles genes vs X
            geom_tile(aes(
                y = genes,
                x = position,
                # if the fill exide min or max values level them to the manual imposed ones
                fill = ifelse(
                    value < min(limits),
                    min(limits),
                    ifelse(value > max(limits), max(limits),
                           value)
                )
            )) +
            # plot labels and breaks
            scale_x_continuous(breaks = breaks, labels = labels) +
            # style
            style +
            theme(
                axis.line.y = element_blank(),
                axis.title = element_blank(),
                legend.position = legend.position,
                legend.title = element_blank(),
                axis.ticks.y = element_blank(),
                plot.title = element_text(
                    hjust = 0.5,
                    family = font,
                    size = tex_size,
                    face = fontface,
                    color = 'black'
                ),
                axis.text.x = element_blank(),
                legend.key.width = unit(0.2, "cm"),
                legend.key.height = unit(0.1, "cm")
            ) +
            # drow box
            geom_rect(
                data = max_min,
                aes(
                    xmin = min_x,
                    ymin = min_y - 0.5,
                    ymax = max_y + 0.5,
                    xmax = max_x
                ),
                size = line_size,
                fill = NA,
                color = 'black'
            ) + 
            #add title
            ggtitle(name) +
            # draw dashed lines 
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
        
        if (!show_name) {
            # if names are give and we want to hide them 
            matrix = matrix + theme(axis.text.y = element_blank())
        }
        # no names case
    } else{
        #plot
        matrix = ggplot(matrix) +
            # plot line (y) vs matrix signal
            geom_tile(aes(
                y = line,
                x = position,
                #saturate fill if outside of user limits
                fill = ifelse(
                    value < min(limits),
                    min(limits),
                    ifelse(value > max(limits), max(limits),
                           value)
                )
            )) +
            # assign breaks and labels 
            scale_x_continuous(breaks = breaks, labels = labels) +
            #add style
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
                ),
                axis.text.x = element_blank(),
                legend.text = element_text(
                    angle = 45,
                    hjust = 1,
                    vjust = 1,
                    family = font,
                    size = tex_size,
                    face = fontface,
                    color = 'black'
                ),
                legend.key.width = unit(0.2, "cm"),
                legend.key.height = unit(0.1, "cm")
            ) +
            # add box around plot
            geom_rect(
                data = max_min,
                aes(
                    xmin = min_x,
                    ymin = min_y - 0.5,
                    ymax = max_y + 0.5,
                    xmax = max_x
                ),
                size = line_size,
                fill = NA,
                color = 'black'
            ) + 
            #add title
            ggtitle(name) +
            #add dashed segments aligned with breaks
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
    }
    
    # number of breaks in the color legend 
    n_fill_breaks = n_fill_breaks - 1
    
    # mid color is no provided 
    if (is.na(mid_color)) {
        
        #find color positions to place a reference number 
        colo_breaks=seq(limits[1], limits[2], sum(abs(limits)) /
                            n_fill_breaks)
        # design the label 
        labels_colors=tibble(color_breks=colo_breaks,more_or_less_label=more_or_less_label)%>%
            mutate(
                labels_colors=case_when(
                    
                    more_or_less_label %in% c('B','L') & color_breks==min(color_breks) ~ paste0('≤',colo_breaks),
                    more_or_less_label %in% c('B','M') & color_breks==max(color_breks) ~ paste0('≥',colo_breaks),
                    T ~ as.character(colo_breaks)
                    
                ))%>%pull(labels_colors)
        # add legend
        matrix = matrix + scale_fill_gradient(
            low = low_color,
            high = high_color,
            limits = limits,
            name = name,
            breaks = colo_breaks,
            labels=labels_colors
        )
        
        # mid color is provided
        } else{
        
        #identify color breaks 
        colo_breaks=seq(limits[1], limits[2], sum(abs(limits)) /
                            n_fill_breaks)
        #create labels 
        labels_colors=tibble(color_breks=colo_breaks,more_or_less_label=more_or_less_label)%>%
            mutate(
                labels_colors=case_when(
                    
                    more_or_less_label %in% c('B','L') & color_breks==min(color_breks) ~ paste0('≤',colo_breaks),
                    more_or_less_label %in% c('B','M') & color_breks==max(color_breks) ~ paste0('≥',colo_breaks),
                    T ~ as.character(colo_breaks)
                    
                ))%>%pull(labels_colors)
        # add scale 
        matrix = matrix + scale_fill_gradient2(
            low = low_color,
            mid = mid_color,
            high = high_color,
            midpoint = 0,
            limits = limits,
            name = name,
            breaks = colo_breaks,
            labels=labels_colors, oob=squish
        )
        
    }
    # return plot 
    return(matrix)
}

#assemble plot
Heat_maps_small = plot_grid(
        heatmap_small(
            matrix = exp_sdr,
            limits = c(0, 500),
            more_or_less_label = 'M',
            line_order = line_order_small,
            name = 'Gro-seq',
            high_color = 'black',
            low_color = 'white',
            gene_list = region_small$name2,
            show_name = T,
            tex_size = tex_size,
            font = font,
            fontface = fontface,
            n_fill_breaks = 2,
            line_size = line_size,
            legend.position = 'top'
        ) ,
    
        heatmap_small(
            matrix = rt_sdr,
            limits = c(0, 1),
            line_order = line_order_small,
            name = 'RT (NT)',
            high_color = 'red',
            low_color = 'green',
            gene_list = region_small$name2,
            show_name = F,
            tex_size = tex_size,
            font = font,
            fontface = fontface,
            n_fill_breaks = 3,
            line_size = line_size,
            legend.position = 'top'
        ),
    
        heatmap_small(
            matrix = mat_Aph_sdr,
            limits = c(-3, 3),
            more_or_less_label = 'B',
            line_order = line_order_small,
            name = 'URI (Aph)',
            high_color = 'blue',
            low_color = 'red',
            mid_color = 'white',
            gene_list = region_small$name2,
            show_name = F,
            tex_size = tex_size,
            font = font,
            fontface = fontface,
            n_fill_breaks = 3,
            line_size = line_size,
            legend.position = 'top'
        ) ,
    
        heatmap_small(
            matrix = mat_AphRO_sdr,
            limits = c(-3, 3),
            more_or_less_label = 'B',
            line_order = line_order_small,
            name = 'URI (ARO)',
            high_color = 'blue',
            low_color = 'red',
            mid_color = 'white',
            gene_list = region_small$name2,
            show_name = F,
            tex_size = tex_size,
            font = font,
            fontface = fontface,
            n_fill_breaks = 3,
            line_size = line_size,
            legend.position = 'top'
        ) ,
    
        heatmap_small(
            matrix = mat_HU_sdr,
            limits = c(-3, 3),
            more_or_less_label = 'B',
            line_order = line_order_small,
            name = 'URI (HU)',
            high_color = 'blue',
            low_color = 'red',
            mid_color = 'white',
            gene_list = region_small$name2,
            show_name = F,
            tex_size = tex_size,
            font = font,
            fontface = fontface,
            n_fill_breaks = 3,
            line_size = line_size,
            legend.position = 'top'
        ) ,
        heatmap_small(
            matrix = mat_ATRiHU_sdr,
            limits = c(-3, 3),
            more_or_less_label = 'B',
            line_order = line_order_small,
            name = 'URI (ATRiHU)',
            high_color = 'blue',
            low_color = 'red',
            mid_color = 'white',
            gene_list = region_small$name2,
            show_name = F,
            tex_size = tex_size,
            font = font,
            fontface = fontface,
            n_fill_breaks = 3,
            line_size = line_size,
            legend.position = 'top'
        ) ,
    nrow = 1,
    rel_widths = c(1.28, 1.06, 1, 1, 1,1)
)


ggsave(
    plot =Heat_maps_small,
    filename = paste0(output_dir,'/Fig4A.pdf'),
    width = 17,
    height = 11,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)
