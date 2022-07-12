#########################################################
#########################################################

# RUN Fig4A.R, Fig4B, Fig4C and FigS3A first


#########################################################
#########################################################
output_dir = '~/Desktop/Figure_paper_to_update/FigS3/FigS3B-D'
system(paste0('mkdir -p ',output_dir))

#add a color for SDR
Colors = c(Colors, 'SDR' = 'red')

# # boxplots for statistics
split_region300 = split(region_300, region_300$split_expression)

exclude_common_elements = function(x, y) {
    # convert a list into a string and removes the part that is in common between the two
    # it is used to remove the bins that are shared between the gene body and the gene body+flanking regions
    y = as.numeric(str_split(
        string =  str_replace(
            string = str_remove(
                pattern = paste(round(x, 4), collapse =  '_'),
                string = paste(round(y, 4), collapse =  '_')
            ),
            pattern = '__',
            replacement = '_'
        ),
        pattern = '_'
    )[[1]])
    return(y[!is.na(y)])
    
}

library(foreach)

# for the granges overlap we require a minimum overlap of 25Kb
min_ov = 25000
#### Fig S3 C

#### Aph
data = foreach(c = names(split_region300), .combine = 'rbind') %do% {
    # add flanking regions to the gene body 100kb both sides
    increased <-
        resize(split_region300[[c]],
               width = width(split_region300[[c]]) + 200000,
               fix = 'center')
    
    # find overlaps for the gene body
    hits1 = findOverlaps(subject = split_region300[[c]],
                         query = URI,
                         minoverlap = min_ov)
    # find overlaps for the gene body+flanking
    hits2 = findOverlaps(subject = increased,
                         query = URI,
                         minoverlap = min_ov)
    
    Gene = split_region300[[c]][subjectHits(hits1)] %>% as_tibble() %>%
        # assign a list of URI to a gene
        mutate(Gene = URI$Aph[queryHits(hits1)]) %>% dplyr::select(name2, Gene, split_expression) %>%
        drop_na() %>%
        group_by(name2, split_expression) %>%
        summarise(Gene = list(Gene))
    
    Flanking = increased[subjectHits(hits2)] %>% as_tibble() %>%
        # assign a list of URI to a gene+flanking
        mutate(Flanking = URI$Aph[queryHits(hits2)]) %>% dplyr::select(name2, Flanking, split_expression) %>%
        drop_na() %>%
        group_by(name2, split_expression) %>%
        summarise(Flanking = list(Flanking))
    
    #join the two df
    inner_join(Gene, Flanking) %>%
        rowwise() %>%
        # remove the sequence of bins in the gene body that is at the center of the gene+flanking
        mutate(
            Flanking = list(exclude_common_elements(
                x = unlist(Gene), y = unlist(Flanking)
            )),
            # calculate mean URI values for both of them
            Gene = mean(unlist(Gene), rm.na = T),
            Flanking = mean(unlist(Flanking), rm.na = T)
        ) %>%
        ungroup() %>%
        #reshape df
        gather(Feature, Value, Gene, Flanking) %>%
        spread(Feature, Value)
    
}

# calculate pvalues using wilcoxon
pval = data %>% group_by(split_expression) %>% summarise(pval = wilcox.test(
    x = Gene ,
    y = Flanking ,
    paired = T,
    alternative = 'less'
)$p.value) %>% ungroup()
# adjust with false discovery rate
pval$qval = p.adjust(pval$pval, method = 'fdr')
# filter signification pvalues
pval = pval[pval$qval < 0.05, ]

p1 = data %>%
    dplyr::rename(group = split_expression) %>%
    gather(Condition, Value, Gene, Flanking) %>%
    group_by(group) %>%
    ggplot() +
    geom_boxplot(aes(x = group, y = Value, fill = Condition), position = position_dodge(width = 0.6),width=0.4) +
    #add adjusted pval
    geom_text(
        data = pval,
        aes(
            x = split_expression,
            y = 3.5,
            label = format(qval, digits = 2, scientific = T)
        ),
        vjust = -0.25
    ) +
    annotate(
        'segment',
        x = (1:5) - 0.25,
        xend = (1:5) + 0.25,
        y = 3.5,
        yend = 3.5
    ) +
    labs(fill = '', x = 'Gene Group', y = 'mean URI in each gene body and relative flanking regions in Aph') +
    scale_fill_manual(values = c('Gene' = '#FCE621', 'Flanking' = '#CCCCCC')) +
    general_theme +
    stat_summary(
        fun.data = function(x) {
            return(c(y = -2, label = length(x)))
        } ,
        geom = "text",
        fun = median,aes(x = group, y = Value, group = Condition), position = position_dodge(width = 0.6)
    )

#save
ggsave(
    plot = p1,
    filename = paste0(output_dir, '/FigS3B.pdf'),
    width = 17,
    height = 14,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)


#### Fig S3 B only flanking (for stat)

# genes group
list_gene_grops = unique(data$split_expression)
# the foreach loop will go though all the unique combinations of gene groups
pval = foreach(a = 1:(length(list_gene_grops) - 1), .combine = 'rbind') %:%
    foreach(b = (a + 1):length(list_gene_grops),
            .combine = 'rbind') %do% {
                data %>%
                    # select the groups of interest
                    filter(split_expression == list_gene_grops[a] |
                               split_expression == list_gene_grops[b]) %>%
                    # select column of interest
                    dplyr::select(split_expression, Flanking) %>%
                    group_by(split_expression) %>%
                    #group into a list
                    summarise(Flanking = list(Flanking), x = 1) %>%
                    ungroup() %>%
                    # assign names
                    mutate(split_expression = ifelse(split_expression == list_gene_grops[a], 'Group1', 'Group2')) %>%
                    # spread values
                    spread(key = split_expression, value = Flanking) %>%
                    # run wilcoxon
                    mutate(
                        pval = wilcox.test(x = unlist(Group1) , y = unlist(Group2))$p.value,
                        Group1 = list_gene_grops[a],
                        Group2 = list_gene_grops[b]
                    ) %>%
                    dplyr::select(-x)
                
            }
# adjust with FDR
pval$qval = p.adjust(pval$pval, method = 'fdr')
# select significant qval
pval = pval[pval$qval < 0.05, ]
# adjust for plotting
pval = pval %>%
    #convert groups in factors
    mutate(
        Group1 = factor(Group1, levels = c('1st', '2nd', '3rd', '4th', 'SDR')),
        Group2 = factor(Group2, levels = c('1st', '2nd', '3rd', '4th', 'SDR')),
        #calculate center of two factors
        Center = (as.numeric(Group1) + as.numeric(Group2)) / 2
    ) %>%
    #arrange position y
    arrange(-(as.numeric(Group1) - as.numeric(Group2))) %>%
    mutate(y = seq(3, 5, length.out = n()))

p2 = data %>%
    #rename column
    dplyr::rename(group = split_expression) %>%
    # mutate groups into factor
    mutate(group = factor(group, levels = c('1st', '2nd', '3rd', '4th', 'SDR'))) %>%
    ggplot() +
    #plot boxplot of flanking regions
    geom_boxplot(aes(x = group, y = Flanking, fill = group),width=0.4) +
    #add adjusted pval
    geom_text(data = pval,
              aes(
                  x = Center,
                  y = y,
                  label = format(qval, digits = 2, scientific = T)
              ),
              vjust = -0.25) +
    geom_segment(data = pval, aes(
        x = Group1 ,
        xend = Group2,
        y = y,
        yend = y
    )) +
    labs(fill = '', x = 'Gene Group', y = 'Average Aph URI in each flanking associated to a gene') +
    scale_fill_manual(values = Colors) + general_theme+
    stat_summary(
        fun.data = function(x) {
            return(c(y = -2, label = length(x)))
        } ,
        geom = "text",
        fun = median,aes(x = group, y = Flanking, group = group), position = position_dodge(width = 0.6)
    )

#save
ggsave(
    plot = p2,
    filename = paste0(output_dir, '/FigS3B_only flanking.pdf'),
    width = 17,
    height = 14,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)


#### Fig S3 C left
hit = findOverlaps(query = URI,
                   subject = NT,
                   minoverlap = min_ov)

URI$RT=NA
### check that there are no duplicates
if (length(subjectHits(hit)) == length(unique(subjectHits(hit))) &
    length(queryHits(hit)) == length(unique(queryHits(hit)))) {
    URI$RT[queryHits(hit)] = NT$s50[subjectHits(hit)]
    
} else{
    print('Error!!!!!')
}


filter_by_RT = function(X, RT, rt_limit = 0.5) {
    # it takes a list of bins wit URI values and matching RT values
    # it filters out bins based on an RT value limit
    RT = unlist(RT)
    X = unlist(X)
    return(list(X[which(RT >= rt_limit)]))
}



# calculates each group separately
data = foreach(c = names(split_region300), .combine = 'rbind') %do% {
    # find overlaps for the gene body
    hits = findOverlaps(subject = split_region300[[c]],
                        query = URI,
                        minoverlap = min_ov)
    
    Gene = split_region300[[c]][subjectHits(hits)] %>% as_tibble() %>%
        #same as before but we keep Aph and ARO value
        mutate(Aph = URI$Aph[queryHits(hits)],
               ARO = URI$AphRO[queryHits(hits)],
               RT = URI$RT[queryHits(hits)]) %>% dplyr::select(name2, Aph, ARO, RT, split_expression) %>%
        drop_na() %>%
        group_by(name2, split_expression) %>%
        summarise(Aph_G = list(Aph),
                  ARO_G = list(ARO),
                  RT_G = list(RT)) %>%
        rowwise() %>%
        # remove the gene body bins from flanking+gene body region
        mutate(
            # filter RT >= 0 (all the bins)
            Aph_G = filter_by_RT(X = Aph_G, RT_G, rt_limit = 0),
            ARO_G = filter_by_RT(X = ARO_G, RT_G, rt_limit = 0),
            #calculate average UR
            Aph_G = mean(unlist(Aph_G), rm.na = T),
            ARO_G = mean(unlist(ARO_G), rm.na = T)
        ) %>%
        ungroup() %>%
        drop_na() %>%
        #reshape data
        dplyr::select(name2, group = split_expression, Aph_G, ARO_G) %>%
        gather(type, URI, Aph_G, ARO_G) %>%
        separate(type,
                 into = c('Treatment', 'Region'),
                 sep = '_') %>%
        mutate(Region = 'Gene body')
    
}

pval = data %>%
    #reshape data
    group_by(group) %>% spread(Treatment, URI) %>%
    #run wilcoxon test per group
    summarise(pval = wilcox.test(
        x = Aph ,
        y = ARO ,
        paired = T,
        alternative = 'less'
    )$p.value) %>% ungroup()
#adjust pvalue
pval$qval = p.adjust(pval$pval, method = 'fdr')

#filter non sig
pval = pval[pval$qval < 0.05, ]

p3 = data %>%
    ggplot() +
    # plot boxplot
    geom_boxplot(aes(x = group, y = URI, fill = Treatment), position = position_dodge(width = 0.6),width=0.4)+
    #add adjusted pvalue
    geom_text(data = pval,
              aes(
                  x = group,
                  y = 3,
                  label = format(qval, digits = 2, scientific = T)
              ),
              vjust = -0.25) +
    # put a segment over all the groups
    annotate(
        'segment',
        x = (1:5) - 0.25,
        xend = (1:5) + 0.25,
        y = 3,
        yend = 3
    ) +
    #labels
    labs(fill = '', x = 'Gene Group', y = 'mean URI in each gene body') +
    scale_fill_manual(values =  c('Aph' = '#95C11F',
                                  'ARO' = '#F9B233')) + general_theme + coord_cartesian(ylim = c(-3, 3.3))+
    stat_summary(
        fun.data = function(x) {
            return(c(y = -3, label = length(x)))
        } ,
        geom = "text",
        fun = median,aes(x = group, y = URI, group = Treatment), position = position_dodge(width = 0.6)
    )
#save
ggsave(
    plot = p3,
    filename = paste0(output_dir, '/FigS3C_left.pdf'),
    width = 17,
    height = 14,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)
#### Fig S3 C right

# calculates each group separately
data_later = foreach(c = names(split_region300), .combine = 'rbind') %do%
{
    # find overlaps for the gene body
    hits = findOverlaps(subject = split_region300[[c]],
                        query = URI,
                        minoverlap = min_ov)
    
    Gene = split_region300[[c]][subjectHits(hits)] %>% as_tibble() %>%
        #same as before but we keep Aph and ARO value
        mutate(Aph = URI$Aph[queryHits(hits)],
               ARO = URI$AphRO[queryHits(hits)],
               RT = URI$RT[queryHits(hits)]) %>% dplyr::select(name2, Aph, ARO, RT, split_expression) %>%
        drop_na() %>%
        group_by(name2, split_expression) %>%
        summarise(Aph_G = list(Aph),
                  ARO_G = list(ARO),
                  RT_G = list(RT)) %>%
        rowwise() %>%
        # remove the gene body bins from flanking+gene body region
        mutate(
            # filter RT >= 0 (all the bins)
            Aph_G = filter_by_RT(X = Aph_G, RT_G, rt_limit = quantile(unlist(RT_G), 0.5)[1]),
            ARO_G = filter_by_RT(X = ARO_G, RT_G, rt_limit = quantile(unlist(RT_G), 0.5)[1]),
            #calculate average UR
            Aph_G = mean(unlist(Aph_G), rm.na = T),
            ARO_G = mean(unlist(ARO_G), rm.na = T)
        ) %>%
        ungroup() %>%
        drop_na() %>%
        #reshape data
        dplyr::select(name2, group = split_expression, Aph_G, ARO_G) %>%
        gather(type, URI, Aph_G, ARO_G) %>%
        separate(type,
                 into = c('Treatment', 'Region'),
                 sep = '_') %>%
        mutate(Region = 'Gene body')
    
}

pval = data_later %>%
    #reshape data
    group_by(group) %>% spread(Treatment, URI) %>%
    #run wilcoxon test per group
    summarise(pval = wilcox.test(
        x = Aph ,
        y = ARO ,
        paired = T,
        alternative = 'less'
    )$p.value) %>% ungroup()
#adjust pvalue
pval$qval = p.adjust(pval$pval, method = 'fdr')

#save pval_later
pval_later_gene=pval

#filter non sig
pval = pval[pval$qval < 0.05, ]


p4 = data_later %>%
    ggplot() +
    # plot boxplot
    geom_boxplot(aes(x = group, y = URI, fill = Treatment), position = position_dodge(width = 0.6),width=0.4) +
    #add adjusted pvalue
    geom_text(data = pval,
              aes(
                  x = group,
                  y = 4,
                  label = format(qval, digits = 2, scientific = T)
              ),
              vjust = -0.25) +
    # put a segment over all the groups
    annotate(
        'segment',
        x = (1:5) - 0.25,
        xend = (1:5) + 0.25,
        y = 4,
        yend = 4
    ) +
    #labels
    labs(fill = '', x = 'Gene Group', y = 'mean URI in each gene body per beans with \nS50 >= 50 percentile of the S50') +
    scale_fill_manual(values =  c('Aph' = '#95C11F',
                                  'ARO' = '#F9B233')) + general_theme+
    stat_summary(
        fun.data = function(x) {
            return(c(y = -3, label = length(x)))
        } ,
        geom = "text",
        fun = median,aes(x = group, y = URI, group = Treatment), position = position_dodge(width = 0.6)
    )
ggsave(
    plot = p4,
    filename = paste0(output_dir, '/FigS3C_right.pdf'),
    width = 17,
    height = 14,
    limitsize = F,
    units = 'cm',
    device = cairo_pdf
)

#figure S3 D

URI=read_tsv('URi_simulation_periodic_1kbXY.tsv')%>%
    filter(Condition %in% c('Aph','AphRO'))
    
S50=read_tsv('S50_simulation_periodic_1kbXY.tsv')%>%
    filter(Condition=='NT')%>%
    select(-Condition)

URI=URI%>%inner_join(S50)%>%
    mutate(s50 = case_when(
                            0.1 <= s50 & s50 < 0.33 ~ 'Early',
                               0.33 <= s50 & s50 < 0.66 ~ 'Mid',
                               0.66 <= s50 & s50 <= 0.9 ~ 'Late'),
    s50=factor(s50,levels = c('Early','Mid', 'Late')))



pval = URI %>%
    spread(Condition,URI)%>%
    drop_na()%>%
    #reshape data
    group_by(s50) %>%
    #run wilcoxon test per group
    summarise(pval = wilcox.test(
        x = Aph ,
        y = AphRO ,
        paired = T
    )$p.value,
    pval_less = wilcox.test(
        x = Aph ,
        y = AphRO ,
        paired = T,
        alternative = 'less'
    )$p.value) %>% ungroup()

pval$qval=p.adjust(pval$pval,method = 'fdr')
pval$qval_less=p.adjust(pval$pval_less,method = 'fdr')

pval=pval%>%
    mutate(qval=ifelse(qval>=0.05,'NS',format(qval, digits = 2, scientific = T)))

p5 = URI %>%
    drop_na()%>%
    dplyr::select('group'=s50,URI,Condition)%>%
    mutate(Condition=ifelse(Condition=='AphRO','ARO',Condition))%>%
    ggplot() +
    # plot boxplot
    geom_boxplot(aes(x = group, y = URI, fill = Condition), position = position_dodge(width = 0.6),width=0.4) +
    #add adjusted pvalue
    geom_text(data = pval,
              aes(
                  x = s50,
                  y = 10,
                  label = qval
              ),
              vjust = -0.25) +
   
    # put a segment over all the groups
    annotate(
        'segment',
        x = (1:3) - 0.25,
        xend = (1:3) + 0.25,
        y = 10,
        yend = 10
    ) +
    #labels
    labs(fill = '', x = 'Replication Timing', y = 'URI') +
    scale_fill_manual(values =  c('Aph' = '#95C11F',
                                  'ARO' = '#F9B233')) + general_theme+
    stat_summary(
        fun.data = function(x) {
            return(c(y = -4, label = length(x)))
        } ,
        geom = "text",
        fun = median,aes(x = group, y = URI, group = Condition), position = position_dodge(width = 0.6)
    )


ggsave(p5,filename = '~/Desktop/Figure_paper_to_update/FigS3/FigS3B-D/FigS3D.pdf',device = cairo_pdf,width = 10,height = 5)

