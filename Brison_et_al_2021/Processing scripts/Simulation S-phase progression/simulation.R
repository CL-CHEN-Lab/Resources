library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

results=paste(c(
    'Forks',
    'Initiations',
    'termination',
    'Activated_origins',
    'Passivly_replicated_origins',
    'replicated_dna',
    'unreplicated_dna',
    'percent' ,
    'cycles',
    'time',
    'simulation',
    'Free_factor',
    'Kod',
    'genome_size_kb',
    'delay_LF',
    'available_origins',
    'total_limiting_factors',
    'total_available_limiting_factors',
    'produced_limiting_factors'),collapse = '\t')

system(paste0('echo ',results))

bin_size=1000
p0=0.012
ref=read_tsv('./hg38.chr.size',col_names = c('chr','size'))%>%
     rowwise()%>%
    mutate(chr=ifelse(!chr %in% c('chrX','chrY'),list(paste0(chr,'_',c(1,2))),list(paste0(chr,'_',c(1)))))%>%
    unnest(chr)%>%
    group_by(chr)%>%mutate(size=round(size/bin_size),
             pos=list(1:size))%>%
    unnest(pos)%>%
    ungroup()%>%
    mutate(n=1:n())%>%
    group_by(chr)%>%
    mutate(origins = ifelse(n %% round((1/p0)*1000/bin_size) ==0,'origins',0),
           borders=ifelse(n==min(n)| n==max(n),'Border',NA))
    
    


                      # reset system
                      border=ref%>%filter(borders=='Border')
                      border=border$n
                      origins=ref%>%filter(origins=='origins')%>%
                          ungroup()
                      
                      origins=origins
                      origins=origins$n
                      
                      origins= origins[!origins %in% border]
                      
                      n_genome=length(ref$n)
                      genome=rep(1,n_genome)
                      n_genome_kb=n_genome*bin_size/1000
                      
                      next_fork_f=0
                      next_fork_r=0
                      memory_forks=0
                      Passive_origins=0
                      Activated_origins=0
                      percent=0
                      Free_factors=0
                      memory_firing=0

                      origins_density=length(origins)/n_genome_kb
                      fork_speed=1.46
                      bin=bin_size/1000
                      referece_length_genome=3000
                      alpha=n_genome_kb/referece_length_genome
                      total_free_factors=round(12*alpha)
                      Kod=(1.2*10^-2)/alpha
                      time_cycle=(bin/fork_speed)
                      delay=15/(time_cycle)
                      total_cycle=0
                      while (percent!=100 & total_cycle < 10000) {
                          
                        total_cycle=total_cycle+1
                                 #The eukaryotic bell-shaped temporal rate of DNA replication origin firing emanates from a balance between origin activation and passivation Jean-Michel Arbona, Alain Arneodo4, Arach Goldar, Olivier Hyrien, Benjamin Audit1 
                                 #limiting factors increas gradually as described in 3D replicon distributions arise from stochastic initiation and domino-like DNA replication progression D. Löb, N. Lengert, V. O. Chagin, M. Reinhart, C. S. Casas-Delucchi, M. C. Cardoso & B. Drossel 
                                Free_factors=Free_factors+round(total_free_factors*
                                                                    ((1-exp(-total_cycle/delay))-(1-exp(-(total_cycle-1)/delay))))
                                
                                 prob_firing=1-(1-Kod)^length(origins)
                                 n_of_firing_origins=runif(trunc(Free_factors/2))
                                 n_of_firing_origins=trunc(sum(prob_firing>n_of_firing_origins))
                                 n_of_firing_origins=ifelse(n_of_firing_origins > trunc(Free_factors/2),trunc(Free_factors/2),n_of_firing_origins)
                                 n_of_firing_origins=ifelse(n_of_firing_origins > length(origins),length(origins),n_of_firing_origins)
                                 
                                 Free_factors=Free_factors-2*n_of_firing_origins
                                 
                          #progression of the fork
                          # previous cycle orgins are replicated (value 3, forks)
                          #previous fork position
                          # previous forks pass to stage 5
                          genome[genome == 4] = 5
                          genome[genome == 3] = 4
                          genome[genome == 2] = 3
                          
                          # id replication foks
                          
                          forks = which(genome == 3)
                          next_fork_f = forks[!forks %in% border] + 1
                          next_fork_r = forks[!forks %in% border] - 1
                          
                          next_fork_r=next_fork_r[genome[next_fork_r]!=4]
                          next_fork_f=next_fork_f[genome[next_fork_f]!=4]
                          
                          next_fork_f = next_fork_f[genome[next_fork_f] == 1]
                          next_fork_r = next_fork_r[genome[next_fork_r] == 1]
                          
                          genome[next_fork_f] = 2
                          genome[next_fork_r] = 2
                          genome[next_fork_f[next_fork_f %in% next_fork_r]]=3
                          # calculate passive origins 
                          Passive_origins=Passive_origins+length(origins)-length(origins[genome[origins]==1])
                          #remove passive origins from the pool
                          origins=origins[genome[origins]==1]
                          
                          if(length(origins) >  n_of_firing_origins  ) {
                              
                              firing = sample(origins,size = n_of_firing_origins)
                              origins = origins[!origins %in% firing]
                              
                          }   else if(length(origins) ==  n_of_firing_origins & length(origins) != 0){
                              firing = origins
                              origins = origins[!origins %in% firing]
                          }else{
                              firing=NULL
                          }
                            
                          genome[firing]  = 2
                          firing_2=firing-1
                          firing_2=firing[firing_2 > 1]
                          firing_2=firing_2[genome[firing_2]==1]
                          genome[firing_2]  = 2
                       
                          origins=origins[genome[origins]==1]
                          Activated_origins=Activated_origins+length(firing)
                          percent = 100*sum(genome != 1) /n_genome
                      
                          Forks =sum(!next_fork_f %in% next_fork_r)+sum(!next_fork_r %in% next_fork_f)+length(firing)+length(firing_2)
                          termination_events=ifelse((memory_forks+length(firing)+length(firing_2)-Forks)<0,0,memory_forks+length(firing)+length(firing_2)-Forks)
                          Free_factors=Free_factors+termination_events
                          
                          results=paste(
                              Forks,
                              length(firing),
                              termination_events,
                              Activated_origins,
                              Passive_origins,
                              sum(genome != 1)*bin_size/1000,
                              sum(genome == 1)*bin_size/1000,
                              percent,
                              total_cycle,
                              time_cycle*total_cycle,
                              args,
                              Free_factors,
                              Kod,
                              n_genome_kb,
                              delay,
                              length(origins),
                              total_free_factors,
                              round(total_free_factors*(1-exp(-total_cycle/delay))),
                              round(total_free_factors*((1-exp(-total_cycle/delay))-(1-exp(-(total_cycle-1)/delay)))),sep = '\t')
                                                          
                          memory_forks=sum(!next_fork_f %in% next_fork_r)+sum(!next_fork_r %in% next_fork_f)+2*length(firing)
                          system(paste0('echo ',results))
             
                    }
                    
                    


