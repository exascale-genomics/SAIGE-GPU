# loading libraries. --------------------------------------
# https://github.com/milesilab/peakdetection/tree/v1.1
if( !requireNamespace("data.table", quietly = T) )
  install.packages("data.table")
if( !requireNamespace("dplyr", quietly = T) )
  install.packages("dplyr")
if( !requireNamespace("ggplot2", quietly = T) )
  install.packages("ggplot2")

library(data.table)
library(dplyr)
library(ggplot2)
library(R.utils)

# declaring variables. ------------------------------------
vars = 6
slw.size       = 300
smooth.intens  = 0.90
extens         = "pcaAdapt" ### nom pr output: extens.var.
type           = "Pval" #"BF" #"Z" # Pval
proba          = 0.95 # to be specified
iter           = 1000
p_val_threshold = 1e-8

## helping functions. -------------------------------------
computing.profile = function( y, dens.cut )
{
  profile = y > dens.cut
  profile[is.na(profile)] = 0
  counter = profile[slw.size:length(profile)] -
    c(0, profile[1:(length(profile)-slw.size)])
  initial.sum = sum( profile[1:(slw.size-1)] )

  return( c( rep(NA, slw.size-1), initial.sum + cumsum(counter) ) )
}

randomization = function( .data, iter, dens.cut )
{
  replicate(
    iter,
    .data$y %>% sample %>%
      computing.profile( dens.cut ) %>% max( na.rm = T )
  )
}

is.defined <- function(sym) {
  sym <- deparse(substitute(sym))
  env <- parent.frame()
  exists(sym, env)
}

# loading datasets. ---------------------------------------
gwas_output_path = "/ccs/home/arodriguez/med112/task0101113/output/HARE_ANC_Run"
phecodes = list.dirs(gwas_output_path, recursive=FALSE, full.names=FALSE)
group = "ASN"
output.path = "/ccs/home/arodriguez/med112/task0101113/output/HARE_ANC_Run/SUBMIT/peaks"
#gwas_path = "/ccs/home/arodriguez/med112/task0101113/output/HARE_ANC_Run/Phe_454_1/ASN/step2/"

for (phecode in phecodes) {
  gwas_path = paste(gwas_output_path, phecode, group, "step2", sep="/")
  if (dir.exists(gwas_path)) {  
    final_stats <- matrix()
    gwas_files = list.files(path = gwas_path, pattern = ".assoc.gpu.txt$")
    dir.create( paste(output.path, "stat/", sep="/"), showWarnings = F )

    for (gwas_f in gwas_files) {
       print(gwas_f)
       gwas = fread(paste(gwas_path, gwas_f, sep="/")) %>% as_tibble()

       ## Peak detection

       dataf2 = 
         gwas %>% 
           mutate( Scaff = chrom, Pos = 1:n(), ID = rsid ) %>% 
           select( Scaff, pos, rsid, chrom, all_of(vars) )

       liss    = matrix( NA, nrow(dataf2), length(vars) )
       signals = list() 
       length(signals) = length(vars)

       ## peak detection. ----------------------------------------
       for( var in length(vars) )
       {
           # loading data. -------------------------------
           varname = names(dataf2)[4+var]
           dataf = 
             dataf2 %>% 
            rename( y = !!varname ) %>% 
            mutate(
                y = case_when( type == "BF" | type == "Z" ~ log10(y),
                       type == "Pval" ~ -log10(y) ))
  
           dens.cut = quantile( dataf$y, p = proba )
           cut.quantile = dataf %>% randomization(iter, dens.cut) %>% max
  
           # computing the density. ----------------------
           dataf = 
               dataf %>% 
               group_by(chrom) %>% 
               mutate(
                  avg    = computing.profile( y, dens.cut ),
                  smooth = ifelse( avg < quantile( avg, p = smooth.intens, na.rm = T ),
                         0, avg )
                  ) %>% ungroup
  
           # outputting. ---------------------------------
           # smoothed peaks. 
           liss[, var] = dataf %>% pull(smooth)
  
           # positions of the peaks. 
           .signals  = dataf %>% filter( avg > cut.quantile ) %>% pull(pos)
           .dsignals = diff(.signals)
  
           impl.signals = NULL
           i = 1
           while( i < length(.signals) )
           {
             j = i
             # strictly lower to avoid inter LG peaks. 
             while( j < length(.dsignals) && .dsignals[j] < slw.size ) j = j + 1
    
             impl.signals = 
                 c( impl.signals, 
                 dataf %>% 
                     slice( .signals[i:(j-1)] ) %>% # ~ filter( Pos %in% .signals )
                     filter( avg == max(avg) ) %>% 
                     pull(pos) %>% median() %>% round())
             i = j + 1
           }
  
           signals[[var]] = impl.signals
  
           ## displaying. ------------------------------------------
           print( paste( round(var/length(vars)*100, 0), " %", sep = "" ) )
       }

       peak.annotation = list()
       cand = list()
       cand2 = list()
       snp.sign = list()
       scaff.sign = list()

       for( var in seq_along(vars) )
       {
           # computing the peak limits. ------------------
           peaklimit = array(NA, c( length(signals[[var]]), 2 ) )
           for( i in seq_along(signals[[var]]) )
           {
               if (is.infinite(which( liss[,var] == 0 &
                     seq_along(liss[,var]) < signals[[var]][i] ) %>% max)) {
	          peaklimit[i, 1] = signals[[var]][i]
               }else {
      	           peaklimit[i, 1] = which( liss[,var] == 0 & 
                     seq_along(liss[,var]) < signals[[var]][i] ) %>% max
               }
               if (is.infinite(which( liss[,var] == 0 &
                     seq_along(liss[,var]) > signals[[var]][i] ) %>% min)) {
	           peaklimit[i, 2] = signals[[var]][i]
               }else {
                   peaklimit[i, 2] = which( liss[,var] == 0 & 
                       seq_along(liss[,var]) > signals[[var]][i] ) %>% min
               }
           }
  
           # annotating the peaks. -----------------------
           .peak.annotation = matrix( NA, nrow(peaklimit), 8 )
           intervals  = c()
           intervals2 = c()
           .snp.sign = c()
           .scaff.sign = c()
           for( i in 1:nrow(peaklimit) )
           {
               interval = peaklimit[i,1]:peaklimit[i,2]
               tmp.data = dataf2[interval, 4+var] %>% pull
    
               .peak.annotation[i,1] = 
               ifelse( type == "BF" || type == "Z", max( tmp.data ), min( tmp.data ) )
    
               .peak.annotation[i,2] = 
               ifelse( type == "BF" || type == "Z", interval[ which.max(tmp.data) ], interval[ which.min(tmp.data) ])
    
               .peak.annotation[i,3] = min(interval)
               .peak.annotation[i,4] = max(interval)
               .peak.annotation[i,5] = max(interval) - min(interval)
    
               if(type == "BF" | type == "Z") {
                   qpeak = which( log10(tmp.data) > quantile( dataf$y, p = proba) )
               } else 
                   qpeak = which( -log10(tmp.data) > quantile( dataf$y, p = proba ) )
    
                   .peak.annotation[i,6] = length( qpeak )
    
                   intervals = c( intervals, peaklimit %>% t %>% c)
                   intervals2 = c(intervals2, peaklimit %>% t %>% c)
               }
  
               colnames( .peak.annotation ) = c("Best", "Pos", "Start", "End", "Wide", "SNPs", "Genes", "Sign" )
               rownames( .peak.annotation ) = dataf2 %>% slice( .peak.annotation[,2] ) %>% pull(ID) 
  
               peak.annotation[[var]] = .peak.annotation
               cand[[var]] = unique(intervals)
               cand2[[var]] = unique(intervals2)
               snp.sign[[var]] = .snp.sign
               scaff.sign[[var]] = .scaff.sign
           }

           E = sort(-log10(runif(nrow(dataf2),0,1)))
  
           run<-0
           for( var in seq_along(vars) )
           {
               if (is.defined(stats)) { 
                   stats = rbind(stats, peak.annotation[[var]])
               } else {
                   stats = peak.annotation[[var]]
               }

               stats = subset(stats,stats[, "Best"] < p_val_threshold) 
               final_stats <- stats
           }
       }
    }
    if( length(final_stats) == 0){
        # do nothing.
    } else {
        stat_file = paste(output.path, paste(phecode, group, "peak", extens, colnames(dataf2)[[4+var]], proba, "txt", sep="."), sep="/")
        write.table(final_stats, stat_file, col.names = T, row.names = T, quote = F, sep = "\t" )
    }
  }
}
```

