library(plotly)
library(ggplot2)
problem = 'ICA'
setwd(paste('//loki/export/mialab/users/rsilva/repos/code/MISA/results/sim/',
            problem, sep = ''
            )
      )

makemynames <- function(mole){
  moles = NULL
  for(kk in 1:length(mole)){
    moles = c(moles, sprintf('%.2g',mole[kk]))
  }
  return(moles)
}

tab1 = NULL
tab2 = NULL
for(Run in 1:10){
  # for(r0 in 1:10){
    fname = paste(problem,
                  'sim_w0_results_rr', sprintf('%03d',Run),
#                   '_ir0', sprintf('%03d',r0),
#                   '_er0', sprintf('%03d',r0),
                  '.csv', sep = '')
    tmp = read.csv(file = fname, head = TRUE, sep = ",")
#     if(Run == 6){
#       print(Run)
#       print(dim(tmp))
#       print(tmp)
#       # print(which(is.na(tmp$SNRdB)))
#     }
    tab1 = rbind(tab1, tmp[1:80,])#
    tab2 = rbind(tab2, tmp[81:180,])#
  # }
}
tab2 <- tab2[!is.na(tab2$SNRdB),] # Removes repeats that did not run
SNR = unique(tab2$SNRdB)
SNRs = makemynames(SNR)
Acond = unique(tab1$Acond)
Aconds = makemynames(Acond)
Alg = levels(tab2$Alg)
Algs = c(Alg[1],'RE+MISA')

# For ggplot2 labeller below:
Acond_names <- as.vector(apply(as.matrix(Aconds),1,function(x) paste('cond(A) = ',x,sep='')))
names(Acond_names) <- Aconds
SNR_names <- as.vector(apply(as.matrix(SNRs),1,function(x) paste('SNR = ',x,'dB',sep='')))
names(SNR_names) <- SNRs
Run_names <- as.vector(apply(as.matrix(1:10),1,function(x) paste('Run ',x,sep='')))
names(Run_names) <- as.character(1:10)

tab1$MISI <- tab1$MISI + .Machine$double.eps

tab1$SNRdB <- factor(tab1$SNRdB, levels = SNR, labels = SNRs)
tab1$Acond <- factor(tab1$Acond, levels = Acond, labels = Aconds)
tab1$Run <- as.factor(tab1$Run)
tab1$Initial <- factor(tab1$Initial, levels = c('Wi','W','W0'))
tab1$Alg <- factor(tab1$Alg, levels = Alg, labels = Algs)

tab2$MISI <- tab2$MISI + .Machine$double.eps

tab2$SNRdB <- factor(tab2$SNRdB, levels = SNR, labels = SNRs)
tab2$Acond <- factor(tab2$Acond, levels = Acond, labels = Aconds)
tab2$Run <- as.factor(tab2$Run)
tab2$Initial <- factor(tab2$Initial, levels = c('Wi','W','W0'))
tab2$Alg <- factor(tab2$Alg, levels = Alg, labels = Algs)

# MISI_lims = c(min(min(tab1$MISI), min(tab2$MISI)), max(max(tab1$MISI), max(tab2$MISI)))
# MISI_lims[1] <- 10^floor(log10(MISI_lims[1]))
# # MISI_lims
# MMSE_lims = c(min(min(tab1$MMSE), min(tab2$MMSE)), max(max(tab1$MMSE), max(tab2$MMSE)))
# MMSE_lims[1] <- 10^floor(log10(MMSE_lims[1]))
# # MMSE_lims
# equal_lims = c(min(MISI_lims,MMSE_lims),max(MISI_lims,MMSE_lims))
# equal_lims
equal_lims = c(0.002, 1.6)

# tab2 <- tab[-which(tab$SNRdB == levels(tab$SNRdB)[1]),]
p <- ggplot(tab1,aes(x=MISI, y=MMSE, colour = Alg))
#p <- p + geom_vline(linetype = 2, colour = 'black', xintercept = 79.6733)#log10(.Machine$double.eps))
#p <- p + geom_path(data = tab[tab$State == 'Mid',])#, aes(linetype = RElam, fill = Replicate))
p <- p + geom_point(aes(shape = Alg, size = Initial))#shape = Initial, alpha = 1))
#p <- p + scale_shape_discrete(solid=F)
p <- p + scale_shape_manual(values=c(2,1))
#p <- p + scale_size_manual(values=c(2,2,3))
p <- p + scale_color_manual(values=c('#e6ab02','#e7298a'))
p <- p + facet_grid(Run~Acond, labeller = labeller(Acond = Acond_names, Run = Run_names))#, margins=F, space = "free")
# print(p)
p <- p + scale_x_log10()
p <- p + scale_y_log10()
p <- p + expand_limits(x = equal_lims, y = equal_lims)
p <- p + ggtitle('SNR = 3dB')
p <- p + theme(aspect.ratio=1)
print(p)
# Save pdf with Device: 6.7 (width) x 11.50 (height)

p <- ggplot(tab2,aes(x=MISI, y=MMSE, colour = Alg))
#p <- p + geom_vline(linetype = 2, colour = 'black', xintercept = 79.6733)#log10(.Machine$double.eps))
#p <- p + geom_path(data = tab[tab$State == 'Mid',])#, aes(linetype = RElam, fill = Replicate))
p <- p + geom_point(aes(shape = Alg, size = Initial))#shape = Initial, alpha = 1))
#p <- p + scale_shape_discrete(solid=F)
p <- p + scale_shape_manual(values=c(2,1))
#p <- p + scale_size_manual(values=c(2,2,3))
p <- p + scale_color_manual(values=c('#e6ab02','#e7298a'))
p <- p + facet_grid(Run~SNRdB, labeller = labeller(SNRdB = SNR_names, Run = Run_names))#, margins=F, space = "free")
# print(p)
p <- p + scale_x_log10()
p <- p + scale_y_log10()
p <- p + expand_limits(x = equal_lims, y = equal_lims)
p <- p + ggtitle('cond(A) = 7')
p <- p + theme(aspect.ratio=1) 
print(p)
# Save pdf with Device: 7.75 (width) x 11.50 (height)

gg <- ggplotly(p)
print(gg)

################################### OLD

T = read.csv(file='pareto_case16_r020_rA000_rr001_branch001_ss002_optmodestrict_bb004.csv',head=TRUE,sep=",")
T$Range <- as.factor(T$Range)
T$opt <- factor(T$opt, levels = c('MISA','none','REP','RET','RET_U0'))
T$state <- as.factor(T$state)
T$initial <- factor(T$initial, levels = c('w0','w0_REP','w0_RET','w0_RET_U0','GT'))
T$MISAabs <- abs(T$MISA)
# T$MISA.scl <- log10(T$MISA - min(T$MISA))
T$MISA.scl <- T$MISA - 79.6733
T <- T[T$MISA.scl >= 0,]
T$MISA.scl <- log10(T$MISA.scl - min(T$MISA.scl))
T$MISA.scl[T$MISA.scl == -Inf] <- log10(.Machine$double.eps)
p <- ggplot(T,aes(x=MISA.scl, y=RE))
p <- p + geom_vline(linetype = 2, color = 'black',xintercept = 79.6733)#log10(.Machine$double.eps))
#p <- p + geom_line(data = T[T$Range == -4,], colour = 'black')
p <- p + geom_point(aes(color = opt, shape = state, size = state), alpha = .5)
p <- p + scale_shape_manual(values=c(15,16,17))
p <- p + scale_size_manual(values=c(4,1,4))
p <- p + scale_y_log10()
p <- p + facet_grid(initial~Range, margins=TRUE)
print(p)
gg <- ggplotly(p)
print(gg)

p <- p + scale_y_log10(breaks = 10^(seq(-9,0,by = 1)))
gg <- ggplotly(p)
print(gg)


T = read.csv(file='pareto_case16_r020.csv',head=TRUE,sep=",")
T2 <- T[T$Range <= .1 & T$Range >=0,]
T2$MISA <- (T2$MISA - T2$MISA[T2$Range == 0])
T2$Range <- as.factor(T2$Range)
p <- ggplot(T2)
p <- p + geom_point(aes(y = T2$MISA, x = as.numeric(row.names(T2)), colour = as.factor(T2$Range)))
#p <- p + scale_y_log10()
gg <- ggplotly(p)
print(gg)
