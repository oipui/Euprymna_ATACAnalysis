
data20_1=read.table("103061_stage20_1_insertsizes.txt")
#change data to positive values: 
data20_1$V1=abs(as.numeric(data20_1$V1)) #column V1 as numeric, convert then to absolute and put into data20_1 column V1

data20_1[1,] #[row,column] to see first row 
data20_1[,1] # to see first column
#histogramm uses matrix, dataframe conversion to matrix with as.matrix
h<-hist(as.matrix(data20_1),breaks = 200, main='Fragmentsizes 103061_Stage20_1')
#show with h enter
plot(h$mids,h$counts,log='y', main ='Fragmentsizes 103061_Stage20_1')


library(magrittr)
library(dplyr)
library(ggplot2)

#97307 stage20 
data97307_st20=read.table("/Users/pui/Documents/Lab_Vienna/ATAC/Insertsizes/pos_97307_Stage20_insertsizes.txt")
pos_data97307_st20<-data97307_st20%>%filter(data97307_st20$V1>0)

fragLenPlot_97307_st20 <- table(pos_data97307_st20) %>% 
    data.frame %>% 
    rename(InsertSize = pos_data97307_st20, Count = Freq) %>% 
    mutate(InsertSize = as.numeric(as.vector(InsertSize)), Count = as.numeric(as.vector(Count)))  %>% 
  ggplot(aes(x = InsertSize, y = Count)) + geom_line()

#visualize 
fragLenPlot_97307_st20 + theme_bw() + ggtitle("97307 Stage20 Insertsizes")

#log scaled 
fragLenPlot_97307_st20 + scale_y_continuous(trans = "log2") + theme_bw() + ggtitle("97307 Stage20 Insertsizes")

#add nucleosome free (< 100bp), mono-nucleosome (180bp-247bp) and di-nucleosome (315-437) markings
fragLenPlot_97307_st20 + geom_vline(xintercept = c(180, 247), colour = "red") + 
  geom_vline(xintercept = c(315, 437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + 
  theme_bw() + ggtitle("97307 Stage20 Insertsizes")

#and the same for log scaled
fragLenPlot_97307_st20 + 
  scale_y_continuous(trans = "log2") + 
  geom_vline(xintercept = c(180, 247), colour = "red") + 
  geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
  geom_vline(xintercept = c(100), colour = "darkgreen") +
  theme_bw() + 
  theme(text = element_text(size=20)) + 
  ggtitle("A' 97307 Stage 20 ") + 
  ylab("Log2(Count)") + 
  xlab("Fragment size") 


#103061 Stage 20_1
data103061_st20_1=read.table("/Users/pui/Documents/Lab_Vienna/ATAC/Insertsizes/pos_103061_stage20_1_insertsizes.txt")
data103061_st20_1$V1=abs(as.numeric(data103061_st20_1$V1))

fragLenPlot_103061_st20_1 <- table(data103061_st20_1) %>% data.frame %>% rename(InsertSize = data103061_st20_1, Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), Count = as.numeric(as.vector(Count)))%>% ggplot(aes(x = InsertSize, y = Count)) + geom_line()

#visualize 
fragLenPlot_103061_st20_1 + theme_bw() + ggtitle("103061 Stage20_1 Insertsizes")

#log scaled 
fragLenPlot_103061_st20_1 + scale_y_continuous(trans = "log2") + theme_bw() +ggtitle("103061 Stage20_1 Insertsizes")

#add nucleosome free (< 100bp), mono-nucleosome (180bp-247bp) and di-nucleosome (315-437) markings
fragLenPlot_103061_st20_1 + geom_vline(xintercept = c(180, 247), colour = "red") + 
  geom_vline(xintercept = c(315, 437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + 
  theme_bw() + ggtitle("103061 Stage 20_1 Insertsizes")

#and the same for log scaled
fragLenPlot_103061_st20_1 + 
  scale_y_continuous(trans = "log2") + 
  geom_vline(xintercept = c(180, 247), colour = "red") + 
  geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
  geom_vline(xintercept = c(100), colour = "darkgreen") +
  theme_bw() + 
  theme(text = element_text(size=20)) + 
  ggtitle("A'' 103061 Stage 20 ") + 
  ylab("Log2(Count)") + 
  xlab("Fragment size") 



#97308 Stage25
data97308_st25=read.table("/Users/pui/Documents/Lab_Vienna/ATAC/Insertsizes/pos_97308_Stage25_insertsizes.txt")

fragLenPlot_97308_st25 <- table(data97308_st25) %>% data.frame %>% rename(InsertSize = data97308_st25, Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), Count = as.numeric(as.vector(Count)))%>% ggplot(aes(x = InsertSize, y = Count)) + geom_line()

#visualize 
fragLenPlot_97308_st25 + theme_bw() + ggtitle("97308 Stage25 Insertsizes")

#log scaled 
fragLenPlot_97308_st25 + scale_y_continuous(trans = "log2") + theme_bw() +ggtitle("97308 Stage25 Insertsizes")

#add nucleosome free (< 100bp), mono-nucleosome (180bp-247bp) and di-nucleosome (315-437) markings
fragLenPlot_97308_st25 + geom_vline(xintercept = c(180, 247), colour = "red") + 
  geom_vline(xintercept = c(315, 437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + 
  theme_bw() + ggtitle("97308 Stage 25 Insertsizes")

#and the same for log scaled
fragLenPlot_97308_st25 + 
  scale_y_continuous(trans = "log2") + 
  geom_vline(xintercept = c(180, 247), colour = "red") + 
  geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
  geom_vline(xintercept = c(100), colour = "darkgreen") +
  theme_bw() + 
  theme(text = element_text(size=20)) + 
  ggtitle("B' 97308 Stage 25 ") + 
  ylab("Log2(Count)") + 
  xlab("Fragment size") 

#103062 Stage 25
data103062_st25=read.table("/Users/pui/Documents/Lab_Vienna/ATAC/Insertsizes/pos_103062_Stage25_insertsizes.txt")

fragLenPlot_103062_st25 <- table(data103062_st25) %>% data.frame %>% rename(InsertSize = data103062_st25, Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), Count = as.numeric(as.vector(Count)))%>% ggplot(aes(x = InsertSize, y = Count)) + geom_line()

#visualize 
fragLenPlot_103062_st25 + theme_bw() + ggtitle("103062 Stage25 Insertsizes")

#log scaled 
fragLenPlot_103062_st25 + scale_y_continuous(trans = "log2") + theme_bw() +ggtitle("97308 Stage25 Insertsizes")

#add nucleosome free (< 100bp), mono-nucleosome (180bp-247bp) and di-nucleosome (315-437) markings
fragLenPlot_103062_st25 + geom_vline(xintercept = c(180, 247), colour = "red") + 
  geom_vline(xintercept = c(315, 437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + 
  theme_bw() + ggtitle("103062 Stage 25 Insertsizes")

#and the same for log scaled
fragLenPlot_103062_st25 +
  scale_y_continuous(trans = "log2") + 
  geom_vline(xintercept = c(180, 247), colour = "red") + 
  geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
  geom_vline(xintercept = c(100), colour = "darkgreen") +
  theme_bw() + 
  theme(text = element_text(size=20)) + 
  ggtitle("B'' 103062 Stage 25 ") + 
  ylab("Log2(Count)") + 
  xlab("Fragment size") 

#97309 Stage29
data97309_st29=read.table("/Users/pui/Documents/Lab_Vienna/ATAC/Insertsizes/pos_97309_Stage29_insertsizes.txt")

fragLenPlot_97309_st29 <- table(data97309_st29) %>% data.frame %>% rename(InsertSize = data97309_st29, Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), Count = as.numeric(as.vector(Count)))%>% ggplot(aes(x = InsertSize, y = Count)) + geom_line()

#visualize 
fragLenPlot_97309_st29 + theme_bw() + ggtitle("97309 Stage 29 Insertsizes")

#log scaled 
fragLenPlot_97309_st29 + scale_y_continuous(trans = "log2") + theme_bw() +ggtitle("97309 Stage 29 Insertsizes")


#add nucleosome free (< 100bp), mono-nucleosome (180bp-247bp) and di-nucleosome (315-437) markings
fragLenPlot_97309_st29 + geom_vline(xintercept = c(180, 247), colour = "red") + 
  geom_vline(xintercept = c(315, 437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + 
  theme_bw() + ggtitle("97309 Stage 29 Insertsizes")

#and the same for log scaled

fragLenPlot_97309_st29 + 
  scale_y_continuous(trans = "log2") + 
  geom_vline(xintercept = c(180, 247), colour = "red") + 
  geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
  geom_vline(xintercept = c(100), colour = "darkgreen") +
  theme_bw() + 
  theme(text = element_text(size=20)) + 
  ggtitle("C' 97309 Stage 29 ") + 
  ylab("Log2(Count)") + 
  xlab("Fragment size") 



#103063 Stage29
data103063_st29=read.table("/Users/pui/Documents/Lab_Vienna/ATAC/Insertsizes/pos_103063_Stage29_insertsizes.txt")

fragLenPlot_103063_st29 <- table(data103063_st29) %>% data.frame %>% rename(InsertSize = data103063_st29, Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), Count = as.numeric(as.vector(Count)))%>% ggplot(aes(x = InsertSize, y = Count)) + geom_line()

#visualize 
fragLenPlot_103063_st29 + theme_bw() + ggtitle("103063 Stage 29 Insertsizes")

#log scaled 
fragLenPlot_103063_st29 + scale_y_continuous(trans = "log2") + theme_bw() +ggtitle("103063 Stage 29 Insertsizes")


#add nucleosome free (< 100bp), mono-nucleosome (180bp-247bp) and di-nucleosome (315-437) markings
fragLenPlot_103063_st29 + geom_vline(xintercept = c(180, 247), colour = "red") + 
  geom_vline(xintercept = c(315, 437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + 
  theme_bw() + ggtitle("103063 Stage 29 Insertsizes")

#and the same for log scaled
fragLenPlot_103063_st29 + 
  scale_y_continuous(trans = "log2") + 
  geom_vline(xintercept = c(180, 247), colour = "red") + 
  geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
  geom_vline(xintercept = c(100), colour = "darkgreen") +
  theme_bw() + 
  theme(text = element_text(size=20)) + 
  ggtitle("C'' 103063 Stage 29 ") + 
  ylab("Log2(Count)") + 
  xlab("Fragment size") 
