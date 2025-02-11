library(ggplot2)

load("Res_all_100var_n2.RData")
CV_arr = array(unlist(CV), dim = c(7, 8, 100)) #7 vs 13
CV_avg = apply(CV_arr, c(1,2), mean)
noquote(paste0(round(t(CV_avg*100),2)[1,c(1,6,3)], "%&"))
noquote(paste0(round(t(CV_avg*100),2)[2,c(1,6,3)], "%&"))
noquote(paste0(round(t(CV_avg*100),2)[3,c(1,6,3)], "%&"))
noquote(paste0(round(t(CV_avg*100),2)[4,c(1,6,3)], "%&"))
noquote(paste0(round(t(CV_avg*100),2)[5,c(1,6,3)], "%&"))
noquote(paste0(round(t(CV_avg*100),2)[6,c(1,6,3)], "%&"))
noquote(paste0(round(t(CV_avg*100),2)[7,c(1,6,3)], "%&"))
noquote(paste0(round(t(CV_avg*100),2)[8,c(1,6,3)], "%&"))

CV_sd = apply(CV_arr, c(1,2), sd)
# alpha = .05
# CV_CI1 = apply(CV_arr, c(1,2), quantile, probs = alpha/2)
# CV_CI2 = apply(CV_arr, c(1,2), quantile, probs = (1-alpha/2))

idx = c(1,3,6)
RES_avg = t(CV_avg[idx,])
RES_sd = t(CV_sd[idx,])
# RES_CI1 = t(CV_CI1[idx,])
# RES_CI2 = t(CV_CI2[idx,])

RES_avg =  as.data.frame(RES_avg)
RES_sd = as.data.frame(RES_sd)
# RES_CI1 = as.data.frame(RES_CI1)
# RES_CI2 = as.data.frame(RES_CI2)
colnames(RES_avg) = colnames(RES_sd) = #colnames(RES_CI1) = colnames(RES_CI2) = 
  c("GMSE :: Linearized", "GMSE :: Monte Carlo", "GMSE :: Bootstrap")

titstu = c("1 Illiterate",
           "2 Literate but no education",
           "3 Primary", 
           "4 Lower secondary",
           "5 Upper secondary",
           "6 Bachelor degree",
           "7 Master degree",
           "8 PhD level")
RES_avg$tit = RES_sd$tit = #RES_CI1$tit = RES_CI2$tit = 
  c(titstu[1:5], "6 Bachelor", "7 Master", "8 Phd")

var_n = c("GMSE :: Monte Carlo", "GMSE :: Linearized", "GMSE :: Bootstrap")#, "GMSE :: P-Bootstrap")
labels_n = c("GMSE :: Monte Carlo", "GMSE :: Linearized", "GMSE :: Bootstrap")#, "GMSE :: P-Bootstrap")

df <- reshape2::melt(RES_avg, id.var='tit')
df_var <- reshape2::melt(RES_sd, id.var='tit')
# df_CI1 <- reshape2::melt(RES_CI1, id.var='tit')
# df_CI2 <- reshape2::melt(RES_CI2, id.var='tit')
df = cbind(df, sd = df_var$value)#, CI1 = df_CI1$value, CI2 = df_CI2$value)

pos = 0
plot_1nb = ggplot(df, aes(as.character(tit), value, 
                          group=as.character(variable),
                          colour=as.character(variable))) + 
  geom_line(aes(linetype=as.character(variable), colour=as.character(variable)), linewidth = .5,
            position = position_dodge(pos), cex = 2) +
  geom_point(aes(colour=as.character(variable)), alpha = 1, cex = 1.5, position = position_dodge(pos)) +
  geom_errorbar(aes(ymin = value-2*sd/10, ymax = value+2*sd/10, colour=as.character(variable),
                    group=as.character(variable)), width = 0.2, position = position_dodge(pos))+
  #scale_linetype_manual("",values=c("solid", "dashed")) +
  scale_linetype_manual(" ",values =c("dashed", "solid", "solid","dotted"),breaks=var_n, labels=labels_n)+
  scale_colour_manual(" ",values=c("red", "#00BFC4", "#E7B800"),breaks=var_n, labels=labels_n)+
  scale_shape_manual(" ",values=c(19, 19, 19),breaks=var_n, labels=labels_n)+
  xlab('Education level') + #xlim(0, 300) + ylim(0, 50)+
  ylab("CV") +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(color = 'grey99'),
        panel.grid.minor = element_line(color = 'grey99')) +
  theme(axis.line = element_line(colour = "black"), axis.title=element_text(size=10)) +
  theme(legend.position="bottom")

plot_1nb



# one box per variety
change_n = function(x){
  rownames(x) = c(labels_n[2], "CV_hat_theta", labels_n[1],
                  "CV_hat_theta_III", "CV_hat_theta", labels_n[3], "CV_hat_theta_III")
  return(x)
}

change_cn = function(x){
  colnames(x) = titstu
  return(x)
}

library(tidyr)
CV = lapply(CV, change_n)
CV_df = data.frame(do.call("rbind", CV))
names(CV_df) = titstu
CV_df$GMSE = rep(rownames(CV[[1]]), 100)
CV_df = CV_df[CV_df$GMSE%in%rownames(CV[[1]])[c(1,3,6)],]

CV_df <- CV_df %>%
  pivot_longer("1 Illiterate":"8 PhD level", names_to = "Education", values_to = "CV")

CV_df$CV = CV_df$CV*100

ggplot(CV_df, aes(x=Education, y=CV, fill=GMSE)) + 
  geom_boxplot(alpha=0.7) +
  theme_classic() +
  #theme(legend.position="right") +
  theme(legend.position = "bottom", legend.title=element_blank()) +
  ylab("CV (%)") + ylim(0, 40) +
  xlab("Attained Education") +
  #scale_color_brewer(palette="Dark2") +
  ggtitle("Estimated CV by education class") +
  scale_fill_manual(values = c("#00BFC4", "#E7B800", "black")) #+
  #scale_fill_brewer(palette="Dark2")#+
  #scale_shape_manual(" ",values=c(19, 19, 19),breaks=var_n, labels=labels_n)+

