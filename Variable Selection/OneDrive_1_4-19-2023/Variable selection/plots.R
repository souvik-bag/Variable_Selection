bdata_boxplot21 <- cbind(r_lasso21$brier,r_alasso21$brier,r_sparse21$brier,r_enet21$brier,r_best21$brier,r_scad21$brier,
                         r_mcp21$brier,r_sis21$brier,r_isis21$brier,r_sivs21$brier,r_boruta21$brier,r_vsurf21$brier,
                         r_rrf21$brier,r_pimp21$brier,r_nta21$brier,"Method")
colnames(bdata_boxplot21) = c("LASSO","ALASSO","SparseStep","ElasticNet","BestSubset","SCAD","MCP","SIS","ISIS","SIVS","Boruta","VSURF", "RRF",
                              "PIMP","NTA","Type")
bdata_boxplot21 <- as.data.frame(bdata_boxplot21)
bdata_boxplot_long21 <- melt(bdata_boxplot21,id = "Type" )  
#sapply(bdata_boxplot_long21,"class")
#bdata_boxplot_long21$value <- as.numeric(bdata_boxplot_long21$value)
colnames(bdata_boxplot_long21) <- c("Type","Method","BrierScore")

sapply(bdata_boxplot_long21,"class")
bdata_boxplot_long21$`BrierScore` <- round(as.numeric(bdata_boxplot_long21$`BrierScore`),3)
b <- ggplot(bdata_boxplot_long21, aes(x = Method, y = BrierScore)) +  # ggplot function
  geom_boxplot() +
  theme_minimal() +theme(axis.text = element_text(size = 15),
                         axis.title = element_text(size = 15),
                         legend.text = element_text(size = 15),
                         legend.title = element_text( size = 15),
                         axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
                         legend.position = "bottom")

bdata_boxplot22 <- cbind(r_lasso22$brier,r_alasso22$brier,r_sparse22$brier,r_enet22$brier,r_best22$brier,r_scad22$brier,
                         r_mcp22$brier,r_sis22$brier,r_isis22$brier,r_sivs22$brier,r_boruta22$brier,r_vsurf22$brier,
                         r_rrf22$brier,r_pimp22$brier,r_nta22$brier,"Method")
colnames(bdata_boxplot22) = c("LASSO","ALASSO","SparseStep","ElasticNet","BestSubset","SCAD","MCP","SIS","ISIS","SIVS","Boruta","VSURF", "RRF",
                              "PIMP","NTA","Type")
bdata_boxplot22 <- as.data.frame(bdata_boxplot22)
bdata_boxplot_long22 <- melt(bdata_boxplot22,id = "Type" )  
#sapply(bdata_boxplot_long22,"class")
#bdata_boxplot_long22$value <- as.numeric(bdata_boxplot_long22$value)
colnames(bdata_boxplot_long22) <- c("Type","Method","BrierScore")

sapply(bdata_boxplot_long22,"class")
bdata_boxplot_long22$`BrierScore` <- round(as.numeric(bdata_boxplot_long22$`BrierScore`),3)
b1<- ggplot(bdata_boxplot_long22, aes(x = Method, y = BrierScore)) +  # ggplot function
  geom_boxplot() +
  theme_minimal() +theme(axis.text = element_text(size = 15),
                         axis.title = element_text(size = 15),
                         legend.text = element_text(size = 15),
                         legend.title = element_text( size = 15),
                         axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
                         legend.position = "bottom")

z1 <- ggpubr::ggarrange(b,b1,ncol = 1, nrow = 2)

ggsave(plot = z1,filename = "selected_fig3.eps",
       width = 300,
       height = 400,
       units = "mm",
       dpi = 800)