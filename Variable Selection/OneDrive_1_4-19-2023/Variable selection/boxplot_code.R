plot22 <- data.frame(r_lasso22$selected,r_alasso22$selected,
                     r_enet22$selected,
                     r_sparse22$selected,r_best22$selected,r_scad22$selected
                     ,r_mcp22$selected,r_sis22$selected,r_isis22$selected,r_sivs22$selected,r_boruta22$selected,r_vsurf22$selected,r_rrf22$selected,
                     r_pimp22$selected,r_nta22$selected)


par(mar=c(7,5,1,1))
boxplot(plot22, names = c("LASSO","Ada_LASSO","Elastic_net","Sparse_step","Best_subset","SCAD","MCP","SIS","ISIS","SIVS","BORUTA",
                         "VSURF","RRF","PIMP","NTA"), las = 2, ylab = "Number of Selected Variables",ylim = c(0,60))
geom_boxplot(plot22)

names_methods = c("LASSO","Adaptive_LASSO","Elastic_net","Sparse_step","Best_subset","SCAD","MCP","SIS","ISIS","SIVS","BORUTA",
                  "VSURF","RRF","PIMP","NTA")

boxplot(plot21~names_methods)


# r_lasso21 = data.frame(r_lasso21)
# r_enet21 = data.frame(r_enet21)
# r_alasso21 = data.frame( r_alasso21)
# r_sivs21 = data.frame(r_sivs21)
# r_sparse21 = data.frame(r_sparse21)
# r_best21 = data.frame(r_best21)
# r_scad21 = data.frame(r_scad21)
# r_mcp21 = data.frame(r_mcp21)
# r_sis21 = data.frame(r_sis21)
# r_isis21 = data.frame(r_isis21)
# r_boruta21 = data.frame(r_boruta21)
# r_vsurf21 = data.frame(r_vsurf21)
# r_rrf21 = data.frame(r_rrf21)
# r_pimp21 = data.frame(r_pimp21)
# r_nta21 =  data.frame(r_nta21)
# 
# 
# 
# colnames(r_lasso21) = colnames(r_enet21) = colnames(r_alasso21) =
#   colnames(r_sivs21) =colnames(r_sparse21) =
#   colnames(r_best21) =colnames(r_scad21) =colnames(r_mcp21) =
#   colnames(r_sis21) =colnames(r_isis21) =colnames(r_boruta21) =
#   colnames(r_vsurf21) =colnames(r_rrf21) =colnames(r_pimp21) = colnames(r_nta21) = 
#   c("selected","Imp%","unimp%","mse","mae","misclas","pre","reca","time","brier")
