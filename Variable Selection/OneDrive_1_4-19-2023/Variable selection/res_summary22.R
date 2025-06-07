
r_lasso22 = data.frame(r_lasso22)
r_enet22 = data.frame(r_enet22)
r_alasso22 = data.frame( r_alasso22)
r_sivs22 = data.frame(r_sivs22)
r_sparse22 = data.frame(r_sparse22)
r_best22 = data.frame(r_best22)
r_scad22 = data.frame(r_scad22)
r_mcp22 = data.frame(r_mcp22)
r_sis22 = data.frame(r_sis22)
r_isis22 = data.frame(r_isis22)
r_boruta22 = data.frame(r_boruta22)
r_vsurf22 = data.frame(r_vsurf22)
r_rrf22 = data.frame(r_rrf22)
r_pimp22 = data.frame(r_pimp22)
r_nta22 =  data.frame(r_nta22)
colnames(r_lasso22) = colnames(r_enet22) = colnames(r_alasso22) =
  colnames(r_sivs22) =colnames(r_sparse22) =
  colnames(r_best22) =colnames(r_scad22) =colnames(r_mcp22) =
  colnames(r_sis22) =colnames(r_isis22) =colnames(r_boruta22) =
  colnames(r_vsurf22) =colnames(r_rrf22) =colnames(r_pimp22) = colnames(r_nta22) = 
  c("selected","Imp%","unimp%","mse","mae","misclas","pre","reca","time","brier")

#    sim = sim + 1


res_summary22 = data.frame(cbind(colMeans(r_lasso22,na.rm = T),colMeans(r_alasso22,na.rm = T),
                               colMeans(r_sparse22,na.rm = T),colMeans(r_enet22,na.rm = T),
                               colMeans(r_best22,na.rm = T),
                               colMeans(r_scad22,na.rm = T),colMeans(r_mcp22,na.rm = T),
                               colMeans(r_sis22,na.rm = T),colMeans(r_isis22,na.rm = T),
                               colMeans(r_sivs22,na.rm = T),
                               colMeans(r_boruta22,na.rm = T),colMeans(r_vsurf22,na.rm = T),
                               colMeans(r_rrf22,na.rm = T),colMeans(r_pimp22,na.rm = T),
                               colMeans(r_nta22,na.rm = T)))
colnames(res_summary22) = c("lasso","alasso","sparse","enet","best","scad","mcp","sis","isis","sivs","boruta","vsurf", "rrf",
                          "pimp","nta")
res_summary22 %>%
  mutate_if(is.numeric,round,digits = 3) %>%
  View()
