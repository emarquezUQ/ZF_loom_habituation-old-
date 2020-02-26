load('final_S60_step1.mat','ZS_s60','rawregressS20','idx_rsq_test_s60short','High_corr_Nb_s60','High_corr_Nb_s60_short','gooodmaps');





%%% this is just to get the means and make the rasterplots for figures 1
%%% and 2.

%%
%%% to get the means of the 4 main clusters and the 3 kinds of fast hab
%%% clusters for figures. 






%%% first for CL4
for j=gooodmaps
idx_temp=find(High_corr_Nb_s60_short==j);

if j==1
mean_CL4_s60.fasthab=mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60));
elseif j==5 
mean_CL4_s60.slopehab=mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60));
elseif j==3
mean_CL4_s60.nonhab=mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60));   
else
mean_CL4_s60.inhib=mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60));
end

end

%figure;plot(mean_CL4_s60.inhib);
% j=4
% idx_temp=find(High_corr_Nb_s60_short==j);
%figure;plot(mean(ZS_s60(idx_rsq_test_s60short(idx_temp),:)));

for j=1:size(rawregressS20,1)
idx_temp=find(High_corr_Nb_s60==j);

if j==4
mean_CL7_s60.fasthab_med=mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60));
mean_CL7_s60_short.fasthab_med=mean_CL7_s60.fasthab_med(1,55:180);
elseif j==6
mean_CL7_s60.sound=mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60));  
mean_CL7_s60_short.sound=mean_CL7_s60.sound(1,55:180);  
elseif j==2
mean_CL7_s60.fasthab_sharp=mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60));
mean_CL7_s60_short.fasthab_sharp=mean_CL7_s60.fasthab_sharp(1,55:180);
elseif j==1
mean_CL7_s60.fasthab_broad=mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60));
mean_CL7_s60_short.fasthab_broad=mean_CL7_s60.fasthab_broad(1,55:180);
elseif j==5
mean_CL7_s60.slopehab=mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60));
mean_CL7_s60_short.slopehab=mean_CL7_s60.slopehab(1,55:180);
elseif j==3
mean_CL7_s60.nonhab=mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60)); 
mean_CL7_s60_short.nonhab=mean_CL7_s60.nonhab(1,55:180);
else
mean_CL7_s60.inhib=mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60));
mean_CL7_s60_short.inhib=mean_CL7_s60.inhib(1,55:180);
end

end

%figure;plot(mean_CL7_s60.nonhab);
%figure;plot(mean_CL7_s60_short.slopehab);

save('means_S60_CL4n7.mat','mean_CL4_s60','mean_CL7_s60','mean_CL7_s60_short','-v7.3');



%%% now getting the raster plots of the 4 main clusters

for j=gooodmaps
idx_temp=find(High_corr_Nb_s60_short==j);

if j==1
figure;imagesc(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),[min(min(mean_CL7_s60.inhib))-0.5 max(max(mean_CL7_s60.fasthab_sharp))]);colormap hot
saveas(gcf,'rasterplot_s60_CL4_fasthab.tif');

elseif j==5 
 figure;imagesc(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),[min(min(mean_CL7_s60.inhib))-0.5 max(max(mean_CL7_s60.fasthab_sharp))]);colormap hot
saveas(gcf,'rasterplot_s60_CL4_slopehab.tif');   

elseif j==3
  figure;imagesc(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),[min(min(mean_CL7_s60.inhib))-0.5 max(max(mean_CL7_s60.fasthab_sharp))]);colormap hot
saveas(gcf,'rasterplot_s60_CL4_nonhab.tif');  

else
  figure;imagesc(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),[min(min(mean_CL7_s60.inhib))-0.5 max(max(mean_CL7_s60.fasthab_sharp))]);colormap hot; 
saveas(gcf,'rasterplot_s60_CL4_inhib.tif');    
    
figure;imagesc(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),[min(min(mean_CL7_s60.inhib))-0.5 max(max(mean_CL7_s60.fasthab_sharp))]);colormap hot; colorbar
saveas(gcf,'rasterplot_s60_CL4_inhib_wbar.tif');   

end

end
