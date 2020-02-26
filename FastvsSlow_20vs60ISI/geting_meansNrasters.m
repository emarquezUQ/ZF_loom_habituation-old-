
cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab

%load('final_F60_step1.mat','ZS_f60','rawregressF60','idx_rsq_test_f60short','High_corr_Nb_f60','High_corr_Nb_f60_short','gooodmaps');

load('final_F60_step1_2.mat','ZS_f60','rawregressF20','idx_rsq_test_f60short','High_corr_Nb_f60','High_corr_Nb_f60_short','gooodmaps');



%%% this is just to get the means and make the rasterplots for figures 1
%%% and 2.

%%
%%% to get the means of the 4 main clusters and the 3 kinds of fast hab
%%% clusters for figures. 


%%%% note!!!!! i need to change the values of j between Fs and S and the
%%%% clusters are not in the same order. 


%%% first for CL4
for j=gooodmaps
idx_temp=find(High_corr_Nb_f60_short==j);

if j==2
mean_CL4_f60.fasthab=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),:));
elseif j==6 
mean_CL4_f60.slopehab=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),:));
elseif j==1
mean_CL4_f60.nonhab=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),:));   
else
mean_CL4_f60.inhib=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),:));
end

end

%figure;plot(mean_CL4_f60.slopehab);
%figure;plot(mean(ZS_f60(idx_rsq_test_f60short(idx_temp),:)));

for j=1:size(rawregressF20,1)
idx_temp=find(High_corr_Nb_f60==j);

if j==2
mean_CL7_f60.fasthab_med=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),:));
mean_CL7_f60_short.fasthab_med=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),55:180));
elseif j==3
mean_CL7_f60.sound=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),:));  
mean_CL7_f60_short.sound=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),55:180));  
elseif j==4
mean_CL7_f60.fasthab_sharp=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),:));
mean_CL7_f60_short.fasthab_sharp=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),55:180));
elseif j==5
mean_CL7_f60.fasthab_broad=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),:));
mean_CL7_f60_short.fasthab_broad=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),55:180));
elseif j==6 
mean_CL7_f60.slopehab=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),:));
mean_CL7_f60_short.slopehab=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),55:180));
elseif j==1
mean_CL7_f60.nonhab=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),:)); 
mean_CL7_f60_short.nonhab=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),55:180)); 
else
mean_CL7_f60.inhib=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),:));
mean_CL7_f60_short.inhib=mean(ZS_f60(idx_rsq_test_f60short(idx_temp),55:180));
end

end

%figure;plot(mean_CL7_f60.fasthab_sharp);
%figure;plot(mean_CL7_f60_short.slopehab);

save('means_F60_CL4n7.mat','mean_CL4_f60','mean_CL7_f60','mean_CL7_f60_short','-v7.3');



%%% now getting the raster plots of the 4 main clusters

for j=gooodmaps
idx_temp=find(High_corr_Nb_f60_short==j);

if j==2
figure;imagesc(ZS_f60(idx_rsq_test_f60short(idx_temp),:),[min(min(mean_CL7_f60.inhib))-0.5 max(max(mean_CL7_f60.fasthab_sharp))]);colormap hot
saveas(gcf,'rasterplot_f60_CL4_fasthab.tif');

elseif j==6 
 figure;imagesc(ZS_f60(idx_rsq_test_f60short(idx_temp),:),[min(min(mean_CL7_f60.inhib))-0.5 max(max(mean_CL7_f60.fasthab_sharp))]);colormap hot
saveas(gcf,'rasterplot_f60_CL4_slopehab.tif');   

elseif j==1
  figure;imagesc(ZS_f60(idx_rsq_test_f60short(idx_temp),:),[min(min(mean_CL7_f60.inhib))-0.5 max(max(mean_CL7_f60.fasthab_sharp))]);colormap hot
saveas(gcf,'rasterplot_f60_CL4_nonhab.tif');  

else
  figure;imagesc(ZS_f60(idx_rsq_test_f60short(idx_temp),:),[min(min(mean_CL7_f60.inhib))-0.5 max(max(mean_CL7_f60.fasthab_sharp))]);colormap hot; 
saveas(gcf,'rasterplot_f60_CL4_inhib.tif');    
    
figure;imagesc(ZS_f60(idx_rsq_test_f60short(idx_temp),:),[min(min(mean_CL7_f60.inhib))-0.5 max(max(mean_CL7_f60.fasthab_sharp))]);colormap hot; colorbar
saveas(gcf,'rasterplot_f60_CL4_inhib_wbar.tif');   

end

end


%%% to checkit worked
names=fieldnames(mean_CL4_f60);
counter=1;
figure;
for i=1:length(names)-1
subplot(4,2,counter);plot(mean_CL4_f60.(names{i}));

counter=counter+1;
end


%%% to checkit worked
names=fieldnames(mean_CL4_f60);

figure;
for i=1:length(names)-1
plot(mean_CL4_f60.(names{i}));
hold on;

end

