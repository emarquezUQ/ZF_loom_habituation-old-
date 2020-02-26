%%% i am trying to look at how many ROIs of the clusters intersect with the
%%% movement ROIs


%idx_temp=intersect(idx_rsq_Mov_cleaned3,clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned);


figure;imagesc(ZS_f20(idx_rsq_Mov_cleaned3,:)); %%% i had to resort the idx of movements. look at the script rasterplots_multisensNmotor...


figure;imagesc(ZS_f20(idx_temp,:));
figure;plot(mean(ZS_f20(idx_temp,:)));

size(idx_temp)
figure;histogram(idx_Fish_f20(idx_temp));%%% for the fish location
 

figure;scatter(BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,1),BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,2),'.');
hold on;
scatter(BrainReg_F20.ROI_temp2(idx_temp,1),BrainReg_F20.ROI_temp2(idx_temp,2),'.');



size(intersect(BrainReg_F20.PerBrainRegions.Pallium.idx,idx_temp))
figure;histogram(idx_Fish_f20(intersect(BrainReg_F20.PerBrainRegions.Pallium.idx,idx_temp)));%%% for the fish location
 

figure;scatter(BrainReg_F20.ROI_temp2(intersect(BrainReg_F20.PerBrainRegions.Pallium.idx,idx_rsq_test_f20short_cleaned),1),BrainReg_F20.ROI_temp2(intersect(BrainReg_F20.PerBrainRegions.Pallium.idx,idx_rsq_test_f20short_cleaned),2),'.');
hold on;
scatter(BrainReg_F20.ROI_temp2(intersect(BrainReg_F20.PerBrainRegions.Pallium.idx,idx_temp),1),BrainReg_F20.ROI_temp2(intersect(BrainReg_F20.PerBrainRegions.Pallium.idx,idx_temp),2),'.');




%%% to look at how many ROIs the movement overlaps with each cluster... 

size(intersect(idx_rsq_Mov_cleaned3,clust_f20_CL4_cleaned.clust_f20_CL4_1_cleaned)) %%% nonhab
size(intersect(idx_rsq_Mov_cleaned3,clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned)) %%% fasthab
size(intersect(idx_rsq_Mov_cleaned3,clust_f20_CL4_cleaned.clust_f20_CL4_6_cleaned)) %%% slopehab

size(intersect(idx_rsq_Mov_cleaned3,clust_f20_CL7_cleaned.clust_f20_CL7_1_cleaned))%%% nonhab
size(intersect(idx_rsq_Mov_cleaned3,clust_f20_CL7_cleaned.clust_f20_CL7_2_cleaned))%%% fasthab med
size(intersect(idx_rsq_Mov_cleaned3,clust_f20_CL7_cleaned.clust_f20_CL7_4_cleaned))%%% fasthab sharp
size(intersect(idx_rsq_Mov_cleaned3,clust_f20_CL7_cleaned.clust_f20_CL7_5_cleaned))%%% fasthab broad
size(intersect(idx_rsq_Mov_cleaned3,clust_f20_CL7_cleaned.clust_f20_CL7_6_cleaned))%%% slopehab



%%
idx_temp=intersect(idx_rsq_Mov_cleaned3,clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned);


figure;scatter(BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,1),BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,2),'.');
hold on;
scatter(BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,1),BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,2),'.');
hold on;
scatter(BrainReg_F20.ROI_temp2(idx_temp,1),BrainReg_F20.ROI_temp2(idx_temp,2),'.');


%%%% for 3d
figure;scatter3(BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,1),BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,2),BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,3),'.');
hold on;
scatter3(BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,1),BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,2),BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,3),'.');
hold on;
scatter3(BrainReg_F20.ROI_temp2(idx_temp,1),BrainReg_F20.ROI_temp2(idx_temp,2),BrainReg_F20.ROI_temp2(idx_temp,3),'.');

%%% it seems that they could be in the subpallium and not in the pallium...
figure;scatter(BrainReg_F20.ROI_temp2(intersect(BrainReg_F20.PerBrainRegions.Telencephalon.idx,idx_rsq_test_f20short_cleaned),1),BrainReg_F20.ROI_temp2(intersect(BrainReg_F20.PerBrainRegions.Telencephalon.idx,idx_rsq_test_f20short_cleaned),2),'.');
hold on;
scatter(BrainReg_F20.ROI_temp2(intersect(BrainReg_F20.PerBrainRegions.Telencephalon.idx,idx_temp),1),BrainReg_F20.ROI_temp2(intersect(BrainReg_F20.PerBrainRegions.Telencephalon.idx,idx_temp),2),'.');





for i=1:length(RegionList)
    regionName=RegionList{i};
    amountperRegion(:,i)=sum(ismember(PerBrainRegions.f20.(regionName).idx,idx_temp));
end
figure;bar(amountperRegion);set(gca,'xticklabel',RegionList);


%%

%%% to look at the subtypes of fasthab. 
idx_temp1=intersect(idx_rsq_Mov_cleaned3,clust_f20_CL7_cleaned.clust_f20_CL7_4_cleaned);
idx_temp2=intersect(idx_rsq_Mov_cleaned3,clust_f20_CL7_cleaned.clust_f20_CL7_2_cleaned);
idx_temp3=intersect(idx_rsq_Mov_cleaned3,clust_f20_CL7_cleaned.clust_f20_CL7_5_cleaned);

 figure;
%scatter(BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,1),BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,2),'.');
% hold on;
% scatter(BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,1),BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,2),'.');
% hold on;
scatter(BrainReg_F20.ROI_temp2(idx_temp1,1),BrainReg_F20.ROI_temp2(idx_temp1,2),'.');
hold on;
scatter(BrainReg_F20.ROI_temp2(idx_temp2,1),BrainReg_F20.ROI_temp2(idx_temp2,2),'.');
hold on;
scatter(BrainReg_F20.ROI_temp2(idx_temp3,1),BrainReg_F20.ROI_temp2(idx_temp3,2),'.');

%%

%%%% to separete the fish that responded with a movement vs the ones that
%%%% didnt. 
idx_temp=intersect(idx_rsq_Mov_cleaned3,clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned);

figure;imagesc(ZS_f20(idx_rsq_Mov_cleaned3,:)); %%% i had to resort the idx of movements. look at the script rasterplots_multisensNmotor...


figure;imagesc(ZS_f20(idx_temp,:));
figure;plot(mean(ZS_f20(idx_temp,:)));

size(idx_temp)
figure;histogram(idx_Fish_f20(idx_temp));%%% for the fish location


fish_1stLoom_resp=unique(idx_Fish_f20(idx_temp));

fishF20=unique(idx_Fish_f20);

fish_non1stLoom_resp=find(~ismember(fishF20,fish_1stLoom_resp));

fish_non1stLoom_resp=fishF20(fish_non1stLoom_resp);




idx_fish_1stLoom_resp=[];

for i=1:length(fish_1stLoom_resp)
idx_fish_1stLoom_resp_temp=find(idx_Fish_f20==fish_1stLoom_resp(i));

idx_fish_1stLoom_resp=vertcat(idx_fish_1stLoom_resp,idx_fish_1stLoom_resp_temp);

end

clear idx_fish_1stLoom_resp_temp



idx_fish_non1stLoom_resp=[];

for i=1:length(fish_non1stLoom_resp)
idx_fish_1stLoom_resp_temp=find(idx_Fish_f20==fish_non1stLoom_resp(i));

idx_fish_non1stLoom_resp=vertcat(idx_fish_non1stLoom_resp,idx_fish_1stLoom_resp_temp);

end

clear idx_fish_non1stLoom_resp_temp

%%% to check if it worked. 
figure;imagesc(ZS_f20(intersect(idx_rsq_Mov_cleaned3,idx_fish_1stLoom_resp),:)); %%% i had to resort the idx of movements. look at the script rasterplots_multisensNmotor...


size(intersect(idx_rsq_Mov_cleaned3,idx_fish_1stLoom_resp)) %%% 1stLoom_resp
size(intersect(idx_rsq_Mov_cleaned3,idx_fish_non1stLoom_resp)) %%% non1stLoom_resp

idx_temp1=intersect(idx_rsq_Mov_cleaned3,idx_fish_1stLoom_resp);
idx_temp2=intersect(idx_rsq_Mov_cleaned3,idx_fish_non1stLoom_resp);

 figure;
%scatter(BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,1),BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,2),'.');
% hold on;
% scatter(BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,1),BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,2),'.');
% hold on;
scatter(BrainReg_F20.ROI_temp2(idx_temp1,1),BrainReg_F20.ROI_temp2(idx_temp1,2),'.');
hold on;
scatter(BrainReg_F20.ROI_temp2(idx_temp2,1),BrainReg_F20.ROI_temp2(idx_temp2,2),'.');



%%%% for 3D
figure;
% figure;scatter3(BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,1),BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,2),BrainReg_F20.ROI_temp2(idx_rsq_test_f20short_cleaned,3),'.');
% hold on;
% scatter3(BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,1),BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,2),BrainReg_F20.ROI_temp2(idx_rsq_Mov_cleaned3,3),'.');
% hold on;
scatter3(BrainReg_F20.ROI_temp2(idx_temp1,1),BrainReg_F20.ROI_temp2(idx_temp1,2),BrainReg_F20.ROI_temp2(idx_temp1,3),'.');
hold on;
scatter3(BrainReg_F20.ROI_temp2(idx_temp2,1),BrainReg_F20.ROI_temp2(idx_temp2,2),BrainReg_F20.ROI_temp2(idx_temp2,3),'.');
