load('All_More_BrainReg2.mat');
f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');


%%%%%% this script is to make a mask of the medial cerebellum and medial
%%%%%% hindbrain to have a closer look at their responses. 

%%% I will start with the medial hindbrain. I will first get together all
%%% the hindbrain ROIs of the 4 datasets together. 

fast_hindB_idx_f20=intersect(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,PerBrainRegions.f20.Hindbrain.idx);
fast_cereb_idx_f20=intersect(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,PerBrainRegions.f20.Cerebellum.idx);

figure;
scatter3(ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,2),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,3),'.');
hold on;
scatter3(ROI_temp2.f20(fast_cereb_idx_f20,1),ROI_temp2.f20(fast_cereb_idx_f20,2),ROI_temp2.f20(fast_cereb_idx_f20,3),'.');
hold on;
scatter3(ROI_temp2.f20(fast_hindB_idx_f20,1),ROI_temp2.f20(fast_hindB_idx_f20,2),ROI_temp2.f20(fast_hindB_idx_f20,3),'.');

%%
%%%% the hindbrain mask we did includes more than the medial strip that we
%%%% are interested on... so I will try with a mask from rombomeres 2-7.

figure;
scatter3(ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,2),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,3),'.');
hold on;
scatter3(ROI_temp2.f20(fast_cereb_idx_f20,1),ROI_temp2.f20(fast_cereb_idx_f20,2),ROI_temp2.f20(fast_cereb_idx_f20,3),'.');
hold on;
scatter3(ROI_temp2.f20(fast_hindB_idx_f20,1),ROI_temp2.f20(fast_hindB_idx_f20,2),ROI_temp2.f20(fast_hindB_idx_f20,3),'.');
hold on;
scatter3(Zbrain_Masks{219,3}(:,1),Zbrain_Masks{219,3}(:,2),Zbrain_Masks{219,3}(:,3),'.');

%%

%%%% here I am getting the rhombomeres' masks and using them to get the fast hab ROIs in each dataset 
hindbrain_masks=load('All_More_BrainReg3.mat');

Rhombs_list=hindbrain_masks.RegionList(106:111);
regnum=[106:111];

figure;
fast_hindB_idx_f20=[];
for i=1:length(Rhombs_list)
fast_hindB_idx_f20_temp=intersect(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_2_cleaned,hindbrain_masks.PerBrainRegions.f20.(strcat('reg',num2str(regnum(i)))).idx);
fast_hindB_idx_f20=vertcat(fast_hindB_idx_f20,fast_hindB_idx_f20_temp);

scatter3(hindbrain_masks.ROI_temp2.f20(fast_hindB_idx_f20_temp,1),hindbrain_masks.ROI_temp2.f20(fast_hindB_idx_f20_temp,2),hindbrain_masks.ROI_temp2.f20(fast_hindB_idx_f20_temp,3),'.');
hold on;   
    
end
figure;
scatter3(hindbrain_masks.ROI_temp2.f20(fast_hindB_idx_f20,1),hindbrain_masks.ROI_temp2.f20(fast_hindB_idx_f20,2),hindbrain_masks.ROI_temp2.f20(fast_hindB_idx_f20,3),'.');

figure;
fast_hindB_idx_f60=[];
for i=1:length(Rhombs_list)
fast_hindB_idx_f60_temp=intersect(f60_cleaned_idxs.clust_f60_CL4_cleaned.clust_f60_CL4_2_cleaned,hindbrain_masks.PerBrainRegions.f60.(strcat('reg',num2str(regnum(i)))).idx);
fast_hindB_idx_f60=vertcat(fast_hindB_idx_f60,fast_hindB_idx_f60_temp);

scatter3(hindbrain_masks.ROI_temp2.f60(fast_hindB_idx_f60_temp,1),hindbrain_masks.ROI_temp2.f60(fast_hindB_idx_f60_temp,2),hindbrain_masks.ROI_temp2.f60(fast_hindB_idx_f60_temp,3),'.');
hold on;   
    
end
figure;
scatter3(hindbrain_masks.ROI_temp2.f60(fast_hindB_idx_f60,1),hindbrain_masks.ROI_temp2.f60(fast_hindB_idx_f60,2),hindbrain_masks.ROI_temp2.f60(fast_hindB_idx_f60,3),'.');


figure;
fast_hindB_idx_s20=[];
for i=1:length(Rhombs_list)
fast_hindB_idx_s20_temp=intersect(s20_cleaned_idxs.clust_s20_CL4_cleaned.clust_s20_CL4_1_cleaned,hindbrain_masks.PerBrainRegions.s20.(strcat('reg',num2str(regnum(i)))).idx);
fast_hindB_idx_s20=vertcat(fast_hindB_idx_s20,fast_hindB_idx_s20_temp);

scatter3(hindbrain_masks.ROI_temp2.s20(fast_hindB_idx_s20_temp,1),hindbrain_masks.ROI_temp2.s20(fast_hindB_idx_s20_temp,2),hindbrain_masks.ROI_temp2.s20(fast_hindB_idx_s20_temp,3),'.');
hold on;   
    
end
figure;
scatter3(hindbrain_masks.ROI_temp2.s20(fast_hindB_idx_s20,1),hindbrain_masks.ROI_temp2.s20(fast_hindB_idx_s20,2),hindbrain_masks.ROI_temp2.s20(fast_hindB_idx_s20,3),'.');


figure;
fast_hindB_idx_s60=[];
for i=1:length(Rhombs_list)
fast_hindB_idx_s60_temp=intersect(s60_cleaned_idxs.clust_s60_CL4_cleaned.clust_s60_CL4_1_cleaned,hindbrain_masks.PerBrainRegions.s60.(strcat('reg',num2str(regnum(i)))).idx);
fast_hindB_idx_s60=vertcat(fast_hindB_idx_s60,fast_hindB_idx_s60_temp);

scatter3(hindbrain_masks.ROI_temp2.s60(fast_hindB_idx_s60_temp,1),hindbrain_masks.ROI_temp2.s60(fast_hindB_idx_s60_temp,2),hindbrain_masks.ROI_temp2.s60(fast_hindB_idx_s60_temp,3),'.');
hold on;   
    
end
figure;
scatter3(hindbrain_masks.ROI_temp2.s60(fast_hindB_idx_s60,1),hindbrain_masks.ROI_temp2.s60(fast_hindB_idx_s60,2),hindbrain_masks.ROI_temp2.s60(fast_hindB_idx_s60,3),'.');


figure;
scatter3(hindbrain_masks.ROI_temp2.f20(fast_hindB_idx_f20,1),hindbrain_masks.ROI_temp2.f20(fast_hindB_idx_f20,2),hindbrain_masks.ROI_temp2.f20(fast_hindB_idx_f20,3),'.');
hold on;
scatter3(hindbrain_masks.ROI_temp2.f60(fast_hindB_idx_f60,1),hindbrain_masks.ROI_temp2.f60(fast_hindB_idx_f60,2),hindbrain_masks.ROI_temp2.f60(fast_hindB_idx_f60,3),'.');
hold on;
scatter3(hindbrain_masks.ROI_temp2.s20(fast_hindB_idx_s20,1),hindbrain_masks.ROI_temp2.s20(fast_hindB_idx_s20,2),hindbrain_masks.ROI_temp2.s20(fast_hindB_idx_s20,3),'.');
hold on;
scatter3(hindbrain_masks.ROI_temp2.s60(fast_hindB_idx_s60,1),hindbrain_masks.ROI_temp2.s60(fast_hindB_idx_s60,2),hindbrain_masks.ROI_temp2.s60(fast_hindB_idx_s60,3),'.');

fast_hindB_coor_all=[];
fast_hindB_coor_all=hindbrain_masks.ROI_temp2.f20(fast_hindB_idx_f20,:);
fast_hindB_coor_all=vertcat(fast_hindB_coor_all,hindbrain_masks.ROI_temp2.f60(fast_hindB_idx_f60,:));
fast_hindB_coor_all=vertcat(fast_hindB_coor_all,hindbrain_masks.ROI_temp2.s20(fast_hindB_idx_s20,:));
fast_hindB_coor_all=vertcat(fast_hindB_coor_all,hindbrain_masks.ROI_temp2.s60(fast_hindB_idx_s60,:));

figure;
scatter3(fast_hindB_coor_all(:,1),fast_hindB_coor_all(:,2),fast_hindB_coor_all(:,3),'.');

%%% I have decided that the selection of the medial hindbrain ROIs should
%%% be from 375(left end)-305(putative midline)-235(right end) in y coor

fast_hindB_coor_all_good=find(fast_hindB_coor_all(:,2)<375 & fast_hindB_coor_all(:,2)>235);
fast_hindB_coor_all_good=fast_hindB_coor_all(fast_hindB_coor_all_good,:);

figure;
scatter3(fast_hindB_coor_all_good(:,1),fast_hindB_coor_all_good(:,2),fast_hindB_coor_all_good(:,3),'.');



%%
%%% this part is to make a polygone

%k=boundary(fast_hindB_coor_all_good,0.5); %%% i am using a shrink factor (0-1) of 0.1.  

k=convhull(fast_hindB_coor_all_good);

%%% this is to make a Patch 3-D polygon with the boundarys
h=trisurf(k,fast_hindB_coor_all_good(:,1),fast_hindB_coor_all_good(:,2),fast_hindB_coor_all_good(:,3),'FaceColor','red','FaceAlpha',0.1);



%%% to reduce the faces of the patch while keeping the same shape. I am
%%% doing it to 20%
h3=reducepatch(h,.2,'verbose');

patch(h3,'EdgeColor','none','FaceAlpha',0.1);

scatter3(h3.vertices(:,1),h3.vertices(:,2),h3.vertices(:,3));

%%

tri_test = delaunayn([h3.vertices(:,1) h3.vertices(:,2) h3.vertices(:,3)]); % Generate delaunay triangulization
tn_test = tsearchn([h3.vertices(:,1) h3.vertices(:,2) h3.vertices(:,3)], tri_test, fast_hindB_coor_all_good); % Determine which triangle point is within
IsInside_test = ~isnan(tn_test); % Convert to logical vector

patch(h3,'EdgeColor','none','FaceAlpha',0.1);
hold on;
scatter3(fast_hindB_coor_all_good(:,1),fast_hindB_coor_all_good(:,2),fast_hindB_coor_all_good(:,3),'.');

notInsdie=find(IsInside_test==0);

patch(h3,'EdgeColor','none','FaceAlpha',0.1);
hold on;
scatter3(fast_hindB_coor_all_good(:,1),fast_hindB_coor_all_good(:,2),fast_hindB_coor_all_good(:,3),'.');
hold on;
scatter3(fast_hindB_coor_all_good(notInsdie,1),fast_hindB_coor_all_good(notInsdie,2),fast_hindB_coor_all_good(notInsdie,3),'.');


%%%% now tryting with the zbrain mask of the rhomboencephalon

tri_test = delaunayn([h3.vertices(:,1) h3.vertices(:,2) h3.vertices(:,3)]); % Generate delaunay triangulization
tn_test = tsearchn([h3.vertices(:,1) h3.vertices(:,2) h3.vertices(:,3)], tri_test, Zbrain_Masks{259,3}); % Determine which triangle point is within
IsInside_test = ~isnan(tn_test); % Convert to logical vector

IsInside_good=find(IsInside_test);

scatter3(Zbrain_Masks{259,3}(IsInside_good,1),Zbrain_Masks{259,3}(IsInside_good,2),Zbrain_Masks{259,3}(IsInside_good,3),'.');

%%% it seems to work!!!! what if I do it with non convex polygon?
%%% it seems that the delaunayn function doesnt respect the concave
%%% boundaries so making the polygon with boundary gives me the same results than convhull


hindB_mask=Zbrain_Masks{259,3}(IsInside_good,:);
scatter3(hindB_mask(:,1),hindB_mask(:,2),hindB_mask(:,3),'.');

%% now for the cerebellum/hindbrain cluster of nonhab/slopehab ROIs


%%%% this part of the script is for the cluster of slope and nonhab neurons
%%%% in the medial cerebellum/hindbran. 

%%%%% to visualize nonhab cells in cerebellum and hindbrain in the 4 datasets

non_hindB_idx_f20=intersect(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_1_cleaned,PerBrainRegions.f20.Hindbrain.idx);
non_cereb_idx_f20=intersect(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_1_cleaned,PerBrainRegions.f20.Cerebellum.idx);

figure;
scatter3(ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,2),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,3),'.');
hold on;
scatter3(ROI_temp2.f20(non_cereb_idx_f20,1),ROI_temp2.f20(non_cereb_idx_f20,2),ROI_temp2.f20(non_cereb_idx_f20,3),'.');
hold on;
scatter3(ROI_temp2.f20(non_hindB_idx_f20,1),ROI_temp2.f20(non_hindB_idx_f20,2),ROI_temp2.f20(non_hindB_idx_f20,3),'.');


%%

non_hindB_idx_f60=intersect(f60_cleaned_idxs.clust_f60_CL4_cleaned.clust_f60_CL4_1_cleaned,PerBrainRegions.f60.Hindbrain.idx);
non_cereb_idx_f60=intersect(f60_cleaned_idxs.clust_f60_CL4_cleaned.clust_f60_CL4_1_cleaned,PerBrainRegions.f60.Cerebellum.idx);

figure;
scatter3(ROI_temp2.f60(f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned,1),ROI_temp2.f60(f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned,2),ROI_temp2.f60(f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned,3),'.');
hold on;
scatter3(ROI_temp2.f60(non_cereb_idx_f60,1),ROI_temp2.f60(non_cereb_idx_f60,2),ROI_temp2.f60(non_cereb_idx_f60,3),'.');
hold on;
scatter3(ROI_temp2.f60(non_hindB_idx_f60,1),ROI_temp2.f60(non_hindB_idx_f60,2),ROI_temp2.f60(non_hindB_idx_f60,3),'.');

%%

non_hindB_idx_s20=intersect(s20_cleaned_idxs.clust_s20_CL4_cleaned.clust_s20_CL4_3_cleaned,PerBrainRegions.s20.Hindbrain.idx);
non_cereb_idx_s20=intersect(s20_cleaned_idxs.clust_s20_CL4_cleaned.clust_s20_CL4_3_cleaned,PerBrainRegions.s20.Cerebellum.idx);

figure;
scatter3(ROI_temp2.s20(s20_cleaned_idxs.idx_rsq_test_s20short_cleaned,1),ROI_temp2.s20(s20_cleaned_idxs.idx_rsq_test_s20short_cleaned,2),ROI_temp2.s20(s20_cleaned_idxs.idx_rsq_test_s20short_cleaned,3),'.');
hold on;
scatter3(ROI_temp2.s20(non_cereb_idx_s20,1),ROI_temp2.s20(non_cereb_idx_s20,2),ROI_temp2.s20(non_cereb_idx_s20,3),'.');
hold on;
scatter3(ROI_temp2.s20(non_hindB_idx_s20,1),ROI_temp2.s20(non_hindB_idx_s20,2),ROI_temp2.s20(non_hindB_idx_s20,3),'.');


%%
non_hindB_idx_s60=intersect(s60_cleaned_idxs.clust_s60_CL4_cleaned.clust_s60_CL4_3_cleaned,PerBrainRegions.s60.Hindbrain.idx);
non_cereb_idx_s60=intersect(s60_cleaned_idxs.clust_s60_CL4_cleaned.clust_s60_CL4_3_cleaned,PerBrainRegions.s60.Cerebellum.idx);

figure;
scatter3(ROI_temp2.s60(s60_cleaned_idxs.idx_rsq_test_s60short_cleaned,1),ROI_temp2.s60(s60_cleaned_idxs.idx_rsq_test_s60short_cleaned,2),ROI_temp2.s60(s60_cleaned_idxs.idx_rsq_test_s60short_cleaned,3),'.');
hold on;
scatter3(ROI_temp2.s60(non_cereb_idx_s60,1),ROI_temp2.s60(non_cereb_idx_s60,2),ROI_temp2.s60(non_cereb_idx_s60,3),'.');
hold on;
scatter3(ROI_temp2.s60(non_hindB_idx_s60,1),ROI_temp2.s60(non_hindB_idx_s60,2),ROI_temp2.s60(non_hindB_idx_s60,3),'.');

%%

%%%%% to visualize slope cells in cerebellum and hindbrain in the 4 datasets

slope_hindB_idx_f20=intersect(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_6_cleaned,PerBrainRegions.f20.Hindbrain.idx);
slope_cereb_idx_f20=intersect(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_6_cleaned,PerBrainRegions.f20.Cerebellum.idx);

figure;
scatter3(ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,2),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,3),'.');
hold on;
scatter3(ROI_temp2.f20(slope_cereb_idx_f20,1),ROI_temp2.f20(slope_cereb_idx_f20,2),ROI_temp2.f20(slope_cereb_idx_f20,3),'.');
hold on;
scatter3(ROI_temp2.f20(slope_hindB_idx_f20,1),ROI_temp2.f20(slope_hindB_idx_f20,2),ROI_temp2.f20(slope_hindB_idx_f20,3),'.');


%%

slope_hindB_idx_f60=intersect(f60_cleaned_idxs.clust_f60_CL4_cleaned.clust_f60_CL4_6_cleaned,PerBrainRegions.f60.Hindbrain.idx);
slope_cereb_idx_f60=intersect(f60_cleaned_idxs.clust_f60_CL4_cleaned.clust_f60_CL4_6_cleaned,PerBrainRegions.f60.Cerebellum.idx);

figure;
scatter3(ROI_temp2.f60(f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned,1),ROI_temp2.f60(f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned,2),ROI_temp2.f60(f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned,3),'.');
hold on;
scatter3(ROI_temp2.f60(slope_cereb_idx_f60,1),ROI_temp2.f60(slope_cereb_idx_f60,2),ROI_temp2.f60(slope_cereb_idx_f60,3),'.');
hold on;
scatter3(ROI_temp2.f60(slope_hindB_idx_f60,1),ROI_temp2.f60(slope_hindB_idx_f60,2),ROI_temp2.f60(slope_hindB_idx_f60,3),'.');

%%

slope_hindB_idx_s20=intersect(s20_cleaned_idxs.clust_s20_CL4_cleaned.clust_s20_CL4_5_cleaned,PerBrainRegions.s20.Hindbrain.idx);
slope_cereb_idx_s20=intersect(s20_cleaned_idxs.clust_s20_CL4_cleaned.clust_s20_CL4_5_cleaned,PerBrainRegions.s20.Cerebellum.idx);

figure;
scatter3(ROI_temp2.s20(s20_cleaned_idxs.idx_rsq_test_s20short_cleaned,1),ROI_temp2.s20(s20_cleaned_idxs.idx_rsq_test_s20short_cleaned,2),ROI_temp2.s20(s20_cleaned_idxs.idx_rsq_test_s20short_cleaned,3),'.');
hold on;
scatter3(ROI_temp2.s20(slope_cereb_idx_s20,1),ROI_temp2.s20(slope_cereb_idx_s20,2),ROI_temp2.s20(slope_cereb_idx_s20,3),'.');
hold on;
scatter3(ROI_temp2.s20(slope_hindB_idx_s20,1),ROI_temp2.s20(slope_hindB_idx_s20,2),ROI_temp2.s20(slope_hindB_idx_s20,3),'.');


%%
slope_hindB_idx_s60=intersect(s60_cleaned_idxs.clust_s60_CL4_cleaned.clust_s60_CL4_5_cleaned,PerBrainRegions.s60.Hindbrain.idx);
slope_cereb_idx_s60=intersect(s60_cleaned_idxs.clust_s60_CL4_cleaned.clust_s60_CL4_5_cleaned,PerBrainRegions.s60.Cerebellum.idx);

figure;
scatter3(ROI_temp2.s60(s60_cleaned_idxs.idx_rsq_test_s60short_cleaned,1),ROI_temp2.s60(s60_cleaned_idxs.idx_rsq_test_s60short_cleaned,2),ROI_temp2.s60(s60_cleaned_idxs.idx_rsq_test_s60short_cleaned,3),'.');
hold on;
scatter3(ROI_temp2.s60(slope_cereb_idx_s60,1),ROI_temp2.s60(slope_cereb_idx_s60,2),ROI_temp2.s60(slope_cereb_idx_s60,3),'.');
hold on;
scatter3(ROI_temp2.s60(slope_hindB_idx_s60,1),ROI_temp2.s60(slope_hindB_idx_s60,2),ROI_temp2.s60(slope_hindB_idx_s60,3),'.');
% hold on;
% scatter3(Zbrain_Masks{219,3}(:,1),Zbrain_Masks{219,3}(:,2),Zbrain_Masks{219,3}(:,3),'.');

%%
%%% putting all the ROIs together in coordinates

CerbHind_clustCoords_all=[];

CerbHind_clustCoords_all=ROI_temp2.f20(non_hindB_idx_f20,:);
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.f20(non_cereb_idx_f20,:));
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.f60(non_hindB_idx_f60,:));
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.f60(non_cereb_idx_f60,:));
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.s20(non_hindB_idx_s20,:));
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.s20(non_cereb_idx_s20,:));
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.s60(non_hindB_idx_s60,:));
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.s60(non_cereb_idx_s60,:));


CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.f20(slope_hindB_idx_f20,:));
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.f20(slope_cereb_idx_f20,:));
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.f60(slope_hindB_idx_f60,:));
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.f60(slope_cereb_idx_f60,:));
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.s20(slope_hindB_idx_s20,:));
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.s20(slope_cereb_idx_s20,:));
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.s60(slope_hindB_idx_s60,:));
CerbHind_clustCoords_all=vertcat(CerbHind_clustCoords_all,ROI_temp2.s60(slope_cereb_idx_s60,:));

figure;
scatter3(CerbHind_clustCoords_all(:,1),CerbHind_clustCoords_all(:,2),CerbHind_clustCoords_all(:,3),'.');


%%% I have decided that the selection of the ROIs should
%%% be from x=785-650; y=400-305; z=>69

CerbHind_clustCoords_good=find(CerbHind_clustCoords_all(:,2)<400 & CerbHind_clustCoords_all(:,2)>305 & CerbHind_clustCoords_all(:,1)<786 & CerbHind_clustCoords_all(:,1)>650 & CerbHind_clustCoords_all(:,3)>69);
CerbHind_clustCoords_good=CerbHind_clustCoords_all(CerbHind_clustCoords_good,:);

figure;
scatter3(ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,2),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,3),'.');
hold on;
scatter3(CerbHind_clustCoords_good(:,1),CerbHind_clustCoords_good(:,2),CerbHind_clustCoords_good(:,3),'.');


%%
%%% this part is to make a polygone

%k=boundary(fast_hindB_coor_all_good,0.5); %%% i am using a shrink factor (0-1) of 0.1.  

k=convhull(CerbHind_clustCoords_good);

%%% this is to make a Patch 3-D polygon with the boundarys
h=trisurf(k,CerbHind_clustCoords_good(:,1),CerbHind_clustCoords_good(:,2),CerbHind_clustCoords_good(:,3),'FaceColor','red','FaceAlpha',0.1);



%%% to reduce the faces of the patch while keeping the same shape. I am
%%% doing it to 20%
h3=reducepatch(h,.2,'verbose');

patch(h3,'EdgeColor','none','FaceAlpha',0.1);

scatter3(h3.vertices(:,1),h3.vertices(:,2),h3.vertices(:,3));

%%

tri_test = delaunayn([h3.vertices(:,1) h3.vertices(:,2) h3.vertices(:,3)]); % Generate delaunay triangulization
tn_test = tsearchn([h3.vertices(:,1) h3.vertices(:,2) h3.vertices(:,3)], tri_test, CerbHind_clustCoords_good); % Determine which triangle point is within
IsInside_test = ~isnan(tn_test); % Convert to logical vector

patch(h3,'EdgeColor','none','FaceAlpha',0.1);
hold on;
scatter3(CerbHind_clustCoords_good(:,1),CerbHind_clustCoords_good(:,2),CerbHind_clustCoords_good(:,3),'.');

notInsdie=[];
notInsdie=find(IsInside_test==0);

patch(h3,'EdgeColor','none','FaceAlpha',0.1);
hold on;
scatter3(CerbHind_clustCoords_good(:,1),CerbHind_clustCoords_good(:,2),CerbHind_clustCoords_good(:,3),'.');
hold on;
scatter3(CerbHind_clustCoords_good(notInsdie,1),CerbHind_clustCoords_good(notInsdie,2),CerbHind_clustCoords_good(notInsdie,3),'.');


%%%% now tryting with the zbrain mask of the rhomboencephalon

tri_test = delaunayn([h3.vertices(:,1) h3.vertices(:,2) h3.vertices(:,3)]); % Generate delaunay triangulization
tn_test = tsearchn([h3.vertices(:,1) h3.vertices(:,2) h3.vertices(:,3)], tri_test, Zbrain_Masks{259,3}); % Determine which triangle point is within
IsInside_test = ~isnan(tn_test); % Convert to logical vector

IsInside_good=find(IsInside_test);

scatter3(Zbrain_Masks{259,3}(IsInside_good,1),Zbrain_Masks{259,3}(IsInside_good,2),Zbrain_Masks{259,3}(IsInside_good,3),'.');

%%% it seems to work!!!! what if I do it with non convex polygon?
%%% it seems that the delaunayn function doesnt respect the concave
%%% boundaries so making the polygon with boundary gives me the same results than convhull


CerbhindB_mask=Zbrain_Masks{259,3}(IsInside_good,:);

figure;
scatter3(ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,1),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,2),ROI_temp2.f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned,3),'.');
hold on;
scatter3(CerbhindB_mask(:,1),CerbhindB_mask(:,2),CerbhindB_mask(:,3),'.');

save('specialMasks.mat','hindB_mask','CerbhindB_mask');

