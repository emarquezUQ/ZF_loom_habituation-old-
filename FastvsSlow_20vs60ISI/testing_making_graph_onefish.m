%%%% this script is to make a graph with the ROIs of one single fish for
%%%% the first, 10th and 11th looms and compare them. I will by
%%%% separating the clusters.


load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx', 'S_trim');

datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);

%load('All_More_BrainReg.mat','PerBrainRegions');
 load('All_More_BrainReg2.mat')
load('zbrain3D.mat')
 
%RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain'};

RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain','Cerebellum','Tegmentum','Habenula',};


%% making the shape of the brain
load('Zbrain_Masks.mat');

%%% this is to make a matrix with the location of the main structures 
Zbrain_brainMask=vertcat(Zbrain_Masks{[76 113 259 274 294],3}); %%% 78 is the eyes
 
Zbrain_brainMask=unique(Zbrain_brainMask,'rows');

figure;scatter(Zbrain_brainMask(:,1),Zbrain_brainMask(:,2),'.');

Zbrain_brainMask2D_bound=boundary(Zbrain_brainMask(:,1),Zbrain_brainMask(:,2),1);

figure;plot(Zbrain_brainMask(Zbrain_brainMask2D_bound,1),Zbrain_brainMask(Zbrain_brainMask2D_bound,2));

Zbrain_brainMask2D=[];
Zbrain_brainMask2D(:,1)=Zbrain_brainMask(Zbrain_brainMask2D_bound,1);
Zbrain_brainMask2D(:,2)=Zbrain_brainMask(Zbrain_brainMask2D_bound,2);

figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');

%% for f20
%%%for f20. fish1

load('final_F20_step1.mat','ZS_f20','idx_Fish_f20');
f20_cleaned_idxs=load('f20_cleaned_idxs.mat');


unique(idx_Fish_f20)

%%% selecting the fish
idx_rsq_fish40=intersect(find(idx_Fish_f20==40),f20_cleaned_idxs.idx_rsq_test_f20short_cleaned);

%%%% this script is to play around with the data to check some of the
%%%% methods: Organizing principles of whole-brain functional connectivity in zebrafish larvae 
%%%% Betzel 2018 biorxiv

%%%I am using the f20 dataset to test
%%% for CL7


clust_f20_CL7_cleaned=f20_cleaned_idxs.clust_f20_CL7_cleaned;

clust_f20_CL7_cleaned_cell={};
clust=fieldnames(clust_f20_CL7_cleaned);
for j=1:size(clustersF,1)
 clust_f20_CL7_cleaned_cell.(clustersF{j,1})=clust_f20_CL7_cleaned.(clust{j});   
end   


%%% is to make functional modules in the zebrafish brain. 
moduleN=[13 22 2 26 10 26 1]; %%% for 100 modules
%moduleN=[6 10 2 13 5 13 1];  %%% for 50 modules
%moduleN=[4 5 2 6 3 5 1];  %%% for 26 modules
Modules=struct;
Modules.Mod_loc=[];
Modules.Mod_clust=[];
Modules.Mod_mean=[];
for clust=[1 2 4 5 6 7]%1:length(fieldnames(clust_f20_CL7_cleaned_cell))

    
    
    
 idx_temp=intersect(idx_rsq_fish40,clust_f20_CL7_cleaned_cell.(clustersF{clust}));   
    
options = statset('UseParallel',1); [idxKmeans_ROIs_f20 Cmap_ROIs_f20]=kmeans(ROI_temp2.f20(idx_temp,:),moduleN(clust),'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

Modules.Mod_loc=vertcat(Modules.Mod_loc,Cmap_ROIs_f20);
Modules.Mod_clust=vertcat(Modules.Mod_clust,ones(size(Cmap_ROIs_f20,1),1)*clust);

temp_module_means=[];
for i=1:size(Cmap_ROIs_f20,1)
    temp_idx=find(idxKmeans_ROIs_f20==i);
    temp_module_means(i,:)=mean(ZS_f20(idx_temp(temp_idx),:));
    
end
Modules.Mod_mean=vertcat(Modules.Mod_mean,temp_module_means);

end


figure;
%patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
%figure; 
%scatter3(Modules.Mod_loc(:,1),Modules.Mod_loc(:,2),Modules.Mod_loc(:,3));
gscatter(Modules.Mod_loc(:,1),Modules.Mod_loc(:,2),Modules.Mod_clust);
view(-90,90);


%%% to check if it worked
% figure;
% for i=1:length(Modules.Mod_loc)
% plot(Modules.Mod_mean(i,:));
% pause(2)
% end

%%% to make a matrix of the corr
R0 = corrcoef((Modules.Mod_mean(:,10:50))'); %%% for pre loom
R1 = corrcoef((Modules.Mod_mean(:,Loomf20_onset_idx(1):Loomf20_onset_idx(1)+40))'); %%% for loom 1
R2 = corrcoef((Modules.Mod_mean(:,Loomf20_onset_idx(2):Loomf20_onset_idx(2)+40))'); %%% for loom 2
R10 = corrcoef((Modules.Mod_mean(:,Loomf20_onset_idx(10):Loomf20_onset_idx(10)+40))'); %%% for loom 10
R11 = corrcoef((Modules.Mod_mean(:,Loomf20_onset_idx(11):Loomf20_onset_idx(11)+40))'); %%% for loom 11

figure;
subplot(1,6,1);imagesc(R0);
subplot(1,6,2);imagesc(R1);
subplot(1,6,3);imagesc(R2);
subplot(1,6,4);imagesc(R10);
subplot(1,6,5);imagesc(R11);
subplot(1,6,6);plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Modules.Mod_loc(:,1),Modules.Mod_loc(:,2),Modules.Mod_clust);view(-90,90);
title('idx_rsq_fish40');


%%%%% to make a Fisher z-transform
F0 = atanh(R0);
F1 = atanh(R1);
F2 = atanh(R2);
F10 = atanh(R10);
F11 = atanh(R11);

figure;
subplot(1,6,1);imagesc(F0);
subplot(1,6,2);imagesc(F1);
subplot(1,6,3);imagesc(F2);
subplot(1,6,4);imagesc(F10);
subplot(1,6,5);imagesc(F11);
subplot(1,6,6);plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Modules.Mod_loc(:,1),Modules.Mod_loc(:,2),Modules.Mod_clust);view(-90,90);
title('idx_rsq_fish40');
%R=horzcat(R1,R10,R11);

R=[];
R(:,1)=nonzeros(triu(R1,1));
R(:,2)=nonzeros(triu(R2,1));
R(:,3)=nonzeros(triu(R10,1));
R(:,4)=nonzeros(triu(R11,1));

figure;histogram(R(:,1));
std_R(1,1)=std(R(:,1));
Qs_R(:,1)= quantile(R(:,1),[0.025 0.25 0.50 0.75 0.975]);

%% to make the graph
n=length(Modules.Mod_loc);

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

%R1(~R1) = nan; % replace zero weights with nan
%weights = nonzeros(tril(R11,-1));
[~,~,weights] = find(tril(R1,-1));

% [row,col] = find(triu(R1,1));
% weights = R(find(triu(R1,1)));


% create the graph object:
%G = graph(row,col,weights,n);
G = graph(s,t,weights,n);

% mark the lines to remove from the graph:
threshold = 0.75; %  minimum correlation to plot
line_to_remove = isnan(weights) | abs(weights)<threshold;
% remove the lines from the graph:
G = G.rmedge(find(line_to_remove)); %#ok<FNDSB>

% plot it:
figure;
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'k';
% color positive in blue and negative in red:
p.EdgeColor = [G.Edges.Weight<0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight>0.'];
% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*1;
axis off


% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = Modules.Mod_loc(:,1);
y = Modules.Mod_loc(:,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix


hold on;
gscatter(Modules.Mod_loc(:,1),Modules.Mod_loc(:,2),Modules.Mod_clust);
view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold off;


%%%% 3d doesnt work cause the network is in 2d
% hold on;
% patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
% hold on;
%figure; 

%% playing with the brain connectivity toolbox


 W_thr10 = threshold_absolute(abs(R10), 0.75);

 W_thr11 = threshold_absolute(abs(R11), 0.75);

[Ppos10 Pneg10] = participation_coef_sign(W_thr10,Modules.Mod_clust);

[Ppos11 Pneg11] = participation_coef_sign(W_thr11,Modules.Mod_clust);


mean(Ppos10)

ratioPpos10_11=Ppos10./Ppos11;
subsPpos10_11=Ppos10-Ppos11;

figure;histogram(subsPpos10_11);

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Modules.Mod_loc(:,1),Modules.Mod_loc(:,2),25,subsPpos10_11,'filled');
view(-90,90);

%%% to check which cluster has the most nodes that increase participation
%%% at recovery. it seems is the slopehab. 
figure;histogram(Modules.Mod_clust(find(subsPpos10_11<-0.5)));

figure;histogram(Modules.Mod_clust(find(subsPpos10_11>0.5)));


%% testing doing graph with all the ROIs of one fish. 

idx_rsq_fish40


R_f40_0 = corrcoef((ZS_f20(idx_rsq_fish40,10:50))'); %%% for pre loom
R_f40_1 = corrcoef((ZS_f20(idx_rsq_fish40,Loomf20_onset_idx(1):Loomf20_onset_idx(1)+40))'); %%% for loom 1
R_f40_2 = corrcoef((ZS_f20(idx_rsq_fish40,Loomf20_onset_idx(2):Loomf20_onset_idx(2)+40))'); %%% for loom 2
R_f40_10 = corrcoef((ZS_f20(idx_rsq_fish40,Loomf20_onset_idx(10):Loomf20_onset_idx(10)+40))'); %%% for loom 10
R_f40_11 = corrcoef((ZS_f20(idx_rsq_fish40,Loomf20_onset_idx(11):Loomf20_onset_idx(11)+40))'); %%% for loom 11


figure;
subplot(1,6,1);imagesc(R_f40_0);
subplot(1,6,2);imagesc(R_f40_1);
subplot(1,6,3);imagesc(R_f40_2);
subplot(1,6,4);imagesc(R_f40_10);
subplot(1,6,5);imagesc(R_f40_11);
%subplot(1,6,6);plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Modules.Mod_loc(:,1),Modules.Mod_loc(:,2),Modules.Mod_clust);view(-90,90);
title('idx_rsq_fish40');

W_thr_f40_10 = threshold_absolute(abs(R_f40_10), 0.75);
W_thr_f40_11 = threshold_absolute(abs(R_f40_11), 0.75);

[X10,Y10,Z10] = adjacency_plot_und(W_thr_f40_10,ROI_temp2.f20(idx_rsq_fish40,:));

[X11,Y11,Z11] = adjacency_plot_und(W_thr_f40_11,ROI_temp2.f20(idx_rsq_fish40,:));

figure;
subplot(1,2,1);patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);hold on; plot3(X10,Y10,Z10);view(-90,90);
subplot(1,2,2);patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);hold on; plot3(X11,Y11,Z11);view(-90,90);

load('f20_r2050_CL5.mat','idxKmeans_ZS','idxKmeans_final','goodmaps');

%%% to quickly check i got the right kmeans. 
figure;plot(mean(ZS_f20(find(idxKmeans_ZS==goodmaps(3)),:)));



%%% now to look at the participation
[Ppos_f40_10 Pneg_f40_10] = participation_coef_sign(W_thr_f40_10,idxKmeans_ZS(idx_rsq_fish40,:));

[Ppos_f40_11 Pneg_f40_11] = participation_coef_sign(W_thr_f40_11,idxKmeans_ZS(idx_rsq_fish40,:));

%%%% an the participation change
ratioPpos_f40_10_11=Ppos_f40_10./Ppos_f40_11;
subsPpos_f40_10_11=Ppos_f40_10-Ppos_f40_11;

figure;histogram(subsPpos_f40_10_11);

figure;
%plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
scatter3(ROI_temp2.f20(idx_rsq_fish40,1),ROI_temp2.f20(idx_rsq_fish40,2),ROI_temp2.f20(idx_rsq_fish40,3),15,Ppos_f40_11,'filled');
view(-90,90); colorbar; colormap('jet');

%%% to check which cluster has the most nodes that increase participation
%%% at recovery. it seems is the slopehab. 
figure;histogram(idxKmeans_ZS(idx_rsq_fish40(find(subsPpos_f40_10_11<-0.3))));

figure;histogram(idxKmeans_ZS(idx_rsq_fish40(find(subsPpos_f40_10_11>0.3))));


