%%%% this script is to make a graph with the ROIs of all the fish of different datasets for
%%%% the first, 2nd, 10th and 11th looms and compare them. I will do it by
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


load('final_F20_step1.mat','ZS_f20','idx_Fish_f20');
f20_cleaned_idxs=load('f20_cleaned_idxs.mat');



%%% selecting the fish
idx_rsq_f20=f20_cleaned_idxs.idx_rsq_test_f20short_cleaned;

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
Modules_f20=struct;
Modules_f20.Mod_loc=[];
Modules_f20.Mod_clust=[];
Modules_f20.Mod_mean=[];
for clust=[1 2 4 5 6]%1:length(fieldnames(clust_f20_CL7_cleaned_cell))

    
    
    
 idx_temp=intersect(idx_rsq_f20,clust_f20_CL7_cleaned_cell.(clustersF{clust}));   
    
options = statset('UseParallel',1); [idxKmeans_ROIs_f20 Cmap_ROIs_f20]=kmeans(ROI_temp2.f20(idx_temp,:),moduleN(clust),'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

Modules_f20.Mod_loc=vertcat(Modules_f20.Mod_loc,Cmap_ROIs_f20);
Modules_f20.Mod_clust=vertcat(Modules_f20.Mod_clust,ones(size(Cmap_ROIs_f20,1),1)*clust);

temp_module_means=[];
for i=1:size(Cmap_ROIs_f20,1)
    temp_idx=find(idxKmeans_ROIs_f20==i);
    temp_module_means(i,:)=mean(ZS_f20(idx_temp(temp_idx),:));
    
end
Modules_f20.Mod_mean=vertcat(Modules_f20.Mod_mean,temp_module_means);

end


figure;patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
%figure; 
%scatter3(Modules.Mod_loc(:,1),Modules.Mod_loc(:,2),Modules.Mod_loc(:,3));
gscatter(Modules_f20.Mod_loc(:,1),Modules_f20.Mod_loc(:,2),Modules_f20.Mod_clust);
view(-90,90);


%%% to check if it worked
% figure;
% for i=1:length(Modules.Mod_loc)
% plot(Modules.Mod_mean(i,:));
% pause(2)
% end

%%% to make a matrix of the corr
R0 = corrcoef((Modules_f20.Mod_mean(:,10:50))');
R1 = corrcoef((Modules_f20.Mod_mean(:,Loomf20_onset_idx(1):Loomf20_onset_idx(1)+40))'); %%% for loom 1
R2 = corrcoef((Modules_f20.Mod_mean(:,Loomf20_onset_idx(2):Loomf20_onset_idx(2)+40))'); %%% for loom 2
R10 = corrcoef((Modules_f20.Mod_mean(:,Loomf20_onset_idx(10):Loomf20_onset_idx(10)+40))'); %%% for loom 10
R11 = corrcoef((Modules_f20.Mod_mean(:,Loomf20_onset_idx(11):Loomf20_onset_idx(11)+40))'); %%% for loom 11


figure;
subplot(1,6,1);imagesc(R0);
subplot(1,6,2);imagesc(R1);
subplot(1,6,3);imagesc(R2);
subplot(1,6,4);imagesc(R10);
subplot(1,6,5);imagesc(R11);
subplot(1,6,6);plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Modules_f20.Mod_loc(:,1),Modules_f20.Mod_loc(:,2),Modules_f20.Mod_clust);view(-90,90);
title('idx_rsq_f20');
%R=horzcat(R1,R10,R11);

% R=[];
% R(:,1)=nonzeros(triu(R1,1));
% R(:,2)=nonzeros(triu(R10,1));
% R(:,3)=nonzeros(triu(R11,1));


%% to make the graph
n=length(Modules_f20.Mod_loc);

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

%R1(~R1) = nan; % replace zero weights with nan
weights = nonzeros(tril(R11,-1));


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
x = Modules_f20.Mod_loc(:,1);
y = Modules_f20.Mod_loc(:,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix


hold on;
gscatter(Modules_f20.Mod_loc(:,1),Modules_f20.Mod_loc(:,2),Modules_f20.Mod_clust);
view(-90,90)
hold off;


%%%% 3d doesnt work cause the network is in 2d
hold on;
patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
%figure; 


%% for f60

load('final_F60_step1_2.mat','ZS_f60','idx_Fish_f60','ZS_short_F60');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');


%%% selecting the fish
idx_rsq_f60=f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned;

%%% this is to get the ZS of the ROIs of each cluster per fish and per brain region of interest. 

clust_f60_CL7_cleaned=f60_cleaned_idxs.clust_f60_CL7_cleaned;

clust_f60_CL7_cleaned_cell={};
clust=fieldnames(clust_f60_CL7_cleaned);
for j=1:size(clustersF,1)
 clust_f60_CL7_cleaned_cell.(clustersF{j,1})=clust_f60_CL7_cleaned.(clust{j});   
end    

fish=unique(idx_Fish_f60);
fish(find(fish==47))=[]; %%% cause I also took out fish 47...





%%%% this script is to play around with the data to check some of the
%%%% methods: Organizing principles of whole-brain functional connectivity in zebrafish larvae 
%%%% Betzel 2018 biorxiv

%%%I am using the f60 dataset to test
%%% for CL7



%%% is to make functional modules in the zebrafish brain. 
moduleN=[13 22 2 26 10 26 1]; %%% for 100 modules
%moduleN=[6 10 2 13 5 13 1];  %%% for 50 modules
%moduleN=[4 5 2 6 3 5 1];  %%% for 26 modules
Modules_f60=struct;
Modules_f60.Mod_loc=[];
Modules_f60.Mod_clust=[];
Modules_f60.Mod_mean=[];
for clust=[1 2 4 5 6]%1:length(fieldnames(clust_f60_CL7_cleaned_cell))

    
 idx_temp=intersect(idx_rsq_f60,clust_f60_CL7_cleaned_cell.(clustersF{clust}));   
    
options = statset('UseParallel',1); [idxKmeans_ROIs_f60 Cmap_ROIs_f60]=kmeans(ROI_temp2.f60(idx_temp,:),moduleN(clust),'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

Modules_f60.Mod_loc=vertcat(Modules_f60.Mod_loc,Cmap_ROIs_f60);
Modules_f60.Mod_clust=vertcat(Modules_f60.Mod_clust,ones(size(Cmap_ROIs_f60,1),1)*clust);

temp_module_means=[];
for i=1:size(Cmap_ROIs_f60,1)
    temp_idx=find(idxKmeans_ROIs_f60==i);
    temp_module_means(i,:)=mean(ZS_f60(idx_temp(temp_idx),:));
    
end
Modules_f60.Mod_mean=vertcat(Modules_f60.Mod_mean,temp_module_means);

end


figure;
%patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
%figure; 
%scatter3(Modules.Mod_loc(:,1),Modules.Mod_loc(:,2),Modules.Mod_loc(:,3));
gscatter(Modules_f60.Mod_loc(:,1),Modules_f60.Mod_loc(:,2),Modules_f60.Mod_clust);
view(-90,90);


%%% to check if it worked
% figure;
% for i=1:length(Modules.Mod_loc)
% plot(Modules.Mod_mean(i,:));
% pause(2)
% end

%%% to make a matrix of the corr
R0 = corrcoef((Modules_f60.Mod_mean(:,10:50))');
R1 = corrcoef((Modules_f60.Mod_mean(:,ZS_short_F60(Loomf20_onset_idx(1):Loomf20_onset_idx(1)+40)))'); %%% for loom 1
R2 = corrcoef((Modules_f60.Mod_mean(:,ZS_short_F60(Loomf20_onset_idx(2):Loomf20_onset_idx(2)+40)))'); %%% for loom 2
R10 = corrcoef((Modules_f60.Mod_mean(:,ZS_short_F60(Loomf20_onset_idx(10):Loomf20_onset_idx(10)+40)))');%%% for loom 10
R11 = corrcoef((Modules_f60.Mod_mean(:,ZS_short_F60(Loomf20_onset_idx(11):Loomf20_onset_idx(11)+40)))');%%% for loom 11


figure;
subplot(1,6,1);imagesc(R0);
subplot(1,6,2);imagesc(R1);
subplot(1,6,3);imagesc(R2);
subplot(1,6,4);imagesc(R10);
subplot(1,6,5);imagesc(R11);
subplot(1,6,6);plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Modules_f60.Mod_loc(:,1),Modules_f60.Mod_loc(:,2),Modules_f60.Mod_clust);view(-90,90);
title('idx_rsq_f60');
%R=horzcat(R1,R10,R11);

% R=[];
% R(:,1)=nonzeros(triu(R1,1));
% R(:,2)=nonzeros(triu(R10,1));
% R(:,3)=nonzeros(triu(R11,1));


%% to make the graph
n=length(Modules_f60.Mod_loc);

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

%R1(~R1) = nan; % replace zero weights with nan
weights = nonzeros(tril(R11,-1));


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
x = Modules_f60.Mod_loc(:,1);
y = Modules_f60.Mod_loc(:,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix


hold on;
gscatter(Modules_f60.Mod_loc(:,1),Modules_f60.Mod_loc(:,2),Modules_f60.Mod_clust);
view(-90,90)
hold off;


%%%% 3d doesnt work cause the network is in 2d
hold on;
patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
%figure; 



%%

%%% this is to get the ZS of the ROIs of each cluster per fish and per brain region of interest. 

clust_f60_CL4_cleaned=f60_cleaned_idxs.clust_f60_CL4_cleaned;

clust_f60_CL4_cleaned_cell={};
clust=fieldnames(clust_f60_CL4_cleaned);
for j=1:size(clustersF,1)
 clust_f60_CL4_cleaned_cell.(clustersF{j,1})=clust_f60_CL4_cleaned.(clust{j});   
end    

fish=unique(idx_Fish_f60);
fish(find(fish==47))=[]; %%% cause I also took out fish 47...


%% for s20

load('final_S20_step1.mat','ZS_s20','idx_Fish_s20');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');


%%% selecting the fish
idx_rsq_s20=s20_cleaned_idxs.idx_rsq_test_s20short_cleaned;

%%%% this script is to play around with the data to check some of the
%%%% methods: Organizing principles of whole-brain functional connectivity in zebrafish larvae 
%%%% Betzel 2018 biorxiv

%%%I am using the s20 dataset to test
%%% for CL7


clust_s20_CL7_cleaned=s20_cleaned_idxs.clust_s20_CL7_cleaned;

clust_s20_CL7_cleaned_cell={};
clust=fieldnames(clust_s20_CL7_cleaned);
for j=1:size(clustersS,1)
 clust_s20_CL7_cleaned_cell.(clustersS{j,1})=clust_s20_CL7_cleaned.(clust{j});   
end   


%%% is to make functional modules in the zebrafish brain. 
moduleN=[13 22 2 26 10 26 1]; %%% for 100 modules
%moduleN=[6 10 2 13 5 13 1];  %%% for 50 modules
%moduleN=[4 5 2 6 3 5 1];  %%% for 26 modules
Modules_s20=struct;
Modules_s20.Mod_loc=[];
Modules_s20.Mod_clust=[];
Modules_s20.Mod_mean=[];
for clust=[1 2 4 5 6]%1:length(fieldnames(clust_s20_CL7_cleaned_cell))

    
 idx_temp=intersect(idx_rsq_s20,clust_s20_CL7_cleaned_cell.(clustersF{clust}));   
    
options = statset('UseParallel',1); [idxKmeans_ROIs_s20 Cmap_ROIs_s20]=kmeans(ROI_temp2.s20(idx_temp,:),moduleN(clust),'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

Modules_s20.Mod_loc=vertcat(Modules_s20.Mod_loc,Cmap_ROIs_s20);
Modules_s20.Mod_clust=vertcat(Modules_s20.Mod_clust,ones(size(Cmap_ROIs_s20,1),1)*clust);

temp_module_means=[];
for i=1:size(Cmap_ROIs_s20,1)
    temp_idx=find(idxKmeans_ROIs_s20==i);
    temp_module_means(i,:)=mean(ZS_s20(idx_temp(temp_idx),:));
    
end
Modules_s20.Mod_mean=vertcat(Modules_s20.Mod_mean,temp_module_means);

end


figure;
%patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
%figure; 
%scatter3(Modules.Mod_loc(:,1),Modules.Mod_loc(:,2),Modules.Mod_loc(:,3));
gscatter(Modules_s20.Mod_loc(:,1),Modules_s20.Mod_loc(:,2),Modules_s20.Mod_clust);
view(-90,90);


%%% to check if it worked
% figure;
% for i=1:length(Modules.Mod_loc)
% plot(Modules.Mod_mean(i,:));
% pause(2)
% end

%%% to make a matrix of the corr
R0 = corrcoef((Modules_s20.Mod_mean(:,10:50))');
R1 = corrcoef((Modules_s20.Mod_mean(:,S_trim(Loomf20_onset_idx(1):Loomf20_onset_idx(1)+40)))'); %%% for loom 1
R2 = corrcoef((Modules_s20.Mod_mean(:,S_trim(Loomf20_onset_idx(2):Loomf20_onset_idx(2)+40)))'); %%% for loom 2
R10 = corrcoef((Modules_s20.Mod_mean(:,S_trim(Loomf20_onset_idx(10):Loomf20_onset_idx(10)+40)))');%%% for loom 10
R11 = corrcoef((Modules_s20.Mod_mean(:,S_trim(Loomf20_onset_idx(11):Loomf20_onset_idx(11)+40)))');%%% for loom 11


figure;
subplot(1,6,1);imagesc(R0);
subplot(1,6,2);imagesc(R1);
subplot(1,6,3);imagesc(R2);
subplot(1,6,4);imagesc(R10);
subplot(1,6,5);imagesc(R11);
subplot(1,6,6);plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Modules_s20.Mod_loc(:,1),Modules_s20.Mod_loc(:,2),Modules_s20.Mod_clust);view(-90,90);
title('idx_rsq_s20');
%R=horzcat(R1,R10,R11);

% R=[];
% R(:,1)=nonzeros(triu(R1,1));
% R(:,2)=nonzeros(triu(R10,1));
% R(:,3)=nonzeros(triu(R11,1));


%% to make the graph
n=length(Modules_s20.Mod_loc);

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

%R1(~R1) = nan; % replace zero weights with nan
weights = nonzeros(tril(R11,-1));


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
x = Modules_s20.Mod_loc(:,1);
y = Modules_s20.Mod_loc(:,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix


hold on;
gscatter(Modules_s20.Mod_loc(:,1),Modules_s20.Mod_loc(:,2),Modules_s20.Mod_clust);
view(-90,90)
hold off;


%%%% 3d doesnt work cause the network is in 2d
hold on;
patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
%figure; 


%% for s60


%%% now for s60

load('final_S60_step1.mat','ZS_s60','idx_Fish_s60','ZS_short_S60');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');


%%

%%% selecting the fish
idx_rsq_s60=s60_cleaned_idxs.idx_rsq_test_s60short_cleaned;

%%% this is to get the ZS of the ROIs of each cluster per fish and per brain region of interest. 

clust_s60_CL7_cleaned=s60_cleaned_idxs.clust_s60_CL7_cleaned;

clust_s60_CL7_cleaned_cell={};
clust=fieldnames(clust_s60_CL7_cleaned);
for j=1:size(clustersS,1)
 clust_s60_CL7_cleaned_cell.(clustersS{j,1})=clust_s60_CL7_cleaned.(clust{j});   
end    

fish=unique(idx_Fish_s60);
fish(find(fish==47))=[]; %%% cause I also took out fish 47...





%%%% this script is to play around with the data to check some of the
%%%% methods: Organizing principles of whole-brain functional connectivity in zebrafish larvae 
%%%% Betzel 2018 biorxiv

%%%I am using the s60 dataset to test
%%% for CL7



%%% is to make functional modules in the zebrafish brain. 
moduleN=[13 22 2 26 10 26 1]; %%% for 100 modules
%moduleN=[6 10 2 13 5 13 1];  %%% for 50 modules
%moduleN=[4 5 2 6 3 5 1];  %%% for 26 modules
Modules_s60=struct;
Modules_s60.Mod_loc=[];
Modules_s60.Mod_clust=[];
Modules_s60.Mod_mean=[];
for clust=[1 2 4 5 6]%1:length(fieldnames(clust_s60_CL7_cleaned_cell))

    
 idx_temp=intersect(idx_rsq_s60,clust_s60_CL7_cleaned_cell.(clustersF{clust}));   
    
options = statset('UseParallel',1); [idxKmeans_ROIs_s60 Cmap_ROIs_s60]=kmeans(ROI_temp2.s60(idx_temp,:),moduleN(clust),'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

Modules_s60.Mod_loc=vertcat(Modules_s60.Mod_loc,Cmap_ROIs_s60);
Modules_s60.Mod_clust=vertcat(Modules_s60.Mod_clust,ones(size(Cmap_ROIs_s60,1),1)*clust);

temp_module_means=[];
for i=1:size(Cmap_ROIs_s60,1)
    temp_idx=find(idxKmeans_ROIs_s60==i);
    temp_module_means(i,:)=mean(ZS_s60(idx_temp(temp_idx),:));
    
end
Modules_s60.Mod_mean=vertcat(Modules_s60.Mod_mean,temp_module_means);

end


figure;
%patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
%figure; 
%scatter3(Modules.Mod_loc(:,1),Modules.Mod_loc(:,2),Modules.Mod_loc(:,3));
gscatter(Modules_s60.Mod_loc(:,1),Modules_s60.Mod_loc(:,2),Modules_s60.Mod_clust);
view(-90,90);


%%% to check if it worked
% figure;
% for i=1:length(Modules.Mod_loc)
% plot(Modules.Mod_mean(i,:));
% pause(2)
% end

%%% to make a matrix of the corr
R0 = corrcoef((Modules_s60.Mod_mean(:,10:50))');
R1 = corrcoef((Modules_s60.Mod_mean(:,ZS_short_S60(S_trim(Loomf20_onset_idx(1):Loomf20_onset_idx(1)+40))))'); %%% for loom 1
R2 = corrcoef((Modules_s60.Mod_mean(:,ZS_short_S60(S_trim(Loomf20_onset_idx(2):Loomf20_onset_idx(2)+40))))'); %%% for loom 2
R10 = corrcoef((Modules_s60.Mod_mean(:,ZS_short_S60(S_trim(Loomf20_onset_idx(10):Loomf20_onset_idx(10)+40))))');%%% for loom 10
R11 = corrcoef((Modules_s60.Mod_mean(:,ZS_short_S60(S_trim(Loomf20_onset_idx(11):Loomf20_onset_idx(11)+40))))');%%% for loom 11


figure;
subplot(1,6,1);imagesc(R0);
subplot(1,6,2);imagesc(R1);
subplot(1,6,3);imagesc(R2);
subplot(1,6,4);imagesc(R10);
subplot(1,6,5);imagesc(R11);
subplot(1,6,6);plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Modules_s60.Mod_loc(:,1),Modules_s60.Mod_loc(:,2),Modules_s60.Mod_clust);view(-90,90);
title('idx_rsq_s60');
%R=horzcat(R1,R10,R11);

% R=[];
% R(:,1)=nonzeros(triu(R1,1));
% R(:,2)=nonzeros(triu(R10,1));
% R(:,3)=nonzeros(triu(R11,1));


%% to make the graph
n=length(Modules_s60.Mod_loc);

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

%R1(~R1) = nan; % replace zero weights with nan
weights = nonzeros(tril(R11,-1));


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
x = Modules_s60.Mod_loc(:,1);
y = Modules_s60.Mod_loc(:,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix


hold on;
gscatter(Modules_s60.Mod_loc(:,1),Modules_s60.Mod_loc(:,2),Modules_s60.Mod_clust);
view(-90,90)
hold off;


%%%% 3d doesnt work cause the network is in 2d
hold on;
patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
%figure; 




