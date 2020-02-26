%%%% this script is to make a graph with the ROIs of one single fish for
%%%% the first, 10th and 11th looms and compare them. I will first test it
%%%% clustering the ROIs locations without separating clusters and then by
%%%% separating the clusters. 


%%% for CL7



load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx', 'S_trim');

datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20');
clustersF=fieldnames(mean_CL4_f20);
clustersS=fieldnames(mean_CL4_s20);

%load('All_More_BrainReg.mat','PerBrainRegions');
 load('All_More_BrainReg2.mat')
load('zbrain3D.mat')
 
%RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain'};

RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain','Cerebellum','Tegmentum','Habenula',};

%% for f20
%%%for f20. fish1

load('final_F20_step1.mat','ZS_f20','idx_Fish_f20');
f20_cleaned_idxs=load('f20_cleaned_idxs.mat');

idx_rsq_fish1=intersect(find(idx_Fish_f20==1),f20_cleaned_idxs.idx_rsq_test_f20short_cleaned);

%%%% this script is to play around with the data to check some of the
%%%% methods: Organizing principles of whole-brain functional connectivity in zebrafish larvae 
%%%% Betzel 2018 biorxiv

%%%I am using the f20 dataset to test

%%% is to make functional modules in the zebrafish brain. 

options = statset('UseParallel',1); [idxKmeans_ROIs_f20 Cmap_ROIs_f20]=kmeans(ROI_temp2.f20(idx_rsq_fish1,:),100,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

figure;patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
%figure; 
scatter3(Cmap_ROIs_f20(:,1),Cmap_ROIs_f20(:,2),Cmap_ROIs_f20(:,3));
view(-90,90);

module_means=[];
for i=1:length(Cmap_ROIs_f20)
    temp_idx=find(idxKmeans_ROIs_f20==i);
    module_means(i,:)=mean(ZS_f20(idx_rsq_fish1(temp_idx),:));
    
end

%%% to check if it worked
figure;
for i=1:length(Cmap_ROIs_f20)
plot(module_means(i,:));
pause(2)
end

%%% to make a matrix of the corr
R = corrcoef(module_means');
%R = corrcoef((module_means(:,Loomf20_onset_idx(1):Loomf20_onset_idx(1)+40))'); %%% for loom 1

%R2=atanh(R); %%% this is to do a fisher transform it doesnt work very well with it... 
%%%% i think is actually to normalize different correlation results from
%%%% different samples. so in this case it doesnt work.

%%% here i am trying to cluster the correlation matrix
options = statset('UseParallel',1); [idxKmeans_modules_corr Cmap_modules_corr]=kmeans(R,7,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');


module_means2=[];
for i=1:max(idxKmeans_modules_corr)
    temp_idx=find(idxKmeans_modules_corr==i);
    module_means2(i,:)=mean(module_means(temp_idx,:));
    
end

%%% to check if it worked
figure;
for i=1:max(idxKmeans_modules_corr)
    subplot(4,2,i);
    plot(module_means2(i,:));

end

scatter(Cmap_ROIs_f20(:,1),Cmap_ROIs_f20(:,2));


%%% this is to plot them per groups

figure;patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
view(-90,90);
%figure;
gscatter(Cmap_ROIs_f20(:,1),Cmap_ROIs_f20(:,2),idxKmeans_modules_corr);

%%%% i can also do it with a loop 
figure;
for i=1:max(idxKmeans_modules_corr)
temp_idx=find(idxKmeans_modules_corr==i);
scatter(Cmap_ROIs_f20(temp_idx,1),Cmap_ROIs_f20(temp_idx,2),'filled');
hold on;
end


%%%% now just trying to see what happens if I kluster the modules 
options = statset('UseParallel',1); [idxKmeans_modules_mean Cmap_modules_mean]=kmeans(module_means,7,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');



%%% to check if it worked
figure;
for i=1:max(idxKmeans_modules_mean)
subplot(4,2,i);
    plot(Cmap_modules_mean(i,:));
end

figure; 
for i=1:max(idxKmeans_modules_mean)
temp_idx=find(idxKmeans_modules_mean==i);
scatter(Cmap_ROIs_f20(temp_idx,1),Cmap_ROIs_f20(temp_idx,2),10,'filled');
hold on;
end


figure;
patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
view(-90,90);
%figure; 

for i=1:max(idxKmeans_modules_mean)
temp_idx=find(idxKmeans_modules_mean==i);
scatter3(Cmap_ROIs_f20(temp_idx,1),Cmap_ROIs_f20(temp_idx,2),Cmap_ROIs_f20(temp_idx,3),10,'filled');
hold on;
end

%%% this is not working
gscatter(Cmap_ROIs_f20(:,1),Cmap_ROIs_f20(:,2),idxKmeans_modules_mean);

%%%% i get similar results... not sure if i can say one is better than the
%%%% other. 

%%% I will try to classify them by the original clusters. 
rawregressF20=load('rawregressF20.mat','rawregress');
rawregressF20 = rawregressF20.('rawregress');
rawregressF20_CN=load('f20_CN_r2050_CL5_extra.mat','rawregress_CN');
rawregressF20_CN = rawregressF20_CN.('rawregress_CN');
rawregressF20(7,:)=rawregressF20_CN(1,:);

%%% to clasify them with a correlation

Correlation_group_f20={};
counter=1;
for i=1:size(rawregressF20,1)
    Correlation_temp=[];
    for idx=1:size(module_means,1)
        temp_corr=corrcoef(rawregressF20(i,:),module_means(idx,:));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_f20{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat_f20=[];
for n=1:size(rawregressF20,1)
Correlation_group_mat_f20(n,:)=Correlation_group_f20{n}(n,:);

end

%%% now if I filter looking for the max correlation
%%% 

High_corr_Nb_f20=zeros(length(Correlation_group_mat_f20),1);
for i=1:length(Correlation_group_mat_f20)
    [~,I]=max(Correlation_group_mat_f20(:,i));
    High_corr_Nb_f20(i,1)=I;
    
end


figure;
for i=1:max(High_corr_Nb_f20)
subplot(4,2,i);
idx=find(High_corr_Nb_f20==i);

    plot(mean(module_means(idx,:)));
end


figure;
patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
view(-90,90);
gscatter(Cmap_ROIs_f20(:,1),Cmap_ROIs_f20(:,2),High_corr_Nb_f20);

%% to make the graph
n=100;

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);
% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';
R(~R) = nan; % replace zero weights with nan
weights = nonzeros(tril(R,-1));
% create the graph object:
G = graph(s,t,weights,n);
% mark the lines to remove from the graph:
threshold = 0.95; %  minimum correlation to plot
line_to_remove = isnan(weights) | abs(weights)<threshold;
% remove the lines from the graph:
G = G.rmedge(find(line_to_remove)); %#ok<FNDSB>
% plot it:
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'k';
% color positive in blue and negative in red:
p.EdgeColor = [G.Edges.Weight<0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight>0.'];
% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*1;
axis off


% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = Cmap_ROIs_f20(:,1);
y = Cmap_ROIs_f20(:,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix


hold on;
for i=1:max(idxKmeans_modules_mean)
temp_idx=find(idxKmeans_modules_mean==i);
scatter(Cmap_ROIs_f20(temp_idx,1),Cmap_ROIs_f20(temp_idx,2),'filled');
hold on;
end


%%%% 3d doesnt work cause the network is in 2d
hold on;
patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
%figure; 

for i=1:max(idxKmeans_modules_mean)
temp_idx=find(idxKmeans_modules_mean==i);
scatter3(Cmap_ROIs_f20(temp_idx,1),Cmap_ROIs_f20(temp_idx,2),Cmap_ROIs_f20(temp_idx,3),5,'filled');
hold on;
end


