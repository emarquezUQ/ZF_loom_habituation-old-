

%%% this script is to analyze the fmr1 loom habituation data with graph theory.
%%% first I will make nodes of the clusters and then make means per genotype.
%%% they I will get the correlation of the responses per loom and make
%%% averages of the fish per genotype. 

%%%% NOTE: the nodes are made base on the
%%% cluster and brain region they belong!!

%%%from the wildtype loomhab dataset
load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx','S_trim');

load('Nodes_N_means_alldatasets.mat','Zbrain_brainMask2D');

load('s20_fmr1_loomhab_CN.mat','ZS_CN');

load('s20_good_NumbersNidx_Plane.mat','idx_Plane','Numbers');

%%% to get the clusters from the Kmeans done to ALL the ROIs
load('s20_postKmeans_CN.mat','idxKmeans_ZS_CN','Cmap_ZS_CN');

load('s20_good_idx_Fish.mat','idx_Fish');
idx_Fish_cat=categorical(idx_Fish);

load('s20_fmr1_loomhab_CN_part2.mat','GoodBetas_ZS_CN_selected','idx_rsq');

load('s20_fmr1_loomhab_CN_part2_High_corr_Nb.mat','High_corr_Nb');

load('fmr1loomhab_BrainRegNclean.mat','PerBrainRegions','RegionList','ROI_temp2','idx_rsq_cleaned');

%%% cleaning the clasification of the clusters. 
idx_clean=ismember(idx_rsq,idx_rsq_cleaned);
idx_clean=find(idx_clean);

High_corr_Nb=High_corr_Nb(idx_clean);

%%%% i also need to generate the list of fish saved in the
 %%%% fmr1loomhab_lists.m file
 
 load('s20_fmr1_loomhab_CN_part3.mat','idx_temp1','idx_temp2','idx_temp3','idx_temp4','idx_temp5');
 
 RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};

 
 
 %%
 
 %%%% to adapt the indexing of the looms.
 
 s_trim2Fmr1=S_trim(1:896);

 
 figure;
 plot(mean(ZS_CN(idx_rsq_cleaned(find(High_corr_Nb==3)),s_trim2Fmr1)));
 hold on;
 scatter(s_trim2Fmr1(Loomf20_onset_idx(1:20)),(zeros(1,20)),'x');
 hold on;
 scatter(S_trim(Loomf20_onset_idx(1:20)),(ones(1,20)),'*'); %%% same thing
 
 
 %% to get the nodes of the clusters
 
 
  %%% I need to order them as in my draft graphs 
 goodorder_clust=[8 7 6 5 4 9 10];
 
 %%%% merging cluster 3 and 7 (the two broad fasthab). 
 
 High_corr_Nb_new=High_corr_Nb;
 High_corr_Nb_new(find(High_corr_Nb==3))=7;
 unique(High_corr_Nb_new);
 
 figure;
 for i=1:7
    plot(mean(ZS_CN(idx_rsq_cleaned(find(High_corr_Nb_new==goodorder_clust(i))),:))); 
    hold on;
    scatter(S_trim(Loomf20_onset_idx(1:20)),(zeros(1,20)),'*');
    pause(2);
    hold off;
 end
 
%moduleN=[8 22 24 20 10 4 8 4]; %%% for 100 modules. in the right order 

Nodes2=struct;
Nodes2.Mod_loc=[];
Nodes2.Mod_clust=[];
Nodes2.Mod_KmeansID={};
Nodes2.Mod_brain=[];
counter=1;
for clust=goodorder_clust
       
 idx_temp=idx_rsq_cleaned(find(High_corr_Nb_new==clust));   
    
 for brain=1:length(RegionList)
     
     
     idx_brain_temp=PerBrainRegions.(RegionList{brain}).idx;
     
     brain_clust_idx=intersect(idx_temp,idx_brain_temp);
 
     if length(brain_clust_idx)<200
       continue 
     elseif length(brain_clust_idx)<500 & length(brain_clust_idx)>200
         moduleN=1;
     elseif length(brain_clust_idx)<1000 & length(brain_clust_idx)>500
     moduleN=2;
     elseif length(brain_clust_idx)<3000 & length(brain_clust_idx)>999
     moduleN=3;
      elseif length(brain_clust_idx)>3000 
     moduleN=4;
     end   
     
options = statset('UseParallel',1); [idxKmeans_ROIs Cmap_ROIs]=kmeans(ROI_temp2(brain_clust_idx,:),moduleN,'Options',options,'Distance','sqeuclidean','Replicates',5,'MaxIter',1000,'Display','final');

Nodes2.Mod_loc=vertcat(Nodes2.Mod_loc,Cmap_ROIs);
Nodes2.Mod_clust=vertcat(Nodes2.Mod_clust,ones(size(Cmap_ROIs,1),1)*clust);
Nodes2.Mod_brain=vertcat(Nodes2.Mod_brain,ones(size(Cmap_ROIs,1),1)*brain);


%nodes=find(Nodes2.Mod_clust==clust);
KID=unique(idxKmeans_ROIs);
for ID=1:length(KID)
 Nodes2.Mod_KmeansID{counter,1}=brain_clust_idx(find(idxKmeans_ROIs==ID));   
counter=counter+1;
end
 end
end 

%%% to check if it would work
clust=3;
idx_temp=idx_rsq_cleaned(find(High_corr_Nb_new==clust)); 
figure;scatter(ROI_temp2(idx_temp,1),ROI_temp2(idx_temp,2));

%%%% i will need to adjust the colors but it seems that it is working. 
figure;
%patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
%figure; 
%scatter3(Nodes.Mod_loc(:,1),ModulNodeses.Mod_loc(:,2),Nodes.Mod_loc(:,3));
gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_brain);
view(-90,90);

figure;
%patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
%figure; 
%scatter3(Nodes.Mod_loc(:,1),ModulNodeses.Mod_loc(:,2),Nodes.Mod_loc(:,3));
gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_clust,'rcbgyrm','.',20,'on');
view(-90,90);

%%
%%% it seems that it worked... now i need to get ROIs of each fish per
%%% genotype.

fish=vertcat(list1,list2,list3,list4);
list5=union(list1,list3);

%%% it seems that fish 201810048 from list1 dont have any ROIs... dont know
%%% why. might need to check if it is not a mistake in the name. i am
%%% getting rid of it on the list. 
find(idx_Fish==201810048);
fish(find(fish==201810048))=[];


for f=1:length(fish)
    tempfish=find(idx_Fish==fish(f));
    
    if ismember(fish(f),list2)
        group='fmr1';
        %ff=find(ismember(list2,fish(f)));
    elseif ismember(fish(f),list4)
        group='control';
        %ff=find(ismember(list4,fish(f)));
    elseif ismember(fish(f),list5)
        group='hets';
        %ff=find(ismember(list5,fish(f)));
    else 
    end
    temp_group=struct;
  for clust=goodorder_clust
      
      temp_idx_clust=intersect(idx_rsq_cleaned(find(High_corr_Nb_new==clust)),tempfish);
     
      
    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes2.Mod_clust==clust);
    for node=(find(Nodes2.Mod_clust==clust))'
    temp_idx=intersect(temp_idx_clust,Nodes2.Mod_KmeansID{node,1});
    
%      figure;plot(mean(ZS_CN(temp_idx,:))); %% to check if it works.
%     figure;
%     plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
%     hold on;
%     scatter(ROI_temp2(temp_idx,1),ROI_temp2(temp_idx,2),'filled');
%     view(-90,90);
    
    
    
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx;  
    
    end

    Nodes2.ROIs_idx.(group).(strcat('fish_',num2str(fish(f)))).(strcat('clust_',num2str(clust)))=idx_ROIs_Node;
    
    
  end
end


 
 %%
 %%% to get the means
  for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes2.ROIs_idx.(group));
 for f=1:length(fish)
     
for clust=goodorder_clust
   

    for node=(find(Nodes2.Mod_clust==clust))'
    
   
    temp_idx=Nodes2.ROIs_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).idx;
%         figure;
%         plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
%         hold on;
%         scatter(ROI_temp2(temp_idx,1),ROI_temp2(temp_idx,2));
%         view(-90,90);
    
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_CN(temp_idx,:);
        plot(temp_mean);
    else
        temp_mean=mean(ZS_CN(temp_idx,:));
        plot(temp_mean);
    end
    
    Nodes2.ROIs_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

  end
  
  
  %%% making means matrices per fish and testing a correlation
  for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes2.ROIs_idx.(group));
for f=1:length(fish)

    temp_matrix=NaN(size(Nodes2.Mod_clust,1),904);
    for clust=goodorder_clust
        
        for node=(find(Nodes2.Mod_clust==clust))'
            
        temp_mean=Nodes2.ROIs_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes2.mean_matrix.(group).(fish{f})=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes2.corr_matrix.(group).(fish{f})=R_temp;
        
        Nodes2.NaNtest.(group).(fish{f})=isnan(R_temp);
        
end
  end

  %%% testing the matrix with all fish to see if there are still gaps
Matrix_mean=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes2.NaNtest.(group).(fish{f}));     
end

Matrix_mean=nanmean(Matrix_mean,3);
  

%%% is it better to use less distance for recluting the ROIs? cause maybe
%%% sometimes one node 'steals' from another one. I did some test and I got
%%% the best results with the low quantile, although very close to the min(d).
%%% avg of NaN for the whole matrix: min=0.4251 and quantile=0.4179 (mean of NaNs)

mean(mean(Matrix_mean))


%% getting the corrmatrix for each loom and each fish
Data_corrMat2=struct;

for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes2.mean_matrix.(group));
    for f=1:length(fish)
    
        temp_mean=Nodes2.mean_matrix.(group).(fish{f});
        
        temp_R={};
        
        for k=1:21
        
        if k==1
            temp_R{1,k}=corrcoef((temp_mean(:,10:45))');
        elseif k==11||k==21 
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28))');
        else
            
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36))'); 
                           
        end
        end
        
    Data_corrMat2.(group).(fish{f}).loomsR=temp_R;
        
    end
             
end

%% making means of each loom per dataset

for g=1:3
     group=groupnames{g,1};
    
    fish=fieldnames(Data_corrMat2.(group));
    
    Mean_corrMat={};
    for k=1:21
        
        temp_Mean_corrMat=[];
        
        for f=1:length(fish)
            
        temp_mat=Data_corrMat2.(group).(fish{f}).loomsR{1,k};
        
        temp_Mean_corrMat=cat(3,temp_Mean_corrMat,temp_mat);
        end
        
        Mean_corrMat{1,k}=nanmean(temp_Mean_corrMat,3);
    end
    Data_corrMat2.(group).Mean_corrMat=Mean_corrMat;
    
end


%%
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_clust);view(-90,90);
title('Model Nodes');

counter=1;
figure;
for g=1:3
     group=groupnames{g,1};
subplot(3,8,counter);imagesc(Data_corrMat2.(group).Mean_corrMat{1,1}); %% for pre loom
subplot(3,8,counter+1);imagesc(Data_corrMat2.(group).Mean_corrMat{1,2}); %% for 1st loom
subplot(3,8,counter+2);imagesc(Data_corrMat2.(group).Mean_corrMat{1,3});
subplot(3,8,counter+3);imagesc(Data_corrMat2.(group).Mean_corrMat{1,4});
subplot(3,8,counter+4);imagesc(Data_corrMat2.(group).Mean_corrMat{1,5});
subplot(3,8,counter+5);imagesc(Data_corrMat2.(group).Mean_corrMat{1,6});
subplot(3,8,counter+6);imagesc(Data_corrMat2.(group).Mean_corrMat{1,11});%% for 10th loom
subplot(3,8,counter+7);imagesc(Data_corrMat2.(group).Mean_corrMat{1,12}); %% for 11th loom
title(group);

counter=counter+8;
end

%% now to plot some graphs. 

for g=1:3
    figure;
    count=1;
    group=groupnames{g,1};
   sgtitle(group) 
   
    for k=[2 3 4 5 6 11 12]
    
    R=Data_corrMat2.(group).Mean_corrMat{1,k};
        
        
n=length(Nodes2.Mod_loc);

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

[~,~,weights] = find(tril(R,-1));

% create the graph object:
G = graph(s,t,weights,n);

% mark the lines to remove from the graph:
threshold = 0.75; %  minimum correlation to plot
line_to_remove = isnan(weights) | abs(weights)<threshold;
% remove the lines from the graph:
G = G.rmedge(find(line_to_remove)); %#ok<FNDSB>

% plot it:
subplot(1,7,count);
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'k';
% color positive in blue and negative in red:
p.EdgeColor = [G.Edges.Weight<0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight>0.'];
% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*1;
axis off

% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = Nodes2.Mod_loc(:,1);
y = Nodes2.Mod_loc(:,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_clust,'rcbgyrm','.',20,'off');

view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold off;

count=count+1

    end
end

%%



%% for the substraction of contros minus fmr1
counter=1;
figure;    
subplot(1,8,counter);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,1},Data_corrMat2.fmr1.Mean_corrMat{1,1})); caxis([-1 1]);colormap('jet');%% for pre loom
subplot(1,8,counter+1);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,2},Data_corrMat2.fmr1.Mean_corrMat{1,2}));caxis([-1 1]); colormap('jet');%% for 1st loom
subplot(1,8,counter+2);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,3},Data_corrMat2.fmr1.Mean_corrMat{1,3}));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+3);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,4},Data_corrMat2.fmr1.Mean_corrMat{1,4}));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+4);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,5},Data_corrMat2.fmr1.Mean_corrMat{1,5}));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+5);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,6},Data_corrMat2.fmr1.Mean_corrMat{1,6}));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+6);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,11},Data_corrMat2.fmr1.Mean_corrMat{1,11}));caxis([-1 1]);colormap('jet');%% for 10th loom
subplot(1,8,counter+7);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,12},Data_corrMat2.fmr1.Mean_corrMat{1,12})); caxis([-1 1]);colormap('jet');%% for 11th loom
sgtitle('controls minus fmr1');


%% to check the participation change from 10 to 11th 

%%% this needs to be done with the Brain Connectivity Toolbox. 
%%% so i need to have it in the path. 

W_thr10 = threshold_absolute(abs(Data_corrMat2.fmr1.Mean_corrMat{1,11}), 0.75);

 W_thr11 = threshold_absolute(abs(Data_corrMat2.fmr1.Mean_corrMat{1,12}), 0.75);

[Ppos10 Pneg10] = participation_coef_sign(W_thr10,Nodes2.Mod_clust);

[Ppos11 Pneg11] = participation_coef_sign(W_thr11,Nodes2.Mod_clust);


mean(Ppos10)

ratioPpos10_11=Ppos10./Ppos11;
subsPpos10_11=Ppos10-Ppos11;

figure;histogram(subsPpos10_11);

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),25,subsPpos10_11,'filled'); colorbar;colormap('jet')
view(-90,90);
title('fmr1');
%%% to check which cluster has the most nodes that increase participation
%%% at recovery. it seems is the slopehab. 
figure;histogram(Nodes2.Mod_clust(find(subsPpos10_11<-0.3)));

figure;histogram(Nodes2.Mod_clust(find(subsPpos10_11>0.3)));

%% checking a few things


%%% to see which ROIs have less representation
counter=1;
figure;
meanProp=[];
NaN_nodes={};
for g=1:3
    
  group=groupnames{g,1};
    
    fish=fieldnames(Nodes2.NaNtest.(group));  
     
    
    Matrix_mean=[];
    temp_NaN_nodes=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes2.NaNtest.(group).(fish{f}));  
 temp=double(diag(Nodes2.NaNtest.(group).(fish{f})));
 temp_NaN_nodes=horzcat(temp_NaN_nodes,temp);
end

%Matrix_mean=(sum(Matrix_mean,3))/length(fish); %% is the same
Matrix_mean=nanmean(Matrix_mean,3);  


subplot(1,3,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),20,(1-mean(Matrix_mean)),'filled'); colorbar;colormap('jet'); caxis([0 1]);
view(-90,90);
title(group);
hold off;

counter=counter+1;

meanProp=horzcat(meanProp,(mean(Matrix_mean)'));

NaN_nodes{g}=temp_NaN_nodes;
end

fish_perNode=[];
meanProp_good=[];
for  g=1:3
    fish_perNode(:,g)=abs(sum(NaN_nodes{g},2)-size(NaN_nodes{g},2));
    meanProp_good(:,g)=1-(sum(NaN_nodes{g},2)/size(NaN_nodes{g},2));
    
end


%%% to see which fish are contributing to specific nodes. 
fish=fieldnames(Nodes2.NaNtest.hets);
fmr1_hindbrain_fish={};
for f=1:length(fish)
    
    temp_mat=Nodes2.NaNtest.hets.(fish{f});
    if temp_mat(26,1)==0 
    fmr1_hindbrain_fish{f}=(fish{f})
    else
    end
end

%%% to see which fish are contributing to specific nodes. 
fish=fieldnames(Nodes2.NaNtest.hets);
fmr1_hindbrain_fish={};
for f=1:length(fish)
    
    temp_mat=Nodes2.NaNtest.hets.(fish{f});
    if temp_mat(26,1)==0 
    fmr1_hindbrain_fish{f}=(fish{f})
    else
    end
end



%% what if i put a threshold on proportion of fish needed to contribute to a node?
%%%% before I run the crosscorrelation. 
%%% with at least 25% of the fish i loose some nodes in the hindbrain although a visual effect still
%%% seems to remain there. 

%discard=find(min((1-meanProp)')<0.25); %%% I could use 0.25 (discard 8 nodes, but leaves some in the hindbrain), 0.33 (discard 18 nodes) , or 0/5 (discard 60 nodes).
discard=find(min(meanProp_good,[],2)<0.25); %%% this is the proper way!!! with 0.25 i discard 4 nodes, 0.33=7 nodes and 0.5=27 nodes
keep=find(ismember([1:89],discard)==0);


%%%% with 0.25, i discarded mostly from contralateral side. 
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_brain(keep));view(-90,90);
title('Model Nodes');
%%% adding numbers
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_brain(:));view(-90,90);
title('Model Nodes');
a = [1:89]'; b = num2str(a); c = cellstr(b);
dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
text(x+dx, y+dy, c);


counter=1;
figure;
for g=1:3
     group=groupnames{g,1};
subplot(3,8,counter);imagesc(Data_corrMat2.(group).Mean_corrMat{1,1}(keep,keep)); %% for pre loom
subplot(3,8,counter+1);imagesc(Data_corrMat2.(group).Mean_corrMat{1,2}(keep,keep)); %% for 1st loom
subplot(3,8,counter+2);imagesc(Data_corrMat2.(group).Mean_corrMat{1,3}(keep,keep));
subplot(3,8,counter+3);imagesc(Data_corrMat2.(group).Mean_corrMat{1,4}(keep,keep));
subplot(3,8,counter+4);imagesc(Data_corrMat2.(group).Mean_corrMat{1,5}(keep,keep));
subplot(3,8,counter+5);imagesc(Data_corrMat2.(group).Mean_corrMat{1,6}(keep,keep));
subplot(3,8,counter+6);imagesc(Data_corrMat2.(group).Mean_corrMat{1,11}(keep,keep));%% for 10th loom
subplot(3,8,counter+7);imagesc(Data_corrMat2.(group).Mean_corrMat{1,12}(keep,keep)); %% for 11th loom
title(group);

counter=counter+8;
end

%% now to plot some graphs. 

for g=1:3
    figure;
    count=1;
    group=groupnames{g,1};
   sgtitle(group) 
   
    for k=[2 3 4 5 6 11 12]
    
    R=Data_corrMat2.(group).Mean_corrMat{1,k}(keep,keep);
        
        
n=length(Nodes2.Mod_loc(keep));

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

[~,~,weights] = find(tril(R,-1));

% create the graph object:
G = graph(s,t,weights,n);

% mark the lines to remove from the graph:
threshold = 0.75; %  minimum correlation to plot
line_to_remove = isnan(weights) | abs(weights)<threshold;
% remove the lines from the graph:
G = G.rmedge(find(line_to_remove)); %#ok<FNDSB>

% plot it:
subplot(1,7,count);
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'k';
% color positive in blue and negative in red:
p.EdgeColor = [G.Edges.Weight<0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight>0.'];
% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*1;
axis off

% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = Nodes2.Mod_loc(keep,1);
y = Nodes2.Mod_loc(keep,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_clust(keep),'grbbggrm','.',20,'off');
hold on;
view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold off;

count=count+1;

    end
end


%% for the substraction of contros minus fmr1
counter=1;
figure;    
subplot(1,8,counter);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,1}(keep,keep),Data_corrMat2.fmr1.Mean_corrMat{1,1}(keep,keep))); caxis([-1 1]);colormap('jet');%% for pre loom
subplot(1,8,counter+1);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,2}(keep,keep),Data_corrMat2.fmr1.Mean_corrMat{1,2}(keep,keep)));caxis([-1 1]); colormap('jet');%% for 1st loom
subplot(1,8,counter+2);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,3}(keep,keep),Data_corrMat2.fmr1.Mean_corrMat{1,3}(keep,keep)));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+3);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,4}(keep,keep),Data_corrMat2.fmr1.Mean_corrMat{1,4}(keep,keep)));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+4);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,5}(keep,keep),Data_corrMat2.fmr1.Mean_corrMat{1,5}(keep,keep)));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+5);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,6}(keep,keep),Data_corrMat2.fmr1.Mean_corrMat{1,6}(keep,keep)));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+6);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,11}(keep,keep),Data_corrMat2.fmr1.Mean_corrMat{1,11}(keep,keep)));caxis([-1 1]);colormap('jet');%% for 10th loom
subplot(1,8,counter+7);imagesc(minus(Data_corrMat2.control.Mean_corrMat{1,12}(keep,keep),Data_corrMat2.fmr1.Mean_corrMat{1,12}(keep,keep))); caxis([-1 1]);colormap('jet');%% for 11th loom
sgtitle('controls minus fmr1');


%% to check the participation change from 10 to 11th 
W_thr10 = threshold_absolute(abs(Data_corrMat2.control.Mean_corrMat{1,11}(keep,keep)), 0.75);

 W_thr11 = threshold_absolute(abs(Data_corrMat2.control.Mean_corrMat{1,12}(keep,keep)), 0.75);

[Ppos10 Pneg10] = participation_coef_sign(W_thr10,Nodes2.Mod_clust(keep));

[Ppos11 Pneg11] = participation_coef_sign(W_thr11,Nodes2.Mod_clust(keep));


mean(Ppos10)

ratioPpos10_11=Ppos10./Ppos11;
subsPpos10_11=Ppos10-Ppos11;

figure;histogram(subsPpos10_11);

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,subsPpos10_11,'filled'); colorbar;colormap('jet')
view(-90,90);
title('control');
%%% to check which cluster has the most nodes that increase participation
%%% at recovery. it seems is the slopehab. 
figure;histogram(Nodes2.Mod_clust(find(subsPpos10_11<-0.3)));

figure;histogram(Nodes2.Mod_clust(find(subsPpos10_11>0.3)));


%%
%%%% it seems we will work with the atribution by the kmeans to the nodes. 

save('NodesNgraphFmr1Loomhab2.mat','Nodes2','goodorder_clust','Data_corrMat2');

