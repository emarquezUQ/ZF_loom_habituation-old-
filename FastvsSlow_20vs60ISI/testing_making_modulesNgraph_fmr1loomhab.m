%%% this script is to analyze the fmr1 loom habituation data with graph theory.
%%% first I will make nodes of the clusters and then make means per genotype.
%%% they I will get the correlation of the responses per loom and make
%%% averages of the fish per genotype. 

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
 goodorder_clust=[3 7 8 6 5 4 9 10];
 
 figure;
 for i=1:8
    plot(mean(ZS_CN(idx_rsq_cleaned(find(High_corr_Nb==goodorder_clust(i))),:))); 
    hold on;
    scatter(S_trim(Loomf20_onset_idx(1:20)),(zeros(1,20)),'*');
    pause(2);
    hold off;
 end
 
 %%% remember that I am skiping the first 2 clusters from the original 10
 %%% ~95000 ROIs from all the clusters so:
%moduleN=[13 22 2(sound) 26 10 26 3]; %%% for 100 modules of the wild type dataset 

 proportion=[];
 for i=1:8
     temp=(length(idx_rsq_cleaned(find(High_corr_Nb==goodorder_clust(i)))))*100/95000;
     temp=ceil(temp);
     proportion=vertcat(proportion,temp);
     
 end
 proportion'
sum(proportion) %%% is giving me 94=[8 21 24 19 10 3 7 2] so i can allocate some arbitrarly 
 
moduleN=[8 22 24 20 10 4 8 4]; %%% for 100 modules. in the right order 

Nodes=struct;
Nodes.Mod_loc=[];
Nodes.Mod_clust=[];
Nodes.Mod_KmeansID={};
for clust=goodorder_clust
       
 idx_temp=idx_rsq_cleaned(find(High_corr_Nb==clust));   
    
options = statset('UseParallel',1); [idxKmeans_ROIs Cmap_ROIs]=kmeans(ROI_temp2(idx_temp,:),moduleN(find(goodorder_clust==clust)),'Options',options,'Distance','sqeuclidean','Replicates',5,'MaxIter',1000,'Display','final');

Nodes.Mod_loc=vertcat(Nodes.Mod_loc,Cmap_ROIs);
Nodes.Mod_clust=vertcat(Nodes.Mod_clust,ones(size(Cmap_ROIs,1),1)*clust);

nodes=find(Nodes.Mod_clust==clust);
KID=unique(idxKmeans_ROIs);
for ID=1:length(KID)
Nodes.Mod_KmeansID{nodes(ID),1}=idx_temp(find(idxKmeans_ROIs==ID));
end

end 

%%% to check if it would work
clust=3;
idx_temp=idx_rsq_cleaned(find(High_corr_Nb==clust)); 
figure;scatter(ROI_temp2(idx_temp,1),ROI_temp2(idx_temp,2));

%%%% i will need to adjust the colors but it seems that it is working. 
figure;
%patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
%figure; 
%scatter3(Nodes.Mod_loc(:,1),ModulNodeses.Mod_loc(:,2),Nodes.Mod_loc(:,3));
gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Nodes.Mod_clust);
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
      
      d=pdist(Nodes.Mod_loc(find(Nodes.Mod_clust==clust),:));
    %temp_min=(min(d)/2)-1; %%% doesnt work that well
    %temp_min=min(d); %%% maybe 2nd best results. close with the quantile... it depends on the fish
    temp_min=quantile(d,[0.025]); %%% mabye the best results. in varies among fish. what i used before for the wildtypes dataset
    %temp_min=quantile(d,[0.025])/2; %%% doesnt work that well
    
    %%%% i did some calculations on the number of NaNs i get with the
    %%%% different distances. maybe the min is very slightly better...
    %%%% although there is not much difference. 
    
      temp_idx_clust=intersect(idx_rsq_cleaned(find(High_corr_Nb==clust)),tempfish);
      %figure;plot(mean(ZS_CN(temp,:))); %% to check if it works.
       
      D = pdist2(ROI_temp2(temp_idx_clust,:),Nodes.Mod_loc);  

    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    [temp_row,~]=find(D(:,node)<temp_min);
    
    %%% in this bit I am tryting to fix the overlaping. 
        if ~isempty(find(ismember(temp_row,row)))
        temp=find(ismember(temp_row,row));
        temp_row(temp)=[];
        row=vertcat(row,temp_row);
 
        else
        row=vertcat(row,temp_row);
        end
    
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx_clust(temp_row);  
    
    end

    Nodes.ROIs_idx.(group).(strcat('fish_',num2str(fish(f)))).(strcat('clust_',num2str(clust)))=idx_ROIs_Node;
    
    
  end
end

 %%% to check if i fixed the overlaping. it seemed that it worked. 
groupnames=fieldnames(Nodes.ROIs_idx);

 for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_idx.(group));
     
 for f=1:length(fish)
     
 for clust=goodorder_clust
   
    idx_ROIs_Node=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    
    temp_idx=Nodes.ROIs_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).idx;
    
    idx_ROIs_Node=vertcat(idx_ROIs_Node,temp_idx);
    end
     
    B = unique(idx_ROIs_Node); % which will give you the unique elements of A in array B
    Ncount = histc(idx_ROIs_Node, B); % this willgive the number of occurences of each unique element
    unique(Ncount)
    
    
end
end
 end


 
 
 %%
 %%% to get the means
  for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_idx.(group));
 for f=1:length(fish)
     
for clust=goodorder_clust
   

    for node=(find(Nodes.Mod_clust==clust))'
    
   
    temp_idx=Nodes.ROIs_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).idx;
%         figure;
%         plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
%         hold on;
%         scatter(ROI_temp2(temp_idx,1),ROI_temp2(temp_idx,2));
%         view(-90,90);
    
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_CN(temp_idx,:);
        %plot(temp_mean);
    else
        temp_mean=mean(ZS_CN(temp_idx,:));
        %plot(temp_mean);
    end
    
    Nodes.ROIs_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

  end
  
  
  %%% making means matrices per fish and testing a correlation
  for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_idx.(group));
for f=1:length(fish)

    temp_matrix=NaN(100,904);
    for clust=goodorder_clust
        
        for node=(find(Nodes.Mod_clust==clust))'
            
        temp_mean=Nodes.ROIs_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes.mean_matrix.(group).(fish{f})=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes.corr_matrix.(group).(fish{f})=R_temp;
        
        Nodes.NaNtest.(group).(fish{f})=isnan(R_temp);
        
end
  end

  %%% testing the matrix with all fish to see if there are still gaps
Matrix_mean=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes.NaNtest.(group).(fish{f}));     
end

Matrix_mean=nanmean(Matrix_mean,3);
  

%%% is it better to use less distance for recluting the ROIs? cause maybe
%%% sometimes one node 'steals' from another one. I did some test and I got
%%% the best results with the low quantile, although very close to the min(d).
%%% avg of NaN for the whole matrix: min=0.4251 and quantile=0.4179 (mean of NaNs)

mean(mean(Matrix_mean))


%% getting the corrmatrix for each loom and each fish
Data_corrMatQ=struct;

for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.mean_matrix.(group));
    for f=1:length(fish)
    
        temp_mean=Nodes.mean_matrix.(group).(fish{f});
        
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
        
    Data_corrMatQ.(group).(fish{f}).loomsR=temp_R;
        
    end
             
end

%% making means of each loom per dataset

for g=1:3
     group=groupnames{g,1};
    
    fish=fieldnames(Data_corrMatQ.(group));
    
    Mean_corrMat={};
    for k=1:21
        
        temp_Mean_corrMat=[];
        
        for f=1:length(fish)
            
        temp_mat=Data_corrMatQ.(group).(fish{f}).loomsR{1,k};
        
        temp_Mean_corrMat=cat(3,temp_Mean_corrMat,temp_mat);
        end
        
        Mean_corrMat{1,k}=nanmean(temp_Mean_corrMat,3);
    end
    Data_corrMatQ.(group).Mean_corrMat=Mean_corrMat;
    
end


%%
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Nodes.Mod_clust);view(-90,90);
title('Model Nodes');

counter=1;
figure;
for g=1:3
     group=groupnames{g,1};
subplot(3,8,counter);imagesc(Data_corrMatQ.(group).Mean_corrMat{1,1}); %% for pre loom
subplot(3,8,counter+1);imagesc(Data_corrMatQ.(group).Mean_corrMat{1,2}); %% for 1st loom
subplot(3,8,counter+2);imagesc(Data_corrMatQ.(group).Mean_corrMat{1,3});
subplot(3,8,counter+3);imagesc(Data_corrMatQ.(group).Mean_corrMat{1,4});
subplot(3,8,counter+4);imagesc(Data_corrMatQ.(group).Mean_corrMat{1,5});
subplot(3,8,counter+5);imagesc(Data_corrMatQ.(group).Mean_corrMat{1,6});
subplot(3,8,counter+6);imagesc(Data_corrMatQ.(group).Mean_corrMat{1,11});%% for 10th loom
subplot(3,8,counter+7);imagesc(Data_corrMatQ.(group).Mean_corrMat{1,12}); %% for 11th loom
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
    
    R=Data_corrMatQ.(group).Mean_corrMat{1,k};
        
        
n=length(Nodes.Mod_loc);

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
x = Nodes.Mod_loc(:,1);
y = Nodes.Mod_loc(:,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Nodes.Mod_clust,'grbbggrm','.',20,'off');

view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold off;

count=count+1

    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getting the idx with Kmeans of the ROIs
%%%% I will attempt this cause I still get some nodes that are not very
%%%% well represented... 

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
    
       
    idx_ROIs_Node=struct;
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    temp_idx=intersect(tempfish,Nodes.Mod_KmeansID{node,1});
    %figure;plot(mean(ZS_CN(temp_idx,:))); %% to check if it works.
    
    
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx;  
    
    end

    Nodes.ROIs_K_idx.(group).(strcat('fish_',num2str(fish(f)))).(strcat('clust_',num2str(clust)))=idx_ROIs_Node;
    
    
  end
end


%%
 %%% to get the means
  for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_K_idx.(group));
 for f=1:length(fish)
     
for clust=goodorder_clust
   

    for node=(find(Nodes.Mod_clust==clust))'
    
   
    temp_idx=Nodes.ROIs_K_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).idx;
%         figure;
%         plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
%         hold on;
%         scatter(ROI_temp2(temp_idx,1),ROI_temp2(temp_idx,2));
%         view(-90,90);
    
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_CN(temp_idx,:);
        %plot(temp_mean);
    else
        temp_mean=mean(ZS_CN(temp_idx,:));
        %plot(temp_mean);
    end
    
    Nodes.ROIs_K_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

  end
  
  
  %%% making means matrices per fish and testing a correlation
  for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_K_idx.(group));
for f=1:length(fish)

    temp_matrix=NaN(100,904);
    for clust=goodorder_clust
        
        for node=(find(Nodes.Mod_clust==clust))'
            
        temp_mean=Nodes.ROIs_K_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes.mean_matrix_K.(group).(fish{f})=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes.corr_matrix_K.(group).(fish{f})=R_temp;
        
        Nodes.NaNtest_K.(group).(fish{f})=isnan(R_temp);
        
end
  end

  %%% testing the matrix with all fish to see if there are still gaps
Matrix_mean=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes.NaNtest_K.(group).(fish{f}));     
end

Matrix_mean=nanmean(Matrix_mean,3);
  

%%% is it better to use less distance for recluting the ROIs? cause maybe
%%% sometimes one node 'steals' from another one. I did some test and I got
%%% the best results with the min(d), although very close to the quantile.
%%% avg of NaN for the whole matrix: min=0.4251 and quantile=0.4179

mean(mean(Matrix_mean))


%% getting the corrmatrix for each loom and each fish
Data_corrMatK=struct;

for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.mean_matrix_K.(group));
    for f=1:length(fish)
    
        temp_mean=Nodes.mean_matrix_K.(group).(fish{f});
        
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
        
    Data_corrMatK.(group).(fish{f}).loomsR=temp_R;
        
    end
             
end

%% making means of each loom per dataset

for g=1:3
     group=groupnames{g,1};
    
    fish=fieldnames(Data_corrMatK.(group));
    
    Mean_corrMat={};
    for k=1:21
        
        temp_Mean_corrMat=[];
        
        for f=1:length(fish)
            
        temp_mat=Data_corrMatK.(group).(fish{f}).loomsR{1,k};
        
        temp_Mean_corrMat=cat(3,temp_Mean_corrMat,temp_mat);
        end
        
        Mean_corrMat{1,k}=nanmean(temp_Mean_corrMat,3);
    end
    Data_corrMatK.(group).Mean_corrMat=Mean_corrMat;
    
end


%%
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Nodes.Mod_clust);view(-90,90);
title('Model Nodes');

counter=1;
figure;
for g=1:3
     group=groupnames{g,1};
subplot(3,8,counter);imagesc(Data_corrMatK.(group).Mean_corrMat{1,1}); %% for pre loom
subplot(3,8,counter+1);imagesc(Data_corrMatK.(group).Mean_corrMat{1,2}); %% for 1st loom
subplot(3,8,counter+2);imagesc(Data_corrMatK.(group).Mean_corrMat{1,3});
subplot(3,8,counter+3);imagesc(Data_corrMatK.(group).Mean_corrMat{1,4});
subplot(3,8,counter+4);imagesc(Data_corrMatK.(group).Mean_corrMat{1,5});
subplot(3,8,counter+5);imagesc(Data_corrMatK.(group).Mean_corrMat{1,6});
subplot(3,8,counter+6);imagesc(Data_corrMatK.(group).Mean_corrMat{1,11});%% for 10th loom
subplot(3,8,counter+7);imagesc(Data_corrMatK.(group).Mean_corrMat{1,12}); %% for 11th loom
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
    
    R=Data_corrMatK.(group).Mean_corrMat{1,k};
        
        
n=length(Nodes.Mod_loc);

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
x = Nodes.Mod_loc(:,1);
y = Nodes.Mod_loc(:,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Nodes.Mod_clust,'grbbggrm','.',20,'off');

view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold off;

count=count+1

    end
end

%% for the substraction of contros minus fmr1
counter=1;
figure;    
subplot(1,8,counter);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,1},Data_corrMatK.fmr1.Mean_corrMat{1,1})); caxis([-1 1]);colormap('jet');%% for pre loom
subplot(1,8,counter+1);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,2},Data_corrMatK.fmr1.Mean_corrMat{1,2}));caxis([-1 1]); colormap('jet');%% for 1st loom
subplot(1,8,counter+2);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,3},Data_corrMatK.fmr1.Mean_corrMat{1,3}));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+3);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,4},Data_corrMatK.fmr1.Mean_corrMat{1,4}));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+4);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,5},Data_corrMatK.fmr1.Mean_corrMat{1,5}));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+5);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,6},Data_corrMatK.fmr1.Mean_corrMat{1,6}));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+6);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,11},Data_corrMatK.fmr1.Mean_corrMat{1,11}));caxis([-1 1]);colormap('jet');%% for 10th loom
subplot(1,8,counter+7);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,12},Data_corrMatK.fmr1.Mean_corrMat{1,12})); caxis([-1 1]);colormap('jet');%% for 11th loom
sgtitle('controls minus fmr1');


%% to check the participation change from 10 to 11th 

%%% this needs to be done with the Brain Connectivity Toolbox. 
%%% so i need to have it in the path. 

W_thr10 = threshold_absolute(abs(Data_corrMatK.fmr1.Mean_corrMat{1,11}), 0.75);

 W_thr11 = threshold_absolute(abs(Data_corrMatK.fmr1.Mean_corrMat{1,12}), 0.75);

[Ppos10 Pneg10] = participation_coef_sign(W_thr10,Nodes.Mod_clust);

[Ppos11 Pneg11] = participation_coef_sign(W_thr11,Nodes.Mod_clust);


mean(Ppos10)

ratioPpos10_11=Ppos10./Ppos11;
subsPpos10_11=Ppos10-Ppos11;

figure;histogram(subsPpos10_11);

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),25,subsPpos10_11,'filled'); colorbar;colormap('jet')
view(-90,90);
title('fmr1');
%%% to check which cluster has the most nodes that increase participation
%%% at recovery. it seems is the slopehab. 
figure;histogram(Nodes.Mod_clust(find(subsPpos10_11<-0.3)));

figure;histogram(Nodes.Mod_clust(find(subsPpos10_11>0.3)));

%% checking a few things


%%% to see which ROIs have less representation
counter=1;
figure;
meanProp=[];
for g=1:3
    
  group=groupnames{g,1};
    
    fish=fieldnames(Nodes.NaNtest_K.(group));  
     
    
    Matrix_mean=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes.NaNtest_K.(group).(fish{f}));     
end

%Matrix_mean=(sum(Matrix_mean,3))/length(fish); %% is the same
Matrix_mean=nanmean(Matrix_mean,3);  


subplot(1,3,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),20,(1-mean(Matrix_mean)),'filled'); colorbar;colormap('jet'); caxis([0 1]);
view(-90,90);
title(group);
hold off;

counter=counter+1;

meanProp=horzcat(meanProp,(mean(Matrix_mean)'));


end

%%% to see which fish are contributing to specific nodes. 
fish=fieldnames(Nodes.NaNtest.fmr1);
fmr1_hindbrain_fish={};
for f=1:length(fish)
    
    temp_mat=Nodes.NaNtest.fmr1.(fish{f});
    if temp_mat(26,1)==0 
    fmr1_hindbrain_fish{f}=(fish{f})
    else
    end
end




%% what if i put a threshold on proportion of fish needed to contribute to a node?
%%%% before I run the crosscorrelation. 
%%% with at least 33% of the fish i loose some nodes in the hindbrain although a visual effect still
%%% seems to remain there. 

discard=find(min((1-meanProp)')<0.25); %%% I could use 0.25 (discard 6 nodes, but keep some in hindbrain), 0.33 (discard 11 nodes) , or 0/5 (discard 29 nodes).
keep=find(ismember([1:100],discard)==0);



figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes.Mod_loc(keep,1),Nodes.Mod_loc(keep,2),Nodes.Mod_clust(keep));view(-90,90);
title('Model Nodes');
%%% adding numbers
% a = [1:100]'; b = num2str(a(keep)); c = cellstr(b);
% dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
% text(x+dx, y+dy, c);


counter=1;
figure;
for g=1:3
     group=groupnames{g,1};
subplot(3,8,counter);imagesc(Data_corrMatK.(group).Mean_corrMat{1,1}(keep,keep)); %% for pre loom
subplot(3,8,counter+1);imagesc(Data_corrMatK.(group).Mean_corrMat{1,2}(keep,keep)); %% for 1st loom
subplot(3,8,counter+2);imagesc(Data_corrMatK.(group).Mean_corrMat{1,3}(keep,keep));
subplot(3,8,counter+3);imagesc(Data_corrMatK.(group).Mean_corrMat{1,4}(keep,keep));
subplot(3,8,counter+4);imagesc(Data_corrMatK.(group).Mean_corrMat{1,5}(keep,keep));
subplot(3,8,counter+5);imagesc(Data_corrMatK.(group).Mean_corrMat{1,6}(keep,keep));
subplot(3,8,counter+6);imagesc(Data_corrMatK.(group).Mean_corrMat{1,11}(keep,keep));%% for 10th loom
subplot(3,8,counter+7);imagesc(Data_corrMatK.(group).Mean_corrMat{1,12}(keep,keep)); %% for 11th loom
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
    
    R=Data_corrMatK.(group).Mean_corrMat{1,k}(keep,keep);
        
        
n=length(Nodes.Mod_loc(keep));

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
x = Nodes.Mod_loc(keep,1);
y = Nodes.Mod_loc(keep,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
gscatter(Nodes.Mod_loc(keep,1),Nodes.Mod_loc(keep,2),Nodes.Mod_clust(keep),'grbbggrm','.',20,'off');
hold on;
view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold off;

count=count+1

    end
end


%% for the substraction of contros minus fmr1
counter=1;
figure;    
subplot(1,8,counter);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,1}(keep,keep),Data_corrMatK.fmr1.Mean_corrMat{1,1}(keep,keep))); caxis([-1 1]);colormap('jet');%% for pre loom
subplot(1,8,counter+1);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,2}(keep,keep),Data_corrMatK.fmr1.Mean_corrMat{1,2}(keep,keep)));caxis([-1 1]); colormap('jet');%% for 1st loom
subplot(1,8,counter+2);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,3}(keep,keep),Data_corrMatK.fmr1.Mean_corrMat{1,3}(keep,keep)));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+3);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,4}(keep,keep),Data_corrMatK.fmr1.Mean_corrMat{1,4}(keep,keep)));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+4);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,5}(keep,keep),Data_corrMatK.fmr1.Mean_corrMat{1,5}(keep,keep)));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+5);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,6}(keep,keep),Data_corrMatK.fmr1.Mean_corrMat{1,6}(keep,keep)));caxis([-1 1]);colormap('jet');
subplot(1,8,counter+6);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,11}(keep,keep),Data_corrMatK.fmr1.Mean_corrMat{1,11}(keep,keep)));caxis([-1 1]);colormap('jet');%% for 10th loom
subplot(1,8,counter+7);imagesc(minus(Data_corrMatK.control.Mean_corrMat{1,12}(keep,keep),Data_corrMatK.fmr1.Mean_corrMat{1,12}(keep,keep))); caxis([-1 1]);colormap('jet');%% for 11th loom
sgtitle('controls minus fmr1');


%% to check the participation change from 10 to 11th 
W_thr10 = threshold_absolute(abs(Data_corrMatK.control.Mean_corrMat{1,11}(keep,keep)), 0.75);

 W_thr11 = threshold_absolute(abs(Data_corrMatK.control.Mean_corrMat{1,12}(keep,keep)), 0.75);

[Ppos10 Pneg10] = participation_coef_sign(W_thr10,Nodes.Mod_clust(keep));

[Ppos11 Pneg11] = participation_coef_sign(W_thr11,Nodes.Mod_clust(keep));


mean(Ppos10)

ratioPpos10_11=Ppos10./Ppos11;
subsPpos10_11=Ppos10-Ppos11;

figure;histogram(subsPpos10_11);

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes.Mod_loc(keep,1),Nodes.Mod_loc(keep,2),25,subsPpos10_11,'filled'); colorbar;colormap('jet')
view(-90,90);
title('control');
%%% to check which cluster has the most nodes that increase participation
%%% at recovery. it seems is the slopehab. 
figure;histogram(Nodes.Mod_clust(find(subsPpos10_11<-0.3)));

figure;histogram(Nodes.Mod_clust(find(subsPpos10_11>0.3)));


%%
%%%% it seems we will work with the atribution by the kmeans to the nodes. 

save('NodesNgraphFmr1Loomhab.mat','Nodes','goodorder_clust','Data_corrMatQ','Data_corrMatK');

