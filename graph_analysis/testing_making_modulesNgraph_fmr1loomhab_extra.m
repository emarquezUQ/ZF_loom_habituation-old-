%% new graph analysis with no filtered ROIs

%%% I will test now to get ROIs from all the fish but not filtered to
%%% secure the participation of most fish. 


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
    %temp_min=quantile(d,[0.025]); %%% mabye the best results. in varies among fish. what i used before for the wildtypes dataset
    %temp_min=quantile(d,[0.025])/2; %%% doesnt work that well
    temp_min=25; %%% I tested with 12... 25 (I got a very good fish representation but the graphs were not special
    
    %%%% i did some calculations on the number of NaNs i get with the
    %%%% different distances. maybe the min is very slightly better...
    %%%% although there is not much difference. 
    
      temp_idx_clust=tempfish;
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

    Nodes.ROIs_idx_nonloom.(group).(strcat('fish_',num2str(fish(f)))).(strcat('clust_',num2str(clust)))=idx_ROIs_Node;
    
    
  end
end

 %%% to check if i fixed the overlaping. it seemed that it worked. 
groupnames=fieldnames(Nodes.ROIs_idx_nonloom);

 for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_idx_nonloom.(group));
     
 for f=1:length(fish)
     
 for clust=goodorder_clust
   
    idx_ROIs_Node=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    
    temp_idx=Nodes.ROIs_idx_nonloom.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).idx;
    
    idx_ROIs_Node=vertcat(idx_ROIs_Node,temp_idx);
    end
     
    B = unique(idx_ROIs_Node); % which will give you the unique elements of A in array B
    Ncount = histc(idx_ROIs_Node, B); % this willgive the number of occurences of each unique element
    unique(Ncount)
    
    
end
end
 end



for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_idx_nonloom.(group));
     
     temp_allfish=[];
 for f=1:length(fish)
     
     temp_Num_ROIs_Node=struct;
     
 for clust=goodorder_clust
   
    
    for node=(find(Nodes.Mod_clust==clust))'
        
        temp_roiN=size(Nodes.ROIs_idx_nonloom.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).idx,1);                  
        
        temp_Num_ROIs_Node.(group).(fish{f})(node,1)=temp_roiN;
    end   
            
 end
 
 temp_allfish=horzcat(temp_allfish,temp_Num_ROIs_Node.(group).(fish{f}));
 
 end
  
 Nodes.NumROIsnonloom.(group)=temp_allfish; 
 
 
end

%%
%%% i then analyze them in prism there were no sig differences in any
%%% nodes between controls and fmr1


%%
 %%% to get the means
  for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_idx_nonloom.(group));
 for f=1:length(fish)
     
for clust=goodorder_clust
   

    for node=(find(Nodes.Mod_clust==clust))'
    
   
    temp_idx=Nodes.ROIs_idx_nonloom.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).idx;
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
    
    Nodes.ROIs_idx_nonloom.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

  end
  
  
  %%% making means matrices per fish and testing a correlation
  for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_idx_nonloom.(group));
for f=1:length(fish)

    temp_matrix=NaN(100,904);
    for clust=goodorder_clust
        
        for node=(find(Nodes.Mod_clust==clust))'
            
        temp_mean=Nodes.ROIs_idx_nonloom.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes.mean_matrix_raw.(group).(fish{f})=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes.corr_matrix_raw.(group).(fish{f})=R_temp;
        
        Nodes.NaNtest_raw.(group).(fish{f})=isnan(R_temp);
        
end
  end

  %%% testing the matrix with all fish to see if there are still gaps
Matrix_mean=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes.NaNtest_raw.(group).(fish{f}));     
end

Matrix_mean=nanmean(Matrix_mean,3);
  

%%% is it better to use less distance for recluting the ROIs? cause maybe
%%% sometimes one node 'steals' from another one. I did some test and I got
%%% the best results with the low quantile, although very close to the min(d).
%%% avg of NaN for the whole matrix: min=0.4251 and quantile=0.4179 (mean of NaNs)

mean(mean(Matrix_mean))


%% getting the corrmatrix for each loom and each fish
Data_corrMatRaw=struct;

for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.mean_matrix_raw.(group));
    for f=1:length(fish)
    
        temp_mean=Nodes.mean_matrix_raw.(group).(fish{f});
        
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
        
    Data_corrMatRaw.(group).(fish{f}).loomsR=temp_R;
        
    end
             
end

%% making means of each loom per dataset

for g=1:3
     group=groupnames{g,1};
    
    fish=fieldnames(Data_corrMatRaw.(group));
    
    Mean_corrMat={};
    for k=1:21
        
        temp_Mean_corrMat=[];
        
        for f=1:length(fish)
            
        temp_mat=Data_corrMatRaw.(group).(fish{f}).loomsR{1,k};
        
        temp_Mean_corrMat=cat(3,temp_Mean_corrMat,temp_mat);
        end
        
        Mean_corrMat{1,k}=nanmean(temp_Mean_corrMat,3);
    end
    Data_corrMatRaw.(group).Mean_corrMat=Mean_corrMat;
    
end


%%
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Nodes.Mod_clust);view(-90,90);
title('Model Nodes');

counter=1;
figure;
for g=1:3
     group=groupnames{g,1};
subplot(3,8,counter);imagesc(Data_corrMatRaw.(group).Mean_corrMat{1,1}); caxis([-1 1]);%% for pre loom
subplot(3,8,counter+1);imagesc(Data_corrMatRaw.(group).Mean_corrMat{1,2});caxis([-1 1]); %% for 1st loom
subplot(3,8,counter+2);imagesc(Data_corrMatRaw.(group).Mean_corrMat{1,3});caxis([-1 1]);
subplot(3,8,counter+3);imagesc(Data_corrMatRaw.(group).Mean_corrMat{1,4});caxis([-1 1]);
subplot(3,8,counter+4);imagesc(Data_corrMatRaw.(group).Mean_corrMat{1,5});caxis([-1 1]);
subplot(3,8,counter+5);imagesc(Data_corrMatRaw.(group).Mean_corrMat{1,6});caxis([-1 1]);
subplot(3,8,counter+6);imagesc(Data_corrMatRaw.(group).Mean_corrMat{1,11});caxis([-1 1]);%% for 10th loom
subplot(3,8,counter+7);imagesc(Data_corrMatRaw.(group).Mean_corrMat{1,12});caxis([-1 1]); %% for 11th loom
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
    
    R=Data_corrMatRaw.(group).Mean_corrMat{1,k};
        
        
n=length(Nodes.Mod_loc);

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

[~,~,weights] = find(tril(R,-1));

% create the graph object:
G = graph(s,t,weights,n);

% mark the lines to remove from the graph:
threshold = 0.5; %  minimum correlation to plot
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


%% checking a few things


%%% to see which ROIs have less representation
counter=1;
meanProp=[];
figure;
for g=1:3
    
  group=groupnames{g,1};
    
    fish=fieldnames(Nodes.NaNtest_raw.(group));  
     
    
    Matrix_mean=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes.NaNtest_raw.(group).(fish{f}));     
end

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

%%
%% new graph analysis with no filtered ROIs first and then correlation top 10.

%%% I will test now to get ROIs from all the fish but not filtered to
%%% secure the participation of most fish. 


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
    %temp_min=quantile(d,[0.025]); %%% mabye the best results. in varies among fish. what i used before for the wildtypes dataset
    %temp_min=quantile(d,[0.025])/2; %%% doesnt work that well
    temp_min=25; %%% I tested with 12... 25 (I got a very good fish representation but the graphs were not special
    
    %%%% i did some calculations on the number of NaNs i get with the
    %%%% different distances. maybe the min is very slightly better...
    %%%% although there is not much difference. 
    
      temp_idx_clust=tempfish;
      %figure;plot(mean(ZS_CN(temp,:))); %% to check if it works.
       
      D = pdist2(ROI_temp2(temp_idx_clust,:),Nodes.Mod_loc);  

    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    [temp_row,~]=find(D(:,node)<temp_min);
    
    %%% this is to filter them with a correlation to the cluster mean
    if length(temp_row)>10
        corr_temp=[];
    for i=1:length(temp_row)    
    temp_corr=corrcoef(Cmap_ZS_CN(GoodBetas_ZS_CN_selected(find(goodorder_clust==clust)),:),ZS_CN(temp_idx_clust(temp_row(i)),:));
    corr_temp(i)=temp_corr(1,2);
    end
    sorted_corr=sort(corr_temp,'descend');
    temp_row=temp_row(find(ismember(corr_temp,sorted_corr(1:10))));
    else
    end
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

    Nodes.ROIs_idx_nonloom_corr.(group).(strcat('fish_',num2str(fish(f)))).(strcat('clust_',num2str(clust)))=idx_ROIs_Node;
    
    
  end
end

 %%% to check if i fixed the overlaping. it seemed that it worked. 
groupnames=fieldnames(Nodes.ROIs_idx_nonloom_corr);

 for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_idx_nonloom_corr.(group));
     
 for f=1:length(fish)
     
 for clust=goodorder_clust
   
    idx_ROIs_Node=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    
    temp_idx=Nodes.ROIs_idx_nonloom_corr.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).idx;
    
    idx_ROIs_Node=vertcat(idx_ROIs_Node,temp_idx);
    end
     
    B = unique(idx_ROIs_Node); % which will give you the unique elements of A in array B
    Ncount = histc(idx_ROIs_Node, B); % this willgive the number of occurences of each unique element
    unique(Ncount)
    
    
end
end
 end



for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_idx_nonloom_corr.(group));
     
     temp_allfish=[];
 for f=1:length(fish)
     
     temp_Num_ROIs_Node=struct;
     
 for clust=goodorder_clust
   
    
    for node=(find(Nodes.Mod_clust==clust))'
        
        temp_roiN=size(Nodes.ROIs_idx_nonloom_corr.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).idx,1);                  
        
        temp_Num_ROIs_Node.(group).(fish{f})(node,1)=temp_roiN;
    end   
            
 end
 
 temp_allfish=horzcat(temp_allfish,temp_Num_ROIs_Node.(group).(fish{f}));
 
 end
  
 Nodes.NumROIsnonloom_corr.(group)=temp_allfish; 
 
 
end

%%
%%% i then analyze them in prism there were no sig differences in any
%%% nodes between controls and fmr1


%%
 %%% to get the means
  for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_idx_nonloom_corr.(group));
 for f=1:length(fish)
     
for clust=goodorder_clust
   

    for node=(find(Nodes.Mod_clust==clust))'
    
   
    temp_idx=Nodes.ROIs_idx_nonloom_corr.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).idx;
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
    
    Nodes.ROIs_idx_nonloom_corr.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

  end
  
  
  %%% making means matrices per fish and testing a correlation
  for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_idx_nonloom_corr.(group));
for f=1:length(fish)

    temp_matrix=NaN(100,904);
    for clust=goodorder_clust
        
        for node=(find(Nodes.Mod_clust==clust))'
            
        temp_mean=Nodes.ROIs_idx_nonloom_corr.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes.mean_matrix_raw_corr.(group).(fish{f})=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes.corr_matrix_raw_corr.(group).(fish{f})=R_temp;
        
        Nodes.NaNtest_raw_corr.(group).(fish{f})=isnan(R_temp);
        
end
  end

  %%% testing the matrix with all fish to see if there are still gaps
Matrix_mean=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes.NaNtest_raw_corr.(group).(fish{f}));     
end

Matrix_mean=nanmean(Matrix_mean,3);
  

%%% is it better to use less distance for recluting the ROIs? cause maybe
%%% sometimes one node 'steals' from another one. I did some test and I got
%%% the best results with the low quantile, although very close to the min(d).
%%% avg of NaN for the whole matrix: min=0.4251 and quantile=0.4179 (mean of NaNs)

mean(mean(Matrix_mean))


%% getting the corrmatrix for each loom and each fish
Data_corrMatRaw_corr=struct;

for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.mean_matrix_raw_corr.(group));
    for f=1:length(fish)
    
        temp_mean=Nodes.mean_matrix_raw_corr.(group).(fish{f});
        
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
        
    Data_corrMatRaw_corr.(group).(fish{f}).loomsR=temp_R;
        
    end
             
end

%% making means of each loom per dataset

for g=1:3
     group=groupnames{g,1};
    
    fish=fieldnames(Data_corrMatRaw_corr.(group));
    
    Mean_corrMat={};
    for k=1:21
        
        temp_Mean_corrMat=[];
        
        for f=1:length(fish)
            
        temp_mat=Data_corrMatRaw_corr.(group).(fish{f}).loomsR{1,k};
        
        temp_Mean_corrMat=cat(3,temp_Mean_corrMat,temp_mat);
        end
        
        Mean_corrMat{1,k}=nanmean(temp_Mean_corrMat,3);
    end
    Data_corrMatRaw_corr.(group).Mean_corrMat=Mean_corrMat;
    
end


%%
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Nodes.Mod_clust);view(-90,90);
title('Model Nodes');

counter=1;
figure;
for g=1:3
     group=groupnames{g,1};
subplot(3,8,counter);imagesc(Data_corrMatRaw_corr.(group).Mean_corrMat{1,1}); caxis([-1 1]);%% for pre loom
subplot(3,8,counter+1);imagesc(Data_corrMatRaw_corr.(group).Mean_corrMat{1,2});caxis([-1 1]); %% for 1st loom
subplot(3,8,counter+2);imagesc(Data_corrMatRaw_corr.(group).Mean_corrMat{1,3});caxis([-1 1]);
subplot(3,8,counter+3);imagesc(Data_corrMatRaw_corr.(group).Mean_corrMat{1,4});caxis([-1 1]);
subplot(3,8,counter+4);imagesc(Data_corrMatRaw_corr.(group).Mean_corrMat{1,5});caxis([-1 1]);
subplot(3,8,counter+5);imagesc(Data_corrMatRaw_corr.(group).Mean_corrMat{1,6});caxis([-1 1]);
subplot(3,8,counter+6);imagesc(Data_corrMatRaw_corr.(group).Mean_corrMat{1,11});caxis([-1 1]);%% for 10th loom
subplot(3,8,counter+7);imagesc(Data_corrMatRaw_corr.(group).Mean_corrMat{1,12});caxis([-1 1]); %% for 11th loom
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
    
    R=Data_corrMatRaw_corr.(group).Mean_corrMat{1,k};
        
        
n=length(Nodes.Mod_loc);

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

[~,~,weights] = find(tril(R,-1));

% create the graph object:
G = graph(s,t,weights,n);

% mark the lines to remove from the graph:
threshold = 0.5; %  minimum correlation to plot
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


%% checking a few things


%%% to see which ROIs have less representation
counter=1;
meanProp=[];
figure;
for g=1:3
    
  group=groupnames{g,1};
    
    fish=fieldnames(Nodes.NaNtest_raw_corr.(group));  
     
    
    Matrix_mean=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes.NaNtest_raw_corr.(group).(fish{f}));     
end

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

%%

%%%% checking the cluster distribution per fish for the fmr1 loomhab
%%%% dataset



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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

ClustDistPerF=struct;
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
          
      temp_idx_clust=intersect(idx_rsq_cleaned(find(High_corr_Nb==clust)),tempfish);
       
    ClustDistPerF.ROIs_idx.(group).(strcat('fish_',num2str(fish(f)))).(strcat('clust_',num2str(clust)))=temp_idx_clust;
    
    
  end
end



%%%% to check if it worked
temp=ClustDistPerF.ROIs_idx.(group).(strcat('fish_',num2str(fish(f)))).(strcat('clust_',num2str(clust)));

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(ROI_temp2(temp,1),ROI_temp2(temp,2),15,'filled');
view(-90,90);


for g=1:3
     group=groupnames{g,1};
     
     for clust=6%goodorder_clust
      
     figure;
     counter=1;
     sgtitle(strcat(group,'_',num2str(clust)));
     fish=fieldnames(ClustDistPerF.ROIs_idx.(group));
              
        for f=1:length(fish)
            
        if length(fish)<12 
        row=3;col=4;
        
        else
        row=4;col=5;    
        end
        subplot(row,col,counter);
        temp=ClustDistPerF.ROIs_idx.(group).(fish{f}).(strcat('clust_',num2str(clust)));    
        plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
        hold on;
        scatter(ROI_temp2(temp,1),ROI_temp2(temp,2),15,'filled');
        view(-90,90); hold off;
        
        counter=counter+1;
        
        
        end   
      end
     
end

%% what if I lower the r2 threshold? 
%%% i cant do it much cause i did the linear regression to a pre filtered
%%% dataset with a correlation of 

load('rsquare_fmr1loomhab.mat');

load('s20_fmr1_loomhab_CN_part2.mat','idx_corr');
 


idx_rsq=idx_corr(idx_rsq_test_s20short);


%%%% unless I rerun my correlation

%%% i will use this clusters to filter the data with a correlation of 0.5

load('s20_postKmeans_CN.mat','idxKmeans_ZS_CN','Cmap_ZS_CN');
load('s20_fmr1_loomhab_CN_part2.mat','GoodBetas_ZS_CN_selected','idx_rsq');


for i=1:length(GoodBetas_ZS_CN_selected)
GoodBetas_regress(i,:)=Cmap_ZS_CN(GoodBetas_ZS_CN_selected(i),:);

end


corrfilter=[];Threshold=0.3;
for i=1:length(GoodBetas_ZS_CN_selected)    
    
    corr_temp=zeros(1,length(ZS_CN));
    parfor jj=1:size(ZS_CN,1)
        temp_corr=corrcoef(GoodBetas_regress(i,:), ZS_CN(jj,:));
        corr_temp(jj)=temp_corr(1,2);
    end    
    %corrfilter(i).ZS=ZS_CN(find(corr_temp>Threshold),:);
    corrfilter(i).idx=find(corr_temp>Threshold);
    %corrfilter(i).mean=mean(corrfilter(i).ZS,1);
    %corrfilter(i).STD=std(corrfilter(i).ZS,1,1);       
end

idx_corr=[];
for i=1:length(GoodBetas_ZS_CN_selected)    

idx_corr=horzcat(idx_corr,corrfilter(i).idx);

end    
   
idx_corr=unique(idx_corr);

%%% to clasify them with a correlation

Correlation_group={};
counter=1;
for i=1:length(GoodBetas_ZS_CN_selected)
    Correlation_temp=[];
    for idx=1:size(ZS_CN(idx_corr),2)
        temp_corr=corrcoef(GoodBetas_regress(i,:),ZS_CN(idx_corr(idx),:));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat=[];
for n=1:size(GoodBetas_ZS_CN_selected,2)
Correlation_group_mat(n,:)=Correlation_group{n}(n,:);

end

%%% now if I filter looking for the max correlation
%%% 

High_corr_Nb_corr=zeros(length(Correlation_group_mat),1);
for i=1:length(Correlation_group_mat)
    [~,I]=max(Correlation_group_mat(:,i));
    High_corr_Nb_corr(i,1)=I;
    
end


%%
fish=vertcat(list1,list2,list3,list4);
list5=union(list1,list3);

%%% it seems that fish 201810048 from list1 dont have any ROIs... dont know
%%% why. might need to check if it is not a mistake in the name. i am
%%% getting rid of it on the list. 
find(idx_Fish==201810048);
fish(find(fish==201810048))=[];

ClustDistPerF2=struct;
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
          
      temp_idx_clust=intersect(idx_corr(find(High_corr_Nb_corr==clust)),tempfish);
       
    ClustDistPerF2.ROIs_idx.(group).(strcat('fish_',num2str(fish(f)))).(strcat('clust_',num2str(clust)))=temp_idx_clust;
    
    
  end
end

%% to plot

for g=1:3
     group=groupnames{g,1};
     
     for clust=5%goodorder_clust(1:3)%% only fasthab
      
     figure;
     counter=1;
     sgtitle(strcat(group,'_',num2str(clust)));
     fish=fieldnames(ClustDistPerF2.ROIs_idx.(group));
              
        for f=1:length(fish)
            
        if length(fish)<12 
        row=3;col=4;
        
        else
        row=4;col=5;    
        end
        subplot(row,col,counter);
        temp=ClustDistPerF2.ROIs_idx.(group).(fish{f}).(strcat('clust_',num2str(clust)));  
        temp2=ClustDistPerF.ROIs_idx.(group).(fish{f}).(strcat('clust_',num2str(clust))); 
        plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
        hold on;
        scatter(ROI_temp2(temp,1),ROI_temp2(temp,2),15,'filled');
        hold on;
        scatter(ROI_temp2(temp2,1),ROI_temp2(temp2,2),15,'filled');
        view(-90,90); hold off;
        
        counter=counter+1;
        
        
        end   
      end
     
end

%%%% checking how many ROIs are there in the hindbrain. per genotype. taking
%%%% all ROIs not only loom respondent. 


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
    
      temp_idx_clust=tempfish;
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

    Nodes.ROIs_idx_nonloom.(group).(strcat('fish_',num2str(fish(f)))).(strcat('clust_',num2str(clust)))=idx_ROIs_Node;
    
    
  end
end

 %%% to check if i fixed the overlaping. it seemed that it worked. 
groupnames=fieldnames(Nodes.ROIs_idx_nonloom);

 for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_idx_nonloom.(group));
     
 for f=1:length(fish)
     
 for clust=goodorder_clust
   
    idx_ROIs_Node=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    
    temp_idx=Nodes.ROIs_idx_nonloom.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).idx;
    
    idx_ROIs_Node=vertcat(idx_ROIs_Node,temp_idx);
    end
     
    B = unique(idx_ROIs_Node); % which will give you the unique elements of A in array B
    Ncount = histc(idx_ROIs_Node, B); % this willgive the number of occurences of each unique element
    unique(Ncount)
    
    
end
end
 end



for g=1:3
     group=groupnames{g,1};
     fish=fieldnames(Nodes.ROIs_idx_nonloom.(group));
     
     temp_allfish=[];
 for f=1:length(fish)
     
     temp_Num_ROIs_Node=struct;
     
 for clust=goodorder_clust
   
    
    for node=(find(Nodes.Mod_clust==clust))'
        
        temp_roiN=size(Nodes.ROIs_idx_nonloom.(group).(fish{f}).(strcat('clust_',num2str(clust))).(strcat('node_',num2str(node))).idx,1);                  
        
        temp_Num_ROIs_Node.(group).(fish{f})(node,1)=temp_roiN;
    end   
            
 end
 
 temp_allfish=horzcat(temp_allfish,temp_Num_ROIs_Node.(group).(fish{f}));
 
 end
  
 Nodes.NumROIsnonloom.(group)=temp_allfish; 
 
 
end


%%%% i then analyse it in prism. It seems that only a couple of nodes have
%%%% significant differents. although among them node 2 which is in the
%%%% hindbrain. 


%%%% in this case I got them based on the minimum quantile distance... but
%%%% there might be a lot of overlap. so i will make a new analysis with
%%%% smaller area. 



