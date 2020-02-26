%%%% this script is to make a correlation in each fish for each of the loom
%%%% resonses based on the nodes from
%%%% testing_making_modules_ROIs_datasets_perfish.m. Then I will make an
%%%% average of these correlation matrices per dataset to be able to
%%%% compare them. 

datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);

load('Nodes_N_means_alldatasets.mat','Nodes','Zbrain_brainMask2D');

load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx');

%% getting the corrmatrix for each loom and each fish
Data_corrMat=struct;

for data=1:4
    datatemp=datasets(data,:);
    
    fish=fieldnames(Nodes.(datatemp).mean_matrix);
    
    for tempfish=1:length(fish)
        temp_mean=Nodes.(datatemp).mean_matrix.(fish{tempfish});
        
        temp_R={};
        
        for k=1:31
        
        if k==1
            temp_R{1,k}=corrcoef((temp_mean(:,10:45))');
        elseif k==11||k==21 ||k==31
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28))');
        else
            
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36))'); 
                           
        end
        end
        
    Data_corrMat.(datatemp).(fish{tempfish}).loomsR=temp_R;
        
    end
             
end


%% making means of each loom per dataset

for data=1:4
    datatemp=datasets(data,:);
    
    fish=fieldnames(Data_corrMat.(datatemp));
    
    Mean_corrMat={};
    for k=1:31
        
        temp_Mean_corrMat=[];
        
        for tempfish=1:length(fish)
            
        temp_mat=Data_corrMat.(datatemp).(fish{tempfish}).loomsR{1,k};
        
        temp_Mean_corrMat=cat(3,temp_Mean_corrMat,temp_mat);
        end
        
        Mean_corrMat{1,k}=nanmean(temp_Mean_corrMat,3);
    end
    Data_corrMat.(datatemp).Mean_corrMat=Mean_corrMat;
    
end




%%
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Nodes.Mod_clust);view(-90,90);
title('Model Nodes');

counter=1;
figure;
for data=1:4
    datatemp=datasets(data,:);
subplot(4,8,counter);imagesc(Data_corrMat.(datatemp).Mean_corrMat{1,1}); %% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat.(datatemp).Mean_corrMat{1,2}); %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat.(datatemp).Mean_corrMat{1,3});
subplot(4,8,counter+3);imagesc(Data_corrMat.(datatemp).Mean_corrMat{1,4});
subplot(4,8,counter+4);imagesc(Data_corrMat.(datatemp).Mean_corrMat{1,5});
subplot(4,8,counter+5);imagesc(Data_corrMat.(datatemp).Mean_corrMat{1,6});
subplot(4,8,counter+6);imagesc(Data_corrMat.(datatemp).Mean_corrMat{1,11});%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat.(datatemp).Mean_corrMat{1,12}); %% for 11th loom
title(datatemp);

counter=counter+8;
end


%% now to plot some graphs. 




for data=1:4
    figure;
    count=1;
    datatemp=datasets(data,:);
    
   
    for k=[2 3 4 5 6 11 12]
    
    R=Data_corrMat.(datatemp).Mean_corrMat{1,k};
        
        
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
gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Nodes.Mod_clust,'rygcbm','.',20,'off');

view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold off;

count=count+1

    end
end


save('graph_loomNdataset.mat','Data_corrMat');


%% checking a few things

for data=1:4
     datatemp=datasets(data,:);
        
     fish=fieldnames(Nodes.(datatemp).corr_matrix);
     
    for tempfish=1:length(fish)
           
    Nodes.NaNtest.(datatemp).(fish{tempfish})=isnan(Nodes.(datatemp).corr_matrix.(fish{tempfish}));
    end
end


%%% to see which ROIs have less representation
counter=1;
figure;
meanProp=[];
for data=1:4
    datatemp=datasets(data,:);
    
    fish=fieldnames(Nodes.NaNtest.(datatemp));
     
    
    Matrix_mean=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes.NaNtest.(datatemp).(fish{f}));     
end

%Matrix_mean=(sum(Matrix_mean,3))/length(fish); %% is the same
Matrix_mean=nanmean(Matrix_mean,3);  


subplot(1,4,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),20,(1-mean(Matrix_mean)),'filled'); colorbar;colormap('jet'); caxis([0 1]);
view(-90,90);
title(datatemp);
hold off;

counter=counter+1;

meanProp=horzcat(meanProp,(mean(Matrix_mean)'));


end
