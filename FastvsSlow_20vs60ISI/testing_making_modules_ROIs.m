%%%% this script is to play around with the data to check some of the
%%%% methods: Organizing principles of whole-brain functional connectivity in zebrafish larvae 
%%%% Betzel 2018 biorxiv

%%%I am using the s20 dataset to test

%%% is to make functional modules in the zebrafish brain. 

options = statset('UseParallel',1); [idxKmeans_ROIs_s20 Cmap_ROIs_s20]=kmeans(ROI_temp2.s20(idx_rsq_test_s20short_cleaned,:),100,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');

patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
%figure; 
scatter3(Cmap_ROIs_s20(:,1),Cmap_ROIs_s20(:,2),Cmap_ROIs_s20(:,3));


module_means=[];
for i=1:length(Cmap_ROIs_s20)
    temp_idx=find(idxKmeans_ROIs_s20==i);
    module_means(i,:)=mean(ZS_s20(idx_rsq_test_s20short_cleaned(temp_idx),:));
    
end

%%% to check if it worked
figure;
for i=1:length(Cmap_ROIs_s20)
plot(module_means(i,:));
pause(2)
end

%%% to make a matrix of the corr
R = corrcoef(module_means');


R2=atanh(R); %%% this is to do a fisher transform it doesnt work very well with it... 
%%%% i think is actually to normalize different correlation results from
%%%% different samples. so in this case it doesnt work.

%%% here i am trying to cluster the correlation matrix
options = statset('UseParallel',1); [idxKmeans_modules_corr Cmap_modules_corr]=kmeans(R,4,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');


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

scatter(Cmap_ROIs_s20(:,1),Cmap_ROIs_s20(:,2));


%%% this is to plot them per groups but i dont know is not working
gscatter(Cmap_ROIs_s20(:,1),Cmap_ROIs_s20(:,2),idxKmeans_modules_corr,'r');

%%%% so i will do it with a loop

figure;
for i=1:max(idxKmeans_modules_corr)
temp_idx=find(idxKmeans_modules_corr==i);
scatter(Cmap_ROIs_s20(temp_idx,1),Cmap_ROIs_s20(temp_idx,2),'filled');
hold on;
end


%%%% now just trying to see what happens if I kluster the modules 
options = statset('UseParallel',1); [idxKmeans_modules_mean Cmap_modules_mean]=kmeans(module_means,4,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');



%%% to check if it worked
figure;
for i=1:max(idxKmeans_modules_mean)
subplot(4,2,i);
    plot(Cmap_modules_mean(i,:));
end

figure; 
for i=1:max(idxKmeans_modules_mean)
temp_idx=find(idxKmeans_modules_mean==i);
scatter(Cmap_ROIs_s20(temp_idx,1),Cmap_ROIs_s20(temp_idx,2),10,'filled');
hold on;
end


figure;
patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
%figure; 

for i=1:max(idxKmeans_modules_mean)
temp_idx=find(idxKmeans_modules_mean==i);
scatter3(Cmap_ROIs_s20(temp_idx,1),Cmap_ROIs_s20(temp_idx,2),Cmap_ROIs_s20(temp_idx,3),10,'filled');
hold on;
end

%%% this is not working
%gscatter(Cmap_ROIs_s20(:,1),Cmap_ROIs_s20(:,2),idxKmeans_modules_mean);

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
x = Cmap_ROIs_s20(:,1);
y = Cmap_ROIs_s20(:,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix


hold on;
for i=1:max(idxKmeans_modules_mean)
temp_idx=find(idxKmeans_modules_mean==i);
scatter(Cmap_ROIs_s20(temp_idx,1),Cmap_ROIs_s20(temp_idx,2),'filled');
hold on;
end


%%%% 3d doesnt work cause the network is in 2d
hold on;
patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
hold on;
%figure; 

for i=1:max(idxKmeans_modules_mean)
temp_idx=find(idxKmeans_modules_mean==i);
scatter3(Cmap_ROIs_s20(temp_idx,1),Cmap_ROIs_s20(temp_idx,2),Cmap_ROIs_s20(temp_idx,3),5,'filled');
hold on;
end
