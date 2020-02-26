
%%%% now between f20 and f60. 



F20vsF60_disN=struct;

for i=1:length(clustersF)

idx_clust=find(strcmp(clustersF, clustersF(i)));   
temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
temp_fieldname2=fieldnames(f60_cleaned_idxs.clust_f60_CL4_cleaned);temp_fieldname2=temp_fieldname2(idx_clust);


%%% to get the distances
D_f20 = pdist2(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),:),ROI_temp2.f60(f60_cleaned_idxs.clust_f60_CL4_cleaned.(char(temp_fieldname2)),:),'euclidean');

M_f20=struct;

for j=1:size(D_f20,1)
        
temp_D=D_f20(j,:); %%%% here I am getting each row to be albe to take way the 0 of the ROIs that is being comapared to itself
temp_D(find(temp_D==0))=mean(temp_D); %%% here I am changing the 0 value to the mean of that row. just to be albe to get the min afterwards without getting the 0
[min_temp idx_temp]=min(temp_D(temp_D(1,:)>0)); %%% to get the nearest ROI and its location
M_f20.min(j)=min_temp;
M_f20.idx(j)=idx_temp;

end

%%% to put all the data together

F20vsF60_disN.distance.(char(clustersF(i))).Dis=D_f20;
F20vsF60_disN.distance.(char(clustersF(i))).M.min=M_f20.min;
F20vsF60_disN.distance.(char(clustersF(i))).M.idx=M_f20.idx;


%%%% to get the descriptive statistics of the distances and corr
F20vsF60_disN.distance.(char(clustersF(i))).M.Min_dis = min(M_f20.min);
F20vsF60_disN.distance.(char(clustersF(i))).M.Max_dis = max(M_f20.min);
F20vsF60_disN.distance.(char(clustersF(i))).M.mean_dis = mean(M_f20.min);
F20vsF60_disN.distance.(char(clustersF(i))).M.med_dis = median(M_f20.min);
F20vsF60_disN.distance.(char(clustersF(i))).M.sd_dis = std(M_f20.min);
F20vsF60_disN.distance.(char(clustersF(i))).M.Qs_dis = quantile(M_f20.min,[0.025 0.25 0.50 0.75 0.975]);

end



%%
%%% now for s20


F20vsS20_disN=struct;

for i=1:length(clustersF)

idx_clust=find(strcmp(clustersF, clustersF(i)));    
idx_clust2=find(strcmp(clustersS, clustersF(i)));   
temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
temp_fieldname2=fieldnames(s20_cleaned_idxs.clust_s20_CL4_cleaned);temp_fieldname2=temp_fieldname2(idx_clust2);


%%% to get the distances
D_f20 = pdist2(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),:),ROI_temp2.s20(s20_cleaned_idxs.clust_s20_CL4_cleaned.(char(temp_fieldname2)),:),'euclidean');

M_f20=struct;

for j=1:size(D_f20,1)
        
temp_D=D_f20(j,:); %%%% here I am getting each row to be albe to take way the 0 of the ROIs that is being comapared to itself
temp_D(find(temp_D==0))=mean(temp_D); %%% here I am changing the 0 value to the mean of that row. just to be albe to get the min afterwards without getting the 0
[min_temp idx_temp]=min(temp_D(temp_D(1,:)>0)); %%% to get the nearest ROI and its location
M_f20.min(j)=min_temp;
M_f20.idx(j)=idx_temp;

end



%%% to put all the data together

F20vsS20_disN.distance.(char(clustersF(i))).Dis=D_f20;
F20vsS20_disN.distance.(char(clustersF(i))).M.min=M_f20.min;
F20vsS20_disN.distance.(char(clustersF(i))).M.idx=M_f20.idx;



%%%% to get the descriptive statistics of the distances and corr
F20vsS20_disN.distance.(char(clustersF(i))).M.Min_dis = min(M_f20.min);
F20vsS20_disN.distance.(char(clustersF(i))).M.Max_dis = max(M_f20.min);
F20vsS20_disN.distance.(char(clustersF(i))).M.mean_dis = mean(M_f20.min);
F20vsS20_disN.distance.(char(clustersF(i))).M.med_dis = median(M_f20.min);
F20vsS20_disN.distance.(char(clustersF(i))).M.sd_dis = std(M_f20.min);
F20vsS20_disN.distance.(char(clustersF(i))).M.Qs_dis = quantile(M_f20.min,[0.025 0.25 0.50 0.75 0.975]);



end


%%

%%% now for s60


F20vsS60_disN=struct;

for i=1:length(clustersF)

idx_clust=find(strcmp(clustersF, clustersF(i)));    
idx_clust2=find(strcmp(clustersS, clustersF(i)));   
temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
temp_fieldname2=fieldnames(s60_cleaned_idxs.clust_s60_CL4_cleaned);temp_fieldname2=temp_fieldname2(idx_clust2);


%%% to get the distances
D_f20 = pdist2(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),:),ROI_temp2.s60(s60_cleaned_idxs.clust_s60_CL4_cleaned.(char(temp_fieldname2)),:),'euclidean');

M_f20=struct;

for j=1:size(D_f20,1)
        
temp_D=D_f20(j,:); %%%% here I am getting each row to be albe to take way the 0 of the ROIs that is being comapared to itself
temp_D(find(temp_D==0))=mean(temp_D); %%% here I am changing the 0 value to the mean of that row. just to be albe to get the min afterwards without getting the 0
[min_temp idx_temp]=min(temp_D(temp_D(1,:)>0)); %%% to get the nearest ROI and its location
M_f20.min(j)=min_temp;
M_f20.idx(j)=idx_temp;

end


%%% to put all the data together

F20vsS60_disN.distance.(char(clustersF(i))).Dis=D_f20;
F20vsS60_disN.distance.(char(clustersF(i))).M.min=M_f20.min;
F20vsS60_disN.distance.(char(clustersF(i))).M.idx=M_f20.idx;


%%%% to get the descriptive statistics of the distances and corr
F20vsS60_disN.distance.(char(clustersF(i))).M.Min_dis = min(M_f20.min);
F20vsS60_disN.distance.(char(clustersF(i))).M.Max_dis = max(M_f20.min);
F20vsS60_disN.distance.(char(clustersF(i))).M.mean_dis = mean(M_f20.min);
F20vsS60_disN.distance.(char(clustersF(i))).M.med_dis = median(M_f20.min);
F20vsS60_disN.distance.(char(clustersF(i))).M.sd_dis = std(M_f20.min);
F20vsS60_disN.distance.(char(clustersF(i))).M.Qs_dis = quantile(M_f20.min,[0.025 0.25 0.50 0.75 0.975]);



end

%%

%%%% now try to find the ROIs that are close in all datasets. 

limit=5;

All_disN=struct;

for i=1:length(clustersF)

    
idx_clust=find(strcmp(clustersF, clustersF(i)));    
idx_clust2=find(strcmp(clustersS, clustersF(i)));   
temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);
temp_fieldname2=fieldnames(s60_cleaned_idxs.clust_s60_CL4_cleaned);temp_fieldname2=temp_fieldname2(idx_clust2);
   
    
    
idx_close_temp=find(F20vsF60_disN.distance.(char(clustersF(i))).M.min<limit);
All_disN.(char(clustersF(i))).idx_f60_close=f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname))(idx_close_temp);


idx_close_temp=find(F20vsS20_disN.distance.(char(clustersF(i))).M.min<limit);
All_disN.(char(clustersF(i))).idx_s20_close=f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname))(idx_close_temp);
    
    
idx_close_temp=find(F20vsS60_disN.distance.(char(clustersF(i))).M.min<limit);
All_disN.(char(clustersF(i))).idx_s60_close=f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname))(idx_close_temp);


end




figure;
for j=1:length(clustersF)
    

scatter(ROI_temp2.f20(All_disN.(char(clustersF(j))).idx_s60_close,1),ROI_temp2.f20(All_disN.(char(clustersF(j))).idx_s60_close,2),10,'filled','MarkerFaceAlpha',.8); colormap('jet'); colorbar;
hold on;
end


All_disN2=struct;

for j=1:length(clustersF)

idx_temp1=intersect(All_disN.(char(clustersF(j))).idx_f60_close,All_disN.(char(clustersF(j))).idx_s20_close);    
    
All_disN2.(char(clustersF(j))).idx=intersect(All_disN.(char(clustersF(j))).idx_s60_close,idx_temp1);

end

figure;
for j=1:length(clustersF)
    

scatter(ROI_temp2.f20(All_disN2.(char(clustersF(j))).idx,1),ROI_temp2.f20(All_disN2.(char(clustersF(j))).idx,2),10,'filled','MarkerFaceAlpha',.8); colormap('jet'); colorbar;
hold on;
end


figure;
for j=1:length(clustersF)
    

scatter3(ROI_temp2.f20(All_disN2.(char(clustersF(j))).idx,1),ROI_temp2.f20(All_disN2.(char(clustersF(j))).idx,2),ROI_temp2.f20(All_disN2.(char(clustersF(j))).idx,3),10,'filled','MarkerFaceAlpha',.8); colormap('jet'); colorbar;
hold on;
end



%% savig

save('distancesNCorr_groups2.mat','All_disN','All_disN2','F20vsF60_disN','F20vsS20_disN','F20vsS60_disN');


%% to graph

%%% to visualize it:


for j=1:length(clustersF)
    
    idx_clust=find(strcmp(clustersF, clustersF(j)));   
    temp_fieldname=fieldnames(f20_cleaned_idxs.clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(idx_clust);

    
    %%%% here i need to change the datasets to visualize
    X=F20_disNcorr.distance.(char(clustersF(j))).M.min';
    Y=F20_disNcorr.correlation.(char(clustersF(j))).Corr';
    numbins = 100;
    
  
    
    [values, centers] = hist3([X Y], [numbins numbins]);
    centers_X = centers{1,1};
    centers_Y = centers{1,2};
    binsize_X = abs(centers_X(2) - centers_X(1)) / 2;
    binsize_Y = abs(centers_Y(2) - centers_Y(1)) / 2;
    bins_X = zeros(numbins, 2);
    bins_Y = zeros(numbins, 2);
    for i = 1:numbins
        bins_X(i, 1) = centers_X(i) - binsize_X;
        bins_X(i, 2) = centers_X(i) + binsize_X;
        bins_Y(i, 1) = centers_Y(i) - binsize_Y;
        bins_Y(i, 2) = centers_Y(i) + binsize_Y;
    end
    scatter_COL = zeros(length(X), 1);
    onepercent = round(length(X) / 100);
    
    %fprintf('Generating colormap...\n');
    
    for i = 1:length(X)
        if (mod(i,onepercent) == 0)
            fprintf('.');
        end            
        last_lower_X = NaN;
        last_higher_X = NaN;
        id_X = NaN;
        c_X = X(i);
        last_lower_X = find(c_X >= bins_X(:,1));
        if (~isempty(last_lower_X))
            last_lower_X = last_lower_X(end);
        else
            last_higher_X = find(c_X <= bins_X(:,2));
            if (~isempty(last_higher_X))
                last_higher_X = last_higher_X(1);
            end
        end
        if (~isnan(last_lower_X))
            id_X = last_lower_X;
        else
            if (~isnan(last_higher_X))
                id_X = last_higher_X;
            end
        end
        last_lower_Y = NaN;
        last_higher_Y = NaN;
        id_Y = NaN;
        c_Y = Y(i);
        last_lower_Y = find(c_Y >= bins_Y(:,1));
        if (~isempty(last_lower_Y))
            last_lower_Y = last_lower_Y(end);
        else
            last_higher_Y = find(c_Y <= bins_Y(:,2));
            if (~isempty(last_higher_Y))
                last_higher_Y = last_higher_Y(1);
            end
        end
        if (~isnan(last_lower_Y))
            id_Y = last_lower_Y;
        else
            if (~isnan(last_higher_Y))
                id_Y = last_higher_Y;
            end
        end
        scatter_COL(i) = values(id_X, id_Y);
    
    end
    
%     fprintf(' Done!\n');
%     
%     fprintf('Plotting...');
%     
figure;
subplot(1,2,1)
    scatter(X, Y, 10, scatter_COL, 'filled'); colormap('jet'); colorbar;
subplot(1,2,2)  

scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),2),20,X,'filled','MarkerFaceAlpha',.8); colormap('jet'); colorbar;

% subplot(1,3,3)  
% scatter(ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),1),ROI_temp2.f20(f20_cleaned_idxs.clust_f20_CL4_cleaned.(char(temp_fieldname)),2),20,Y,'filled','MarkerFaceAlpha',.8); colormap('jet'); colorbar;
    
end