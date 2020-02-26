
%%% this script is to try to locate the the connections that recover from
%%% the 19 hab-sensitive nodes but based on which brain areas they are
%%% connecting. 

RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};

 
%%
%%%% this part is with the substraction of 10th minus 11th OF the f20 matrices
%%%% that passed the 0.75 threshold

big_subs_f20_10th_11th_habNodes

half_big_subs_f20_10th_11th_habNodes=triu(big_subs_f20_10th_11th_habNodes,1);

%%% or with a matrix based on the connections in both moments
%%% So I will find the non zero values in both matrices (10th and 11th) and
%%% get their values from the non-thresholded matrices to do the substraction. 
raw_10th=abs(Data_corrMat2.f20.Mean_corrMat{1,11}(keep,keep));
raw_11th=abs(Data_corrMat2.f20.Mean_corrMat{1,12}(keep,keep));

new_10th=MatAll_corrected.f20.loom10.Mat;
new_11th=MatAll_corrected.f20.loom11.Mat;

new_10th(find(MatAll_corrected.f20.loom11.Mat))=raw_10th(find(MatAll_corrected.f20.loom11.Mat));

new_11th(find(MatAll_corrected.f20.loom10.Mat))=raw_11th(find(MatAll_corrected.f20.loom10.Mat));

test_mat=new_10th-new_11th;

big_subs_f20_10th_11th_habNodes2=zeros(size(MatAll_corrected.f20.loom11.Mat));

big_subs_f20_10th_11th_habNodes2(:,hab_nodes_Par_f20_f60)=test_mat(:,hab_nodes_Par_f20_f60);

big_subs_f20_10th_11th_habNodes2(hab_nodes_Par_f20_f60,:)=test_mat(hab_nodes_Par_f20_f60,:);

half_big_subs_f20_10th_11th_habNodes=triu(big_subs_f20_10th_11th_habNodes2,1);

%%

figure;hist(half_big_subs_f20_10th_11th_habNodes);
figure;imagesc(half_big_subs_f20_10th_11th_habNodes);


%%% how many and which brain areas are included in the 19 nodes? 7...  this
%%% means I will have a lot of possible pairs... although some might not be
%%% that intereting... 
unique(Nodes2.Mod_brain(keep(hab_nodes_Par_f20_f60)));
length(unique(Nodes2.Mod_brain(keep(hab_nodes_Par_f20_f60))));


%%% to how many brain areas they connect? to 8... all of them except
%%% subpallium. that will be a lot... but maybe we can select the top ones
%%% later... 
unique(Nodes2.Mod_brain(keep(find(sum(half_big_subs_f20_10th_11th_habNodes)))))



%%%% so what i will try to do is generate a cell matrix with as rows the
%%%% brain areas to which the 19 nodes belong and columns the brain areas
%%%% they connect to. inside each cell I will have a histocount

edges=[-1:0.095:1];

subs_brain_conn_count={};
for brain=unique(Nodes2.Mod_brain(keep(hab_nodes_Par_f20_f60)))'
    temp_brain=find(Nodes2.Mod_brain(keep)==brain);
    temp_brain_mat=half_big_subs_f20_10th_11th_habNodes(intersect(temp_brain,hab_nodes_Par_f20_f60),:);
    
    if size(temp_brain_mat,1)==1
        temp_brain_vec=temp_brain_mat;
        
    else 
        temp_brain_vec=sum(temp_brain_mat);
    end
    
    for brain2=unique(Nodes2.Mod_brain(keep(find(temp_brain_vec))))'
       temp_brain2=find(Nodes2.Mod_brain(keep)==brain2);
    [subs_brain_conn_count{brain,brain2},~]=histcounts(nonzeros(temp_brain_mat(:,temp_brain2)),edges);
    end
    
end


%%% stacking them on top of each other. As expected most of them are
%%% negative... as I selected nodes thare were more sensitive to recovery. 
subs_brain_conn_count_mat=[];
brain_id_conn_count=[];
for i=1:9
    for j=1:9
        
    if isempty(subs_brain_conn_count{i,j})
        continue
    else
        temp=subs_brain_conn_count{i,j};
        subs_brain_conn_count_mat=vertcat(subs_brain_conn_count_mat,temp);
        temp2=[i j];
        brain_id_conn_count=vertcat(brain_id_conn_count,temp2);
    end
    end
end

brain_conn_labels={};
for i=1:length(brain_id_conn_count)
    temp=brain_id_conn_count(i,:);
brain_conn_labels{i}=strcat(RegionList(temp(1)),'-',RegionList(temp(2)));
end

%%% making figures
figure;imagesc(subs_brain_conn_count_mat);colormap(inferno);colorbar;%caxis([0 8]);
xticks([1 11 21])
xticklabels({'-1','0','1'})
yticks([1:39]);
yticklabels(string(brain_conn_labels));

figure;plot(sum(subs_brain_conn_count_mat));


figure;spy(subs_brain_conn_count_mat);
xticks([1 11 21])
xticklabels({'-1','0','1'})
yticks([1:length(brain_id_conn_count)]);
yticklabels(string(brain_conn_labels));


%%
%%%% this part is with the substraction of 10th minus 11th OF the matrices
%%%% that had a substraction already between f20-f60!!!
big_subs_f20_f60_habNodes


figure;hist(triu(big_subs_f20_f60_habNodes,1));
figure;imagesc(triu(big_subs_f20_f60_habNodes,1));

half_big_subs_f20_f60_habNodes=triu(big_subs_f20_f60_habNodes,1);

%%% how many and which brain areas are included in the 19 nodes? 7...  this
%%% means I will have a lot of possible pairs... although some might not be
%%% that intereting... 
unique(Nodes2.Mod_brain(keep(hab_nodes_Par_f20_f60)));
length(unique(Nodes2.Mod_brain(keep(hab_nodes_Par_f20_f60))));


%%% to how many brain areas they connect? to 8... all of them except
%%% subpallium. that will be a lot... but maybe we can select the top ones
%%% later... 
unique(Nodes2.Mod_brain(keep(find(sum(half_big_subs_f20_f60_habNodes)))))



%%%% so what i will try to do is generate a cell matrix with as rows the
%%%% brain areas to which the 19 nodes belong and columns the brain areas
%%%% they connect to. inside each cell I will have a histocount

edges=[-1:0.095:1];

subs_brain_conn_count={};
for brain=unique(Nodes2.Mod_brain(keep(hab_nodes_Par_f20_f60)))'
    temp_brain=find(Nodes2.Mod_brain(keep)==brain);
    temp_brain_mat=half_big_subs_f20_f60_habNodes(intersect(temp_brain,hab_nodes_Par_f20_f60),:);
    
    if size(temp_brain_mat,1)==1
        temp_brain_vec=temp_brain_mat;
        
    else 
        temp_brain_vec=sum(temp_brain_mat);
    end
    
    for brain2=unique(Nodes2.Mod_brain(keep(find(temp_brain_vec))))'
       temp_brain2=find(Nodes2.Mod_brain(keep)==brain2);
    [subs_brain_conn_count{brain,brain2},~]=histcounts(nonzeros(temp_brain_mat(:,temp_brain2)),edges);
    end
    
end


%%% stacking them on top of each other. As expected most of them are
%%% negative... as I selected nodes thare were more sensitive to recovery. 
subs_brain_conn_count_mat=[];
brain_id_conn_count=[];
for i=1:9
    for j=1:9
        
    if isempty(subs_brain_conn_count{i,j})
        continue
    else
        temp=subs_brain_conn_count{i,j};
        subs_brain_conn_count_mat=vertcat(subs_brain_conn_count_mat,temp);
        temp2=[i j];
        brain_id_conn_count=vertcat(brain_id_conn_count,temp2);
    end
    end
end

brain_conn_labels={};
for i=1:39
    temp=brain_id_conn_count(i,:);
brain_conn_labels{i}=strcat(RegionList(temp(1)),'-',RegionList(temp(2)));
end

%%% making figures
figure;imagesc(subs_brain_conn_count_mat);colormap(inferno);colorbar;caxis([0 8]);
xticks([1 11 21])
xticklabels({'-1','0','1'})
yticks([1:39]);
yticklabels(string(brain_conn_labels));

figure;plot(sum(subs_brain_conn_count_mat));




conserved_conn=[];
for i=1:size(subs_brain_conn_count_mat,1)
if sum(subs_brain_conn_count_mat(i,11:end))>0

    conserved_conn=vertcat(conserved_conn,i);
    
else  
end

end


%%% there are brain connections that have both positive and negative values
%%% wich ones didnt change or are positive?
brain_id_conn_count(conserved_conn,:)


%%% wich ones recovered? 
connections=1:39;
recovered_conn = setdiff(connections,conserved_conn);

brain_id_conn_count(recovered_conn,:)


%% with the 99 nodes
%%% so what if I use the substraction matrix of the 99 nodes of f20 from the 10th to the 11th loom. 


figure;hist(triu(MatAll_corrected.f20.loom10.Mat-MatAll_corrected.f20.loom11.Mat,1));
figure;imagesc(triu(MatAll_corrected.f20.loom10.Mat-MatAll_corrected.f20.loom11.Mat,1));

%%%% this one gives too clear cut because of the 0.75 threshold. 
%half_subs_f20_10th_11th=triu(MatAll_corrected.f20.loom10.Mat-MatAll_corrected.f20.loom11.Mat,1);

%%% So I will find the non zero values in both matrices (10th and 11th) and
%%% get their values from the non-thresholded matrices to do the substraction. 
raw_10th=abs(Data_corrMat2.f20.Mean_corrMat{1,11}(keep,keep));
raw_11th=abs(Data_corrMat2.f20.Mean_corrMat{1,12}(keep,keep));

new_10th=MatAll_corrected.f20.loom10.Mat;
new_11th=MatAll_corrected.f20.loom11.Mat;

new_10th(find(MatAll_corrected.f20.loom11.Mat))=raw_10th(find(MatAll_corrected.f20.loom11.Mat));

new_11th(find(MatAll_corrected.f20.loom10.Mat))=raw_11th(find(MatAll_corrected.f20.loom10.Mat));

%%% now the substraction
half_subs_f20_10th_11th=triu(new_10th-new_11th,1);

%%% to how many brain areas they connect? to 8... all of them except
%%% subpallium. that will be a lot... but maybe we can select the top ones
%%% later... 
unique(Nodes2.Mod_brain(keep(find(sum(half_subs_f20_10th_11th)))))



%%%% so what i will try to do is generate a cell matrix with as rows the
%%%% brain areas to which the 99 nodes belong and columns the brain areas
%%%% they connect to. inside each cell I will have a histocount

edges=[-1:0.095:1];

subs_brain_conn_count_all={};
for brain=unique(Nodes2.Mod_brain(keep))'
    temp_brain=find(Nodes2.Mod_brain(keep)==brain);
    temp_brain_mat=half_subs_f20_10th_11th(temp_brain,:);
    
    if size(temp_brain_mat,1)==1
        temp_brain_vec=temp_brain_mat;
        
    else 
        temp_brain_vec=sum(temp_brain_mat);
    end
    
    for brain2=unique(Nodes2.Mod_brain(keep(find(temp_brain_vec))))'
       temp_brain2=find(Nodes2.Mod_brain(keep)==brain2);
    [subs_brain_conn_count_all{brain,brain2},~]=histcounts(nonzeros(temp_brain_mat(:,temp_brain2)),edges);
    
    end
    
end

%%% to merge the same combinations
subs_brain_conn_count_all_cleaned={};
subs_brain_conn_count_all_cleaned{9,9}=[];
for i=1:length(subs_brain_conn_count_all)

    for j=1:length(subs_brain_conn_count_all)
       
        if i==j
            subs_brain_conn_count_all_cleaned{i,j}=subs_brain_conn_count_all{i,j};
         elseif ~isempty(subs_brain_conn_count_all{i,j}) & isempty(subs_brain_conn_count_all{j,i})
             subs_brain_conn_count_all_cleaned{i,j}=subs_brain_conn_count_all{i,j};
         elseif ~isempty(subs_brain_conn_count_all{i,j}) & ~isempty(subs_brain_conn_count_all{j,i})
             if isempty(subs_brain_conn_count_all_cleaned{i,j}) & isempty(subs_brain_conn_count_all_cleaned{j,i})
             subs_brain_conn_count_all_cleaned{i,j}=subs_brain_conn_count_all{i,j}+subs_brain_conn_count_all{j,i};
             else
             subs_brain_conn_count_all_cleaned{i,j}=[];  
             end
         else
            subs_brain_conn_count_all_cleaned{i,j}=[];
        end
    
    end
end


%%% stacking them on top of each other. As expected most of them are
%%% negative... as I selected nodes thare were more sensitive to recovery. 
subs_brain_conn_count_mat2=[];
brain_id_conn_count_all=[];
for i=1:9
    for j=1:9
        
    if isempty(subs_brain_conn_count_all_cleaned{i,j})
        continue
    else
        temp=subs_brain_conn_count_all_cleaned{i,j};
        subs_brain_conn_count_mat2=vertcat(subs_brain_conn_count_mat2,temp);
        temp2=[i j];
        brain_id_conn_count_all=vertcat(brain_id_conn_count_all,temp2);
    end
    end
end


brain_conn_labels_all={};
for i=1:size(brain_id_conn_count_all,1)
    temp=brain_id_conn_count_all(i,:);
brain_conn_labels_all{i}=strcat(RegionList(temp(1)),'-',RegionList(temp(2)));
end


figure;imagesc(subs_brain_conn_count_mat2);
xticks([1 11 21])
xticklabels({'-1','0','1'})
yticks([1:size(brain_id_conn_count_all,1)]);
yticklabels(string(brain_conn_labels_all));


figure;plot(sum(subs_brain_conn_count_mat2));

%% now with the functional clusters



%%% this script is to try to locate the the connections that recover from
%%% the 19 hab-sensitive nodes but based on which functional clusters they are
%%% connecting. 

%%% I will be merging the fast-hab clusters

short_clust=Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60));
short_clust(find(short_clust==4))=2;
short_clust(find(short_clust==5))=2;

short_clust2=Nodes2.Mod_clust(keep);
short_clust2(find(short_clust2==4))=2;
short_clust2(find(short_clust2==5))=2;


figure;hist(triu(big_subs_f20_f60_habNodes,1));
figure;imagesc(triu(big_subs_f20_f60_habNodes,1));

half_big_subs_f20_f60_habNodes=triu(big_subs_f20_f60_habNodes,1);

%%% how many and which clust areas are included in the 19 nodes? 7...  this
%%% means I will have a lot of possible pairs... although some might not be
%%% that intereting... 



unique(short_clust);
length(unique(short_clust));


%%% to how many clust areas they connect? to 4... all of them. 
unique(short_clust2(find(sum(half_big_subs_f20_f60_habNodes))))



%%%% so what i will try to do is generate a cell matrix with as rows the
%%%% clust areas to which the 19 nodes belong and columns the clust areas
%%%% they connect to. inside each cell I will have a histocount

edges=[-1:0.095:1];

subs_clust_conn_count={};
for clust=unique(short_clust)'
    temp_clust=find(short_clust2==clust);
    temp_clust_mat=half_big_subs_f20_f60_habNodes(intersect(temp_clust,hab_nodes_Par_f20_f60),:);
    
    if size(temp_clust_mat,1)==1
        temp_clust_vec=temp_clust_mat;
        
    else 
        temp_clust_vec=sum(temp_clust_mat);
    end
    
    for clust2=unique(short_clust2(find(temp_clust_vec)))'
       temp_clust2=find(short_clust2==clust2);
    [subs_clust_conn_count{clust,clust2},~]=histcounts(nonzeros(temp_clust_mat(:,temp_clust2)),edges);
    end
    
end


%%% stacking them on top of each other. 
subs_clust_conn_count_mat=[];
clust_id_conn_count=[];
for i=1:size(subs_clust_conn_count,1)
    for j=1:size(subs_clust_conn_count,2)
        
    if isempty(subs_clust_conn_count{i,j})
        continue
    else
        temp=subs_clust_conn_count{i,j};
        subs_clust_conn_count_mat=vertcat(subs_clust_conn_count_mat,temp);
        temp2=[i j];
        clust_id_conn_count=vertcat(clust_id_conn_count,temp2);
    end
    end
end

clust_conn_labels={};
for i=1:length(clust_id_conn_count)
    temp=clust_id_conn_count(i,:);
clust_conn_labels{i}=strcat(clustersF(temp(1)),'-',clustersF(temp(2)));
end

%%% making figures
figure;imagesc(subs_clust_conn_count_mat);colormap(inferno);colorbar;
xticks([1 11 21])
xticklabels({'-1','0','1'})
yticks([1:39]);
yticklabels(string(clust_conn_labels));

figure;plot(sum(subs_clust_conn_count_mat));




conserved_conn=[];
for i=1:size(subs_clust_conn_count_mat,1)
if sum(subs_clust_conn_count_mat(i,11:end))>0

    conserved_conn=vertcat(conserved_conn,i);
    
else  
end

end


%%% there are clust connections that have both positive and negative values
%%% wich ones didnt change or are positive?
clust_id_conn_count(conserved_conn,:)


%%% wich ones recovered? 
connections=1:7;
recovered_conn = setdiff(connections,conserved_conn);

clust_id_conn_count(recovered_conn,:)


%%
%%% so what if I use the substraction matrix of the 11th loom of f20 minus
%%% f60? is seems harderd to interpret. 


figure;hist(triu(MatAll_corrected.f20.loom10.Mat-MatAll_corrected.f20.loom11.Mat,1));
figure;imagesc(triu(MatAll_corrected.f20.loom10.Mat-MatAll_corrected.f20.loom11.Mat,1));

%%%% this one gives too clear cut because of the 0.75 threshold. 
%half_subs_f20_10th_11th=triu(MatAll_corrected.f20.loom10.Mat-MatAll_corrected.f20.loom11.Mat,1);

%%% So I will find the non zero values in both matrices (10th and 11th) and
%%% get their values from the non-thresholded matrices to do the substraction. 
raw_10th=abs(Data_corrMat2.f20.Mean_corrMat{1,11}(keep,keep));
raw_11th=abs(Data_corrMat2.f20.Mean_corrMat{1,12}(keep,keep));

new_10th=MatAll_corrected.f20.loom10.Mat;
new_11th=MatAll_corrected.f20.loom11.Mat;

new_10th(find(MatAll_corrected.f20.loom11.Mat))=raw_10th(find(MatAll_corrected.f20.loom11.Mat));

new_11th(find(MatAll_corrected.f20.loom10.Mat))=raw_11th(find(MatAll_corrected.f20.loom10.Mat));

%%% now the substraction
half_subs_f20_10th_11th=triu(new_10th-new_11th,1);

%%% to how many clust areas they connect? to 8... all of them except
%%% subpallium. that will be a lot... but maybe we can select the top ones
%%% later... 
unique(short_clust2(find(sum(half_subs_f20_10th_11th))))



%%%% so what i will try to do is generate a cell matrix with as rows the
%%%% clust areas to which the 19 nodes belong and columns the clust areas
%%%% they connect to. inside each cell I will have a histocount

edges=[-1:0.095:1];

subs_clust_conn_count_all={};
for clust=unique(short_clust2)'
    temp_clust=find(short_clust2==clust);
    temp_clust_mat=half_subs_f20_10th_11th(temp_clust,:);
    
    if size(temp_clust_mat,1)==1
        temp_clust_vec=temp_clust_mat;
        
    else 
        temp_clust_vec=sum(temp_clust_mat);
    end
    
    for clust2=unique(short_clust2(find(temp_clust_vec)))'
       temp_clust2=find(short_clust2==clust2);
    [subs_clust_conn_count_all{clust,clust2},~]=histcounts(nonzeros(temp_clust_mat(:,temp_clust2)),edges);
    
    end
    
end

%%% to merge the same combinations
% subs_clust_conn_count_all_cleaned{9,9}=[];
% for i=1:size(subs_clust_conn_count_all,1)
% 
%     for j=1:size(subs_clust_conn_count_all,2)
%        
%         if i==j
%             subs_clust_conn_count_all_cleaned{i,j}=subs_clust_conn_count_all{i,j};
%          elseif ~isempty(subs_clust_conn_count_all{i,j}) & isempty(subs_clust_conn_count_all{j,i})
%              subs_clust_conn_count_all_cleaned{i,j}=subs_clust_conn_count_all{i,j};
%          elseif ~isempty(subs_clust_conn_count_all{i,j}) & ~isempty(subs_clust_conn_count_all{j,i})
%              if isempty(subs_clust_conn_count_all_cleaned{i,j}) & isempty(subs_clust_conn_count_all_cleaned{j,i})
%              subs_clust_conn_count_all_cleaned{i,j}=subs_clust_conn_count_all{i,j}+subs_clust_conn_count_all{j,i};
%              else
%              subs_clust_conn_count_all_cleaned{i,j}=[];  
%              end
%          else
%             subs_clust_conn_count_all_cleaned{i,j}=[];
%         end
%     
%     end
% end


%%% stacking them on top of each other. As expected most of them are
%%% negative... as I selected nodes thare were more sensitive to recovery. 
subs_clust_conn_count_mat2=[];
clust_id_conn_count_all=[];
for i=1:6
    for j=1:7
        
    if isempty(subs_clust_conn_count_all{i,j})
        continue
    else
        temp=subs_clust_conn_count_all{i,j};
        subs_clust_conn_count_mat2=vertcat(subs_clust_conn_count_mat2,temp);
        temp2=[i j];
        clust_id_conn_count_all=vertcat(clust_id_conn_count_all,temp2);
    end
    end
end


clust_conn_labels_all={};
for i=1:size(clust_id_conn_count_all,1)
    temp=clust_id_conn_count_all(i,:);
clust_conn_labels_all{i}=strcat(clustersF(temp(1)),'-',clustersF(temp(2)));
end


figure;imagesc(subs_clust_conn_count_mat2);
xticks([1 11 21])
xticklabels({'-1','0','1'})
yticks([1:size(clust_id_conn_count_all,1)]);
yticklabels(string(clust_conn_labels_all));


figure;plot(sum(subs_clust_conn_count_mat2));
