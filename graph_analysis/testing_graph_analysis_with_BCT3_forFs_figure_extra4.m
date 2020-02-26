
%%%% this script is to try to represent the difference between the 10th and 11th loom connections based on the brain
%%%% areas and functional cluster. 


%% I will try first with the 19 nodes

%%% brain regions
RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};

%%% merging the fasthab clusters
short_clust2=Nodes2.Mod_clust(keep);
short_clust2(find(short_clust2==4))=2;
short_clust2(find(short_clust2==5))=2;




%%%% this part is with the substraction of 10th minus 11th OF the f20 matrices
%%%% that passed the 0.75 threshold

%half_big_subs_f20_10th_11th_habNodes=triu(big_subs_f20_10th_11th_habNodes,1);

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


%%%% so what i will try is to get all the connections, their value, which
%%%% parts of the brain they are connecting and which functional clusters
%%%% the nodes belong to and make a matrix with all this information. 


subs_10th_11th_connection_data=[];
% brain_id_conn_count_all=[];
for i=1:99
    temp=half_big_subs_f20_10th_11th_habNodes(i,:);
    temp2=find(temp);
    if ~isempty(temp2)
        for j=temp2
            
            subs=half_big_subs_f20_10th_11th_habNodes(i,j);
            brain1=Nodes2.Mod_brain(keep(i));
            brain2=Nodes2.Mod_brain(keep(j));
            clust1=short_clust2(i);
            clust2=short_clust2(j);
            
            brain=sort([brain1 brain2]);
            brain=num2str(brain);
            brain= brain(~isspace(brain));
            brain= str2num(brain);
            
            y=0;
            
            clust=sort([clust1 clust2]);
            clust=num2str(clust);
            clust= clust(~isspace(clust));
            clust= str2num(clust);    
                                  
            temp3=[subs y brain clust];
            
            subs_10th_11th_connection_data=vertcat(subs_10th_11th_connection_data,temp3);
            
%             temp_brain=[brain1 brain2];
%             brain_id_conn_count_all=vertcat(brain_id_conn_count_all,temp_brain);
            
        end       
    else
    end
end


brain_conn=unique(subs_10th_11th_connection_data(:,3))';
counter=1;
for i=brain_conn
    
    temp=find(subs_10th_11th_connection_data(:,3)==i);
    
    subs_10th_11th_connection_data(temp,2)=counter;
    
    counter=counter+1;
end

[~,temp_idx] = sort(subs_10th_11th_connection_data(:,3)); % sort just the first column
subs_10th_11th_connection_data = subs_10th_11th_connection_data(temp_idx,:);   % sort the whole matrix using the sort indices

clust_conn=unique(subs_10th_11th_connection_data(:,4))';


brain_conn_labels_all={};
for i=1:size(brain_conn,2)
    temp=brain_conn(i);
    temp=num2str(temp);
    
brain_conn_labels_all{i}=strcat(RegionList(str2num(temp(1))),'-',RegionList(str2num(temp(2))));
end



%%% to plot it
% figure;
% for i=clust_conn
%     temp=find(subs_10th_11th_connection_data(:,4)==i);
%     
%     scatter(subs_10th_11th_connection_data(temp,1),subs_10th_11th_connection_data(temp,2),25,'filled');xlim([-0.8 0.2]);
%     yticks([1:size(brain_conn_labels_all,2)]);
%     yticklabels(string(brain_conn_labels_all));
%     
%     hold on;
% end

%%% to plot it better
figure;
gscatter(subs_10th_11th_connection_data(:,1),subs_10th_11th_connection_data(:,2),subs_10th_11th_connection_data(:,4),'rymkgcb','.',25);xlim([-0.8 0.2]);
yticks([1:size(brain_conn_labels_all,2)]);
yticklabels(string(brain_conn_labels_all));

%% now for the 99 nodes

%%% looks very good too. we might stick with this one as there is more
%%% info. 

half_subs_f20_10th_11th


subs_10th_11th_connection_data2=[];
% brain_id_conn_count_all=[];
for i=1:99
    temp=half_subs_f20_10th_11th(i,:);
    temp2=find(temp);
    if ~isempty(temp2)
        for j=temp2
            
            subs=half_subs_f20_10th_11th(i,j);
            brain1=Nodes2.Mod_brain(keep(i));
            brain2=Nodes2.Mod_brain(keep(j));
            clust1=short_clust2(i);
            clust2=short_clust2(j);
            
            brain=sort([brain1 brain2]);
            brain=num2str(brain);
            brain= brain(~isspace(brain));
            brain= str2num(brain);
            
            y=0;
            
            clust=sort([clust1 clust2]);
            clust=num2str(clust);
            clust= clust(~isspace(clust));
            clust= str2num(clust);    
                                  
            temp3=[subs y brain clust];
            
            subs_10th_11th_connection_data2=vertcat(subs_10th_11th_connection_data2,temp3);
            
%             temp_brain=[brain1 brain2];
%             brain_id_conn_count_all=vertcat(brain_id_conn_count_all,temp_brain);
            
        end       
    else
    end
end


brain_conn2=unique(subs_10th_11th_connection_data2(:,3))';
brain_conn2=sort(brain_conn2,'descend');
counter=1;
for i=brain_conn2
    
    temp=find(subs_10th_11th_connection_data2(:,3)==i);
    
    subs_10th_11th_connection_data2(temp,2)=counter;
    
    counter=counter+1;
end

[~,temp_idx] = sort(subs_10th_11th_connection_data2(:,3)); % sort just the first column
subs_10th_11th_connection_data2 = subs_10th_11th_connection_data2(temp_idx,:);   % sort the whole matrix using the sort indices

clust_conn2=unique(subs_10th_11th_connection_data2(:,4))';


brain_conn_labels_all2={};
for i=1:size(brain_conn2,2)
    temp=brain_conn2(i);
    temp=num2str(temp);
    
brain_conn_labels_all2{i}=strcat(RegionList(str2num(temp(1))),'-',RegionList(str2num(temp(2))));
end



%%% to plot it better
figure;
gscatter(subs_10th_11th_connection_data2(:,1),subs_10th_11th_connection_data2(:,2),subs_10th_11th_connection_data2(:,4),'rymkgcb','.',25);xlim([-0.8 0.2]);
yticks([1:size(brain_conn_labels_all2,2)]);
yticklabels(string(brain_conn_labels_all2));


%%% to get the data for prism
subs_10th_11th_connection_data2_cell={};
for i=1:length(clust_conn2)
    
temp=find(subs_10th_11th_connection_data2(:,4)==clust_conn2(i));
subs_10th_11th_connection_data2_cell{i}=subs_10th_11th_connection_data2(temp,:);

end

%% now for the 99 nodes but with less brain regions

%%% only keeping Pallium Thalamus Tectum Tegmentum Hindbrain

%brain_to_discard=[4 5 8];

conn_to_discard=[];
for i=1:length(subs_10th_11th_connection_data2)
   temp=num2str(subs_10th_11th_connection_data2(i,3)); 
    
    if str2num(temp(1))==4 | str2num(temp(1))==5 | str2num(temp(1))==8 | str2num(temp(2))==4 | str2num(temp(2))==5 | str2num(temp(2))==8
    
        conn_to_discard(i)=1;
    else
        conn_to_discard(i)=0;
    end
    
end

idx_conn_discard=find(conn_to_discard);

subs_10th_11th_connection_data2_short=subs_10th_11th_connection_data2;
subs_10th_11th_connection_data2_short(idx_conn_discard,:)=[];

brain_conn2_short=unique(subs_10th_11th_connection_data2_short(:,3))';
brain_conn2_short=sort(brain_conn2_short,'descend');
counter=1;
for i=brain_conn2_short
    
    temp=find(subs_10th_11th_connection_data2_short(:,3)==i);
    
    subs_10th_11th_connection_data2_short(temp,2)=counter;
    
    counter=counter+1;
end


brain_conn_labels_all2_short={};
for i=1:size(brain_conn2_short,2)
    temp=brain_conn2_short(i);
    temp=num2str(temp);
    
brain_conn_labels_all2_short{i}=strcat(RegionList(str2num(temp(1))),'-',RegionList(str2num(temp(2))));
end



%%% to plot it better
figure;
gscatter(subs_10th_11th_connection_data2_short(:,1),subs_10th_11th_connection_data2_short(:,2),subs_10th_11th_connection_data2_short(:,4),'rymkgcb','.',25);xlim([-0.8 0.6]);
yticks([1:size(brain_conn_labels_all2_short,2)]);
yticklabels(string(brain_conn_labels_all2_short));

%%% to get the data for prism
subs_10th_11th_connection_data2_short_cell={};
for i=1:length(clust_conn2)
    
temp=find(subs_10th_11th_connection_data2_short(:,4)==clust_conn2(i));
subs_10th_11th_connection_data2_short_cell{i}=subs_10th_11th_connection_data2_short(temp,:);

end