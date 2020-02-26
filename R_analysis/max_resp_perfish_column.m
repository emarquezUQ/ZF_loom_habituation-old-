

%%%%% geting the max values for each fish and each cluster. put them in a
%%%%% column per cluster so i can paste them into my excel template to then
%%%%% use GLMMs to analyse them. at the moment is only the tectum responses

%%%%% I need to mind the order of the clusters!!
%%%%% (nonhab,fasthab,slopehab,inhib) for the Fs
%%%%% (fasthab,nonhab,slopehab,inhib) for the Ss



load('All_means_maxResp_CL4_normalized.mat'); %%%% find in which folder it is

MaxColumn_fisnNclust_1stBlock={};
S_order=[2 1 3 4];
for clust=1:4
MaxColumn_fisnNclust_1stBlock{1,clust}=[];
for fish=1:size(Max_resp_f20_perfish{1,clust},1)
    temp=[];
    temp(:,1)=Max_resp_f20_perfish{1,clust}(fish,1:10);
    MaxColumn_fisnNclust_1stBlock{1,clust}=vertcat(MaxColumn_fisnNclust_1stBlock{1,clust},temp);
end

for fish=1:size(Max_resp_f60_perfish{1,clust},1)
    temp=[];
    temp(:,1)=Max_resp_f60_perfish{1,clust}(fish,1:10);
    MaxColumn_fisnNclust_1stBlock{1,clust}=vertcat(MaxColumn_fisnNclust_1stBlock{1,clust},temp);
end

for fish=1:size(Max_resp_s20_perfish{1,S_order(clust)},1)
    temp=[];
    temp(:,1)=Max_resp_s20_perfish{1,S_order(clust)}(fish,1:10);
    MaxColumn_fisnNclust_1stBlock{1,clust}=vertcat(MaxColumn_fisnNclust_1stBlock{1,clust},temp);
end

for fish=1:size(Max_resp_s60_perfish{1,S_order(clust)},1)
    temp=[];
    temp(:,1)=Max_resp_s60_perfish{1,S_order(clust)}(fish,1:10);
    MaxColumn_fisnNclust_1stBlock{1,clust}=vertcat(MaxColumn_fisnNclust_1stBlock{1,clust},temp);
end

end

%% second block


MaxColumn_fisnNclust_2ndBlock={};

for clust=1:4
MaxColumn_fisnNclust_2ndBlock{1,clust}=[];
for fish=1:size(Max_resp_f20_perfish{1,clust},1)
    temp=[];
    temp(:,1)=Max_resp_f20_perfish{1,clust}(fish,11:20);
    MaxColumn_fisnNclust_2ndBlock{1,clust}=vertcat(MaxColumn_fisnNclust_2ndBlock{1,clust},temp);
end

for fish=1:size(Max_resp_f60_perfish{1,clust},1)
    temp=[];
    temp(:,1)=Max_resp_f60_perfish{1,clust}(fish,11:20);
    MaxColumn_fisnNclust_2ndBlock{1,clust}=vertcat(MaxColumn_fisnNclust_2ndBlock{1,clust},temp);
end

for fish=1:size(Max_resp_s20_perfish{1,S_order(clust)},1)
    temp=[];
    temp(:,1)=Max_resp_s20_perfish{1,S_order(clust)}(fish,11:20);
    MaxColumn_fisnNclust_2ndBlock{1,clust}=vertcat(MaxColumn_fisnNclust_2ndBlock{1,clust},temp);
end

for fish=1:size(Max_resp_s60_perfish{1,S_order(clust)},1)
    temp=[];
    temp(:,1)=Max_resp_s60_perfish{1,S_order(clust)}(fish,11:20);
    MaxColumn_fisnNclust_2ndBlock{1,clust}=vertcat(MaxColumn_fisnNclust_2ndBlock{1,clust},temp);
end

end

