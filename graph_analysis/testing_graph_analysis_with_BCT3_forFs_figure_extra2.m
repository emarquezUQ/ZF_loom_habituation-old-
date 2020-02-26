%%% this script is to look at how the connections from the 19 hab-sensitive
%%% nodes from the previous script change after the stimuli presentation.
%%% Although I also did it for the 99 nodes.

%%% I had a few weird results that seem to be due to multiple divisions
%%% that ended up making many NaN and Inf values. the final result is
%%% interesting, specially for the blue to red and blue to blue connections in the 19 nodes
%%% which have a very strong recovery.However, the connections to from
%%% green to blue or green to red didnt increased dramatically. 


%%%%% I am not getting the same values comparisons between the same
%%%%% clusters. for example i dont get the same values for fasthab vs
%%%%% nonhab than nonhab vs fasthab. I think it is because I am 
%%%%% looking at the connections of the 19 to others but not only with in
%%%%% the 19. so there could be connections to other nodes from the same
%%%%% cluster that are not taking into account afterwords. 

%%% merging the green clusters
short_clust=Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60));
short_clust(find(short_clust==4))=2;
short_clust(find(short_clust==5))=2;

ratio_connect_19_short=struct;

for data=1
    
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));
    
    temp_deg1=MatAll_corrected.f20.(loom{1}).NodeConnClustDeg(hab_nodes_Par_f20_f60,:);
    
    %%%% merging fasthab and getting rid of inhib
    temp_deg1(:,1)=temp_deg1(:,1)+temp_deg1(:,2)+temp_deg1(:,3);
    temp_deg1(:,[2 3 6])=[];
    
   for k=1:length(loom)
       
    temp_deg2=MatAll_corrected.f20.(loom{k}).NodeConnClustDeg(hab_nodes_Par_f20_f60,:);
    
    %%%% merging fasthab and getting rid of inhib
    temp_deg2(:,1)=temp_deg2(:,1)+temp_deg2(:,2)+temp_deg2(:,3);
    temp_deg2(:,[2 3 6])=[];
       
       
     temp_ratio=temp_deg2./temp_deg1; %%%% is generating NaNs when it divides by 0 and Inf when 1/0
     
    for clust=[2 6 1]   
    temp_idx=find(short_clust==clust);
    
    ratio_connect_19_short.(loom{k}).(clustersF{clust}).ratios=temp_ratio(temp_idx,:);
    
    temp_ratio2=temp_ratio(temp_idx,:);
    
    temp_ratio2(find(temp_ratio2==Inf))=NaN;
    
    if size(temp_ratio2,1)==1
    
        ratio_connect_19_short.(loom{k}).(clustersF{clust}).means=temp_ratio2;
    else
    ratio_connect_19_short.(loom{k}).(clustersF{clust}).means=nanmean(temp_ratio2);  
    end
    
    end
    
   end
end


ratio_connect_19_short_all=[];
for data=1
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for clust=[2 6 1] 
    
 for  i=1:3 
     temp=[];
    for k=1:length(loom)
     temp(1,k)=ratio_connect_19_short.(loom{k}).(clustersF{clust}).means(i); 
    end 
    ratio_connect_19_short_all=vertcat(ratio_connect_19_short_all,temp);
 end 
 
end
end

%%

%%%% checking if my comment from the top is true by doing it for all the
%%%% nodes. no, i am not getting the same values and I am not sure why... 
%%%% most of the values are kind of similar but not the same. it could be because
%%%% of NaNs and Inf values... the problem is that for some of the looms
%%%% the values do differ significantly.

%%% so it seems that is because of some of the divisions... cause i did it
%%% just summing the degrees and then i did got the same values when
%%% comparing for example green vs blue and blue vs green. see below

%%% merging the green clusters
short_clust2=Nodes2.Mod_clust(keep);
short_clust2(find(short_clust2==4))=2;
short_clust2(find(short_clust2==5))=2;

ratio_connect_99_short2=struct;

for data=1
    
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));
    
    temp_deg1=MatAll_corrected.f20.(loom{1}).NodeConnClustDeg;
    
    %%%% merging fasthab and getting rid of inhib
    temp_deg1(:,1)=temp_deg1(:,1)+temp_deg1(:,2)+temp_deg1(:,3);
    temp_deg1(:,[2 3 6])=[];
    
   for k=1:length(loom)
       
    temp_deg2=MatAll_corrected.f20.(loom{k}).NodeConnClustDeg;
    
    %%%% merging fasthab and getting rid of inhib
    temp_deg2(:,1)=temp_deg2(:,1)+temp_deg2(:,2)+temp_deg2(:,3);
    temp_deg2(:,[2 3 6])=[];
       
       
     temp_ratio=temp_deg2./temp_deg1; %%%% is generating NaNs when it divides by 0 and Inf when 1/0
     
    for clust=[2 6 1]   
    temp_idx=find(short_clust2==clust);
    
    ratio_connect_99_short2.(loom{k}).(clustersF{clust}).ratios=temp_ratio(temp_idx,:);
    
    temp_ratio2=temp_ratio(temp_idx,:);
    
    temp_ratio2(find(temp_ratio2==Inf))=NaN;
    
    if size(temp_ratio2,1)==1
    
        ratio_connect_99_short2.(loom{k}).(clustersF{clust}).means=temp_ratio2;
    else
    ratio_connect_99_short2.(loom{k}).(clustersF{clust}).means=nanmean(temp_ratio2);  
    end
    
    end
    
   end
end


ratio_connect_99_short_all2=[];
for data=1
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for clust=[2 6 1] 
    
 for  i=1:3 
     temp=[];
    for k=1:length(loom)
     temp(1,k)=ratio_connect_99_short2.(loom{k}).(clustersF{clust}).means(i); 
    end 
    ratio_connect_99_short_all2=vertcat(ratio_connect_99_short_all2,temp);
 end 
 
end
end




%%


%% testing just summing the degrees. 

%%% here I do get the same values. so it most be in some of the divisions
%%% that I get the divergence. 

ratio_connect_99_short3=struct;

for data=1
    
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));
    
    temp_deg1=MatAll_corrected.f20.(loom{1}).NodeConnClustDeg;
    
    %%%% merging fasthab and getting rid of inhib
    temp_deg1(:,1)=temp_deg1(:,1)+temp_deg1(:,2)+temp_deg1(:,3);
    temp_deg1(:,[2 3 6])=[];
    
   for k=1:length(loom)
       
    temp_deg2=MatAll_corrected.f20.(loom{k}).NodeConnClustDeg;
    
    %%%% merging fasthab and getting rid of inhib
    temp_deg2(:,1)=temp_deg2(:,1)+temp_deg2(:,2)+temp_deg2(:,3);
    temp_deg2(:,[2 3 6])=[];
       
       
     test_deg=temp_deg2;%./temp_deg1; %%%% is generating NaNs when it divides by 0 and Inf when 1/0
     
%      temp_ratio(find(temp_ratio==Inf))=NaN;
%      temp_ratio(find(isnan(temp_ratio)))=0;
     
    for clust=[2 6 1]   
    temp_idx=find(short_clust2==clust);
    
    ratio_connect_99_short3.(loom{k}).(clustersF{clust}).ratios=test_deg(temp_idx,:);
    
    temp_ratio2=test_deg(temp_idx,:);
    
%     temp_ratio2(find(temp_ratio2==Inf))=NaN;
    
    if size(temp_ratio2,1)==1
    
        ratio_connect_99_short3.(loom{k}).(clustersF{clust}).means=temp_ratio2;
    else
    ratio_connect_99_short3.(loom{k}).(clustersF{clust}).means=nansum(temp_ratio2);  
    end
    
    end
    
   end
end


ratio_connect_99_short_all3=[];
for data=1
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for clust=[2 6 1] 
    
 for  i=1:3 
     temp=[];
    for k=1:length(loom)
     temp(1,k)=ratio_connect_99_short3.(loom{k}).(clustersF{clust}).means(i); 
    end 
    ratio_connect_99_short_all3=vertcat(ratio_connect_99_short_all3,temp);
 end 
 
end
end

%% so I will try doing the normalization at the end. 
ratio_connect_99_short_all4=[];
for i=1:size(ratio_connect_99_short_all3,1)
ratio_connect_99_short_all4(i,:)=ratio_connect_99_short_all3(i,:)/ratio_connect_99_short_all3(i,1);

end



%% so now back again with the main 19...



ratio_connect_19_short2=struct;

for data=1
    
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));
    
    temp_deg1=MatAll_corrected.f20.(loom{1}).NodeConnClustDeg(hab_nodes_Par_f20_f60,:);
    
    %%%% merging fasthab and getting rid of inhib
    temp_deg1(:,1)=temp_deg1(:,1)+temp_deg1(:,2)+temp_deg1(:,3);
    temp_deg1(:,[2 3 6])=[];
    
   for k=1:length(loom)
       
    temp_deg2=MatAll_corrected.f20.(loom{k}).NodeConnClustDeg(hab_nodes_Par_f20_f60,:);
    
    %%%% merging fasthab and getting rid of inhib
    temp_deg2(:,1)=temp_deg2(:,1)+temp_deg2(:,2)+temp_deg2(:,3);
    temp_deg2(:,[2 3 6])=[];
       
       
     test_deg=temp_deg2;%./temp_deg1; %%%% is generating NaNs when it divides by 0 and Inf when 1/0
     
%      temp_ratio(find(temp_ratio==Inf))=NaN;
%      temp_ratio(find(isnan(temp_ratio)))=0;
     
    for clust=[2 6 1]   
    temp_idx=find(short_clust==clust);
    
    ratio_connect_19_short2.(loom{k}).(clustersF{clust}).ratios=test_deg(temp_idx,:);
    
    temp_ratio2=test_deg(temp_idx,:);
    
%     temp_ratio2(find(temp_ratio2==Inf))=NaN;
    
    if size(temp_ratio2,1)==1
    
        ratio_connect_19_short2.(loom{k}).(clustersF{clust}).means=temp_ratio2;
    else
    ratio_connect_19_short2.(loom{k}).(clustersF{clust}).means=nansum(temp_ratio2);  
    end
    
    end
    
   end
end


ratio_connect_19_short_all2=[];
for data=1
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for clust=[2 6 1] 
    
 for  i=1:3 
     temp=[];
    for k=1:length(loom)
     temp(1,k)=ratio_connect_19_short2.(loom{k}).(clustersF{clust}).means(i); 
    end 
    ratio_connect_19_short_all2=vertcat(ratio_connect_19_short_all2,temp);
 end 
 
end
end

%% so I will try doing the normalization at the end. 
ratio_connect_19_short_all3=[];
for i=1:size(ratio_connect_19_short_all2,1)
ratio_connect_19_short_all3(i,:)=ratio_connect_19_short_all2(i,:)/ratio_connect_19_short_all2(i,1);

end


figure;
for i=1:size(ratio_connect_19_short_all3,1)
    title('19 nodes connections')
    plot(ratio_connect_19_short_all3(i,:));
    hold on;
end
