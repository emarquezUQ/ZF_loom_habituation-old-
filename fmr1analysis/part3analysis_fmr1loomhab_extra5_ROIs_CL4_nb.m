
%%%% this script is to check the number of ROIs of each functional cluster
%%%% in the different genotypes as it seems that fmr1 has less fasthab
%%%% ROIs. 

%%%from the wildtype loomhab dataset. to plot it if i need to
load('Nodes_N_means_alldatasets.mat','Zbrain_brainMask2D');


%%
load('s20_good_idx_Fish.mat','idx_Fish');
idx_Fish_cat=categorical(idx_Fish);

load('s20_fmr1_loomhab_CN_part2.mat','GoodBetas_ZS_CN_selected','idx_rsq');

load('s20_fmr1_loomhab_CN_part2_High_corr_Nb.mat','High_corr_Nb');

load('fmr1loomhab_BrainRegNclean.mat','PerBrainRegions','RegionList','ROI_temp2','idx_rsq_cleaned');

%% cleaning the clasification of the clusters. 
idx_clean=ismember(idx_rsq,idx_rsq_cleaned);
idx_clean=find(idx_clean);

High_corr_Nb=High_corr_Nb(idx_clean);

%%%% i also need to generate the list of fish saved in the
 %%%% fmr1loomhab_lists.m file
 
 load('s20_fmr1_loomhab_CN_part3.mat','idx_temp1','idx_temp2','idx_temp3','idx_temp4','idx_temp5');
 
 RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};



 load('fmr1_ROIs_with_WT_clusters.mat','idx_rsq_good','High_corr_Nb_s20','rawregress','fmr1_wt_clust_corr');
 
 
 %load('NodesNgraphFmr1Loomhab4.mat','Nodes4','goodorder_clust','Data_corrMat4','keep','discard');

 
 groupnames=fieldnames(Nodes4.ROIs_idx);
 
%%% this is what I found on the testing_making_modulesNgraph_fmr1loomhab4.m
%%% script after I got the coords for Unity. 

figure;
counter=1;
for g=[3 1 2]
     group=groupnames{g,1};
    
     if g==1
         temp_group_idx=idx_temp1;%idx_temp5;
     elseif g==2
         temp_group_idx=idx_temp2;
     else
         temp_group_idx=idx_temp4;
     end
         
fast_temp_idx1=[];
fast_temp_idx2=[];
for i=1:6

    subplot(3,6,counter);
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
    hold on;
    
    if i==1
        temp_idx1=find(High_corr_Nb_s20==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('fasthabSharp_',group,'_ROIs:',num2str(length(temp_coor))));  
        
    elseif i==2
        temp_idx1=find(High_corr_Nb_s20==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('fasthabMed_',group,'_ROIs:',num2str(length(temp_coor))));    
        
    elseif i==3
        temp_idx1=find(High_corr_Nb_s20==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('fasthabBroad_',group,'_ROIs:',num2str(length(temp_coor))));    
    
    elseif i==4
        temp_idx1=find(High_corr_Nb_s20==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('slopehab_',group,'_ROIs:',num2str(length(temp_coor))));        
    
    elseif i==5
        temp_idx1=find(High_corr_Nb_s20==i);        
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('nonhab_',group,'_ROIs:',num2str(length(temp_coor))));       
           
    elseif i==6
       temp_idx1=find(High_corr_Nb_s20==i);
        temp_idx2=intersect(temp_group_idx,idx_rsq_good(temp_idx1));
        temp_coor=ROI_temp2(temp_idx2,:);
        scatter(temp_coor(:,1),temp_coor(:,2),5,'filled');view(-90,90);
        title(strcat('inhib_',group,'_ROIs:',num2str(length(temp_coor))));              
          
    end
               
    counter=counter+1;
    
end       
        
end




%%

%%%% now, to get the number of ROIs per cluster (with fasthab subtypes and
%%%% without). i will make an average with the fish of each genotype. 



fish=vertcat(list1,list2,list3,list4);
list5=union(list1,list3);

%%% it seems that fish 201810048 from list1 dont have any ROIs... dont know
%%% why. might need to check if it is not a mistake in the name. i am
%%% getting rid of it on the list. 
find(idx_Fish==201810048);
fish(find(fish==201810048))=[];

%%%% Note. 'fish' list (41 fish) and unique(idx_Fish) (49) are different because the list has only the fish that got genotyped.  

Nb_ROIs_WT_clusters=struct;

for g=1:3    
    group=groupnames{g,1};
    
    if g==1 %hets
    idx_genotype=idx_temp5;%idx_temp5;  %%% temp5 is with all hets but I could also try with half to have similar numbers
    elseif g==2 % fmr1
     idx_genotype=idx_temp2;  
    elseif g==3 % controls
      idx_genotype=idx_temp4;  
    end
    
    for b=1:length(RegionList)
    
    idx_brain=PerBrainRegions.(RegionList{b}).idx;
        
    idx_temp=intersect(idx_brain,idx_genotype);
    
        for c=1:length(unique(High_corr_Nb_s20))
        
        idx_clust=idx_rsq_good(find(High_corr_Nb_s20==c));    
            
        idx_temp_2=intersect(idx_temp,idx_clust);
        
        if isempty(idx_temp_2)
            Nb_ROIs_WT_clusters.CL6.(group).(RegionList{b}){c}=[];
        else    
        temp_fish=unique(idx_Fish(idx_temp_2));  %%% I am not counting the fish that dont have ROIs. should I?
            for f=1:length(temp_fish) 
           
            idx_fish=intersect(idx_temp_2,find(idx_Fish==temp_fish(f)));
    
            Nb_ROIs_WT_clusters.CL6.(group).(RegionList{b}){c}(f,:)=length(idx_fish);
               
            end
        end
        end
                
        for c=[1 4 5 6]
        
         if c==1   
        idx_clust=idx_rsq_good(find(High_corr_Nb_s20==1 | High_corr_Nb_s20==2 | High_corr_Nb_s20==3 ));    
         else
          idx_clust=idx_rsq_good(find(High_corr_Nb_s20==c));   
         end
        
        idx_temp_2=intersect(idx_temp,idx_clust);
        
        if isempty(idx_temp_2)
            Nb_ROIs_WT_clusters.CL4.(group).(RegionList{b}){c}=[];
        else    
        temp_fish=unique(idx_Fish(idx_temp_2));  %%% I am not counting the fish that dont have ROIs. should I?
            for f=1:length(temp_fish) 
           
            idx_fish=intersect(idx_temp_2,find(idx_Fish==temp_fish(f)));
    
            Nb_ROIs_WT_clusters.CL4.(group).(RegionList{b}){c}(f,:)=length(idx_fish);
               
            end
        end
        end
                      
    end
 end



for g=1:3    
    group=groupnames{g,1};
    
    if g==1
      num=length(unique(idx_Fish(idx_temp5)));  %%% temp5 is with all hets but I could also try with half to have similar numbers
    elseif g==2
        num=length(unique(idx_Fish(idx_temp2)));
    else
        num=length(unique(idx_Fish(idx_temp4)));
    end
    
    for b=1:length(RegionList)
     
    for c=1:length(unique(High_corr_Nb_s20))
    %temp=[];
    temp=zeros(num,1);%%%% here i am taking into account the total number of fish per genotype to get the means and sd    
    temp(1:length(Nb_ROIs_WT_clusters.CL4.(group).(RegionList{b}){c}),1)=Nb_ROIs_WT_clusters.CL4.(group).(RegionList{b}){c};
    
    Nb_ROIs_WT_clusters.CL4_good.(group).(RegionList{b}){c}=temp;
        
    end   
    end   
end




for c=[1 4 5 6]
   figure;
   temp_mat={};
   for g=[3 2 1]    
    group=groupnames{g,1};
    temp=[];
   for b=1:length(RegionList)
       
   temp(1,b)=mean(Nb_ROIs_WT_clusters.CL4_good.(group).(RegionList{b}){c}); 
   temp_mat{g}(b,:)=Nb_ROIs_WT_clusters.CL4_good.(group).(RegionList{b}){c};
   end
   
   %plot(temp);
    
   %hold on; 
   end 
end

save('Nb_ROIs_WT_clusters.mat','Nb_ROIs_WT_clusters');


%% to check in the whole brain and not per brain region


fish=vertcat(list1,list2,list3,list4);
list5=union(list1,list3);

%%% it seems that fish 201810048 from list1 dont have any ROIs... dont know
%%% why. might need to check if it is not a mistake in the name. i am
%%% getting rid of it on the list. 
find(idx_Fish==201810048);
fish(find(fish==201810048))=[];


Nb_ROIs_WT_clusters_allbrain=struct;

for g=1:3    
    group=groupnames{g,1};
    
    if g==1 %hets
      num=length(unique(idx_Fish(idx_temp5)));  %%% temp5 is with all hets but I could also try with half to have similar numbers
      idx_genotype=idx_temp5;%idx_temp5;  %%% temp5 is with all hets but I could also try with half to have similar numbers
    
    elseif g==2 % fmr1
        num=length(unique(idx_Fish(idx_temp2)));
        idx_genotype=idx_temp2;  
    elseif g==3  % controls
        num=length(unique(idx_Fish(idx_temp4)));
         idx_genotype=idx_temp4;
    end
    
    
     
    for c=[1 4 5 6]
        
        idx_clust=idx_rsq_good(find(High_corr_Nb_s20==c));    
            
        idx_temp_2=intersect(idx_genotype,idx_clust);
    %temp=[];
    
    
    if isempty(idx_temp_2)
            Nb_ROIs_WT_clusters_allbrain.CL4_good.(group){c}=[];
        else    
        temp_fish=unique(idx_Fish(idx_temp_2));  %%% I am not counting the fish that dont have ROIs. should I?
        temp=zeros(num,1);%%%% here i am taking into account the total number of fish per genotype to get the means and sd    
               
        for f=1:length(temp_fish) 
           
            idx_fish=intersect(idx_temp_2,find(idx_Fish==temp_fish(f)));
    
             temp(f,1)=length(idx_fish);
                          
        end
            
        Nb_ROIs_WT_clusters_allbrain.CL4_good.(group){c}=temp;
        
        end
            
    end   
       
end

%%% for the ROIs that just passed the rsq threshold from any cluster and
%%% all brain ROIs (i did this before but with the 10 original clusters I
%%% think). 
for g=1:3    
    group=groupnames{g,1};
    
    if g==1 %hets
      num=length(unique(idx_Fish(idx_temp5)));  %%% temp5 is with all hets but I could also try with half to have similar numbers
      idx_genotype=idx_temp5;%idx_temp5;  %%% temp5 is with all hets but I could also try with half to have similar numbers
    
    elseif g==2 % fmr1
        num=length(unique(idx_Fish(idx_temp2)));
        idx_genotype=idx_temp2;  
    elseif g==3  % controls
        num=length(unique(idx_Fish(idx_temp4)));
         idx_genotype=idx_temp4;
    end
     
    
    idx_temp_2=intersect(idx_genotype,idx_rsq_good);
    
        temp_fish=unique(idx_Fish(idx_genotype));  %%% I am not counting the fish that dont have ROIs. should I?
        temp1=zeros(num,1);%%%% here i am taking into account the total number of fish per genotype to get the means and sd    
        temp2=zeros(num,1);%%%% here i am taking into account the total number of fish per genotype to get the means and sd    
                
        for f=1:length(temp_fish) 
           
            idx_fish1=intersect(idx_genotype,find(idx_Fish==temp_fish(f)));    
             temp1(f,1)=length(idx_fish1);
             
             idx_fish2=intersect(idx_temp_2,find(idx_Fish==temp_fish(f)));    
             temp2(f,1)=length(idx_fish2);
                          
        end
            
        Nb_ROIs_WT_clusters_allbrain.allROIs.(group){1}=temp1;
        Nb_ROIs_WT_clusters_allbrain.rsq_ROIs.(group){1}=temp2;
        
        end
            
 
