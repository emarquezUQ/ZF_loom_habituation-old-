


%%%%% this script is to find the possible boundaries of gamma and omega in the fmr1 dataset. 

%%%% the maximum gamma could be when having a number of communities that is
%%%% 2/3 of the nodes = 60, half the number of nodes = 45  or a quarter = 23. the minim gamma
%%%% could be when there are less than 4 comunities as this are my
%%%% functional clusters

%%%% the maximum omega could be when there is no change across the time
%%%% dimension in at least one node. Or when a minimum portion of the nodes
%%%% fail to detect (change community) at the first and 11th looms. the min
%%%% omega could be when there is a change at every time point or just 0.1 value



load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

keepFmr1=load('graphs_fmr1Exp.mat','keep');
keepFmr1=keepFmr1.keep;

%%


%%%% just checking how many nodes I have per brain region and cluster

%%% the clustID tags are not in the right order... I will change it
oldclusttag=[4 2 5 6 1 7];
ClustID=zeros(size(Nodes.Nod_clustID(keep)));
for i=1:length(unique(Nodes.Nod_clustID(keep)))
    
    idx_temp=find(Nodes.Nod_clustID(keep)==oldclusttag(i));
    ClustID(idx_temp)=i;
   
end
 
clustnames={'fasthab1','fasthab2','fasthab3','modhab','weakhab','inhib'};

% figure;
% subplot(1,2,1);histogram(Nodes.Nod_brainID(keep)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);
% subplot(1,2,2);histogram(ClustID); xticks([1:6]); xticklabels(clustnames);xtickangle(45);


%%%% with fasthab merged
oldclusttag=[4 2 5 6 1 7];
ClustID_CL4=zeros(size(Nodes.Nod_clustID(keep)));
for i=1:length(unique(Nodes.Nod_clustID(keep)))
    
    idx_temp=find(Nodes.Nod_clustID(keep)==oldclusttag(i));
    
    if i<4    
    ClustID_CL4(idx_temp)=1;
    elseif i==4
   ClustID_CL4(idx_temp)=2;
   elseif i==5
   ClustID_CL4(idx_temp)=3;
   elseif i==6
   ClustID_CL4(idx_temp)=4;
    end
end
 
clustnames_CL4={'fasthab','modhab','weakhab','inhib'};




%%
%%%% now, I will be testing this using the results from the consensus attempt

load('S_cons_testing_GnO_measures6.mat');


%% first, just looking at the general values
%%%% in general, it seems that fmr1 has higher flexibility values

counter=1;
figure;
for group=[3 1 2]
   subplot(3,3,counter);imagesc(FlexMat.(groupnames{group}));caxis([0 1]);
   title(strcat('flex/',groupnames{group},'/med=',num2str(median(FlexMat.(groupnames{group})(:)))));
   
   counter=counter+1;
end
for group=[3 1 2]
   subplot(3,3,counter);imagesc(CoheMat.(groupnames{group}));caxis([0 1]);
   title(strcat('cohe/',groupnames{group},'/med=',num2str(median(CoheMat.(groupnames{group})(:)))));
   
   counter=counter+1;
end
for group=[3 1 2]
   subplot(3,3,counter);imagesc(PromMat.(groupnames{group}));caxis([0 1]);
   title(strcat('prom/',groupnames{group},'/med=',num2str(median(PromMat.(groupnames{group})(:)))));
   
   counter=counter+1;
end


counter=1;
test=[];
for group=[3 1 2]
    test(:,counter)=FlexMat.(groupnames{group})(:);
    
    counter=counter+1;
end


figure;boxplot(test);

%     [p, tbl, stats]=anova1(test);
%     figure;
%     [c, m]=multcompare(stats,'CType','bonferroni');
    
   %%%%OR
   for group=[1 2 3]  %%% cause now they are in the right order
   kstest(test(:,group))
   end 
   
   %%%% they are not normal... so  kruskalwallis? or Friedman? I think it
   %%%% would be friedman cause is not replicates on the same conditions as
   %%%% gamma and omega change
        %[p, tbl, stats]=kruskalwallis(test);
        [p, tbl, stats]=friedman(test,1);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');
    
    
    %%%%% would it be better to do a ranksum comparison? it looks at the
    %%%%% differences in the median. I think this is more accurate. In this
    %%%%% case the difference between WTvsfmr1 is also sig. 
    
%     [p, tbl, stats]=ranksum(test(:,1),groups_flexibilityTest(:,3));
%     p
    



%%

%%% looking at the concensus dynamic community detections
counter=1;
figure;
for i=1:25
    
    for j=1:20
        
   subplot(25,20,counter);imagesc(Big_FMR1_OPT{i,j}.control.S_cons);
   counter=counter+1;
    end
    
end

%%%% looking at a random raw sample of the community detection
counter=1;
figure;
for i=1:25
    
    for j=1:20
      sample=randperm(100,1);  
   subplot(25,20,counter);imagesc(Big_FMR1_OPT{i,j}.control.S_test(:,:,sample));
   counter=counter+1;
    end
    
end

%% stablishing limits

MAT_O_limits=[];
MAT_G_limits=[];
for i=1:25
    
    for j=1:20
        
        temp_mat=Big_FMR1_OPT{i,j}.control.S_cons;
        nodes=size(temp_mat,1);
        times=size(temp_mat,2);
        
        %%% testing for number of communities
        temp=[];
        for t=1:times
            temp(t)=max(temp_mat(:,t));           
        end
        
        if max(temp)>60  %%% I am testing with a quarter too
            
           MAT_G_limits(i,j)=1;
         elseif min(temp)<4  %%% 4 as is the number of my hab clusters... 
             MAT_G_limits(i,j)=-1;
        else
            MAT_G_limits(i,j)=0;
        end
            
        %%%% testing changes in time (min and max omega)
        temp=[];
        temp_lim=[];
        for n=1:nodes
            
            for t=2:length(temp_mat(n,:)) %%%% i am taking from timepoint 2 to compare with previous timepoint
                if temp_mat(n,t)==temp_mat(n,t-1)
                    temp(n,t)=1;
                else
                    temp(n,t)=0;
                end
                               
            end
            
            if temp(n,2)==0 && temp(n,12)==0 %%% to check the first and recovery looms
                temp_lim(n)=0;
            else
                temp_lim(n)=1;
            end
        end
        
        if sum(temp_lim)>60 %%%% i am not taking into account the 1 node cause it never changes community as is the reference node for being the first    
        MAT_O_limits(i,j)=1; 
        else
        MAT_O_limits(i,j)=0;    
        end
    
        
    end
    
end

figure;imagesc(MAT_G_limits);
figure;imagesc(MAT_O_limits);

%%%% to check the differences in flexibility with the
%%%% checking_GnO_measures.m script

good=intersect(find(MAT_G_limits==0),find(MAT_O_limits==0));

goodMat=zeros(size(MAT_G_limits));
goodMat(good)=1;
figure;imagesc(goodMat);

%[g,o]=find(goodMat);


%%

%%%%% now getting the null model differences 

load('QoptFMR1_All2.mat');



%%


Allthebest_fmr1=[];
Allverybest_fmr1=[];
for group=1:3
    temp_dif=QoptFMR1_All.(groupnames{group}).meanQopt2-QoptFMR1_All.(groupnames{group}).meanQopt2_null2;
    
    %figure;imagesc(temp_dif);
           
temp_relative_var=-(QoptFMR1_All.(groupnames{group}).varQopt2-max(QoptFMR1_All.(groupnames{group}).varQopt2(:)));
figure;imagesc(temp_relative_var);

thegood=temp_dif.*temp_relative_var;

%figure;imagesc(thegood);

thebest=thegood;thebest(thegood<0)=0;

figure;imagesc(thebest);

Allthebest_fmr1=cat(3,Allthebest_fmr1,thebest);

%%%% I will also look at the higher 1 or 2SD

x=mean(thebest(:));
y=std(thebest(:));

good_idx=find(thebest>(x+2*y));
%good_idx=find(thebest>x);

very_best=zeros(size(thebest));

very_best(good_idx)=1;

Allverybest_fmr1=cat(3,Allverybest_fmr1,very_best);

%figure;imagesc(very_best);
      
end

%%
Allthebest_all=Allthebest_fmr1;

Opt=mean(Allthebest_all,3);

figure;
imagesc(Opt);

%%%% to look at the parts with more difference with the null model
% 
% x=mean(Opt(:));
% y=std(Opt(:));
% 
% %good_idx=find(Opt>x);
% %good_idx=find(Opt>(x+1*y));
% good_idx=find(Opt>(x+2*y));
% 
% Opt_best=zeros(size(Opt));
% 
% Opt_best(good_idx)=1;
% 
% figure;imagesc(Opt_best);

%[g,o]=find(Opt_best);



%%
%%%%%% to look at the difference with the null model but based on limits
%figure;imagesc(Opt);

 Opt_new2=zeros(25,20);
 Opt_new2(find(goodMat))=Opt(find(goodMat));
 
 
figure;imagesc(Opt_new2);


 x=mean(Opt_new2(find(Opt_new2)));
 y=std(Opt_new2(find(Opt_new2)));

 good_idx=find(Opt_new2>x);
% good_idx=find(Opt_new2>(x+1*y));
%good_idx=find(Opt_new2>(x+2*y));

Opt_best2=zeros(size(Opt_new2));

Opt_best2(good_idx)=1;

figure;imagesc(Opt_best2);

[g,o]=find(Opt_best2);




counter=1;
test=struct;
for i=1:length(o)
%%%% comparing flexibility


 %for g=7
 %for o=9
    
    for group=1:3   
    
    test.(groupnames{group})(counter,:)=Big_FMR1_OPT{g(i),o(i)}.(groupnames{group}).flex;
    end
    counter=counter+1;
 %end
end


%%
%groups_flexibilityTest=NaN(90,3);
groups_flexibilityTest=NaN(length(o),3);
counter=1;
for group=[3 1 2]
    groups_flexibilityTest(1:length(o),counter)=mean(test.(groupnames{group})'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes. 
    %groups_flexibilityTest(1:90,counter)=mean(test.(groupnames{group}));
    %groups_flexibilityTest(1:90,counter)=test.(groupnames{group});
counter=counter+1;
end

figure;boxplot(groups_flexibilityTest);


%%%% they are not normal... so  kruskalwallis? or Friedman? I think it
   %%%% would be friedman cause is not replicates on the same conditions as
   %%%% gamma and omega change
    
%     [p, tbl,
%     stats]=ranksum(groups_flexibilityTest(:,1),groups_flexibilityTest(:,3));
%      %%%%% the ranksum test also gives sig. 
%     p
    
    [p, tbl, stats]=friedman(groups_flexibilityTest);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');

    
 %%
    
 
%%% to check which ones is worth plotting per brain region
Flex_perBrain=struct;
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    temp_idx=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);
    temp_groups=[];
       
        counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Flex_perBrain.(RegionList{brain})=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    figure;title(RegionList{brain});
    [c, m]=multcompare(stats,'CType','bonferroni');
end

%%% 




%%%%%%%% now by cluster with the CL4
Flex_perClust4=struct;
for clust=unique(ClustID_CL4(keepFmr1))'
    temp_idx=find(ClustID_CL4(keepFmr1)==clust);
    temp_groups=[];
    
        counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Flex_perClust4.(clustnames_CL4{clust})=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    [c, m]=multcompare(stats,'CType','bonferroni');
end

%%%%% for the subtypes in the OT

Flex_OT_CL4=struct;
brain=6;
temp_idx1=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);

for clust=1:3 %% I am not doing inhib cause there is only 1 node
    temp_idx=find(ClustID_CL4(keepFmr1)==clust);
    temp_idx=intersect(temp_idx,temp_idx1);
    temp_groups=[];
    counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Flex_OT_CL4.(clustnames_CL4{clust})=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    [c, m]=multcompare(stats,'CType','bonferroni');
end



    %% making figures
    %%%% for this one I need to add the cbrewer
    
    groups_flexibilityTest=NaN(90,3);
%groups_flexibilityTest=NaN(length(o),3);
counter=1;
for group=[3 1 2]
    %groups_flexibilityTest(1:length(o),counter)=mean(test.(groupnames{group})'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes. 
    groups_flexibilityTest(1:90,counter)=mean(test.(groupnames{group}));
    %groups_flexibilityTest(1:90,counter)=test.(groupnames{group});
counter=counter+1;
end
    
flex_dif=groups_flexibilityTest(:,1)-groups_flexibilityTest(:,3);
%figure;histogram(flex_dif);

figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
%   hold on;
%   gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),30,flex_dif,'filled');%colormap('RdBu');%caxis([-0.6 0.6]);%colorbar;
view(-90,90);
%title('top 1SD dif flex');
hold off;


%%

%% now cohesion 


counter=1;
test2=struct;
for i=1:length(o)
% for g=7:10
% for o=8:10
    
    for group=1:3   
    
    test2.(groupnames{group})(counter,:)=Big_FMR1_OPT{g(i),o(i)}.(groupnames{group}).node_cohesion;
    end
    counter=counter+1;
% end
% end
end

%groups_CoheTest=NaN(90,3);
groups_CoheTest=NaN(length(o),3);
counter=1;
for group=[3 1 2]
    groups_CoheTest(1:length(o),counter)=mean(test2.(groupnames{group})');
counter=counter+1;
end

figure;boxplot(groups_CoheTest);

    [p, tbl, stats]=friedman(groups_CoheTest);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');

    
    
    
%%% to check which ones is worth plotting per brain region
Cohe_perBrain=struct;
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    temp_idx=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);
    temp_groups=[];
       
        counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test2.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Cohe_perBrain.(RegionList{brain})=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    figure;title(RegionList{brain});
    [c, m]=multcompare(stats,'CType','bonferroni');
end

%%% 


%%%%%%%% now by cluster with the CL4
Cohe_perClust4=struct;
for clust=unique(ClustID_CL4(keepFmr1))'
    temp_idx=find(ClustID_CL4(keepFmr1)==clust);
    temp_groups=[];
    
        counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test2.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Cohe_perClust4.(clustnames_CL4{clust})=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    [c, m]=multcompare(stats,'CType','bonferroni');
end

%%%%% for the subtypes in the OT

Cohe_OT_CL4=struct;
brain=6;
temp_idx1=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);

for clust=1:3 %% I am not doing inhib cause there is only 1 node
    temp_idx=find(ClustID_CL4(keepFmr1)==clust);
    temp_idx=intersect(temp_idx,temp_idx1);
    temp_groups=[];
    counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test2.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Cohe_OT_CL4.(clustnames_CL4{clust})=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    [c, m]=multcompare(stats,'CType','bonferroni');
end


    %% making figures
    %%%% for this one I need to add the cbrewer
    
    groups_CoheTest=NaN(90,3);
%groups_flexibilityTest=NaN(length(o),3);
counter=1;
for group=[3 1 2]
    %groups_flexibilityTest(1:length(o),counter)=mean(test.(groupnames{group})'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes. 
    groups_CoheTest(1:90,counter)=mean(test2.(groupnames{group}));
    %groups_flexibilityTest(1:90,counter)=test.(groupnames{group});
counter=counter+1;
end
    
cohe_dif=groups_CoheTest(:,1)-groups_CoheTest(:,3);
%figure;histogram(flex_dif);

figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
%   hold on;
%   gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),30,cohe_dif,'filled');%colormap('RdBu');%caxis([-0.6 0.6]);%colorbar;
view(-90,90);
%title('top 1SD dif flex');
hold off;

    

    
%% and promiscuity
 

counter=1;
test3=struct;
for i=1:length(o)
% for g=7:10
% for o=8:10
    
    for group=1:3   
    
    test3.(groupnames{group})(counter,:)=Big_FMR1_OPT{g(i),o(i)}.(groupnames{group}).P;
    end
    counter=counter+1;
% end
end


groups_PromTest=NaN(length(o),3);
counter=1;
for group=[3 1 2]
    groups_PromTest(1:length(o),counter)=mean(test3.(groupnames{group})');
counter=counter+1;
end

figure;boxplot(groups_PromTest);

    [p, tbl, stats]=friedman(groups_PromTest);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');

    
    
       
%%% to check which ones is worth plotting per brain region
Prom_perBrain=struct;
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    temp_idx=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);
    temp_groups=[];
       
        counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test3.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Prom_perBrain.(RegionList{brain})=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    figure;title(RegionList{brain});
    [c, m]=multcompare(stats,'CType','bonferroni');
end

%%% 


%%%%%%%% now by cluster with the CL4
Prom_perClust4=struct;
for clust=unique(ClustID_CL4(keepFmr1))'
    temp_idx=find(ClustID_CL4(keepFmr1)==clust);
    temp_groups=[];
    
        counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test3.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Prom_perClust4.(clustnames_CL4{clust})=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    [c, m]=multcompare(stats,'CType','bonferroni');
end

%%%%% for the subtypes in the OT

Prom_OT_CL4=struct;
brain=6;
temp_idx1=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);

for clust=1:3 %% I am not doing inhib cause there is only 1 node
    temp_idx=find(ClustID_CL4(keepFmr1)==clust);
    temp_idx=intersect(temp_idx,temp_idx1);
    temp_groups=[];
    counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test3.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Prom_OT_CL4.(clustnames_CL4{clust})=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    [c, m]=multcompare(stats,'CType','bonferroni');
end


    %% making figures
    %%%% for this one I need to add the cbrewer
    
    groups_PromTest=NaN(90,3);
%groups_flexibilityTest=NaN(length(o),3);
counter=1;
for group=[3 1 2]
    %groups_flexibilityTest(1:length(o),counter)=mean(test.(groupnames{group})'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes. 
    groups_PromTest(1:90,counter)=mean(test3.(groupnames{group}));
    %groups_flexibilityTest(1:90,counter)=test.(groupnames{group});
counter=counter+1;
end
    
Prom_dif=groups_PromTest(:,1)-groups_PromTest(:,3);
%figure;histogram(flex_dif);

figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
%   hold on;
%   gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),30,Prom_dif,'filled');%colormap('RdBu');%caxis([-0.6 0.6]);%colorbar;
view(-90,90);
%title('top 1SD dif flex');
hold off;

    
    
%% checking also the number of communities:


%%% both hets and fmr1 seem to have a slower decay in number of comunities

CommNumMats=struct;
for group=[3 1 2]
temp=[];
for g=1:25
    for o=1:20
    CommNumMats.(groupnames{group})(g,o)=length(unique(Big_FMR1_OPT{g,o}.(groupnames{group}).S_cons(:)));
    end
    
end

end

%% comparing
counter=1;
for group=[3 1 2]
  
    temp(:,counter)=CommNumMats.(groupnames{group})(:);
    counter=counter+1;
end


    [p, tbl, stats]=friedman(temp);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');
    
   
%% making a figure
for group=[3 1 2]
max(CommNumMats.(groupnames{group})(:))
min(CommNumMats.(groupnames{group})(:))
end 
counter=1;
figure;
for group=[3 1 2]
   subplot(1,3,counter);imagesc(CommNumMats.(groupnames{group}));caxis([0 90]);
   title(strcat('#comm/',groupnames{group},'/med=',num2str(median(CommNumMats.(groupnames{group})(:)))));
   
   counter=counter+1;
end



%% # of communities by time in 3D

%%% taking into account the time (loom presentation
%%% not really sure it will be useful


CommNumMats3D=struct;
for group=[3 1 2]
temp=[];
for g=1:25
    for o=1:20
    for loom=1:21
        CommNumMats3D.(groupnames{group})(g,o,loom)=length(unique(Big_FMR1_OPT{g,o}.(groupnames{group}).S_cons(:,loom)));
    end
    end
end

end



%% comparing
counter=1;
for group=[3 1 2]
  
    temp(:,counter)=CommNumMats3D.(groupnames{group})(:);
    counter=counter+1;
end


    [p, tbl, stats]=friedman(temp);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');


    
%%%%%% trying as 3 different multiplots

for group=[3 1 2]
   counter=1;
    figure;
   for loom=1:size(CommNumMats3D.(groupnames{group}),3)
    subplot(3,7,counter);imagesc(CommNumMats3D.(groupnames{group})(:,:,loom));caxis([0 90]);
    counter=counter+1;
   end
end


%%%%% there are some small diferences only. probably nothing relevant. 

%%



%% ploting specific community detections

%%%% as an example, picking one that is included in our analysis

%% dynamic network community detection


figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
sgtitle(strcat('g=',num2str(1.6),'/','o=',num2str(0.9)));
counter=1;
for group=[3 1 2]
S_good=Big_FMR1_OPT{16,9}.(groupnames{group}).S_cons;

low=min(min(S_good));
high=max(max(S_good));

 subplot(1,3,counter);imagesc(S_good);colormap('jet');caxis([low high]);colorbar;
  title(groupnames{group}); 
 counter=counter+1;
    
 end   
% counter=counter+1;
% for data=1:4
% S_good=S_cons_FnS_All.S_cons.(datasets(data,:)).S_cons;
% 
% low=min(min(S_good));
% high=max(max(S_good));
% 
%  subplot(2,4,counter);imagesc(S_good);colormap('jet');caxis([low high]);colorbar;
%  title(datasets(data,:)); 
%  counter=counter+1;
%     
% end   

saveas(gcf,'communities_fmr1_g16o09.svg')
close 

%% making a figure with the flex,cohe and prom matrices of all datasets


GnO_FnS=load('S_cons_testing_GnO_measures4.mat','FlexMat','CoheMat','PromMat');

counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
for data=1:4
   subplot(3,7,counter);imagesc(GnO_FnS.FlexMat.(datasets(data,:)));caxis([0 1]);colormap(inferno);colorbar;
   tempMean=mean(GnO_FnS.FlexMat.(datasets(data,:))(:));
   tempSD=std(GnO_FnS.FlexMat.(datasets(data,:))(:));
   title(strcat('flex/',datasets(data,:),'/M=',num2str(round(tempMean,2)),'/SD=',num2str(round(tempSD,2))));
   
   counter=counter+1;
end
for group=[3 1 2]
   subplot(3,7,counter);imagesc(FlexMat.(groupnames{group}));caxis([0 1]);colormap(inferno);colorbar;
   tempMean=mean(FlexMat.(groupnames{group})(:));
   tempSD=std(FlexMat.(groupnames{group})(:));
   title(strcat('flex/',groupnames{group},'/M=',num2str(round(tempMean,2)),'/SD=',num2str(round(tempSD,2))));
   
   counter=counter+1;
end
for data=1:4
   subplot(3,7,counter);imagesc(GnO_FnS.CoheMat.(datasets(data,:)));caxis([0 1]);colormap(inferno);colorbar;
   tempMean=mean(GnO_FnS.CoheMat.(datasets(data,:))(:));
   tempSD=std(GnO_FnS.CoheMat.(datasets(data,:))(:));
   title(strcat('cohe/',datasets(data,:),'/M=',num2str(round(tempMean,2)),'/SD=',num2str(round(tempSD,2))));
   
   counter=counter+1;
end
for group=[3 1 2]
   subplot(3,7,counter);imagesc(CoheMat.(groupnames{group}));caxis([0 1]);colormap(inferno);colorbar;
   tempMean=mean(CoheMat.(groupnames{group})(:));
   tempSD=std(CoheMat.(groupnames{group})(:));
   title(strcat('Cohe/',groupnames{group},'/M=',num2str(round(tempMean,2)),'/SD=',num2str(round(tempSD,2))));
   
   counter=counter+1;
end
for data=1:4
   subplot(3,7,counter);imagesc(GnO_FnS.PromMat.(datasets(data,:)));caxis([0 1]);colormap(inferno);colorbar;
   tempMean=mean(GnO_FnS.PromMat.(datasets(data,:))(:));
   tempSD=std(GnO_FnS.PromMat.(datasets(data,:))(:));
   title(strcat('prom/',datasets(data,:),'/M=',num2str(round(tempMean,2)),'/SD=',num2str(round(tempSD,2))));
   
   counter=counter+1;
end
for group=[3 1 2]
   subplot(3,7,counter);imagesc(PromMat.(groupnames{group}));caxis([0 1]);colormap(inferno);colorbar;
   tempMean=mean(PromMat.(groupnames{group})(:));
   tempSD=std(PromMat.(groupnames{group})(:));
   title(strcat('prom/',groupnames{group},'/M=',num2str(round(tempMean,2)),'/SD=',num2str(round(tempSD,2))));
   
   counter=counter+1;
end

%saveas(gcf,'Measures_All.svg');
saveas(gcf,'Measures_All_wcolorbar.svg');

