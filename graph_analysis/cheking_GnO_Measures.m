

%%%%% this script is a coninuation of
%%%%% testing_GenLouvain_optimazing_All_wTemporalNull.m to explore the
%%%%% dynamic community measures based on the optimizatin of gamma and
%%%%% omega. 


%%% Note: g and o values are calculated in
%%% testing_GenLouvain_optimazing_All_wTemporalNull.m code, unless I choose
%%% specific ones. 
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
groups_flexibilityTest=NaN(90,3);
counter=1;
for group=[3 1 2]
    groups_flexibilityTest(1:90,counter)=mean(test.(groupnames{group}));
    %groups_flexibilityTest(1:90,counter)=test.(groupnames{group});
counter=counter+1;
end

figure;boxplot(groups_flexibilityTest);

    [p, tbl, stats]=anova1(groups_flexibilityTest);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');

    
    
    %%
    
 
%%% to check which ones is worth plotting
Flex_perBrain=struct;
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    temp_idx=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);
    temp_groups=[];
    for group=1:3 %%% keep in mind that I changed the order. now is controls, hets and fmr1
    temp_groups(:,group)=groups_flexibilityTest(temp_idx,group);   
    end
    
    Flex_perBrain.(RegionList{brain})=temp_groups;
    
    [p, tbl, stats]=anova1(temp_groups);
    figure;title(RegionList{brain});
    [c, m]=multcompare(stats,'CType','bonferroni');
end

%%% 

%%%%%%%% now by cluster with the CL4

for clust=unique(ClustID_CL4(keepFmr1))'
    temp_idx=find(ClustID_CL4(keepFmr1)==clust);
    temp_groups=[];
    for group=1:3 %%% keep in mind that I changed the order
    temp_groups(:,group)=groups_flexibilityTest(temp_idx,group);   
    end
    
    Flex_perClust4.(clustnames_CL4{clust})=temp_groups;
    
    [p, tbl, stats]=anova1(temp_groups);
    figure;title(clustnames_CL4{clust});
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
    for group=1:3 %%% keep in mind that I changed the order
    temp_groups(:,group)=groups_flexibilityTest(temp_idx,group);   
    end
    
    Flex_OT_CL4.(clustnames_CL4{clust})=temp_groups;
    
    [p, tbl, stats]=anova1(temp_groups);
    [c, m]=multcompare(stats,'CType','bonferroni');
end


    %%
    
    
flex_dif=groups_flexibilityTest(:,1)-groups_flexibilityTest(:,3);
%figure;histogram(flex_dif);

figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
%   hold on;
%   gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),30,flex_dif,'filled');colormap('RdBu');%caxis([-0.6 0.6]);%colorbar;
view(-90,90);
%title('top 1SD dif flex');
hold off;

%%

%%%% now cohesion 
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

groups_CoheTest=NaN(90,3);
counter=1;
for group=[3 1 2]
    groups_CoheTest(1:90,counter)=mean(test2.(groupnames{group}));
counter=counter+1;
end

figure;boxplot(groups_CoheTest);

    [p, tbl, stats]=anova1(groups_CoheTest);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');

    
    
 %%%% and promiscuity
 

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


groups_PromTest=NaN(90,3);
counter=1;
for group=[3 1 2]
    groups_PromTest(1:90,counter)=mean(test3.(groupnames{group}));
counter=counter+1;
end

figure;boxplot(groups_PromTest);

    [p, tbl, stats]=anova1(groups_PromTest);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');

    
%% checking also the number of communities:


%%% both hets and fmr1 seem to have a slower decay in number of comunities

CommNumMats=struct;
for group=[3 1 2]
temp=[];
for g=1:16
    for o=1:10
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


    [p, tbl, stats]=anova1(temp);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');


%% making a figure
for group=[3 1 2]
max(CommNumMats.(groupnames{group})(:))
min(CommNumMats.(groupnames{group})(:))
end

%%%% max 90 and min 3

counter=1;
figure;
for group=[3 1 2]
   subplot(1,3,counter);imagesc(CommNumMats.(groupnames{group}));caxis([0 90]);
   title(strcat('#Comm/',groupnames{group}));
   
   counter=counter+1;
end


%% # of communities by time in 3D


%%% both hets and fmr1 seem to have a slower decay in number of comunities

CommNumMats3D=struct;
for group=[3 1 2]
temp=[];
for g=1:16
    for o=1:10
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


    [p, tbl, stats]=anova1(temp);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');


%% making a figure



%%%% mmm this is not working... not sure how I could visualize it

for group=[3 1 2]
max(CommNumMats3D.(groupnames{group})(:))
min(CommNumMats3D.(groupnames{group})(:))
end

%%%% max 90 and min 3



counter=1;
figure;
for group=[3 1 2]
   subplot(3,3,counter);surf(CommNumMats3D.(groupnames{group}));caxis([0 90]);
   title(strcat('#Comm/',groupnames{group}));
   
   counter=counter+1;
end
for group=[3 1 2]
   subplot(3,3,counter);surf(CommNumMats3D.(groupnames{group}));caxis([0 90]);
   title(strcat('#Comm/',groupnames{group}));
   
   counter=counter+1;
end
for group=[3 1 2]
   subplot(3,3,counter);surf(CommNumMats3D.(groupnames{group}));caxis([0 90]);
   title(strcat('#Comm/',groupnames{group}));
   
   counter=counter+1;
end

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



%% ploting specific community detections

%% dynamic network community detection


figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
sgtitle(strcat('g=',num2str(0.8),'/','o=',num2str(0.9)));
counter=1;
for group=[3 1 2]
S_good=Big_FMR1_OPT{8,9}.(groupnames{group}).S_cons;

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

saveas(gcf,'communities_fmr1_g08o09.svg')
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