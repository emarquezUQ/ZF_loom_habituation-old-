
%%%%% this script is to make the panels for the multilayer community
%%%%% detection figure. I used g=1.4 o=0.8.

load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

keepFmr1=load('graphs_fmr1Exp.mat','keep');
keepFmr1=keepFmr1.keep;

S_cons_FMR1=load('S_cons_FMR1_G14O07.mat');
S_cons_FnS_All=load('S_cons_FnS_All_G14O07.mat');

cbrewer()

[RdBu]=cbrewer('div','RdBu',101);


%%
%%%% with fasthab subtypes
oldclusttag=[4 2 5 6 1 7];
ClustID=zeros(size(Nodes.Nod_clustID(keep)));
for i=1:length(unique(Nodes.Nod_clustID(keep)))
    
    idx_temp=find(Nodes.Nod_clustID(keep)==oldclusttag(i));
    ClustID(idx_temp)=i;
   
end
 
clustnames={'fasthab1','fasthab2','fasthab3','modhab','weakhab','inhib'};

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


%% dynamic network community detection
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
sgtitle(strcat('g=',num2str(1.4),'/','o=',num2str(0.8)));
counter=1;
for group=[3 1 2]
S_good=S_cons_FMR1.S_cons.(groupnames{group}).S_cons;

low=min(min(S_good));
high=max(max(S_good));

 subplot(2,4,counter);imagesc(S_good);colormap('jet');caxis([low high]);colorbar;
  title(groupnames{group}); 
 counter=counter+1;
    
end   
counter=counter+1;
for data=1:4
S_good=S_cons_FnS_All.S_cons.(datasets(data,:)).S_cons;

low=min(min(S_good));
high=max(max(S_good));

 subplot(2,4,counter);imagesc(S_good);colormap('jet');caxis([low high]);colorbar;
 title(datasets(data,:)); 
 counter=counter+1;
    
end   

saveas(gcf,'communities_All2.svg')
close 

%% getting the values for flexibility plot of all datasets

groups_flexibility=NaN(99,7);
counter=1;
for group=[3 1 2]
    groups_flexibility(1:90,counter)=S_cons_FMR1.S_cons.(groupnames{group}).flex;
counter=counter+1;
end
for data=1:4
    groups_flexibility(:,counter)=S_cons_FnS_All.S_cons.(datasets(data,:)).flex;
counter=counter+1;
end

figure;boxplot(groups_flexibility);title('flexibility');ylim([0 1]);

%%%% just checking
[p, tbl, stats]=anova1(groups_flexibility);
figure;
[c, m]=multcompare(stats,'CType','bonferroni');

%%% then I took the groups_flexibility to graphpad to graph it

%% getting the values for cohesion plot of all datasets

groups_cohesion=NaN(99,7);
counter=1;
for group=[3 1 2]
    groups_cohesion(1:90,counter)=S_cons_FMR1.S_cons.(groupnames{group}).node_cohesion;
counter=counter+1;
end
for data=1:4
    groups_cohesion(:,counter)=S_cons_FnS_All.S_cons.(datasets(data,:)).node_cohesion;
counter=counter+1;
end

figure;boxplot(groups_cohesion);title('cohesion');ylim([0 1]);

%%%% just checking
[p, tbl, stats]=anova1(groups_cohesion);
figure;
[c, m]=multcompare(stats,'CType','bonferroni');

%%% then I took the groups_cohesion to graphpad to graph it


%% getting the values for promiscuity plot of all datasets

groups_prom=NaN(99,7);
counter=1;
for group=[3 1 2]
    groups_prom(1:90,counter)=S_cons_FMR1.S_cons.(groupnames{group}).P;
counter=counter+1;
end
for data=1:4
    groups_prom(:,counter)=S_cons_FnS_All.S_cons.(datasets(data,:)).P;
counter=counter+1;
end

figure;boxplot(groups_prom);title('promiscuity');ylim([0 1]);

%%%% just checking
[p, tbl, stats]=anova1(groups_prom);
figure;
[c, m]=multcompare(stats,'CType','bonferroni');

%%% then I took the groups_prom to graphpad to graph it

%% brain region and clusters comparison. For flexibility

%%% to check which ones is worth plotting
Flex_perBrain=struct;
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    temp_idx=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);
    temp_groups=[];
    for group=1:3 %%% keep in mind that I changed the order. now is controls, hets and fmr1
    temp_groups(:,group)=groups_flexibility(temp_idx,group);   
    end
    
    Flex_perBrain.(RegionList{brain})=temp_groups;
    
    [p, tbl, stats]=anova1(temp_groups);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');
end

%%% Hindbrain and tectum are among the sig ones. and maybe thalamus to show
%%% is not sig? I will do it in graphpad. the data is in Flex_perBrain

%%%%%%%% now by cluster with the CL4

for clust=unique(ClustID_CL4(keepFmr1))'
    temp_idx=find(ClustID_CL4(keepFmr1)==clust);
    temp_groups=[];
    for group=1:3 %%% keep in mind that I changed the order
    temp_groups(:,group)=groups_flexibility(temp_idx,group);   
    end
    
    Flex_perClust4.(clustnames_CL4{clust})=temp_groups;
    
    [p, tbl, stats]=anova1(temp_groups);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');
end


%% relative difference in flex plot


flex_dif=S_cons_FMR1.S_cons.(groupnames{3}).flex-S_cons_FMR1.S_cons.(groupnames{2}).flex;
%figure;histogram(flex_dif);

figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
%   hold on;
%   gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),30,flex_dif,'filled');colormap(RdBu);caxis([-0.3 0.3]);colorbar;
view(-90,90);
%title('top 1SD dif flex');
hold off;

saveas(gcf,'Flex_WTvsfmr1_2.svg')
%saveas(gcf,'Flex_WTvsfmr1_wcolorbar2.svg')
close 

% top_P_dif=find(flex_dif>0.1);
% 
% figure;
% subplot(1,2,1);histogram(NodesFmr1.Nod_brainID(keepFmr1)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);hold on;histogram(NodesFmr1.Nod_brainID(keepFmr1(top_P_dif)));
% subplot(1,2,2);histogram(ClustID(keepFmr1)); xticks([1:6]); xticklabels(clustnames);xtickangle(45);hold on;histogram(ClustID(keepFmr1(top_P_dif))); 

