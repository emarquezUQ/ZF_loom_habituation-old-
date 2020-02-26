
%%%%% acording to the comparisons with a nodal null model (see
%%%%% QoptFMR1_WT.mat from testing_GenLouvain_optimazingFMR1_controls.m
%%%%% script). some of the optimal values could be g=0.7 o=0.5; g=1.2 o=0.5; 
%%%%% g=1.4 o=0.8. The last one is the one that convinces me the most. 
%%%%% i also tried with g=1.2 o=0.5 but the results are much more
%%%%% homogenic. 


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

figure;
subplot(1,2,1);histogram(Nodes.Nod_brainID(keep)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);
subplot(1,2,2);histogram(ClustID); xticks([1:6]); xticklabels(clustnames);xtickangle(45);

%% testing values


%%
for group=1:3
A={};

for loom=1:21
  
   temp=Data_corrMat4.(groupnames{group}).Mean_corrMat{1,loom}(keepFmr1,keepFmr1);
   temp(isnan(temp))=0;
    A{loom}=temp;  
end
clear temp

%%


gamma = 0.9;
omega = 0.8; %%% this makes a big influence... 

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
S = reshape(S,N,T);

%% checking results

%figure;imagesc(S);colormap('jet');caxis([1 5]);


%%


time=[1 2 6 11 12];

low=min(min(S));
high=max(max(S));
c=jet(high);
%c=c(1:5,:);

%figure;imagesc(S);colormap(c);colorbar;

  figure;
  
  subplot(3,6,1);imagesc(S);colormap('jet');caxis([low high]);
  subplot(3,6,7);imagesc(S);colormap('jet');caxis([low high]);
  subplot(3,6,13);imagesc(S);colormap('jet');caxis([low high]);
  temp=[];
  for i=1:5
  subplot(3,6,i+1);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
%   gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],21,'off'); 
%     hold on;
    scatter(Nodes.Nod_coor(keep(keepFmr1),1),Nodes.Nod_coor(keep(keepFmr1),2),20,S(:,time(i)),'filled');colormap('jet');caxis([low high]);
    view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  

    
subplot(3,6,i+7);histogram(Nodes.Nod_brainID(keep(keepFmr1))); xticks([1:9]); xticklabels(RegionList);xtickangle(45);
for j=1:max(max(S))
    hold on;histogram(Nodes.Nod_brainID(keep(keepFmr1(find(S(:,time(i))==j)))),'DisplayStyle','stairs','EdgeColor',c(j,:));
end

subplot(3,6,i+13);histogram(ClustID(keepFmr1)); xticks([1:6]); xticklabels(clustnames);xtickangle(45);
for j=1:max(max(S))
hold on;histogram(ClustID(keepFmr1(find(S(:,time(i))==j))),'DisplayStyle','stairs','EdgeColor',c(j,:)); 
end
      
  end
  
end


%% making 1000 iterations per group to calculate a representative modularity 
S_cons=struct;
%%
for group=1:3
A={};

for loom=1:21
  
   temp=Data_corrMat4.(groupnames{group}).Mean_corrMat{1,loom}(keepFmr1,keepFmr1);
   temp(isnan(temp))=0;
    A{loom}=temp;  
end
clear temp

%%

%%%% making a multidimensional structure where to store things
S_test=[];

for test=1:1000
     
gamma = 0.9;
omega = 0.8; %%% this makes a big influence... 

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
S = reshape(S,N,T);

S_test=cat(3,S_test,S);
  
  
end


%% trying consensus_iterative 


S_good=[];
for i=1:21
   C=S_test(:,i,:); 
   C=squeeze(C);
   [S2, Q2, X_new3, qpc] = consensus_iterative(C');
    S_good(:,i)=S2(i,:);
end

S_cons.(groupnames{group}).S_test=S_test;
S_cons.(groupnames{group}).S_cons=S_good;


low=min(min(S_good));
high=max(max(S_good));

time=[1 2 6 11 12];

  figure;
  
  subplot(1,7,1);imagesc(S_good);colormap('jet');caxis([low high]);
  temp=[];
  for i=1:5
  subplot(1,7,i+1);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
%    gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],21,'off'); 
%      hold on;
    scatter(Nodes.Nod_coor(keepFmr1,1),Nodes.Nod_coor(keepFmr1,2),20,S_good(:,time(i)),'filled');colormap('jet');caxis([low high]);
    view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  
     
  end
  subplot(1,7,7);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
   gscatter(Nodes.Nod_coor(keepFmr1,1),Nodes.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],20,'off');     
     view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  
    
    sgtitle(strcat('g=',num2str(gamma),'/','o=',num2str(omega)));
end



%% # of communities

%%% both hets and fmr1 seem to have a slower decay in number of comunities

figure;
for group=[3 1 2]
temp=[];
for loom=1:21
temp(loom)=length(unique(S_cons.(groupnames{group}).S_cons(:,loom)));
end
plot(temp);
hold on;

end
hold off;
clear temp
%% looking at flexibility and other measures


%% testing flexibility

for group=1:3

flex=flexibility(S_cons.(groupnames{group}).S_cons','temp'); %%% i need to mind the orientation of the matrix to do it properly

S_cons.(groupnames{group}).flex=flex;

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_clustID(keepFmr1,:),'rgggbm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),20,flex,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title('flexibility');
hold off;


clustnames={'strhab1','strhab2','strhab3','modhab','weakhab','inhib'};

    
figure;
subplot(1,2,1);boxplot(flex,NodesFmr1.Nod_brainID(keepFmr1)),xticklabels(RegionList);xtickangle(45);ylim([0 1]);
subplot(1,2,2);boxplot(flex,ClustID(keepFmr1)),xticklabels(clustnames);xtickangle(45);ylim([0 1]);

end

%% cohesion strength and related

options.figureFlag	= 0;
options.colormap	= 'jet';

for group=1:3

[Cij,node_cohesion,node_disjoint,node_flexibility,strength_cohesion,commChanges,commCohesion,commDisjoint,commIndex] = calc_node_cohesion(S_cons.(groupnames{group}).S_cons,options);

S_cons.(groupnames{group}).node_cohesion=node_cohesion;

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
% hold on;
% gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_clustID(keepFmr1,:),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),20,node_cohesion,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title('cohesion');
hold off;


figure;
subplot(1,2,1);boxplot(node_cohesion,NodesFmr1.Nod_brainID(keepFmr1)),xticklabels(RegionList);xtickangle(45);ylim([0 1]);
subplot(1,2,2);boxplot(node_cohesion,ClustID(keepFmr1)),xticklabels(clustnames);xtickangle(45);ylim([0 1]);

end
%%

for group=1:3

P = promiscuity(S_cons.(groupnames{group}).S_cons');  %%% i need to mind the orientation of the matrix to do it properly

S_cons.(groupnames{group}).P=P;

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
% hold on;
% gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_clustID(keepFmr1,:),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),20,P,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title('promiscuity');
hold off;
  
figure;
subplot(1,2,1);boxplot(P,NodesFmr1.Nod_brainID(keepFmr1)),xticklabels(RegionList);xtickangle(45);ylim([0 1]);
subplot(1,2,2);boxplot(P,ClustID(keepFmr1)),xticklabels(clustnames);xtickangle(45);ylim([0 1]);

end

%%
counter=1;
figure;
for group=[3 1 2]
subplot(3,3,counter);

plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),20,S_cons.(groupnames{group}).flex,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title(strcat('flexibility','/',groupnames{group}));
hold off;

counter=counter+1;
end
for group=[3 1 2]
subplot(3,3,counter);

plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),20,S_cons.(groupnames{group}).node_cohesion,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title(strcat('cohesion','/',groupnames{group}));
hold off;

counter=counter+1;
end
for group=[3 1 2]
subplot(3,3,counter);

plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),20,S_cons.(groupnames{group}).P,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title(strcat('promiscuity','/',groupnames{group}));
hold off;

counter=counter+1;
end

%%


%%
groups_flexibility=[];
counter=1;
for group=[3 1 2]

    groups_flexibility(:,counter)=S_cons.(groupnames{group}).flex;

counter=counter+1;
end

figure;boxplot(groups_flexibility);title('flexibility');



groups_cohesion=[];
counter=1;
for group=[3 1 2]

    groups_cohesion(:,counter)=S_cons.(groupnames{group}).node_cohesion;

counter=counter+1;
end

figure;boxplot(groups_cohesion);title('cohesion');



groups_prom=[];
counter=1;
for group=[3 1 2]

    groups_prom(:,counter)=S_cons.(groupnames{group}).P;

counter=counter+1;
end

figure;boxplot(groups_prom);title('promiscuity');


figure;
subplot(1,3,1);boxplot(groups_flexibility);title('flexibility');ylim([0 1]);
subplot(1,3,2);boxplot(groups_cohesion);title('cohesion');ylim([0 1]);
subplot(1,3,3);boxplot(groups_prom);title('promiscuity');ylim([0 1]);

%%%%% comparing statistical differences

%%%% it seems that both hets and wt are sig different from fmr1!!

[p, tbl, stats]=anova1(groups_flexibility);
figure;
[c, m]=multcompare(stats,'CType','bonferroni');

[p, tbl, stats]=anova1(groups_cohesion);
figure;
[c, m]=multcompare(stats,'CType','bonferroni');

[p, tbl, stats]=anova1(groups_prom);
figure;
[c, m]=multcompare(stats,'CType','bonferroni');


%% comparing per brain region and functional cluster

%%%%%%% first with Brain regions
%%%% it seems that subpallium, pretectum, tectum and Hindbrain have sig
%%%% differences between fmr1 and WT. tegmentum has sig between fmr1 and
%%%% hets

for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    temp_idx=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);
    temp_groups=[];
    for group=1:3 %%% keep in mind that I changed the order
    temp_groups(:,group)=groups_flexibility(temp_idx,group);   
    end
    [p, tbl, stats]=anova1(temp_groups);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');
end

%%%%%%%% now by cluster
%%%%% it seems that the 3 fast hab clusters are sig. then the weakly
%%%%% habituating and inhibited are too. but hets is also sig lower than WT
%%%%% which is an odd result. for the moderately hab cluster fmr1 is sig
%%%%% different from hets but not WT. 

for clust=unique(ClustID(keepFmr1))'
    temp_idx=find(ClustID(keepFmr1)==clust);
    temp_groups=[];
    for group=1:3 %%% keep in mind that I changed the order
    temp_groups(:,group)=groups_flexibility(temp_idx,group);   
    end
    [p, tbl, stats]=anova1(temp_groups);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');
end

%%%% from this analysis it seems that the flexibility in green and red
%%%% clusters in the midbrain (pretectum and tectum) are the ones most
%%%% affected. 

%% testing some more things

%%% checking relationship between dynamic variables
for group=1:3

figure;
subplot(1,3,1);scatter(S_cons.(groupnames{group}).flex,S_cons.(groupnames{group}).node_cohesion,'filled');xlim([0 1]);ylim([0 1]);title('flexVScohe');
subplot(1,3,2);scatter(S_cons.(groupnames{group}).flex,S_cons.(groupnames{group}).P,'filled');xlim([0 1]);ylim([0 1]);title('flexVSprom');
subplot(1,3,3);scatter(S_cons.(groupnames{group}).node_cohesion,S_cons.(groupnames{group}).P,'filled');xlim([0 1]);ylim([0 1]);title('coheVSprom');

end


%%%%% looking for the top flexibility
for group=1:3

threshold=mean(S_cons.(groupnames{group}).flex)+1*std(S_cons.(groupnames{group}).flex);

idx_top=find(S_cons.(groupnames{group}).flex>threshold);

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
  hold on;
  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1(idx_top),1),NodesFmr1.Nod_coor(keepFmr1(idx_top),2),20,'k','filled');
view(-90,90);
title('top 1SD flexibility');
hold off;


figure;
subplot(1,2,1);histogram(NodesFmr1.Nod_brainID(keepFmr1)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);hold on;histogram(NodesFmr1.Nod_brainID(keepFmr1(idx_top)));
subplot(1,2,2);histogram(ClustID(keepFmr1)); xticks([1:6]); xticklabels(clustnames);xtickangle(45);hold on;histogram(ClustID(keepFmr1(idx_top))); 

end


%%%%% looking for the top cohesion
for group=1:3

threshold=mean(S_cons.(groupnames{group}).node_cohesion)+1*std(S_cons.(groupnames{group}).node_cohesion);

idx_top=find(S_cons.(groupnames{group}).node_cohesion>threshold);

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
  hold on;
  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1(idx_top),1),NodesFmr1.Nod_coor(keepFmr1(idx_top),2),20,'k','filled');
view(-90,90);
title('top 1SD cohesion');
hold off;


figure;
subplot(1,2,1);histogram(NodesFmr1.Nod_brainID(keepFmr1)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);hold on;histogram(NodesFmr1.Nod_brainID(keepFmr1(idx_top)));
subplot(1,2,2);histogram(ClustID(keepFmr1)); xticks([1:6]); xticklabels(clustnames);xtickangle(45);hold on;histogram(ClustID(keepFmr1(idx_top))); 

end

%%%%% looking for the top promiscuity
for group=1:3

threshold=mean(S_cons.(groupnames{group}).P)+1*std(S_cons.(groupnames{group}).P);

idx_top2=find(S_cons.(groupnames{group}).P>threshold);

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
  hold on;
  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1(idx_top2),1),NodesFmr1.Nod_coor(keepFmr1(idx_top2),2),20,'k','filled');
view(-90,90);
title('top 1SD prom');
hold off;


figure;
subplot(1,2,1);histogram(NodesFmr1.Nod_brainID(keepFmr1)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);hold on;histogram(NodesFmr1.Nod_brainID(keepFmr1(idx_top2)));
subplot(1,2,2);histogram(ClustID(keepFmr1)); xticks([1:6]); xticklabels(clustnames);xtickangle(45);hold on;histogram(ClustID(keepFmr1(idx_top2))); 

end



%% relative change


%%%%%% flexibility

flex_dif=S_cons.(groupnames{3}).flex-S_cons.(groupnames{2}).flex;
figure;histogram(flex_dif);
x=mean(flex_dif);
y=std(flex_dif);

top_flex_dif=find(flex_dif<(x-1*y));


figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
  hold on;
  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1(top_flex_dif),1),NodesFmr1.Nod_coor(keepFmr1(top_flex_dif),2),20,'k','filled');
view(-90,90);
title('top 1SD dif flex');
hold off;

figure;
subplot(1,2,1);histogram(NodesFmr1.Nod_brainID(keepFmr1)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);hold on;histogram(NodesFmr1.Nod_brainID(keepFmr1(top_flex_dif)));
subplot(1,2,2);histogram(ClustID(keepFmr1)); xticks([1:6]); xticklabels(clustnames);xtickangle(45);hold on;histogram(ClustID(keepFmr1(top_flex_dif))); 


%%%%%% cohesion

cohe_dif=S_cons.(groupnames{3}).node_cohesion-S_cons.(groupnames{2}).node_cohesion;
figure;histogram(cohe_dif);
x=mean(cohe_dif);
y=std(cohe_dif);

top_cohe_dif=find(cohe_dif<(x-1*y));


figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
  hold on;
  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1(top_cohe_dif),1),NodesFmr1.Nod_coor(keepFmr1(top_cohe_dif),2),20,'k','filled');
view(-90,90);
title('top 1SD dif cohe');
hold off;

figure;
subplot(1,2,1);histogram(NodesFmr1.Nod_brainID(keepFmr1)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);hold on;histogram(NodesFmr1.Nod_brainID(keepFmr1(top_cohe_dif)));
subplot(1,2,2);histogram(ClustID(keepFmr1)); xticks([1:6]); xticklabels(clustnames);xtickangle(45);hold on;histogram(ClustID(keepFmr1(top_cohe_dif))); 


%%%%%% promiscuity
P_dif=S_cons.(groupnames{3}).P-S_cons.(groupnames{2}).P;
figure;histogram(P_dif);
x=mean(P_dif);
y=std(P_dif);

top_P_dif=find(P_dif<(x-1*y));


figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
  hold on;
  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1(top_P_dif),1),NodesFmr1.Nod_coor(keepFmr1(top_P_dif),2),20,'k','filled');
view(-90,90);
title('top 1SD dif prom');
hold off;

figure;
subplot(1,2,1);histogram(NodesFmr1.Nod_brainID(keepFmr1)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);hold on;histogram(NodesFmr1.Nod_brainID(keepFmr1(top_P_dif)));
subplot(1,2,2);histogram(ClustID(keepFmr1)); xticks([1:6]); xticklabels(clustnames);xtickangle(45);hold on;histogram(ClustID(keepFmr1(top_P_dif))); 


%%%%% 


figure;
subplot(1,3,1);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),20,flex_dif,'filled'); colormap('cool'); caxis([-0.3 0.3]);colorbar;
view(-90,90);
title(strcat('flexibility','/wt-fmr1'));
hold off;
subplot(1,3,2);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),20,cohe_dif,'filled'); colormap('cool'); caxis([-0.3 0.3]);colorbar;
view(-90,90);
title(strcat('cohesion','/wt-fmr1'));
hold off;
subplot(1,3,3);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),20,P_dif,'filled'); colormap('cool'); caxis([-0.3 0.3]);colorbar;
view(-90,90);
title(strcat('promiscuity','/wt-fmr1'));
hold off;




%% checking flexibility per block

%%%% not sure what to make of it. in general the flexibility is higher
%%%% always for fmr1. so this could be a phenotype. but is hard to conclude
%%%% on the dynamics. 

for group=1:3

flex1=flexibility(S_cons.(groupnames{group}).S_cons(:,2:11)','temp'); %%% i need to mind the orientation of the matrix to do it properly

flex2=flexibility(S_cons.(groupnames{group}).S_cons(:,12:21)','temp');

S_cons.(groupnames{group}).flex1=flex1;
S_cons.(groupnames{group}).flex2=flex2;

mean(flex1)
mean(flex2)
end

%% saving

save('S_cons_FMR1_G09O08.mat','S_cons');

%%

%%%%% flexibility per blocks again but
%%%%% in 4 blocks of 5 looms

block_flex=[];
for group=1:3

flex1=flexibility(S_cons.(groupnames{group}).S_cons(:,2:6)','temp'); %%% i need to mind the orientation of the matrix to do it properly

flex2=flexibility(S_cons.(groupnames{group}).S_cons(:,7:11)','temp');

flex3=flexibility(S_cons.(groupnames{group}).S_cons(:,12:16)','temp'); %%% i need to mind the orientation of the matrix to do it properly

flex4=flexibility(S_cons.(groupnames{group}).S_cons(:,17:21)','temp');

block_flex(group,1)=mean(flex1);
block_flex(group,2)=mean(flex2);
block_flex(group,3)=mean(flex3);
block_flex(group,4)=mean(flex4);

end

figure;
for group=[3 1 2]
   plot(block_flex(group,:));
  hold on;  
end

