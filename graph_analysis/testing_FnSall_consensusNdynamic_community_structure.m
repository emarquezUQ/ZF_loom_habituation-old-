
%%%%% this script is to thest with a specific gamma and omega combination
%%%%% the dynamic community structure in the FnS dataset

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

%% 

S_cons=struct;

for data=1:4

A={};

for loom=1:31
   temp=Data_corrMat2.(datasets(data,:)).Mean_corrMat{1,loom}(keep,keep);
   temp(isnan(temp))=0;
    A{loom}=temp;  
end
clear temp


%%%% making a multidimensional structure where to store things
S_test=[];

for test=1:1000
     
gamma = 0.9;
omega = 0.9; %%% this makes a big influence... 

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
for i=1:31
   C=S_test(:,i,:); 
   C=squeeze(C);
   [S2, Q2, X_new3, qpc] = consensus_iterative(C');
    S_good(:,i)=S2(i,:);
end


S_cons.(datasets(data,:)).S_test=S_test;
S_cons.(datasets(data,:)).S_cons=S_good;


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
    scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,S_good(:,time(i)),'filled');colormap('jet');caxis([low high]);
    view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  
     
  end
  subplot(1,7,7);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
   gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],20,'off');     
     view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  
    
end

%%



for data=1:4

flex=flexibility(S_cons.(datasets(data,:)).S_cons','temp'); %%% i need to mind the orientation of the matrix to do it properly

S_cons.(datasets(data,:)).flex=flex;

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),Nodes.Nod_clustID(keep,:),'rgggbm',[],21); 
hold on;
scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,flex,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title('flexibility');
hold off;


clustnames={'strhab1','strhab2','strhab3','modhab','weakhab','inhib'};

    
figure;
subplot(1,2,1);boxplot(flex,Nodes.Nod_brainID(keep)),xticklabels(RegionList);xtickangle(45);ylim([0 1]);
subplot(1,2,2);boxplot(flex,ClustID),xticklabels(clustnames);xtickangle(45);ylim([0 1]);

end

%% cohesion strength and related

options.figureFlag	= 0;
options.colormap	= 'jet';

for data=1:4

[Cij,node_cohesion,node_disjoint,node_flexibility,strength_cohesion,commChanges,commCohesion,commDisjoint,commIndex] = calc_node_cohesion(S_cons.(datasets(data,:)).S_cons,options);

S_cons.(datasets(data,:)).node_cohesion=node_cohesion;

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
% hold on;
% gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),Nodes.Nod_clustID(keep,:),'gggbrm',[],21); 
hold on;
scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,node_cohesion,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title('cohesion');
hold off;


figure;
subplot(1,2,1);boxplot(node_cohesion,Nodes.Nod_brainID(keep)),xticklabels(RegionList);xtickangle(45);ylim([0 1]);
subplot(1,2,2);boxplot(node_cohesion,ClustID),xticklabels(clustnames);xtickangle(45);ylim([0 1]);

end
%%

for data=1:4

P = promiscuity(S_cons.(datasets(data,:)).S_cons');  %%% i need to mind the orientation of the matrix to do it properly

S_cons.(datasets(data,:)).P=P;

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
% hold on;
% gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),Nodes.Nod_clustID(keep,:),'gggbrm',[],21); 
hold on;
scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,P,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title('promiscuity');
hold off;
  
figure;
subplot(1,2,1);boxplot(P,Nodes.Nod_brainID(keep)),xticklabels(RegionList);xtickangle(45);ylim([0 1]);
subplot(1,2,2);boxplot(P,ClustID),xticklabels(clustnames);xtickangle(45);ylim([0 1]);

end



save('S_cons_FnS_All_G09O09.mat','S_cons');
%%


counter=1;
figure;
for data=1:4
subplot(3,4,counter);

plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
hold on;
scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,S_cons.(datasets(data,:)).flex,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title(strcat('flexibility','/',datasets(data,:)));
hold off;

counter=counter+1;
end
for data=1:4
subplot(3,4,counter);

plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
hold on;
scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,S_cons.(datasets(data,:)).node_cohesion,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title(strcat('cohesion','/',datasets(data,:)));
hold off;

counter=counter+1;
end
for data=1:4
subplot(3,4,counter);

plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
hold on;
scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,S_cons.(datasets(data,:)).P,'filled'); colormap('cool'); caxis([0 1]);colorbar;
view(-90,90);
title(strcat('promiscuity','/',datasets(data,:)));
hold off;

counter=counter+1;
end


%%


%% comparing groups with stats


groups_flexibility=[];
counter=1;
for data=1:4

    groups_flexibility(:,counter)=S_cons.(datasets(data,:)).flex;

counter=counter+1;
end

groups_cohesion=[];
counter=1;
for data=1:4

    groups_cohesion(:,counter)=S_cons.(datasets(data,:)).node_cohesion;

counter=counter+1;
end

groups_prom=[];
counter=1;
for data=1:4

    groups_prom(:,counter)=S_cons.(datasets(data,:)).P;

counter=counter+1;
end


figure;
subplot(1,3,1);boxplot(groups_flexibility);title('flexibility');ylim([0 1]);
subplot(1,3,2);boxplot(groups_cohesion);title('cohesion');ylim([0 1]);
subplot(1,3,3);boxplot(groups_prom);title('promiscuity');ylim([0 1]);


[p, tbl, stats]=anova1(groups_flexibility);
figure;
[c, m]=multcompare(stats,'CType','bonferroni');


[p, tbl, stats]=anova1(groups_cohesion);
figure;
[c, m]=multcompare(stats,'CType','bonferroni');


[p, tbl, stats]=anova1(groups_prom);
figure;
[c, m]=multcompare(stats,'CType','bonferroni');


[p, tbl, stats]=anova1(groups_prom);
figure;
[c2, m2]=multcompare(stats); %%% to check if it is different without bonferroni. 

