


%%%%% this script is to see if there is a trend in flexibility, cohesion
%%%%% and promiscuity. 


%%%%% I will expand the range of gammas (2.5) and omegas (2) but reduce the
%%%%% repetitions of the community detection to a 100. dani said it was
%%%%% fine

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

Big_FMR1_OPT={};
FlexMat=struct;
CoheMat=struct;
PromMat=struct;


%%
Allgamma = 0.1:0.1:2.5;
Allomega = 0.1:0.1:2;


%for g=1%:16
for g=1:25
for o=1:20

gamma = Allgamma(g);
omega = Allomega(o);   
    
%% making 100 iterations per group to calculate a representative modularity 
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

for test=1:100
     
%gamma = 0.9;
%omega = 0.7; %%% this makes a big influence... 

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



%% testing flexibility

flex=flexibility(S_cons.(groupnames{group}).S_cons','temp'); %%% i need to mind the orientation of the matrix to do it properly

S_cons.(groupnames{group}).flex=flex;

%% cohesion strength and related

options.figureFlag	= 0;
options.colormap	= 'jet';

[Cij,node_cohesion,node_disjoint,node_flexibility,strength_cohesion,commChanges,commCohesion,commDisjoint,commIndex] = calc_node_cohesion(S_cons.(groupnames{group}).S_cons,options);

S_cons.(groupnames{group}).node_cohesion=node_cohesion;

%% promiscuity

P = promiscuity(S_cons.(groupnames{group}).S_cons');  %%% i need to mind the orientation of the matrix to do it properly

S_cons.(groupnames{group}).P=P;

end

%%
groups_flexibility=[];
counter=1;
for group=[3 1 2]

    groups_flexibility(:,counter)=S_cons.(groupnames{group}).flex;

counter=counter+1;
end

groups_cohesion=[];
counter=1;
for group=[3 1 2]

    groups_cohesion(:,counter)=S_cons.(groupnames{group}).node_cohesion;

counter=counter+1;
end

groups_prom=[];
counter=1;
for group=[3 1 2]

    groups_prom(:,counter)=S_cons.(groupnames{group}).P;

counter=counter+1;
end


%% relative change


%%%%%% flexibility

flex_dif=S_cons.(groupnames{3}).flex-S_cons.(groupnames{2}).flex;

%%%%%% cohesion

cohe_dif=S_cons.(groupnames{3}).node_cohesion-S_cons.(groupnames{2}).node_cohesion;

%%%%%% promiscuity
P_dif=S_cons.(groupnames{3}).P-S_cons.(groupnames{2}).P;

%%%%% 

FlexMat.(groupnames{3})(g,o)=mean(S_cons.(groupnames{3}).flex);
FlexMat.(groupnames{1})(g,o)=mean(S_cons.(groupnames{1}).flex);
FlexMat.(groupnames{2})(g,o)=mean(S_cons.(groupnames{2}).flex);

CoheMat.(groupnames{3})(g,o)=mean(S_cons.(groupnames{3}).node_cohesion);
CoheMat.(groupnames{1})(g,o)=mean(S_cons.(groupnames{1}).node_cohesion);
CoheMat.(groupnames{2})(g,o)=mean(S_cons.(groupnames{2}).node_cohesion);

PromMat.(groupnames{3})(g,o)=mean(S_cons.(groupnames{3}).P);
PromMat.(groupnames{1})(g,o)=mean(S_cons.(groupnames{1}).P);
PromMat.(groupnames{2})(g,o)=mean(S_cons.(groupnames{2}).P);


%%
name=strcat('g',num2str(10*gamma),'_o',num2str(10*omega));

Big_FMR1_OPT{g,o}=S_cons;


end
end


%% saving

save('S_cons_testing_GnO_measures6.mat','Big_FMR1_OPT','FlexMat','CoheMat','PromMat','Allgamma','Allomega');



%%

counter=1;
figure;
for group=[3 1 2]
   subplot(3,3,counter);imagesc(FlexMat.(groupnames{group}));caxis([0 1]);
   title(strcat('flex/',groupnames{group}));
   
   counter=counter+1;
end
for group=[3 1 2]
   subplot(3,3,counter);imagesc(CoheMat.(groupnames{group}));caxis([0 1]);
   title(strcat('cohe/',groupnames{group}));
   
   counter=counter+1;
end
for group=[3 1 2]
   subplot(3,3,counter);imagesc(PromMat.(groupnames{group}));caxis([0 1]);
   title(strcat('prom/',groupnames{group}));
   
   counter=counter+1;
end


%%


counter=1;
test=[];
for group=[3 1 2]
    test(:,counter)=FlexMat.(groupnames{group})(:);
    
    counter=counter+1;
end


figure;boxplot(test);

    [p, tbl, stats]=anova1(test);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');
    
   %%%%OR
   for group=[1 2 3]  %%% cause now they are in the right order
   kstest(test(:,group))
   end 
   
   %%%% they are not normal... so  kruskalwallis? or Friedman? I think it
   %%%% would be friedman cause is not replicates on the same conditions
        %[p, tbl, stats]=kruskalwallis(test);
        [p, tbl, stats]=friedman(test,1);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');

%%


%%%%%% NOTE: I read that the multiple comparison test in R for friedman
%%%%%% tests are not accurate because they are sensitive to the order of
%%%%%% the variables... that is better to test this in R. I tried and I had
%%%%%% very similar results. 

%%%% comparing flexibility
counter=1;
test=struct;
for g=1:25
for o=1:20
    
    for group=1:3   
    
    test.(groupnames{group})(counter,:)=Big_FMR1_OPT{g,o}.(groupnames{group}).flex
    end
    counter=counter+1;
end
end


groups_flexibilityTest=NaN(90,3);
counter=1;
for group=[3 1 2]
    groups_flexibilityTest(1:90,counter)=mean(test.(groupnames{group}));
counter=counter+1;
end

figure;boxplot(groups_flexibilityTest);

    [p, tbl, stats]=kruskalwallis(groups_flexibilityTest);
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
    
    [p, tbl, stats]=kruskalwallis(temp_groups);
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
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),30,flex_dif,'filled');colormap('parula');%caxis([-0.6 0.6]);%colorbar;
view(-90,90);
%title('top 1SD dif flex');
hold off;

%%

%%%% now cohesion 
counter=1;
test2=struct;
for g=1:16
for o=1:10
    
    for group=1:3   
    
    test2.(groupnames{group})(counter,:)=Big_FMR1_OPT{g,o}.(groupnames{group}).node_cohesion
    end
    counter=counter+1;
end
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
test2=struct;
for g=1:16
for o=1:10
    
    for group=1:3   
    
    test2.(groupnames{group})(counter,:)=Big_FMR1_OPT{g,o}.(groupnames{group}).P
    end
    counter=counter+1;
end
end


groups_PromTest=NaN(90,3);
counter=1;
for group=[3 1 2]
    groups_PromTest(1:90,counter)=mean(test2.(groupnames{group}));
counter=counter+1;
end

figure;boxplot(groups_PromTest);

    [p, tbl, stats]=anova1(groups_PromTest);
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');

    
%%

x=mean(FlexMat.control(:));
y=std(FlexMat.control(:));

temp_idx=find(FlexMat.control<x+y & FlexMat.control>x-y);

Test=zeros(size(FlexMat.control));
Test(temp_idx)=1;
figure;imagesc(Test);



%%

figure;
for group=[3 1 2]
plot(mean(FlexMat.(groupnames{group})));
hold on;
end


