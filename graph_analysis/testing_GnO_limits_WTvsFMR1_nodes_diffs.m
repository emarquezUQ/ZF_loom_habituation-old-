
%%%%% this script is to explore the differences per individual nodes in
%%%%% the community detection measures. is a continuation of
%%%%% testing_GnO_limits_N_Qopt_FMR1.m and related to
%%%%% figures_n_stats_commun_measures.m

%% exploring the WTvsFMR1 diff per nodes


figure;
scatter3(flex_dif,cohe_dif,Prom_dif,10,'filled');

commMeasures=[];

commMeasures(:,1)=flex_dif;
commMeasures(:,2)=cohe_dif;
commMeasures(:,3)=Prom_dif;

STATS_FMR1_all.difs_measures=commMeasures;

%%%% raster plot of the differences sorted by cluster
figure;
imagesc(commMeasures);colormap(RdBu);caxis([-0.4 0.4]);%colorbar;

%%%% raster plot of the differences sorted by brain region
[B,I] = sort(NodesFmr1.Nod_brainID(keepFmr1));
figure;
imagesc(commMeasures(I,:));colormap(RdBu);caxis([-0.4 0.4]);%colorbar;



%%%% what if I cluster them?
eva = evalclusters(commMeasures,'kmeans','Silhouette','Klist',[1:8],'Distance','sqEuclidean');

figure;plot(eva)

kmeans_idx=kmeans(commMeasures,5,'Distance','sqEuclidean');

%%%% the pca didnt seem very useful
%  [coeff_Mes,score_Mes,~,~,explained_Mes,~] = pca(commMeasures);
%  figure;plot(cumsum(explained_Mes)); 

%%%% the tsne didnt seem very useful
Y=tsne(commMeasures,'Exaggeration',20);

figure;
gscatter(Y(:,1),Y(:,2),kmeans_idx);

%%%%% ploting identities of nodes 


figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
%   hold on;
   gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),kmeans_idx,'rygbm',[],21); 
% hold on;
% scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),30,cohe_dif,'filled');colormap(RdBu);caxis([-0.4 0.4]);colorbar;
view(-90,90);
%title('top 1SD dif flex');
hold off;

figure;
imagesc(kmeans_idx);

%%

figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
   hold on;
   gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
view(-90,90);
%title('top 1SD dif flex');
hold off;


figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
   hold on;
   gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1)); 
view(-90,90);
%title('top 1SD dif flex');
hold off;

figure;
imagesc(NodesFmr1.Nod_brainID(keepFmr1(I)));
figure;
imagesc(NodesFmr1.Nod_brainID(keepFmr1));
figure;
imagesc(NodesFmr1.Nod_clustID(keepFmr1));


%%

Flex_perNode=struct;
Flex_NodeTest=[];
for i=1:90
temp_groups=[];


counter=1;
for group=[3 2] %%% only controls and fmr1
    temp_groups(1:length(o),counter)=test.(groupnames{group})(:,i)';
counter=counter+1;
end


Flex_NodeTest(i,1)=median(temp_groups(:,1))-median(temp_groups(:,2));

%figure;boxplot(temp_groups);

Flex_perNode.(strcat('Node',num2str(i))).groups=temp_groups;

    [p, tbl, stats]=friedman(temp_groups,1,'off');
       
    Flex_perNode.(strcat('Node',num2str(i))).friedman.p=p;
    Flex_perNode.(strcat('Node',num2str(i))).friedman.p=tbl;
    Flex_perNode.(strcat('Node',num2str(i))).friedman.p=stats;
    
    %figure;
    [c, m]=multcompare(stats,'CType','bonferroni','display','off');

    Flex_perNode.(strcat('Node',num2str(i))).multcomp.c=c;
    Flex_perNode.(strcat('Node',num2str(i))).multcomp.m=m;
       
    Flex_NodeTest(i,2)=c(6);
    
    if c(6) < 0.05/90  %%%% with bonferroni correction.        
    Flex_NodeTest(i,3)=1;
    else
    Flex_NodeTest(i,3)=0;
    end
    
end


%%

Cohe_perNode=struct;
Cohe_NodeTest=[];
for i=1:90
temp_groups=[];


counter=1;
for group=[3 2] %%% only controls and fmr1
    temp_groups(1:length(o),counter)=test2.(groupnames{group})(:,i)';
counter=counter+1;
end


Cohe_NodeTest(i,1)=median(temp_groups(:,1))-median(temp_groups(:,2));

%figure;boxplot(temp_groups);

Cohe_perNode.(strcat('Node',num2str(i))).groups=temp_groups;

    [p, tbl, stats]=friedman(temp_groups,1,'off');
       
    Cohe_perNode.(strcat('Node',num2str(i))).friedman.p=p;
    Cohe_perNode.(strcat('Node',num2str(i))).friedman.p=tbl;
    Cohe_perNode.(strcat('Node',num2str(i))).friedman.p=stats;
    
    %figure;
    [c, m]=multcompare(stats,'CType','bonferroni','display','off');

    Cohe_perNode.(strcat('Node',num2str(i))).multcomp.c=c;
    Cohe_perNode.(strcat('Node',num2str(i))).multcomp.m=m;
       
    Cohe_NodeTest(i,2)=c(6);
    
    if c(6) < 0.05/90  %%%% with bonferroni correction.        
    Cohe_NodeTest(i,3)=1;
    else
    Cohe_NodeTest(i,3)=0;
    end
    
end


%%

Prom_perNode=struct;
Prom_NodeTest=[];
for i=1:90
temp_groups=[];


counter=1;
for group=[3 2] %%% only controls and fmr1
    temp_groups(1:length(o),counter)=test3.(groupnames{group})(:,i)';
counter=counter+1;
end


Prom_NodeTest(i,1)=median(temp_groups(:,1))-median(temp_groups(:,2));

%figure;boxplot(temp_groups);

Prom_perNode.(strcat('Node',num2str(i))).groups=temp_groups;

    [p, tbl, stats]=friedman(temp_groups,1,'off');
       
    Prom_perNode.(strcat('Node',num2str(i))).friedman.p=p;
    Prom_perNode.(strcat('Node',num2str(i))).friedman.p=tbl;
    Prom_perNode.(strcat('Node',num2str(i))).friedman.p=stats;
    
    %figure;
    [c, m]=multcompare(stats,'CType','bonferroni','display','off');

    Prom_perNode.(strcat('Node',num2str(i))).multcomp.c=c;
    Prom_perNode.(strcat('Node',num2str(i))).multcomp.m=m;
       
    Prom_NodeTest(i,2)=c(6);
    
    if c(6) < 0.05/90  %%%% with bonferroni correction.       
    Prom_NodeTest(i,3)=1;
    else
    Prom_NodeTest(i,3)=0;
    end
    
end

commMeasures2(:,1)=Flex_NodeTest(:,1);
commMeasures2(:,2)=Cohe_NodeTest(:,1); 
commMeasures2(:,3)=Prom_NodeTest(:,1);

figure;
imagesc(commMeasures2);colormap(RdBu);caxis([-0.4 0.4]);%colorbar;



commMeasuresPBonfe(:,1)=Flex_NodeTest(:,3);
commMeasuresPBonfe(:,2)=Cohe_NodeTest(:,3); 
commMeasuresPBonfe(:,3)=Prom_NodeTest(:,3);

figure;
imagesc(commMeasuresPBonfe);


%%

%%% flexibility, significant nodes
figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1(find(commMeasuresPBonfe(:,1))),1),NodesFmr1.Nod_coor(keepFmr1(find(commMeasuresPBonfe(:,1))),2),30,commMeasures2(find(commMeasuresPBonfe(:,1)),1),'filled');colormap(RdBu);caxis([-0.6 0.6]);colorbar;
view(-90,90);
%title('top 1SD dif flex');

%%% cohesion, significant nodes
figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1(find(commMeasuresPBonfe(:,2))),1),NodesFmr1.Nod_coor(keepFmr1(find(commMeasuresPBonfe(:,2))),2),30,commMeasures2(find(commMeasuresPBonfe(:,2)),2),'filled');colormap(RdBu);caxis([-0.6 0.6]);colorbar;
view(-90,90);
%title('top 1SD dif flex');

%%% promiscuity, significant nodes
figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1(find(commMeasuresPBonfe(:,3))),1),NodesFmr1.Nod_coor(keepFmr1(find(commMeasuresPBonfe(:,3))),2),30,commMeasures2(find(commMeasuresPBonfe(:,3)),3),'filled');colormap(RdBu);caxis([-0.2 0.2]);colorbar;
view(-90,90);
%title('top 1SD dif flex');


