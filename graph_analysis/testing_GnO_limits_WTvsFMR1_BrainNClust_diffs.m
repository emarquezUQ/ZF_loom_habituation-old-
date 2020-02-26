
%%%%% this script is to explore the differences per brain region and per 
%%%%% cluster type in the community detection measures. is a continuation of
%%%%% testing_GnO_limits_N_Qopt_FMR1.m and related to
%%%%% figures_n_stats_commun_measures.m

STATS_FMR1_BrainNClust=struct;

%% getting the combinations

BrainNClustTest={};
%Flex_perBrainNClust=struct;

testNum=1;
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    temp_idx1=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);
    for clust=1:4 
    temp_idx2=find(ClustID_CL4(keepFmr1)==clust);
    temp_idx=intersect(temp_idx2,temp_idx1);
    
    if isempty(temp_idx)
        continue 
        
    else
    
    BrainNClustTest{testNum}=strcat(RegionList{brain},'/',clustnames_CL4{clust});
    
    testNum=testNum+1;
    
    end
    end
end


%% for flexibility. 

Flex_BrainNClustTest=[];
%Flex_perBrainNClust=struct;

testNum=1;
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    temp_idx1=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);
    for clust=1:4 
    temp_idx2=find(ClustID_CL4(keepFmr1)==clust);
    temp_idx=intersect(temp_idx2,temp_idx1);
    
    if isempty(temp_idx)
        continue 
        
    else
    
    temp_groups=[];
       
        counter=1;
        for group=[3 2] %%% only WT vs fmr1
            temp_groups(1:length(o),counter)=mean(test.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    %Flex_perBrainNClust.(RegionList{brain}).(clustnames_CL4{clust})=temp_groups;
    
    STATS_FMR1_BrainNClust.Flex_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).groups=temp_groups;
    
    %figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups,1,'off');
    
    STATS_FMR1_BrainNClust.Flex_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).friedman.p=p;
    STATS_FMR1_BrainNClust.Flex_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).friedman.tbl=tbl;
    STATS_FMR1_BrainNClust.Flex_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).friedman.stats=stats;
    
    %figure;title(RegionList{brain});
    [c, m]=multcompare(stats,'CType','bonferroni','display','off');
        
    STATS_FMR1_BrainNClust.Flex_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).multcomp.c=c;
    STATS_FMR1_BrainNClust.Flex_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).multcomp.m=m;
    
    
    Flex_BrainNClustTest(1,testNum)=median(temp_groups(:,1))-median(temp_groups(:,2));
    Flex_BrainNClustTest(2,testNum)=c(6);
    if c(6) < 0.05/26  %%%% with bonferroni correction.        
    Flex_BrainNClustTest(3,testNum)=1;
    else
    Flex_BrainNClustTest(3,testNum)=0;
    end
    
    testNum=testNum+1;
    
    end
    end
end

%% for cohesion



Cohe_BrainNClustTest=[];
%Flex_perBrainNClust=struct;

testNum=1;
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    temp_idx1=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);
    for clust=1:4 
    temp_idx2=find(ClustID_CL4(keepFmr1)==clust);
    temp_idx=intersect(temp_idx2,temp_idx1);
    
    if isempty(temp_idx)
        continue 
        
    else
    
    temp_groups=[];
       
        counter=1;
        for group=[3 2] %%% only WT vs fmr1
            temp_groups(1:length(o),counter)=mean(test2.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    %Flex_perBrainNClust.(RegionList{brain}).(clustnames_CL4{clust})=temp_groups;
    
    STATS_FMR1_BrainNClust.Cohe_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).groups=temp_groups;
    
    %figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups,1,'off');
    
    STATS_FMR1_BrainNClust.Cohe_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).friedman.p=p;
    STATS_FMR1_BrainNClust.Cohe_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).friedman.tbl=tbl;
    STATS_FMR1_BrainNClust.Cohe_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).friedman.stats=stats;
    
    %figure;title(RegionList{brain});
    [c, m]=multcompare(stats,'CType','bonferroni','display','off');
        
    STATS_FMR1_BrainNClust.Cohe_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).multcomp.c=c;
    STATS_FMR1_BrainNClust.Cohe_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).multcomp.m=m;
    
    
    Cohe_BrainNClustTest(1,testNum)=median(temp_groups(:,1))-median(temp_groups(:,2));
    Cohe_BrainNClustTest(2,testNum)=c(6);
    if c(6) < 0.05/26  %%%% with bonferroni correction.        
    Cohe_BrainNClustTest(3,testNum)=1;
    else
    Cohe_BrainNClustTest(3,testNum)=0;
    end
    
    testNum=testNum+1;
    
    end
    end
end


%% for promiscuity


Prom_BrainNClustTest=[];
%Flex_perBrainNClust=struct;

testNum=1;
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    temp_idx1=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);
    for clust=1:4 
    temp_idx2=find(ClustID_CL4(keepFmr1)==clust);
    temp_idx=intersect(temp_idx2,temp_idx1);
    
    if isempty(temp_idx)
        continue 
        
    else
    
    temp_groups=[];
       
        counter=1;
        for group=[3 2] %%% only WT vs fmr1
            temp_groups(1:length(o),counter)=mean(test3.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    %Flex_perBrainNClust.(RegionList{brain}).(clustnames_CL4{clust})=temp_groups;
    
    STATS_FMR1_BrainNClust.Prom_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).groups=temp_groups;
    
    %figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups,1,'off');
    
    STATS_FMR1_BrainNClust.Prom_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).friedman.p=p;
    STATS_FMR1_BrainNClust.Prom_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).friedman.tbl=tbl;
    STATS_FMR1_BrainNClust.Prom_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).friedman.stats=stats;
    
    %figure;title(RegionList{brain});
    [c, m]=multcompare(stats,'CType','bonferroni','display','off');
        
    STATS_FMR1_BrainNClust.Prom_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).multcomp.c=c;
    STATS_FMR1_BrainNClust.Prom_perBrain.(RegionList{brain}).(clustnames_CL4{clust}).multcomp.m=m;
    
    
    Prom_BrainNClustTest(1,testNum)=median(temp_groups(:,1))-median(temp_groups(:,2));
    Prom_BrainNClustTest(2,testNum)=c(6);
    if c(6) < 0.05/26  %%%% with bonferroni correction.        
    Prom_BrainNClustTest(3,testNum)=1;
    else
    Prom_BrainNClustTest(3,testNum)=0;
    end
    
    testNum=testNum+1;
    
    end
    end
end



