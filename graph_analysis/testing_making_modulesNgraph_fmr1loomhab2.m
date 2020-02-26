



%%%%% this script is to try to further analyze the graph fmr1 data using
%%%%% the Brain Connectivity Toolbox. this is the second attempt, now with
%%%%% the nodes also sorted by brain regions. 

%%% i need to have it on the path. 




%% checking the corr. distributions
%%%% first just having a look at the corr coeficient values

for g=1:3
    figure;
count=1;
    count=1;
    group=groupnames{g,1};
   sgtitle(group) 
   
    for k=[2 3 4 5 6 11 12]
    subplot(1,7,count)
        
    R=Data_corrMat2.(group).Mean_corrMat{1,k};
        

    [~,~,weights] = find(tril(R,-1));
    %quantile(weights,[0.025 0.25 0.50 0.75 0.975])
    
    count=count+1;
    histogram(weights,20)
    end
    
end
%%



%%% setting the diagonal to 0. 
%%% looms 1,2,3,4,5,10 and 11
MatAll=struct;
for g=1:3
    group=groupnames{g,1};
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,4,5,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat2.(group).Mean_corrMat{1,moment(m)}),0.75);
     MatAll.(group).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end




%% density
%


%kden = density_und(CIJ);


for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)

kden=density_und(MatAll.(group).(loom{i}).Mat);

MatAll.(group).(loom{i}).kden=kden;
end
end

figure;
temp=[];
for g=[3 2]  %%%% only controls and fmr1 and in that order. 
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)
    
    temp(1,i)=MatAll.(group).(loom{i}).kden;
end
plot(temp);
hold on;

end


%%
%%%% I will look at the number of degrees and strength.
%%%% I could maybe use them to see how they decrease through habituation or
%%%% increase at revovery


%%%% degrees and strengths. I am doing absolute ones. very similar results than degrees. 

for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)

deg=degrees_und(MatAll.(group).(loom{i}).Mat);
str=strengths_und_sign(abs(MatAll.(group).(loom{i}).Mat));

MatAll.(group).(loom{i}).deg=deg;
MatAll.(group).(loom{i}).str=str;
end
end

%%%% to plot histograms
figure;
for i=1:7
subplot(1,7,i);
for g=1:3
group=groupnames{g,1};
histogram(MatAll.(group).(loom{i}).deg,20);

hold on;
end
end

%%%% non-parametrical test to check if they are different. 
h = ranksum(MatAll.control.loom1.deg,MatAll.fmr1.loom1.deg) %%% it seems that the differences for the first loom are sig. but i should test for normality
h = ranksum(MatAll.control.loom2.deg,MatAll.fmr1.loom2.deg) %%% nonsig
h = ranksum(MatAll.control.loom3.deg,MatAll.fmr1.loom3.deg) %%% nonsig
h = ranksum(MatAll.control.loom10.deg,MatAll.fmr1.loom10.deg) %%% sig
h = ranksum(MatAll.control.loom11.deg,MatAll.fmr1.loom11.deg) %%% nonsig


%%%% to plot with brains
figure;
counter=1;

for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)
    if counter==1|counter==8|counter==15|counter==22
    low=0;high=100;
    else
    low=0;high=50; 
    end
    
  subplot(3,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),25,MatAll.(group).(loom{i}).deg,'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 title(group);
counter=counter+1;
end
end

%%%% to substract the controls minus fmr1
figure;
counter=1;
for g=3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)
    
    low=-30;high=30;
    
subplot(1,7,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); 
hold on; scatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),25,(MatAll.(group).(loom{i}).deg-MatAll.(groupnames{g-1,1}).(loom{i}).deg),'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end
end

%% betweeness centrality. 
%%% i need to first to a connection-length matrix
%%% I tried and the results look kind of weird... not sure what to make of
%%% them. so i drop them. 



%% participation and gateway coef

%P = participation_coef(X,modulesID);

%P=participation_coef(Matctrl_1,Nodes2.Mod_clust);
%[Gpos,~]=gateway_coef_sign(Matctrl_1,Nodes2.Mod_clust,1);

for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)

temp_mat=MatAll.(group).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes2.Mod_clust);
[Gpos,~]=gateway_coef_sign(MatAll.(group).(loom{i}).Mat,Nodes2.Mod_clust,1);

MatAll.(group).(loom{i}).P=P;
MatAll.(group).(loom{i}).Gpos=Gpos;
end
end

%%% to plot it on brains
figure;
counter=1;

for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)
    
    low=0;high=0.8;
     
  subplot(3,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); title(group);
  hold on; scatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),25,MatAll.(group).(loom{i}).P,'filled');colormap('jet');view(-90,90);%caxis([low high]);%colorbar; %
 
counter=counter+1;
end
end


%%%% to substract the controls minus fmr1
figure;
counter=1;
for g=3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)
    
    low=-1;high=1;
    
subplot(1,7,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); 
hold on; scatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),25,(MatAll.(group).(loom{i}).P-MatAll.(groupnames{g-1,1}).(loom{i}).P),'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end
end




%% core vs periphery
% in the function guidelines it says that its supposed to have directed values...
% so i dont know if we can use this, i tried it anyway. 


%CoreP     = core_periphery_dir(W);

for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)
temp_mat=MatAll.(group).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
CoreP=core_periphery_dir(temp_mat);

MatAll.(group).(loom{i}).CoreP=CoreP;
end
end

%%% to plot it on brains
figure;
counter=1;

for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)
    
    low=0;high=1;
     
  subplot(3,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); title(group);
  hold on; scatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),25,MatAll.(group).(loom{i}).CoreP,'filled');colormap('jet');view(-90,90);%caxis([low high]);%colorbar; %
 
counter=counter+1;
end
end

%%%% to substract the controls minus fmr1
figure;
counter=1;
for g=3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)
    
    low=-1;high=1;
    
subplot(1,7,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); 
hold on; scatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),25,(MatAll.(group).(loom{i}).CoreP-MatAll.(groupnames{g-1,1}).(loom{i}).CoreP),'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end
end


%% clustering coef
% i tried it... it might be interesting to see how the connections are lost
% during habituation but I dont think it gives us much more information. 
%%% I could calculate averages and compare between ctrl and fmr1

%Ccoef = clustering_coef_wu(W);
for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)
temp_mat=MatAll.(group).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
Ccoef=clustering_coef_wu(temp_mat);

MatAll.(group).(loom{i}).Ccoef=Ccoef;
end
end

%%% to test particular combinations. 
data=1;
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));
i=1
h = ranksum(MatAll.control.(loom{i}).Ccoef,MatAll.fmr1.(loom{i}).Ccoef) %%% 

%%% to plot it on brains
figure;
counter=1;

for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)
    
    low=0;high=1;
     
  subplot(3,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);title(group);
  hold on; scatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),25,MatAll.(group).(loom{i}).Ccoef,'filled');colormap('jet');view(-90,90);%caxis([low high]);%colorbar; %
 
counter=counter+1;
end
end



%% characteristic path length
%  i first need to get a connection-length matrix, then a distance matrix. I dont think the charpath will
%  really give me something interesting but it can be used to calculate teh
%  small-worldness. I could also compare it between datasets


for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)
temp_mat=MatAll.(group).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
LMat = weight_conversion(temp_mat, 'lengths');
DMat = distance_wei(LMat);
[CharP,~,~,~,~]=charpath(DMat,0,0);

MatAll.(group).(loom{i}).LMat=LMat;
MatAll.(group).(loom{i}).DMat=DMat;
MatAll.(group).(loom{i}).CharP=CharP;
end
end

%%
%%% calculating small-worldness as burgstaller et al 2019 bioRxiv did
%%% although I think I am doing this wrong... I think I should calculate
%%% all the values per fish. 

%%%% so far I can see some differences between genotypes and how the
%%%% differences drop with the stimulus but I am not really sure what to do
%%%% with it. 

for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)

    SW=mean(MatAll.(group).(loom{i}).Ccoef)/MatAll.(group).(loom{i}).CharP;

MatAll.(group).(loom{i}).SW=SW;

end
end

figure;
temp=[];
for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll.(group));

for i=1:length(loom)
    
    temp(1,i)=MatAll.(group).(loom{i}).SW;
end
plot(temp);
hold on;

end

%% modularity and modules based on correlation
%%% it gives some intresting results but they are not the same as the
%%% HierarchicalConsensus method. it is much fasther though. 

%Ci = modularity_und(W);
for g=1:3
    
    group=groupnames{g,1}; 
   
   loom=fieldnames(MatAll.(group));

for i=1:length(loom)

    temp_mat=MatAll.(group).(loom{i}).Mat;
    temp_mat(isnan(temp_mat))=0;
    
mod=modularity_und(temp_mat);

MatAll.(group).(loom{i}).mod=mod;
end
end

for g=1:3
    figure;
    group=groupnames{g,1}; 
   
   loom=fieldnames(MatAll.(group));
for i=1:length(loom)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); hold on; gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),MatAll.(group).(loom{i}).mod);view(-90,90);  %colorbar;colormap('jet');caxis([-1 1]);
title(strcat(groupnames{g,1},'_',loom{i}));

pause(3);

end
end

 figure;
 counter=1;
for g=1:3
   
    group=groupnames{g,1}; 
   
   loom=fieldnames(MatAll.(group));
for i=[1 2 3 6 7]
 subplot(3,5,counter);   
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); hold on; gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),MatAll.(group).(loom{i}).mod);view(-90,90);  %colorbar;colormap('jet');caxis([-1 1]);
title(strcat(groupnames{g,1},'_',loom{i}));

counter=counter+1;
end
end

%%%%% HierarchicalConsensus. Based on Betzel 2018 bioRxiv

%%%%%%%%%%%%  NOTE

%%%% I havent tested this yet.

%%% this gave me good results. but... see below
[gamma_min,gamma_max]=gammaRange(Data_corrMat2.fmr1.Mean_corrMat{1,2});
S = fixedResSamples(Data_corrMat2.fmr1.Mean_corrMat{1,2}, 10000, 'Gamma', 1.2);
[Sc, Tree] = hierarchicalConsensus(S,0.5);
[Sall, thresholds] = allPartitions(Sc, Tree);
drawHierarchy(Sc, Tree)
C = coclassificationMatrix(S);
figure;consensusPlot(C, Sc, Tree)


%%% this onw is not bad at all. it is actually very good! many interesting
%%% connections! I tried with different sizes of ensamble and it seems that
%%% I always got the same results from 1000000 to 150000. but not with
%%% 100000. the maximum number of modules I got was 9.
S = eventSamples(Matfmr1_1, 1000000); %%%% it seems that it works too with 500000, 250000 and 150000. It doesnt work with 100000.
[gamma_min,gamma_max]=gammaRange(Matfmr1_1);
%S = fixedResSamples(Data_corrMat2.fmr1.Mean_corrMat{1,2}, 10000, 'Gamma', 1.2);
%S = exponentialSamples(Matfmr1_1, 10000);
[Sc, Tree] = hierarchicalConsensus(S,0.05);
[Sall, thresholds] = allPartitions(Sc, Tree);
drawHierarchy(Sc, Tree)
C = coclassificationMatrix(S);
consensusPlot(C, Sc, Tree)

unique(Sall(:,1))
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); hold on; gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Sall(:,1));view(-90,90);  %colorbar;colormap('jet');caxis([-1 1]);
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); hold on; gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Sall(:,2));view(-90,90);  %colorbar;colormap('jet');caxis([-1 1]);
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); hold on; gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_clust);view(-90,90);  %colorbar;colormap('jet');caxis([-1 1]);

%%% for fmr1
S4 = eventSamples(Matfmr1_1, 150000);
[gamma_min,gamma_max]=gammaRange(Matfmr1_1);
%S = fixedResSamples(Data_corrMat2.fmr1.Mean_corrMat{1,2}, 10000, 'Gamma', 1.2);
%S = exponentialSamples(Matfmr1_1, 10000);
[Sc4, Tree4] = hierarchicalConsensus(S4,0.05);
[Sall4, thresholds4] = allPartitions(Sc4, Tree4);
%drawHierarchy(Sc4, Tree4)
C4 = coclassificationMatrix(S4);
figure;consensusPlot(C4, Sc4, Tree4)

figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); hold on; gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Sall4(:,1));view(-90,90);  %colorbar;colormap('jet');caxis([-1 1]);


%%% now with ctrls. they have very different gamma range...
S2 = eventSamples(Matctrl_1, 150000);  %%%% I tried with 150000 and 500000 and got the same results. 
[gamma_min,gamma_max]=gammaRange(Matctrl_1);
%S = fixedResSamples(Data_corrMat2.fmr1.Mean_corrMat{1,2}, 10000, 'Gamma', 1.2);
%S = exponentialSamples(Matfmr1_1, 10000);
[Sc2, Tree2] = hierarchicalConsensus(S2,0.05);
[Sall2, thresholds2] = allPartitions(Sc2, Tree2);
%drawHierarchy(Sc2, Tree2)
C2 = coclassificationMatrix(S2);
figure;consensusPlot(C2, Sc2, Tree2)

figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); hold on; gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Sall2(:,1));view(-90,90);  %colorbar;colormap('jet');caxis([-1 1]);




%%% for all of them. 

%%%%% NOTE: something wrong happend with this loop. it didnt work well. i got the
%%%%% same results for all the genotypes and looms. 
%%%%% I found the error. I forgot to change the name of "Sc" so i got the same one every time. I
%%%%% fixed it. And I corrected it with a loop bellow. 


for g=1:3
    
    group=groupnames{g,1}; 
   
   loom=fieldnames(MatAll.(group));

for i=1:length(loom)
  
    temp_mat=MatAll.(group).(loom{i}).Mat;
    temp_mat(isnan(temp_mat))=0;

    S = eventSamples(temp_mat, 150000); % 150000 seems to work fine
    %[gamma_min,gamma_max]=gammaRange(temp_mat);
    %S = fixedResSamples(temp_mat, 10000, 'Gamma', 1.2);
    %S = exponentialSamples(temp_mat, 10000);
    [Sc, Tree] = hierarchicalConsensus(S,0.05);
    [Sall, thresholds] = allPartitions(Sc, Tree);
    %drawHierarchy(Sc2, Tree2)
    C = coclassificationMatrix(S);
    % figure;consensusPlot(C2, Sc2, Tree2)

    % figure;
    % plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
    % hold on;
    % gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Sall2(:,1),[],'.',25);view(-90,90);
    
    
    
MatAll.(group).(loom{i}).S=S;
MatAll.(group).(loom{i}).Sc=Sc;
MatAll.(group).(loom{i}).Tree=Tree;
MatAll.(group).(loom{i}).Sall=Sall;
MatAll.(group).(loom{i}).thresholds=thresholds;
MatAll.(group).(loom{i}).C=C;

end
end


for g=1:3
    
    group=groupnames{g,1}; 
   
   loom=fieldnames(MatAll.(group));

for i=1:length(loom)
     
    [Sall, thresholds] = allPartitions(MatAll.(group).(loom{i}).Sc, MatAll.(group).(loom{i}).Tree);
        
MatAll.(group).(loom{i}).Sall=Sall;
MatAll.(group).(loom{i}).thresholds=thresholds;


end
end



%%% to plot it. 

for g=1:3
 figure;   
    group=groupnames{g,1}; 
loom=fieldnames(MatAll.(group));
for i=1:length(loom)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); hold on; gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),MatAll.(group).(loom{i}).Sall(:,1));view(-90,90);  %colorbar;colormap('jet');caxis([-1 1]);
title(strcat(groupnames{g,1},'_',loom{i}));

pause(3);

end
end


figure;
 counter=1;
for g=1:3
   
    group=groupnames{g,1}; 
   
   loom=fieldnames(MatAll.(group));
for i=[1 2 3 6 7]
 subplot(3,5,counter);   
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); hold on; gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),MatAll.(group).(loom{i}).Sall(:,1));view(-90,90);  %colorbar;colormap('jet');caxis([-1 1]);
title(strcat(groupnames{g,1},'_',loom{i}));

counter=counter+1;
end
end


save('graph_analysis_fmr1loomhab.mat', 'MatAll');

%% trying a substraction of the graphs

%%%%% this script is to try to do a substracted graph for the fmr1loomhab
%%%%% dataset. for looms 1-5, 10 and 11th. The idea is to colorcode the
%%%%% connections based on the difference in strenght (correlation). 




%%  first substracting the matrices. 


ctrl_fmr1_subs_mat=struct;   
    
for k=[1:21]
  
    temp_mat1=Data_corrMat2.control.Mean_corrMat{1,k};
    temp_mat1(isnan(temp_mat1))=0;
    
    temp_mat2=Data_corrMat2.fmr1.Mean_corrMat{1,k};
    temp_mat2(isnan(temp_mat2))=0;
    
    temp_mat=temp_mat1-temp_mat2;
    
    ctrl_fmr1_subs_mat.Subs_Mean_corrMat{1,k}=temp_mat;
    
end


%%
%%%% here i am trying to get the connections from my previous filtering
  
for k=[1:21]
  
    temp_mat1=Data_corrMat2.control.Mean_corrMat{1,k};
    temp_mat1(isnan(temp_mat1))=0;
    
    temp_mat1_idx=find(threshold_absolute(abs(temp_mat1),0.75));
    
    temp_mat2=Data_corrMat2.fmr1.Mean_corrMat{1,k};
    temp_mat2(isnan(temp_mat2))=0;
    
    temp_mat2_idx=find(threshold_absolute(abs(temp_mat2),0.75));
    
    temp_mat=temp_mat1-temp_mat2;
    
    temp_mat_idx=union(temp_mat1_idx,temp_mat2_idx);
    
    temp_mat_cleaned=zeros(size(temp_mat));
    temp_mat_cleaned(temp_mat_idx)=temp_mat(temp_mat_idx);
    
    ctrl_fmr1_subs_mat.Subs_Mean_corrMat_cleaned{1,k}=temp_mat_cleaned;
    
end
    
    %% now to plot some graphs. 
   %figure;
   %count=1;
    for k=[2 3 4 5 6 11 12]
    figure;
    set(gcf, 'Position',  [100, 100, 700, 900])
    
    
    %R=ctrl_fmr1_subs_mat.Subs_Mean_corrMat{1,k};
    R=ctrl_fmr1_subs_mat.Subs_Mean_corrMat_cleaned{1,k};
    %R=Data_corrMat2.control.Mean_corrMat{1,k};   
        
n=length(Nodes2.Mod_loc);

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

%[~,~,weights] = find(tril(R,-1));
weights = nonzeros(tril(R,-1));

%quantile(abs(weights),[0.025 0.25 0.50 0.75 0.975])


% create the graph object:
%G = graph(s,t,weights,n);
G = graph(R);

% mark the lines to remove from the graph:
threshold = 0; %  minimum correlation to plot
line_to_remove = isnan(weights) | abs(weights)<threshold;
% remove the lines from the graph:
G = G.rmedge(find(line_to_remove)); %#ok<FNDSB>

% plot it:
%subplot(1,7,count);
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'k';
% color positive in blue and negative in red:
  colormap jet;caxis([-1 1]);%colorbar;
  p.EdgeCData=G.Edges.Weight;
%p.EdgeColor = [G.Edges.Weight>0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight<0.'];


% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*3;
axis off

% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = Nodes2.Mod_loc(:,1);
y = Nodes2.Mod_loc(:,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_clust,'rbbggrm','.',20,'off');

view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold off;


saveas(gcf,strcat('subsCtrlFmr1_',num2str(k),'.svg'));
%count=count+1;

    end


