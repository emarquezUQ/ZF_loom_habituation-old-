

%%%%% this script is to try to further analyze the graph of the original loomhab data using
%%%%% the Brain Connectivity Toolbox



load('Nodes_N_means_alldatasets.mat','Nodes','ROI_temp2_all','Zbrain_brainMask2D');

load('graph_loomNdataset.mat','Data_corrMat');

datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);


%%% testing the loop works
for data=1:4
     datatemp=datasets(data,:);
        
     fish=fieldnames(Nodes.(datatemp).corr_matrix);
     
    for tempfish=1:length(fish)
           
    Nodes.NaNtest.(datatemp).(fish{tempfish})
    end
end

%% checking the corr. distributions
%%%% first just having a look at the corr coeficient values

for data=1:4
    figure;
count=1;
    datatemp=datasets(data,:);
   sgtitle(datatemp) 
   
    for k=[2 3 4 5 6 11 12]
    subplot(1,7,count)
        
    R=Data_corrMat.(datatemp).Mean_corrMat{1,k};
        
    [~,~,weights] = find(tril(R,-1));
    %quantile(weights,[0.025 0.25 0.50 0.75 0.975])
    
    count=count+1;
    histogram(weights,20)
    end
    
end



%%
%%%% I will look at density, the number of degrees and strength.
%%%% I could maybe use them to see how they decrease through habituation or
%%%% increase at revovery

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll=struct;
for data=1:4
    datatemp=datasets(data,:);
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat.(datatemp).Mean_corrMat{1,moment(m)}),0.75);
     MatAll.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end

%% density
%


%kden = density_und(CIJ);

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)

temp_kden=density_und(MatAll.(datatemp).(loom{i}).Mat);

MatAll.(datatemp).(loom{i}).kden=temp_kden;
end
end

figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

%%%% degrees and Strength

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)

deg=degrees_und(MatAll.(datatemp).(loom{i}).Mat);
str=strengths_und_sign(abs(MatAll.(datatemp).(loom{i}).Mat));

MatAll.(datatemp).(loom{i}).deg=deg;
MatAll.(datatemp).(loom{i}).str=str;
end
end

%%%% to plot histograms
figure;
for i=1:7
subplot(1,7,i);
for data=1:4
datatemp=datasets(data,:);
histogram(MatAll.(datatemp).(loom{i}).str,20);

hold on;
end
end

%%%% to plot with brains
figure;
counter=1;

for data=1:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)
    if counter==1|counter==8|counter==15|counter==22
    low=0;high=100;
    else
    low=0;high=100; 
    end
    
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),25,MatAll.(datatemp).(loom{i}).deg,'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end
end

%%%% there are many combinations for substractions... depending on what we
%%%% would like to see. i will first test f20 minus f60 and then s20-s60


figure;
counter=1;

for data=3
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)
    
    low=-50;high=50;
    
subplot(1,7,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); 
hold on; scatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),25,(MatAll.(datatemp).(loom{i}).deg-MatAll.(datasets(data+1,:)).(loom{i}).deg),'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end
end

%% participation and gateway coef

%P = participation_coef(X,modulesID);

%P=participation_coef(Matctrl_1,Nodes.Mod_clust);
%[Gpos,~]=gateway_coef_sign(Matctrl_1,Nodes.Mod_clust,1);


for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)

temp_mat=MatAll.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes.Mod_clust);
[Gpos,~]=gateway_coef_sign(MatAll.(datatemp).(loom{i}).Mat,Nodes.Mod_clust,1);

MatAll.(datatemp).(loom{i}).P=P;
MatAll.(datatemp).(loom{i}).Gpos=Gpos;
end
end

%%% to plot it on brains
figure;
counter=1;

for data=1:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)
    
    low=0;high=0.8;
     
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),25,MatAll.(datatemp).(loom{i}).P,'filled');colormap('jet');view(-90,90);%caxis([low high]);%colorbar; %
 
counter=counter+1;
end
end

%% core vs periphery
% in the function guidelines it says that its supposed to have directed values...
% so i dont know if we can use this, i tried it anyway. 


%CoreP     = core_periphery_dir(W);

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)
temp_mat=MatAll.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
CoreP=core_periphery_dir(temp_mat);

MatAll.(datatemp).(loom{i}).CoreP=CoreP;
end
end

%%% to plot it on brains
figure;
counter=1;

for data=1:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)
    
    low=0;high=1;
     
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),25,MatAll.(datatemp).(loom{i}).CoreP,'filled');colormap('jet');view(-90,90);%caxis([low high]);%colorbar; %
 
counter=counter+1;
end
end


%% clustering coef
% i tried it... it might be interesting to see how the connections are lost
% during habituation but I dont think it gives us much more information. 
%%% I could calculate averages and compare between datasets

%Ccoef = clustering_coef_wu(W);
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)
temp_mat=MatAll.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
Ccoef=clustering_coef_wu(temp_mat);

MatAll.(datatemp).(loom{i}).Ccoef=Ccoef;
end
end

%%% to test particular combinations. 
data=1;
datatemp=datasets(data,:);
loom=fieldnames(MatAll.(datatemp));
i=1
h = ranksum(MatAll.f20.(loom{i}).Ccoef,MatAll.s20.(loom{i}).Ccoef) %%% 

%%% to plot it on brains
figure;
counter=1;

for data=1:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)
    
    low=0;high=1;
     
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),25,MatAll.(datatemp).(loom{i}).Ccoef,'filled');colormap('jet');view(-90,90);%caxis([low high]);%colorbar; %
 
counter=counter+1;
end
end


%% characteristic path length
%  i first need to get a connection-length matrix, then a distance matrix. I dont think the charpath will
%  really give me something interesting but it can be used to calculate teh
%  small-worldness. I could also compare it between datasets


for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)
temp_mat=MatAll.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
LMat = weight_conversion(temp_mat, 'lengths');
DMat = distance_wei(LMat);
[CharP,~,~,~,~]=charpath(DMat,0,0);

MatAll.(datatemp).(loom{i}).LMat=LMat;
MatAll.(datatemp).(loom{i}).DMat=DMat;
MatAll.(datatemp).(loom{i}).CharP=CharP;
end
end

%%
%%% calculating small-worldness as burgstaller et al 2019 bioRxiv did
%%% although I think I am doing this wrong... I think I should calculate
%%% all the values per fish. 

%%%% so far I can see some differences between genotypes and how the
%%%% differences drop with the stimulus but I am not really sure what to do
%%%% with it. 

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)

    SW=mean(MatAll.(datatemp).(loom{i}).Ccoef)/MatAll.(datatemp).(loom{i}).CharP;

MatAll.(datatemp).(loom{i}).SW=SW;

end
end

figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll.(datatemp).(loom{i}).SW;
end
plot(temp);
hold on;

end

%% modularity and modules based on correlation

% modularity and modules based on correlation
%%% it gives some intresting results but they are not the same as the
%%% HierarchicalConsensus method. it is much fasther though. 


%Ci = modularity_und(W);
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)
temp_mat=MatAll.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;    
    
mod=modularity_und(temp_mat);

MatAll.(datatemp).(loom{i}).mod=mod;
end
end

figure;
data=1;
datatemp=datasets(data,:);
loom=fieldnames(MatAll.(datatemp));
for i=1:length(loom)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); hold on; gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),MatAll.(datatemp).(loom{i}).mod);view(-90,90);  %colorbar;colormap('jet');caxis([-1 1]);
pause(3);
end

%%%%% HierarchicalConsensus. Based on Betzel 2018 bioRxiv


%%%% this part crashed overnight... i need to do it again. 

%%%% I will test with the one that seemed to have worked for the
%%%% fmr1ldataset
%%% I will get some nans to 0 cause I have in some matrices. 
%%% first with f20
temp_mat=MatAll.f20.loom1.Mat;
temp_mat(isnan(temp_mat))=0;

S2 = eventSamples(temp_mat, 150000); % 500000=5 modules. Not sure it was enought though... with 150000 i got the same results. 5 modules
[gamma_min,gamma_max]=gammaRange(temp_mat);
%S = fixedResSamples(temp_mat, 10000, 'Gamma', 1.2);
%S = exponentialSamples(temp_mat, 10000);
[Sc2, Tree2] = hierarchicalConsensus(S2,0.05);
[Sall2, thresholds2] = allPartitions(Sc2, Tree2);
%drawHierarchy(Sc2, Tree2)
C2 = coclassificationMatrix(S2);
figure;consensusPlot(C2, Sc2, Tree2)

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on;
gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Sall2(:,1),[],'.',25);view(-90,90);
 

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on;
gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Nodes.Mod_clust(:,1),[],'.',25);view(-90,90);


%%% now with f60
temp_mat=MatAll.f60.loom1.Mat;
temp_mat(isnan(temp_mat))=0;

S3 = eventSamples(temp_mat, 150000); % 150000=10 modules.
[gamma_min,gamma_max]=gammaRange(temp_mat);
%S = fixedResSamples(temp_mat, 10000, 'Gamma', 1.2);
%S = exponentialSamples(temp_mat, 10000);
[Sc3, Tree3] = hierarchicalConsensus(S3,0.05);
[Sall3, thresholds3] = allPartitions(Sc3, Tree3);
%drawHierarchy(Sc3, Tree3)
C3 = coclassificationMatrix(S3);
figure;consensusPlot(C3, Sc3, Tree3)

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on;
gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Sall3(:,1),[],'.',25);view(-90,90);

%%% now with s20
temp_mat=MatAll.s20.loom1.Mat;
temp_mat(isnan(temp_mat))=0;

S4 = eventSamples(temp_mat, 150000); % 150000=2 modules. with 500000=2 modules. same ones.
[gamma_min,gamma_max]=gammaRange(temp_mat);
%S = fixedResSamples(temp_mat, 10000, 'Gamma', 1.2);
%S = exponentialSamples(temp_mat, 10000);
[Sc4, Tree4] = hierarchicalConsensus(S4,0.05);
[Sall4, thresholds4] = allPartitions(Sc4, Tree4);
%drawHierarchy(Sc4, Tree4)
C4 = coclassificationMatrix(S4);
figure;consensusPlot(C4, Sc4, Tree4)

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on;
gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Sall4(:,1),[],'.',25);view(-90,90);


%%% now with s60
temp_mat=MatAll.s60.loom1.Mat;
temp_mat(isnan(temp_mat))=0;

S5 = eventSamples(temp_mat, 150000); % 150000=4 modules.
[gamma_min,gamma_max]=gammaRange(temp_mat);
%S = fixedResSamples(temp_mat, 10000, 'Gamma', 1.2);
%S = exponentialSamples(temp_mat, 10000);
[Sc5, Tree5] = hierarchicalConsensus(S5,0.05);
[Sall5, thresholds5] = allPartitions(Sc5, Tree5);
%drawHierarchy(Sc5, Tree5)
C5 = coclassificationMatrix(S5);
figure;consensusPlot(C5, Sc5, Tree5)

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on;
gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Sall5(:,1),[],'.',25);view(-90,90);



%%% for all of them. 

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll.(datatemp));

for i=1:length(loom)
  
    temp_mat=MatAll.(datatemp).(loom{i}).Mat;
    temp_mat(isnan(temp_mat))=0;

    S = eventSamples(temp_mat, 150000); % 150000 seems to work fine
    %[gamma_min,gamma_max]=gammaRange(temp_mat);
    %S = fixedResSamples(temp_mat, 10000, 'Gamma', 1.2);
    %S = exponentialSamples(temp_mat, 10000);
    [Sc, Tree] = hierarchicalConsensus(S,0.05);
    [Sall, thresholds] = allPartitions(Sc2, Tree);
    %drawHierarchy(Sc2, Tree2)
    C = coclassificationMatrix(S);
    % figure;consensusPlot(C2, Sc2, Tree2)

    % figure;
    % plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
    % hold on;
    % gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Sall2(:,1),[],'.',25);view(-90,90);
    
    
    
MatAll.(datatemp).(loom{i}).S=S;
MatAll.(datatemp).(loom{i}).Sc=Sc;
MatAll.(datatemp).(loom{i}).Tree=Tree;
MatAll.(datatemp).(loom{i}).Sall=Sall;
MatAll.(datatemp).(loom{i}).thresholds=thresholds;
MatAll.(datatemp).(loom{i}).thresholds=thresholds;
MatAll.(datatemp).(loom{i}).C=C;

end
end



