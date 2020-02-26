
%%%% this script is to generate a null model of using the WT nodes from the
%%%% fmr1 experiment. 


%%% i am using the aaft function from : https://au.mathworks.com/matlabcentral/fileexchange/16062-test-of-non-linearity
%%% which is used in this paper: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1369647
%%% and in turn bases its aaft in this one: D. Kugiumtzis, Surrogate data test for nonlinearity including monotonic transformations, Phys Rev E, vol. 62, 1, 2000


load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

datasets=['f20'; 'f60'; 's20'; 's60'];

%% doing the null graph
%%% I am first taking the timelines of the nodes to use them as sample to
%%% generate a surrogates with aaft
%%% as some nodes are empty I wil do an average of all of them 

%wtFish=fieldnames(NodesFmr1.wt.NodeMats); %%% to grab the WT from the fmr1
%dataset

f20Fish=fieldnames(Nodes.f20.NodeMats);

avgWTmat=[];
for i=1:length(f20Fish)
  
  tempMat=Nodes.f20.NodeMats.(f20Fish{i});  
  avgWTmat=cat(3,avgWTmat,tempMat);  
    
end

avgWTmat=nanmean(avgWTmat,3);

avgWTmat=avgWTmat(keep,:);

figure;imagesc(avgWTmat);
%%
%%%surrMat=aaft(avgWTmat,99); %%% it is not working... it doenst take the
%%%matrix

figure;plot(mean(avgWTmat));

surrMat=aaft(mean(avgWTmat),99);

surrMat=surrMat';

figure;imagesc(surrMat);

%%% this works... but not sure if it is representative... 
%%
%%% maybe I should do submatrices of each functional cluster and then stack them. 
good_surrMat=[];
for i=1:max(NodesFmr1.Nod_clustID)

    tempMean=mean(avgWTmat(find(NodesFmr1.Nod_clustID==i),:)); 
    figure;plot(tempMean)
   temp_surrMat=aaft(tempMean,length(find(NodesFmr1.Nod_clustID==i)));
   
   good_surrMat=vertcat(good_surrMat,temp_surrMat');
   
end


figure;imagesc(good_surrMat);
%%
%%% or maybe even do a surrogate timeline for each node...

very_good_surrMat=[];
for i=1:length(NodesFmr1.Nod_clustID)

    tempMean=avgWTmat(i,:); 
   %figure;plot(tempMean)
   temp_surrMat=aaft(tempMean,1);
   
   very_good_surrMat=vertcat(very_good_surrMat,temp_surrMat');
   
end


figure;imagesc(very_good_surrMat);

%%
%%%%% finally, I will try doing a surrogate for each loom

%%% I need to do it per loom. as previously
%%%from the wildtype loomhab dataset
load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx','S_trim');

super_good_surrMat=avgWTmat;

 for node=1:size(avgWTmat,1)              
        
     temp_mean=avgWTmat(node,:);
     
        for k=1:31
        
        if k==1
            temp_mean(:,10:45)=aaft(temp_mean(:,10:45),1);
        elseif k==11||k==21 ||k==31
            temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28)= aaft(temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28),1);
        else
            
            temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36)= aaft(temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36),1); 
                           
        end
        end
        
        super_good_surrMat(node,:)=temp_mean;
 end

figure;imagesc(super_good_surrMat);


%% making the graph

%%% just testing
%R_temp=corrcoef(good_surrMat');




Data_corrMat_Null=struct;

temp_mean=very_good_surrMat;
                
        temp_R={};
        
        for k=1:21
        
        if k==1
            temp_R{1,k}=corrcoef((temp_mean(:,10:45))');
        elseif k==11||k==21 
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28))');
        else
            
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36))'); 
                           
        end
        end
        
    Data_corrMat_Null.loomsR=temp_R;
        
%% doing graph analysis

%%%% plotting the matrices i used in the paper for f20


counter=1;
figure;
for data=1%:4
     datatemp=datasets(data,:);
subplot(1,8,counter);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,1}(keep,keep)); caxis([-1 1])%% for pre loom
subplot(1,8,counter+1);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep,keep));caxis([-1 1]) %% for 1st loom
subplot(1,8,counter+2);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,3}(keep,keep));caxis([-1 1])
subplot(1,8,counter+3);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,4}(keep,keep));caxis([-1 1])
subplot(1,8,counter+4);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,5}(keep,keep));caxis([-1 1])
subplot(1,8,counter+5);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,6}(keep,keep));caxis([-1 1])
subplot(1,8,counter+6);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11}(keep,keep));caxis([-1 1])%% for 10th loom
subplot(1,8,counter+7);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep,keep));caxis([-1 1]) %% for 11th loom
title(datatemp);

counter=counter+8;
end




 %%%% plotting the matrices of the null model
figure;

counter=1;
subplot(1,8,counter);imagesc(Data_corrMat_Null.loomsR{1,1}); caxis([-1 1])%% for pre loom
subplot(1,8,counter+1);imagesc(Data_corrMat_Null.loomsR{1,2});caxis([-1 1]) %% for 1st loom
subplot(1,8,counter+2);imagesc(Data_corrMat_Null.loomsR{1,3});caxis([-1 1])
subplot(1,8,counter+3);imagesc(Data_corrMat_Null.loomsR{1,4});caxis([-1 1])
subplot(1,8,counter+4);imagesc(Data_corrMat_Null.loomsR{1,5});caxis([-1 1])
subplot(1,8,counter+5);imagesc(Data_corrMat_Null.loomsR{1,6});caxis([-1 1])
subplot(1,8,counter+6);imagesc(Data_corrMat_Null.loomsR{1,11});caxis([-1 1])%% for 10th loom
subplot(1,8,counter+7);imagesc(Data_corrMat_Null.loomsR{1,12});caxis([-1 1]) %% for 11th loom
title('null model 1');

%% comparing matrices

%%% I will compare the matrices that I got with the null model with the
%%% wild type results from the fmr1 experiment. 

%%% using euclidean distance... this didnt seem to be that informative
%%% although is interesting to see how the WT from the fmr1 experiment and
%%% the s20 group have similar distances vs the null model. it also works
%%% with the f20 data instead of the WT from the fmr1 experiment

for i=1:21
D = pdist2(Data_corrMat_Null.loomsR{1,i},Data_corrMat4.control.Mean_corrMat{1,i},'euclidean');
Mat_dist_NullvsWT(i)=mean(mean(D));
end


%%% also s20 wt dataset
for i=1:21
D = pdist2(Data_corrMat_Null.loomsR{1,i},Data_corrMat2.s20.Mean_corrMat{1,i}(keep,keep),'euclidean'); %%% the keep is from the FnS dataset
Mat_dist_NullvsS20(i)=mean(mean(D));
end


figure;plot(Mat_dist_NullvsWT);hold on; plot(Mat_dist_NullvsS20);


%%%% testing doing multiple null models to play around with things. 

for j=1:10

    very_good_surrMat2=[];
for i=1:length(NodesFmr1.Nod_clustID)

    tempMean=avgWTmat(i,:); 
   %figure;plot(tempMean)
   temp_surrMat=aaft(tempMean,1);
   
   very_good_surrMat2=vertcat(very_good_surrMat2,temp_surrMat');
   
end
    
    
Data_corrMat_Null2=struct;

temp_mean=very_good_surrMat2;
                
        temp_R={};
        
        for k=1:21
        
        if k==1
            temp_R{1,k}=corrcoef((temp_mean(:,10:45))');
        elseif k==11||k==21 
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28))');
        else
            
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36))'); 
                           
        end
        end
        
    Data_corrMat_Null2.loomsR=temp_R;
    
    
    Mat_dist_NullvsNull2=[];
for i=1:21
D = pdist2(Data_corrMat_Null.loomsR{1,i},Data_corrMat_Null2.loomsR{1,i},'euclidean');
Mat_dist_NullvsNull2(i)=mean(mean(D));
end

hold on; plot(Mat_dist_NullvsNull2); %%%% this is a follow up from the previous section to plot on top of the other two responses

end



%% BCT analysis
%%%%% now further analyze the null model graph using
%%%%% the Brain Connectivity Toolbox


%%%% first just having a look at the corr coeficient values


    figure;
count=1;
   
    for k=[2 3 4 5 6 11 12]
    subplot(1,7,count)
       
    R=Data_corrMat_Null.loomsR{1,k};
        
    [~,~,weights] = find(tril(R,-1));
    %quantile(weights,[0.025 0.25 0.50 0.75 0.975])
    
    count=count+1;
    histogram(weights,20)
    end
    
%% with the BCT

%%
%%%% I will look at density, the number of degrees and strength.
%%%% I could maybe use them to see how they decrease through habituation or
%%%% increase at revovery

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected_null=struct;

    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat_Null.loomsR{1,moment(m)}),0.75);
     MatAll_corrected_null.(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end


%% density



loom=fieldnames(MatAll_corrected_null);

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected_null.(loom{i}).Mat);

MatAll_corrected_null.(loom{i}).kden=temp_kden;
end


figure;
temp=[];

loom=fieldnames(MatAll_corrected_null);

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected_null.(loom{i}).kden;
end
plot(temp);
hold on;


%%%% degrees and Strength



for i=1:length(loom)

deg=degrees_und(MatAll_corrected_null.(loom{i}).Mat);
str=strengths_und_sign(abs(MatAll_corrected_null.(loom{i}).Mat));

MatAll_corrected_null.(loom{i}).deg=deg;
MatAll_corrected_null.(loom{i}).str=str;
end


%%%% to plot histograms
figure;
for i=1:7
subplot(1,7,i);


histogram(MatAll_corrected_null.(loom{i}).deg,20);

hold on;
end


%%%% to plot with brains
figure;
counter=1;

    loom=fieldnames(MatAll_corrected_null);

for i=1:length(loom)
    if counter==1|counter==8|counter==15|counter==22
    low=0;high=100;
    else
    low=0;high=50; 
    end
    
  subplot(1,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected_null.(loom{i}).deg,'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end


%% participation and gateway coef

%P = participation_coef(X,modulesID);

%P=participation_coef(Matctrl_1,Nodes2.Mod_clust);
%[Gpos,~]=gateway_coef_sign(Matctrl_1,Nodes2.Mod_clust,1);


loom=fieldnames(MatAll_corrected_null);

for i=1:length(loom)

temp_mat=MatAll_corrected_null.(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes.Nod_clustID(keep));
[Gpos,~]=gateway_coef_sign(MatAll_corrected_null.(loom{i}).Mat,Nodes.Nod_clustID(keep),1);

MatAll_corrected_null.(loom{i}).P=P;
MatAll_corrected_null.(loom{i}).Gpos=Gpos;
end


%%% to plot it on brains
figure;
counter=1;

    loom=fieldnames(MatAll_corrected_null);

for i=1:length(loom)
    
    low=0;high=0.8;
     
  subplot(1,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected_null.(loom{i}).P,'filled');colormap('jet');view(-90,90);caxis([low high]);%colorbar; %
 
counter=counter+1;
end

%%
 %% now to plot some graphs. 
    

 figure;
 count=1;
    for k=[2 3 4 5 6 11 12]
   % figure;
    %set(gcf, 'Position',  [200, 200, 700, 900]);
    set(gcf, 'Position',  [200, 200, 1200, 900]);
    
    R=Data_corrMat_Null.loomsR{1,k};
    
    
    
n=length(Nodes.Nod_coor(keep));

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
threshold = 0.75; %  minimum correlation to plot
%threshold = 0; %  minimum correlation to plot
line_to_remove = isnan(weights) | abs(weights)<threshold;
% remove the lines from the graph:
G = G.rmedge(find(line_to_remove)); %#ok<FNDSB>


subplot(1,7,count);
% plot it:
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'k';
% color positive in blue and negative in red:
  colormap jet;caxis([-1 1]);%colorbar;
  %colormap jet;caxis([0 3]);%colorbar;
  p.EdgeCData=G.Edges.Weight;
%p.EdgeColor = [G.Edges.Weight>0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight<0.']; %% red high, blue low
%p.EdgeColor = [G.Edges.Weight<0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight>0.'];

% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*3;
%p.LineWidth = 1;
axis off

% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = Nodes.Nod_coor(keep,1);
y = Nodes.Nod_coor(keep,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
%gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_clust,'rgggbm','.',20,'off');
gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),Nodes.Nod_clustID(keep),'gggbrm','.',20,'off');

view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold off;

%saveas(gcf,strcat('subsf20vsf60_',num2str(k),'.svg'));

count=count+1;

    end

    %end



%% checking the FnS datasets to compare


load('graphs_FnS_all.mat')

%% checking the corr. distributions
%%%% first just having a look at the corr coeficient values

for data=1:4
    figure;
count=1;
    datatemp=datasets(data,:);
   sgtitle(datatemp) 
   
    for k=[2 3 4 5 6 11 12]
    subplot(1,7,count)
        
    R=Data_corrMat2.(datatemp).Mean_corrMat{1,k}(keep,keep);
        
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
MatAll_corrected=struct;
for data=1:4
    datatemp=datasets(data,:);
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat2.(datatemp).Mean_corrMat{1,moment(m)}(keep,keep)),0.75);
     MatAll_corrected.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end

%% density
%


%kden = density_und(CIJ);

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected.(datatemp).(loom{i}).Mat);

MatAll_corrected.(datatemp).(loom{i}).kden=temp_kden;
end
end

figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

%%%% degrees and Strength

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

deg=degrees_und(MatAll_corrected.(datatemp).(loom{i}).Mat);
str=strengths_und_sign(abs(MatAll_corrected.(datatemp).(loom{i}).Mat));

MatAll_corrected.(datatemp).(loom{i}).deg=deg;
MatAll_corrected.(datatemp).(loom{i}).str=str;
end
end

%%%% to plot histograms
figure;
for i=1:7
subplot(1,7,i);
for data=1:4
datatemp=datasets(data,:);
histogram(MatAll_corrected.(datatemp).(loom{i}).deg,20);

hold on;
end
end

%%%% to plot with brains
figure;
counter=1;

for data=1:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    if counter==1|counter==8|counter==15|counter==22
    low=0;high=100;
    else
    low=0;high=50; 
    end
    
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected.(datatemp).(loom{i}).deg,'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end
end

%% participation and gateway coef

%P = participation_coef(X,modulesID);

%P=participation_coef(Matctrl_1,Nodes2.Mod_clust);
%[Gpos,~]=gateway_coef_sign(Matctrl_1,Nodes2.Mod_clust,1);


for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

temp_mat=MatAll_corrected.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes.Nod_clustID(keep));
[Gpos,~]=gateway_coef_sign(MatAll_corrected.(datatemp).(loom{i}).Mat,Nodes.Nod_clustID(keep),1);

MatAll_corrected.(datatemp).(loom{i}).P=P;
MatAll_corrected.(datatemp).(loom{i}).Gpos=Gpos;
end
end

%%% to plot it on brains
figure;
counter=1;

for data=1:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    low=0;high=0.8;
     
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected.(datatemp).(loom{i}).P,'filled');colormap('jet');view(-90,90);caxis([low high]);%colorbar; %
 
counter=counter+1;
end
end


%% what if i do it again but with avgWTmat? like a positive control

Data_corrMat_avgWTmat=struct;

temp_mean=avgWTmat;
                
        temp_R={};
        
        for k=1:21
        
        if k==1
            temp_R{1,k}=corrcoef((temp_mean(:,10:45))');
        elseif k==11||k==21 
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28))');
        else
            
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36))'); 
                           
        end
        end
        
    Data_corrMat_avgWTmat.loomsR=temp_R;
        
%% doing graph analysis

figure;

counter=1;
subplot(1,8,counter);imagesc(Data_corrMat_avgWTmat.loomsR{1,1}); caxis([-1 1])%% for pre loom
subplot(1,8,counter+1);imagesc(Data_corrMat_avgWTmat.loomsR{1,2});caxis([-1 1]) %% for 1st loom
subplot(1,8,counter+2);imagesc(Data_corrMat_avgWTmat.loomsR{1,3});caxis([-1 1])
subplot(1,8,counter+3);imagesc(Data_corrMat_avgWTmat.loomsR{1,4});caxis([-1 1])
subplot(1,8,counter+4);imagesc(Data_corrMat_avgWTmat.loomsR{1,5});caxis([-1 1])
subplot(1,8,counter+5);imagesc(Data_corrMat_avgWTmat.loomsR{1,6});caxis([-1 1])
subplot(1,8,counter+6);imagesc(Data_corrMat_avgWTmat.loomsR{1,11});caxis([-1 1])%% for 10th loom
subplot(1,8,counter+7);imagesc(Data_corrMat_avgWTmat.loomsR{1,12});caxis([-1 1]) %% for 11th loom
title('avg model');

%%
%%%%% this script is to try to further analyze the graph of the original loomhab data using
%%%%% the Brain Connectivity Toolbox


%%%% first just having a look at the corr coeficient values


    figure;
count=1;
   
    for k=[2 3 4 5 6 11 12]
    subplot(1,7,count)
       
    R=Data_corrMat_avgWTmat.loomsR{1,k};
        
    [~,~,weights] = find(tril(R,-1));
    %quantile(weights,[0.025 0.25 0.50 0.75 0.975])
    
    count=count+1;
    histogram(weights,20)
    end
    
%% with the BCT

%%
%%%% I will look at density, the number of degrees and strength.
%%%% I could maybe use them to see how they decrease through habituation or
%%%% increase at revovery

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected_avg=struct;

    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat_avgWTmat.loomsR{1,moment(m)}),0.75);
     MatAll_corrected_avg.(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end


%% density



loom=fieldnames(MatAll_corrected_avg);

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected_avg.(loom{i}).Mat);

MatAll_corrected_avg.(loom{i}).kden=temp_kden;
end


figure;
temp=[];

loom=fieldnames(MatAll_corrected_avg);

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected_avg.(loom{i}).kden;
end
plot(temp);
hold on;


%%%% degrees and Strength



for i=1:length(loom)

deg=degrees_und(MatAll_corrected_avg.(loom{i}).Mat);
str=strengths_und_sign(abs(MatAll_corrected_avg.(loom{i}).Mat));

MatAll_corrected_avg.(loom{i}).deg=deg;
MatAll_corrected_avg.(loom{i}).str=str;
end


%%%% to plot histograms
figure;
for i=1:7
subplot(1,7,i);


histogram(MatAll_corrected_avg.(loom{i}).deg,20);

hold on;
end


%%%% to plot with brains
figure;
counter=1;

    loom=fieldnames(MatAll_corrected_avg);

for i=1:length(loom)
    if counter==1|counter==8|counter==15|counter==22
    low=0;high=100;
    else
    low=0;high=50; 
    end
    
  subplot(1,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected_avg.(loom{i}).deg,'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end


%% participation and gateway coef

%P = participation_coef(X,modulesID);

%P=participation_coef(Matctrl_1,Nodes2.Mod_clust);
%[Gpos,~]=gateway_coef_sign(Matctrl_1,Nodes2.Mod_clust,1);


loom=fieldnames(MatAll_corrected_avg);

for i=1:length(loom)

temp_mat=MatAll_corrected_avg.(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes.Nod_clustID(keep));
[Gpos,~]=gateway_coef_sign(MatAll_corrected_avg.(loom{i}).Mat,Nodes.Nod_clustID(keep),1);

MatAll_corrected_avg.(loom{i}).P=P;
MatAll_corrected_avg.(loom{i}).Gpos=Gpos;
end


%%% to plot it on brains
figure;
counter=1;

    loom=fieldnames(MatAll_corrected_avg);

for i=1:length(loom)
    
    low=0;high=0.8;
     
  subplot(1,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected_avg.(loom{i}).P,'filled');colormap('jet');view(-90,90);caxis([low high]);%colorbar; %
 
counter=counter+1;
end

%% having a look at the null model 2
%%%% surrogate timelines just at the loom timepoints


Data_corrMat_Null2=struct;

temp_mean=super_good_surrMat;
                
        temp_R={};
        
        for k=1:21
        
        if k==1
            temp_R{1,k}=corrcoef((temp_mean(:,10:45))');
        elseif k==11||k==21 
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28))');
        else
            
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36))'); 
                           
        end
        end
        
    Data_corrMat_Null2.loomsR=temp_R;
        
%% doing graph analysis

 %%%% plotting the matrices of the null model
figure;

counter=1;
subplot(1,8,counter);imagesc(Data_corrMat_Null2.loomsR{1,1}); caxis([-1 1])%% for pre loom
subplot(1,8,counter+1);imagesc(Data_corrMat_Null2.loomsR{1,2});caxis([-1 1]) %% for 1st loom
subplot(1,8,counter+2);imagesc(Data_corrMat_Null2.loomsR{1,3});caxis([-1 1])
subplot(1,8,counter+3);imagesc(Data_corrMat_Null2.loomsR{1,4});caxis([-1 1])
subplot(1,8,counter+4);imagesc(Data_corrMat_Null2.loomsR{1,5});caxis([-1 1])
subplot(1,8,counter+5);imagesc(Data_corrMat_Null2.loomsR{1,6});caxis([-1 1])
subplot(1,8,counter+6);imagesc(Data_corrMat_Null2.loomsR{1,11});caxis([-1 1])%% for 10th loom
subplot(1,8,counter+7);imagesc(Data_corrMat_Null2.loomsR{1,12});caxis([-1 1]) %% for 11th loom
title('null model 2');

%% with the BCT

%%
%%%% I will look at density, the number of degrees and strength.
%%%% I could maybe use them to see how they decrease through habituation or
%%%% increase at revovery

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected_null2=struct;

    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat_Null2.loomsR{1,moment(m)}),0.75);
     MatAll_corrected_null2.(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end


%% density



loom=fieldnames(MatAll_corrected_null2);

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected_null2.(loom{i}).Mat);

MatAll_corrected_null2.(loom{i}).kden=temp_kden;
end


figure;
temp=[];

loom=fieldnames(MatAll_corrected_null2);

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected_null2.(loom{i}).kden;
end
plot(temp);
%hold on;


%%%% degrees and Strength

for i=1:length(loom)

deg=degrees_und(MatAll_corrected_null2.(loom{i}).Mat);
str=strengths_und_sign(abs(MatAll_corrected_null2.(loom{i}).Mat));

MatAll_corrected_null2.(loom{i}).deg=deg;
MatAll_corrected_null2.(loom{i}).str=str;
end


%%%% to plot histograms
figure;
for i=1:7
subplot(1,7,i);


histogram(MatAll_corrected_null2.(loom{i}).deg,20);

hold on;
end


%% making some summary figures

%%%%%%%%%%%% the nodes timelines
figure;
subplot(1,3,1);imagesc(avgWTmat);title('avg model');
subplot(1,3,2);imagesc(very_good_surrMat);title('null model 1');
subplot(1,3,3);imagesc(super_good_surrMat);title('null model 2');


%%%%%%%%%%%%% matrices: 

%%% f20 used in the paper
counter=1;
figure;
for data=1%:4
     datatemp=datasets(data,:);
subplot(4,8,counter);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,1}(keep,keep)); caxis([-1 1])%% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep,keep));caxis([-1 1]) %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,3}(keep,keep));caxis([-1 1])
subplot(4,8,counter+3);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,4}(keep,keep));caxis([-1 1])
subplot(4,8,counter+4);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,5}(keep,keep));caxis([-1 1])
subplot(4,8,counter+5);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,6}(keep,keep));caxis([-1 1])
subplot(4,8,counter+6);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11}(keep,keep));caxis([-1 1])%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep,keep));caxis([-1 1]) %% for 11th loom
title(datatemp);

%counter=counter+8;
end

%%% average time lines
counter=counter+8;
subplot(4,8,counter);imagesc(Data_corrMat_avgWTmat.loomsR{1,1}); caxis([-1 1])%% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat_avgWTmat.loomsR{1,2});caxis([-1 1]) %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat_avgWTmat.loomsR{1,3});caxis([-1 1])
subplot(4,8,counter+3);imagesc(Data_corrMat_avgWTmat.loomsR{1,4});caxis([-1 1])
subplot(4,8,counter+4);imagesc(Data_corrMat_avgWTmat.loomsR{1,5});caxis([-1 1])
subplot(4,8,counter+5);imagesc(Data_corrMat_avgWTmat.loomsR{1,6});caxis([-1 1])
subplot(4,8,counter+6);imagesc(Data_corrMat_avgWTmat.loomsR{1,11});caxis([-1 1])%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat_avgWTmat.loomsR{1,12});caxis([-1 1]) %% for 11th loom
title('avg model');


%%% null model 1

counter=counter+8;
subplot(4,8,counter);imagesc(Data_corrMat_Null.loomsR{1,1}); caxis([-1 1])%% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat_Null.loomsR{1,2});caxis([-1 1]) %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat_Null.loomsR{1,3});caxis([-1 1])
subplot(4,8,counter+3);imagesc(Data_corrMat_Null.loomsR{1,4});caxis([-1 1])
subplot(4,8,counter+4);imagesc(Data_corrMat_Null.loomsR{1,5});caxis([-1 1])
subplot(4,8,counter+5);imagesc(Data_corrMat_Null.loomsR{1,6});caxis([-1 1])
subplot(4,8,counter+6);imagesc(Data_corrMat_Null.loomsR{1,11});caxis([-1 1])%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat_Null.loomsR{1,12});caxis([-1 1]) %% for 11th loom
title('null model 1');


%%% average time lines

counter=counter+8;
subplot(4,8,counter);imagesc(Data_corrMat_Null2.loomsR{1,1}); caxis([-1 1])%% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat_Null2.loomsR{1,2});caxis([-1 1]) %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat_Null2.loomsR{1,3});caxis([-1 1])
subplot(4,8,counter+3);imagesc(Data_corrMat_Null2.loomsR{1,4});caxis([-1 1])
subplot(4,8,counter+4);imagesc(Data_corrMat_Null2.loomsR{1,5});caxis([-1 1])
subplot(4,8,counter+5);imagesc(Data_corrMat_Null2.loomsR{1,6});caxis([-1 1])
subplot(4,8,counter+6);imagesc(Data_corrMat_Null2.loomsR{1,11});caxis([-1 1])%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat_Null2.loomsR{1,12});caxis([-1 1]) %% for 11th loom
title('null model 2');

%%%%%%%%% density

figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

%%% avg model
temp=[];
for i=1:length(loom)   
    temp(1,i)=MatAll_corrected_avg.(loom{i}).kden;
end
plot(temp);
hold on;

%%% null 1
temp=[];
for i=1:length(loom)   
    temp(1,i)=MatAll_corrected_null.(loom{i}).kden;
end
plot(temp);
hold on;

%%% null 3
temp=[];
for i=1:length(loom)   
    temp(1,i)=MatAll_corrected_null2.(loom{i}).kden;
end
plot(temp);
legend('f20','f60','s20','s60','avg model','null model 1','null model 2');
hold off;


%%%%%%%%%%%% degrees histograms per loom
edges=[0:5:100];
figure;

subplot(1,4,1);
for i=1:7
 data=1;
datatemp=datasets(data,:);
histogram(MatAll_corrected.(datatemp).(loom{i}).deg,edges);ylim([0 100]);title('f20');
hold on;
end

subplot(1,4,2);
for i=1:7
histogram(MatAll_corrected_avg.(loom{i}).deg,edges);ylim([0 100]);title('avg model');
hold on;
end

subplot(1,4,3);
for i=1:7
histogram(MatAll_corrected_null.(loom{i}).deg,edges);ylim([0 100]);title('null model 1');
hold on;
end

subplot(1,4,4);
for i=1:7
histogram(MatAll_corrected_null2.(loom{i}).deg,edges);ylim([0 100]);title('null model 2');
hold on;
end

%%%%%%%%%%% brain plots


%%%% to plot with brains
figure;

   low=0;high=75;
    

counter=1;
for i=1:length(loom)
       
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected.f20.(loom{i}).deg,'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
  title(strcat('f20','/',(loom{i})));  
counter=counter+1;
end


for i=1:length(loom)
    
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected_avg.(loom{i}).deg,'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 title(strcat('avg model','/',(loom{i})));
counter=counter+1;
end


for i=1:length(loom)
         
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected_null.(loom{i}).deg,'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 title(strcat('null model 1','/',(loom{i})));
counter=counter+1;
end


for i=1:length(loom)   
       
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected_null2.(loom{i}).deg,'filled');colormap('jet');caxis([low high]);view(-90,90);%colorbar; 
 title(strcat('null model 2','/',(loom{i})));
counter=counter+1;
end
