%%% this script is to make a graph of the connections that the habituation
%%% sensitive nodes based on participation change have. 

%%% for f20 10th minus 11th loom
test_mat=MatAll_corrected.f20.loom10.Mat-MatAll_corrected.f20.loom11.Mat;


big_subs_f20_10th_11th_habNodes=zeros(size(MatAll_corrected.f20.loom11.Mat));

big_subs_f20_10th_11th_habNodes(:,hab_nodes_Par_f20_f60)=test_mat(:,hab_nodes_Par_f20_f60);

big_subs_f20_10th_11th_habNodes(hab_nodes_Par_f20_f60,:)=test_mat(hab_nodes_Par_f20_f60,:);



%%% for the substracted f20-f60 matrices and then 10th minus 11th loom
test_mat=f20_f60_subs_mat.Subs_Mean_corrMat_cleaned{1,11}-f20_f60_subs_mat.Subs_Mean_corrMat_cleaned{1,12};


big_subs_f20_f60_habNodes=zeros(size(f20_f60_subs_mat.Subs_Mean_corrMat_cleaned{1,12}));

big_subs_f20_f60_habNodes(:,hab_nodes_Par_f20_f60)=test_mat(:,hab_nodes_Par_f20_f60);

big_subs_f20_f60_habNodes(hab_nodes_Par_f20_f60,:)=test_mat(hab_nodes_Par_f20_f60,:);



 figure;
    set(gcf, 'Position',  [200, 200, 700, 900]);
   
       
    %R=mini_subs_f20_f60_habNodes; 
    %R=f20_f60_subs_mat.Subs_Mean_corrMat_cleaned{1,11}-f20_f60_subs_mat.Subs_Mean_corrMat_cleaned{1,12};
    %R=big_subs_f20_f60_habNodes; 
    R=big_subs_f20_10th_11th_habNodes; 
    
n=length(Nodes2.Mod_loc(keep));

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

G = graph(R);

% plot it:
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'w';
p.MarkerSize = 2;
% color positive in blue and negative in red:
  colormap(RdBu) ;caxis([-1 1]);colorbar;
  %colormap jet;caxis([0 3]);%colorbar;
  p.EdgeCData=G.Edges.Weight;
%p.EdgeColor = [G.Edges.Weight>0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight<0.']; %% red high, blue low
%p.EdgeColor = [G.Edges.Weight<0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight>0.'];

% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*1;
%p.LineWidth = 1;
%axis off

p.NodeLabel=[];
p.EdgeAlpha=0.75;
% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = Nodes2.Mod_loc(keep,1);
y = Nodes2.Mod_loc(keep,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
%gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_clust,'rgggbm','.',20,'off');
gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60),1),Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60),2),Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60)),'rgggbm','.',25,'off');

view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'w');xlim([300 1350])
set(gca,'Color','k');%%% to make the background black
hold off;


%%

figure; set(gcf, 'Position',  [200, 200, 700, 900]);
gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_clust(keep),'wwwwww','.',15,'off');
hold on;
gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60),1),Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60),2),Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60)),'rgggbm','.',25,'off');
view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'w');xlim([300 1350])
set(gca,'Color','k');%%% to make the background black
hold off;

figure; set(gcf, 'Position',  [200, 200, 700, 900]);
gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_clust(keep),'kkkkkk','.',15,'off');
hold on;
gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60),1),Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60),2),Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60)),'rgggbm','.',25,'off');
view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([300 1350])
hold off;



%%

%%%% in this script I am trying to localize the nodes that had the biggest
%%%% change of participation to then look how their connectivity changed
%%%% from loom 10 to 11. 

%% if i want to check their corr substraction... is not the same than participation though
%[B,I] = sort(f20_f60_subs_mat.Subs_Mean_corrMat_cleaned{1,12}(hab_nodes_Par_f20_f60,11),'descend')


%% based on the top difference in Par at the recovery

[B,I] = sort(subs_Par_perLoom_f20_f60(hab_nodes_Par_f20_f60,11),'descend')

top_change=I(1:5);

figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_change)),1),Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_change)),2),Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60(top_change))),'rgggb','.',30,'on');xlim([300 1400]);view(-90,90);%colorbar; %

temp_clust=[];
counter=1;
for i=unique(Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60)))'
temp=find(Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60(I)))==i);

temp_clust(counter)=temp(1);

counter=counter+1;
end

top_change_clust=I(temp_clust);

%top_change_clust=I([14 2 10 1 11]);

figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_change_clust)),1),Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_change_clust)),2),Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60(top_change_clust))),'rgggb','.',30,'on');xlim([300 1400]);view(-90,90);%colorbar; %
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_change_clust)),1),Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_change_clust)),2),Nodes2.Mod_brain(keep(hab_nodes_Par_f20_f60(top_change_clust))),[],'.',30,'on');xlim([300 1400]);view(-90,90);%colorbar; %


figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_change_clust)),1),Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_change_clust)),2),Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60(top_change_clust))),[],'.',30,'on');xlim([300 1400]);view(-90,90);%colorbar; %


%%%% getting the degrees and the proportions to each cluster. 

%%%% first to put it on the right order
Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60(top_change_clust)))
%%% the nonhab needs to be moved to the end and fast hab first
top_change_clust=[1 5 13 17 18];

figure;
subplot(1,2,1);bar(MatAll_corrected.f20.loom10.NodeConnClustDeg(hab_nodes_Par_f20_f60(top_change_clust),:),'stacked');
subplot(1,2,2);bar(MatAll_corrected.f20.loom11.NodeConnClustDeg(hab_nodes_Par_f20_f60(top_change_clust),:),'stacked');

%%%% not bad... but I loose the inner responses among the blue cluster. i
%%%% will check with a node from them tectum. 


chosen_clust=[1 3 13 15 18]; %%% 1=tectum sharp-fasthab; 3=thalamus med-fasthab; 13=tectum broad-fasthab; 15=tectum slopehab and 18=pallium nonhab. 


figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(chosen_clust)),1),Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(chosen_clust)),2),Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60(chosen_clust))),'rgggb','.',30,'on');xlim([300 1400]);view(-90,90);%colorbar; %
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(chosen_clust)),1),Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(chosen_clust)),2),Nodes2.Mod_brain(keep(hab_nodes_Par_f20_f60(chosen_clust))),[],'.',30,'on');xlim([300 1400]);view(-90,90);%colorbar; %



figure;
subplot(1,2,1);H=bar(MatAll_corrected.f20.loom10.NodeConnClustDeg(hab_nodes_Par_f20_f60(chosen_clust),:),'stacked');
H(1).FaceColor = 'y';
H(2).FaceColor = 'g';
H(3).FaceColor = 'c';
H(4).FaceColor = 'b';
H(5).FaceColor = 'r';
H(6).FaceColor = 'm';
subplot(1,2,2);H=bar(MatAll_corrected.f20.loom11.NodeConnClustDeg(hab_nodes_Par_f20_f60(chosen_clust),:),'stacked');
H(1).FaceColor = 'y';
H(2).FaceColor = 'g';
H(3).FaceColor = 'c';
H(4).FaceColor = 'b';
H(5).FaceColor = 'r';
H(6).FaceColor = 'm';

chosen_clust_deg_loom10_11_f20=[];

for i=1:length(chosen_clust)
    
    temp1=MatAll_corrected.f20.loom10.NodeConnClustDeg(hab_nodes_Par_f20_f60(chosen_clust(i)),:);
    temp2=MatAll_corrected.f20.loom11.NodeConnClustDeg(hab_nodes_Par_f20_f60(chosen_clust(i)),:);
    
    temp=vertcat(temp1,temp2);
    
    chosen_clust_deg_loom10_11_f20=vertcat(chosen_clust_deg_loom10_11_f20,temp);
end


figure;
H=bar(chosen_clust_deg_loom10_11_f20,'stacked');
%%% with colors for fasthab
H(1).FaceColor = 'y';
H(2).FaceColor = 'g';
H(3).FaceColor = 'c';
H(4).FaceColor = 'b';
H(5).FaceColor = 'r';
H(6).FaceColor = 'm';


figure;
H=bar(chosen_clust_deg_loom10_11_f20,'stacked');
%%% all green for fasthab
H(1).FaceColor = 'g';
H(2).FaceColor = 'g';
H(3).FaceColor = 'g';
H(4).FaceColor = 'b';
H(5).FaceColor = 'r';
H(6).FaceColor = 'm';

%%  based on the top Par at the recovery for f20
%%% mmm not that interesting


[B,I] = sort(MatAll_corrected2.f20.loom11.P(hab_nodes_Par_f20_f60),'descend')

top_Par=I(1:5);

figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_Par)),1),Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_change)),2),Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60(top_change))),'rgggb','.',30,'on');xlim([300 1400]);view(-90,90);%colorbar; %

temp_clust=[];
counter=1;
for i=unique(Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60)))'
temp=find(Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60(I)))==i);

temp_clust(counter)=temp(1);

counter=counter+1;
end

top_Par_clust=I(temp_clust);


figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_Par_clust)),1),Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_Par_clust)),2),Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60(top_Par_clust))),'rgggb','.',30,'on');xlim([300 1400]);view(-90,90);%colorbar; %
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_Par_clust)),1),Nodes2.Mod_loc(keep(hab_nodes_Par_f20_f60(top_Par_clust)),2),Nodes2.Mod_brain(keep(hab_nodes_Par_f20_f60(top_Par_clust))),[],'.',30,'on');xlim([300 1400]);view(-90,90);%colorbar; %


%% side notes

%%% there are nonhab ROIs in the tegmentum??
%%% it seems so!! and I double checked in the original datasets and got the
%%% same distribution whe adding all the datasets together

hab_nodes_idx=[];
for i=1
  temp=Nodes2.Mod_KmeansID{keep(hab_nodes_Par_f20_f60(top_Par_clust(i)))};
    hab_nodes_idx=vertcat(hab_nodes_idx,temp);
end

figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(ROI_temp2_all(hab_nodes_idx,1),ROI_temp2_all(hab_nodes_idx,2));view(-90,90);%colorbar; %

%%%% what about slopehabs??
%%%% also, and they are also in the 4 datasets!!!


%% checking the 19 changing connections. 

%%% edges normalized to amount in loom 1. I get Inf values because
%%% sometimes i have divisions 1/0 when there was no connection in the
%%% first loom. I am ignoring them at the moment transforming them as NaN. 
%%% i am doing it per cluster 

ratio_connect_19=struct;

for data=1
    
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));
    
    temp_deg1=MatAll_corrected.f20.(loom{1}).NodeConnClustDeg(hab_nodes_Par_f20_f60,:);
        
   for k=1:length(loom)
       
     temp_ratio=(MatAll_corrected.f20.(loom{k}).NodeConnClustDeg(hab_nodes_Par_f20_f60,:))./temp_deg1; %%%% is generating NaNs when it divides by 0 
     
    for clust=[4 2 5 6 1 7]   
    temp_idx=find(Nodes2.Mod_clust(keep(hab_nodes_Par_f20_f60))==clust);
    
    ratio_connect_19.(loom{k}).(clustersF{clust}).ratios=temp_ratio(temp_idx,:);
    
    temp_ratio2=temp_ratio(temp_idx,:);
    
    temp_ratio2(find(temp_ratio2==Inf))=NaN;
    
    if size(temp_ratio2,1)==1
    
        ratio_connect_19.(loom{k}).(clustersF{clust}).means=temp_ratio2;
    else
    ratio_connect_19.(loom{k}).(clustersF{clust}).means=nanmean(temp_ratio2);  
    end
    
    end
    
   end
end


ratio_connect_19_all=[];
for data=1
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for clust=[4 2 5 6 1 7] 
    
 for  i=1:6 
     temp=[];
    for k=1:length(loom)
     temp(1,k)=ratio_connect_19.(loom{k}).(clustersF{clust}).means(i); 
    end 
    ratio_connect_19_all=vertcat(ratio_connect_19_all,temp);
 end 
 
end
end
 
%%%% is a bit of a mess... I will merge the fasthab and get rid of the
%%%% inhibited cluster cause it only give me NaNs
