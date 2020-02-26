


%%%%% this script is to make the graphs that will go on the paper. 

%%% this is for the fmr1 dataset. 

%%% and using the nodes and cluster from the wild type dataset
load('NodesNgraphFmr1Loomhab2.mat');


load('Nodes_N_means_alldatasets2.mat','Zbrain_brainMask2D');


cbrewer()

[RdYlBu]=cbrewer('div','RdYlBu',101);
[RdBu]=cbrewer('div','RdBu',101);
[RdBu2]=cbrewer('div','RdBu',1001);
[PRGn]=cbrewer('div','PRGn',101);
[PiYG]=cbrewer('div','PiYG',101);

%%


groupnames=fieldnames(Nodes2.ROIs_idx);


%% checking a few things


%%% to see which ROIs have less representation
counter=1;
figure;
meanProp=[];
NaN_nodes={};
for g=1:3
    
  group=groupnames{g,1};
    
    fish=fieldnames(Nodes2.NaNtest.(group));  
     
    
    Matrix_mean=[];
    temp_NaN_nodes=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes2.NaNtest.(group).(fish{f}));  
 temp=double(diag(Nodes2.NaNtest.(group).(fish{f})));
 temp_NaN_nodes=horzcat(temp_NaN_nodes,temp);
end

%Matrix_mean=(sum(Matrix_mean,3))/length(fish); %% is the same
Matrix_mean=nanmean(Matrix_mean,3);  


subplot(1,3,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),20,(1-mean(Matrix_mean)),'filled'); colorbar;colormap('jet'); caxis([0 1]);
view(-90,90);
title(group);
hold off;

counter=counter+1;

meanProp=horzcat(meanProp,(mean(Matrix_mean)'));

NaN_nodes{g}=temp_NaN_nodes;
end

fish_perNode=[];
meanProp_good=[];
for  g=1:3
    fish_perNode(:,g)=abs(sum(NaN_nodes{g},2)-size(NaN_nodes{g},2));
    meanProp_good(:,g)=1-(sum(NaN_nodes{g},2)/size(NaN_nodes{g},2));
    
end

%% what if i put a threshold on proportion of fish needed to contribute to a node?
%%%% before I run the crosscorrelation. 
%%% with at least 25% of the fish i loose some nodes in the hindbrain although a visual effect still
%%% seems to remain there. 

%discard=find(min((1-meanProp)')<0.25); %%% I could use 0.25 (discard 8 nodes, but leaves some in the hindbrain), 0.33 (discard 18 nodes) , or 0/5 (discard 60 nodes).
discard=find(min(meanProp_good,[],2)<0.25); %%% this is the proper way!!! with 0.25 i discard 4 nodes, 0.33=7 nodes and 0.5=27 nodes
keep=find(ismember([1:89],discard)==0);


%%%% with 0.25, i discarded mostly from contralateral side. 
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_brain(keep));view(-90,90);
title('Model Nodes');
%%% adding numbers
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_brain(:));view(-90,90);
title('Model Nodes');



%%

%%% first, the matrices of the 1st, 10th, 11th and loom similar to 11th. I
%%% need to find a way to get it in a cuantitative way.


%%% probably a correlation is the easiest? 
similar_to_rec=struct;
for g=1:3
     group=groupnames{g,1};
     
     temp_mat1=Data_corrMat2.(group).Mean_corrMat{1,12}(keep,keep);
    temp_mat1(isnan(temp_mat1))=0;

    temp=[];
    temp2=[];
    
for k=1:11
    temp_mat2=Data_corrMat2.(group).Mean_corrMat{1,k}(keep,keep);
    temp_mat2(isnan(temp_mat2))=0;
    temp_corr=corrcoef(temp_mat1,temp_mat2);
    temp(k)=temp_corr(1,2);
    temp_mean=mean(mean(temp_mat2));
    temp2(k)=temp_mean;
end
similar_to_rec.(group).corr=temp;
similar_to_rec.(group).Mean_corr=temp2;
end

%%%% to see how the correlations drop during habituation. is more of the
%%%% same but could be useful... 
figure;
for g=1:3
     group=groupnames{g,1};
plot(similar_to_rec.(group).Mean_corr);
hold on;
end

%%%% so, it seems that for all of them the most similar one to recovery is
%%%% on the second loom... so it might not be very relevant to show it in the
%%%% figure. or maybe to show that that is not a difference?
for g=1:3
     group=groupnames{g,1};
    
    [M,I] = max(similar_to_rec.(group).corr);
    %similar_to_rec.(group).corr
    I-1 %%% i take one out cause the first one is the pre loom moment. so this result is the actual loom
end


%% for contros vs fmr1

%% first the matrices
%%%%
counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 g=3; %% controls
     group=groupnames{g,1};
subplot(2,4,counter);imagesc(Data_corrMat2.(group).Mean_corrMat{1,2}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for first loom
subplot(2,4,counter+1);imagesc(Data_corrMat2.(group).Mean_corrMat{1,11}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for 10th loom
subplot(2,4,counter+2);imagesc(Data_corrMat2.(group).Mean_corrMat{1,12}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for 11th loom
subplot(2,4,counter+3);imagesc(Data_corrMat2.(group).Mean_corrMat{1,3}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for the 2nd loom
title(group);
g=2; %% fmr1
     group=groupnames{g,1};
subplot(2,4,counter+4);imagesc(Data_corrMat2.(group).Mean_corrMat{1,2}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for first loom
subplot(2,4,counter+5);imagesc(Data_corrMat2.(group).Mean_corrMat{1,11}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for 10th loom
subplot(2,4,counter+6);imagesc(Data_corrMat2.(group).Mean_corrMat{1,12}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for 11th loom
subplot(2,4,counter+7);imagesc(Data_corrMat2.(group).Mean_corrMat{1,3}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for the 2nd loom
title(group);
%counter=counter+8;




%%% and now with the flitered matrix above 0.75

MatAll_corrected2=struct;
for g=1:3
    group=groupnames{g,1};
    moment=[2 3 4 5 6 7 8 9 10 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 6 7 8 9 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat2.(group).Mean_corrMat{1,moment(m)}(keep,keep)),0.75);
     MatAll_corrected2.(group).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end

%%%% checking if the correlation falls in the same place...
%%% probably a correlation is the easiest? 
similar_to_rec2=struct;
for g=1:3
     group=groupnames{g,1};
     loom=fieldnames(MatAll_corrected2.(group));
     temp_mat1=MatAll_corrected2.(group).(loom{11}).Mat;
    temp_mat1(isnan(temp_mat1))=0;

    temp=[];
    temp2=[];
for k=1:10
    temp_mat2=MatAll_corrected2.(group).(loom{k}).Mat;
    temp_mat2(isnan(temp_mat2))=0;
    temp_corr=corrcoef(temp_mat1,temp_mat2);
    temp(k)=temp_corr(1,2);
    temp_mean=mean(mean(temp_mat2));
    temp2(k)=temp_mean;
end
similar_to_rec2.(group).corr=temp;
similar_to_rec2.(group).Mean_corr=temp2;
end

%%%% to see how the correlations drop during habituation. is more of the
%%%% same but could be useful... 
figure;
for g=1:3
     group=groupnames{g,1};
plot(similar_to_rec2.(group).Mean_corr);
hold on;
end

%%%% all in the second... so more of the same
for g=1:3
     group=groupnames{g,1};
    
    [M,I] = max(similar_to_rec2.(group).corr);
    %similar_to_rec.(group).corr
    I %%% 
end

%%%%%
counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
  g=3; %% controls
     group=groupnames{g,1};
     loom=fieldnames(MatAll_corrected2.(group));
subplot(2,4,counter);imagesc(MatAll_corrected2.(group).(loom{1}).Mat); pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar;%% for pre loom
subplot(2,4,counter+1);imagesc(MatAll_corrected2.(group).(loom{10}).Mat);pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar;%% for 10th loom
subplot(2,4,counter+2);imagesc(MatAll_corrected2.(group).(loom{11}).Mat);pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar; %% for 11th loom
subplot(2,4,counter+3);imagesc(MatAll_corrected2.(group).(loom{4}).Mat);pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar; %% for the 4th loom
title(group);
 g=2; %% fmr1
     group=groupnames{g,1};
     loom=fieldnames(MatAll_corrected2.(group));
subplot(2,4,counter+4);imagesc(MatAll_corrected2.(group).(loom{1}).Mat); pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar;%% for pre loom
subplot(2,4,counter+5);imagesc(MatAll_corrected2.(group).(loom{10}).Mat);pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar;%% for 10th loom
subplot(2,4,counter+6);imagesc(MatAll_corrected2.(group).(loom{11}).Mat);pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar; %% for 11th loom
subplot(2,4,counter+7);imagesc(MatAll_corrected2.(group).(loom{6}).Mat);pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar; %% for the 6th loom
title(group);


%% making some matrices substraction of substractions

counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 g=3; %% controls
     group=groupnames{g,1};
     loom=fieldnames(MatAll_corrected2.(group));
subplot(2,3,counter);imagesc(MatAll_corrected2.(group).(loom{1}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{1}).Mat); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for pre loom
subplot(2,3,counter+1);imagesc(MatAll_corrected2.(group).(loom{10}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{10}).Mat);pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for 10th loom
subplot(2,3,counter+2);imagesc(MatAll_corrected2.(group).(loom{11}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{11}).Mat);pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for 11th loom
title(group);


%%%% to try to put the empty places with black boxes. 

sub1=MatAll_corrected2.(group).(loom{1}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{1}).Mat;
sub10=MatAll_corrected2.(group).(loom{10}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{10}).Mat;
sub11=MatAll_corrected2.(group).(loom{11}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{11}).Mat;

max(max(sub1))
min(min(sub11))

sub1(find(sub1==0)) = -1;
sub10(find(sub10==0)) = -1;
sub11(find(sub11==0)) = -1;

ddd=[0 0 0;RdBu];

counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 g=3; %%controls
     group=groupnames{g,1};
     loom=fieldnames(MatAll_corrected2.(group));
subplot(2,3,counter);imagesc(sub1); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);colorbar;%% for pre loom
subplot(2,3,counter+1);imagesc(sub10);pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);colorbar;%% for 10th loom
subplot(2,3,counter+2);imagesc(sub11);pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);colorbar; %% for 11th loom
title('control minus fmr1');


%% substraction of connections
%%%% I am trying control minus fmr1 first. 


%%
%%%% here i am trying to get the connections from my previous filtering
  
for k=[1:21]
  
    temp_mat1=Data_corrMat2.control.Mean_corrMat{1,k}(keep,keep);
    temp_mat1(isnan(temp_mat1))=0;
    
    temp_mat1_idx=find(threshold_absolute(abs(temp_mat1),0.75));
    
    temp_mat2=Data_corrMat2.fmr1.Mean_corrMat{1,k}(keep,keep);
    temp_mat2(isnan(temp_mat2))=0;
    
    temp_mat2_idx=find(threshold_absolute(abs(temp_mat2),0.75));
    
    temp_mat=temp_mat1-temp_mat2;
    
    temp_mat_idx=union(temp_mat1_idx,temp_mat2_idx);
    
    temp_mat_cleaned=zeros(size(temp_mat));
    temp_mat_cleaned(temp_mat_idx)=temp_mat(temp_mat_idx);
    
    ctrl_fmr1_subs_mat.Subs_Mean_corrMat_cleaned{1,k}=temp_mat_cleaned;
    
end

    
    %% now to plot some graphs. 
    
%     figure;
%      count=1;
%     for g=1:3
%  group=groupnames{g,1};
 %figure;set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%sgtitle(group);
 count=1;
    for k=[2 3 4 5 6 11 12]
    figure;
    set(gcf, 'Position',  [200, 200, 700, 900]);
    %set(gcf, 'Position',  [200, 200, 1200, 900]);
    
    
    R=ctrl_fmr1_subs_mat.Subs_Mean_corrMat_cleaned{1,k}; 
    
    
n=length(Nodes2.Mod_loc(keep));

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
%threshold = 0; %  minimum correlation to plot
%threshold = 0; %  minimum correlation to plot
%line_to_remove = isnan(weights) | abs(weights)<threshold;
% remove the lines from the graph:
%G = G.rmedge(find(line_to_remove)); %#ok<FNDSB>


%subplot(1,3,count);
% plot it:
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'k';
% color positive in blue and negative in red:
  colormap(RdBu) ;caxis([-1 1]);colorbar;
  %colormap jet;caxis([0 3]);%colorbar;
  p.EdgeCData=G.Edges.Weight;
%p.EdgeColor = [G.Edges.Weight>0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight<0.']; %% red high, blue low
%p.EdgeColor = [G.Edges.Weight<0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight>0.'];

% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*3;
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
gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_clust(keep),'gggbrm','.',15,'off');


view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'w');xlim([300 1350])
set(gca,'Color','k');%%% to make the background black
hold off;

%saveas(gcf,strcat('subsf20vsf60_',num2str(k),'.svg'));

count=count+1;

    end

    %end

%% now degrees, strenth and participation

%%%% find ways to show it in a nice way. 


for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll_corrected2.(group));

for i=1:length(loom)

deg=degrees_und(MatAll_corrected2.(group).(loom{i}).Mat);
str=strengths_und_sign(abs(MatAll_corrected2.(group).(loom{i}).Mat));

MatAll_corrected2.(group).(loom{i}).deg=deg;
MatAll_corrected2.(group).(loom{i}).str=str;
end
end


%%%% to plot with brains
figure;
counter=1;

for g=[3 2]
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));

for i=[1:5 10 11]
    if counter==1|counter==8|counter==15|counter==22
    low=0;high=100;
    else
    low=0;high=50; 
    end
    
  subplot(2,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,MatAll_corrected2.(group).(loom{i}).deg,'filled');colormap(inferno);caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end
end



%%%% there are many combinations for substractions... depending on what we
%%%% would like to see. i will first test f20 minus f60 and then s20-s60
%%% not sure if including them.

figure;
counter=1;

for g=3
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));

for i=[1:5 10 11]
    
    low=-50;high=50;
    
subplot(1,7,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); 
hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,(MatAll_corrected2.(group).(loom{i}).deg-MatAll_corrected2.(groupnames{g-1,1}).(loom{i}).deg),'filled');colormap(RdBu);caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end
end

%%%% degrees as rasterplot

figure;
counter=1;
for g=[3 2]
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));
temp=[];
for i=1:length(loom)

    temp(:,i)=(MatAll_corrected2.(group).(loom{i}).str)';
    

end 

subplot(1,2,counter);
imagesc(temp);caxis([0 25]);colormap(inferno); colorbar;
counter=counter+1;
end

%%% raster substraction
figure;
counter=1;
for g=3
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));
temp=[];
for i=1:length(loom)

    temp(:,i)=(MatAll_corrected2.(group).(loom{i}).deg-MatAll_corrected2.(groupnames{g-1,1}).(loom{i}).deg)';
    

end 

imagesc(temp);caxis([-25 25]);colormap(RdBu); colorbar;
counter=counter+1;
end


%%%% the nodes that are the most senstive in fmr1.
hab_nodes_deg_ctrl_fmr1=find((temp(:,2)<-5 | temp(:,3)<-5) & temp(:,11)<-5);
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(Nodes2.Mod_loc(keep(hab_nodes_deg_ctrl_fmr1),1),Nodes2.Mod_loc(keep(hab_nodes_deg_ctrl_fmr1),2),Nodes2.Mod_clust(keep(hab_nodes_deg_ctrl_fmr1)));view(-90,90);%colorbar; %
 figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(Nodes2.Mod_loc(keep(hab_nodes_deg_ctrl_fmr1),1),Nodes2.Mod_loc(keep(hab_nodes_deg_ctrl_fmr1),2),Nodes2.Mod_brain(keep(hab_nodes_deg_ctrl_fmr1)));view(-90,90);%colorbar; %
 


%%% degrees or strength as timelines for all datasets.
%%% they drop very quickly and kind of hit the floor... not that interesting to show. 
figure;

for g=1:3
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));
temp=[];
for i=1:length(loom)
    
    temp(:,i)=nanmean(MatAll_corrected2.(group).(loom{i}).deg);
    
end 

plot(temp);
hold on;
end


%%%% participation

for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll_corrected2.(group));

for i=1:length(loom)

temp_mat=MatAll_corrected2.(group).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes2.Mod_clust(keep));
[Gpos,~]=gateway_coef_sign(MatAll_corrected2.(group).(loom{i}).Mat,Nodes2.Mod_clust(keep),1);

MatAll_corrected2.(group).(loom{i}).P=P;
MatAll_corrected2.(group).(loom{i}).Gpos=Gpos;
end
end

%%%% interesting... the fmr1 have more participation in loom 2 and recovery
%%% in the brains
figure;
counter=1;

for g=[3 2]
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));

for i=[1:5 10 11]
    
    low=0;high=0.8;
     
  subplot(2,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,MatAll_corrected2.(group).(loom{i}).P,'filled');colormap(inferno);view(-90,90);caxis([low high]);%colorbar; %
 title(group);
counter=counter+1;
end
end

%%% raster
figure;
counter=1;
for g=[3 2]
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));
temp=[];
for i=1:length(loom)

    temp(:,i)=(MatAll_corrected2.(group).(loom{i}).P)';
    

end 

subplot(1,2,counter);
imagesc(temp);caxis([0 0.8]);colormap(inferno); colorbar;
counter=counter+1;
end

%%% raster substraction
figure;
counter=1;
subs_Par_perLoom_f20_f60=[];
for g=3
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));

for i=1:length(loom)

    subs_Par_perLoom_f20_f60(:,i)=(MatAll_corrected2.(group).(loom{i}).P-MatAll_corrected2.(groupnames{g-1,1}).(loom{i}).P)';
    

end 

imagesc(subs_Par_perLoom_f20_f60);caxis([-0.8 0.8]);colormap(RdBu); colorbar;
counter=counter+1;
end


%%%%  the nodes that are the most senstive to the for fmr1. 
hab_nodes_Par_ctrl_fmr1=find((subs_Par_perLoom_f20_f60(:,3)<-0.1 | subs_Par_perLoom_f20_f60(:,4)<-0.1) & subs_Par_perLoom_f20_f60(:,11)<-0.1);
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_ctrl_fmr1),1),Nodes2.Mod_loc(keep(hab_nodes_Par_ctrl_fmr1),2),Nodes2.Mod_clust(keep(hab_nodes_Par_ctrl_fmr1)));view(-90,90);%colorbar; %
 figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
hold on; gscatter(Nodes2.Mod_loc(keep(hab_nodes_Par_ctrl_fmr1),1),Nodes2.Mod_loc(keep(hab_nodes_Par_ctrl_fmr1),2),Nodes2.Mod_brain(keep(hab_nodes_Par_ctrl_fmr1)));view(-90,90);%colorbar; %



%%% participation as timelines for all datasets 
%%% interesting fmr1 stays at the top... and hets are the same as controls.
figure;
for g=1:3
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));
temp=[];
for i=1:length(loom)
    
    temp(:,i)=mean(MatAll_corrected2.(group).(loom{i}).P(hab_nodes_Par_ctrl_fmr1));
    
end 

plot(temp);
hold on;
end


%%%%% other graphs I could do...

%%%% clustering coef

for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll_corrected2.(group));

for i=1:length(loom)
temp_mat=MatAll_corrected2.(group).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
Ccoef=clustering_coef_wu(temp_mat);

MatAll_corrected2.(group).(loom{i}).Ccoef=Ccoef;
end
end


%%% maybe there is a bit more for fmr1... but is not big differences
%%% to plot it on brains
figure;
counter=1;

for g=[3 2]
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));

for i=[1:5 10 11]
    
    low=0;high=1;
     
  subplot(2,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,MatAll_corrected2.(group).(loom{i}).Ccoef,'filled');colormap(inferno);view(-90,90);caxis([low high]);%colorbar; %
 
counter=counter+1;
end
end


%%% raster
figure;
counter=1;
for g=[3 2]
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));
temp=[];
for i=1:length(loom)

    temp(:,i)=(MatAll_corrected2.(group).(loom{i}).Ccoef)';
    

end 

subplot(1,2,counter);
imagesc(temp);caxis([0 1]);colormap(inferno); colorbar;
counter=counter+1;
end

%%% timelines for al datasets
%%% more of the same... 
figure;
for g=1:3
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));
temp=[];
for i=1:length(loom)
    
    temp(:,i)=mean(MatAll_corrected2.(group).(loom{i}).Ccoef);
    
end 

plot(temp);
hold on;
end


%%% to plot it on brains... making an average. 
figure;
counter=1;
Temp=[];
for g=1:3
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));

    temp=[];
for i=1:length(loom)
temp(:,i)=MatAll_corrected2.(group).(loom{i}).Ccoef;


end  
temp=mean(temp');
Temp(data,:)=temp;

    low=0;high=1;
     
  subplot(1,4,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,temp,'filled');colormap(inferno);view(-90,90);%caxis([low high]);%colorbar; %
 
counter=counter+1;

end

Temp=mean(Temp);
low=0;high=1;
     
  figure;
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,Temp,'filled');colormap(inferno);view(-90,90);%caxis([low high]);%colorbar; %
 
 main_visual_idx2=find(Temp>0.5);


