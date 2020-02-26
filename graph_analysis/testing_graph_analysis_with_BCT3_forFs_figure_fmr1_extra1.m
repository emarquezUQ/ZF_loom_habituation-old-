
 
%%  trying to make circular graph

%%% trying to make circular graph
g=3;
    group=groupnames{g,1};
%%% substractions of ctrl-fmr1
sub1=MatAll_corrected2.(group).(loom{1}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{1}).Mat;
sub2=MatAll_corrected2.(group).(loom{2}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{2}).Mat;
sub3=MatAll_corrected2.(group).(loom{3}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{3}).Mat;
sub4=MatAll_corrected2.(group).(loom{4}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{4}).Mat;
sub5=MatAll_corrected2.(group).(loom{5}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{5}).Mat;
sub10=MatAll_corrected2.(group).(loom{10}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{10}).Mat;
sub11=MatAll_corrected2.(group).(loom{11}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{11}).Mat;


%%% first with okomarov-schemaball

figure;schemaball(sub11);

schemaball(sub11,num2str(Nodes4.Mod_brain(keep)),[1,0,0;0 0 1],[]);
 
    h = schemaball;
set(h.l(~isnan(h.l)), 'LineWidth',1.2)
set(h.s, 'MarkerEdgeColor','red','LineWidth',2,'SizeData',100)
set(h.t, 'EdgeColor','white','LineWidth',1)


%%%%%



%%% and with paul-kassebaum-mathworks-circularGraph

%%% not great... there are some good things but I can only use binary data.
%%% although maybe if you play with it more. 

%%% it cannot take nans or negative values... 
sub11(find(isnan(sub11)))= 0;
% thresh = 0.75;
% sub11(sub11 >  thresh) = 1;
% sub11(sub11 <= thresh) = 0;

figure;circularGraph(abs(sub11));

%%%% this is with some help with gilles:

[Blues]=cbrewer('seq','Blues',length(sub11));
[Reds]=cbrewer('seq','Reds',length(sub11));

%%%% so the way we found to put both negative and positive values is to do
%%%% the absolutes values and overlay the two graphs. or to put them
%%%% separated. 
sub_temp=sub11;sub_temp(sub11<0)=0;
figure;
circularGraph(sub_temp/5,'Colormap',Blues);
hold on;
sub_temp=sub11;sub_temp(sub11>0)=0;
%figure;
circularGraph(abs(sub_temp/5),'Colormap',Reds);

loom=fieldnames(MatAll_corrected2.(group));
ctrl_11=MatAll_corrected2.(group).(loom{11}).Mat;
ctrl_11(find(isnan(ctrl_11)))= 0;

fmr1_11=MatAll_corrected2.(groupnames{g-1,1}).(loom{11}).Mat;
fmr1_11(find(isnan(fmr1_11)))= 0;

figure;
circularGraph(ctrl_11,'Colormap',Blues);
figure;
circularGraph(fmr1_11,'Colormap',Reds);

%%%% after showing it to Ethan. he likes the order per functional cluster,
%%%% a graph for each genotype and the loom graphs not the substractions. 


%%% creating a functional cluster colormap
funct_clust_colors=zeros(length(Nodes4.Mod_clust(keep)),3);
for i=unique(Nodes4.Mod_clust(keep))'
    
    temp_node=find(Nodes4.Mod_clust(keep)==i);
    
    if i<4    
    funct_clust_colors(temp_node,:)= repmat([0 1 0],size(temp_node,1),1);  
    elseif i==4      
     funct_clust_colors(temp_node,:)= repmat([0 0 1],size(temp_node,1),1);          
    elseif i==5
     funct_clust_colors(temp_node,:)= repmat([1 0 0],size(temp_node,1),1);      
    elseif i==6
     funct_clust_colors(temp_node,:)= repmat([1 0 1],size(temp_node,1),1);                     
    end   
    
end


black_nodes=zeros(length(Nodes4.Mod_clust(keep)),3);

figure;
circularGraph(ctrl_11,'Colormap',funct_clust_colors);
figure;
circularGraph(fmr1_11,'Colormap',funct_clust_colors);



figure;
circularGraph(ctrl_11,'ColorMapNode',funct_clust_colors,'ColorMapEdges',Blues);
figure;
circularGraph(fmr1_11,'ColorMapNode',funct_clust_colors,'ColorMapEdges',Reds);


%%% trying to make a colormap
[Blues]=cbrewer('seq','Blues',length(nonzeros(ctrl_11)));
[Reds]=cbrewer('seq','Reds',length(nonzeros(fmr1_11)));

BLUE=repmat([0 0 1],length(nonzeros(ctrl_11)),1);
RED=repmat([1 0 0],length(nonzeros(fmr1_11)),1);

BLUE2=zeros(size(BLUE));
BLUE2(:,1)=nonzeros(ctrl_11);
RED2=zeros(size(RED));
RED2(:,1)=nonzeros(fmr1_11);

figure;h=circularGraph(fmr1_11/5,'ColorMapNode',funct_clust_colors,'ColorMapEdges',Reds);


[~,~,bin] = histcounts(nonzeros(ctrl_11),20);
[Blues]=cbrewer('seq','Blues',20);
ctrl_11_blues=Blues(bin,:);

[BluesDark]=cbrewer('seq','Blues',40);
BluesDark=BluesDark(21:40,:);
ctrl_11_blues_dark=BluesDark(bin,:);

[~,~,bin] = histcounts(nonzeros(fmr1_11),20);
[Reds]=cbrewer('seq','Reds',20);
fmr1_11_reds=Reds(bin,:);

[RedsDark]=cbrewer('seq','Reds',40);
RedsDark=RedsDark(21:40,:);
fmr1_11_reds_dark=RedsDark(bin,:);


%%%% I have been changing this in the original script to add other
%%%% parameters like the edges colors. so some parts wont work anymore...

%%%% i dont manage to change the filling of the markers... so i am plutting
%%%% on top of them another node. 

figure;h=circularGraph(fmr1_11,'ColorMapNode',funct_clust_colors,'ColorMapEdges',fmr1_11_reds_dark);

positions=[];
for i=1:length(h.Node)
positions(i,:)=h.Node(i).Position
end

hold on;
scatter(positions(:,1),positions(:,2),30,funct_clust_colors,'filled');


%%%% making brain labels
RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};

RegionList_av={'Pal','Sp','Th','Hb','Pt','Tec','Tg','Cb','HB'};


Brain_label={};
for i=unique(Nodes4.Mod_brain(keep))'
    
    temp_node=find(Nodes4.Mod_brain(keep)==i);
    
     for j=temp_node'
    Brain_label{j}= RegionList{i}(1:4); 
     end
    
end

Brain_label2={};
for i=unique(Nodes4.Mod_brain(keep))'
    
    temp_node=find(Nodes4.Mod_brain(keep)==i);
    
     for j=temp_node'
    Brain_label2{j}= RegionList_av{i}; 
     end
    
end

figure;h=circularGraph(ctrl_11,'ColorMapNode',black_nodes,'ColorMapEdges',ctrl_11_blues_dark,'Label',Brain_label2);

positions=[];
for i=1:length(h.Node)
positions(i,:)=h.Node(i).Position
end

hold on;
scatter(positions(:,1),positions(:,2),30,funct_clust_colors,'filled');
colormap(BluesDark);
cbh=colorbar;
set(cbh,'XTickLabel',{'0.75','0.8','0.85','0.9','0.95','1'})

figure;h=circularGraph(fmr1_11,'ColorMapNode',black_nodes,'ColorMapEdges',fmr1_11_reds_dark,'Label',Brain_label2);

positions=[];
for i=1:length(h.Node)
positions(i,:)=h.Node(i).Position;
end

hold on;
scatter(positions(:,1),positions(:,2),30,funct_clust_colors,'filled');
colormap(RedsDark);
cbh=colorbar;
set(cbh,'XTickLabel',{'0.75','0.8','0.85','0.9','0.95','1'})


%%% 

%%%% I will try to see how they look order by brain regions. 

%%%% first I need to make it symetrical



%%%% I will try to see how they look order by brain regions. 

[B I]=sort(Nodes4.Mod_brain(keep));

alt_ctrl_11=zeros(size(ctrl_11));

for i=1:length(I)
for j=1:length(I)    

    alt_ctrl_11(i,j)=ctrl_11(I(i),I(j));

end
end


alt_fmr1_11=zeros(size(fmr1_11));

for i=1:length(I)
for j=1:length(I)    

    alt_fmr1_11(i,j)=fmr1_11(I(i),I(j));

end
end

%%%% they look interesting... although the heapmap is probably wrong
figure;h=circularGraph(alt_ctrl_11,'ColorMapNode',black_nodes,'ColorMapEdges',ctrl_11_blues_dark2,'Label',Brain_label(I));

positions=[];
for i=1:length(h.Node)
positions(i,:)=h.Node(i).Position;
end

hold on;
scatter(positions(:,1),positions(:,2),30,funct_clust_colors(I,:),'filled');


figure;h=circularGraph(alt_fmr1_11,'ColorMapNode',black_nodes,'ColorMapEdges',fmr1_11_reds_dark2,'Label',Brain_label(I));

positions=[];
for i=1:length(h.Node)
positions(i,:)=h.Node(i).Position;
end

hold on;
scatter(positions(:,1),positions(:,2),30,funct_clust_colors(I,:),'filled');


%%
%%%% first I need to make it symetrical
%%% finding the midline. I think 310 is pretty good. 
figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([300 1350])
hold on;
gscatter(Nodes4.Mod_loc(keep,1),Nodes4.Mod_loc(keep,2),Nodes4.Mod_brain(keep));
view(-90,90)
plot([500:1200],(ones(length(500:1200))*310),'k');


%%% finding left and right nodes. 

leftH_nodes=find(Nodes4.Mod_loc(keep,2)>310);
rightH_nodes=find(Nodes4.Mod_loc(keep,2)<310);

%%% now i need to see how many i have per brain region in each side and
%%% make it equal
for i=1:length(leftH_nodes)
leftH_nodes(i,2)=Nodes4.Mod_brain(keep(leftH_nodes(i)));
leftH_nodes(i,3)=Nodes4.Mod_clust(keep(leftH_nodes(i)));
end

for i=1:length(rightH_nodes)
rightH_nodes(i,2)=Nodes4.Mod_brain(keep(rightH_nodes(i)));
rightH_nodes(i,3)=Nodes4.Mod_clust(keep(rightH_nodes(i)));
end

Nb_nodes_perBrain=[];
for i=unique(Nodes4.Mod_brain(keep))'

    temp=length(find(leftH_nodes(:,2)==i));
    if isempty(temp)
        temp=0;
    end
    Nb_nodes_perBrain(1,i)=temp;
    
    
     temp=length(find(rightH_nodes(:,2)==i));
    if isempty(temp)
        temp=0;
    end
    Nb_nodes_perBrain(2,i)=temp;
    
end

 nodes_needed=Nb_nodes_perBrain(1,:)- Nb_nodes_perBrain(2,:);
 
 %%%  re-ordering them. left goes ascending and right goes descending
length(leftH_nodes)
[~,I]=sort(leftH_nodes(:,2));
leftH_nodes=leftH_nodes(I,:);


 length(rightH_nodes)
 [~,I]=sort(rightH_nodes(:,2),'descend');
 rightH_nodes=rightH_nodes(I,:);
 
 %%% adding empty nodes to the specific brain regions
new_rightH_nodes=[];
for i=sort(unique(Nodes4.Mod_brain(keep)),'descend')'
    
    add=nodes_needed(i);
    temp1=rightH_nodes(find(rightH_nodes(:,2)==i),:);  
           
    new_rightH_nodes=vertcat(new_rightH_nodes,temp1);
    
    temp2=repmat([0 i 10],add,1);

    new_rightH_nodes=vertcat(new_rightH_nodes,temp2);
    
end



new_brain_order=vertcat(leftH_nodes,new_rightH_nodes);

big_ctrl_11=zeros(length(new_brain_order),length(new_brain_order));
for i=1:length(new_brain_order)
    
   for j=1:length(new_brain_order) 
    
       if new_brain_order(i,1)==0 | new_brain_order(j,1)==0
       
           big_ctrl_11(i,j)=0;
       else
           
           big_ctrl_11(i,j)=ctrl_11(new_brain_order(i,1),new_brain_order(j,1));
           
       end
end

end


big_fmr1_11=zeros(length(new_brain_order),length(new_brain_order));
for i=1:length(new_brain_order)
    
   for j=1:length(new_brain_order) 
    
       if new_brain_order(i,1)==0 | new_brain_order(j,1)==0
       
           big_fmr1_11(i,j)=0;
       else
           
           big_fmr1_11(i,j)=fmr1_11(new_brain_order(i,1),new_brain_order(j,1));
           
       end
end

end


black_nodes_big=zeros(length(new_brain_order),3); %%% for the nodes edges

%%%% making brain labels
Brain_label_big={};
for i=unique(new_brain_order(:,2))'
    
    temp_node=find(new_brain_order(:,2)==i);
    
     for j=temp_node'
    Brain_label_big{j}= RegionList{i}(1:4); 
     end
    
end

%%% creating a functional cluster colormap
funct_clust_colors_big=zeros(length(new_brain_order),3);
for i=unique(new_brain_order(:,3))'
    
    temp_node=find(new_brain_order(:,3)==i);
    
    if i<4    
    funct_clust_colors_big(temp_node,:)= repmat([0 1 0],size(temp_node,1),1);  
    elseif i==4      
     funct_clust_colors_big(temp_node,:)= repmat([0 0 1],size(temp_node,1),1);          
    elseif i==5
     funct_clust_colors_big(temp_node,:)= repmat([1 0 0],size(temp_node,1),1);      
    elseif i==6
     funct_clust_colors_big(temp_node,:)= repmat([1 0 1],size(temp_node,1),1);  
    elseif i==10
     funct_clust_colors_big(temp_node,:)= repmat([0 0 0],size(temp_node,1),1);   
    end   
    
end

%%%% for the edges colormap
[D,edges,bin] = histcounts(nonzeros(big_ctrl_11),20);
[BluesDark]=cbrewer('seq','Blues',40);
BluesDark=BluesDark(21:40,:);
ctrl_11_blues_dark2=BluesDark(bin,:);

[~,~,bin] = histcounts(nonzeros(big_fmr1_11),20);
[RedsDark]=cbrewer('seq','Reds',40);
RedsDark=RedsDark(21:40,:);
fmr1_11_reds_dark2=RedsDark(bin,:);



figure;h=circularGraph(big_ctrl_11,'ColorMapNode',black_nodes_big,'ColorMapEdges',ctrl_11_blues_dark2,'Label',Brain_label_big);

positions=[];
for i=1:length(h.Node)
positions(i,:)=h.Node(i).Position;
end

hold on;
scatter(positions(:,1),positions(:,2),30,funct_clust_colors_big,'filled');



figure;h=circularGraph(big_fmr1_11,'ColorMapNode',black_nodes_big,'ColorMapEdges',fmr1_11_reds_dark2,'Label',Brain_label_big);

positions=[];
for i=1:length(h.Node)
positions(i,:)=h.Node(i).Position;
end

hold on;
scatter(positions(:,1),positions(:,2),30,funct_clust_colors_big,'filled');

%% making the circos figures brain location for looms 1-3,10 and 11


edges2=0.75:0.0125:1;
for L=[1]%[1 2 3 10 11]
    
    %%% for ctrls
   for g= [3 2] 
    group=groupnames{g,1}; 
    loom=fieldnames(MatAll_corrected2.(group));
    
    temp_mat=MatAll_corrected2.(group).(loom{L}).Mat;
    temp_mat(find(isnan(temp_mat)))= 0;
    
    big_temp_mat=zeros(length(new_brain_order),length(new_brain_order));
    for ii=1:length(new_brain_order)    
        for j=1:length(new_brain_order)     
            if new_brain_order(ii,1)==0 | new_brain_order(j,1)==0     
           big_temp_mat(ii,j)=0;
            else          
           big_temp_mat(ii,j)=temp_mat(new_brain_order(ii,1),new_brain_order(j,1));         
            end
        end
    end
    
    [~,~,bin] = histcounts(nonzeros(big_temp_mat),edges2);
    
    if g==3
    [TempColor]=cbrewer('seq','Blues',40);
    else
     [TempColor]=cbrewer('seq','Reds',40);   
    end
    
    TempColor=TempColor(21:40,:);
    TempColor_dark=[];
    TempColor_dark=TempColor(bin,:);
    
    %figure('Renderer', 'painters', 'Position', [100 100 650 650]);
    figure('Position', [100 100 650 650]);
    h=circularGraph(big_temp_mat,'ColorMapNode',black_nodes_big,'ColorMapEdges',TempColor_dark,'Label',Brain_label_big);

    positions=[];
        for p=1:length(h.Node)
        positions(p,:)=h.Node(p).Position;
        end

    hold on;
    scatter(positions(:,1),positions(:,2),30,funct_clust_colors_big,'filled');
    hold off;
    
     %set(gcf, 'Renderer','painters');
    %saveas(gcf,strcat('circularGraph_brain_',group,'_',loom{L},'.svg'));
    %saveas(gcf,strcat('circularGraph_brain_',group,'_',loom{L},'.emf'));
    
   end
end


%%

%%% for the graphs based on functional clustering again... just to see how
%%% they look. 




edges2=0.75:0.0125:1;
for L=[2 3 10 11]%[1 2 3 10 11]
    
    %%% for ctrls
   for g= [3 2] 
    group=groupnames{g,1}; 
    loom=fieldnames(MatAll_corrected2.(group));
    
    temp_mat=MatAll_corrected2.(group).(loom{L}).Mat;
    temp_mat(find(isnan(temp_mat)))= 0;
         
    [~,~,bin] = histcounts(nonzeros(temp_mat),edges2);
    
    if g==3
    [TempColor]=cbrewer('seq','Blues',40);
    else
     [TempColor]=cbrewer('seq','Reds',40);   
    end
    
    TempColor=TempColor(21:40,:);
    TempColor_dark=[];
    TempColor_dark=TempColor(bin,:);
    
    figure('Renderer', 'painters', 'Position', [100 100 550 550]);
    %figure('Position', [100 100 550 550]);
    h=circularGraph(temp_mat,'ColorMapNode',black_nodes,'ColorMapEdges',TempColor_dark,'Label',Brain_label2);

    positions=[];
        for p=1:length(h.Node)
        positions(p,:)=h.Node(p).Position;
        end

    hold on;
    scatter(positions(:,1),positions(:,2),30,funct_clust_colors,'filled');
    hold off;
    
    %set(gcf, 'Renderer','painters');
    saveas(gcf,strcat('circularGraph_cluster_',group,'_',loom{L},'_av','.svg'));
    %saveas(gcf,strcat('circularGraph_cluster_',group,'_',loom{L},'_av','.emf'));
   end
end
