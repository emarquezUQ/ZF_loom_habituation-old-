
%%%%% making circulargraphs without the names on the edges. based on Dani's
%%%%% comment. I do the first part of 
%%%%% testing_graph_analysis_with_BTC3_forFs_figure_fmr1.m and from
%%%%% testing_graph_analysis_with_BTC3_forFs_figure_fmr1_extra1.m

%%% then, I had to comment the part on the circulargraph code to take away
%%% the names... didnt work. but I could do a empty label cell array.

emptyLabel={};
for i=1:90
emptyLabel{1,i}=' ';
end

%%% it worked

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
    h=circularGraph(temp_mat,'ColorMapNode',black_nodes,'ColorMapEdges',TempColor_dark,'Label',emptyLabel);
    %h=circularGraph(temp_mat,'ColorMapNode',black_nodes,'ColorMapEdges',TempColor_dark);

    
    positions=[];
        for p=1:length(h.Node)
        positions(p,:)=h.Node(p).Position;
        end

    hold on;
    scatter(positions(:,1),positions(:,2),30,funct_clust_colors,'filled');
    hold off;
    
    %set(gcf, 'Renderer','painters');
    %saveas(gcf,strcat('circularGraph_cluster_',group,'_',loom{L},'_av_nolabel','.svg'));
    saveas(gcf,strcat('circularGraph_cluster_',group,'_',loom{L},'_av_nolabel','.emf'));
   end
end



