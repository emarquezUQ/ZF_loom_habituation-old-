%%%% this script is to further analyse and plot the distribution of max
%%%% responses to the looms per cluster. It continues from things done in
%%%% the part3analysis_fmr1loomhab.m script.

%%%% one of the things I will try to do is to normalize the area under the
%%%% curve of the distributions so we can compare them even if they have
%%%% different number of ROIs (as fmr1 have less loom respondent ROIs).


load('s20_fmr1_loomhab_CN_part3.mat','GoodClust_goodmembers_full','idx_temp1','idx_temp2','idx_temp3','idx_temp4','idx_temp5');

edges=[0.25:0.25:15];

%%


%%% this part is to get some of usual clusters and get the means of the
%%% distribution to plot them in the right order. 

group = fieldnames(GoodClust_goodmembers_full(i).MaxCountPerFish);
counter=1;

 figure;
for i=[3 7 8 6 5 4 9 10] %%% I am not taking the first 2 clusters. 
   
for loom=1:6
subplot(8,6,counter);
    for k=1:3
     
     plot(edges,GoodClust_goodmembers_full(i).MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean); 
     hold on;
  
    end
     hold off;
     counter=counter+1;
    
end
end


%%

group = fieldnames(GoodClust_goodmembers_full(i).MaxCountPerFish);
counter=1;

figure;
for i=[3 7 8 6 5 4 9 10] %%% I am not taking the first 2 clusters. 
   
for loom=1:6
subplot(8,6,counter);
    for k=1:3
     
     temp_area=trapz(edges,GoodClust_goodmembers_full(i).MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean);
     temp_mean=GoodClust_goodmembers_full(i).MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean;
     temp_mean_norm=temp_mean./temp_area;
     %trapz(edges, temp_mean_norm) %%% to check if it work. the answer
     %should be 1.
     if k==1
       alpha=15;  
     elseif k==2
         alpha=0;
     else
         alpha=10;
     end
     
     plot(edges,temp_mean_norm,'linewidth',1.5,'Color',[0 0 0]+0.05*alpha);
     hold on;
     
    end   
    
    hold off;
     counter=counter+1;
    
end
end


%%
%%% just for the first 3 looms and the 11th

group = fieldnames(GoodClust_goodmembers_full(i).MaxCountPerFish);
counter=1;

figure;
for i=[3 7 8 6 5 4 9 10] %%% I am not taking the first 2 clusters. 
   
for loom=[1 2 3 6]
subplot(8,4,counter);
    for k=1:3
     
     temp_area=trapz(edges,GoodClust_goodmembers_full(i).MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean);
     temp_mean=GoodClust_goodmembers_full(i).MaxCountPerFish.(group{k}).(strcat('loom',num2str(loom))).mean;
     temp_mean_norm=temp_mean./temp_area;
     %trapz(edges, temp_mean_norm) %%% to check if it work. the answer
     %should be 1.
     if k==1
       alpha=15;  
     elseif k==2
         alpha=0;
     else
         alpha=10;
     end
     
     plot(edges,temp_mean_norm,'linewidth',1.5,'Color',[0 0 0]+0.05*alpha);
     hold on;
     
    end   
    
    hold off;
     counter=counter+1;
    
end
end
