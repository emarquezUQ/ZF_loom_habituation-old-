%%% this part is to analyze the tail movements


%%% I will make a new folder to save the graphs
%destdirectory = strcat('/QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab','/f20_mov_graphs_HPC','/');
destdirectory1 = strcat('C:/Emm_temp/FvsS_20vs60_analysis','/f20_mov_graphs','/');
mkdir(destdirectory1);   %create the directory


%%%%%%%this is to look at the tail movements
%%%only works with slooms fish. above fish20
load allfish_f20looms_tailmov




allfish_f20looms_tailmov=cell2mat(allfish_f20looms_tailmov(2:end,:));


%%%first to generate the regressors for the movements. 
Fish_movements={}; %%%first we make a structure to put the data in
for fish=1:size(allfish_f20looms_tailmov,2) %%%we make a variable 'fish' as long as the columns from the Behavioural_responses variable which has the data copied from excel
    for i=1:2 %%%% Note: here I had a loop to work with the different strenght responses but in the new analysis i am not doing that.
        %%% cause I have resposnes from 1-4 depending on the strenght of the moment.
        
        temp=find(allfish_f20looms_tailmov(:,fish)==i); %%%to find the type of responses in turn. 
        temp=round(temp/5);%%%% i this is to round it to the frame rate that we use for the SPIM videos. the tail movies have a 10 Hz and the SPIM is 2Hz. 
        temp(temp<1)=[];
        temp(temp>(length(allfish_f20looms_tailmov)-1))=[];        
        Fish_movements{fish,i}=temp;
        
        
    end
    
    %%% this is loop is to locate the moments the loom is ON.
    %%% but I need to make a vector with the momennts the loom is on. 
    if fish==13 %%at the moment is 13 (for f20)
        temp=find(allfish_f20looms_tailmov(:,fish)==2);
        temp=round(temp/5);
        temp(temp<1)=[];
        temp(temp>(length(allfish_f20looms_tailmov)-1))=[];        
        Fish_movements{fish,1}=temp;
    end
end

%%%here is for the lenght of the SPIM movies

%ZS=zeros(1344,1); ZS=ZS'; %% for f20
%ZS=zeros(3504,1); ZS=ZS'; %% for f60


%%%this is to put Gcamp spikes to the fish movements. 
Movements=zeros(size(allfish_f20looms_tailmov,2)-1,1,size(ZS_f20,2)); %%%to mae a 3D array the size of the number of fish (18 at the moment), the type of movement (1) and the lenght of the SPIM movies (depending on the variable). 
%GCaMP6=[5.13796058542217,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502]';
GCaMP6=[0,0.5,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
for fish=1:size(allfish_f20looms_tailmov,2)-1
    for movement=1  %movement=1:4
        idx=Fish_movements{fish,movement};%%%%to take the timepoints of the movements.  
        if idx
            for i=idx'
                Movements(fish,movement,i:i+size(GCaMP6,1)-1)=GCaMP6/6;%%% to place the Gcamp spikes in the places where the movements occur
            end
        end
    end
end
Movements=Movements(:,:,1:size(ZS_f20,2));

%%

%%% this is the fish list. but I want find a way to put them on top of
%%% every fish behavioural gcamp trace

% Fish_with_behav_f20=[1 7 13 17 21 25 29 34 38 40 44 48];
% 
% Fish_with_behav_f60=[3 4 8 9 15 19 26 32 37 41 47];
% 
% Fish_with_behav_s20=[5 10 14 18 22 24 28 33 36 43 45];
% 
% Fish_with_behav_s60=[6 11 12 16 20 23 27 30 31 35 39 42 46];



%%%%I think this is to correlate the Gcamp spikes generated with the movements to the
%%%%Gcamps of the actual brain imaging of each fish. 
Correlation_movement={};
counter=1;
for fish=Fish_with_behav_f20
    temp_ZS=ZS_f20(find(idx_Fish_f20==fish),:);Correlation=[];
    for idx=1:size(temp_ZS,1)
        for movement=1  %:4
            temp=corrcoef(squeeze(squeeze(Movements(counter,movement,:))),temp_ZS(idx,:));
            Correlation(movement,idx)=temp(1,2);            
        end
    end
    Correlation_movement{counter}=Correlation;
    counter=counter+1;
end


%%% to show where when all the looms are presented


loomtimes=zeros(1,size(ZS_f20,2)); %%%to mae a 3D array the size of the number of fish (18 at the moment), the type of movement (1) and the lenght of the SPIM movies (depending on the variable). 

idx=Fish_movements{13,1};%%%%to take the timepoints of the movements.  
       
  for i=idx'
   loomtimes(1,i)=0.5;%%% to place the Gcamp spikes in the places where the movements occur
  end
        
 


figure;plot(loomtimes);




%%%%% this is to plot the mean of the ROIs that correlate to the looms and
%%%%% the movements above the chosen Threshold. the plot has 4 subplots, each for a kind of movement. 
%%%in this part of the script and some others further down I added an if
%%%loop to make it work when there are less fish than in the
%%%Behavioural_responses matrix.
Threshold=0.5;


for fish=1:length(Fish_with_behav_f20)
    if isempty(Correlation_movement{fish})~=1 % idx_Fish_f20==Fish_with_behav_f20(fish)
    Fighandle=figure;
    temp_ZS=ZS_f20(find(idx_Fish_f20==Fish_with_behav_f20(fish)),:);
    %for 
        i=1
        Correlation=Correlation_movement{fish};
        %subplot(2,2,i);
        plot(mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1));
        hold on
        Movements2=squeeze(Movements(:,i,:));
        Movements2=squeeze(Movements2);
        plot(Movements2(fish,:));
        plot(loomtimes);
        
        %mdl=stepwiselm(Stimuli',mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
        mdl=fitlm(squeeze(squeeze(Movements(counter,movement,:)))',mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1)); %%% I changed it to fitlm cause is faster
        title(strcat('lin_reg to movs : ',num2str([mdl.Rsquared.Adjusted]),' nb of ROIs : ', num2str(length(find(Correlation(i,:)>Threshold)))));
    
    %end
    else
    end
    %print(Fighandle,strcat(destdirectory1,'f20_movs_Fish',num2str(Fish_with_behav_f20(fish))),'-dpng','-r0');
    %close all;
    
end


%%%this is to find the movements that overlap while the looms are ON. 
Correlation_movements_loomsb={};
for fish=1:length(Fish_with_behav_f20)   %%%for the amount of fish recorded. the -1 is because the last one are the looms 
    Mov_temp=Fish_movements{fish,1}; %%%we take the index of their movements.
    %for i=1  %%%2:4  %%%cause we only have 1 movement.
        %Mov_temp=vertcat(Mov_temp,Fish_movements{fish,i}); %%%to do a vertical concatenation of the indexes of the movements acording to the kind of movement. 
    %end
    Correlation_movements_loomsb{fish}=intersect([Fish_movements{13,2}],Mov_temp);%%%%find the movements that happend when the looms was ON. 
end

%%%this is to put Gcamp spikes on the movements that happend while the
%%%looms were ON.
Movements_looms=zeros(size(allfish_f20looms_tailmov,2)-1,size(ZS_f20,2)); %%%to mae a 3D array the size of the number of fish (12 at the moment), the type of movement (1) and the lenght of the SPIM movies (depending on the variable). 

GCaMP6=[0,0.5,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
for fish=1:length(Fish_with_behav_f20)
    if Correlation_movements_loomsb{fish}
        for i=Correlation_movements_loomsb{fish}'
            Movements_looms(fish,i:i+size(GCaMP6,1)-1)=GCaMP6/6; %%% I divide by 6 cause i think the original gcapm is too big
        end
    end
end
Movements_looms=Movements_looms(:,1:size(ZS_f20,2));


%%%%I think this is to correlate the Gcamp spikes generated with the movements that happend while the looms were ON to the
%%%%Gcamps of the actual brain imaging of each fish. 
Correlation_movement_loom={};
counter=1;
for fish=Fish_with_behav_f20
    temp_ZS=ZS_f20(find(idx_Fish_f20==fish),:);Correlation=[];
    for idx=1:size(temp_ZS,1)
        temp=corrcoef(squeeze(squeeze(Movements_looms(counter,:))),temp_ZS(idx,:));
        Correlation(idx)=temp(1,2);        
    end
    Correlation_movement_loom{counter}=Correlation;
    counter=counter+1;
end


%%%%and this is to plot the means of the traces of the ROIs that correlated
%%%%with the loom movements. It makes a figure with a subplot of each fish.
%%%%it also adds a red rectangle where the looms was shown. 
figure;%Threshold=0.5;
counter=1;counter2=1;xplot=floor(sqrt(length(Fish_with_behav_f20)));yplot=ceil(length(Fish_with_behav_f20)/xplot); 
for fish=1:length(Fish_with_behav_f20)    
    temp_ZS=ZS_f20(find(idx_Fish_f20==Fish_with_behav_f20(fish)),:);    
    Corr_temp=Correlation_movement_loom{fish};
    subplot(xplot,yplot,fish);plot(mean(temp_ZS(find(Corr_temp>Threshold),:),1));
%%% this is to put a rectangle where the looms appear but
%%% i need to modify the times for this new dataset
%     for i=1:10  
%         rectangle('FaceColor','r','Position',[50+((i-1)*60) -1 10 0.25]); %%%%this is for the rectangles. it was at 40 but I changed it to 50 cause the looms starts 5 seconds after.
%     end
    
    title(strcat(' nb of ROIs : ', num2str(length(find(Corr_temp>Threshold)))));
end

%%%this is for a raster plot for each fish with the traces of the ROIs that
%%%correlated to the loom movements. 
figure;
xplot=floor(sqrt(length(Fish_with_behav_f20)));yplot=ceil(length(Fish_with_behav_f20)/xplot);
for fish=1:length(Fish_with_behav_f20)    
    temp_ZS=ZS_f20(find(idx_Fish_f20==Fish_with_behav_f20(fish)),:);    
    Corr_temp=Correlation_movement_loom{fish};
    subplot(xplot,yplot,fish);imagesc(temp_ZS(find(Corr_temp>Threshold),:),[-0.5 3]);    
end

 
%%%%and this is to plot the means of the traces of the ROIs that correlated
%%%%with the movements. It makes a figure with a subplot of each fish. 
 
%for 
    i=1 %%%:4
    figure; 
    xplot=floor(sqrt(length(Fish_with_behav_f20)));yplot=ceil(length(Fish_with_behav_f20)/xplot);
    for fish=1:length(Fish_with_behav_f20)%fish=2
        if isempty(Correlation_movement{fish})~=1 % idx_Fish_f20==Fish_with_behav_f20(fish)
        temp_ZS=ZS_f20(find(idx_Fish_f20==Fish_with_behav_f20(fish)),:);        
        Correlation=Correlation_movement{fish};
        subplot(xplot,yplot,fish);plot(mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1));
        hold on
        Movements2=squeeze(Movements(:,i,:));
        Movements2=squeeze(Movements2);
        plot(Movements2(fish,:))
        plot(loomtimes);
        hold off
        
        %mdl=stepwiselm(Stimuli',mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
        %title(num2str([mdl.Rsquared.Adjusted]));
        title(strcat(' nb of ROIs : ', num2str(length(find(Correlation(i,:)>Threshold)))));
        else
        end
    end
%end

%%%this is for a raster plot for each fish with the traces of the ROIs that
%%%correlated to the movements. we take the ZS traces from the ROIs that had a correlation coeficient higher than 0.5. 
 
%for 
    i=1 %%%:4
    Fighandle=figure;
    xplot=floor(sqrt(length(Fish_with_behav_f20)));yplot=ceil(length(Fish_with_behav_f20)/xplot);
    for fish=1:length(Fish_with_behav_f20)
        if isempty(Correlation_movement{fish})~=1 % idx_Fish_f20==Fish_with_behav_f20(fish)
        temp_ZS=ZS_f20(find(idx_Fish_f20==Fish_with_behav_f20(fish)),:);        
        Correlation=Correlation_movement{fish};
        subplot(xplot,yplot,fish);imagesc(temp_ZS(find(Correlation(i,:)>Threshold),:),[-0.5 3]);       
        else
        end  
       
    end
%end

 print(Fighandle,strcat(destdirectory1,'f20_movs_raster_all'),'-dpng','-r0');
    %close all;




%%%this is for a raster plot for each fish with the traces of the ROIs that
%%%correlated to the loom movements but also had a correlation higher than 0.5. before it was with the ones that  passed the crieterion of
%%%coefficient significance and are above the thereshold the r2 value of the linear regression we
%%%chose above. and we add red rectangles when the looms were ON.
figure;%Threshold=0.5;
xplot=floor(sqrt(length(Fish_with_behav_f20)));yplot=ceil(length(Fish_with_behav_f20)/xplot);
for fish=1:length(Fish_with_behav_f20)    
    idx_fish=find(idx_Fish_f20==Fish_with_behav_f20(fish)); %%%this is to take the index of the ROIs of the fish in turn
    temp_ZS=ZS_f20(intersect(idx_coef_rsq,idx_fish),:); %%%this is to take the filtered above (sig beta and above the r2 value that we want)ZS of the ROIs of the fish in turn   
    Corr_temp=Correlation_movement_loom{fish};
    [~,idx_interecpt]=intersect(idx_fish,idx_coef_rsq); %%%this is to get the index of the filtered above (sig beta and above the r2 value that we want)ZS of the ROIs of the fish in turn    
    subplot(xplot,yplot,fish);imagesc(temp_ZS(find(Corr_temp(idx_interecpt)>Threshold),:),[-0.5 3]);
%%% this is for the rectagles when the loom appear but I need to adapt it
%%% to this dataset
%     for j=1:10
%         rectangle('FaceColor','r','Position',[50+((j-1)*60) -1 10 25]); %%%this is  to put a rectangle when the looms appear. Gilles had it at 40 but it should be at 50 so i changed it. 
%     end
end

%%%this is for a raster plot for each fish with the traces of the ROIs that
%%%correlated to the type of movements but also passed the crieterion of
%%%of a simple correlation above the thereshold. and we add red rectangles when the looms were ON.

%for
    i=1  %%%:4
    figure;
    xplot=floor(sqrt(length(Fish_with_behav_f20)));yplot=ceil(length(Fish_with_behav_f20)/xplot);
    for fish=1:length(Fish_with_behav_f20)
        if isempty(Correlation_movement{fish})~=1 % idx_Fish_f20==Fish_with_behav_f20(fish)
        idx_fish=find(idx_Fish_f20==Fish_with_behav_f20(fish));
        temp_ZS=ZS_f20(find(idx_Fish_f20==Fish_with_behav_f20(fish)),:);        
        Correlation=Correlation_movement{fish};
        [~,idx_interecpt]=intersect(idx_fish,idx_corr);    
        subplot(xplot,yplot,fish);imagesc(temp_ZS(find(Correlation(i,idx_interecpt)>Threshold),:),[-0.5 3]);       
%         for j=1:10
%             rectangle('FaceColor','r','Position',[50+((j-1)*60) -1 10 25]); %%%this is  to put a rectangle when the looms appear. Gilles had it at 40 but it should be at 50 so i changed it.
%         end
        else
        end
    end
%end

%%
Threshold2=0.25;
%%%this is to get the idx of the ROIs that correlated to the loom movements

idx_LoomMov={};
counter=1;
for fish=1:length(Fish_with_behav_f20)    
    temp_idx=find(idx_Fish_f20==Fish_with_behav_f20(fish)); idx_LoomMov_temp=[];   
    Corr_temp=Correlation_movement_loom{fish};
    idx_LoomMov_temp=temp_idx(find(Corr_temp>Threshold2),:);
    idx_LoomMov{counter}=idx_LoomMov_temp;
    counter=counter+1;
end

clear idx_LoomMov_temp temp_idx Corr_temp counter

 idx_LoomMov=vertcat(idx_LoomMov{:});


%%% this is for a multigraph of the ROIs that correlated to the loom
%%% movements for individual fish. is not very informative...
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(Fish_with_behav_f20,2);counter=1;
for fish=1:length(Fish_with_behav_f20)
    temp_idx=find(idx_Fish_f20==Fish_with_behav_f20(fish)); idx_LoomMov_temp=[];   
    Corr_temp=Correlation_movement_loom{fish};
    idx_LoomMov_temp=temp_idx(find(Corr_temp>Threshold),:);
    %idx_temp=find(High_corr_Nb_all==i);
    subplot(rows,4,counter);plot(mean(ZS_f20(idx_LoomMov_temp,:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f20(idx_LoomMov_temp,:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f20(idx_LoomMov_temp)); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f20(idx_LoomMov_temp));%%% for the fish location
    counter=counter+4;
end

%saveas(gcf,'multigraph_LoomMov_CN_F20','jpg');


%%%this is to get the idx of the ROIs that correlated to the movements
%%%(regardless of the loom)


idx_Mov={};
counter=1;
for fish=1:length(Fish_with_behav_f20)    
    temp_idx=find(idx_Fish_f20==Fish_with_behav_f20(fish)); idx_Mov_temp=[];   
    Corr_temp=Correlation_movement{fish};
    idx_Mov_temp=temp_idx(find(Corr_temp>Threshold),:);
    idx_Mov{counter}=idx_Mov_temp;
    counter=counter+1;
end

clear idx_Mov_temp temp_idx Corr_temp counter

 idx_Mov=vertcat(idx_Mov{:});

 
 %%% this is for a multigraph of the ROIs that correlated to the loom
%%% movements for individual fish. is not very informative...

 

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(Fish_with_behav_f20,2);counter=1;
for fish=1:length(Fish_with_behav_f20)
    temp_idx=find(idx_Fish_f20==Fish_with_behav_f20(fish)); idx_Mov_temp=[];   
    Corr_temp=Correlation_movement{fish};
    idx_Mov_temp=temp_idx(find(Corr_temp>Threshold),:);
    %idx_temp=find(High_corr_Nb_all==i);
    subplot(rows,4,counter);plot(mean(ZS_f20(idx_Mov_temp,:),1)); %%%to plot the mean
    hold on %%% this bit is also to plot the mov regressor
%         Movements2=squeeze(Movements(:,i,:)); %%% done before
%         Movements2=squeeze(Movements2);
        plot(Movements2(fish,:))
        hold off
    subplot(rows,4,counter+1);imagesc(ZS_f20(idx_Mov_temp,:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f20(idx_Mov_temp)); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f20(idx_Mov_temp));%%% for the fish location
    counter=counter+4;
end

%saveas(gcf,'multigraph_Mov_CN_F20','jpg');


idx_good_LoomMov=intersect(idx_LoomMov,idx_Mov);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% this part is to make the Kmeans figures with the localization of the selected clusters in the brain  %%%%%%%%%

temp=[];
counter=1;
temp{counter}=idx_LoomMov;    
counter=counter+1;    
temp{counter}=idx_Mov; 
counter=counter+1;    
temp{counter}=idx_good_LoomMov; 


Numbers(1)=0; %%% to change the first value of Numbers to 0 (cause it was not relevant)
%colors = [1,0,0;0,1,0;0,0,1;1,0.600000000000000,0;1,0,0.7];
colors = distinguishable_colors(length(temp),[1 1 1; 0 0 0]); %%%here we use a script from Matlab (downloaded, and needs to be in the folder) to generate colors
% colors = [0         0    1.0000
%          0    0.5000    1.0000
%     1.0000         0         0
%     1.0000    0.1034    0.7241
%     1.0000    0.5000    0.3000
%          0    0.7000    0.2000
%     0.5000    0.5000         0
%          0    0.5000    0.5000];
colors = colors*256;


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;counter2=1;xplot=floor(sqrt(length(temp)));yplot=ceil(length(temp)/xplot); rows=length(temp);
for i=1:length(temp)
    idx_temp=temp{i};
    subplot(rows,4,counter);plot(mean(ZS_f20(idx_temp,:),1),'color',colors(counter2,:)/256); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f20(idx_temp,:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f20(idx_temp)); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f20(idx_temp));%%% for the fish location
    %counter=counter+1;
    counter2=counter2+1
     counter=counter+4;
end

%saveas(gcf,'multigraph_LoomMov_Mov_CN_F20_R04','jpg');
 print(Fighandle,strcat(destdirectory1,'multigraph_LoomMov_Mov_CN_F20_R04'),'-dpng','-r0');
    %close all;


idx=1
    filename=MatFiles(idx).name;
%destdirectory = strcat('/QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab','/Kmeans_f20_r2050_CL50_HPC',filename(8:11),'/');
destdirectory = strcat('C:/Emm_temp/FvsS_20vs60_analysis','/Kmeans_f20_CN_good_LoomMov_R04',filename(8:11),'/');
mkdir(destdirectory);   %create the directory
clearvars idx filename

for idx=1:length(MatFiles)
    filename=MatFiles(idx).name;%%%to put the respective name of the files (in this case the slices)
    ROIsNb=[];ClusterNb=[];%%% make the variables we are going to use
    %for k = 1 : length(temp)
    %%%%%this is to locate in every plane which ROIs and cluster can be
    %%%%%found
    for k = 3 %1 : length(temp)
        tempROIsNb=find([temp{k}]<=Numbers(idx+1) & [temp{k}]>Numbers(idx)); %%%to put the data filtered of the selected clusters on a new variable. but i dont understand the sentece... why not just ask for Numbers(idx)?
        if tempROIsNb            
            ROIsNb=[ROIsNb; temp{k}(tempROIsNb)];%%%%
            %temp{k}(tempROIsNb)=[];
            ClusterNb=[ClusterNb; repmat(k,length(tempROIsNb),1)];
        end
    end
    
    %%%% and this part is to make and image using the _mean.tif image as
    %%%% base, where we will add the ROIs located before and colored base
    %%%% on the number of the cluster
    if ROIsNb
        imagename=regexp(filename,'_output_analysis','split');
        %imagename=regexp(imagename,'_output_analysis_matlab2.mat','split');
        imagename=strcat(imagename{1},'_mean.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*64;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=MatFiles(idx).ROI;       
        ROIsNb=ROIsNb-Numbers(idx);
        ROIs=ROIs(:,ROIsNb);
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
            for j=1:3
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
        %image3(:,:,3)=image;
        %%%finally we save the Image
        name=strcat(destdirectory,'_Kmeans_f20_CN_good_LoomMov_R04',imagename(4:end));
    imwrite(image3,name,'tif');
    end
    %image3=uint8(image3);

end


clearvars idx i tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster


%%
%%% i am going to look where are the neuros that that responded when
%%% movements were happening during looms

%%% not sure if this next step is necessary

%%%% to correlate the Gcamp spikes generated with the movements that happend while the looms were ON to the
%%%%Gcamps of the actual brain imaging of each fish.I am taking away the
%%%%first 3 looms. for f20 it would be from 180 onwards
Correlation_movement_loom_post3={};
counter=1;
for fish=Fish_with_behav_f20
    temp_ZS=ZS_f20(find(idx_Fish_f20==fish),180:end);Correlation=[];
    for idx=1:size(temp_ZS,1)
        temp=corrcoef(squeeze(squeeze(Movements_looms(counter,180:end))),temp_ZS(idx,:));
        Correlation(idx)=temp(1,2);        
    end
    Correlation_movement_loom_post3{counter}=Correlation;
    counter=counter+1;
end

%%%%and this is to plot the means of the traces of the ROIs that correlated
%%%%with the loom movements. It makes a figure with a subplot of each fish.
%%%%it also adds a red rectangle where the looms was shown. 
figure;%Threshold=0.5;
counter=1;counter2=1;xplot=floor(sqrt(length(Fish_with_behav_f20)));yplot=ceil(length(Fish_with_behav_f20)/xplot); 
for fish=1:length(Fish_with_behav_f20)    
    temp_ZS=ZS_f20(find(idx_Fish_f20==Fish_with_behav_f20(fish)),:);    
    Corr_temp=Correlation_movement_loom_post3{fish};
    subplot(xplot,yplot,fish);plot(mean(temp_ZS(find(Corr_temp>Threshold),:),1));
%%% this is to put a rectangle where the looms appear but
%%% i need to modify the times for this new dataset
%     for i=1:10  
%         rectangle('FaceColor','r','Position',[50+((i-1)*60) -1 10 0.25]); %%%%this is for the rectangles. it was at 40 but I changed it to 50 cause the looms starts 5 seconds after.
%     end
    
    title(strcat(' nb of ROIs : ', num2str(length(find(Corr_temp>Threshold)))));
end


%%%this is for a raster plot for each fish with the traces of the ROIs that
%%%correlated to the loom movements but also had a correlation higher than 0.5. before it was with the ones that  passed the crieterion of
%%%coefficient significance and are above the thereshold the r2 value of the linear regression we
%%%chose above. and we add red rectangles when the looms were ON.
figure;%Threshold=0.5;
xplot=floor(sqrt(length(Fish_with_behav_f20)));yplot=ceil(length(Fish_with_behav_f20)/xplot);
for fish=1:length(Fish_with_behav_f20)    
    idx_fish=find(idx_Fish_f20==Fish_with_behav_f20(fish)); %%%this is to take the index of the ROIs of the fish in turn
    temp_ZS=ZS_f20(intersect(idx_coef_rsq,idx_fish),:); %%%this is to take the filtered above (sig beta and above the r2 value that we want)ZS of the ROIs of the fish in turn   
    Corr_temp=Correlation_movement_loom_post3{fish};
    [~,idx_interecpt]=intersect(idx_fish,idx_coef_rsq); %%%this is to get the index of the filtered above (sig beta and above the r2 value that we want)ZS of the ROIs of the fish in turn    
    subplot(xplot,yplot,fish);imagesc(temp_ZS(find(Corr_temp(idx_interecpt)>Threshold),:),[-0.5 3]);
%%% this is for the rectagles when the loom appear but I need to adapt it
%%% to this dataset
%     for j=1:10
%         rectangle('FaceColor','r','Position',[50+((j-1)*60) -1 10 25]); %%%this is  to put a rectangle when the looms appear. Gilles had it at 40 but it should be at 50 so i changed it. 
%     end
end


%%%this is to get the idx of the ROIs that correlated to the loom movements

idx_LoomMov_post3={};
counter=1;
for fish=1:length(Fish_with_behav_f20)    
    temp_idx=find(idx_Fish_f20==Fish_with_behav_f20(fish)); idx_LoomMov_temp=[];   
    Corr_temp=Correlation_movement_loom_post3{fish};
    idx_LoomMov_temp=temp_idx(find(Corr_temp>Threshold2),:);
    idx_LoomMov_post3{counter}=idx_LoomMov_temp;
    counter=counter+1;
end

clear idx_LoomMov_temp temp_idx Corr_temp counter

 idx_LoomMov_post3=vertcat(idx_LoomMov_post3{:});


%%
%%% this is if i try the same but with the whole mov regressors. 


%%%%I think this is to correlate the Gcamp spikes generated with the movements to the
%%%%Gcamps of the actual brain imaging of each fish. 
Correlation_movement_post3={};
counter=1;
for fish=Fish_with_behav_f20
    temp_ZS=ZS_f20(find(idx_Fish_f20==fish),180:end);Correlation=[];
    for idx=1:size(temp_ZS,1)
        for movement=1  %:4
            temp=corrcoef(squeeze(squeeze(Movements(counter,movement,180:end))),temp_ZS(idx,:));
            Correlation(movement,idx)=temp(1,2);            
        end
    end
    Correlation_movement_post3{counter}=Correlation;
    counter=counter+1;
end

%%%%and this is to plot the means of the traces of the ROIs that correlated
%%%%with the movements. It makes a figure with a subplot of each fish. 
 
%for 
    i=1 %%%:4
    figure; 
    xplot=floor(sqrt(length(Fish_with_behav_f20)));yplot=ceil(length(Fish_with_behav_f20)/xplot);
    for fish=1:length(Fish_with_behav_f20)%fish=2
        if isempty(Correlation_movement_post3{fish})~=1 % idx_Fish_f20==Fish_with_behav_f20(fish)
        temp_ZS=ZS_f20(find(idx_Fish_f20==Fish_with_behav_f20(fish)),:);        
        Correlation=Correlation_movement_post3{fish};
        subplot(xplot,yplot,fish);plot(mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1));
        hold on
%       Movements2=squeeze(Movements(:,i,:)); %%% this was already done before
%       Movements2=squeeze(Movements2);
        plot(Movements2(fish,:))
        plot(loomtimes);
        hold off
        
        %mdl=stepwiselm(Stimuli',mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
        %title(num2str([mdl.Rsquared.Adjusted]));
        title(strcat(' nb of ROIs : ', num2str(length(find(Correlation(i,:)>Threshold)))));
        else
        end
    end
%end


%%%this is for a raster plot for each fish with the traces of the ROIs that
%%%correlated to the movements but also had a correlation higher than 0.5. before it was with the ones that  passed the crieterion of
%%%coefficient significance and are above the thereshold the r2 value of the linear regression we
%%%chose above. and we add red rectangles when the looms were ON.
figure;%Threshold=0.5;
xplot=floor(sqrt(length(Fish_with_behav_f20)));yplot=ceil(length(Fish_with_behav_f20)/xplot);
for fish=1:length(Fish_with_behav_f20)    
    idx_fish=find(idx_Fish_f20==Fish_with_behav_f20(fish)); %%%this is to take the index of the ROIs of the fish in turn
    temp_ZS=ZS_f20(intersect(idx_rsq_test,idx_fish),:); %%%this is to take the filtered above (sig beta and above the r2 value that we want)ZS of the ROIs of the fish in turn   
    Corr_temp=Correlation_movement_post3{fish};
    [~,idx_interecpt]=intersect(idx_fish,idx_rsq_test); %%%this is to get the index of the filtered above (sig beta and above the r2 value that we want)ZS of the ROIs of the fish in turn. i am testing with 0.3r2    
    subplot(xplot,yplot,fish);imagesc(temp_ZS(find(Corr_temp(idx_interecpt)>Threshold),:),[-0.5 3]);
%%% this is for the rectagles when the loom appear but I need to adapt it
%%% to this dataset
%     for j=1:10
%         rectangle('FaceColor','r','Position',[50+((j-1)*60) -1 10 25]); %%%this is  to put a rectangle when the looms appear. Gilles had it at 40 but it should be at 50 so i changed it. 
%     end
end


%%%this is to get the idx of the ROIs that correlated to the movements
%%%(regardless of the loom)

idx_Mov_post3={};
counter=1;
for fish=1:length(Fish_with_behav_f20)    
    temp_idx=find(idx_Fish_f20==Fish_with_behav_f20(fish)); idx_Mov_temp=[];   
    Corr_temp=Correlation_movement_post3{fish};
    idx_Mov_temp=temp_idx(find(Corr_temp>Threshold),:);
    idx_Mov_post3{counter}=idx_Mov_temp;
    counter=counter+1;
end

clear idx_Mov_temp temp_idx Corr_temp counter

 idx_Mov_post3=vertcat(idx_Mov_post3{:});

idx_good_LoomMov_post3=intersect(idx_LoomMov_post3,idx_Mov_post3);

%%% this is to get the idx of filtered by rsq and coef ROIs but that
%%% correlated to movements post 3 looms. 



idx_filtered_good_LoomMov_post3=intersect(idx_rsq_test,idx_Mov_post3); %%% i am using a thereshold of the regression of 0.3 to see if it improves.

figure;imagesc(ZS_f20(idx_filtered_good_LoomMov_post3,:),[-0.5 4]);colormap hot

figure;imagesc(ZS_f20(idx_good_LoomMov_post3,:),[-0.5 4]);colormap hot

figure;imagesc(ZS_f20(idx_Mov_post3,:),[-0.5 4]); colormap hot

figure;imagesc(ZS_f20(idx_Mov,:),[-0.5 4]); colormap hot

%%
%%% testing with a linear regression!!

%%% i will check how it looks with a linear regression with the movements that happend while the looms were ON to the
%%%%Gcamps of the actual brain imaging of each fish.I am taking away the
%%%%first 3 looms. for f20 it would be from 180 onwards

%%% it doenst look very good... 

ModelResults_movPerFish={};
counter=1;
rsquare_loom_mov={};
for fish= Fish_with_behav_f20
   if ismember(fish,unique(idx_Fish_f20))  
    temp_ZS=ZS_f20(find(idx_Fish_f20==fish),180:end);ModelResults_temp=[];
    
        parfor i=1:size(temp_ZS,1)  %%%parfor is to do it in parallel
    
        %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
        %%vector Y as a response variable and the columns of the matrix X as
        %%%predictor variables, performs stepwise regression, and returns the
        %%final result as the linear model LM.
        mdl=fitlm(squeeze(squeeze(Movements_looms(counter,180:end)))',temp_ZS(i,:));
    
        %%%this is to put the results in the ModelResulsts variable
        ModelResults_temp(i).coef=mdl.Coefficients;
        ModelResults_temp(i).MSE=mdl.MSE;
        ModelResults_temp(i).Fitted=mdl.Fitted;
        ModelResults_temp(i).rsquared=mdl.Rsquared.Adjusted;
        end
        
        ModelResults_movPerFish{counter}=ModelResults_temp;
       
        rsquare_loom_mov{counter}=[ModelResults_temp.rsquared];
        
        counter=counter+1;
        
   else
   end
end

rsquare_loom_mov=horzcat(rsquare_loom_mov{:});
rsquare_loom_mov=rsquare_loom_mov';

figure;histogram(rsquare_loom_mov);

%%%not working due to bad indexing
idx_rsq_move=find(rsquare_loom_mov>0.1 & rsquare_loom_mov<1); %%%then select the rsquare that are between 0.3 and 1
figure; imagesc(ZS_f20(idx_rsq_move,:), [-0.5 4]);colormap hot


idx_rsq_LoomMov={};
rsquare_loom_loomMov={};
counter=1;
for fish= Fish_with_behav_f20
   if ismember(fish,unique(idx_Fish_f20))  
    temp_idx=find(idx_Fish_f20==fish); idx_Mov_temp=[];   
    rsq_temp=[ModelResults_movPerFish{1,counter}.rsquared];
    idx_Mov_temp=temp_idx(find(rsq_temp>0.1 & rsq_temp<1),:);
    idx_rsq_LoomMov{counter}=idx_Mov_temp;
    rsquare_loom_loomMov{counter}=rsq_temp;
    counter=counter+1;
   else
   end    
end

clear idx_Mov_temp temp_idx Corr_temp counter

idx_rsq_LoomMov=vertcat(idx_rsq_LoomMov{:});
figure; imagesc(ZS_f20(idx_rsq_LoomMov,:), [-0.5 4]);colormap hot





%%% i will check how it looks with a linear regression with the all the movements to the
%%%%Gcamps of the actual brain imaging of each fish.I am taking away the
%%%%first 3 looms. for f20 it would be from 180 onwards




%%% I might aslo try to select only the
%%% parts where the movement responses could be happening. 

ModelResults_allmovPerFish={};
counter=1;
rsquare_loom_allmov={};
for fish= Fish_with_behav_f20
   if ismember(fish,unique(idx_Fish_f20))  
    temp_ZS=ZS_f20(find(idx_Fish_f20==fish),:);ModelResults_temp=[];
    %temp_ZS=ZS_f20(find(idx_Fish_f20==fish),:);ModelResults_temp=[];
        parfor i=1:size(temp_ZS,1)  %%%parfor is to do it in parallel
    
        %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
        %%vector Y as a response variable and the columns of the matrix X as
        %%%predictor variables, performs stepwise regression, and returns the
        %%final result as the linear model LM.
        
        mdl=fitlm(squeeze(squeeze(Movements(find(Fish_with_behav_f20==fish),:)))',temp_ZS(i,:));
         
        %mdl=fitlm(rawregress_CN',temp_ZS(i,:)); %%to test if the code is working
         
        %%% this one is to take just the timepoints where the looms should
        %%% be on... i might need to fidle with it to add extra timepoints
        %%% before and after.
        %mdl=fitlm(squeeze(squeeze(Movements(find(Fish_with_behav_f20==fish),find(Movements2(find(Fish_with_behav_f20==fish),:)))))',temp_ZS(i,find(Movements2(find(Fish_with_behav_f20==fish),:))));
    
        
        %%%this is to put the results in the ModelResulsts variable
        ModelResults_temp(i).coef=mdl.Coefficients;
        ModelResults_temp(i).MSE=mdl.MSE;
        ModelResults_temp(i).Fitted=mdl.Fitted;
        ModelResults_temp(i).rsquared=mdl.Rsquared.Adjusted;
        end
        
        ModelResults_allmovPerFish{counter}=ModelResults_temp;
       
        rsquare_loom_allmov{counter}=[ModelResults_temp.rsquared];
        
         counter=counter+1;
        
   else
   end
end

%%% this was to test some things. cause it work for individual fish
%rsquare_loom_allmov=[];
%rsquare_loom_allmov = [ModelResults_allmovPerFish{1,11}.rsquared]; %%% to
%select from a specific fish
%idx_rsq_allmove=find(rsquare_loom_allmov>0.15 & rsquare_loom_allmov<1); %%%then select the rsquare that are between 0.15 and 1
%figure; imagesc(temp_ZS(idx_rsq_allmove,:), [-0.5 4]);colormap hot

%%% this didnt work due to the indexing
% rsquare_loom_allmov=horzcat(rsquare_loom_allmov{:});
% figure;histogram(rsquare_loom_allmov);
% idx_rsq_allmove=find(rsquare_loom_allmov>0.15 & rsquare_loom_allmov<1); %%%then select the rsquare that are between 0.3 and 1
% figure; imagesc(ZS_f20(idx_rsq_allmove,:), [-0.5 4]);colormap hot


%%%this is to get the idx of the ROIs that passed the rsq2 thereshold that we choose to the movements
%%%(regardless of the loom). it solves the problem of the indexing. 
%%% i think 0.2 is a good threshold


idx_rsq_Mov={};
rsquare_loom_allmov2={};
counter=1;
for fish= Fish_with_behav_f20
   if ismember(fish,unique(idx_Fish_f20))  
    temp_idx=find(idx_Fish_f20==fish); idx_Mov_temp=[];   
    rsq_temp=[ModelResults_allmovPerFish{1,counter}.rsquared];
    idx_Mov_temp=temp_idx(find(rsq_temp>0.2 & rsq_temp<1),:);
    idx_rsq_Mov{counter}=idx_Mov_temp;
    rsquare_loom_allmov2{counter}=rsq_temp;
    counter=counter+1;
   else
   end    
end

clear idx_Mov_temp temp_idx rsq_temp counter

idx_rsq_Mov=vertcat(idx_rsq_Mov{:});
figure; imagesc(ZS_f20(idx_rsq_Mov,:), [-0.5 4]);colormap hot


rsquare_loom_allmov2=horzcat(rsquare_loom_allmov2{:});

figure;histogram(rsquare_loom_allmov2);

%%%% this doesnt work because of the indexing
% idx_rsq_allmove2=find(rsquare_loom_allmov2>0.15 & rsquare_loom_allmov2<1); %%%then select the rsquare that are between 0.3 and 1
% figure; imagesc(ZS_f20(idx_rsq_allmove2,:), [-0.5 4]);colormap hot


%%%%and this is to plot the means of the traces of the ROIs that had an r2
%%%%value above the threshold set above.
%%%% It makes a figure with a subplot of each fish. 
 
%%% i still need to adjust it

   
    counter=1;
    figure; 
    xplot=floor(sqrt(length(Fish_with_behav_f20)));yplot=ceil(length(Fish_with_behav_f20)/xplot);
    for fish=Fish_with_behav_f20
        if ismember(fish,unique(idx_Fish_f20)) %%isempty(Correlation_movement_post3{fish})~=1 % idx_Fish_f20==Fish_with_behav_f20(fish)
        temp_ZS=ZS_f20(find(idx_Fish_f20==fish),:);        
        rsq_temp=[ModelResults_allmovPerFish{1,counter}.rsquared];
        subplot(xplot,yplot,counter);plot(mean(temp_ZS(find(rsq_temp>0.1 & rsq_temp<1),:),1));
        hold on
%       Movements2=squeeze(Movements(:,i,:)); %%% this was already done before
%       Movements2=squeeze(Movements2);
        plot(Movements2(find(Fish_with_behav_f20==fish),:))
        plot(loomtimes);
        hold off
        
        %mdl=stepwiselm(Stimuli',mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
        %title(num2str([mdl.Rsquared.Adjusted]));
        title(strcat(' nb of ROIs : ', num2str(length(find(rsq_temp>0.1 & rsq_temp<1)))));
        
         counter=counter+1;
        else
        end
    end

    clear temp_ZS rsq_temp counter

    
    %%
   
    
    %%% this is to look for the ROIs of the movements that also correlated
    %%% to looms... it doesnt look good. 
   idx_rsq_Mov_N_Loom_corr=intersect(idx_LoomMov_post3,idx_rsq_Mov);
   figure; imagesc(ZS_f20(idx_rsq_Mov_N_Loom_corr,:), [-0.5 4]);colormap hot
   
   
    %%% this is to look for the ROIs of the movements that also part of the loom responding cells... it doesnt look good. 
   idx_rsq_Mov_N_rsq_test=intersect(idx_rsq_test,idx_rsq_Mov);
   figure; imagesc(ZS_f20(idx_rsq_Mov_N_rsq_test,:), [-0.5 4]);colormap hot
   
   
   temp2=[];
counter=1;
temp2{counter}=idx_rsq_Mov;    
counter=counter+1;    
temp2{counter}=idx_rsq_Mov_N_Loom_corr; 
counter=counter+1;    
temp2{counter}=idx_rsq_Mov_N_rsq_test; 
   

Numbers(1)=0; %%% to change the first value of Numbers to 0 (cause it was not relevant)
colors2 = distinguishable_colors(length(temp2),[1 1 1; 0 0 0]); %%%here we use a script from Matlab (downloaded, and needs to be in the folder) to generate colors
colors2 = colors2*256;
 
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;counter2=1;xplot=floor(sqrt(length(temp2)));yplot=ceil(length(temp2)/xplot); rows=length(temp2);
for i=1:length(temp2)
    idx_temp=temp2{i};
    subplot(rows,4,counter);plot(mean(ZS_f20(idx_temp,:),1),'color',colors(counter2,:)/256); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f20(idx_temp,:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f20(idx_temp)); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f20(idx_temp));%%% for the fish location
    %counter=counter+1;
    counter2=counter2+1
     counter=counter+4;
end


    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%% this part is to make the Kmeans figures with the localization of the selected clusters in the brain  %%%%%%%%%

temp=[];
counter=1;
temp{counter}=idx_LoomMov;    
counter=counter+1;    
temp{counter}=idx_Mov; 
counter=counter+1;    
temp{counter}=idx_good_LoomMov; 
counter=counter+1;    
temp{counter}=idx_LoomMov_post3; 
counter=counter+1;    
temp{counter}=idx_Mov_post3; 
counter=counter+1;    
temp{counter}=idx_good_LoomMov_post3; 
counter=counter+1;    
temp{counter}=idx_filtered_good_LoomMov_post3; 

Numbers(1)=0; %%% to change the first value of Numbers to 0 (cause it was not relevant)
%colors = [1,0,0;0,1,0;0,0,1;1,0.600000000000000,0;1,0,0.7];
colors = distinguishable_colors(length(temp),[1 1 1; 0 0 0]); %%%here we use a script from Matlab (downloaded, and needs to be in the folder) to generate colors
% colors = [0         0    1.0000
%          0    0.5000    1.0000
%     1.0000         0         0
%     1.0000    0.1034    0.7241
%     1.0000    0.5000    0.3000
%          0    0.7000    0.2000
%     0.5000    0.5000         0
%          0    0.5000    0.5000];
colors = colors*256;


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;counter2=1;xplot=floor(sqrt(length(temp)));yplot=ceil(length(temp)/xplot); rows=length(temp);
for i=1:length(temp)
    idx_temp=temp{i};
    subplot(rows,4,counter);plot(mean(ZS_f20(idx_temp,:),1),'color',colors(counter2,:)/256); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f20(idx_temp,:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f20(idx_temp)); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f20(idx_temp));%%% for the fish location
    %counter=counter+1;
    counter2=counter2+1
     counter=counter+4;
end

%saveas(gcf,'multigraph_LoomMov_Mov_CN_F20_R04','jpg');
 %print(Fighandle,strcat(destdirectory1,'multigraph_LoomMov_Mov_post3_CN_F20_R05'),'-dpng','-r0');
    %close all;


idx=1
    filename=MatFiles(idx).name;
%destdirectory = strcat('/QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab','/Kmeans_f20_r2050_CL50_HPC',filename(8:11),'/');
destdirectory = strcat('C:/Emm_temp/FvsS_20vs60_analysis','/Kmeans_f20_CN_rsq_LoomMov_R02_testingthings',filename(8:11),'/');
mkdir(destdirectory);   %create the directory
clearvars idx filename

for idx=1:length(MatFiles)
    filename=MatFiles(idx).name;%%%to put the respective name of the files (in this case the slices)
    ROIsNb=[];ClusterNb=[];%%% make the variables we are going to use
    %for k = 1 : length(temp)
    %%%%%this is to locate in every plane which ROIs and cluster can be
    %%%%%found
    for k = 1 : length(temp2) %%% to select the one i want
        tempROIsNb=find([temp2{k}]<=Numbers(idx+1) & [temp2{k}]>Numbers(idx)); %%%to put the data filtered of the selected clusters on a new variable. but i dont understand the sentece... why not just ask for Numbers(idx)?
        if tempROIsNb            
            ROIsNb=[ROIsNb; temp2{k}(tempROIsNb)];%%%%
            %temp{k}(tempROIsNb)=[];
            ClusterNb=[ClusterNb; repmat(k,length(tempROIsNb),1)];
        end
    end
    
    %%%% and this part is to make and image using the _mean.tif image as
    %%%% base, where we will add the ROIs located before and colored base
    %%%% on the number of the cluster
    if ROIsNb
        imagename=regexp(filename,'_output_analysis','split');
        %imagename=regexp(imagename,'_output_analysis_matlab2.mat','split');
        imagename=strcat(imagename{1},'_mean.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*64;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=MatFiles(idx).ROI;       
        ROIsNb=ROIsNb-Numbers(idx);
        ROIs=ROIs(:,ROIsNb);
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
            for j=1:3
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
        %image3(:,:,3)=image;
        %%%finally we save the Image
        name=strcat(destdirectory,'_Kmeans_f20_CN_rsq_LoomMov_R02_testingthings',imagename(4:end));
    imwrite(image3,name,'tif');
    end
    %image3=uint8(image3);

end


%clearvars idx i tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster
