

load('f20_CN_r2050_CL5_extra.mat','MatFiles','Numbers','ZS_CN','goodmaps_CN','idx_Fish','rawregress_CN','idxKmeans_ZS_CN','rsquare_loom');

idx_rsq_test=find(rsquare_loom>0.3 & rsquare_loom<1);


Fish_list_f20=num2str(unique(idx_Fish));


%%%% this is for testing with only f20 fish

%%%% now, this is the part where i used the warped datasets to get the
%%%% brain locations of my ROIs.

CSV_Files=dir('_2Warped*_Resized.csv');
%CSV_Files=dir('_2Warped1_*.csv'); %% to test with first one
ROIs=struct();truth=[]; counter=1;
for i=1:length(CSV_Files)
      
    Fishname=regexp(CSV_Files(i).name,'_2Warped(\d+)_Resized.csv','tokens');Fishname=Fishname{1}{1};    
    
    if ismember(str2num(Fishname),str2num(Fish_list_f20))
    
    temp_warp=csvread(CSV_Files(i).name,1);    
    ROIs(counter).name=Fishname;    
    ROIs(counter).coord=temp_warp(:,1:3);
    ROIs(counter).idx=temp_warp(:,5);
    truth(counter)=size(temp_warp,1)==sum(idx_Fish==str2num(Fishname)) %%% important!! i need to have an idx_Fish with all the fish.
    
     counter=counter+1;
    else
    end
      
end

%%% this is to check for errors. it looks if the number of ROIs in the csv
%%% files is the same than the goodnumber ROIs of each fish. If there are
%%% errors they need to be fixed before to move on. 
for i=1:length(truth)
    if ~truth(i)
        strcat('Fish',num2str(ROIs(i).name),': error = ',num2str(length(ROIs(i).coord)-sum(idx_Fish==str2num(ROIs(i).name))))
    end    
end

%%% this is to get the ROIs of all the fish (or the ones i used the warped files off) and put them all together.
i=1;ROI_pool=ROIs(i).coord;
for i=2:length(ROIs)
    ROI_pool=[ROI_pool; ROIs(i).coord];
end


%%%to get the individual fish nambes or number
idx_Fish_name={};
for i=1:length(MatFiles)
    name=strcat(MatFiles(i).name);
    [Fish,~]=regexp(name,'fish(\d+)_','tokens','match');Fish=Fish{1}; %%%to get the number of the fish 
    if iscell(Fish)
        Fish=Fish{1};
    end
    idx_Fish_name{i}=Fish;
end
clearvars i Fish  name 
Fish_list=unique(idx_Fish_name);


%%% note: %%% i need to get the GoodNumbers again too

Sort_ROIs=[];temp_nb=0;%truth=[];
MatFiles_names={MatFiles.name};
for fish_nb=1:length(Fish_list_f20)
    %temp_warp=num2str(Fish_list(fish_nb)); %Gilles
    
    temp_warp=str2num(Fish_list_f20);
    temp_warp=temp_warp(fish_nb);
	IndexC=strfind({MatFiles.name}, strcat('_fish',num2str(temp_warp),'_')); 
    
    fish_name=strcat('_fish',num2str(temp_warp),'_'); %%% I added this to get the fish name
    
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    for file_nb=1:length(MatFiles_fish)
        if MatFiles_fish(file_nb)==1
            numbersForROIs=[1 MatFiles(1).GoodNumber];
            
        else
            numbersForROIs=[MatFiles(MatFiles_fish(file_nb)-1).GoodNumber+1 MatFiles(MatFiles_fish(file_nb)).GoodNumber]; 
        end
        if ismember(numbersForROIs,Sort_ROIs)
            
            %%% note: %%% i need to get the fish_name again
            fish_name
            break
        end
        Sort_ROIs=[Sort_ROIs numbersForROIs(1):1:numbersForROIs(2)];        
    end    
    if ~length(Sort_ROIs)-temp_nb==sum(idx_Fish==temp_warp)
        %~length(Sort_ROIs)-temp_nb==sum(idx_Fish==str2num(cell2mat(Fish_list_small(fish_nb)))) %%% i added the str2num and cell2mat
        fish_name
        break
    end
    temp_nb=length(Sort_ROIs);
end
clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb

%%% I need to get the Zbrain_Masks

load('Zbrain_Masks.mat');


ROI_fish(Sort_ROIs,:)=ROI_pool; %%% Sort_ROIs and ROI_pool need, to have the same amount of elements .doesnt work... different numbers. ask gilles. I think i know, it is cause i just have the warped images of f20 and s60 but not from the rest. i changed that but still doesnt work...
ROI_temp=round(ROI_fish);           %%% here i am rounding the values of the coordenates cause in the masks they are integers 
ROI_temp = ROI_temp(any(ROI_temp,2),:);%%% i need to adjust cause the sort_ROIs worked using all the Goodnumbers from all the fish.

ROI_temp2=ROI_temp;
ROI_temp2(:,1)=ROI_temp(:,2); 
ROI_temp2(:,2)=ROI_temp(:,1);  


%%% checking that the rotation worked

figure;scatter(ROI_temp(find(idx_Fish==1),1),ROI_temp(find(idx_Fish==1),2)); 

%%% so it seems they are looking up...


figure;scatter(Zbrain_Masks{294,3}(:,1),Zbrain_Masks{294,3}(:,2)); %%% this is telencephalon in Zbrain
%%% and I need them to look to the right. 

figure;scatter(ROI_temp2(find(idx_Fish==1),1),ROI_temp2(find(idx_Fish==1),2)); %%% this is my rotation
%%%as I have it at the moment is looking to the right. 

%%% chekcing in 3D
figure;scatter3(ROI_temp2(find(idx_Fish==1),1),ROI_temp2(find(idx_Fish==1),2),ROI_temp2(find(idx_Fish==1),3)); %%% this is my rotation



%%% now checking on the z axis that my ROIs are not outside the range. 

%%% zbrain OT
figure;histogram(Zbrain_Masks{105,3}(:,3));
%%% so the OT reaches till a bit more than 130. the maximum in zbrain is
%%% 137


max(Zbrain_Masks{105,3}(:,3))


%%% one of my fish
figure;histogram(ROI_temp2(find(idx_Fish==1),3));

%%
%%%% NOTE: I am having super high numbers on the Z axis... but Gilles thinks it should be
%%%% enough to divide it by 2

ROI_fish(Sort_ROIs,:)=ROI_pool; %%% Sort_ROIs and ROI_pool need, to have the same amount of elements .doesnt work... different numbers. ask gilles. I think i know, it is cause i just have the warped images of f20 and s60 but not from the rest. i changed that but still doesnt work...

ROI_fish(:,3)=ROI_fish(:,3)/2;

ROI_temp=round(ROI_fish);           %%% here i am rounding the values of the coordenates cause in the masks they are integers 
ROI_temp = ROI_temp(any(ROI_temp,2),:);%%% i need to adjust cause the sort_ROIs worked using all the Goodnumbers from all the fish.

ROI_temp2=ROI_temp;
ROI_temp2(:,1)=ROI_temp(:,2); 
ROI_temp2(:,2)=ROI_temp(:,1);  

%%

%%% now i can play again with the data. but is not working... 

PerBrainRegions=struct();
RegionList={'Thalamus','Cerebellum','Semicircularis','Telencephalon','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain'};

RegionList2={'Telencephalon','Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain','Msystem','coeruleus'};

progressbar;
for i=1:length(RegionList2)
    progressbar(i/length(RegionList2));
    regionName=RegionList2{i};
    if strcmp(regionName,'Telencephalon')
        Mask=Zbrain_Masks{294,3};
    elseif strcmp(regionName,'Hindbrain')
        Hindbrain_Mask=Zbrain_Masks{259,3};
        Mask=Zbrain_Masks{131,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove cerebellum
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Zbrain_Masks{295,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove MON
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Hindbrain_Mask;
    elseif strcmp(regionName,'Msystem')
        Mask=[];
        Mask_temp=[];
        Msystem_masks=[184 186 187];
        for j=1:3
           %Mask_temp=Zbrain_Masks{Msystem_masks(j),3};
           Mask=vertcat(Mask,Zbrain_Masks{Msystem_masks(j),3});
        end
    
        clear Mask_temp
    else
        Mask=[];
        IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
        IndexC=find(not(cellfun('isempty', IndexC)));
        for j=IndexC
            if isempty(Mask)
                Mask=Zbrain_Masks{j,3};
            else
                Mask=vertcat(Mask,Zbrain_Masks{j,3});
            end
        end
    end
    Mask=unique(Mask,'rows');
    IsInBrainRegion=ismember(ROI_temp2,Mask,'rows');
    PerBrainRegions.(regionName).idx=find(IsInBrainRegion==1);    
end

idxKmeans_final=zeros(size(idxKmeans_ZS_CN));
idxKmeans_final(idx_rsq_test)=idxKmeans_ZS_CN(idx_rsq_test);

idxKmeans2_coef_rsq_CN=idxKmeans_ZS_CN(idx_rsq_test);


%%% note, in many of this figures the indexing is wrong. I need to change
%%% it. 

figure;imagesc(ZS_CN(intersect(PerBrainRegions.Telencephalon.idx,idx_rsq_test),:), [-0.5 4]);colormap hot

figure;imagesc(ZS_CN(intersect(PerBrainRegions.Thalamus.idx,idx_rsq_test),:), [-0.5 4]);colormap hot

figure;imagesc(ZS_CN(intersect(PerBrainRegions.Pretectum.idx,idx_rsq_test),:), [-0.5 4]);colormap hot

figure;imagesc(ZS_CN(intersect(PerBrainRegions.Tectum.idx,idx_rsq_test),:), [-0.5 4]);colormap hot


figure;imagesc(ZS_CN(intersect(PerBrainRegions.Telencephalon.idx,find(idxKmeans_final==47)),:), [-0.5 4]);colormap hot
figure;imagesc(ZS_CN(intersect(PerBrainRegions.Telencephalon.idx,find(idxKmeans_final==45)),:), [-0.5 4]);colormap hot

figure;imagesc(ZS_CN(intersect(PerBrainRegions.Tectum.idx,find(idxKmeans_final==47)),:), [-0.5 4]);colormap hot
figure;imagesc(ZS_CN(intersect(PerBrainRegions.Tectum.idx,find(idxKmeans_final==45)),:), [-0.5 4]);colormap hot


figure;plot(mean(ZS_CN(intersect(PerBrainRegions.Telencephalon.idx,find(idxKmeans_final==47))),:));
figure;plot(mean(ZS_CN(intersect(PerBrainRegions.Telencephalon.idx,find(idxKmeans_final==45))),:));

figure;plot(mean(ZS_CN(intersect(PerBrainRegions.Tectum.idx,find(idxKmeans_final==47)),:)));
figure;plot(mean(ZS_CN(intersect(PerBrainRegions.Tectum.idx,find(idxKmeans_final==45)),:)));


figure;imagesc(ZS_CN(intersect(PerBrainRegions.Telencephalon.idx,find(idxKmeans_final==22))),:), [-0.5 4]);colormap hot
figure;imagesc(ZS_CN(intersect(PerBrainRegions.Tectum.idx,find(idxKmeans_final==22))),:), [-0.5 4]);colormap hot

figure;plot(mean(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Telencephalon.idx,find(idxKmeans2_coef_rsq_CN==22))),:)));
figure;plot(mean(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Tectum.idx,find(idxKmeans2_coef_rsq_CN==22))),:)));



figure;plot(mean(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Thalamus.idx,find(idxKmeans2_coef_rsq_CN==47))),:)));
figure;plot(mean(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Tectum.idx,find(idxKmeans2_coef_rsq_CN==47))),:)));

figure;plot(mean(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Thalamus.idx,find(idxKmeans2_coef_rsq_CN==45))),:)));



%%% to look at sharp fast-hab
figure;
for i=1:length(RegionList2)
    regionName=RegionList2{i};
    plot(mean(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.(regionName).idx,find(idxKmeans2_coef_rsq_CN==47))),:)));
hold on;
end
legend(RegionList2);
hold off;

%%% to look at broad fast-hab
figure;
for i=1:length(RegionList2)
    regionName=RegionList2{i};
    plot(mean(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.(regionName).idx,find(idxKmeans2_coef_rsq_CN==45))),:)));
hold on;


end
legend(RegionList2);
hold off;

%%% all from different regions... i dont know why but it looks too
%%% similar...
figure;
for i=1:length(RegionList2)
    regionName=RegionList2{i};
    plot(mean(ZS_CN(intersect(PerBrainRegions.(regionName).idx,idx_rsq_test),:)));
hold on;


end
legend(RegionList2);
hold off;



for i=1:length(RegionList2)
    regionName=RegionList2{i};
    amountperRegion(:,i)=sum(ismember(PerBrainRegions.(regionName).idx,find(idxKmeans_final==22)));
end
figure;bar(amountperRegion);set(gca,'xticklabel',RegionList2);

%%
RegionList3={'Pallium','Subpallium','Thalamus','Tectum','Hindbrain'};

for i=1:length(RegionList3)
    regionName=RegionList3{i};
    amountperRegion(:,i)=length(intersect(PerBrainRegions.(regionName).idx,find(idxKmeans2_coef_rsq_CN==22)));
end
figure;bar(amountperRegion);set(gca,'xticklabel',RegionList3);



figure;
for i=1:length(RegionList3)
    regionName=RegionList3{i};
    plot(mean(ZS_CN(intersect(PerBrainRegions.(regionName).idx,idx_rsq_test),:)));
hold on;


end
legend(RegionList3);
hold off;


