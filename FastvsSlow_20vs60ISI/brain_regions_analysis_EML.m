%%% this code is based on the one Gilles used to analyse Rebecca_aud data
%%% but i need to adjust it

%%% is to look at the regions of the brain in the zbrain atlas and the ROIs
%%% from our experiments

%%
%%% checking that the new generated ROIs for ANTs have only the relevent ROIS

%%% it seems to be working well except for fish 9... (where i seem to loose 1 ROI. not sure why). apart from that i have an error of -1
%%% in all the fish but thats because gilles was taking away the first
%%% lineonf the warped files (for which this code was done in the first
%%% place). 

%%% I will need to first get a list of idx_Fish with all ROIs and all fish
%%% and then I will check if the new version of makeing ROIs for ants is
%%% working.

%%% i need to run this before to get the idx_Fish of all the fish.


%%%the following is to create 1 column variables with the number of the
%%%fish and the slice number
Numbers_all=[1 [MatFiles.GoodNumber]]; %%%to take make a vector with the GoodNumber field
counter=1;
idx_Plane_all=nan(length(length(Fitness)),1);%%%% to make an empty (with nans) one column variable the size of ZS
idx_Fish_all=nan(length(length(Fitness)),1);%%%% to make an empty (with nans) one column variable the size of ZS
name=strcat(MatFiles(1).name);%%%to get the name of the files (is actually to create the variable name before the loop)
for i=1:length(MatFiles) %%%%to take slices one by one	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'Slice(\d+)_','tokens','match');Plane=str2num(Plane{1}{1}); %%%to get the number of the plane    
    idx_Plane_all(Numbers_all(i):Numbers_all(i+1))=Plane; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Plane   
    [Fish,~]=regexp(name,'fish(\d+)_','tokens','match');Fish=str2num(Fish{1}{1}); %%%to get the number of the fish 
    idx_Fish_all(Numbers_all(i):Numbers_all(i+1))=Fish; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Fish
end
clearvars i Fish Plane name counter %%%to get rid of vairables we will not use anymore




%%% and this part was to check for the new ROIs for ANTs

CSV_Files2=dir('_ROIsFish*.csv');
%CSV_Files2=dir('_ROIsFish1.csv'); %% to test with first one
%M = csvread('_ROIsFish1.csv');
ROIs_new=struct();truth_2=[];
for i=1:length(CSV_Files2)
    temp_fishROI=csvread(CSV_Files2(i).name,1);   
    Fishname=regexp(CSV_Files2(i).name,'_ROIsFish(\d+).csv','tokens');Fishname=Fishname{1}{1};    
    ROIs_new(i).name=Fishname;    
    %ROIs_new(i).coord=temp_fishROI(:,1:3);
    ROIs_new(i).idx=temp_fishROI(:,5);
    truth_2(i)=size(temp_fishROI,1)==sum(idx_Fish_all==str2num(Fishname))
end
%clearvars i temp CSV_Files Fishname

%%% this is to check for errors. it looks if the number of ROIs in the csv
%%% files is the same than the goodnumber ROIs of each fish. If there are
%%% errors they need to be fixed before to move on. 
for i=1:length(truth_2)
    if ~truth_2(i)
        strcat('Fish',num2str(ROIs_new(i).name),': error = ',num2str(length(ROIs_new(i).idx)-sum(idx_Fish_all==str2num(ROIs_new(i).name))))
    end    
end

clearvars i temp_fishROI CSV_Files2 Fishname ROIs_new



%%

%%%% now, this is the part where i used the warped datasets to get the
%%%% brain locations of my ROIs. 




CSV_Files=dir('_2Warped*_Resized.csv');
%CSV_Files=dir('_2Warped1_*.csv'); %% to test with first one
ROIs=struct();truth=[];
for i=1:length(CSV_Files)
    temp_warp=csvread(CSV_Files(i).name,1);   
    Fishname=regexp(CSV_Files(i).name,'_2Warped(\d+)_Resized.csv','tokens');Fishname=Fishname{1}{1};    
    
    ROIs(i).name=Fishname;    
    ROIs(i).coord=temp_warp(:,1:3);
    ROIs(i).idx=temp_warp(:,5);
    truth(i)=size(temp_warp,1)==sum(idx_Fish_all==str2num(Fishname)) %%% important!! i need to have an idx_Fish with all the fish.
end
clearvars i temp CSV_Files Fishname

%%% this is to check for errors. it looks if the number of ROIs in the csv
%%% files is the same than the goodnumber ROIs of each fish. If there are
%%% errors they need to be fixed before to move on. 
for i=1:length(truth)
    if ~truth(i)
        strcat('Fish',num2str(ROIs(i).name),': error = ',num2str(length(ROIs(i).coord)-sum(idx_Fish_all==str2num(ROIs(i).name))))
    end    
end


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


%%% this is to get the ROIs of all the fish (or the ones i used the warped files off) and put them all together.
i=1;ROI_pool=ROIs(i).coord;
for i=2:length(ROIs)
    ROI_pool=[ROI_pool; ROIs(i).coord];
end



%%% this is not relvant for me... or maybe some of it?

i=1;Fish=ROIs(i).name;Fish=str2num(Fish);
%ROI_genotype=[];
 ROI_FishNb=[];
 for i=1:length(ROIs)
     temp=zeros(1,length(ROIs(i).coord));
     Fish=ROIs(i).name;Fish=str2num(Fish);
     ROI_FishNb=[ROI_FishNb temp+Fish];
%     if ismember(Fish,WT_fishNB)
%         temp=temp+1;
%         ROI_genotype=[ROI_genotype temp];
%     elseif ismember(Fish,FMR_fishNB)
%         temp=temp+2;
%         ROI_genotype=[ROI_genotype temp];
%     else
%         ROI_genotype=[ROI_genotype temp];
%     end
end


%%% to get the data
MatFiles=dir('*analysis_matlab.mat');


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

%%
MatFiles=dir('*analysis_matlab.mat'); %%to get the files
name=strcat(MatFiles(1).name); %%%to get the name of the files
%Calcium=load(name, 'DenoisedTraces'); %%to load only the DenoisedTraces from the file, the raw data was denoised by the CNMF (The Cluster Analysis tool calculates clusters based on a Constrained non-negative matrix factorization (NMF) clustering method.)
%Calcium=Calcium.DenoisedTraces; %%%% <-- take the field called DenoisedTraces from the Calcium structure and make it the new Calcium
Fitness=load(name, 'idx_components');%%to load only the idx_components from the file, they are based on what a Gcamp spike should be and they will filter the true spikes in our data
Fitness=Fitness.idx_components+1; %%%% <-- take the field called idx_components from the Fitness structure and make it the new Fitness but why +1?? Because python indexing starts at 0 ant matlab at 1
%GoodCalcium=Calcium(Fitness,:); 


MatFiles(1).GoodNumber=length(Fitness); %%%% <-- Create a field inside MatFilesCalcium called GoodNumber the size of Fitness.
for i = 2:length(MatFiles) %%%%to take the slices one by one starting by the second one cause we already did this with the first one
    %%%% we are going to do the same thing that before but for all the
    %%%% slices
name=strcat(MatFiles(i).name);%%%%to take the name of the slice in turn
%C=load(name, 'DenoisedTraces');%%to load only the DenoisedTraces from the file
%C=C.DenoisedTraces;%%%% <-- take the field called DenoisedTraces from the C structure and make it the new C
%     if i==3
%         C=[C(:,1) C(:,1) C(:,1:58)];
%     end
% S=load(name, 'Spikes');
% S=S.Spikes;
%N=load(name, 'Noise');
%N=N.Noise;
F=load(name, 'idx_components');
F=F.idx_components+1;%%%because indexing in pythong is from 0 and matlab is at 1
%D=load(name, 'dFonF');
%D=D.dFonF;
%GC=C(F,:);
%GS=S(F,:);
%GD=D(F,:);

    %Noise=vertcat(Noise,N);
    %GN=N(F,:);
    %Calcium=vertcat(Calcium,C);
    %DF=vertcat(DF,D);
    %Spikes=vertcat(Spikes,S);
    Fitness=horzcat(Fitness,F);
    %GoodCalcium=vertcat(GoodCalcium,GC); %The fish 20+ are longer
    %GoodNoise=vertcat(GoodNoise,GN);
    %GoodDF=vertcat(GoodDF,GD);
    %GoodSpikes=vertcat(GoodSpikes,GS);

%MatFiles(i).number=size(Calcium,1);
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N; %%%to get rid of vairables we will not use anymore




%%



%%% in this part is important that i can get the files according to the
%%% name they have... i might need to play with some things to make it work.


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

ROI_fish(Sort_ROIs,:)=ROI_pool; %%% Sort_ROIs and ROI_pool need, to have the same amount of elements .doesnt work... different numbers. ask gilles. I think i know, it is cause i just have the warped images of f20 and s60 but not from the rest. i changed that but still doesnt work...
ROI_temp=round(ROI_fish);           %%% here i am rounding the values of the coordenates cause in the masks they are integers 
ROI_temp(:,1)=round(ROI_fish(:,2)); %%% this is to make the x axis the new y axis
ROI_temp(:,2)=round(ROI_fish(:,1));  %%% this is to make the y axis the new x axis
%ROI_genotype(Sort_ROIs)=ROI_genotype;
ROI_FishNb(Sort_ROIs)=ROI_FishNb;%%%???? what is happpening here? this variable doesnt appear before


%%% I need to get the Zbrain_Masks

load('Zbrain_Masks.mat');

%%% I also need to rotate my ROIs cause my fish are pointing down but the
%%% Zbrain mask is pointing to the right (so a 90deg rotation). I will try
%%% first by just inverting the x and y columns in the ROIs. It seems that
%%% that was done before already in Rebeccas code.

%%% i also need to round them cause the mask values are in integers!!!


PerBrainRegions=struct();
RegionList={'Thalamus','Cerebellum','Semicircularis','Telencephalon','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain'};
progressbar;
for i=1:length(RegionList)
    progressbar(i/length(RegionList));
    regionName=RegionList{i};
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
    IsInBrainRegion=ismember(ROI_temp,Mask,'rows');
    PerBrainRegions.(regionName).idx=find(IsInBrainRegion==1);    
end

%%
%%% and this is to look at stuff. I problably dont need all of it
 
% for i=1:length(RegionList)
%     regionName=RegionList{i};
%     %PerBrainRegions.(regionName).genotype=ROI_genotype(PerBrainRegions.(regionName).idx);
%     PerBrainRegions.(regionName).FishNb=ROI_FishNb(PerBrainRegions.(regionName).idx);
%     %PerBrainRegions.(regionName).ZS_WT=ZS2(PerBrainRegions.(regionName).idx(PerBrainRegions.(regionName).genotype==1),:);
%     %PerBrainRegions.(regionName).ZS_FMR=ZS2(PerBrainRegions.(regionName).idx(PerBrainRegions.(regionName).genotype==2),:);
% end



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


%%% now checking on the z axis that my ROIs are not outside the range. 

%%% zbrain OT
figure;histogram(Zbrain_Masks{105,3}(:,3));
%%% so the OT reaches till a bit more than 130. the maximum in zbrain is
%%% 137


max(Zbrain_Masks{105,3}(:,3))


%%% one of my fish
figure;histogram(ROI_temp2(find(idx_Fish==1),3));
%%% but I am getting much higher values... from -9 to 185.
%%% all my fish go from -22 to 190. but it seems to vary from fish to fish.
%%% so i will need to adjust it to each fish. 

%%% I will test a formula for rescaling. 

min(ROI_temp2(:,3))
max(ROI_temp2(find(idx_Fish==1),3))
min(ROI_temp2(find(idx_Fish==1),3))

old_Zfish1=ROI_temp2(find(idx_Fish==1),3);
new_Zfish1=(((137-0)*(old_Zfish1-min(old_Zfish1)))/(max(old_Zfish1)-min(old_Zfish1)))+0; 

%%% or:

new_Zfish1=(old_Zfish1-min(old_Zfish1))/(max(old_Zfish1)-min(old_Zfish1))*(137-0)+0; 

figure;histogram(new_Zfish1);

min(new_Zfish1)
max(new_Zfish1)

%%% it worked!! so now a loop for each fish

%%% NOTE: my fish dont go as deep as the zbrain stack... so i will try to
%%% fit it from 45-135 instead. it seems that it didnt work... maybe cause
%%% is already taken into account after the warping. I will try again with
%%% 0-137


ROI_fish2=ROI_fish;
ROI_fish2 = ROI_fish2(any(ROI_fish2,2),:);%%% i need to adjust cause the sort_ROIs worked using all the Goodnumbers from all the fish.


for i=1:length(Fish_list_f20)
    
    temp_fish=str2num(Fish_list_f20);
    temp_fish=temp_fish(i);
    temp_Zvalues=ROI_fish(find(idx_Fish==temp_fish),3);
    new_Zvalues=(temp_Zvalues-min(temp_Zvalues))/(max(temp_Zvalues)-min(temp_Zvalues))*(137-0)+0;
    ROI_fish2(find(idx_Fish==temp_fish),3)=new_Zvalues;
end


clear temp_fish temp_Zvalues new_Zvalues


ROI_temp2=round(ROI_fish2);
ROI_temp2(:,1)=ROI_fish2(:,2); 
ROI_temp2(:,2)=ROI_fish2(:,1); 

%%% checking rotation and Z values


figure;scatter(ROI_temp2(find(idx_Fish==1),1),ROI_temp2(find(idx_Fish==1),2)); %%% this is my rotation
%%%as I have it at the moment is looking to the right. 

figure;histogram(ROI_temp2(:,3));

min(ROI_temp2(:,3))
max(ROI_temp2(:,3))

%%

%%% now i can play again with the data. but is not working... why?? stuck
%%% here. 

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

idxKmeans2_coef_rsq_CN=idxKmeans_ZS_CN(idx_rsq_test);

figure;imagesc(ZS_CN(intersect(PerBrainRegions.Telencephalon.idx,idx_rsq_test),:), [-0.5 4]);colormap hot

figure;imagesc(ZS_CN(intersect(PerBrainRegions.Thalamus.idx,idx_rsq_test),:), [-0.5 4]);colormap hot

figure;imagesc(ZS_CN(intersect(PerBrainRegions.Pretectum.idx,idx_rsq_test),:), [-0.5 4]);colormap hot

figure;imagesc(ZS_CN(intersect(PerBrainRegions.Tectum.idx,idx_rsq_test),:), [-0.5 4]);colormap hot


figure;imagesc(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Telencephalon.idx,find(idxKmeans2_coef_rsq_CN==47))),:), [-0.5 4]);colormap hot
figure;imagesc(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Telencephalon.idx,find(idxKmeans2_coef_rsq_CN==45))),:), [-0.5 4]);colormap hot

figure;imagesc(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Tectum.idx,find(idxKmeans2_coef_rsq_CN==47))),:), [-0.5 4]);colormap hot
figure;imagesc(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Tectum.idx,find(idxKmeans2_coef_rsq_CN==45))),:), [-0.5 4]);colormap hot



figure;plot(mean(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Telencephalon.idx,find(idxKmeans2_coef_rsq_CN==47))),:)));
figure;plot(mean(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Telencephalon.idx,find(idxKmeans2_coef_rsq_CN==45))),:)));

figure;plot(mean(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Tectum.idx,find(idxKmeans2_coef_rsq_CN==47))),:)));
figure;plot(mean(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Tectum.idx,find(idxKmeans2_coef_rsq_CN==45))),:)));


figure;imagesc(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Telencephalon.idx,find(idxKmeans2_coef_rsq_CN==22))),:), [-0.5 4]);colormap hot
figure;imagesc(ZS_CN(idx_rsq_test(intersect(PerBrainRegions.Tectum.idx,find(idxKmeans2_coef_rsq_CN==22))),:), [-0.5 4]);colormap hot

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
    amountperRegion(:,i)=length(intersect(PerBrainRegions.(regionName).idx,find(idxKmeans2_coef_rsq_CN==22)));
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




