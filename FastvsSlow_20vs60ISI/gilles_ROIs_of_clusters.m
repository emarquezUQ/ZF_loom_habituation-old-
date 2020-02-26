CSV_Files=dir('_2WarpedLong*.csv');

ROIs=struct();

for i=1:length(CSV_Files);

    temp=csvread(CSV_Files(i).name,0);   

    Fishname=regexp(CSV_Files(i).name,'_2WarpedLong(\d+).csv','tokens');Fishname=Fishname{1}{1};    

    ROIs(i).name=Fishname;    

    ROIs(i).coord=temp(:,1:3);

    ROIs(i).idx=temp(:,5);

    size(temp,1)==sum(idx_Fish==str2num(Fishname))

end

clearvars i temp CSV_Files Fishname



i=1;ROI_pool=ROIs(i).coord;

for i=2:length(ROIs);

    ROI_pool=[ROI_pool; ROIs(i).coord];

end



Sort_ROIs=[];temp_nb=0;

MatFiles_names={MatFiles.name};

for fish_nb=1:length(Fish_list)

    if fish_nb <3

        fish_name=strcat('f',num2str(Fish_list(fish_nb)),'-');

        IndexC=strfind(MatFiles_names,fish_name);

    else

        fish_name=strcat(num2str(Fish_list(fish_nb)),'_');

        IndexC=regexp(MatFiles_names,strcat('^',fish_name));

    end            

    MatFiles_fish = find(not(cellfun('isempty', IndexC)));

    for file_nb=1:length(MatFiles_fish)

        if MatFiles_fish(file_nb)==1

            numbersForROIs=[1 MatFiles(1).GoodNumber];

        else

            numbersForROIs=[MatFiles(MatFiles_fish(file_nb)-1).GoodNumber+1 MatFiles(MatFiles_fish(file_nb)).GoodNumber];

        end

        if ismember(numbersForROIs,Sort_ROIs)

            fish_name

            break

        end

        Sort_ROIs=[Sort_ROIs numbersForROIs(1):1:numbersForROIs(2)];        

    end    

    if ~length(Sort_ROIs)-temp_nb==sum(idx_Fish==(Fish_list(fish_nb)))

        fish_name

        break

    end

    temp_nb=length(Sort_ROIs);

end

clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb



ROI_fish(Sort_ROIs,:)=ROI_pool;

clearvars ROI_pool

 

%ROTATE ROI_fish IF NEEDED !!!!!!!!!!!

 

 

for i=1:length(GoodBetas_select)

    idx_temp=GoodClusters_goodmembers(i).idx;    %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    CSV_temp=ROI_rotated(idx_temp,:);

    CSV_temp(:,3)=CSV_temp(:,3);

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_clust_Basic',num2str(i),'b.csv');

    csvwrite(filename,CSV_temp);

end
