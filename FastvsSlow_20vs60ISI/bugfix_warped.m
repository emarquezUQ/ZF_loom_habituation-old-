

list=str2num(Fish_list_f20);
CSV_Files=dir('_2Warped*_Resized.csv');
%CSV_Files=dir('_2Warped1_*.csv'); %% to test with first one
ROIs=struct();truth=[]; counter=1;
for i=1:length(CSV_Files)
      
    Fishname=regexp(CSV_Files(i).name,'_2Warped(\d+)_Resized.csv','tokens');Fishname=Fishname{1}{1};    
    
    for j=1:length(Fish_list_f20)
     
    if str2num(Fishname)==list(j)
    %if ismember(str2num(Fishname),str2num(Fish_list_f20))
    
    temp_warp=csvread(CSV_Files(i).name,1);    
    ROIs(j).name=Fishname;    
    ROIs(j).coord=temp_warp(:,1:3);
    ROIs(j).idx=temp_warp(:,5);
    truth(j)=size(temp_warp,1)==sum(idx_Fish_f20==str2num(Fishname)) %%% important!! i need to have an idx_Fish with all the fish.
    
     counter=counter+1;
    else
    end
    end  
end



