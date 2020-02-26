
%%%% this script is for checking the participation at different correlation thresholds


figure;
counter=1;
for threshold=[0.05:0.05:0.95]

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected=struct;
for data=1:4
    datatemp=datasets(data,:);
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat2.(datatemp).Mean_corrMat{1,moment(m)}(keep,keep)),threshold);
     MatAll_corrected.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end

%% participation
%

for data=1%:4
    datatemp=datasets(data,:);
             
        loom=fieldnames(MatAll_corrected.(datatemp));
        
        for i=1:length(loom)
            
            temp_mat=MatAll_corrected.(datatemp).(loom{i}).Mat;
            temp_mat(isnan(temp_mat))=0;
            P=participation_coef(temp_mat,ClustID);
            Ccoef=clustering_coef_wu(temp_mat);
            
            MatAll_corrected.(datatemp).(loom{i}).Ccoef=Ccoef;
            MatAll_corrected.(datatemp).(loom{i}).P=P;
            
        end
end
    
temp=[];
for data=1%:4
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected.(datatemp));
    
    for i=1:length(loom)
        
        temp(:,i)=MatAll_corrected.(datatemp).(loom{i}).Ccoef;
    end
    
    temp_mean(counter,:)=nanmean(temp);
    
    plot(nanmean(temp));
    hold on;
        
end
counter=counter+1;

end
plot(nanmean(temp_mean),'k','LineWidth',3);



