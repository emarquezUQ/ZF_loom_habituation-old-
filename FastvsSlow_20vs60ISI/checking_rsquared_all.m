%%% this script is just to check at the distribution of the rsquare values
%%% in the 4 datasets to confirm that 0.3 is a good filter. more or less + 2SD

%% for f20

ModelResults_f20=load('final_F20_step1.mat','ModelResults_shortF20_all');



rsq_test_f20_all=[];
for j=1:7
    temp_rsq=[ModelResults_f20.ModelResults_shortF20_all{1,j}.rsquared];
      
rsq_test_f20_all(:,j)=temp_rsq;

end


rsquare_loom_f20=[];
for i=1:length(rsq_test_f20_all)
    [M,I]=max(rsq_test_f20_all(i,:));
    rsquare_loom_f20(i,1)=M;
    
end



figure;histogram(rsquare_loom_f20);
2*std(rsquare_loom_f20) %% the result is 0.2567! 

clearvars -except rsq_test_f20_all rsquare_loom_f20 



%% for f60

load('final_F60_step1_2.mat','ModelResults_shortF60_all');
load('final_F60_step1_2.mat','idx_Fish_f60');

rsq_test_f60_all=[];
for j=1:7
    temp_rsq=[ModelResults_shortF60_all{1,j}.rsquared];
      
rsq_test_f60_all(:,j)=temp_rsq;

end


rsquare_loom_f60=[];
for i=1:length(rsq_test_f60_all)
    [M,I]=max(rsq_test_f60_all(i,:));
    rsquare_loom_f60(i,1)=M;
    
end


figure;histogram(rsquare_loom_f60);
2*std(rsquare_loom_f60) %% the result is 0.2919! so far so goood



clearvars -except rsq_test_f20_all rsquare_loom_f20 rsq_test_f60_all rsquare_loom_f60

save('rsquare_fast.mat');

%% for s20


load('final_S20_step1.mat','ModelResults_shortS20_all');
%load('final_S20_step1_2.mat','idx_Fish_s20');

rsq_test_s20_all=[];
for j=1:7
    temp_rsq=[ModelResults_shortS20_all{1,j}.rsquared];
      
rsq_test_s20_all(:,j)=temp_rsq;

end


rsquare_loom_s20=[];
for i=1:length(rsq_test_s20_all)
    [M,I]=max(rsq_test_s20_all(i,:));
    rsquare_loom_s20(i,1)=M;
    
end


figure;histogram(rsquare_loom_s20);
2*std(rsquare_loom_s20) %% the result is 0.2789! one more to go!


clearvars -except rsq_test_f20_all rsquare_loom_f20 rsq_test_f60_all rsquare_loom_f60 rsq_test_s20_all rsquare_loom_s20

%% for s60

load('final_S60_step1.mat','ModelResults_shortS60_all');
%load('final_S60_step1.mat','idx_Fish_s60');

rsq_test_s60_all=[];
for j=1:7
    temp_rsq=[ModelResults_shortS60_all{1,j}.rsquared];
      
rsq_test_s60_all(:,j)=temp_rsq;

end


rsquare_loom_s60=[];
for i=1:length(rsq_test_s60_all)
    [M,I]=max(rsq_test_s60_all(i,:));
    rsquare_loom_s60(i,1)=M;
    
end


figure;histogram(rsquare_loom_s60);
2*std(rsquare_loom_s60) %% the result is 0.2405 so they all fall under 0.3! excellent!


clearvars -except rsq_test_f20_all rsquare_loom_f20 rsq_test_f60_all rsquare_loom_f60 rsq_test_s20_all rsquare_loom_s20 rsq_test_s60_all rsquare_loom_s60

save('rsquare_slow.mat');
