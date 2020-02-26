%%% i will be trying to get the raw clusters again by puting more clusters
 
 cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab
 
 load('s60_postKmeans_CN')
 
 %%%here i am doing a kmeans to get some raw regressors again but with more
 %%%clusters
[idxKmeans_ZS_CN_200CL Cmap_ZS_CN_200CL]=kmeans(ZS_CN,200,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
%[Model_ZS_CN_200CL,GoodBetas_ZS_CN2_200CL]=Test_Regress(Cmap_ZS_CN_200CL,Cmap_ZS_CN_200CL,idxKmeans_ZS_CN_200CL,0.1);%%%here we do another linear regression and we select the ones with an r2 value above 0.3

saveas(gcf,'200rawclusters_CN_s60_job','jpg');

counter2=0;
for i=1:4
    Fighandle=figure;counter=1;
    set(Fighandle, 'Position', [100, 100, 1500, 1000]);
    for j=1:50
        subplot(7,8,counter);plot(Cmap_ZS_CN_200CL(j+counter2,:));
        counter=counter+1;
    end
    counter2=counter2+50;
    print(Fighandle,strcat('200rawclusters_CN_s60_job_',num2str(i)),'-dpng','-r0');
    close all;
end

 
 save('s60_postKmeans_CN_200CL_job.mat','-v7.3');