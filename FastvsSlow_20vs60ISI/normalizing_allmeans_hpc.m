

%%% this script to compare the means of the 4 datasets. FvsS and 20vs60
cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab


load('means_F20_CL4n7.mat')

load('means_F60_CL4n7.mat')

load('means_S20_CL4n7.mat')

load('means_S60_CL4n7.mat')



figure;

plot(mean_CL7_s20.nonhab(1,1:1344));

hold on;

plot(mean_CL7_f20.nonhab(1,1:1344));



S_trim=[1:448 453:901 906:1352];



figure;

plot(mean_CL7_s20.nonhab(1,S_trim));

hold on;

plot(mean_CL7_f20.nonhab(1,1:1344));

hold on;

plot(mean_CL7_s60.nonhab(1,S_trim));

hold on;

plot(mean_CL7_f60.nonhab(1,1:1344));



 

norm_s20_nonhab=(mean_CL7_s20.nonhab(1,S_trim) - min(mean_CL7_s20.nonhab(1,S_trim))) / ( max(mean_CL7_s20.nonhab(1,S_trim)) - min(mean_CL7_s20.nonhab(1,S_trim)) );

norm_f20_nonhab=(mean_CL7_f20.nonhab(1,1:1344) - min(mean_CL7_f20.nonhab(1,1:1344))) / ( max(mean_CL7_f20.nonhab(1,1:1344)) - min(mean_CL7_f20.nonhab(1,1:1344)) );

norm_s60_nonhab=(mean_CL7_s60.nonhab(1,S_trim) - min(mean_CL7_s60.nonhab(1,S_trim))) / ( max(mean_CL7_s60.nonhab(1,S_trim)) - min(mean_CL7_s60.nonhab(1,S_trim)) );

norm_f60_nonhab=(mean_CL7_f60.nonhab(1,1:1344) - min(mean_CL7_f60.nonhab(1,1:1344))) / ( max(mean_CL7_f60.nonhab(1,1:1344)) - min(mean_CL7_f60.nonhab(1,1:1344)) );



 

 

figure;

plot(norm_s20_nonhab);

hold on;

plot(norm_f20_nonhab);

hold on;

plot(norm_s60_nonhab);

hold on;

plot(norm_f60_nonhab);



%%%% it works! now I need to do if for all the clusters... I will try with

%%%% a loop



 

datasets=['f20'; 'f60'; 's20'; 's60'];

clusters=fieldnames(mean_CL7_f20);


counter=1;
All_means_CL7_normalized={};

figure;

for i=1:4

    datasetname=strcat('mean_CL7_',datasets(i,:));

    tempVar=eval(datasetname);

    

    if ~isempty(regexp(datasetname,'f'))

        for j=1:7

      
            tempClust=tempVar.(clusters{j,1});

            

             All_means_CL7_normalized.(datasets(i,:)).(clusters{j,1})=(tempClust(1,1:1344) - min(tempClust(1,1:1344))) / ( max(tempClust(1,1:1344)) - min(tempClust(1,1:1344)) );



        

             subplot(4,7,counter);plot(All_means_CL7_normalized.(datasets(i,:)).(clusters{j,1}));

             counter=counter+1

        end

    

    else

        for j=1:7

       

            tempClust=tempVar.(clusters{j,1});

            

             All_means_CL7_normalized.(datasets(i,:)).(clusters{j,1})=(tempClust(1,S_trim) - min(tempClust(1,S_trim))) / ( max(tempClust(1,S_trim)) - min(tempClust(1,S_trim)) );



        

             subplot(4,7,counter);plot(All_means_CL7_normalized.(datasets(i,:)).(clusters{j,1}));

        counter=counter+1

        end

    end

    

end



clear tempVar tempClust datasetname



 
All_means_CL7_normalized=struct2cell(All_means_CL7_normalized);


        for j=1:7
figure;
      for i=1:4
           plot(All_means_CL7_normalized{i,1}.(clusters{j,1}))
            
hold on;
      end 

    

        end

    

%%
%%% what if I dont normalize it...

%%% it gets biased cause the 60 movies are much longer and have periods
%%% where nothing is happening. so the variance changes and makes the first
%%% spike to be stronger... 


%%
%%% I will try to find the peaks to look at the mean and SD to be able to
%%% compare them. 
%%% NOTE: If I do this I need to get the means per fish!!

    
