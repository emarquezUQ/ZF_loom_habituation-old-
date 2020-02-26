

%%% this script to compare the means of the 4 datasets. FvsS and 20vs60

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




counter=1;
All_means_CL7_normalized={};
figure;
for i=1:4
    datasetname=strcat('mean_CL7_',datasets(i,:));
    tempVar=eval(datasetname);
    
    if ~isempty(regexp(datasetname,'f'))
        for j=1:7
      
            tempClust=tempVar.(clusters{j,1});
         
             All_means_CL7_normalized.(datasets(i,:)).(clusters{j,1})=(tempClust(1,1:1344) - tempClust(1,60)) / ( max(tempClust(1,1:1344)) - tempClust(1,60) );    

             subplot(4,7,counter);plot(All_means_CL7_normalized.(datasets(i,:)).(clusters{j,1}));

             counter=counter+1
        end
    else
        for j=1:7     

            tempClust=tempVar.(clusters{j,1});

             All_means_CL7_normalized.(datasets(i,:)).(clusters{j,1})=(tempClust(1,S_trim) - tempClust(1,60)) / ( max(tempClust(1,S_trim)) - tempClust(1,60) );   

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

%%% now to normalize it base on the max of the first loom and the onset of
%%% every loom

%%% it doesnt seem to change much and there some weird glitches... 
%%% also, it seems that i need to do it in 2 steps. first normalize it all
%%% to the first loom and then substract the drift. 

%%% firs i need the onsets and the periods of each loom to normalize


%%% to take the loom onsets based on the movie

%%% for f20
Loomf20_onset=zeros(6720,1);
startpoints=[30,254,478]; %%in seconds
loom_times=[0,22,40,60,78,100,120,140,158,180]; %%in seconds
for p=1:3
for k=1:10
    Loomf20_onset(startpoints(p)*10+loom_times(k)*10)=2;
end
end

figure; plot(Loomf20_onset);

%%

Loomf20_onset_idx=find(Loomf20_onset==2);
Loomf20_onset_idx=round(Loomf20_onset_idx/5);
Loomf20_onset_idx(Loomf20_onset_idx<1)=[];
Loomf20_onset_idx(Loomf20_onset_idx>(length(allfish_f20looms_tailmov)-1))=[];        
  
%%
%%% this is to try to get the idx to normalize with each loom. 


loom_moments={};
for i=1:30
    if i==1||i==11||i==21
    loom_moments{i}=(Loomf20_onset_idx(i)-59):Loomf20_onset_idx(i+1);
    elseif i==10||i==20||i==30
      loom_moments{i}=Loomf20_onset_idx(i)+1:Loomf20_onset_idx(i)+28;  
    else
     loom_moments{i}=Loomf20_onset_idx(i)+1:Loomf20_onset_idx(i+1);   
    end
end

%%% to check
for i=1:30
temp(1,i)=length(loom_moments{i});

end
sum(temp(1,:))



%%% now to normalize it base on the max of the first loom and the onset of
%%% every loom

%%% it doesnt seem to change much and there some weird glitches... 
%%% also, it seems that i need to do it in 2 steps. first normalize it all
%%% to the first loom and then substract the drift. 

counter=1;
All_means_CL7_normalized2={};
figure;
for i=1:4
    datasetname=strcat('mean_CL7_',datasets(i,:));
    tempVar=eval(datasetname);
    
    if ~isempty(regexp(datasetname,'f'))
        for j=1:7
         
            tempClust=tempVar.(clusters{j,1});
            
         tempClust =(tempClust(1,1:1344) - tempClust(1,60)) / ( max(tempClust(1,1:1344)) - tempClust(1,60) );    
  
%             for k=1:30
%               tempClust(1,loom_moments{1,k})=(tempClust(1,loom_moments{1,k}) - tempClust(1,Loomf20_onset_idx(k)));    
%               
%               
%             end
            



             All_means_CL7_normalized2.(datasets(i,:)).(clusters{j,1})=  tempClust;  

             subplot(4,7,counter);plot(All_means_CL7_normalized2.(datasets(i,:)).(clusters{j,1}));

             counter=counter+1
        end
    else
        for j=1:7     

            tempClust=tempVar.(clusters{j,1});
            tempClust=tempClust(1,S_trim)
            
            tempClust=(tempClust(1,1:1344) - tempClust(1,60)) / ( max(tempClust(1,1:1344)) - tempClust(1,60) );   

%              for k=1:30
%               tempClust(1,loom_moments{1,k})=(tempClust(1,loom_moments{1,k}) - tempClust(1,Loomf20_onset_idx(k)));    
%                 
%               
%              end
            




              All_means_CL7_normalized2.(datasets(i,:)).(clusters{j,1})=  tempClust;
              
             subplot(4,7,counter);plot(All_means_CL7_normalized2.(datasets(i,:)).(clusters{j,1}));
        
             counter=counter+1

        end
    end  

end

clear tempVar tempClust datasetname



 
All_means_CL7_normalized2=struct2cell(All_means_CL7_normalized2);


        for j=1:7
figure;
      for i=1:4
           plot(All_means_CL7_normalized2{i,1}.(clusters{j,1}))
            
hold on;
      end 

    

        end





        
        
        
%%
%%% I will try to find the peaks to look at the mean and SD to be able to
%%% compare them. 
%%% NOTE: If I do this I need to get the means per fish!!

%%% this is to try a method to get only the peaks of the loom responses and
%%% the their difference to the onset of each loom all this compared with
%%% the first loom response. 

%%% i will do it by playing around with the means dataset as dommy data. it
%%% should also hint at possible differences. 


   all_matx={};
        for j=1:7
      for i=1:4        
       Matx(i,:)= All_means_CL7_normalized{i,1}.(clusters{j,1}); 
      
      end 
      
       all_matx{j}=Matx;
        end

Max_resp={};
for j=1:7
    temp_Max_resp=[]; 
for i=1:4
for k=1:30

   temp_Max_resp(i,k)= max(all_matx{j}(i,loom_moments{1,k}))-min(all_matx{j}(i,Loomf20_onset_idx(k)+1));
end


temp_Max_resp(i,:)=temp_Max_resp(i,:)/temp_Max_resp(i,1);

end
Max_resp{j}=temp_Max_resp;
end


for j=1:7
figure;
for i=1:4

 plot(Max_resp{j}(i,:));    
 %plot(Max_resp{j}(i,:),'.'); 
%bar(Max_resp(i,:));
hold on;
end
end

 
