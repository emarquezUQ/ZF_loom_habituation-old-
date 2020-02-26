
%%%%% this script is to find the possible boundaries of gamma and omega. 

%%%% the maximum gamma could be when having a number of communities that is
%%%% half the number of nodes = 45  or a quarter = 23 and a minimum gamma of 4 (based on our
%%%% functional clusters). 

%%%% the maximum omega could be when there is no change across the time
%%%% dimension in at least one node. the min omega could be when there is a
%%%% change at every time point

%%%% I will be testing this using the results from the consensus attempt

load('S_cons_testing_GnO_measures2.mat');

%%

%%% looking at the concensus
counter=1;
figure;
for i=1:16
    
    for j=1:10
        
   subplot(16,10,counter);imagesc(Big_FMR1_OPT{i,j}.control.S_cons);
   counter=counter+1;
    end
    
end

%%%% looking at a random raw 
counter=1;
figure;
for i=1:16
    
    for j=1:10
      sample=randperm(1000,1);  
   subplot(16,10,counter);imagesc(Big_FMR1_OPT{i,j}.control.S_test(:,:,sample));
   counter=counter+1;
    end
    
end

%%
MAT_O_limits=[];
MAT_G_limits=[];
for i=1:16
    
    for j=1:10
        
        temp_mat=Big_FMR1_OPT{i,j}.control.S_cons;
        nodes=size(temp_mat,1);
        times=size(temp_mat,2);
        
        %%% testing for number of communities
        temp=[];
        for t=1:times
            temp(t)=max(temp_mat(:,t));           
        end
        
        if max(temp)>90  %%% I am testing with a quarter too
            
           MAT_G_limits(i,j)=1;
         elseif min(temp)<4  %%% 4 as is the number of my hab clusters... 
             MAT_G_limits(i,j)=-1;
        else
            MAT_G_limits(i,j)=0;
        end
            
        %%%% testing changes in time (min and max omega)
        temp=[];
        for n=1:nodes
            
            for t=2:length(temp_mat(n,:)) %%%% i am taking from timepoint 2 to compare with previous timepoint
                if temp_mat(n,t)==temp_mat(n,t-1)
                    temp(n,t)=1;
                else
                    temp(n,t)=0;
                end
                               
            end
            
            if length(unique(temp(n,2:end)))<2 %%% from t=2 to start because the first one will have only 0s
                temp_lim(n)=1;
            else
                temp_lim(n)=0;
            end
        end
        
        if sum(temp_lim(2:end))>23 %%%% i am not taking into account the 1 node cause it never changes community as is the reference node for being the first    
        MAT_O_limits(i,j)=1; 
        else
        MAT_O_limits(i,j)=0;    
        end
    
    end
    
end



%%%% to check the differences in flexibility with the
%%%% checking_GnO_measures.m script

good=intersect(find(MAT_G_limits==0),find(MAT_O_limits==0));

goodMat=zeros(size(MAT_G_limits));
goodMat(good)=1;
figure;imagesc(goodMat);

[g,o]=find(goodMat);

%% with sample raw

%%%% looking at a random raw. is probably not as accurate as using the representative partitions but just checking 

%%
MAT_O_limits_raw=[];
MAT_G_limits_raw=[];
sample=randperm(1000,1);
for i=1:16
    
    for j=1:10
        
        temp_mat=Big_FMR1_OPT{i,j}.control.S_test(:,:,sample);
        nodes=size(temp_mat,1);
        times=size(temp_mat,2);
        
        %%% testing for number of communities
        temp=[];
        for t=1:times
            temp(t)=max(temp_mat(:,t));           
        end
        
        if max(temp)>90
            
           MAT_G_limits_raw(i,j)=1;
%         elseif min(temp)<4
%             MAT_G_limits_raw(i,j)=-1;
        else
            MAT_G_limits_raw(i,j)=0;
        end
            
        %%%% testing changes in time (min and max omega)
        temp=[];
        for n=1:nodes
            
            for t=2:length(temp_mat(n,:)) %%%% i am taking from timepoint 2 to compare with previous timepoint
                if temp_mat(n,t)==temp_mat(n,t-1)
                    temp(n,t)=1;
                else
                    temp(n,t)=0;
                end
                               
            end
            
            if length(unique(temp(n,2:end)))<2 %%% from t=2 to start because the first one will have only 0s
                temp_lim(n)=1;
            else
                temp_lim(n)=0;
            end
        end
        
        
        if sum(temp_lim(2:end))>23 %%%% i am not taking into account the 1 node cause it never changes community as is the reference node for being the first    
        MAT_O_limits_raw(i,j)=1; 
        else
        MAT_O_limits_raw(i,j)=0;    
        end
    
    end
    
end

%%%% to check the differences in flexibility with the
%%%% checking_GnO_measures.m script

good=intersect(find(MAT_G_limits_raw==0),find(MAT_O_limits_raw==0));

goodMat=zeros(size(MAT_G_limits_raw));
goodMat(good)=1;
figure;imagesc(goodMat);

[g,o]=find(goodMat);
