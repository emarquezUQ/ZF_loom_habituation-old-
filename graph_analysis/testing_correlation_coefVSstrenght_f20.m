
%%%% this scrit is to check the correlation between the correlation coef
%%%% and the intensity of the response in the nodes. this is to try to show
%%%% to reviewer one that we are not including noise
%%%% i will look at the max response at each timepoint for each fish. first
%%%% in the f20 dataset. 


%%%% Note: we finally discarded this at it doesnt seem necessary. 

%%
load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

datasets=['f20'; 'f60'; 's20'; 's60'];


%%% I need to do it per loom. as previously
%%%from the wildtype loomhab dataset
load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx','S_trim');

%%

Data_corrMat_Max=struct;

for data=1:4
    datatemp=datasets(data,:);
    
    fish=fieldnames(Nodes.(datatemp).NodeMats);
    
    for tempfish=1:length(fish)
        temp_mean=Nodes.(datatemp).NodeMats.(fish{tempfish});
        
        temp_R={};
        temp_Max_mat={};
        for k=1:31
        
        if k==1
            temp_R{1,k}=max((temp_mean(:,10:45))');
        elseif k==11||k==21 ||k==31
            temp_R{1,k}= max((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28))');
        else           
            temp_R{1,k}= max((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36))');                           
        end
        
        
        for i=1:length(temp_R{1,k})
        for j=1:length(temp_R{1,k})
            temp_Max_mat{1,k}(i,j)=temp_R{1,k}(j);
            
        end
        end
        
        end
        
    Data_corrMat_Max.(datatemp).(fish{tempfish}).Max=temp_R;
    
    Data_corrMat_Max.(datatemp).(fish{tempfish}).loomsR=temp_Max_mat;
   
    
    end            
end

%%

%% ploting it

%%%% it looks terrible... hard to show the point. but we thought that is
%%%% not needed anymore

figure;
for data=1%:4
    datatemp=datasets(data,:);
    
    fish=fieldnames(Nodes.(datatemp).NodeMats);
    
    for tempfish=1:length(fish)
              
    for loom=1:31
        
    temp1=abs(Data_corrMat2.(datatemp).(fish{tempfish}).loomsR{1,loom});   
    temp2=Data_corrMat_Max.(datatemp).(fish{tempfish}).loomsR{1,loom};

    [~,~,temp1] = find(tril(temp1,-1));
    [~,~,temp2] = find(tril(temp2,-1));
    
    scatter(temp2,temp1,'.');
    
    hold on;
    end
    end            
end


%% gilles plots, per loom
figure;
for data=1%:4
    datatemp=datasets(data,:);
    
    fish=fieldnames(Nodes.(datatemp).NodeMats);
    
    for tempfish=1:length(fish)
        
        for loom=1:12
            
            temp1=abs(Data_corrMat2.(datatemp).(fish{tempfish}).loomsR{1,loom});
            temp2=Data_corrMat_Max.(datatemp).(fish{tempfish}).loomsR{1,loom};
            
            [~,~,temp1] = find(tril(temp1,-1));
            [~,~,temp2] = find(tril(temp2,-1));
            subplot(4,3,loom);
            scatter(temp2,temp1,'.');axis([0 15 0 1])            
            
            hold on;
        end
    end
end


