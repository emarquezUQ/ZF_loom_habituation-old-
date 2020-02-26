


%%%% this script is to do the cross-validation test. I will take one fish
%%%% out and regenerate the mean corr-matrices and show that the general
%%%% results still hold. now with the fmr1 dataset

%%% part of the script is based on how gilles did it for the fmr1 dataset

%%
load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

datasets=['f20'; 'f60'; 's20'; 's60'];


keepFmr1=load('graphs_fmr1Exp.mat','keep');
keepFmr1=keepFmr1.keep;



%%
cross_val_fmr1=struct;
for group=1:3

names = fieldnames(Data_corrMat4.(groupnames{group}));
CorrMatrices_mean2=zeros(21,length(names)-1,90,90);
for loom=1:21        
    for fish_rem_nb=1:length(names)-1
        temp=nan(length(names)-2,90,90);    
        counter=1;
        for fish_nb=1:length(names)-1
            if fish_nb ~= fish_rem_nb
                fish_name=names(fish_nb);
                temp(counter,:,:)=Data_corrMat4.(groupnames{group}).(fish_name{1}).loomsR{1,loom}(keepFmr1,keepFmr1);
                counter=counter+1;            
            end
        end
        CorrMatrices_mean2(loom,fish_rem_nb,:,:)=squeeze(nanmean(temp,1));
        
    end       
end
cross_val_fmr1.(groupnames{group}).CorrMatrices_mean2=CorrMatrices_mean2;
end

%% to show it

%%%% 5 random group of matrices without a fish at random
for group=1:3
names = fieldnames(Data_corrMat4.(groupnames{group}));
figure;
counter=1;
for fish=randperm(length(names)-1,5)
for loom=[1 2 3 4 5 6 11 12]

    subplot(5,8,counter);imagesc(squeeze(squeeze(cross_val_fmr1.(groupnames{group}).CorrMatrices_mean2(loom,fish,:,:))));caxis([-1 1])
    
    counter=counter+1;
end
end

end
%%%% changes in density. plot the density for fmr1 experiment and then as scatter plot
%%%% the individual densities of the substracted means

%%%%% getting the thresholed matrices (I need the BCT toolbox)


%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected=struct;
for group=1:3
    datatemp=groupnames{group};
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,4,5,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat4.(groupnames{group}).Mean_corrMat{1,moment(m)}(keepFmr1,keepFmr1)),0.75);
     MatAll_corrected.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end

%% density
%


%kden = density_und(CIJ);

for group=1:3
datatemp=groupnames{group};
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected.(datatemp).(loom{i}).Mat);

MatAll_corrected.(datatemp).(loom{i}).kden=temp_kden;
end
end

figure;
temp=[];
for group=1:3
datatemp=groupnames{group};
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end


for group=1:3

names = fieldnames(Data_corrMat4.(groupnames{group}));

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected_crossval=[];

for fish=1:length(names)-1

    moment=[1 2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    %loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(squeeze(squeeze(cross_val_fmr1.(groupnames{group}).CorrMatrices_mean2(moment(m),fish,:,:)))),0.75);
     MatAll_corrected_crossval(m,fish,:,:)=Mat;
     end

end
cross_val_fmr1.(groupnames{group}).MatAll_corrected_crossval=MatAll_corrected_crossval;

end



%%%% calculating the density

for group=1:3

names = fieldnames(Data_corrMat4.(groupnames{group}));

crossval_density=[];
for fish=1:length(names)-1
for m=1:length(moment)

temp_kden=density_und(squeeze(squeeze(cross_val_fmr1.(groupnames{group}).MatAll_corrected_crossval(m,fish,:,:))));

crossval_density(m,fish)=temp_kden;
end
end
cross_val_fmr1.(groupnames{group}).crossval_density=crossval_density;
end



figure;
temp=[];
for group=1:3
datatemp=groupnames{group};
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

group=2; %%% chose the one to look at
names = fieldnames(Data_corrMat4.(groupnames{group}));
for fish=1:length(names)-1
for m=1:length(moment)-1    
scatter(m,cross_val_fmr1.(groupnames{group}).crossval_density(m+1,fish),'filled');
hold on;
end
end
hold off;


%%% getting the values to do stats. 

density_Groups_fmr1={};

for m=1:length(moment)-1
    temp=nan(20,3);
    counter=1;
for group=[3 1 2] 
names = fieldnames(Data_corrMat4.(groupnames{group}));
for fish=1:length(names)-1
    
temp(fish,counter)=cross_val_fmr1.(groupnames{group}).crossval_density(m+1,fish);

end
counter=counter+1;
end
density_Groups_fmr1{m}=temp;
end


%%

%%% now participation


for group=1:3
datatemp=groupnames{group};
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)


temp_mat=MatAll_corrected.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes.Nod_clustID(keepFmr1));

MatAll_corrected.(datatemp).(loom{i}).P=P;
end
end

figure;
temp=[];
for group=1:3
datatemp=groupnames{group};
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=mean(MatAll_corrected.(datatemp).(loom{i}).P);
end
plot(temp);
hold on;

end

%%%% calculating the participation for the crossval

for group=1:3

names = fieldnames(Data_corrMat4.(groupnames{group}));

crossval_P=[];
for fish=1:length(names)-1
for m=1:length(moment)

    
temp_mat=squeeze(squeeze(cross_val_fmr1.(groupnames{group}).MatAll_corrected_crossval(m,fish,:,:)));
temp_mat(isnan(temp_mat))=0;
    
temp_P=participation_coef(temp_mat,Nodes.Nod_clustID(keepFmr1));

crossval_P(m,fish)=mean(temp_P);
end
end
cross_val_fmr1.(groupnames{group}).crossval_P=crossval_P;
end



figure;
temp=[];
for group=1:3
datatemp=groupnames{group};
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=mean(MatAll_corrected.(datatemp).(loom{i}).P);
end
plot(temp);
hold on;

end

group=2; %%% chose the one to look at
names = fieldnames(Data_corrMat4.(groupnames{group}));
for fish=1:length(names)-1
for m=1:length(moment)-1    
scatter(m,cross_val_fmr1.(groupnames{group}).crossval_P(m+1,fish),'filled');
hold on;
end
end
hold off;


%%% getting the values to do stats. 

Particip_Groups_fmr1={};

for m=1:length(moment)-1
    temp=nan(20,3);
    counter=1;
for group=[3 1 2] 
names = fieldnames(Data_corrMat4.(groupnames{group}));
for fish=1:length(names)-1
    
temp(fish,counter)=cross_val_fmr1.(groupnames{group}).crossval_P(m+1,fish);

end
counter=counter+1;
end
Particip_Groups_fmr1{m}=temp;
end


%% making figures

figure;
counter=1;
for group=[3 1 2]
temp=[];

    for m=1:length(moment)-1
       temp(1,m)=nanmean(density_Groups_fmr1{m}(:,counter));
        
    end  
plot(temp);
legend('WT','Hets','fmr1');
hold on;

counter=counter+1;
end
hold off
%saveas(gcf,'density_crossvall_fmr1.svg');

figure;
counter=1;
for group=[3 1 2]
temp=[];

    for m=1:length(moment)-1
       temp(1,m)=nanmean(Particip_Groups_fmr1{m}(:,counter));
        
    end  
plot(temp);
legend('WT','Hets','fmr1');
hold on;

counter=counter+1;
end

%saveas(gcf,'Part_crossvall_fmr1.svg');


%%
%save('crossvalidation_1fishout_fmr1_results.mat');


