
%%%% this script is to do the cross-validation test. I will take one fish
%%%% out and regenerate the mean corr-matrices and show that the general
%%%% results still hold.

%%% part of the script is based on how gilles did it for the fmr1 dataset

%%
load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

datasets=['f20'; 'f60'; 's20'; 's60'];


%%
for data=1%:4

names = fieldnames(Data_corrMat2.(datasets(data,:)));
CorrMatrices_mean2=zeros(31,length(names)-1,99,99);
for loom=1:31        
    for fish_rem_nb=1:length(names)-1
        temp=nan(length(names)-2,99,99);    
        counter=1;
        for fish_nb=1:length(names)-1
            if fish_nb ~= fish_rem_nb
                fish_name=names(fish_nb);
                temp(counter,:,:)=Data_corrMat2.(datasets(data,:)).(fish_name{1}).loomsR{1,loom}(keep,keep);
                counter=counter+1;            
            end
        end
        CorrMatrices_mean2(loom,fish_rem_nb,:,:)=squeeze(nanmean(temp,1));
    end       
end

end

%% to show it

%%%% 5 random group of matrices without a fish at random

figure;
counter=1;
for fish=randperm(length(names)-1,5)
for loom=[1 2 3 4 5 6 11 12]

    subplot(5,8,counter);imagesc(squeeze(squeeze(CorrMatrices_mean2(loom,fish,:,:))));caxis([-1 1])
    
    counter=counter+1;
end
end


%%%% changes in density. plot the density for f20 and then as scatter plot
%%%% the individual densities of the substracted means

%%%%% getting the thresholed matrices (I need the BCT toolbox)


%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected=struct;
for data=1:4
    datatemp=datasets(data,:);
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat2.(datatemp).Mean_corrMat{1,moment(m)}(keep,keep)),0.75);
     MatAll_corrected.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end

%% density
%


%kden = density_und(CIJ);

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected.(datatemp).(loom{i}).Mat);

MatAll_corrected.(datatemp).(loom{i}).kden=temp_kden;
end
end

figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end



%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected_crossval=[];

for fish=1:length(names)-1

    moment=[1 2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    %loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(squeeze(squeeze(CorrMatrices_mean2(moment(m),fish,:,:)))),0.75);
     MatAll_corrected_crossval(m,fish,:,:)=Mat;
     end

end

%%%% calculating the density


crossval_density=[];
for fish=1:length(names)-1
for m=1:length(moment)

temp_kden=density_und(squeeze(squeeze(MatAll_corrected_crossval(m,fish,:,:))));

crossval_density(m,fish)=temp_kden;
end
end


figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

for fish=1:length(names)-1
for m=1:length(moment)-1    
scatter(m,crossval_density(m+1,fish),'filled');
hold on;
end
end
hold off;

%% multilayer categorical network
%%%% a multilayer categorical network in some sample looms where the layers
%%%% are from different substracted means. we would expect that the
%%%% comunities would remain on every layer. 

%%%% i guess I would need to do it with the optimal gamma and omega... but
%%%% i dont know it yet. waiting for Dani's advice. i will try with the
%%%% usual gamma=1 and omega=0.5

 


%%% the clustID tags are not in the right order... I will change it
oldclusttag=[4 2 5 6 1 7];
ClustID=zeros(size(Nodes.Nod_clustID(keep)));
for i=1:length(unique(Nodes.Nod_clustID(keep)))
    
    idx_temp=find(Nodes.Nod_clustID(keep)==oldclusttag(i));
    ClustID(idx_temp)=i;
   
end
 
clustnames={'fasthab1','fasthab2','fasthab3','modhab','weakhab','inhib'};

figure;
subplot(1,2,1);histogram(Nodes.Nod_brainID(keep)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);
subplot(1,2,2);histogram(ClustID); xticks([1:6]); xticklabels(clustnames);xtickangle(45);


for moment=[1 2 3 6 11 12]

A={};

for fish=1:length(names)-1
   temp=squeeze(squeeze(CorrMatrices_mean2(moment,fish,:,:)));
   temp(isnan(temp))=0;
    A{fish}=temp;  
end
clear temp

%%

 
 gamma = 1;
 omega = 0.5; %%% this makes a big influence on how it changes through time... 
 
% N=length(A{1});
% T=length(A);
% B=spalloc(N*T,N*T,N*N*T+2*N*T);
% twomu=0;
% for s=1:T
%     k=sum(A{s});
%     twom=sum(k);
%     twomu=twomu+twom;
%     indx=[1:N]+(s-1)*N;
%     B(indx,indx)=A{s}-gamma*k'*k/twom;
% end
% twomu=twomu+2*omega*N*(T-1);
% B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
% %[S,Q] = genlouvain(B);
% [S,Q,nb_it] = iterated_genlouvain(B);
% Q = Q/twomu;
% S = reshape(S,N,T);

%%

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,(N+T)*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+T*omega*N*(T-1);
all2all = N*[(-T+1):-1,1:(T-1)];
B = B + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
[S,Q,nb_it_f] = iterated_genlouvain(B);
Q = Q/twomu;
S = reshape(S,N,T);



%% checking results

%figure;imagesc(S);colormap('jet');caxis([1 5]);


%%


time=[1 2 3 6 11]; %%% in this case time represents individual fish

low=min(min(S));
high=max(max(S));
c=jet(high);
%c=c(1:5,:);

%figure;imagesc(S);colormap(c);colorbar;

  figure;
  
  subplot(3,6,1);imagesc(S);colormap('jet');caxis([low high]);
  subplot(3,6,7);imagesc(S);colormap('jet');caxis([low high]);
  subplot(3,6,13);imagesc(S);colormap('jet');caxis([low high]);
  temp=[];
  for i=1:5
  subplot(3,6,i+1);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
%   gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],21,'off'); 
%     hold on;
    scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,S(:,time(i)),'filled');colormap('jet');caxis([low high]);
    view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  

    
subplot(3,6,i+7);histogram(Nodes.Nod_brainID(keep)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);
for j=1:max(max(S))
    hold on;histogram(Nodes.Nod_brainID(keep(find(S(:,time(i))==j))),'DisplayStyle','stairs','EdgeColor',c(j,:));
end

subplot(3,6,i+13);histogram(ClustID); xticks([1:6]); xticklabels(clustnames);xtickangle(45);
for j=1:max(max(S))
hold on;histogram(ClustID(find(S(:,time(i))==j)),'DisplayStyle','stairs','EdgeColor',c(j,:)); 
end
      
  end
  
end


%%
%%%% I could do a participation substraction to the original matrix...
%%%% probalby not necessary 

%%

save('crossvalidation_1fishout_f20_results.mat');


%% testing something

%%%% i noticed something interesting... the communities per loom across all
%%%% the cross validation matrices follow an interesting pattern. First
%%%% they seem to agregate by brain region, then they all merge with the
%%%% first loom to then gradually separete to form something similar to the
%%%% main clusters. 

%%%% I will try to represent it in figure. 


Comm_minusFish=[];
counter=1;
for moment=[1 2 3 6 11 12]

A={};

for fish=1:length(names)-1
   temp=squeeze(squeeze(CorrMatrices_mean2(moment,fish,:,:)));
   temp(isnan(temp))=0;
    A{fish}=temp;  
end
clear temp

%%

 gamma = 1.05;
 omega = 0.5; %%% this makes a big influence on how it changes through time... 
 
%%

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,(N+T)*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+T*omega*N*(T-1);
all2all = N*[(-T+1):-1,1:(T-1)];
B = B + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
[S,Q,nb_it_f] = iterated_genlouvain(B);
Q = Q/twomu;
S = reshape(S,N,T);

Comm_minusFish=cat(3,Comm_minusFish,S);

time=[1 2 3 6 11]; %%% in this case time represents individual fish

low=min(min(S));
high=max(max(S));
c=jet(high);
%c=c(1:5,:);

%figure;imagesc(S);colormap(c);colorbar;

  figure;
  
  subplot(3,6,1);imagesc(S);colormap('jet');caxis([low high]);
  subplot(3,6,7);imagesc(S);colormap('jet');caxis([low high]);
  subplot(3,6,13);imagesc(S);colormap('jet');caxis([low high]);
  temp=[];
  for i=1:5
  subplot(3,6,i+1);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
%   gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],21,'off'); 
%     hold on;
    scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,S(:,time(i)),'filled');colormap('jet');caxis([low high]);
    view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  

    
subplot(3,6,i+7);histogram(Nodes.Nod_brainID(keep)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);
for j=1:max(max(S))
    hold on;histogram(Nodes.Nod_brainID(keep(find(S(:,time(i))==j))),'DisplayStyle','stairs','EdgeColor',c(j,:));
end

subplot(3,6,i+13);histogram(ClustID); xticks([1:6]); xticklabels(clustnames);xtickangle(45);
for j=1:max(max(S))
hold on;histogram(ClustID(find(S(:,time(i))==j)),'DisplayStyle','stairs','EdgeColor',c(j,:)); 
end
      
  end
  
end

%save('Comm_minusFish_test_g10o05.mat','Comm_minusFish','CorrMatrices_mean2','gamma','omega');
%save('Comm_minusFish_test_g105o05.mat','Comm_minusFish','CorrMatrices_mean2','gamma','omega');

%%%%% trying to put it together
%load('Comm_minusFish_test_g105o05.mat','Comm_minusFish','gamma','omega');
  figure;
  sgtitle(strcat('gamma=',num2str(gamma),'omega=',num2str(omega)));
moment=[1 2 3 6 11 12];
for i=1:length(moment)
name=strcat('moment_',num2str(moment(i)));

S_good=Comm_minusFish(:,:,i);


low=min(min(S_good));
high=max(max(S_good));
c=jet(high);
time=randperm(11,1);%[1 2 6 8 9];%%%% in this case this is actually the fish
  
  subplot(4,6,i);imagesc(S_good);colormap('jet');caxis([low high]);
  temp=[];  
  subplot(4,6,i+6);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
    scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,S_good(:,time),'filled');colormap('jet');caxis([low high]);
    view(-90,90); 
    hold off;  
     
%   subplot(5,6,i+12);
%   plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%   hold on;
%    gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],20,'off');     
%      view(-90,90);
%     %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
%     hold off;  
      
subplot(4,6,i+12);histogram(Nodes.Nod_brainID(keep)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);
for j=1:max(max(S_good))
    hold on;histogram(Nodes.Nod_brainID(keep(find(S_good(:,time)==j))),'DisplayStyle','stairs','EdgeColor',c(j,:));
end

subplot(4,6,i+18);histogram(ClustID); xticks([1:6]); xticklabels(clustnames);xtickangle(45);
for j=1:max(max(S_good))
hold on;histogram(ClustID(find(S_good(:,time)==j)),'DisplayStyle','stairs','EdgeColor',c(j,:)); 
end
   
end


%%
%%%% what if I try the consensus?


S_cons=struct;

for moment=[1 2 3 6 11 12]

A={};

for fish=1:length(names)-1
   temp=squeeze(squeeze(CorrMatrices_mean2(moment,fish,:,:)));
   temp(isnan(temp))=0;
    A{fish}=temp;  
end
clear temp


%%%% making a multidimensional structure where to store things
S_test=[];

for test=1:1000
     
gamma = 1;
omega = 0.5;  

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,(N+T)*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+T*omega*N*(T-1);
all2all = N*[(-T+1):-1,1:(T-1)];
B = B + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);

[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
S = reshape(S,N,T);

S_test=cat(3,S_test,S);
  
  
end


%% trying consensus_iterative 


S_good=[];
for i=1:size(S_test,2)
   C=S_test(:,i,:); 
   C=squeeze(C);
   [S2, Q2, X_new3, qpc] = consensus_iterative(C');
    S_good(:,i)=S2(i,:);
end

name=strcat('moment_',num2str(moment));

S_cons.(name).S_test=S_test;
S_cons.(name).S_cons=S_good;


low=min(min(S_good));
high=max(max(S_good));

time=[1 2 6 8 9];%%%% in this case this is actually the fish

  figure;
  
  subplot(1,7,1);imagesc(S_good);colormap('jet');caxis([low high]);
  temp=[];
  for i=1:5
  subplot(1,7,i+1);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
%    gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],21,'off'); 
%      hold on;
    scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,S_good(:,time(i)),'filled');colormap('jet');caxis([low high]);
    view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  
     
  end
  subplot(1,7,7);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
   gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],20,'off');     
     view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  
    
end

%save('Cons_Comm_minusFish_test_g09o08.mat','S_cons','CorrMatrices_mean2','gamma','omega');
%save('Cons_Comm_minusFish_test_g11o01.mat','S_cons','CorrMatrices_mean2','gamma','omega');
%save('Cons_Comm_minusFish_test_g10o05.mat','S_cons','CorrMatrices_mean2','gamma','omega');

%%%%% trying to put it together
load('Cons_Comm_minusFish_test_g10o05.mat','S_cons','gamma','omega');
  figure;
moment=[1 2 3 6 11 12];
for i=1:length(moment)
name=strcat('moment_',num2str(moment(i)));

S_good=S_cons.(name).S_cons;


low=min(min(S_good));
high=max(max(S_good));
c=jet(high);
time=randperm(11,1);%[1 2 6 8 9];%%%% in this case this is actually the fish
  
  subplot(4,6,i);imagesc(S_good);colormap('jet');caxis([low high]);
  temp=[];  
  subplot(4,6,i+6);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
  hold on;
    scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),20,S_good(:,time),'filled');colormap('jet');caxis([low high]);
    view(-90,90); 
    hold off;  
     
%   subplot(5,6,i+12);
%   plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%   hold on;
%    gscatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),ClustID,'gggbrm',[],20,'off');     
%      view(-90,90);
%     %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
%     hold off;  
      
subplot(4,6,i+12);histogram(Nodes.Nod_brainID(keep)); xticks([1:9]); xticklabels(RegionList);xtickangle(45);
for j=1:max(max(S_good))
    hold on;histogram(Nodes.Nod_brainID(keep(find(S_good(:,time)==j))),'DisplayStyle','stairs','EdgeColor',c(j,:));
end

subplot(4,6,i+18);histogram(ClustID); xticks([1:6]); xticklabels(clustnames);xtickangle(45);
for j=1:max(max(S_good))
hold on;histogram(ClustID(find(S_good(:,time)==j)),'DisplayStyle','stairs','EdgeColor',c(j,:)); 
end
   
end

