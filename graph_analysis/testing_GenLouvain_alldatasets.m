

load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

keepFmr1=load('graphs_fmr1Exp.mat','keep');
keepFmr1=keepFmr1.keep;

%%


%%%% just checking how many nodes I have per brain region and cluster

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


%% doing it for a multilayer matrix

%%% first for the averaged matrix for f20 in time 

%%% making the cell array

for data=1:4

A={};

for loom=1:31
   temp=Data_corrMat2.(datasets(data,:)).Mean_corrMat{1,loom}(keep,keep);
   temp(isnan(temp))=0;
    A{loom}=temp;  
end
clear temp

%%


gamma = 0.9;
omega = 0.8; %%% this makes a big influence... 

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
%[S,Q] = genlouvain(B);
[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu
S = reshape(S,N,T);

%% checking results

%figure;imagesc(S);colormap('jet');caxis([1 5]);


%%


time=[1 2 6 11 12];

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


for group=1:3
A={};

for loom=1:21
  
   temp=Data_corrMat4.(groupnames{group}).Mean_corrMat{1,loom}(keepFmr1,keepFmr1);
   temp(isnan(temp))=0;
    A{loom}=temp;  
end
clear temp

%%


gamma = 1;
omega = 0.8; %%% this makes a big influence... 

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
%[S,Q] = genlouvain(B);
[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu
S = reshape(S,N,T);

%% checking results

%figure;imagesc(S);colormap('jet');caxis([1 5]);


%%


time=[1 2 6 11 12];

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
    scatter(Nodes.Nod_coor(keep(keepFmr1),1),Nodes.Nod_coor(keep(keepFmr1),2),20,S(:,time(i)),'filled');colormap('jet');caxis([low high]);
    view(-90,90);
    %title('o=0.8/test=',num2str(test),'/time=',num2str(time(i)));
    hold off;  

    
subplot(3,6,i+7);histogram(Nodes.Nod_brainID(keep(keepFmr1))); xticks([1:9]); xticklabels(RegionList);xtickangle(45);
for j=1:max(max(S))
    hold on;histogram(Nodes.Nod_brainID(keep(keepFmr1(find(S(:,time(i))==j)))),'DisplayStyle','stairs','EdgeColor',c(j,:));
end

subplot(3,6,i+13);histogram(ClustID(keepFmr1)); xticks([1:6]); xticklabels(clustnames);xtickangle(45);
for j=1:max(max(S))
hold on;histogram(ClustID(keepFmr1(find(S(:,time(i))==j))),'DisplayStyle','stairs','EdgeColor',c(j,:)); 
end
      
  end
  
end
