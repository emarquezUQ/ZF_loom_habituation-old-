

%%%%%% this script is to try to find the optimal gamma and omega for the
%%%%%% control group (to use as reference for the others). I will try to
%%%%%% find the values that give a higher difference with a nodal null
%%%%%% model and has low variance. 

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



%%

group=3; %%% for controls
A={};

for loom=1:21
  
   temp=Data_corrMat4.(groupnames{group}).Mean_corrMat{1,loom}(keepFmr1,keepFmr1);
   temp(isnan(temp))=0;
    A{loom}=temp;  
end
clear temp

%%
Qopt2=[];

gamma = 0.1:0.1:3.7;
omega = 0.1:0.1:3.7;

for test=1:100
Qmat2=[];
for i=1:length(gamma)   
parfor j=1:length(omega)

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma(i)*k'*k/twom;
end
twomu=twomu+2*omega(j)*N*(T-1);
B = B + omega(j)*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
%S = reshape(S,N,T);

Qmat2(i,j)=Q;
   
end
end

Qopt2=cat(3,Qopt2,Qmat2);

end

meanQopt2=mean(Qopt2,3);
varQopt2=var(Qopt2,0,3);

figure;
subplot(1,2,1);imagesc(meanQopt2); colorbar;
subplot(1,2,2);imagesc(varQopt2); colorbar;


%% nodal null model

randOrder = randperm(99,99);

%unique(randOrder)


An_null={};

for loom=1:21
    
    randOrder = randperm(99,99);
    
   temp=Data_corrMat4.(groupnames{group}).Mean_corrMat{1,loom}(keepFmr1,keepFmr1);
   temp(isnan(temp))=0;
    An_null{loom}=temp;  
end
clear temp



Qopt2_null=[];

gamma = 0.1:0.1:3.7;
omega = 0.1:0.1:3.7;

for test=1:100
Qmat2_null=[];
for i=1:length(gamma)   
parfor j=1:length(omega)

N=length(An_null{1});
T=length(An_null);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(An_null{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=An_null{s}-gamma(i)*k'*k/twom;
end
twomu=twomu+2*omega(j)*N*(T-1);
B = B + omega(j)*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
%S = reshape(S,N,T);

Qmat2_null(i,j)=Q;
   
end
end

Qopt2_null=cat(3,Qopt2_null,Qmat2_null);

end

meanQopt2_null=mean(Qopt2_null,3);
varQopt2_null=var(Qopt2_null,0,3);

figure;
subplot(1,2,1);imagesc(meanQopt2_null); colorbar;
subplot(1,2,2);imagesc(varQopt2_null); colorbar;


figure;imagesc(meanQopt2-meanQopt2_null);colorbar;

save('QoptFMR1_WT.mat','A','Qmat2','Qopt2','meanQopt2','varQopt2','An_null','randOrder','Qmat2_null','Qopt2_null','meanQopt2_null','varQopt2_null');

%%

%%% trying to find the minimum variance of Q and maximum difference in Q
%%% with the null model:

%%% visually

meanQopt2_dif=meanQopt2-meanQopt2_null;
low=min(min(meanQopt2_dif));
high=max(max(meanQopt2_dif));

meanQopt2_dif_norm=meanQopt2_dif-low;
meanQopt2_dif_norm=meanQopt2_dif_norm/(high-low);


varQopt2_norm=varQopt2;
low=min(min(varQopt2_norm));
high=max(max(varQopt2_norm));

varQopt2_norm=varQopt2_norm-low;
varQopt2_norm=varQopt2_norm/(high-low);

meanQopt2_dif_norm-varQopt2_norm;

%%%% with stats


x=mean(mean(meanQopt2_dif));

y=std(mean(meanQopt2_dif));

good_idx=find(meanQopt2_dif>(x+2*y));

good=zeros(size(meanQopt2_dif));

good(good_idx)=1;

figure;imagesc(good);

x2=mean(mean(varQopt2));

y2=std(mean(varQopt2));

good_idx=find(varQopt2<(x2-1*y2));

good2=zeros(size(varQopt2));

good2(good_idx)=1;

figure;imagesc(good2);

thegood=zeros(size(good));
thegood(intersect(find(good),find(good2)))=1;

figure;imagesc(thegood);
