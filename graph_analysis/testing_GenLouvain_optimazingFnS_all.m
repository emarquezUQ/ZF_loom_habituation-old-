

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

QoptAll=[];

for data=1:4

A={};

for loom=1:31
   temp=Data_corrMat2.(datasets(data,:)).Mean_corrMat{1,loom}(keep,keep);
   temp(isnan(temp))=0;
    A{loom}=temp;  
end
clear temp


%%
Qopt2=[];

gamma = 0.1:0.1:1.6;
omega = 0.1:0.1:1;

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

for loom=1:31
    
    randOrder = randperm(99,99);
    
   
   temp=Data_corrMat2.(datasets(data,:)).Mean_corrMat{1,loom}(keep,keep);
   temp(isnan(temp))=0;
    An_null{loom}=temp;  
end
clear temp



Qopt2_null=[];

gamma = 0.1:0.1:1.6;
omega = 0.1:0.1:1;

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



QoptAll.(datasets(data,:)).A=A;
QoptAll.(datasets(data,:)).Qopt2=Qopt2;
QoptAll.(datasets(data,:)).meanQopt2=meanQopt2;
QoptAll.(datasets(data,:)).varQopt2=varQopt2;

QoptAll.(datasets(data,:)).An_null=An_null;
QoptAll.(datasets(data,:)).Qopt2_null=Qopt2_null;
QoptAll.(datasets(data,:)).meanQopt2_null=meanQopt2_null;
QoptAll.(datasets(data,:)).varQopt2_null=varQopt2_null;
end


save('QoptAll.mat','QoptAll');
