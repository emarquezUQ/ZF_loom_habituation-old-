
%%%%% testing multilayer nodal and temporal null models
%%%% it seems that the nodal model makes very different resutls. the
%%%% temporal model differs a bit but less... the differences are more
%%%% clear a low omega values

load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')


A={};

for loom=1:31
   temp=Data_corrMat2.f20.Mean_corrMat{1,loom}(keep,keep);
   temp(isnan(temp))=0;
    A{loom}=temp;  
end
clear temp


gamma =1;
omega = 2;

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
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu
S = reshape(S,N,T);


%% checking results

figure;imagesc(S);colormap('jet');%caxis([1 5]);

%%
%%%% making a null nodal model of A

randOrder = randperm(99,99);

%unique(randOrder)


An_null={};

for loom=1:31
    
    randOrder = randperm(99,99);
    
   temp=Data_corrMat2.f20.Mean_corrMat{1,loom}(keep(randOrder),keep(randOrder));
   temp(isnan(temp))=0;
    An_null{loom}=temp;  
end
clear temp



gamma = 1;
omega = 2;

N=length(An_null{1});
T=length(An_null);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(An_null{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=An_null{s}-gamma*k'*k/twom;
end
twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
S_null = reshape(S,N,T);


figure;
imagesc(S_null);colormap('jet');%caxis([1 5]);



%% trying also with a temporal null model



%%%%% this is changing the order of the layers


randOrder2 = randperm(31,31);

%unique(randOrder2)

At_null={};

for loom=1:31
    
   temp=Data_corrMat2.f20.Mean_corrMat{1,randOrder2(loom)}(keep,keep);
   temp(isnan(temp))=0;
    At_null{loom}=temp;  
end
clear temp



gamma = 1;
omega = 2;

N=length(At_null{1});
T=length(At_null);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(At_null{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=At_null{s}-gamma*k'*k/twom;
end
twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
S_null2 = reshape(S,N,T);


figure;
imagesc(S_null2);colormap('jet');%caxis([1 5]);
