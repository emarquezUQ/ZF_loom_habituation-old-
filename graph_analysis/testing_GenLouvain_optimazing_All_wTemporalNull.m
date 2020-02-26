%%%% trying to optimize with a temporal null model for all datasets 


load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

keepFmr1=load('graphs_fmr1Exp.mat','keep');
keepFmr1=keepFmr1.keep;


%%
%%%% with fasthab subtypes
oldclusttag=[4 2 5 6 1 7];
ClustID=zeros(size(Nodes.Nod_clustID(keep)));
for i=1:length(unique(Nodes.Nod_clustID(keep)))
    
    idx_temp=find(Nodes.Nod_clustID(keep)==oldclusttag(i));
    ClustID(idx_temp)=i;
   
end
 
clustnames={'fasthab1','fasthab2','fasthab3','modhab','weakhab','inhib'};

%%%% with fasthab merged
oldclusttag=[4 2 5 6 1 7];
ClustID_CL4=zeros(size(Nodes.Nod_clustID(keep)));
for i=1:length(unique(Nodes.Nod_clustID(keep)))
    
    idx_temp=find(Nodes.Nod_clustID(keep)==oldclusttag(i));
    
    if i<4    
    ClustID_CL4(idx_temp)=1;
    elseif i==4
   ClustID_CL4(idx_temp)=2;
   elseif i==5
   ClustID_CL4(idx_temp)=3;
   elseif i==6
   ClustID_CL4(idx_temp)=4;
    end
end
 
clustnames_CL4={'fasthab','modhab','weakhab','inhib'};


%%%% for the fmr1 dataset:

group=3; %%% for controls
A={};

for loom=1:21
  
   temp=Data_corrMat4.(groupnames{group}).Mean_corrMat{1,loom}(keepFmr1,keepFmr1);
   temp(isnan(temp))=0;
    A{loom}=temp;  
end
clear temp


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

%% temporal null model


randOrder2 = randperm(21,21);

%unique(randOrder2)


At_null={};

for loom=1:21
    temp=Data_corrMat4.(groupnames{group}).Mean_corrMat{1,randOrder2(loom)}(keepFmr1,keepFmr1);  
   temp(isnan(temp))=0;
    At_null{loom}=temp;  
end
clear temp



Qopt2_null2=[];

gamma = 0.1:0.1:1.6;
omega = 0.1:0.1:1;

for test=1:100
Qmat2_null2=[];
for i=1:length(gamma)   
parfor j=1:length(omega)

N=length(At_null{1});
T=length(At_null);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(At_null{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=At_null{s}-gamma(i)*k'*k/twom;
end
twomu=twomu+2*omega(j)*N*(T-1);
B = B + omega(j)*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
%S = reshape(S,N,T);

Qmat2_null2(i,j)=Q;
   
end
end

Qopt2_null2=cat(3,Qopt2_null2,Qmat2_null2);

end

meanQopt2_null2=mean(Qopt2_null2,3);
varQopt2_null2=var(Qopt2_null2,0,3);

figure;
subplot(1,2,1);imagesc(meanQopt2_null2); colorbar;
subplot(1,2,2);imagesc(varQopt2_null2); colorbar;


figure;imagesc(meanQopt2-meanQopt2_null2);colorbar;

%%


meanQopt2_dif2=meanQopt2-meanQopt2_null2;


%% gilles method
figure;imagesc(meanQopt2_dif2);
figure;imagesc(varQopt2);

relative_var=-(varQopt2-max(varQopt2(:)));
figure;imagesc(relative_var);

thegood=meanQopt2_dif2.*relative_var;

figure;imagesc(thegood);

thebest_WT=thegood;thebest_WT(thegood<0)=0;figure;imagesc(thebest_WT);

figure;imagesc(thebest_WT);
%%%% I will also look at the higher 1 or 2SD


x=mean(thebest_WT(:));
y=std(thebest_WT(:));

good_idx=find(thebest_WT>(x+1*y));
%good_idx=find(thebest>x);

good_fmr1_WT=zeros(size(thebest_WT));

good_fmr1_WT(good_idx)=1;
figure;imagesc(good_fmr1_WT);


%% for the FnS dataset

%%%% I already have the normal ones an the nodal null model
load('QoptAll.mat','QoptAll');

randOrder2 = randperm(31,31);
for data=1:4


%% temporal null model

%unique(randOrder)

At_null={};

for loom=1:31
        
   temp=Data_corrMat2.(datasets(data,:)).Mean_corrMat{1,randOrder2(loom)}(keep,keep);
   temp(isnan(temp))=0;
    At_null{loom}=temp;  
end
clear temp



Qopt2_null=[];

gamma = 0.1:0.1:1.6;
omega = 0.1:0.1:1;

for test=1:100
Qmat2_null=[];
for i=1:length(gamma)   
parfor j=1:length(omega)

N=length(At_null{1});
T=length(At_null);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(At_null{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=At_null{s}-gamma(i)*k'*k/twom;
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


QoptAll.(datasets(data,:)).At_null=At_null;
QoptAll.(datasets(data,:)).Qopt2_null_T=Qopt2_null;
QoptAll.(datasets(data,:)).meanQopt2_null_T=meanQopt2_null;
QoptAll.(datasets(data,:)).varQopt2_null_T=varQopt2_null;
end


%save('QoptAll_null_Temp.mat','QoptAll','gamma','omega');
%save('QoptAll_null_Temp.mat');

load('QoptAll_null_Temp.mat');

Allthebest=[];
Allverybest=[];
for data=1:4
    temp_dif=QoptAll.(datasets(data,:)).meanQopt2-QoptAll.(datasets(data,:)).meanQopt2_null_T;
    
    figure;imagesc(temp_dif);
           
temp_relative_var=-(QoptAll.(datasets(data,:)).varQopt2-max(QoptAll.(datasets(data,:)).varQopt2(:)));
figure;imagesc(temp_relative_var);

thegood=temp_dif.*temp_relative_var;

%figure;imagesc(thegood);

thebest=thegood;thebest(thegood<0)=0;

figure;imagesc(thebest);

Allthebest=cat(3,Allthebest,thebest);

%%%% I will also look at the higher 1 or 2SD

x=mean(thebest(:));
y=std(thebest(:));

good_idx=find(thebest>(x+1*y));
%good_idx=find(thebest>x);

very_best=zeros(size(thebest));

very_best(good_idx)=1;

Allverybest=cat(3,Allverybest,Allverybest);

figure;imagesc(very_best);
      
end

%%

load('QoptFMR1_All.mat');


Allthebest_fmr1=[];
Allverybest_fmr1=[];
for group=1:3
    temp_dif=QoptFMR1_All.(groupnames{group}).meanQopt2-QoptFMR1_All.(groupnames{group}).meanQopt2_null2;
    
    %figure;imagesc(temp_dif);
           
temp_relative_var=-(QoptFMR1_All.(groupnames{group}).varQopt2-max(QoptFMR1_All.(groupnames{group}).varQopt2(:)));
figure;imagesc(temp_relative_var);

thegood=temp_dif.*temp_relative_var;

%figure;imagesc(thegood);

thebest=thegood;thebest(thegood<0)=0;

figure;imagesc(thebest);

Allthebest_fmr1=cat(3,Allthebest_fmr1,thebest);

%%%% I will also look at the higher 1 or 2SD

x=mean(thebest(:));
y=std(thebest(:));

good_idx=find(thebest>(x+1*y));
%good_idx=find(thebest>x);

very_best=zeros(size(thebest));

very_best(good_idx)=1;

Allverybest_fmr1=cat(3,Allverybest_fmr1,very_best);

figure;imagesc(very_best);
      
end


%%

Allthebest_all=cat(3,Allthebest,Allthebest_fmr1);

%Allthebest_all=Allthebest_fmr1;

Opt=mean(Allthebest_all,3);

figure;
imagesc(Opt);

x=mean(Opt(:));
y=std(Opt(:));

%good_idx=find(Opt>x);
%good_idx=find(Opt>(x+1*y));
good_idx=find(Opt>(x+2*y));

Opt_best=zeros(size(Opt));

Opt_best(good_idx)=1;

figure;imagesc(Opt_best);

[g,o]=find(Opt_best);

%% 

%%%% to take the points that for each dataset passed the threshold that I
%%%% choose. 

all_good_idx=[];
for i=1:size(Allthebest_all,3)
    
    temp=Allthebest_all(:,:,i);
    x=mean(temp(:));
    y=std(temp(:));

    good_idx=find(temp>(x+1*y));
    %good_idx=find(Opt>x);

    all_good_idx=union(all_good_idx,good_idx);
    temp_best=zeros(size(temp));

    temp_best(good_idx)=1;

    figure;imagesc(temp_best);
end

temp_best=zeros(size(temp));

    temp_best(all_good_idx)=1;

    figure;imagesc(temp_best);
    
    
    [g,o]=find(temp_best);
    
    %% making a figure
    
    
 %%%% with f20 Qdif, f20 Q variance, the product of per element multiplication and then the average of all datasets
 
 figure;
 subplot(1,4,1);
 data=1;
 temp_dif=QoptAll.(datasets(data,:)).meanQopt2-QoptAll.(datasets(data,:)).meanQopt2_null_T;
 imagesc(temp_dif);colormap(inferno);colorbar;
 
 subplot(1,4,2);
 imagesc(QoptAll.(datasets(data,:)).varQopt2);colormap(inferno);colorbar;
 
 subplot(1,4,3);
temp_relative_var=-(QoptAll.(datasets(data,:)).varQopt2-max(QoptAll.(datasets(data,:)).varQopt2(:)));
thegood=temp_dif.*temp_relative_var;
thebest=thegood;thebest(thegood<0)=0;
imagesc(thebest);colormap(inferno);colorbar;

subplot(1,4,4);
Opt=mean(Allthebest_all,3);
imagesc(Opt);colormap(inferno);colorbar;

 %%%%% adjust shape and save as SVG   