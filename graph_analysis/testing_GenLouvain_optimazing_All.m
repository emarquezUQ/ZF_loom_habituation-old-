%%%% this script is to look for a good gamma and omega values for all the
%%%% datasets... see if we can find a good number. 

load('QoptFMR1_WT.mat');


meanQopt2_dif2=meanQopt2-meanQopt2_null;
meanQopt2_dif2short=meanQopt2_dif2(1:16,1:10);

varQopt2short=varQopt2(1:16,1:10);

%% gilles method
figure;imagesc(meanQopt2_dif2short);
figure;imagesc(varQopt2short);

relative_var=-(varQopt2short-max(varQopt2short(:)));
figure;imagesc(relative_var);

thegood=meanQopt2_dif2short.*relative_var;

figure;imagesc(thegood);

thebest=thegood;thebest(thegood<0)=0;figure;imagesc(thebest);

figure;imagesc(thebest);
%%%% I will also look at the higher 1 or 2SD


x=mean(thebest(:));
y=std(thebest(:));

good_idx=find(thebest>(x+1*y));
%good_idx=find(thebest>x);

good_fmr1_WT=zeros(size(thebest));

good_fmr1_WT(good_idx)=1;
figure;imagesc(good_fmr1_WT);

%% now for the Fns datasets

load('QoptAll.mat');
datasets=fieldnames(QoptAll);
thebest_FnS=[];
good1SD=[];
good2SD=[];
for data=1:4
    temp_dif=QoptAll.(datasets{data}).meanQopt2-QoptAll.(datasets{data}).meanQopt2_null;
    temp_var=QoptAll.(datasets{data}).varQopt2;
    
     %figure;imagesc(temp_dif);
    
    temp_relative_var=-(temp_var-max(temp_var(:)));
    
    temp_thegood=temp_dif.*temp_relative_var;

    figure;imagesc(temp_thegood);

    temp_thebest=temp_thegood;temp_thebest(temp_thegood<0)=0;figure;imagesc(temp_thebest);
    
    thebest_FnS=cat(3,thebest_FnS,temp_thebest);
    
    x=mean(temp_thebest(:));
    y=std(temp_thebest(:));

    good_idx=find(temp_thebest>(x+1*y));
    good_idx2=find(temp_thebest>(x+2*y));

    temp_good=zeros(size(temp_thebest));

    temp_good(good_idx)=1;
    figure;imagesc(temp_good);
    
    good1SD=vertcat(good1SD,good_idx);
    good2SD=vertcat(good2SD,good_idx2);
    
end


thebest_All=cat(3,thebest_FnS,thebest);

meanBest=mean(thebest_All,3);

figure;imagesc(meanBest);

figure;imagesc(meanBest(10:15,3:8));

%%

good1SD=unique(good1SD);

good2SD=unique(good2SD);

 temp_good=zeros(size(meanBest));

    temp_good(union(good1SD,good_idx))=1;
    figure;imagesc(temp_good);
    
   