



%%%%% this script is to test histograms of the substraction matrices. 
%%%%% I tried previously with the f20-f60 to test but I dont think I used
%%%%% the right approaches. So what I did is to also use the
%%%%% crossvalidation_1fishout_fmr1.m script to be able to do
%%%%% stats on the # of edges above 0.75 with the multiple matrices

load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

datasets=['f20'; 'f60'; 's20'; 's60'];
keepFmr1=load('graphs_fmr1Exp.mat','keep');
keepFmr1=keepFmr1.keep;
%cbrewer()

[RdBu]=cbrewer('div','RdBu',101);


%% first the matrices
%%%% to plot the 1st-3rd, 10th and 11th 
counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 group=3;
     grouptemp=groupnames{group,:};
subplot(2,5,counter);imagesc(Data_corrMat4.(grouptemp).Mean_corrMat{1,2}(keepFmr1,keepFmr1)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for first loom
subplot(2,5,counter+1);imagesc(Data_corrMat4.(grouptemp).Mean_corrMat{1,3}(keepFmr1,keepFmr1));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for the 2nd loom
subplot(2,5,counter+2);imagesc(Data_corrMat4.(grouptemp).Mean_corrMat{1,4}(keepFmr1,keepFmr1));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for the 3rd loom
subplot(2,5,counter+3);imagesc(Data_corrMat4.(grouptemp).Mean_corrMat{1,11}(keepFmr1,keepFmr1));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for 10th loom
subplot(2,5,counter+4);imagesc(Data_corrMat4.(grouptemp).Mean_corrMat{1,12}(keepFmr1,keepFmr1));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for 11th loom
title(grouptemp);
group=2;
     grouptemp=groupnames{group,:};
subplot(2,5,counter+5);imagesc(Data_corrMat4.(grouptemp).Mean_corrMat{1,2}(keepFmr1,keepFmr1)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for first loom
subplot(2,5,counter+6);imagesc(Data_corrMat4.(grouptemp).Mean_corrMat{1,3}(keepFmr1,keepFmr1));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for the 2nd loom
subplot(2,5,counter+7);imagesc(Data_corrMat4.(grouptemp).Mean_corrMat{1,4}(keepFmr1,keepFmr1));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for the 3rd loom
subplot(2,5,counter+8);imagesc(Data_corrMat4.(grouptemp).Mean_corrMat{1,11}(keepFmr1,keepFmr1));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for 10th loom
subplot(2,5,counter+9);imagesc(Data_corrMat4.(grouptemp).Mean_corrMat{1,12}(keepFmr1,keepFmr1));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for 11th loom
title(grouptemp);

%%%%% how doest it look if is not thresholded? 
sub1_weigths=Data_corrMat4.(grouptemp).Mean_corrMat{1,2}(keepFmr1,keepFmr1)-Data_corrMat4.(groupnames{group-1,:}).Mean_corrMat{1,2}(keepFmr1,keepFmr1);
sub10_weigths=Data_corrMat4.(grouptemp).Mean_corrMat{1,11}(keepFmr1,keepFmr1)-Data_corrMat4.(groupnames{group-1,:}).Mean_corrMat{1,11}(keepFmr1,keepFmr1);
sub11_weigths=Data_corrMat4.(grouptemp).Mean_corrMat{1,12}(keepFmr1,keepFmr1)-Data_corrMat4.(groupnames{group-1,:}).Mean_corrMat{1,12}(keepFmr1,keepFmr1);

%%%% is a bit caothic... to many grouppoints
figure;
subplot(1,3,1);histogram(sub1_weigths);
subplot(1,3,2);histogram(sub10_weigths);
subplot(1,3,3);histogram(sub11_weigths);

nanmean(sub1_weigths(:))
nanmean(sub10_weigths(:))
nanmean(sub11_weigths(:))





%%%% this is to get the thresholded group
MatAll_corrected2=struct;
for group=1:3
    grouptemp=groupnames{group,:};
    moment=[2 3 4 5 6 7 8 9 10 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 6 7 8 9 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat4.(grouptemp).Mean_corrMat{1,moment(m)}(keepFmr1,keepFmr1)),0.75);
     MatAll_corrected2.(grouptemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end



counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 group=3;
     grouptemp=groupnames{group,:};
     loom=fieldnames(MatAll_corrected2.(grouptemp));
subplot(2,3,counter);imagesc(MatAll_corrected2.(grouptemp).(loom{1}).Mat-MatAll_corrected2.(groupnames{group-1,:}).(loom{1}).Mat); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for pre loom
subplot(2,3,counter+1);imagesc(MatAll_corrected2.(grouptemp).(loom{10}).Mat-MatAll_corrected2.(groupnames{group-1,:}).(loom{10}).Mat);pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for 10th loom
subplot(2,3,counter+2);imagesc(MatAll_corrected2.(grouptemp).(loom{11}).Mat-MatAll_corrected2.(groupnames{group-1,:}).(loom{11}).Mat);pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for 11th loom
title(grouptemp);


sub1=MatAll_corrected2.(grouptemp).(loom{1}).Mat-MatAll_corrected2.(groupnames{group-1,:}).(loom{1}).Mat;
sub2=MatAll_corrected2.(grouptemp).(loom{2}).Mat-MatAll_corrected2.(groupnames{group-1,:}).(loom{2}).Mat;
sub3=MatAll_corrected2.(grouptemp).(loom{3}).Mat-MatAll_corrected2.(groupnames{group-1,:}).(loom{3}).Mat;
sub10=MatAll_corrected2.(grouptemp).(loom{10}).Mat-MatAll_corrected2.(groupnames{group-1,:}).(loom{10}).Mat;
sub11=MatAll_corrected2.(grouptemp).(loom{11}).Mat-MatAll_corrected2.(groupnames{group-1,:}).(loom{11}).Mat;


figure;
subplot(1,5,1);histogram(sub1(find(sub1(:))));
subplot(1,5,2);histogram(sub2(find(sub2(:))));
subplot(1,5,3);histogram(sub3(find(sub3(:))));
subplot(1,5,4);histogram(sub10(find(sub10(:))));
subplot(1,5,5);histogram(sub11(find(sub11(:))));

nanmedian(sub1(find(sub1(:))))
nanmedian(sub1(find(sub2(:))))
nanmedian(sub1(find(sub3(:))))
nanmedian(sub10(find(sub10(:))))
nanmedian(sub11(find(sub11(:))))

%%
%%%%% i will try to represent it first with the substraction, as is the way
%%%%% we are presenting it in the schematic graph on the paper
edges=[-1:0.2:1];
figure;set(gcf, 'Position',  [100, 100, 1400, 300]);
subplot(1,5,1);histogram(sub1(find(sub1>0)),edges,'FaceColor','b');hold on;histogram(sub1(find(sub1<0)),edges,'FaceColor','r');title('1st Loom');
subplot(1,5,2);histogram(sub2(find(sub2>0)),edges,'FaceColor','b');hold on;histogram(sub2(find(sub2<0)),edges,'FaceColor','r');title('2nd Loom');
subplot(1,5,3);histogram(sub3(find(sub3>0)),edges,'FaceColor','b');hold on;histogram(sub3(find(sub3<0)),edges,'FaceColor','r');title('3rd Loom');
subplot(1,5,4);histogram(sub10(find(sub10>0)),edges,'FaceColor','b');hold on;histogram(sub10(find(sub10<0)),edges,'FaceColor','r');title('10th Loom');
subplot(1,5,5);histogram(sub11(find(sub11>0)),edges,'FaceColor','b');hold on;histogram(sub11(find(sub11<0)),edges,'FaceColor','r');title('11th Loom');

saveas(gcf,'histogram_WTminusFmr1_edges.svg');

%%
%%%% and now just total count of edges above 0.75
highcorrWTvsFmr1(1,1)=length(find(MatAll_corrected2.(grouptemp).(loom{1}).Mat));
highcorrWTvsFmr1(1,2)=length(find(MatAll_corrected2.(groupnames{group-1,:}).(loom{1}).Mat));
highcorrWTvsFmr1(2,1)=length(find(MatAll_corrected2.(grouptemp).(loom{2}).Mat));
highcorrWTvsFmr1(2,2)=length(find(MatAll_corrected2.(groupnames{group-1,:}).(loom{2}).Mat));
highcorrWTvsFmr1(3,1)=length(find(MatAll_corrected2.(grouptemp).(loom{3}).Mat));
highcorrWTvsFmr1(3,2)=length(find(MatAll_corrected2.(groupnames{group-1,:}).(loom{3}).Mat));
highcorrWTvsFmr1(4,1)=length(find(MatAll_corrected2.(grouptemp).(loom{10}).Mat));
highcorrWTvsFmr1(4,2)=length(find(MatAll_corrected2.(groupnames{group-1,:}).(loom{10}).Mat));
highcorrWTvsFmr1(5,1)=length(find(MatAll_corrected2.(grouptemp).(loom{11}).Mat));
highcorrWTvsFmr1(5,2)=length(find(MatAll_corrected2.(groupnames{group-1,:}).(loom{11}).Mat));


figure;set(gcf, 'Position',  [100, 100, 500, 300]);
bar(highcorrWTvsFmr1);breakyaxis([1000 3500]) %%% althought it looks weird, kind of bended. I might just pass the values to Prism

%%%% i broke the axis with
%%%% MikeCF (2020). Break Y Axis (https://www.mathworks.com/matlabcentral/fileexchange/45760-break-y-axis), MATLAB Central File Exchange. Retrieved February 15, 2020.

%saveas('bargraph_WTminusFmr1_edges.svg')

%% with cross-validation values

%%%% I will quickly try to get the cross-validation matrices to be able to
%%%% do stats on this. 


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


for group=1:3

names = fieldnames(Data_corrMat4.(groupnames{group}));

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected_crossval=[];

for fish=1:length(names)-1

    moment=[1 2 3 4 5 6 11 12]; %%% looms 1,2,3,4,5,10 and 11 (cause 1 is pre loom)
    %loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(squeeze(squeeze(cross_val_fmr1.(groupnames{group}).CorrMatrices_mean2(moment(m),fish,:,:)))),0.75);
     MatAll_corrected_crossval(m,fish,:,:)=Mat;
     end

end
cross_val_fmr1.(groupnames{group}).MatAll_corrected_crossval=MatAll_corrected_crossval;

end

%%


%%%% getting the number of high corr >0.75

for group=1:3

names = fieldnames(Data_corrMat4.(groupnames{group}));

crossval_highcorr=[];
for fish=1:length(names)-1
for m=1:length(moment)

temp_highcorr=length(find(squeeze(squeeze(cross_val_fmr1.(groupnames{group}).MatAll_corrected_crossval(m,fish,:,:)))));

crossval_highcorr(m,fish)=temp_highcorr;
end
end
cross_val_fmr1.(groupnames{group}).crossval_highcorrNum=crossval_highcorr;
end


