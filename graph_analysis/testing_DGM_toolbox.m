
%%%%% this script is to play around with the toolbox from Ann Sizemore (Dynamic Graph Metrics)
%%%%% i download it from github: https://github.com/asizemore/Dynamic-Graph-Metrics
% Reference: Ann E. Sizemore and Danielle S. Bassett, "Dynamic Graph 
% Metrics: Tutorial, Toolbox, and Tale." Submitted. (2017)

test3DMat=[];

for k=1:21

    test3DMat=cat(3,test3DMat,Data_corrMat_Null2.loomsR{1,k});

end

%help arrayToContactSeq

test3DMat=arrayToContactSeq(test3DMat,0);

%timeMat=timeAggregate_bin(test3DMat);



%%
[tadj,threshold]=thresholdMatDensity(Data_corrMat_Null2.loomsR{1,2},0.05);
%%%% note: this function doesnt take into account absolute values as it
%%%% is... so negative correlations are lost. 

threshold

%%
Nc=[0 0 0];

Ec=[0 0 1];



figure;plotArcNetwork(tadj,Nc,Ec);



%% trying to put it together...

   
    
 test3DMat=[];

for k=21:21 %%% i am starting at the 3rd loom or other looms to test

    Mat = threshold_absolute(abs(Data_corrMat2.f20.Mean_corrMat{1,k}(keep,keep)),0.75); 
    test3DMat=cat(3,test3DMat,Mat);

end   

test3DMat=arrayToContactSeq(test3DMat,0);

timeMat=timeAggregate_bin(test3DMat);


Nc=[0 0 0];

Ec=[0 0 1];

figure;plotArcNetwork(timeMat,Nc,Ec);
