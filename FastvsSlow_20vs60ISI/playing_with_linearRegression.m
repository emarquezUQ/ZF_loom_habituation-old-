
figure;
for i=1:8
subplot(3,3,i);
plot(rawregress_CN(i,:));
end


figure;
scatter(rawregress_CN(2,:),rawregress_CN(3,:)); lsline;

corr1=corrcoef(rawregress_CN(2,:),rawregress_CN(3,:));

figure;
scatter(rawregress_CN(2,:),rawregress_CN(5,:));lsline;

corr2=corrcoef(rawregress_CN(2,:),rawregress_CN(5,:));


figure;
scatter(rawregress_CN(2,:),rawregress_CN(7,:));lsline;

corr3=corrcoef(rawregress_CN(2,:),rawregress_CN(7,:));

figure;
for i=1:8
subplot(3,3,i);
histogram(rawregress_CN(i,:));
end

normality=[];
for i=1:8
normality(i) = kstest(rawregress_CN(i,:));
end


figure;
for i=1:8
subplot(3,3,i);
scatter(rawregress_CN(2,:),rawregress_CN(i,:));lsline;
end

corr_all=[];
for i=1:8
    corr_temp=corrcoef(rawregress_CN(2,:),rawregress_CN(i,:));
    
    
    corr_all(i)=corr_temp(1,2);
end

corr_all2= corr_all.^2; 



%%
ModelResults=[];
  for i=1:8
    mdl=fitlm(rawregress_CN(2,:)',rawregress_CN(i,:),'linear');
    
    
    ModelResults(i).coef=mdl.Coefficients;
    ModelResults(i).MSE=mdl.MSE;
    ModelResults(i).Fitted=mdl.Fitted;
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
    ModelResults(i).Residuals=mdl.Residuals.Raw;
  end

  figure;
   for i=1:8
    subplot(3,3,i);
    histogram(ModelResults(i).Residuals);
   end
  
   
   normality_res=[];
   for i=1:8
    normality_res(i) = kstest(ModelResults(i).Residuals);
   end
   