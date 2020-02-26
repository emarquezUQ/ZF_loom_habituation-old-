
load('Max_response_perBrain_all_CL4_corrected.mat');


for brain=1:length(RegionList)

figure('Position',[0 100 1200 300]);
counter=1;
for clust=[2 3 1]
subplot(1,3,counter);

if clust==2
    c=[0 1 0];
elseif clust==3
    c=[0 0 1];
else
    c=[1 0 0];
end

if ismember(clustersF{clust},fieldnames(Max_resp_f20_perfishNbrain2.(RegionList{brain})))
plot(nanmean(Max_resp_f20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})),'Color',c*0.75);
else
end
hold on;
if ismember(clustersF{clust},fieldnames(Max_resp_f60_perfishNbrain2.(RegionList{brain})))
plot(nanmean(Max_resp_f60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})),'Color',c*0.75,'LineStyle',':');
else
end
hold on;
if ismember(clustersF{clust},fieldnames(Max_resp_s20_perfishNbrain2.(RegionList{brain})))
plot(nanmean(Max_resp_s20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})),'Color',c);
else
end
hold on;
if ismember(clustersF{clust},fieldnames(Max_resp_s60_perfishNbrain2.(RegionList{brain})))
plot(nanmean(Max_resp_s60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})),'Color',c,'LineStyle',':');
else
end

xlabel('Stimuli');ylabel('Normalized response');%title(clustersF{clust});

v = get(gca,'Position');
set(gca,'Position',[v(1) v(2)*1.5 v(3:4)])

counter=counter+1;

end
suptitle(RegionList{brain});


saveas(gcf,strcat('MaxPerFish_',RegionList{brain},'.svg'));
end

