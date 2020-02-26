%%% this script is to make 3d figures with the parts of the brain of Zbrain
%%% atlas and to plot on top of them the ROIs of the experiments. 

load('Zbrain_Masks.mat');

%%% this is to make a matrix with the location of the main structures 
Zbrain_brainMask=vertcat(Zbrain_Masks{[76 113 259 274 294],3}); %%% 78 is the eyes
 
Zbrain_brainMask=unique(Zbrain_brainMask,'rows');

figure;scatter(Zbrain_brainMask(:,1),Zbrain_brainMask(:,2),'.');
figure;scatter3(Zbrain_brainMask(:,1),Zbrain_brainMask(:,2),Zbrain_brainMask(:,3),'.');

%%%% this is to find the points that are in the boundarys of the brain.
k=boundary(Zbrain_brainMask,1); %%% i am using a shrink factor (0-1) of 1.  

%%% this is to make a Patch 3-D polygon with the boundarys
h=trisurf(k,Zbrain_brainMask(:,1),Zbrain_brainMask(:,2),Zbrain_brainMask(:,3),'FaceColor','red','FaceAlpha',0.1);

%%% to reduce the faces of the patch while keeping the same shape. I am
%%% doing it to 20%
h3=reducepatch(h,.2,'verbose');
brain3D=h3;

%%% to make a Brain Patch of the reduced patch from before
patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
view(3);

%%% and to test if i can plot things inside it. 
hold on;
scatter3(Zbrain_Masks{282,3}(:,1),Zbrain_Masks{282,3}(:,2),Zbrain_Masks{282,3}(:,3),'filled');


%%% its working!!

%%

%%% now I will try to make shapes of structures inside the brain. 

load('All_More_BrainReg2.mat');


BrainRegions3D=struct();

%progressbar;
for i=1:length(RegionList)
    %progressbar(i/length(RegionList));
    regionName=RegionList{i}%;
    if strcmp(regionName,'Telencephalon')
        Mask=Zbrain_Masks{294,3};
    elseif strcmp(regionName,'Hindbrain')
        Hindbrain_Mask=Zbrain_Masks{259,3};
        Mask=Zbrain_Masks{131,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove cerebellum
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Zbrain_Masks{295,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove MON
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Hindbrain_Mask;
    elseif strcmp(regionName,'Msystem')
        Mask=[];
        Mask_temp=[];
        Msystem_masks=[184 186 187];
        for j=1:3
           %Mask_temp=Zbrain_Masks{Msystem_masks(j),3};
           Mask=vertcat(Mask,Zbrain_Masks{Msystem_masks(j),3});
        end
    
        clear Mask_temp
    else
        Mask=[];
        IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
        IndexC=find(not(cellfun('isempty', IndexC)));
        for j=IndexC
            if isempty(Mask)
                Mask=Zbrain_Masks{j,3};
            else
                Mask=vertcat(Mask,Zbrain_Masks{j,3});
            end
        end
    end
    Mask=unique(Mask,'rows');
    
    %%%% this is to find the points that are in the boundarys of the brain.
    k_temp=boundary(Mask,1); %%% i am using a shrink factor (0-1) of 1.  

    %%% this is to make a Patch 3-D polygon with the boundarys
    h_temp=trisurf(k_temp,Mask(:,1),Mask(:,2),Mask(:,3),'FaceColor','red','FaceAlpha',0.1);

    %%% to reduce the faces of the patch while keeping the same shape. I am
    %%% doing it to 20%
    h2_temp=reducepatch(h_temp,.2,'verbose');
    
    
    if strcmp(regionName,'(M1)')
        regionName='M1';
    elseif strcmp(regionName,'(M2)')
        regionName='M2';
    end
    BrainRegions3D.(regionName(~isspace(regionName))).poly=h2_temp;    
end


%%% to test if it works
patch(h3,'EdgeColor','none','FaceAlpha',0.1);
hold on;
patch(BrainRegions3D.Pallium.poly,'EdgeColor','none','FaceColor','red','FaceAlpha',0.1);

%%% it works!!

save('zbrain3D.mat','brain3D','BrainRegions3D');


%%
%%% now, to put selected areas with different colors. 

RegionList_short={'Pallium','Dorsal Thalamus','Ventral Thalamus','(M2)','Cerebellum','Tectum','Tegmentum','Habenula','Pretectum','Msystem'};

colors = distinguishable_colors(length(RegionList_short),[1 1 1; 0 0 0]);
colors = colors*256;

patch(brain3D,'EdgeColor','none','FaceAlpha',0.025);
for i=1:length(RegionList_short)
  hold on;  
  regionName=RegionList_short{i};
  if strcmp(regionName,'(M1)')
        regionName='M1';
    elseif strcmp(regionName,'(M2)')
        regionName='M2';
    end
  patch(BrainRegions3D.(regionName(~isspace(regionName))).poly,'EdgeColor','none','FaceColor',colors(i,:)/256,'FaceAlpha',0.05);
 legend(horzcat('Brain',RegionList_short));
end  

%%% and to add the ROIs
%load('f20_cleaned_idxs.mat') %% for f20

colors2=['r','g','b','m'];

for i=1:4
hold on;
temp_fieldname=fieldnames(clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(i);
scatter3(ROI_temp2.f20(clust_f20_CL4_cleaned.(char(temp_fieldname)),1),ROI_temp2.f20(clust_f20_CL4_cleaned.(char(temp_fieldname)),2),ROI_temp2.f20(clust_f20_CL4_cleaned.(char(temp_fieldname)),3),20,'filled',colors2(i),'MarkerFaceAlpha',0.2);
end


%%
%%% trying with specific structures. Its very cool!!!!
RegionList_single={'Subpallium'};
regionName=RegionList_single{1};
  if strcmp(regionName,'(M1)')
        regionName='M1';
    elseif strcmp(regionName,'(M2)')
        regionName='M2';
    end
  patch(BrainRegions3D.(regionName(~isspace(regionName))).poly,'EdgeColor','none','FaceColor',colors(1,:)/256,'FaceAlpha',0.05);
 legend(horzcat(RegionList_single));
for i=1:4
hold on;
temp_fieldname=fieldnames(clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(i);
scatter3(ROI_temp2.f20(intersect(PerBrainRegions.f20.(regionName(~isspace(regionName))).idx,clust_f20_CL4_cleaned.(char(temp_fieldname))),1),ROI_temp2.f20(intersect(PerBrainRegions.f20.(regionName(~isspace(regionName))).idx,clust_f20_CL4_cleaned.(char(temp_fieldname))),2),ROI_temp2.f20(intersect(PerBrainRegions.f20.(regionName(~isspace(regionName))).idx,clust_f20_CL4_cleaned.(char(temp_fieldname))),3),20,'filled',colors2(i),'MarkerFaceAlpha',1);
end

%%
%%% multiple specific structures

%%
%%% now, to put selected areas with different colors. 

RegionList_short={'Dorsal Thalamus','Ventral Thalamus','(M2)','Tegmentum','Pretectum','Msystem'};

colors = distinguishable_colors(length(RegionList_short),[1 1 1; 0 0 0]);
colors = colors*256;

patch(brain3D,'EdgeColor','none','FaceAlpha',0.025);
for i=1:length(RegionList_short)
  hold on;  
  regionName=RegionList_short{i};
  if strcmp(regionName,'(M1)')
        regionName='M1';
    elseif strcmp(regionName,'(M2)')
        regionName='M2';
    end
  patch(BrainRegions3D.(regionName(~isspace(regionName))).poly,'EdgeColor','none','FaceColor',colors(i,:)/256,'FaceAlpha',0.2);
 legend(horzcat('Brain',RegionList_short));
end  

%%% and to add the ROIs
%load('f20_cleaned_idxs.mat') %% for f20

colors2=['r','g','b','m'];

for i=1:length(RegionList_short)
    hold on;  
  regionName=RegionList_short{i};
  if strcmp(regionName,'(M1)')
        regionName='M1';
    elseif strcmp(regionName,'(M2)')
        regionName='M2';
    end
for j=2%:4
hold on;
temp_fieldname=fieldnames(clust_f20_CL4_cleaned);temp_fieldname=temp_fieldname(j);
scatter3(ROI_temp2.f20(intersect(PerBrainRegions.f20.(regionName(~isspace(regionName))).idx,clust_f20_CL4_cleaned.(char(temp_fieldname))),1),ROI_temp2.f20(intersect(PerBrainRegions.f20.(regionName(~isspace(regionName))).idx,clust_f20_CL4_cleaned.(char(temp_fieldname))),2),ROI_temp2.f20(intersect(PerBrainRegions.f20.(regionName(~isspace(regionName))).idx,clust_f20_CL4_cleaned.(char(temp_fieldname))),3),20,'filled',colors2(j),'MarkerFaceAlpha',1);
end
end
