
%%%% this script is to collect the idxs of more brain regions, in
%%%% particular the hindbrain regions to try to find the area that switches
%%%% from slopehab to nonhab depending on the ISI. 


load('Zbrain_Masks.mat');
BrainReg_F20=load('BrainReg_F20.mat');
BrainReg_F60=load('BrainReg_F60.mat');
BrainReg_S20=load('BrainReg_S20.mat');
BrainReg_S60=load('BrainReg_S60.mat');

ROI_temp2.f20=BrainReg_F20.ROI_temp2;
ROI_temp2.f60=BrainReg_F60.ROI_temp2;
ROI_temp2.s20=BrainReg_S20.ROI_temp2;
ROI_temp2.s60=BrainReg_S60.ROI_temp2;

datasets=['f20'; 'f60'; 's20'; 's60'];

PerBrainRegions=struct();

RegionList={};

for i=find(contains(Zbrain_Masks(:,1),'Rhombencephalon'))

RegionList=(Zbrain_Masks(i,2));

end

RegionList(146)=[]; %%% this is the one with the whole Rhombencephalon

%RegionList={'Telencephalon','Pallium','Subpallium','Caudal Hypothalamus','Dorsal Thalamus','Ventral Thalamus','Eminentia Thalami','(M1)','(M2)','AF4','AF6','AF8','Thalamus','Cerebellum','Semicircularis','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain','Msystem'};

%RegionList2={'Telencephalon','Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain','Msystem','coeruleus'};





%%

%%% this is to look at brain regions

for data=1:4
progressbar;
for i=1:length(RegionList)
    progressbar(i/length(RegionList));
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
    IsInBrainRegion=ismember(ROI_temp2.(datasets(data,:)),Mask,'rows');
    if strcmp(regionName,'(M1)')
        regionName='M1';
    elseif strcmp(regionName,'(M2)')
        regionName='M2';
    end
    PerBrainRegions.(datasets(data,:)).(strcat('reg',num2str(i))).idx=find(IsInBrainRegion==1);    
end

end


save('All_More_BrainReg3.mat','PerBrainRegions','RegionList','ROI_temp2','-v7.3');


