% %% Adding the correct toolboxes 
% pharm_path(0)
% %% Load the fMRI and plot a slices
% 
% [F_FileName,F_PathName] = uigetfile('*', 'Select the fMRI');
% fMRI_Orig = load_untouch_nii(strcat([F_PathName F_FileName]))
% % fMRI = load_untouch_nii('filtered_func_data');
% 
% [FileName,PathName] = uigetfile('*', 'Select the mask');
% Mask = load_untouch_nii(strcat([PathName FileName]))
% 
% figure
% Z_Slice = fMRI_Orig.img(:,:,19,20) ;
% mesh(Z_Slice); view([0 80])
% xlabel('X')
% ylabel('Y')


%% Converting the 4-D matirx to a 2-D matrix
cfg = []
cfg.direction = 42 ;   % 4D to 2D conversion
fMRI = FourDTwoDConvert(cfg, fMRI3, Mask) 

%% Converting the 2-D matirx to a 4-D matrix
% cfg = []
% cfg.direction = 24 ;   % 2D to 4D conversion
% fMRI_saved = FourDTwoDConvert(cfg, fMRI, Mask) 

%% Dummy fMRI for fieldtrip
clear fMRI_data ;
fMRI_data.trial{1} = fMRI.Time_Series ; 
fMRI_data.time{1} = 0:2.2:(size(fMRI.Time_Series, 2))*2.2 -2.2 ; 
fMRI_data.fsample = 1/2.2 ; 

for I = 1:size(fMRI.Time_Series,1)
    fMRI_data.label{I} = strcat('vox',num2str(I));
end

%%
cfg_filt = [] ;
cfg_filt.bpfilter = 'yes';
cfg = [] ; 
cfg.direction = 24 ;
fMRI_Bands = [ 0.01  0.02 ] ; %Define bands. For example NBack is 1/60 = 0.016

% cfg_filt.bpfreq = [0.001 .0033] ;
% [fMRI_data_Filtered] = ft_preprocessing(cfg_filt, fMRI_data) ;
% fMRI.Time_Series = fMRI_data_Filtered.trial{1}(:,:) ;
% fMRI_ppp1 = FourDTwoDConvert(cfg, fMRI, Mask) 

cfg_filt.bpfreq = fMRI_Bands(1,:) ;
[fMRI_data_Filtered] = ft_preprocessing(cfg_filt, fMRI_data) ;
fMRI.Time_Series = fMRI_data_Filtered.trial{1}(:,:) ;
fMRI_1 = FourDTwoDConvert(cfg, fMRI, Mask) 

if   size(fMRI_Bands,1) > 1
    cfg_filt.bpfreq = fMRI_Bands(2,:) ;
    [fMRI_data_Filtered] = ft_preprocessing(cfg_filt, fMRI_data) ;
    fMRI.Time_Series = fMRI_data_Filtered.trial{1}(:,:)  ;
    fMRI_2 = FourDTwoDConvert(cfg, fMRI, Mask) 
end

if   size(fMRI_Bands,1) > 2 
    cfg_filt.bpfreq = fMRI_Bands(3,:) ;
    [fMRI_data_Filtered] = ft_preprocessing(cfg_filt, fMRI_data) ;
    fMRI.Time_Series = fMRI_data_Filtered.trial{1}(:,:)  ;
    fMRI_3 = FourDTwoDConvert(cfg, fMRI, Mask) 
end

if   size(fMRI_Bands,1) > 3 
    cfg_filt.bpfreq = fMRI_Bands(4,:) ;
    [fMRI_data_Filtered] = ft_preprocessing(cfg_filt, fMRI_data) ;
    fMRI.Time_Series = fMRI_data_Filtered.trial{1}(:,:)  ;
    fMRI_4 = FourDTwoDConvert(cfg, fMRI, Mask) 
end
%% To create multi-band image
fMRI_saved = fMRI ;
ZSize = size(fMRI.img,3) ; 
fMRI_saved.hdr.dime.dim(1,4) = ZSize*size(fMRI_Bands,1) ;

Mask2 = Mask ;
ZSize = size(Mask.img,3) ; 
Mask2.hdr.dime.dim(1,4) = ZSize*size(fMRI_Bands,1) ;

if     size(fMRI_Bands,1) == 1
    fMRI_saved.img(:,:, 1          :   ZSize, :)   = fMRI_1.img ; 
    
elseif size(fMRI_Bands,1) == 2
    fMRI_saved.img(:,:, 1          :   ZSize, :)   = fMRI_1.img ; 
    fMRI_saved.img(:,:, 1 +  ZSize : 2*ZSize, :)   = fMRI_2.img ;
    
    Mask2.img(:,:,ZSize   + 1:ZSize*2, :) = Mask.img ;
elseif size(fMRI_Bands,1) == 3
    fMRI_saved.img(:,:, 1          :   ZSize, :)   = fMRI_1.img ; 
    fMRI_saved.img(:,:, 1 +  ZSize : 2*ZSize, :)   = fMRI_2.img ;
    fMRI_saved.img(:,:, 1 +2*ZSize : 3*ZSize, :)   = fMRI_3.img ;

    Mask2.img(:,:,ZSize   + 1:ZSize*2, :) = Mask.img ;
    Mask2.img(:,:,ZSize*2 + 1:ZSize*3, :) = Mask.img ;
elseif size(fMRI_Bands,1) == 4
    fMRI_saved.img(:,:, 1          :   ZSize, :)   = fMRI_1.img ; 
    fMRI_saved.img(:,:, 1 +  ZSize : 2*ZSize, :)   = fMRI_2.img ;
    fMRI_saved.img(:,:, 1 +2*ZSize : 3*ZSize, :)   = fMRI_3.img ;
    fMRI_saved.img(:,:, 1 +3*ZSize : 4*ZSize, :)   = fMRI_4.img ;
    
    Mask2.img(:,:,ZSize   + 1:ZSize*2, :) = Mask.img ;
    Mask2.img(:,:,ZSize*2 + 1:ZSize*3, :) = Mask.img ;
    Mask2.img(:,:,ZSize*3 + 1:ZSize*4, :) = Mask.img ; 
end

%%
save_untouch_nii(fMRI_saved,  'SimulatedSource_028_MultiBand') % Bands are concatenated in the Z dimension
%save_untouch_nii(Mask2,  'Multi_Band_Mask')

