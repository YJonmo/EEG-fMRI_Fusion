

%% Adding the correct toolboxes 
% pharm_path(0)
% CustomPath_SourceReconst()
%% EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG EEG 
    %% Loading the Source space and sensor space data     
    [FileName,PathName] = uigetfile('*.mat', 'Select the source space EEG data (e.g., KET1)');
    load(strcat([PathName FileName]))
%     [FileName,PathName] = uigetfile('*.mat', 'Select the source space EEG data (e.g., KET1)');
%     load(strcat([PathName FileName]))
    %% Resample the previously cleaned EEG 
    Data = data_reconst ;
    cfg2 = [] ;
    cfg2.lpfilter = 'yes' ; 
    cfg2.lpfreq = [20] ;
    [Data] = ft_preprocessing(cfg2, Data) ;
    cfg2 = [];
    cfg2.resamplefs = 75;
    data_short_resamp = ft_resampledata(cfg2, Data)
    if isfield(Data, 'sampleinfo')
        sampleinfo = ceil(Data.sampleinfo/(Data.fsample/cfg2.resamplefs));
        data_short_resamp.sampleinfo = sampleinfo;
    end
    
    %% Rereferencing
    ft_progress('init', 'text', 'Rereferencing the EEG data ...');
    cfg2 = [];
    cfg2.implicitref   = 'FCz' ;
    cfg2.reref = 'yes';
    cfg2.refchannel = 'all';
    [data_short_resamp] = ft_preprocessing(cfg2, data_short_resamp);

    %% Loading the forward model
   
    [FileName,PathName] = uigetfile('*.mat', 'Select a forward model >><< Or press cancel to use the BEM template model');
    if FileName == 0
        cfg = [] ; 
        %ForwardModel = CreatT1BEM(SensorData);
        load BrainCap1010
        cfg.elec = elec;
        ForwardModel = Create_TemplateBEM(cfg, data_short_resamp);
        ForwardModel_Orig = ForwardModel;
    else
        load(strcat([PathName FileName]))
        ForwardModel_Orig = ForwardModel ; 
    end
    
%% Calling the source-space ICA_V3 or source-space ICA_Simplified (version 3).	This version performs band pass filtring on the source space data
    % Match the forward model structure to the EEG data
    ForwardModel = MatchForward2Data(data_short_resamp, ForwardModel)

    cfg = [] ; 
    cfg.bandpass = [  1 3.8; 4 8; 8.2 12;  12.2 22; 23 74]
    % cfg.reducerankBy = [5; 6; 7; 8] ; 
    % cfg.rankThresh = [20 ; 20 ; 20; 20; 20] ;
    % cfg.rankThreshTotl = [20] ;
    % cfg.rankKaiser = [] ; 
    cfg.bandRanks = [24;]% 8; 8]% 12; ] ; 
    % cfg.bandpass = [ ] ;
    %cfg.reref = 'FCz' ; 
    %cfg.trials = 1:50 ;
    cfg.headmodel = ForwardModel.headmodel;
    cfg.grid =  ForwardModel.grid;
    cfg.channel = ForwardModel.elec.label; 
    cfg.elec = ForwardModel.elec;
    %cfg.resampleFs = 200;
    %cfg.method = 'sobi' ; % Other methods 'runica', 'fastica', 'binica', 'pca', 'svd', 'jader', 'varimax', 'dss', 'cca', 'sobi', 'white'
    [SourceSpaceStuff] = Source_Space_ICA_V3(cfg, data_short_resamp );
    %[SourceSpaceStuff] = Source_Space_ICA_Simplified(cfg, SensorData);

    No_Vox = size(SourceSpaceStuff.SpatialICs_Maps(:,1),1) ;

    %% Check the components first

%     %% Choose the brain signals and reduce the size of the source space components
%     cfg = [] ;
%     %GoodComponents = [] ;
%     cfg.channel = GoodComponents ;
%     cfg.channel = [  8 ] ;
%     SourceSpaceStuff2 = SourceSpaceStuff ; 
%     SourceSpaceStuff2.TemporalPCs = ft_selectdata(cfg, SourceSpaceStuff2.TemporalPCs);
%     SourceSpaceStuff2.TemporalICs = ft_selectdata(cfg, SourceSpaceStuff2.TemporalICs);
%     SourceSpaceStuff2.SpatialICs = SourceSpaceStuff2.SpatialICs(:, cfg.channel) ;
%     SourceSpaceStuff2.SpatialICs_Maps = SourceSpaceStuff2.SpatialICs_Maps(:, cfg.channel) ;
%     SourceSpaceStuff2.SpatialPCs_Maps = SourceSpaceStuff2.SpatialPCs_Maps(:, cfg.channel) ;
%     SourceSpaceStuff2.MixingMatrix = SourceSpaceStuff2.MixingMatrix(:, cfg.channel) ;   
    
    %% Or choose all the components
        SourceSpaceStuff2 =SourceSpaceStuff ;
    %% Concatenate the trials (from 2.2s trial length trials back to their largest possible trial lengths)
    cfg = [];
    data_sections = ft_combinetrials(cfg, SourceSpaceStuff2.TemporalICs);
    %% Asign the spatial maps 
    for Comp_Index = 1:size(SourceSpaceStuff2.SpatialICs_Maps,2)
        data_sections.Spatial(:,Comp_Index) = SourceSpaceStuff2.SpatialICs_Maps(:,Comp_Index)/max(abs(SourceSpaceStuff2.SpatialICs_Maps(:,Comp_Index))) ;
    end
%     %% Reduce data to just the comps we want
%     cfg = [];
%     cfg.channel = 'all';
%     %cfg.channel = GoodComponents ;
%     data_ICcont = ft_selectdata(cfg, data_sections);
    
    
    
    %% Calculate the power envelope
    filtabs = data_sections;
    for i = 1 : length(filtabs.trial)
        filtabs.trial{i} = abs(filtabs.trial{i}); %Envelope
    end
    
    %% Generation of standard HRF/convolution (run in FSL)

    %Create hrf at the correct sampling interval - 0 input effectively removes the second gamma function from this
    %[hrf] = doubleGammaHrf(1/filtabs.fsample,[4 16],[1 .5],.2,20);%,tp,beta,rt,len)
    %[hrf] = doubleGammaHrf(1/filtabs.fsample,[4 16],[1 1],.0,20);%,tp,beta,rt,len)
    [hrf] = doubleGammaHrf(1/filtabs.fsample,[6 16],[1 1],.0,20);%,tp,beta,rt,len)
    hrf_time =( 0 :1/filtabs.fsample:20)';
    hrf_time(end) = [];
    figure; 
    plot(hrf_time, hrf);
    title('Standard HRF');
    %% 
%     %Now convolve power envelope with the hrf
%     filthrf = filtabs;
%     filthrf = rmfield(filthrf, 'trial') ;
%     i = 1 ;
%     for i = 1 : length(filtabs.trial)
%         filthrf.trial{i} = zeros(size(filtabs.trial{i},1), size(conv(filtabs.trial{i}(1,:),hrf),2)) ; 
%         for Comp_Index = 1:size(filtabs.trial{1},1)
%             filthrf.trial{i}(Comp_Index,:) = conv(filtabs.trial{i}(Comp_Index,:),hrf) *1; %Envelope
%         end
%     end
%     %plot power envelop plus convolved version
%     figure
%     hold on
%     plot(real(filtabs.trial{1}(1,:)), 'g')
%     plot(real(filthrf.trial{1}(1,:)), 'r') 

    %% 
    %NTRs = 437;
    NTRs = 243;  % Here 243 was for NBack

    %This demonstrates warm up effects in the convolution and suggests that restarting the convolution for every bad trial means losing a lot of data to warm-up effects.
    %function "nanconv" on MATLAB website might be an option on the entire time series. Here I have just used linear interpolation to 
    %glue trials together ..... looks ugly but things will get downsampled and convolved later......

    %Approach - Glue all the sections together, interpolate (bad trials)
    MaxSampleIndex = filtabs.sampleinfo(end,2);
    Tau = filtabs.time{1}(1,2) - filtabs.time{1}(1,1 ) 
    timee  = [1:MaxSampleIndex ] ;
    timee = timee*Tau ;
    SamplesPerTR = MaxSampleIndex./ NTRs;
    rawvec = NaN([1 MaxSampleIndex]);
    %%
    filthrf2 = filtabs;
    filthrf2 = rmfield(filthrf2, 'trial') ;
    filthrf2 = rmfield(filthrf2, 'time') ;
    filthrf2 = rmfield(filthrf2, 'sampleinfo') ;
    i = 1 ;
    %%
    filthrf2.trial{1} = zeros(size(filtabs.trial{1},1), MaxSampleIndex) ; 
    for Comp_Index = 1:size(filtabs.trial{1},1) 
        if length(filtabs.trial) > 1
            for i = 1: length(filtabs.trial) - 1
                rawvec(filtabs.sampleinfo(i,1):filtabs.sampleinfo(i,2)) = filtabs.trial{i}(Comp_Index,:);
                Means{i}(1,2) = mean(rawvec(filtabs.sampleinfo(i,2) - 2.0*filtabs.fsample : filtabs.sampleinfo(i,2))) ; 
                Means{i}(1,1) = mean(rawvec(filtabs.sampleinfo(i,1): filtabs.sampleinfo(i,1) + 2.0*filtabs.fsample)) ; 
                Vars{i}(1,2)  = mean(rawvec(filtabs.sampleinfo(i,2) - 2.0*filtabs.fsample : filtabs.sampleinfo(i,2))) ; 
                Vars{i}(1,1)  = mean(rawvec(filtabs.sampleinfo(i,1): filtabs.sampleinfo(i,1) + 2.0*filtabs.fsample)) ; 
                Size_NAN = size(filtabs.sampleinfo(i,2)+1:filtabs.sampleinfo(i+1,1),2) ; 
                rawvec(filtabs.sampleinfo(i,2)+1:filtabs.sampleinfo(i+1,1)) = abs(randn(1,Size_NAN)*mean(Means{i})) ;%+ linspace(Means{i}(1,1), Means{i}(1,2), Size_NAN) ;
            end
            i = i + 1
            rawvec(filtabs.sampleinfo(i,1):filtabs.sampleinfo(i,2)) = filtabs.trial{i}(Comp_Index,:);
            Means{i}(1,2) = mean(rawvec(filtabs.sampleinfo(i,2) - 2.0*filtabs.fsample : filtabs.sampleinfo(i,2))) ; 
            Means{i}(1,1) = mean(rawvec(filtabs.sampleinfo(i,1): filtabs.sampleinfo(i,1) + 2.0*filtabs.fsample)) ; 
            Vars{i}(1,2)  = mean(rawvec(filtabs.sampleinfo(i,2) - 2.0*filtabs.fsample : filtabs.sampleinfo(i,2))) ; 
            Vars{i}(1,1)  = mean(rawvec(filtabs.sampleinfo(i,1): filtabs.sampleinfo(i,1) + 2.0*filtabs.fsample)) ; 
            Size_NAN = size(filtabs.sampleinfo(i,2)+1:length(rawvec),2) ; 
            rawvec(filtabs.sampleinfo(i,2)+1:end) = abs(randn(1,Size_NAN)*mean(Means{i})) ;
        else
            i = 1 
            rawvec(filtabs.sampleinfo(i,1):filtabs.sampleinfo(i,2)) = filtabs.trial{i}(Comp_Index,:);
            Means{i}(1,2) = mean(rawvec(filtabs.sampleinfo(i,2) - 2.0*filtabs.fsample : filtabs.sampleinfo(i,2))) ; 
            Means{i}(1,1) = mean(rawvec(filtabs.sampleinfo(i,1): filtabs.sampleinfo(i,1) + 2.0*filtabs.fsample)) ; 
            Vars{i}(1,2)  = mean(rawvec(filtabs.sampleinfo(i,2) - 2.0*filtabs.fsample : filtabs.sampleinfo(i,2))) ; 
            Vars{i}(1,1)  = mean(rawvec(filtabs.sampleinfo(i,1): filtabs.sampleinfo(i,1) + 2.0*filtabs.fsample)) ; 
            Size_NAN = size(filtabs.sampleinfo(i,2)+1:length(rawvec),2) ; 
            rawvec(filtabs.sampleinfo(i,2)+1:end) = abs(randn(1,Size_NAN)*mean(Means{i})) ;
        end

        X = find(~isnan(rawvec));
        Xq = find(isnan(rawvec));
        V = rawvec(X);
        Vq = interp1(X,V,Xq, 'linear'); %rough but ok since will downsample
        interpvec = rawvec;
        interpvec(Xq) = Vq;
        filthrf2.trial{1}(Comp_Index,:) = interpvec ; 
    end
    %% 
    filthrf2.time{1} = timee ;
    filthrf2.sampleinfo(1) = 1 ;
    filthrf2.sampleinfo(2) = filtabs.sampleinfo(end,2) ;

    filthrf3 = trial2continuous(filthrf2) ; 

%     X = find(~isnan(rawvec));
%     Xq = find(isnan(rawvec));
%     V = rawvec(X);
%     Vq = interp1(X,V,Xq, 'linear'); %rough but ok since will downsample
% 
%     interpvec = rawvec;
%     interpvec(Xq) = Vq;
% 
% 
%     %Plot the interpolation result
%     figure
%     %inds = [Xq(1)  Xq(1)+ 10000]; %Plot a section of data around a missing trial
%     %plot(timee(inds(1):inds(2)),interpvec(inds(1):inds(2)), 'r')
%     plot(timee,interpvec,'g')
%     hold on
%     %plot(rawvec(inds(1):inds(2)))
%     plot(timee, rawvec)
% 
%     title('Quick linear interpolation of missing trial');

    %% Now convolve with the hrf
    postconv = conv(interpvec, hrf);

    figure
    hold on
    plot(timee, interpvec, 'b')
    %plot(timee, postconv(125*10:end-125*10),'r')
    timee2 = [Tau:Tau:10        timee+10         Tau+10+timee(end):Tau:20+timee(end) - Tau] ;
    plot(timee2, postconv,'r')

    
    i = 1 ;
    filtconv.trial{1} = zeros(size(filtabs.trial{i},1), length(postconv)) ; 
    for Comp_Index = 1:size(filtabs.trial{1},1)
        filtconv.trial{1}(Comp_Index,:) = conv(filthrf3.trial{1}(Comp_Index,:) , hrf);
    end
    filtconv.time{1} = timee2 ;
    filtconv.label = filthrf3.label ; 
    filtconv.fsample = filthrf3.fsample ; 
    filtconv.Spatial = filthrf3.Spatial ;    
    
    %% Reject the first  11 s and the last 22 s of the EEG
    filtconv2 = filtconv ;
    filtconv2.trial{1} = filtconv.trial{1}(:,filtconv.fsample*11  : end-filtconv.fsample*22) ; 
    filtconv2.time{1} = filtconv.time{1}(1,1:size(filtconv2.trial{1},2)) ; 

    %% Down sample the data
    % postconv2 = postconv;
    cfg = [];
    cfg.resamplefs = 0.454545454545455
    %sampleinfo = ceil(Data.sampleinfo/(Data.fsample/cfg.resamplefs));
    %trialinfo  =      Data.trialinfo;
    filtconv2 = ft_resampledata(cfg, filtconv2)
    %Data.trialinfo = trialinfo; 
    %Data.sampleinfo = sampleinfo;
    % 
    % %% Now lets downsample into TR's ready for GLM
    % postconv(MaxSampleIndex+1:end)=[]; %remove extra part of convolution 
    % lowreg = decimate(postconv,(SamplesPerTR));   
    % lowreg_nohrf = decimate(interpvec, (SamplesPerTR));

    %% Band pass the convolve data 
    cfg_filt = [] ;
    cfg_filt.bpfilter = 'yes';
    
    cfg = [] ; 
    cfg.direction = 24 ;
    EEG_data = filtconv2 ; 
    EEG_HRFed_Bands = [ 0.01  0.02 ] ; %THese should be the same as what were defined in fMRI_Multiband.m
    No_Freq_Bands = size(EEG_HRFed_Bands,1) ;  
    
    
    Size_temp = size(EEG_data.trial{1},1) ; 
    EEG_total = EEG_data ; 
    EEG_total.trial{1} = zeros(Size_temp*No_Freq_Bands, size(EEG_data.trial{1},2)) ;  
    
    
    cfg_filt.bpfreq = EEG_HRFed_Bands(1,:) ;
    [EEG_data_Filtered] = ft_preprocessing(cfg_filt, EEG_data) ;
    EEG1.Time_Series = EEG_data_Filtered.trial{1}(:,:) ;
    
    EEG_total.trial{1}(1 :Size_temp             , :)   = EEG1.Time_Series ;
    EEG_total.Spatial = [EEG_data.Spatial] ;
    
    if  size(EEG_HRFed_Bands,1) > 1
        cfg_filt.bpfreq = EEG_HRFed_Bands(2,:) ;
        [EEG_data_Filtered] = ft_preprocessing(cfg_filt, EEG_data) ;
        EEG2.Time_Series = EEG_data_Filtered.trial{1}(:,:)  ;
        EEG_total.trial{1}(Size_temp   + 1:Size_temp*2, :) = EEG2.Time_Series ;
        EEG_total.Spatial = [EEG_data.Spatial EEG_data.Spatial] ;
    end
    if  size(EEG_HRFed_Bands,1) > 2
        cfg_filt.bpfreq = EEG_HRFed_Bands(3,:) ;
        [EEG_data_Filtered] = ft_preprocessing(cfg_filt, EEG_data) ;
        EEG3.Time_Series = EEG_data_Filtered.trial{1}(:,:)  ;
        EEG_total.trial{1}(Size_temp*2 + 1:Size_temp*3, :) = EEG3.Time_Series ;
        EEG_total.Spatial = [EEG_data.Spatial EEG_data.Spatial EEG_data.Spatial] ;
    end
    if  size(EEG_HRFed_Bands,1) > 3
        cfg_filt.bpfreq = EEG_HRFed_Bands(4,:) ;
        [EEG_data_Filtered] = ft_preprocessing(cfg_filt, EEG_data) ;
        EEG4.Time_Series = EEG_data_Filtered.trial{1}(:,:)  ;
        EEG_total.trial{1}(Size_temp*3 + 1:Size_temp*4, :) = EEG4.Time_Series ;
        EEG_total.Spatial = [EEG_data.Spatial EEG_data.Spatial EEG_data.Spatial EEG_data.Spatial] ;
    end
  
    for I = 1:size(EEG_total.trial{1},1)
        EEG_total.label{I} = strcat('Sig',num2str(I));
    end
    
    %% Creating the time course of the paradigm ----____----____---- 
    paradigm = [] ;
    paradigm.time = filthrf3.time ; 
    paradigm.trial{1} = ones(No_Freq_Bands, size(filthrf3.time{1},2)) ; 
    paradigm.fsample = filthrf3.fsample ; 
    
    for I = 1:No_Freq_Bands
        paradigm.label{I} = strcat('Pardig',num2str(I));
    end
 
    for sample_Index = 1:16
        Start_temp = (sample_Index-1)*60*paradigm.fsample ; 
        paradigm.trial{1}(:, Start_temp + 1 : Start_temp + 30*paradigm.fsample ) = 0 ; 
    end
    
    paradigm2 = paradigm ; 
    paradigm2.trial{1} = zeros(size(paradigm.trial{1},1), length(postconv)) ;

    for I = 1:No_Freq_Bands
        paradigm2.trial{1}(I,:) = conv( paradigm.trial{1}(I,:), hrf);
    end

    paradigm2.trial{1} = paradigm2.trial{1}(:,paradigm2.fsample*11  : end-paradigm2.fsample*22) ; 
    paradigm2.time{1} = paradigm2.time{1}(1,1:size(paradigm2.trial{1},2)) ; 
        
    cfg = [];
    cfg.resamplefs = 0.454545454545455
    paradigm2 = ft_resampledata(cfg, paradigm2)
     
    paradigm3 = paradigm2 ;
    
    cfg_filt = [] ;
    cfg_filt.bpfilter = 'yes';
%     cfg_filt.bpfreq = [0.009 .052] ;
    cfg_filt.bpfreq = EEG_HRFed_Bands(1,:) ;
    
    [paradigm_filtered] = ft_preprocessing(cfg_filt, paradigm2) ;
    paradigm3.trial{1}(1,:) = paradigm_filtered.trial{1}(1,:) ; 
    
%      cfg_filt.bpfreq = [0.02 .05] ;
        cfg_filt.bpfreq = EEG_HRFed_Bands(2,:) ;
     [paradigm_filtered] = ft_preprocessing(cfg_filt, paradigm2) ;
     paradigm3.trial{1}(2,:) = paradigm_filtered.trial{1}(1,:) ;

%     cfg_filt.bpfreq = [0.06 .11] ;
%     [paradigm_filtered] = ft_preprocessing(cfg_filt, paradigm2) ;
%     paradigm3.trial{1}(3,:) = paradigm_filtered.trial{1}(1,:) ; 


%% fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI fMRI

    %% To produce the time-series of the ICA of the MELODIC
    % Load the ICA file and then load the original data 
    [FileName,PathName] = uigetfile('*', 'Select the fMRI_filtered_func');
    fMRI_Orig = load_untouch_nii(strcat([PathName FileName]))
    [FileName,PathName] = uigetfile('*', 'Select the melodic_IC');
    fMRI_ica_maps = load_untouch_nii(strcat([PathName FileName]))
    [FileName,PathName] = uigetfile('*', 'Select the Mask');
    Mask = load_untouch_nii(strcat([PathName FileName]))
    
    
    %% Converting the 4-D matirx to a 2-D matrix
    cfg = []
    cfg.direction = 42 ;
    fMRI = FourDTwoDConvert(cfg, fMRI_Orig, Mask) 
    fMRI_ica_maps2 = FourDTwoDConvert(cfg, fMRI_ica_maps, Mask) 
    
    %% Now multiply the data with the ICA weights to get the time-series of the components
    ICA_Time_Series_fMRI = (fMRI_ica_maps2.Time_Series'*fMRI.Time_Series) ; 

    
    %% Dummy fMRI for fieldtrip
    clear fMRI_data ;
    fMRI_data.trial{1} = ICA_Time_Series_fMRI ;
    % fMRI_data.trial{1} = fMRI.Time_Series ; 
    fMRI_data.time{1} = 0:2.2:(size(fMRI.Time_Series, 2))*2.2 -2.2 ; 
    fMRI_data.fsample = 1/2.2 ; 

    for I = 1:size(ICA_Time_Series_fMRI,1)
        fMRI_data.label{I} = strcat('vox',num2str(I));
    end
    

    %% Prepare the structure for the fusion  
    
%     cfg2 = [] ;
%     cfg2.bpfilter = 'yes' ; 
%     cfg2.bpfreq = [0.009 .11] ;
%     [fMRI_data] = ft_preprocessing(cfg2, fMRI_data) ;
%     
    
    fMRI1 = fMRI_data ;
    fMRI1.Spatial = fMRI_ica_maps2.Time_Series ; 
    
    
    
    %% fMRI Spatial map normalization
     
    for Comp_Index = 1:size(fMRI1.Spatial,2)
        fMRI1.Spatial(:,Comp_Index) = fMRI1.Spatial(:,Comp_Index)/max(abs(fMRI1.Spatial(:,Comp_Index))) ;
    end
    
    %% Rejection of the insignificant bands
    fMRI_temp = fMRI_ica_maps2 ; 
    fMRI_temp.Time_Series = fMRI1.Spatial ;
    cfg = []
    cfg.direction = 24 ;
    fMRI_temp = FourDTwoDConvert(cfg, fMRI_temp, Mask) 

    Zsize = size(fMRI_temp.img, 3)/No_Freq_Bands ;
    for Comp_Index = 1:size(fMRI1.Spatial,2)
        [row,col,v] = ind2sub(size(fMRI_temp.img(:,:,:,Comp_Index)), find(abs(fMRI_temp.img(:,:,:,Comp_Index)) == 1));
        if  v <= size(fMRI_temp.img, 3)/No_Freq_Bands                                    % First band
            fMRI_temp.img(:,:,Zsize+1:end, Comp_Index) = 0 ;
        elseif (v <= Zsize*2)                                                      % Second band
            fMRI_temp.img(:,:,[1:Zsize Zsize*2+1:end], Comp_Index) = 0 ;
        elseif (v <= Zsize*3)                                                      % Third band                        
            fMRI_temp.img(:,:, [1:Zsize*2], Comp_Index) = 0 ;
%           fMRI_temp.img(:,:, [1:Zsize*2 Zsize*3+1:end], Comp_Index) = 0 ;
        end
    end
    
    cfg = []
    cfg.direction = 42 ;
    fMRI_temp = FourDTwoDConvert(cfg, fMRI_temp, Mask) 
    %%
    fMRI1.Spatial = fMRI_temp.Time_Series ; 
    
    
    %% Reject the first  11 s of the fMRI

    fMRI1.trial{1} = fMRI1.trial{1}(:,fMRI1.fsample*11 + 1 : end) ; 
    fMRI1.time{1} =  fMRI1.time{1}(1,1:size(fMRI1.trial{1},2)) ;      
    
%% Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion Fusion 

    clear fMRI2 ;
    fMRI2 = fMRI1 ; 
    cfg2 = [];
    cfg2.resamplefs = fMRI2.fsample*10;
    
    fMRI2 = ft_resampledata(cfg2, fMRI2)    

    %EEG2 = filtconv2 ;
    %EEG2 = EEG_PrVPo2 ;
    
    EEG2 = EEG_total ; 
    EEG2.Spatial = EEG_total.Spatial ; 
%     EEG2.freq_analysis.powspctrm = EEG_total.freq_analysis.powspctrm ; 
    
%         EEG2 = Multi_ICA_EEG.TemporalICs ; 
%         EEG2.Spatial = Multi_ICA_EEG.SpatialICs ; 
%         EEG2.freq_analysis.powspctrm = Multi_ICA_EEG.freq_analysis.powspctrm ; 

    cfg2 = [];
    cfg2.resamplefs = EEG2.fsample*10;
    EEG2 = ft_resampledata(cfg2, EEG2)
    % EEG2.Spatial  = Multi_ICA_EEG.SpatialICs ; 

    cfg2 = [];
    cfg2.resamplefs = paradigm3.fsample*10;
    paradigm3_upsam = ft_resampledata(cfg2, paradigm3) ;    
    
    %%
    cfg_fus = [] ; 
    cfg_fus.method = 'sobi' ;    % Pref use SOBI or SVD !!  - Other methods 'runica', 'fastica', 'binica', 'pca', 'svd', 'jader', 'varimax', 'dss', 'cca', 'sobi', 'white'
    cfg_fus.sobi_delay = 20 ; %If it has beeen upsampled (e.g. by 10 then this 10 represents about 1 TR
    cfg_fus.rank_Total = [25] ;
    cfg_fus.fMRI{1} = fMRI2 ;
    cfg_fus.EEG{1}  = EEG2 ; 
    cfg_fus.paradigm = paradigm3_upsam ; 
    %cfg_fus.rank_auto = [] ; 
    FusedStuff = fuse_group(cfg_fus) ; 

    %% Plot the time courses
    cfg = [];
    Data = EEG2;
    cfg.viewmode = 'vertical';
    cfg.channelcolormap = [0 0 0];  %Above two lines plot in blue (events will be black)
    cfg.ploteventlabels = 'colorvalue';
    Limm = .0010000101
    cfg.ylim = [-Limm Limm]; %sets yscale
    cfg.blocksize = 541.2;
    cfg.continuous              = 'yes';
    ft_databrowser(cfg,Data);

    %% Freq analysis
    

    cfg = [];
    cfg.length  = 5;
    cfg.overlap = 0;
    data_cut = ft_redefinetrial(cfg, SourceSpaceStuff2.TemporalICs)

    
    cfg_FFT = [];
    cfg_FFT.method = 'mtmfft';
    cfg_FFT.taper = 'hanning';
    cfg_FFT.foilim = [.001 0.2];
    cfg_FFT.foilim = [1 22];
    cfg_FFT.pad = 'nextpow2' ;
    %cfg_FFT.pad = 400 ;
    freq_SSICs = ft_freqanalysis(cfg_FFT, data_cut);  

    %% Plot Source-space Spatial ICs via 3D scatter plot and the frequency of the component
    Curren_comp = 24 ; 
    % Frequency plot of the components

    FigHandle = figure('Position', [1000, 1400, 550, 170]);
    set(gcf, 'color', [1 1 1])
    plot(freq_SSICs.freq,freq_SSICs.powspctrm(Curren_comp,:))
    xlabel('Frequencies Hz')
    ylabel('Power')

    % Correct the frequency analysis
%     EEG2.freq_analysis.powspctrm = FusedStuff.Mixing_Spatial_EEG'*EEG2.freq_analysis.powspctrm ;

    % Plot EEG Source-space Spatial ICs via 3D scatter plot and the frequency of the component
%    Curren_comp = 15
%     figure ; plot(freq_analysis.freq, EEG2.freq_analysis.powspctrm(Curren_comp,:))
%     xlabel('Frequencies Hz')
%     ylabel('Power')
    %I = I +1 ;
    % Frequency plot of the components

    simple3Dmapping(EEG2.Spatial(:,Curren_comp), ForwardModel.grid.pos(ForwardModel.grid.inside,:))
    %simple3Dmapping(Multi_ICA_EEG.SpatialICs(:,Curren_comp), ForwardModel.grid.pos(ForwardModel.grid.inside,:))
    %simple3Dmapping(SourceSpaceStuff2.SpatialICs_Maps(:,Curren_comp), ForwardModel.grid.pos(ForwardModel.grid.inside,:))
    % simple3Dmapping(FusedStuff.SpatialICs_EEG(:,Curren_comp), ForwardModel.grid.pos(ForwardModel.grid.inside,:))

    set(gca,'XLim',[-90 90],'YLim',[-110 100],'ZLim',[-90 90])
    title(strcat(['Component ', num2str(Curren_comp)]))
    %% Group maps of the EEG tomos 
%     cfg = [] ;    %% Resample the previously cleaned EEG 
    Data = data_reconst ;
    cfg2 = [] ;
    cfg2.lpfilter = 'yes' ; 
    cfg2.lpfreq = [20] ;
    [Data] = ft_preprocessing(cfg2, Data) ;
    cfg2 = [];
    cfg2.resamplefs = 75;
    data_short_resamp = ft_resampledata(cfg2, Data)
    if isfield(Data, 'sampleinfo')
        sampleinfo = ceil(Data.sampleinfo/(Data.fsample/cfg2.resamplefs));
        data_short_resamp.sampleinfo = sampleinfo;
    end
    
    %% Rereferencing
    ft_progress('init', 'text', 'Rereferencing the EEG data ...');
    cfg2 = [];
    cfg2.implicitref   = 'FCz' ;
    cfg2.reref = 'yes';
    cfg2.refchannel = 'all';
    [data_short_resamp] = ft_preprocessing(cfg2, data_short_resamp);

%     cfg.Comps = [21] ; 
%     SourceICA_plotmaps(cfg, SourceSpaceStuff, ForwardModel)
    
    %% Map on the MRI image
    %Curren_comp = 14 ; 

    Inside_count = 0;
    source = ForwardModel.grid;
    for Vox_Index = 1 : size(ForwardModel.grid.inside,1)
        if  ForwardModel.grid.inside(Vox_Index) == 1
            Inside_count = Inside_count + 1;
            source.avg.Map(Vox_Index) = EEG2.Spatial(Inside_count,Curren_comp);
        else
            source.avg.Map(Vox_Index) = 0;
        end
    end

    cfg              = [];
    cfg.parameter    = 'Map';
    cfg.interpmethod = 'nearest';
    source_diff_int  = ft_sourceinterpolate(cfg, source, ForwardModel.mri);

    %And plot the results
    source_diff_int.coordsys = 'mni';
    cfg = [];
    cfg.funparameter    = 'Map';
    cfg.funcolorlim   = 'maxabs';
    %cfg.funcolorlim   = [-10 10];
    cfg.method = 'slice';
   % cfg.method = 'ortho';
    % cfg.atlas = [FIELDTRIPDIR '/template/atlas/aal/ROI_MNI_V4.nii'];
    ft_sourceplot(cfg,source_diff_int);
    
     %% Create the nii (nifti for FSL) from the 3D maps of the components
%     nii = make_nii(source_diff_int.Map) ;
%     save_nii(nii,strcat(['Comp_', num2str(Curren_comp),  '_PLA_NBk_018.nii']))
%     nii = make_nii(fMRI_ica_maps2.Temp_Map) ;
%     save_nii(nii,strcat(['Comp_', num2str(Curren_comp),  '_PLA_NBk_018.nii']))   
    %% Create the nii (nifti for FSL) from the 3D maps of the components
    cfg = [];
    cfg.filetype = 'nifti'
    cfg.datatype = 'float';
    cfg.parameter = 'Map'

    cfg.filename = strcat(['Comp_', '5Hz',  '_Simulated.nii']);
    ft_volumewrite(cfg, source_diff_int);
%     

%% Create the nii MRI of the individual
% cfg = [];
% cfg.filetype = 'nifti'
% cfg.datatype = 'float';
% cfg.parameter = 'anatomy'
% cfg.filename = strcat(['038',  '_PLA_NBk_MRI.nii']);
% ft_volumewrite(cfg, ForwardModel.mri);

% 
% cfg = [];
% cfg.filetype = 'nifti'
% cfg.datatype = 'float';
% cfg.parameter = 'Map_temp'
% cfg.filename = strcat(['018',  'fMRI_IC1.nii']);
% ft_volumewrite(cfg, fMRI_ica_maps2);
    %% Convert the 2D iamges of spatial maps into the 4D structure
    % These maps are the corresponding fMRI maps of fused data and are
    % different from Multi_ICA_fMRI.SpatialICs. Some of these maps may not
    % mean anything. Only by using the FusedStuff.SpatialICs you can
    % understand how much the fMRI has contributed to every fused data
    % and if contribution is high then the map is meaningful. Otherwise it 
    % would be random mixture of several Multi_ICA_fMRI.SpatialICs maps. 
    cfg = []
    cfg.direction = 24 ;
    %fMRI_saved = fMRI2 ; 
    %fMRI.Time_Series = Multi_ICA_fMRI.SpatialICs ;  
    % fMRI.Time_Series = fMRI1.Spatial ;  
    fMRI.Time_Series = FusedStuff.SpatialICs_fMRI ;  
    %fMRI.Time_Series = fMRI_ica_maps2.Time_Series*Multi_ICA_fMRI.SpatialICs ; 
    fMRI.Time_Series = (fMRI.Time_Series./(max(max(fMRI.Time_Series)))).*10 ;
    fMRI_saved = FourDTwoDConvert(cfg, fMRI, Mask) 
    fMRI_saved.hdr.dime.dim(1,5) = size(fMRI_saved.Time_Series,2) ;
    
    %% Convert the 2D iamges of spatial maps into the 4D structure and sum up the images of different bands
    % These maps are the corresponding fMRI maps of fused data and are
    % different from Multi_ICA_fMRI.SpatialICs. Some of these maps may not
    % mean anything. Only by using the FusedStuff.SpatialICs you can
    % understand how much the fMRI has contributed to every fused data
    % and if contribution is high then the map is meaningful. Otherwise it 
    % would be random mixture of several Multi_ICA_fMRI.SpatialICs maps.  
    ZSize = size(fMRI_Orig.img,3)/No_Freq_Bands ; 
    cfg.direction = 24 ;
    fMRI_Temp.Time_Series = FusedStuff.SpatialICs_fMRI ;  
    fMRI_Temp = FourDTwoDConvert(cfg, fMRI, Mask) 


    for I_I = 2:No_Freq_Bands
        fMRI_Temp.img(:,:,1:ZSize,:) = fMRI_Temp.img(:,:,(ZSize)*(I_I-1) + 1 : ZSize*I_I , :) + fMRI_Temp.img(:,:,1:ZSize,:) ;  
    end
    fMRI_Temp.img = fMRI_Temp.img(:,:,1:ZSize,:) ;
    fMRI_saved = fMRI_Temp;
    fMRI_saved.hdr.dime.dim(1,4) = ZSize ;
    fMRI_saved.hdr.dime.dim(1,5) = size(fMRI_saved.img,4) ;
    fMRI_saved = rmfield(fMRI_saved, 'Time_Series');
    fMRI_saved = rmfield(fMRI_saved, 'Time_Series_MeanRem');
    fMRI_saved = rmfield(fMRI_saved, 'Coordinate');
    clear fMRI_Temp ;
  
    
    %% Save the fMIR maps
    save_untouch_nii(fMRI_saved,  'Fused_fMRI_Component')

