%% 
clearvars;



%% Define Paths & Settings
addpath(genpath('dep'))
Par.Paths.out_path='/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/bstrasser/Temp3/LukisScriptTest2';
Par.Paths.csi_path{1} = '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/bstrasser/Projects/Project9_ImplementRecoInICE/Step5_MultiCenterStudy/LargeData_d3hj/MeasAndLogData/UCSF/Volunteer5_20250912/meas_MID00033_FID33731_csi_fid_ViennaCrt_v1_01.dat';
% Par.Paths.csi_path{1} = '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Measurement_Data/specBONN_ptx/meas_MID00176_FID18878_csi_fid_ViennaCr_2D_cp.dat';
Par.Paths.T1wImage='/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/bstrasser/Projects/Project9_ImplementRecoInICE/Step5_MultiCenterStudy/LargeData_d3hj/MeasAndLogData/UCSF/Volunteer5_20250912/MP2RAGE_0_7MM_TR4500_UNI_IMAGES_0006';
Par.Paths.T1wImage_AntiNoise='/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/bstrasser/Projects/Project9_ImplementRecoInICE/Step5_MultiCenterStudy/LargeData_d3hj/MeasAndLogData/UCSF/Volunteer5_20250912/MP2RAGE_0_7MM_TR4500_INV2_0007';
bet_path='/usr/local/fsl/bin/bet';

Par.Flags.NuisRem_WALINET_flag = 0;
Par.Flags.noisedecorrelation_flag = 1;
Par.Flags.hamming_flag = 1;
Par.Flags.hamming_z_flag = 1;
Par.Settings.LipidDecon_L2BetaCorrFactor=0.0; %0.2 default

Par.Settings.GradientDelayPerTempInt_x = [12.42, 12.38, 10.14];
Par.Settings.GradientDelayPerTempInt_y = [10.27,10.75,8.99];
Par.Settings.hamming_factor = 100;
Settings.NonCartReco.DensComp.Normalize_flag = false;
Settings.NonCartReco.DensComp.AutoScale_flag = false;



Settings.ICE_flag = 0;
Settings.NonCartReco.DensComp.ICE_flag = 0;


%% Read Data

[csi, image, NoiseData,coilcompressscan] = io_ReadAndReshapeSiemensData(Par.Paths.csi_path{1});


%% BET MASKING from mp2rage
run create_brain_mask.m

%% average
if ~csi.Par.dicom_flag
    image = op_AverageMRData(image);        % Currently does nothing
    csi = op_AverageMRData(csi);        % Currently does nothing
end


%% Noise-Decorrelate Data
if ~csi.Par.dicom_flag 

    if(exist('NoiseCorrMatStruct','var'))
        if(exist('NoiseData','var'))
            Settings.CreateSyntheticNoise=1;
        end
        csi = op_PerformNoiseDecorrelation(csi,NoiseCorrMatStruct,Settings);
        if(exist('image','var'))
            Settings.CreateSyntheticNoise = 0;
            image = op_PerformNoiseDecorrelation(image,NoiseCorrMatStruct,Settings);
        end
    end
end


%% Do CoilCompression
if ~csi.Par.dicom_flag
    
    if(exist('coilcompressscan','var') && isfield(coilcompressscan,'CoilCompScan') && isfield(coilcompressscan.CoilCompScan,'Data') && isfield(csi.Par,'coilcompression_to_numb_coils') && csi.Par.coilcompression_to_numb_coils<csi.Par.total_channel_no_measured && csi.Par.coilcompression_to_numb_coils>0)
        Text=sprintf('\n\n Do coil compression to %d virtual coils!',csi.Par.coilcompression_to_numb_coils);
        disp(Text)
    
        %set water data equals to image from iMUsical PRescan
        for ii = 1:size(csi.Data,2); wuff.Data{ii} = reshape(coilcompressscan.CoilCompScan.Data{ii},[prod(size_MultiDims(coilcompressscan.CoilCompScan.Data{ii},1:5)) size(coilcompressscan.CoilCompScan.Data{ii},6)]); end
        water.Data = cat(1,wuff.Data{:});
        clear wuff
        
        K=csi.Par.coilcompression_to_numb_coils;
        [rest NbCoil] = size(water.Data);
    
    %     fprintf('\n\nCoil compression: going from ',num2str(NbCoil),' to ',num2str(K),' coils.');
    
    %      Data_kc=permute(reshape(sqz(sum(Water_cast(:,:,:,1:mrsiReconParams.NbPtForWaterPhAmp),4)),[NbCoil Na*Ns]),[2 1]);
         % Data_kc=permute(reshape(sqz(mrsiReconParams.SENSE),[NbCoil M*N]),[2 1]);
         % Data_kc=permute(reshape(sqz(sum(mrsiData_cast(:,:,:,1:mrsiReconParams.NbPtForWaterPhAmp),4)),[NbCoil NbAng*NbTilt]),[2 1]);
    
         [~,~,V_cK]=svd(water.Data,0);
         V_cK=double(V_cK(:,1:K));
    
          %Water data
    %      Temp=zeros([rest K],'single');
    %      for NewCoil=1:K
    %          Temp(:,NewCoil)=single(water.Data)*V_cK(:,NewCoil);
    %      end
    %      image.Data = Temp;
        
          for j=1:size(image.Data,2)
             Temp{j}=zeros(  [size_MultiDims(image.Data{j},1:5) K],'double');
             Temp2=reshape(image.Data{j},[prod(size_MultiDims(image.Data{j},1:5)) size(image.Data{j},6)]);
             for NewCoil=1:K
                 Temp{j}(:,:,:,:,:,NewCoil)=reshape(double(Temp2)*V_cK(:,NewCoil),[size_MultiDims(image.Data{j},1:5)]);             
             end
             clear Temp2
             image.Par.DataSize{j}(6)=K;
             image.RecoPar.DataSize{j}(6)=K;
         end
            image.Data = Temp;
            clear Temp
            clear Temp2
            %Main Data
         
         for j=1:size(csi.Data,2)
             Temp{j}=zeros(  [size_MultiDims(csi.Data{j},1:5) K],'double');
             Temp2=reshape(csi.Data{j},[prod(size_MultiDims(csi.Data{j},1:5)) size(csi.Data{j},6)]);
             for NewCoil=1:K
                 Temp{j}(:,:,:,:,:,NewCoil)=reshape(double(Temp2)*V_cK(:,NewCoil),[size_MultiDims(csi.Data{j},1:5)]);
             end
             clear Temp2
             csi.Par.DataSize{j}(6)=K;
             csi.RecoPar.DataSize{j}(6)=K;
         end
        csi.Data = Temp;
        clear Temp
        clear Temp2
        
    %     for j=1:size(csi.Data,2)
    %          Temp{j}=zeros(  [size_MultiDims(csi.Data{j},1:5) K],'double');
    %          Temp2=reshape(csi.Data{j},[prod(size_MultiDims(csi.Data{j},1:5)) size(csi.Data{j},6)]);
             for NewCoil=1:K
                 Temp(:,NewCoil)=double(NoiseData.Data)*V_cK(:,NewCoil);
             end
             NoiseData.Par.DataSize(2)=K;
    %      end
        NoiseData.Data = Temp;
            clear Temp
     
        
        image.RecoPar.total_channel_no_measured=K;
        image.RecoPar.total_channel_no_reco=K;
        image.Par.total_channel_no_measured=K;
        image.Par.total_channel_no_reco=K;
        csi.RecoPar.total_channel_no_measured=K;
        csi.RecoPar.total_channel_no_reco=K;
        csi.Par.total_channel_no_measured=K;
        csi.Par.total_channel_no_reco=K;
        NoiseData.RecoPar.total_channel_no_measured=K;
        NoiseData.RecoPar.total_channel_no_reco=K;
        NoiseData.Par.total_channel_no_measured=K;
        NoiseData.Par.total_channel_no_reco=K;
        
    
    
    
        Par.CSI.total_channel_no=K;
    
    end
end



%% Calc NoiseCorrMat
if ~csi.Par.dicom_flag
    
    if(isfield(NoiseData,'Data') && numel(NoiseData.Data) > 1 && Par.Flags.noisedecorrelation_flag)
        [NoiseCorrMatStruct,Dummy] = op_CalcNoiseCorrMat(NoiseData);
        NoiseData = Dummy.NoiseData; clear Dummy;
    end
    if(~isfield(image,'Data'))
        clear image;
    end

end



%% FT
if ~csi.Par.dicom_flag
    
    Settings.NonCartReco.DensComp.Method = 'ConcentricRingTrajectory_Theoretical';   
    Settings.NonCartReco.DensComp.ApplyHammingFilter_flag = Par.Flags.hamming_flag;
    Settings.NonCartReco.DensComp.ApplyHammingFilter_z_flag = Par.Flags.hamming_z_flag;
    Settings.NonCartReco.ConjInkSpace_flag = true; 
    Settings.NonCartReco.Phaseroll_flag = true;
    Settings.NonCartReco.ConjIniSpace_flag = false;    
    Settings.NonCartReco.FlipDim = 1;
    Settings.NonCartReco.CircularSFTFoV_flag = false; 
    Settings.PreWhitenData_flag = 0;
    
    if(exist('image','var'))
        if(image.Par.SpatialSpectralEncoding_flag) 
    
            if(isfield(Par.Settings,'GradientDelayPerAngInt_x') && isfield(Par.Settings,'GradientDelayPerAngInt_y'))
                Settings.ReadInTraj.GradDelayPerAngInt_x_us = Par.Settings.GradientDelayPerAngInt_x;
                Settings.ReadInTraj.GradDelayPerAngInt_y_us = Par.Settings.GradientDelayPerAngInt_y;
            
            elseif(isfield(Par.Settings,'GradientDelayPerTempInt_x') && isfield(Par.Settings,'GradientDelayPerTempInt_y'))
                Settings.ReadInTraj.GradDelayPerTempInt_x_us = Par.Settings.GradientDelayPerTempInt_x;
                Settings.ReadInTraj.GradDelayPerTempInt_y_us = Par.Settings.GradientDelayPerTempInt_y;
            else
                Settings.ReadInTraj.GradDelayPerTempInt_x_us = 0;
                Settings.ReadInTraj.GradDelayPerTempInt_y_us = 0;
            end            
            image = op_ReconstructMRData(image,struct(),Settings);
    
        end
    
    end
    
        
    if(exist('image','var'))
        if(size(image.Data,4) > 4)
            image.Data = image.Data(:,:,:,5,:);	% change to 4
        else
            image.Data = image.Data(:,:,:,1,:);
        end
    end


% Reco MRSI Data
    csi = op_ReconstructMRData(csi,struct(),Settings);

end




%% PRE COIL COMBI BILGIC LIPID SUPP

% DEBUG DELETE ME LATER
% mask_bak = mask; mask = MaskShrinkOrGrow(mask,5,1,1);
% run ./NuisanceRemoval_HSVD.m
% mask = mask_bak; clear mask_bak;
% DEBUG END
if ~csi.Par.dicom_flag
    
    Info=csi.Par;
    
    if(Par.Settings.LipidDecon_L2BetaCorrFactor>0.0)
        
        csi_bak = csi;
        clear csi;
        csi = csi_bak.Data;
        csi_bak.Data = [];
        % csi = permute(csi,[5 1 2 3 4 6 7]);
        
        Scalle = 1.1047e+09/norm(csi(:));
        csi=csi*Scalle; % It scales csi data similarly like with Lukas's pipeline --> it makes the L2 work with the same factor    
    
        fprintf('\n\nPerform Bilgic Lipid Decontamination: channelwise & sensitivity-weighted\n')
    
        mkdir([tmp '/LipidDecontamination'])
        mkdir([tmp '/LipidMaskNii'])
    
        x_dim=size(csi,1);
        y_dim=size(csi,2);
        z_dim=size(csi,3);
    
        lipid_mask_total=zeros(x_dim,y_dim,z_dim);
    
        f_range=size(csi,4);
    %     f_range_start=ceil(f_range*0.82);
    %     f_range_end=ceil(f_range*0.96);
        
        
        ppmRange = [1.65, 0.35];
        ppm = compute_chemshift_vector(csi_bak.RecoPar);
        f_range_start = FindClosestIndex(ppm,ppmRange,[],'min');
        f_range_end = f_range_start{2};
        f_range_start = f_range_start{1};
        
        Text=sprintf('\n\nVectorsize Points of the csi go from 0 to %d!',f_range);
        disp(Text)
        Text=sprintf('\nRemoving Lipids in the Range from %d to %d:\n',f_range_start,f_range_end);
        disp(Text)
    
        for nCha = 1:size(csi,5)     
    
            csi_temp = csi(:,:,:,:,nCha); 
            csi_lip = fftshift(fft(csi_temp(:,:,:,:),[],4),4);
    
            for a=1:x_dim
                for b=1:y_dim
                    for c=1:z_dim
                           ftspectra(a,b,c,:)=abs((fftshift(fft(csi_temp(a,b,c,:)))));
                           ftspectramax(a,b,c)=max(ftspectra(a,b,c,f_range_start:f_range_end));
                    end
                end
            end
    
            max_sig=0.30*max(ftspectramax(:));
            sens=ftspectramax;
            sens(sens<max_sig)=0;
            sens(sens~=0)=1;
            lipid_mask=sens;
    
            struct_elm=strel('disk',1,0);
            tumor_mask=imerode(mask,struct_elm);
            lipid_mask=lipid_mask-tumor_mask;
            lipid_mask(lipid_mask~=1)=0;
    
            struct_elm=strel('disk',6,0);
            corner_mask=imdilate(mask,struct_elm);
            corner_mask(corner_mask~=1)=2;
            corner_mask(corner_mask~=2)=0;
            corner_mask(corner_mask~=0)=1;
    
            lipid_mask=lipid_mask-corner_mask;
            lipid_mask(lipid_mask~=1)=0;
    
            lipid_mask_total=lipid_mask_total+lipid_mask;
    
            save(sprintf('%s/LipidMaskNii/lipid_mask_chn_%d.mat',tmp,nCha),'lipid_mask');  
    
            mkdir([sprintf('%s/LipidDecontamination/Channel_%d', tmp,nCha)])
    
    %         for Slc = 1:size(csi_temp,3)
    %             Lipid_fig = figure('visible','on');
    %             imagesc(squeeze(lipid_mask(:,:,Slc)),[0 1])
    %             colorbar;      
    %             saveas(Lipid_fig,sprintf('%s/LipidDecontamination/Channel_%d/lipid_mask_slc_%d', tmp,nCha,Slc),'epsc2')
    %             %saveas(Lipid_fig,sprintf('%s/LipidDecontamination/Channel_%d/lipid_mask_slc_%d', tmp,nCha,Slc),'fig')            
    %             close(Lipid_fig)
    % 
    %             Brain_fig = figure('visible','on');
    %             imagesc(squeeze(abs(csi(:,:,Slc,1,nCha))),[0 max_sig])
    %             colorbar;      
    %             saveas(Brain_fig,sprintf('%s/LipidDecontamination/Channel_%d/brain_mask_slc_%d', tmp,nCha,Slc),'epsc2')
    %             %saveas(Brain_fig,sprintf('%s/LipidDecontamination/Channel_%d/brain_mask_slc_%d', tmp,nCha,Slc),'fig')
    %             close(Brain_fig)
    %         end
    
            fprintf('Decontaminating channel %d\n', nCha)    
    
            for Slc = 1:size(csi_temp,3)
    
                % Prepare slices: Get slice csi data, Extract Relevant Lipids, Get slice mask
                csi_slc = fftshift(fft(squeeze_single_dim(csi_temp(:,:,Slc,:),3),[],3),3);
    
                    % Get L2 Parameter
                    param.beta =  Par.Settings.LipidDecon_L2BetaCorrFactor * 10^-15; 
    
                    % Get L2 Lipids
                    if(1)%Par.CSI.ThreeD_flag || 1)
                        param.Lipid = [];
                        for SlcAlias = 1:size(csi_temp,3)
                            param.Lipid = cat(1, param.Lipid, get_LipidBasis(squeeze(csi_lip(:,:,SlcAlias,:)),lipid_mask(:,:,SlcAlias)));
                        end
    %                 else
    %                     param.Lipid = get_LipidBasis(csi_slc,lipid_mask(:,:,Slc));
                    end
    
                    param.Lipid = transpose(param.Lipid);
    
    %                 % Get L2 Brainmask
    %                 if(Par.Flags.InterpolateCSIResolution_flag && ~Par.Settings.InterpolateCSIResolution_InkSpace)
    %                     param.Bmask = mask_BefInterpol(:,:,Slc);
    %                 else
                        param.Bmask = mask(:,:,Slc);            
    %                 end
    
                % L2 Regularization
                csi_slc = LipidDecon_L2(csi_slc,param);
                csi(:,:,Slc,:,nCha) = ifft(fftshift(csi_slc,3),[],3);
    
            end
    
            clear csi_slc csi_lip sens Slc SlcAlias param
    
        end
    
        save(sprintf('%s/LipidMaskNii/lipid_mask_total.mat',tmp),'lipid_mask_total'); 
    
        mkdir([sprintf('%s/LipidDecontamination/Total', tmp)])    
    
    %     for Slc = 1:size(csi_temp,3)
    %             Lipid_fig_total = figure('visible','on');
    %             imagesc(squeeze(lipid_mask_total(:,:,Slc)),[0 nCha])
    %             colorbar;      
    %             saveas(Lipid_fig_total,sprintf('%s/LipidDecontamination/Total/lipid_mask_total_slc_%d', tmp,Slc),'epsc2')
    % %             saveas(Lipid_fig_total,sprintf('%s/LipidDecontamination/Total/lipid_mask_total_slc_%d', tmp,Slc),'fig')            
    %             close(Lipid_fig_total)
    %     end    
    
        csi_bak.Data = double(csi)/Scalle;
        csi = csi_bak;
        clear csi_bak;
        
        
    end   
    
    
    fprintf('\n\nAfter L2');  
end



%% Coil Combine MRSI Data
if ~csi.Par.dicom_flag

    weights = image; 
    weights.Data = conj(weights.Data(:,:,:,1,:));
    csi = op_CoilCombineData(csi,weights);
    
    % Rescale csi
    csi.Data = csi.Data * 10^5;
    if(isfield(csi,'NoiseData'))
        csi.NoiseData = csi.NoiseData * 10^5;
    end
    
    size_csi = size(csi.Data);
    % csi_struct=csi;
end



%% LCModel Fitting
% Info: This works only if LCModel is installed on your computer

Paths.LCM_Path = '/home/bstrasser/.lcmodel/bin/lcmodel';
Paths.TmpDir = '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/bstrasser/Temp3/LukisScriptTest1/Tmp';
Paths.BasisSetFile = '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/lab/Basis_Sets/LCModelOutput_MM_Param_whole_AMARES_phased/fid_1.300000ms.basis';


% % FittingMask = mask;
% % FOR DEBUGGING:
FittingMask = zeros(size_MultiDims(csi.Data,1:3));
FittingMask(29:34,29:34,1) = 1;
% % FOR DEBUGGING END


Settings.AllowNonEmptyDirs = true;
ControlFile = './LCModel_Control_Volunteers_4_2_to_1_8.m';

Paths.OutputDir = '/ceph/mri.meduniwien.ac.at/departments/radiology/mrsbrain/home/bstrasser/Temp3/LukisScriptTest2/LCModelOutDir';
[csi_FitResMaps,csi_FitResSpec] = op_PerformLCModelFitting(csi,FittingMask,Paths,ControlFile,Settings);

mkdir([Par.Paths.out_path '/maps'])
save([Par.Paths.out_path '/maps/AllMaps.mat'],'csi_FitResMaps','csi_FitResSpec')


%% Create Nifti Files from LCM-Fitted Maps
% NOT YET IMPLEMENTED


%% Save Data

clearvars -except Par csi weights image image_FullFID image_VC mask g_FactorMap NoiseScalingMatrix_spectral NoiseScalingMatrix_time IsWatRef NoiseCorrMatStruct NoiseData

save([Par.Paths.out_path '/CombinedCSI.mat'], '-v7.3')

