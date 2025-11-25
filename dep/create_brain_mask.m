
%%% MAKE FOV MASK (not VOI mask for now)

tmp=[pwd '/temp'];
if ~exist(tmp, 'dir')
    mkdir(tmp);
else
    rmdir(tmp,'s');
    mkdir(tmp);
end
%% Stepsize (Voxel Size) and Orientations
csi.Par.StepRead = -csi.Par.FoV_Read(1) / csi.Par.nFreqEnc;		% Coordinate system is reversed in minc with respect to DICOM
csi.Par.StepPhase = -csi.Par.FoV_Phase(1) / csi.Par.nPhasEnc;    
if(csi.Par.nPartEnc == 1)
    csi.Par.FoV_Partition(1) = csi.Par.VoI_Partition(1);
end
csi.Par.StepSlice = csi.Par.FoV_Partition(1) / (csi.Par.nPartEnc * csi.Par.nSLC); 

% Rename Position Fields
[csi.Par.POS_X] = csi.Par.Pos_Sag;
[csi.Par.POS_Y] = csi.Par.Pos_Cor;
[csi.Par.POS_Z] = csi.Par.Pos_Tra;
% csi.Par = rmfield(csi.Par,{'Pos_Sag','Pos_Cor','Pos_Tra'});

% Only for CRT?   
csi.Par.InPlaneRotation=csi.Par.InPlaneRotation;
csi.Par.InPlaneRotation_VOI=csi.Par.InPlaneRotation_VOI;

% Compute direction cosine from x,y and z components of slice normal vector
[csi.Par.PhaseNormalVector, csi.Par.ReadNormalVector] = compute_dircos([csi.Par.SliceNormalVector_x(1) csi.Par.SliceNormalVector_y(1) csi.Par.SliceNormalVector_z(1)],csi.Par.InPlaneRotation);
csi.Par.SliceNormalVector = [csi.Par.SliceNormalVector_x(1) csi.Par.SliceNormalVector_y(1) csi.Par.SliceNormalVector_z(1)];    

MinusVec1 = [-1 -1 1];        % MinusVec are here only for "tuning" signs of the final rotation matrix (RotMat). The values of MinusVecs are 
MinusVec2 = [1 1 -1];         % empirical, thus there might exist a dataset, which will have signs of RotMat uncorrect. However, this setup
MinusVec3 = [-1 -1 +1];       % works for all tested datasets. The uncorrect signs of direction cosines might be caused by the fact, that the attribute
						  % Image Orientation (Patient) tag(0020,0037) is not used in computation of dircos. Instead the Siemens private Tag (0029, 1020) 
						  % is used


csi.Par.ReadNormalVector = csi.Par.ReadNormalVector .* MinusVec1;
csi.Par.PhaseNormalVector = csi.Par.PhaseNormalVector .* MinusVec2;
csi.Par.SliceNormalVector = csi.Par.SliceNormalVector .* MinusVec3;

% Create rotation matrix
RotMat = cat(1,csi.Par.ReadNormalVector,csi.Par.PhaseNormalVector,csi.Par.SliceNormalVector);

% Reverse x- and y- coordinates due to the reversed coordinate system of minc with respect to DICOM
Pos = [-csi.Par.POS_X(1),-csi.Par.POS_Y(1),csi.Par.POS_Z(1)];

% Convert position from DICOM world coordinates to MINC start values 
% This approach is probably prone to extreme rotation of FOV and is not
% universal, however for the tested datasets it yielded correct results
Pos_Minc = RotMat * transpose(Pos); 
% The following line is old, and I think wrong. What I think happened is: Michal fixed the shift
% in the z-direction by half a voxel with the line below (effectively this line calculates
% Pos_z - FoV_z/2 + FoV_z/(2*N_z). Then I figured out that the x- and y-positions have to be
% shifted by half a voxel, and thought by analogy also the z-dimension has to be shifted, not
% knowing that Michal did that already with the FoVHalf. Now it should be fixed:
% The FoVHalf is defined "normally" also for z, and the half-voxel shift is done in the 3D-case.
% For checks: See git commits #1383, #004f, #a9be
% 	FoVHalf = [csi.Par.FoV_Read(1)/2 csi.Par.FoV_Phase(1)/2 -csi.Par.FoV_Partition(1)/csi.Par.nPartEnc*(csi.Par.nPartEnc-1)/2]; 
FoVHalf = [csi.Par.FoV_Read(1)/2 csi.Par.FoV_Phase(1)/2 -csi.Par.FoV_Partition(1)/2];  
Pos_Minc = transpose(Pos_Minc) + FoVHalf;

% Get from Center of Voxel (DICOM) to corner of voxel (minc) by subtracting half the voxel 
csi.Par.POS_X_FirstVoxel = Pos_Minc(1) + csi.Par.StepRead/2;         % Be aware that StepRead and StepPhase are reversed and thus the sum is effectively a subtraction.
csi.Par.POS_Y_FirstVoxel = Pos_Minc(2) + csi.Par.StepPhase/2;
csi.Par.POS_Z_FirstVoxel = Pos_Minc(3);


if(csi.Par.ThreeD_flag)
    csi.Par.POS_Z_FirstVoxel = csi.Par.POS_Z_FirstVoxel + csi.Par.StepSlice/2;     
end

%% UNIX code


duerer_orig = load('durer');
duerer_orig = duerer_orig.X;
duerer = imresize(duerer_orig,[csi.Par.nPhasEnc,csi.Par.nFreqEnc],'bicubic');
duerer = repmat(duerer,[1 1 csi.Par.nPartEnc*csi.Par.nSLC]);
duerer(:,:,1:2:end) = flip(duerer(:,:,1:2:end),2);                % Flip every second slice, so that we can distinguish between the slices in the mincfile
csi_template_fid = fopen([tmp '/csi_template.raw'],'w');
fwrite(csi_template_fid,duerer,'float');
fclose(csi_template_fid);

csi_path_allfiles = dir( fullfile(Par.Paths.T1wImage,'*.*') );
csi_path_allfiles = natsortfiles(csi_path_allfiles);
csi_path_allfiles = {csi_path_allfiles.name}';
DicomFoundVec = ~cellfun(@isempty,regexpi(csi_path_allfiles,'\.IMA|\.dcm'));
csi_path_allfiles = csi_path_allfiles(DicomFoundVec);
csi_path_allfiles = strcat(Par.Paths.T1wImage,'/',csi_path_allfiles);

if ~isempty(csi_path_allfiles)
    [bla,bla2] = unix(['dcm2mnc ' Par.Paths.T1wImage ' -dname "" ' tmp ' -fname magnitude']);
    [bla,bla2] = unix(['dcm2mnc ' Par.Paths.T1wImage_AntiNoise ' -dname "" ' tmp ' -fname antinoise']);
else
%DCM
    [bla,bla2] = unix(['dcm2niix -o ' tmp ' -z y -f magnitude ' Par.Paths.T1wImage]);
    % [bla,bla2] = unix(['dcm2niix -o ' tmp ' -z y -f antinoise ' Par.Paths.T1wImage_AntiNoise]);
    [bla,bla2] = unix(['gunzip ' tmp '/magnitude.nii.gz']);
    % [bla,bla2] = unix(['gunzip ' tmp '/antinoise.nii.gz']);
    [bla,bla2] = unix(['nii2mnc -quiet ' tmp '/magnitude.nii ' tmp '/magnitude.mnc']);
    % [bla,bla2] = unix(['nii2mnc -quiet ' tmp '/antinoise.nii ' tmp '/antinoise.mnc']);
end
% Normal
bashstring = ['rawtominc -clobber -float ' tmp '/csi_template.mnc -input ' tmp '/csi_template.raw'];
bashstring = sprintf('%s -xstep %8.6f -ystep %8.6f -zstep %8.6f',bashstring,csi.Par.StepRead, csi.Par.StepPhase, csi.Par.StepSlice);
bashstring = sprintf('%s -xstart %8.6f -ystart %8.6f -zstart %8.6f',bashstring, min(csi.Par.POS_X_FirstVoxel), min(csi.Par.POS_Y_FirstVoxel), min(csi.Par.POS_Z_FirstVoxel));
bashstring = sprintf('%s',bashstring);
bashstring = sprintf('%s -xdircos %8.6f %8.6f %8.6f -ydircos %8.6f %8.6f %8.6f -zdircos %8.6f %8.6f %8.6f', bashstring, csi.Par.ReadNormalVector, ...
csi.Par.PhaseNormalVector, csi.Par.SliceNormalVector );            
bashstring = sprintf('%s %d %d %d',bashstring,csi.Par.nPartEnc*csi.Par.nSLC,csi.Par.nPhasEnc,csi.Par.nFreqEnc);

% bashstring_image = ['rawtominc -clobber -float ' tmp '/mag_template.mnc -input ' tmp '/csi_template.raw'];
% bashstring_image = sprintf('%s -xstep %8.6f -ystep %8.6f -zstep %8.6f',bashstring_image, csi.Par.StepRead, csi.Par.StepPhase, csi.Par.StepSlice);    
% bashstring_image = sprintf('%s -xstart %8.6f -ystart %8.6f -zstart %8.6f',bashstring_image, min(csi.Par.POS_X_FirstVoxel), min(csi.Par.POS_Y_FirstVoxel), min(csi.Par.POS_Z_FirstVoxel));        
% bashstring_image = sprintf('%s -xdircos %8.6f %8.6f %8.6f -ydircos %8.6f %8.6f %8.6f -zdircos %8.6f %8.6f %8.6f',bashstring_image,csi.Par.ReadNormalVector, ...
% csi.Par.PhaseNormalVector,csi.Par.SliceNormalVector );            
% bashstring_image = sprintf('%s %d %d %d',bashstring_image,csi.Par.nPartEnc*csi.Par.nSLC,csi.Par.nPhasEnc,csi.Par.nFreqEnc);

fiddy = fopen([tmp '/CreateMincTemplates.sh'],'w+');
fprintf(fiddy, '%s\n%s', bashstring);%, bashstring_image);
fclose(fiddy);
fileattrib([tmp '/CreateMincTemplates.sh'],'+x','gu');   % Executable for group and user.

[bla,bla2] = unix(['bash ' tmp '/CreateMincTemplates.sh']);

if ~isempty(csi_path_allfiles)
    [bla,bla2] = unix(['max_magnitude=`mincstats -quiet -max ' tmp '/antinoise.mnc`; lower_threshold=$(echo "scale=6 ; ${max_magnitude}/84" | bc); mincmath -clobber -segment -const2 $lower_threshold $max_magnitude ' tmp '/antinoise.mnc ' tmp '/antinoise2.mnc']);
    [bla,bla2] = unix(['mincmath -clobber -mult ' tmp '/magnitude.mnc ' tmp '/antinoise2.mnc ' tmp '/magnitude2.mnc']);
    [bla,bla2] = unix(['rm ' tmp '/magnitude.mnc']);
    [bla,bla2] = unix(['rm ' tmp '/antinoise.mnc']);
    [bla,bla2] = unix(['mv ' tmp '/magnitude2.mnc ' tmp '/magnitude.mnc']);
    [bla,bla2] = unix(['mv ' tmp '/antinoise2.mnc ' tmp '/antinoise.mnc']);
end

[bla,bla2] = unix(['mnc2nii -quiet ' tmp '/magnitude.mnc ' tmp '/nii_magnitude.nii']);
[bla,bla2] = unix([bet_path ' ' tmp '/nii_magnitude ' tmp '/brain -f 0.5 -g 0.0 -n -m']);
[bla,bla2] = unix([bet_path ' ' tmp '/nii_magnitude ' tmp '/lipid -f 0.5 -g 0.0 -n -A']);
[bla,bla2] = unix(['gunzip ' tmp '/brain_mask.nii.gz']);
[bla,bla2] = unix(['gunzip ' tmp '/lipid_skull_mask.nii.gz']);
[bla,bla2] = unix(['gunzip ' tmp '/lipid_outskin_mask.nii.gz']);
[bla,bla2] = unix(['gunzip ' tmp '/lipid_outskull_mask.nii.gz']);
[bla,bla2] = unix(['gunzip ' tmp '/lipid_inskull_mask.nii.gz']);
[bla,bla2] = unix(['nii2mnc -quiet ' tmp '/brain_mask.nii ' tmp '/mask_brain_unres.mnc']);
[bla,bla2] = unix(['nii2mnc -quiet ' tmp '/lipid_outskin_mask.nii ' tmp '/OUTSKIN.mnc >/dev/null']);
[bla,bla2] = unix(['nii2mnc -quiet ' tmp '/lipid_outskull_mask.nii ' tmp '/OUTSKULL.mnc >/dev/null']);
[bla,bla2] = unix(['nii2mnc -quiet ' tmp '/lipid_inskull_mask.nii ' tmp '/INSKULL.mnc >/dev/null']);
%[bla,bla2] = unix(['rm ' tmp '/mask_brain_unres.mnc']);
[bla,bla2] = unix(['rm ' tmp '/mask_lipid_unres.mnc']);
[bla,bla2] = unix(['mincmath -sub ' tmp '/OUTSKIN.mnc ' tmp '/mask_brain_unres.mnc ' tmp '/mask_lipid_unres.mnc >/dev/null']);

% [bla,bla2] = unix(['minctoraw ' tmp '/mask_brain_unres.mnc -nonormalize -float > ' tmp '/mask_brain_unres.raw']);
% [bla,bla2] = unix(['rawtominc -float -clobber -like ' tmp '/mag_template.mnc -input ' tmp '/mask_brain_unres.raw ' tmp '/mask_brain_unres2.mnc']);
[bla,bla2] = unix(['mincresample -clobber -nearest_neighbour -like ' tmp '/csi_template.mnc ' tmp '/mask_brain_unres.mnc ' tmp '/mask_brain2.mnc']);
[bla,bla2] = unix(['minctoraw ' tmp '/mask_brain2.mnc -nonormalize -float > ' tmp '/mask_brain.raw']);
[bla,bla2] = unix(['mincresample -clobber -nearest_neighbour -like ' tmp '/csi_template.mnc ' tmp '/mask_lipid_unres.mnc ' tmp '/mask_lipid2.mnc']);
[bla,bla2] = unix(['minctoraw ' tmp '/mask_lipid2.mnc -nonormalize -float > ' tmp '/mask_lipid.raw']);

[bla,bla2] = unix(['rawtominc -float -clobber -like ' tmp '/csi_template.mnc -input ' tmp '/mask_brain.raw ' tmp '/mask_brain.mnc']);
[bla,bla2] = unix(['rawtominc -float -clobber -like ' tmp '/csi_template.mnc -input ' tmp '/mask_lipid.raw ' tmp '/mask_lipid.mnc']);

if(exist([tmp '/mask_brain.raw'],'file'))
        fid_mask = fopen([tmp '/mask_brain.raw'],'r');
        mask = reshape(fread(fid_mask, 'float'), [csi.Par.nPhasEnc csi.Par.nPhasEnc csi.Par.nPartEnc]);
        fclose(fid_mask);
end
if(exist([tmp '/mask_lipid.raw'],'file'))
        fid_mask = fopen([tmp '/mask_lipid.raw'],'r');
        lipid_mask_BET = reshape(fread(fid_mask, 'float'), [csi.Par.nPhasEnc csi.Par.nPhasEnc csi.Par.nPartEnc]);
        fclose(fid_mask);
end



%% Clean Up

if (exist(tmp, 'dir'))
    rmdir(tmp,'s');
end
