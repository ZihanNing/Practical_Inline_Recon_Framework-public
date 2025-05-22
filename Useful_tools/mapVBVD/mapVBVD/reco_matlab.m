function  [output_img] = reco_matlab(varargin)

%% Input handling:

% addpath(genpath(fullfile(getenv('USERPROFILE'),'Documents','MATLAB','RESOURCES')));
path_to_add = genpath(fullfile('..','MATLAB','RESOURCES'));
addpath(path_to_add);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% default input parameters:
opt.raw_data_dir = '.';
opt.raw_data_file = '*.dat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% user-defined input parameters:
idx_arg = 1;
if (nargin > 0)
   while (idx_arg<=nargin && ischar(varargin{idx_arg}) && ~isempty(varargin{idx_arg}) && ~isempty(regexp(varargin{idx_arg}, '^-', 'once')))  
      
      if (strcmpi(varargin{idx_arg}(2:end), 'raw_data_dir'))
          opt.raw_data_dir = varargin{idx_arg+1};
          idx_arg = idx_arg + 2;
          continue;
      end
      
      if (strcmpi(varargin{idx_arg}(2:end), 'raw_data_file'))
          opt.raw_data_file = varargin{idx_arg+1};
          idx_arg = idx_arg + 2;
          continue;
      end
      
      
      %%% This stage should not be reached:
      error ('\nUnrecognized option %s ! Aborting.\n', varargin{idx_arg});
      
   end
end



%% Raw data processing and fft:
temp = dir(fullfile(opt.raw_data_dir,opt.raw_data_file));
if isempty(temp)
    error('\nRaw data file ''%s'' not found in ''%s''.\n',opt.raw_data_file,opt.raw_data_dir);
end
% Pick last (most recent) file
opt.raw_data_file = temp(end);

fprintf('\n\nExtracting all data from file:\n\t%s\n\n',opt.raw_data_file.name);
twix_full = mapVBVD(fullfile(opt.raw_data_dir,opt.raw_data_file.name));
twix = twix_full{end};

twix.image.flagDoAverage = true;
twix.image.flagRemoveOS  = true;


% data dimensions: RO x CHA x PE x PAR:
data = twix.image(:,:,:,:);

%%% account for asymmetric echo:
if (twix.hdr.MeasYaps.sKSpace.ucAsymmetricEchoMode  && size(data,1) < twix.hdr.MeasYaps.sKSpace.lBaseResolution)
    
    data = recoAsymmetricEcho(data,twix.hdr,'conjugate_synth');
    
else
    
    % fft in col and lin dimension:
    fft_dims = [1 3 4];
    for f = fft_dims
        data = ifftshift(ifft(fftshift(data,f),[],f),f);
    end
    
end

img = data; 


%% Building combined image:

%%%% Complex combined data is the complex sum of receive channels 
imgCpx = squeeze(sum(img,2));
%%%% Magnitude --> SoS
imgMag = squeeze(sqrt(sum(abs(img).^2,2)));



%% Output and cleanup:

output_img = imgMag .* exp(1i*angle(imgCpx));
rmpath(path_to_add);


