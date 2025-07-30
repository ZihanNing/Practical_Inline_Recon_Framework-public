function setPath
% Add all subdirectories to path
%
% (c) Corey Baron

  % Delete any parallel pools, because the workers won't get the new path otherwise
  if(exist('gcp'))
    poolobj = gcp('nocreate');
    delete(poolobj);
    clear poolobj
  end

  topPath = mfilename('fullpath');
  topPath = topPath(1:end-length(mfilename));

  % Add paths
  addpath([topPath, 'bview']);
  addpath([topPath, 'computeHarmonics']);
  addpath([topPath, 'demos']);
  addpath([topPath, 'findDel']);
  addpath([topPath, 'general']);
  addpath([topPath, 'iterativeSolvers']);
  addpath([topPath, 'operators']);
  addpath([topPath, 'simulation']);
  addpath([topPath, 'trajectory']);
  addpath([topPath, 'unitTests']);

  
   % Set matMRI version
  global matMRIVersion
  matMRIVersion = 1.07;
end