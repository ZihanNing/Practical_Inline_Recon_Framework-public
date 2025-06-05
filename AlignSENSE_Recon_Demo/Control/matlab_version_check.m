function matlab_version_check
    fprintf('  MATLAB Version: %s\n', version);
    fprintf('  MATLABROOT: %s\n', getenv('MATLABROOT'));
    fprintf('  PATH: %s\n', getenv('PATH'));
    fprintf('  LD_LIBRARY_PATH: %s\n', getenv('LD_LIBRARY_PATH'));
end