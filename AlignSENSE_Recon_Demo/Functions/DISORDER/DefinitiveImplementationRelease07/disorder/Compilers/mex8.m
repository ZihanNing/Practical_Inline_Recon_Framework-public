function mex8()

% need to be in the current directory for mexcuda
oldpath = pwd;
%newpath = '/usr/local/MATLAB/R2015a/extern/include';%fileparts(mfilename('fullpath'));
%cd(newpath);

% if the mexcuda fails, we are stuck - rethrow error
try
    mex_all_compile();
    cd(oldpath)
catch ME
    cd(oldpath)
    rethrow(ME)
end

function mex_all_compile()

% clean
%delete mexGPU8SVD.mex*

%% default mexcuda (cuda 7.5)

%mexcuda mexGPU8SVD.cu -I/usr/local/cuda/include -L/usr/local/cuda/lib64 LDFLAGS='"$LDFLAGS -Wl,--no-as-needed"' -ldl -lcusparse -lcublas_static -lcusparse_static -lculibos

% not used any more
%mexcuda csrsort.cu -I/usr/local/cuda/include -L/usr/local/cuda/lib64 LDFLAGS='"$LDFLAGS -Wl,--no-as-needed"' -ldl -lcusparse -lcublas_static -lcusparse_static -lculibos


%% cuda 8: csrmv.cu (cusparseScsrmv_mp) doesn't work
if true
    files={'../Method/fft_kernel'};
    
    %files = {'mexGPU8SVD','Utilities'};
    
    %matFolder='/usr/local/MATLAB/R2016b/';
    matFolder='/usr/local/MATLAB/R2015a/';
    %archi='-gencode arch=compute_52,code=sm_52';
    cudaFolder='/usr/local/cuda-8.0/targets/x86_64-linux/lib/';
    %matFolder='/usr/local/MATLAB/R2016a/';
    archi='';
    %cudaFolder='/usr/local/cuda-8.0/lib64/';
    
    for k = 1:numel(files)

        cmd1 = ['/usr/local/cuda-8.0/bin/nvcc -c -Xptxas="-v" --compiler-options=-D_GNU_SOURCE,-DMATLAB_MEX_FILE' ...
               ' -I"/usr/local/cuda-8.0/include"' ...
               ' -I"' matFolder 'extern/include" -I"' matFolder 'simulink/include"' ...
               ' -I"' matFolder 'toolbox/distcomp/gpu/extern/include/" -I"' matFolder '/extern/include/"' ...
               ' ' archi '-O3 --use_fast_math --default-stream per-thread --ptxas-options=-v' ...
               ' -std=c++11 --compiler-options=-ansi,-fexceptions,-fPIC,-fno-omit-frame-pointer,-pthread ' ...
               ' ' files{k} '.cu -o ' files{k} '.o'];
           cmd1
        
        cmd2 = ['/usr/bin/g++ -pthread -Wl,--no-undefined -Wl,--no-as-needed -shared -O3' ...
               ' -Wl,--version-script,"' matFolder 'extern/lib/glnxa64/mexFunction.map"' ...
               ' ' strcat(files{k},'.o') ' -ldl' ... 
               ' ' cudaFolder 'libcusparse.so' ... % -lcusparse
               ' ' cudaFolder 'libcusolver.so' ... % -lcusolver
               ' ' cudaFolder 'libcusolver_static.a' ... % -lcusolver_static
               ' ' cudaFolder 'libcufft.so' ... % -lcusolver
               ' ' cudaFolder 'libcufft_static.a' ... % -lcusolver_static
               ' ' cudaFolder 'libcublas_static.a' ... % -lcublas_static
               ' ' cudaFolder 'libcusparse_static.a' ... % -lcusparse_static
               ' ' cudaFolder 'libculibos.a' ... % -lculibos'
               ' -L/usr/local/cuda-8.0/lib64 -Wl,-rpath-link,' matFolder 'bin/glnxa64' ...
               ' -L"' matFolder 'bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -lmwgpu' ...
               ' ' cudaFolder 'libcudart.so' ... % /usr/local/MATLAB/R2016b/bin/glnxa64/libcudart.so.7.5 
               ' -o ' files{k} '.mexa64'];
           cmd2
        
        cmd3 = ['rm -f ' files{k} '.o'];

        disp([files{k} '.cu'])
        if system(cmd1); error('%s failed step 1',files{k}); end
        if system(cmd2); error('%s failed step 2',files{k}); end
        %if system(cmd3); error('%s failed step 3',files{k}); end
        disp('MEX completed successfully.')
        
    end
end