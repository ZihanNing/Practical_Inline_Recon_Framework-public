addpath(genpath('/home/lcg13/Work/DeepNetworks'));

x=readNIICustomAbs('/home/lcg13/Data/rawDestin/mi_23022018_1032032_7_2_sf7zt2mbzax0sensemi8s1V4_Aq.nii');

N=size(x);
x=resampling(x,N/2,[],2*ones(1,3));
x=repmat(x,[1 1 1 4 2]);
vis4DArray(x(:,:,:,:,1),[],6)
%y=augmentSliceRotation(x,1,pi/8);
%y=augmentSliceTranslation(x,1,[10 10]);
%y=augmentSliceBiasField(x,1,[4 4],0.2);
%y=augmentSlicePermuting(x,1,4);
%y=augmentTranslation(x,[10 10 10]);
%y=augmentShear(x,3,1);
%y=augmentScale(x);
%y=augmentRotation(x,3,2*pi);
%y=augmentFlip(x,1);
y=augmentBiasField(x,[3 3 3],0.1);


vis4DArray(y(:,:,:,:,1),[],6)
