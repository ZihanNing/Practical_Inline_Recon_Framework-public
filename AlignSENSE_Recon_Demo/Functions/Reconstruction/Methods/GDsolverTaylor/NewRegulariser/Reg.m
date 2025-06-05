clc
clear all
close all

cd /home/ybr19/Reconstruction/SynthesizeData
addpath(genpath(cd));

%%% load data
studies_20201002_DISORDER_MotionFree;
addpath(genpath('/home/ybr19/Data/2020-10-02_DISORDER_MotionFree'))

% Select studies to reconstruction
id_path= 1; id_file= [1]; id_ref = [2]; %first scan with second lowres
[pathIn, refIn, fileIn] = extract_studies(pathIn, refIn, fileIn, id_path, id_ref, id_file);
ss = load(strcat(pathIn{1},'An-Ve',filesep,fileIn{1}{1},'MotCorr001.mat'));
ssD = load('D.mat');

%%
xGT = ss.rec.d; clear ss;
DGT = ssD.D; clear ssD
DGT = permute(DGT,[3 2 1 4]);

N=size(xGT);
ssD = load('D.mat');

W = ssD.W_delin_deSH; clear ssD
DGT2 = DGT;for i=1:2; DGT2 = dynInd(DGT2,i,4, reshape(gather(W(i,:)), N));end
DGT2 = permute(DGT2,[3 2 1 4]);

plot_(xGT,dynInd(DGT2,1,4),dynInd(DGT2,2,4),[],[-15 15],[],1,0)

%%% create a 3D cartesian grid
rangex = (1:N(1)) - ( floor(N(1)/2)); rangex = rangex / (N(1)/2);
rangey = (1:N(2)) - ( floor(N(2)/2)); rangey = rangey / (N(2)/2);
rangez = (1:N(3)) - ( floor(N(3)/2)); rangez = rangez / (N(3)/2);

rangex= pi * rangex/2;
rangey= pi * rangey/2;
rangez= pi * rangez/2;

[Kx, Ky, Kz] = ndgrid(rangex, rangey, rangez);

D1=dynInd(DGT2,1,4);
D2=dynInd(DGT2,2,4);

S1 = ifftshift(ifftn(Kx.* fftn(fftshift(D1))));
S2 = ifftshift(ifftn(-Ky .*  fftn(fftshift(D2))));

plot_(xGT, S1, S2,[],[-10 10],[-10 10],1,0)
