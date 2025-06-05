clc
clear all
close all

cd /home/ybr19/Reconstruction/SynthesizeData
addpath(genpath(cd));

%%% load data
studies_20201116_DISORDER_LowMotion;
addpath(genpath('/home/ybr19/Data/2020-11-16_DISORDER_LowMotion'))

% Select studies to reconstruction
id_path= 1; id_file= [3]; id_ref = [2]; %first scan with second lowres
[pathIn, refIn, fileIn] = extract_studies(pathIn, refIn, fileIn, id_path, id_ref, id_file);
ss = load(strcat(pathIn{1},'An-Ve',filesep,fileIn{1}{1},'_DephCorr001.mat'));

%% load data
x = ss.rec.d; %clear ss;
x = dynInd(x,1,4);
D = ss.rec.E.Dl.D; %clear ssD
%D = permute(D,[3 2 1 4 5 6]);

plot_(x,dynInd(D,1,6),dynInd(D,2,6),[],[-10 10],[],1,0)

%% Calculate the spectral terms
%%% create a 3D cartesian grid
N=size(x);
rangex = (1:N(1)) - ( floor(N(1)/2)); rangex = rangex / (N(1)/2);
rangey = (1:N(2)) - ( floor(N(2)/2)); rangey = rangey / (N(2)/2);
rangez = (1:N(3)) - ( floor(N(3)/2)); rangez = rangez / (N(3)/2);

rangex= pi * rangex;
rangey= pi * rangey;
rangez= pi * rangez;

[Kx, Ky, Kz] = ndgrid(rangex, rangey, rangez);
%plot_(Kx, Ky, Kz,[],[],[0 5],1,0)

D1=dynInd(D,1,6);
D2=dynInd(D,2,6);

S1 = ifftshift(ifftn(Ky .* fftn(fftshift(D1))));
S2 = ifftshift(ifftn(-Kx .*  fftn(fftshift(D2))));

plot_(x, abs(S1), abs(S2),[],[-10 10],[0 5],1,0)

grad1 =  real( ifftshift( ifftn( conj(Ky).* (Ky .* fftn(fftshift(D1)) + Kx .* fftn(fftshift(D2)))) ) );
grad2 =  real( ifftshift( ifftn( conj(Kx).* (Ky .* fftn(fftshift(D1)) + Kx .* fftn(fftshift(D2)))) ) );

rr = [-10 10]; rr_error = [];
plot_(x, grad1, grad2,[],rr, rr_error,1,0)


%% Do quick GD optimisation to s
close all
Drecon = ss.rec.E.Dl.D; %clear ssD
D = ss.rec.E.Dl.D;

niter = 30;
lambda = 0.5;

EnPrevious = loss(Drecon, D, lambda,x);
alpha_list = 10.^(-1* [2:1.5:35]);

EnList = EnPrevious;
for n = 1:niter
    fprintf('Iteration %d\n', n);
    grad = gradientD(Drecon, D, lambda,x);

    for alpha = alpha_list
        
       update = -alpha * grad;
       En = loss(Drecon+update, D, lambda,x);
       
       if En < EnPrevious
           EnPrevious = En;
           Drecon = Drecon + update;
           break
       end
    end
    
    if alpha == alpha_list(end)
        disp('no step size found');
        break;
    end
    EnList = [EnList En];
    
    figure(100)
    plot(1:length(EnList), EnList)
    
end

plot_(x,dynInd(D,1,6),dynInd(Drecon,1,6),[],[-8 8],[],1,0)
plot_(x,dynInd(D,2,6),dynInd(Drecon,2,6),[],[-8 8],[],1,0)
