function [rawdata,loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17(measfile, dispopt)
%function [rawdata,loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17(measfile, dispopt)
% For Idea. Reads the 2D measurement data got from Siemens scanners following
% (mstart mate->3->10->2(Rawdata2Disk ON));
% I assume the file extesion .out.This mfunction extracts imaging information from the
% mini MDH in the .out file. No .asc file is needed.
%
% The .out file has
%     1) an uint32 indicating hdr_len (only from VB17);
%     2) an common header of hdr_len including the uint32 followed by all ADC data (only from VB15);
%     3) each ADC acquisition is oversampled by 2 and has a 128 byte header;
%     4) each pixel consists of the read and imaginary parts, which are float (4 bytes);
%     5) 384/1024? bytes are put at the end.
%
% dispopt, - 'on' or 'off'. The default is 'off'. The FFTed images for slices of the first 3D
% is displayed only for your reference. if your data were segmented, PF, ReadOutOffcentre, or
% arranged in the special way, you might not see the correct images. You have to turn this
% option off, and work on the rawdata output from this function for the right raw line order.
%
% Output:
% rawdata, - raw data in matrix; oversmapled (x2) in complex (A+iB), indexed by loopcounters and
% ulChannelId;
% loopcounters,- all the sMDH counters and ChannelId corresponding to all the ADC data lines;
% All dimension are in the same order as saved in sMDH.sLC;
% sMDH, - the MDH of the last ADC acquision except that extended sMDH.sLC is set to maximum
% loop values.
%
% EXAMPLE:
%           [rawdata loopcounters sMDH lc_names] = ReadSiemensMeasVB17_idea('Z:\n4\pkg\MeasurementData\FLASH\meas.dat');
%
% Maolin Qiu YALE 9-1-2011 (maolin.qiu(at)yale(dot)edu)

if nargin < 1 | ~exist(measfile)
    help ReadSiemensMeasVB17_idea;
    [filename pathname] = uigetfile( ...
       {'*.dat';'*.out';'*.*'}, ...
        'Pick a Siemens MESUREMENT file');
    if ~filename & ~pathname
        disp(['You selected no file.']);
        return;
    else
        measfile = fullfile(pathname, filename);
    end
end

if nargin < 2
    dispopt = 'off';
end

% measfile

[pa na ex] = fileparts(measfile);
measfile = fullfile(pa, [na ex]);
if exist(measfile, 'file')
    disp(['Measurement data file: ' measfile]);
else
    disp(['Measurement data file does not exist: ' measfile]);
    return;
end

%% init output

data = []; dimensions = []; loopcounters = []; sMDH = []; lc_names = [];

%% Constants used in sMDH

MDH_NUMBEROFEVALINFOMASK   = 2;
MDH_NUMBEROFICEPROGRAMPARA = 4;
MDH_FREEHDRPARA            = 4;

IDX_DIM                    = 30;
lc_names = {'001:sMDH.sLC.ushLine'
       '002:sMDH.sLC.ushAcquisition'
       '003:sMDH.sLC.ushSlice'
       '004:sMDH.sLC.ushPartition'
       '005:sMDH.sLC.ushEcho'
       '006:sMDH.sLC.ushPhase'
       '007:sMDH.sLC.ushRepetition'
       '008:sMDH.sLC.ushSet'
       '009:sMDH.sLC.ushSeg'
       '010:sMDH.sLC.ushIda'
       '011:sMDH.sLC.ushIdb'
       '012:sMDH.sLC.ushIdc'
       '013:sMDH.sLC.ushIdd'
       '014:sMDH.sLC.ushIde'
       '015:sMDH.sLC.ulChannelId'
       '016:sMDH.ushSamplesInScan'
       '017:sMDH.ushKSpaceCentreColumn'
       '018:sMDH.fReadOutOffcentre'
       '019:sMDH.aulEvalInfoMask(1)'
       '020:sMDH.aulEvalInfoMask(2)'
       '021:'
       '022:'
       '023:sMDH.ushKSpaceCentreLineNo'
       '024:sMDH.sCutOff.ushPre'
       '025:sMDH.sCutOff.ushPost'
       '026:sMDH.ushKSpaceCentrePartitionNo' %26
       '027:sMDH.ulDMALength' %27
       '028:sMDH.lMeasUID' %28
       '029:sMDH.ulScanCounter' %29
       '030:sMDH.ulTimeStamp'};
%%--------------------------------------------------------------------------%%
%% Definition of loop counter structure                                     %%
%% Note: any changes of this structure affect the corresponding swapping    %%
%%       method of the measurement data header proxy class (MdhProxy)       %%
%%--------------------------------------------------------------------------%%

sLoopCounter = struct( ...
  'ushLine',0,...                  %% unsigned short  line index                   %%
  'ushAcquisition',0,...           %% unsigned short  acquisition index            %%
  'ushSlice',0,...                 %% unsigned short  slice index                  %%
  'ushPartition',0,...             %% unsigned short  partition index              %%
  'ushEcho',0,...                  %% unsigned short  echo index                   %%
  'ushPhase',0,...                 %% unsigned short  phase index                  %%
  'ushRepetition',0,...            %% unsigned short  measurement repeat index     %%
  'ushSet',0,...                   %% unsigned short  set index                    %%
  'ushSeg',0,...                   %% unsigned short  segment index  (for TSE)     %%
  'ushIda',0,...                   %% unsigned short  IceDimension a index         %%
  'ushIdb',0,...                   %% unsigned short  IceDimension b index         %%
  'ushIdc',0,...                   %% unsigned short  IceDimension c index         %%
  'ushIdd',0,...                   %% unsigned short  IceDimension d index         %%
  'ushIde',0 ...                   %% unsigned short  IceDimension e index         %%
);                                 %% sizeof : 28 byte             %%

%%--------------------------------------------------------------------------%%
%%  Definition of slice vectors                                             %%
%%--------------------------------------------------------------------------%%

sVector = struct( ...
  'flSag',0.0,...       %% float
  'flCor',0.0,...       %% float
  'flTra',0.0 ...       %% float
);

sSliceData = struct( ...
  'sSlicePosVec',sVector,...                   %% slice position vector               %%
  'aflQuaternion',zeros(1,4) ...               %% float rotation matrix as quaternion %%
);                                              %% sizeof : 28 byte                    %%

%%--------------------------------------------------------------------------%%
%%  Definition of cut-off data                                              %%
%%--------------------------------------------------------------------------%%

sCutOffData = struct( ...
  'ushPre',0,...               %% unsigned short  write ushPre zeros at line start %%
  'ushPost',0 ...              %% unsigned short  write ushPost zeros at line end  %%
);

%%--------------------------------------------------------------------------%%
%%  Definition of measurement data header                                   %%
%%--------------------------------------------------------------------------%%

sMDH = struct( ...
  'ulDMALength',0,...                                       %% unsigned long  DMA length [bytes] must be                        4 bytes %% first parameter
  'lMeasUID',0,...                                          %% long           measurement user ID                               4
  'ulScanCounter',0,...                                     %% unsigned long  scan counter [1...]                               4
  'ulTimeStamp',0,...                                       %% unsigned long  time stamp [2.5 ms ticks since 00:00]             4
  'ulPMUTimeStamp',0,...                                    %% unsigned long  PMU time stamp [2.5 ms ticks since last trigger]  4
  'aulEvalInfoMask',zeros(1,MDH_NUMBEROFEVALINFOMASK),...   %% unsigned long  evaluation info mask field                        8
  'ushSamplesInScan',0,...                                  %% unsigned short # of samples acquired in scan                     2
  'ushUsedChannels',0,...                                   %% unsigned short # of channels used in scan                        2   =32
  'sLC',sLoopCounter,...                                    %% loop counters                                                    28  =60
  'sCutOff',sCutOffData,...                                 %% cut-off values                                                   4
  'ushKSpaceCentreColumn',0,...                             %% unsigned short centre of echo                                    2
  'ushDummy',0,...                                          %% unsigned short for swapping                                      2
  'fReadOutOffcentre',0.0,...                               %% float          ReadOut offcenter value                           4
  'ulTimeSinceLastRF',0,...                                 %% unsigned long  Sequence time stamp since last RF pulse           4
  'ushKSpaceCentreLineNo',0,...                             %% unsigned short number of K-space centre line                     2
  'ushKSpaceCentrePartitionNo',0,...                        %% unsigned short number of K-space centre partition                2
  'aushIceProgramPara',zeros(1,MDH_NUMBEROFICEPROGRAMPARA),... %% unsigned short free parameter for IceProgram                  8  =88
  'aushFreePara',zeros(1,MDH_FREEHDRPARA),...               %% unsigned short free parameter                                    4 * 2 =   8
  'sSD',sSliceData,...                                      %% Slice Data                                                       28 =124
  'ulChannelId',0 ...                                       %% unsigned long	 channel Id must be the last parameter            4
);                                                          %% total length: 32 * 32 Bit (128 Byte)                             128
%% MDH_H %%

%% read header information

fid = fopen(measfile,'r');
hdr_len = fread(fid, 1, 'uint32'); %%..maolin modified according to Florian.Knoll at IDEA Board 4/30/2009
disp(['Leangth of the Header : ' num2str(hdr_len)]);
fseek(fid, hdr_len, 'bof'); %%..maolin modified according to Florian.Knoll at IDEA Board 4/30/2009
[tmp Nbytes] = fread(fid,inf,'uchar'); tmp = []; 
FileSize = ftell(fid); 
disp(['Original File Size : ' num2str(FileSize)]);
disp(['Size of measurement data without header is (in bytes): ' num2str(Nbytes) '. This will take long if it is HUGE ... ...']);
%% frewind(fid);
%% fseek(fid, 32, -1); %%..maolin modified according to Florian.Knoll at IDEA Board 4/30/2009
fseek(fid, hdr_len, 'bof'); %%..maolin modified according to Florian.Knoll at IDEA Board 4/30/2009
LastADCData = 0;
lcounter = 1;

while ~feof(fid) & ~LastADCData
try
   sMDH.ulDMALength                             = fread(fid, 1, 'uint32');      % 4
   sMDH.lMeasUID                                = fread(fid, 1,  'int32');      % 8
   sMDH.ulScanCounter                           = fread(fid, 1, 'uint32');      % 12
   sMDH.ulTimeStamp                             = fread(fid, 1, 'uint32');      % 16
   sMDH.ulPMUTimeStamp                          = fread(fid, 1, 'uint32');      % 20
   for i = 1:MDH_NUMBEROFEVALINFOMASK       % 2
        sMDH.aulEvalInfoMask(i)                 = fread(fid, 1, 'uint32');      % 20 + 2 * 4 = 28
   end
   sMDH.ushSamplesInScan                        = fread(fid, 1, 'uint16');      % 30
   sMDH.ushUsedChannels                         = fread(fid, 1, 'uint16');      % 32
   K = sMDH.ushUsedChannels;
   sMDH.sLC.ushLine                             = fread(fid, 1, 'uint16');
   sMDH.sLC.ushAcquisition                      = fread(fid, 1, 'uint16');
   sMDH.sLC.ushSlice                            = fread(fid, 1, 'uint16');
   sMDH.sLC.ushPartition                        = fread(fid, 1, 'uint16');
   sMDH.sLC.ushEcho                             = fread(fid, 1, 'uint16');
   sMDH.sLC.ushPhase                            = fread(fid, 1, 'uint16');
   sMDH.sLC.ushRepetition                       = fread(fid, 1, 'uint16');
   sMDH.sLC.ushSet                              = fread(fid, 1, 'uint16');
   sMDH.sLC.ushSeg                              = fread(fid, 1, 'uint16');
   sMDH.sLC.ushIda                              = fread(fid, 1, 'uint16');
   sMDH.sLC.ushIdb                              = fread(fid, 1, 'uint16');
   sMDH.sLC.ushIdc                              = fread(fid, 1, 'uint16');
   sMDH.sLC.ushIdd                              = fread(fid, 1, 'uint16');
   sMDH.sLC.ushIde                              = fread(fid, 1, 'uint16');      % 32 + 14 * 2 = 60
   sMDH.sCutOff.ushPre                          = fread(fid, 1, 'uint16');
   sMDH.sCutOff.ushPost                         = fread(fid, 1, 'uint16');      % 60 + 2 * 2 = 64
   sMDH.ushKSpaceCentreColumn                   = fread(fid, 1, 'uint16');
   sMDH.ushDummy                                = fread(fid, 1, 'uint16');      % 64 + 2 * 2 = 68
   sMDH.fReadOutOffcentre                       = fread(fid, 1, 'float');       % 68 + 4 = 72
   sMDH.ulTimeSinceLastRF                       = fread(fid, 1, 'uint32');
   sMDH.ushKSpaceCentreLineNo                   = fread(fid, 1, 'uint16');
   sMDH.ushKSpaceCentrePartitionNo              = fread(fid, 1, 'uint16');      % 72 + 4 + 2 + 2 = 80
   for i = 1:MDH_NUMBEROFICEPROGRAMPARA    % 4
        sMDH.aushIceProgramPara(i)              = fread(fid, 1, 'uint16');      % 80 + 4 * 2 = 88
   end
   for i = 1:MDH_FREEHDRPARA  % 4
        sMDH.aushFreePara                       = fread(fid, 1, 'uint16');      % 88 + 4 * 2 = 96
   end
   sMDH.sSD.sVector.flSag                       = fread(fid, 1, 'float');
   sMDH.sSD.sVector.flCor                       = fread(fid, 1, 'float');
   sMDH.sSD.sVector.flTra                       = fread(fid, 1, 'float');       % 96 + 3 * 4 = 108
   for i = 1:4
        sMDH.aflQuaternion(i)                   = fread(fid, 1, 'float');       % 108 + 4 * 4 = 124
   end
   sMDH.ulChannelId                             = mod(fread(fid, 1, 'uint32'),K);      % 124 + 4 = 128 OK!

   if lcounter == 1
       aADC = -2*ones(1,IDX_DIM+sMDH.ushSamplesInScan*2); % the first 14 entries are the indices of this ADC in terms of sLoopCounter
       aADC(1)  = sMDH.sLC.ushLine;
       aADC(2)  = sMDH.sLC.ushAcquisition;
       aADC(3)  = sMDH.sLC.ushSlice;
       aADC(4)  = sMDH.sLC.ushPartition;
       aADC(5)  = sMDH.sLC.ushEcho;
       aADC(6)  = sMDH.sLC.ushPhase;
       aADC(7)  = sMDH.sLC.ushRepetition;
       aADC(8)  = sMDH.sLC.ushSet;
       aADC(9)  = sMDH.sLC.ushSeg;
       aADC(10) = sMDH.sLC.ushIda;
       aADC(11) = sMDH.sLC.ushIdb;
       aADC(12) = sMDH.sLC.ushIdc;
       aADC(13) = sMDH.sLC.ushIdd;
       aADC(14) = sMDH.sLC.ushIde;
       aADC(15) = rem(sMDH.ulChannelId, K);
       aADC(16) = sMDH.ushSamplesInScan;
       aADC(17) = sMDH.ushKSpaceCentreColumn;
       aADC(18) = sMDH.fReadOutOffcentre;
       aADC(19:20) = sMDH.aulEvalInfoMask;
       aADC(23) = sMDH.ushKSpaceCentreLineNo;
       aADC(24) = sMDH.sCutOff.ushPre;
       aADC(25) = sMDH.sCutOff.ushPost;
       aADC(26) = sMDH.ushKSpaceCentrePartitionNo;
       aADC(27) = sMDH.ulDMALength;
       aADC(28) = sMDH.lMeasUID;
       aADC(29) = sMDH.ulScanCounter;
       aADC(30) = sMDH.ulTimeStamp;
       % ... ...
       aADC(IDX_DIM+1:end) = fread(fid, sMDH.ushSamplesInScan*2, 'float'); % 4 bytes in each float
       LADC = 128 + sMDH.ushSamplesInScan*2*4;   % in bytes now
       NADC = floor(Nbytes/LADC)*2;
       data = zeros(NADC, size(aADC,2)); % Optimize to speed up!
       data(lcounter, :) = aADC;
   else
       aADC = -2*ones(1,IDX_DIM+sMDH.ushSamplesInScan*2); % the first 14 entries are the indices of this ADC in terms of sLoopCounter
       aADC(1)  = sMDH.sLC.ushLine;
       aADC(2)  = sMDH.sLC.ushAcquisition;
       aADC(3)  = sMDH.sLC.ushSlice;
       aADC(4)  = sMDH.sLC.ushPartition;
       aADC(5)  = sMDH.sLC.ushEcho;
       aADC(6)  = sMDH.sLC.ushPhase;
       aADC(7)  = sMDH.sLC.ushRepetition;
       aADC(8)  = sMDH.sLC.ushSet;
       aADC(9)  = sMDH.sLC.ushSeg;
       aADC(10) = sMDH.sLC.ushIda;
       aADC(11) = sMDH.sLC.ushIdb;
       aADC(12) = sMDH.sLC.ushIdc;
       aADC(13) = sMDH.sLC.ushIdd;
       aADC(14) = sMDH.sLC.ushIde;
       aADC(15) = mod(sMDH.ulChannelId, K);
       aADC(16) = sMDH.ushSamplesInScan;
       aADC(17) = sMDH.ushKSpaceCentreColumn;
       aADC(18) = sMDH.fReadOutOffcentre;
       aADC(19:20) = sMDH.aulEvalInfoMask;
       aADC(23) = sMDH.ushKSpaceCentreLineNo;
       aADC(24) = sMDH.sCutOff.ushPre;
       aADC(25) = sMDH.sCutOff.ushPost;
       aADC(26) = sMDH.ushKSpaceCentrePartitionNo;
       aADC(27) = sMDH.ulDMALength;
       aADC(28) = sMDH.lMeasUID;
       aADC(29) = sMDH.ulScanCounter;
       aADC(30) = sMDH.ulTimeStamp;
       % ... ...
       aADC(IDX_DIM+1:end) = fread(fid, sMDH.ushSamplesInScan*2, 'float');

       if size(aADC,2) == size(data,2) % readouts of the same length
           data(lcounter,:) = aADC;
       else % readouts of a different length
%           if size(aADC,2) == IDX_DIM + 32  % I regard this as a STOP sign so discard last N acquisitions of 32 (52-20) floats added by VB
%               disp(['WARNING: Discarded scan lcounter = ' num2str(lcounter) ' size(aADC,2) = ' num2str(size(aADC,2)-IDX_DIM) ' that is different from the previously detected ' num2str(size(data,2)-IDX_DIM)]);
%  		       break;
%           else % adjust the center of this readout by cutting short or padding zeros at the beginning of the readout
               data(lcounter,1:IDX_DIM) = aADC(1:IDX_DIM); % loopcounters
               delta = size(data,2)-size(aADC,2);
               if delta >=0
                  data(lcounter,delta+IDX_DIM+1:end) = aADC(IDX_DIM+1:end);
                else
                  data(lcounter,IDX_DIM+1:end) = aADC(-delta+IDX_DIM+1:end);
               end
%           end
       end

   end

   % since for some reasons, several readouts are added at the end let keep the max
   if lcounter == 1
       maxSampleInScan = sMDH.ushSamplesInScan;
   else
       maxSampleInScan = max(maxSampleInScan, sMDH.ushSamplesInScan);
   end

   curr_pos = ftell(fid);
   if curr_pos >= hdr_len + Nbytes, LastADCData = 1; end % exclude the last 384 bytes of tail
   lcounter = lcounter + 1;
catch ME
    disp(['WARNING: Incomplete Measurements Data!']);
    break; % for incomplete data
end

end % of while

fclose(fid);

% put the max maxSampleInScan in sMDH
sMDH.ushSamplesInScan = maxSampleInScan;

lcounter = lcounter - 1;

loopcounters(1:lcounter,1:15) = round(data(1:lcounter,1:15)+1); % all start from 1, unused -1
loopcounters(1:lcounter,16:IDX_DIM) = data(1:lcounter,16:IDX_DIM); % other info about the data, unused -1
dimensions = max(loopcounters);

% Let sMDH have the maximum values of all loop counters

sMDH.sLC.ushLine = dimensions(1)  ;
sMDH.sLC.ushAcquisition = dimensions(2)  ;
sMDH.sLC.ushSlice = dimensions(3)  ;
sMDH.sLC.ushPartition = dimensions(4)  ;
sMDH.sLC.ushEcho = dimensions(5)  ;
sMDH.sLC.ushPhase = dimensions(6)  ;
sMDH.sLC.ushRepetition = dimensions(7)  ;
sMDH.sLC.ushSet = dimensions(8)  ;
sMDH.sLC.ushSeg = dimensions(9)  ;
sMDH.sLC.ushIda = dimensions(10) ;
sMDH.sLC.ushIdb = dimensions(11) ;
sMDH.sLC.ushIdc = dimensions(12) ;
sMDH.sLC.ushIdd = dimensions(13) ;
sMDH.sLC.ushIde = dimensions(14) ;
sMDH.sLC.ulChannelId = dimensions(15) ;
% ... ...

% let it bear the number of channels used
sMDH.ulChannelId = dimensions(15) ;

rawdata = complex(data(1:lcounter,IDX_DIM+1:2:end),data(1:lcounter,IDX_DIM+2:2:end));

%%%%%%%%%%%%%%%%%%%%%%%% show you some information %%%%%%%%%%%%%%%%%%%%%

disp('--------------------------SUMMARY OF RAW DATA-----------------------');
disp('MiniHeader Structure :');
disp(sMDH); 
disp('Max Value of Loop Counters :');
loopcounters_max = sMDH.sLC; disp(loopcounters_max);
disp(['Max line index : ' num2str(sMDH.sLC.ushLine) ', Total number of readouts : ', num2str(size(loopcounters, 1))]);
for jj = 1:15
  eval(['lcv = ' lc_names{jj,:}(5:end) ';'])
  if lcv > 1
      if jj ==1 %% since there are often many lines, give only the several first and last lines FYI
          disp(['repeats for each line -- loopcounters(' num2str(jj) ') or (' lc_names{jj,:} ') : --']);
          for kk = 1:10
              disp(['  No. ' num2str(kk) ' = ' num2str(size(find(loopcounters(:, jj)==kk), 1)) ' repeats']);
          end 
          disp('    ... ...');
          for kk = sMDH.sLC.ushLine-10:sMDH.sLC.ushLine
              disp(['  No. ' num2str(kk) ' = ' num2str(size(find(loopcounters(:, jj)==kk), 1)) ' repeats']);
          end
      else
          disp(['lines for each -- loopcounters(' num2str(jj) ') or (' lc_names{jj,:} ') : --']);
          for kk = 1:lcv
              disp(['  No. ' num2str(kk) ' = ' num2str(size(find(loopcounters(:, jj)==kk), 1)) ' lines']);
          end
      end
  end
end

disp(['Information is retrieved using Matlab function : ']);
disp(['[rawdata,loopcounters,sMDH,lc_names] = ReadSiemensMeasVB17_idea(measfile, dispopt);']);
disp(['See loopcounters for more information. Lines in rawdata are arranged']);
disp(['in the same order as in loopcounters']);
disp(['For details in MATLAB type help ReadSiemensMeasVB17_idea']);
disp('------------------- END OF DATA SUMMARY -------------------------');

%% try to show some way to rearrange the raw data and reconstrut the images
%% if not correct, following the examples below do it yourself AND ... ...
%% YOU CAN DO IT!!

if strcmp(dispopt, 'on')

my_ver = version;
vernum = str2num(my_ver(1:3));
%% A WIP interface works only on Matlab 7.8 and above :(
if vernum >= 7.799999999999999999 & exist('ViewSiemensMeas.m','file') & exist('ViewSiemensMeas.fig','file')
  help ViewSiemensMeas;
% ViewSiemensMeas(rawdata, loopcounters, sMDH); 
  return;
end
what_to_see = 0;

for lc = 1:loopcounters_max.ushLine
  switch what_to_see
  case 1
    lin_idx = find(loopcounters(:,1)   == lc & loopcounters(:,2)   == 1 & loopcounters(:,3)   == 1 & loopcounters(:,5)   == 1  &  loopcounters(:,15)  == 1); 
  case 2
    lin_idx = find(loopcounters(:,1)  == lc & loopcounters(:, 5)  == 1); % to see Echo = 1
  case 3
    lin_idx = find(loopcounters(:,1)  == lc & loopcounters(:,15)  ==  8); % to see Ch = 8
  otherwise
    lin_idx = find(loopcounters(:,1)  == lc & loopcounters(:,3)  == 1);
  end
  araw = rawdata(lin_idx,:);
  if size(araw, 1) == 1
      raw(lc, :) = araw;
  else
      raw(lc, :) = mean(araw, 1);
  end
end

raw_size = size(raw);
aimg = fftshift(fft2(fftshift(raw)));
mag =   abs(aimg); % magnitude image
pha = angle(aimg); % phase image

atitle = ['Slice:' num2str(1)];
figure('Name', ['K-SPACE, PHASE, AND MAGNITUDE (' atitle ')']); 

subplot(3, 1, 1), imagesc(abs(raw)); axis image; colormap(gray);
subplot(3, 1, 2), imagesc(pha); axis image; colormap(gray);
subplot(3, 1, 3),  imagesc(mag); axis image; colormap(gray);


end % dispopt

