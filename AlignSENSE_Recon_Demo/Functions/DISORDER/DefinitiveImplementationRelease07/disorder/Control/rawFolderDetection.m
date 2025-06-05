function pathIn=rawFolderDetection(path,rootIn,infer)

%RAWFOLDERDETECTION   Detects the location of a given raw data path on
%pnraw
%   PATHIN=RAWFOLDERDETECTION(PATH,{ROOTIN})
%   * PATH is the relative path to raw data
%   * {ROOTIN} is the folder where the raw data has been mounted. It
%   defaults to Data/rawSource (relative to the user's home folder)
%   * INFER tries to infer the existance of a given path by using
%   validation and date
%   ** PATHIN is the absolute path to the raw data
%

%BASIC PATH TO RAW DATA
if nargin<2;rootIn=fullfile(filesep,'home','lcg13','Data','rawSource');end
if nargin<3 || isempty(infer);infer=0;end

if infer
    pathSpl=strsplit(path,'/');
    path=pathSpl{1};
end

%SEVERAL OPTIONS FOR RAW DATA PLACEMENT IN ORDER OF PRECEDENCE
pathInV{1}=fullfile(rootIn,'pnraw','raw-nnu');
pathInV{2}=fullfile(rootIn,'archive-dhcp-rawdata','archive-nnu');
pathInV{3}=fullfile(rootIn,'archive-rawdata','archive-nnu');
pathInV{4}=fullfile(rootIn,'pnraw','raw-ingenia');
pathInV{5}=fullfile(rootIn,'archive-dhcp-rawdata','archive-ingenia');
pathInV{6}=fullfile(rootIn,'archive-rawdata','archive-ingenia');

pathInV{7}=fullfile(rootIn,'raw-nnu');
pathInV{8}=fullfile(rootIn);

pathInV{9}=fullfile('/isi01','pnraw','raw-nnu');
pathInV{10}=fullfile('/isi01','archive-dhcp-rawdata','archive-nnu');
pathInV{11}=fullfile('/isi01','archive-rawdata','archive-nnu');
pathInV{12}=fullfile('/isi01','pnraw','raw-ingenia');
pathInV{13}=fullfile('/isi01','archive-dhcp-rawdata','archive-ingenia');
pathInV{14}=fullfile('/isi01','archive-rawdata','archive-ingenia');

pathInV{15}=fullfile(rootIn,'..','DISORDER/PHANTOM');
pathInV{16}=fullfile(rootIn,'..','DISORDER/VIVO');

for n=1:length(pathInV);pathInV{n}=fullfile(pathInV{n},path);end

if ~infer
    pathIn=[];
    for n=1:length(pathInV)
        if exist(pathInV{n},'dir');pathIn=pathInV{n};break;end
    end
    assert(~isempty(pathIn),'Unable to find raw data folder %s',path);
else
    pathIn=0;
    for n=1:length(pathInV)
        if exist(pathInV{n},'dir')
            a=dir(pathInV{n});
            a={a.name};
            for m=1:length(a)
                if contains(a{m},pathSpl{2})
                    pathIn=1;
                    break
                end
            end
        end
        if pathIn;break;end
    end
end
