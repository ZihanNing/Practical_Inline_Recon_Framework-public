
addpath(genpath(fileparts(mfilename('fullpath'))));


inpPath='/home/lcg13/Data/pnrawDe/ReconstructionsRelease03/INF_Priv/';
fileName{1}=sprintf('%s/dHCPRelease03.csv',inpPath);%Reconstructions spreadsheet
fileName{2}=sprintf('%s/All_Scans_with_completion_2019-02-20.xlsx',inpPath);%RecCAP spreadsheet

M=cell(1,2);
N=cell(1,2);
for n=1:2
    M{n}=readtable(fileName{n});
    N{n}=size(M{n});
end
%Extract relevant fields
M{1}=M{1}(:,[2 3 4]);
M{2}=M{2}(:,[1 2 5 10]);

%Remove neonatal rows field
neonatalRows=false(N{2}(1),1);
for m=1:N{2}(1)
    if strcmp(M{2}{m,2}{1},'Neonatal Scan');neonatalRows(m)=true;end
end
M{2}=M{2}(neonatalRows,[1 4 3]);

%Only CC numbers and dates
for s=1:N{1}(1)
    M{1}{s,1}{1}=M{1}{s,1}{1}(1:11);
    M{1}{s,3}{1}=datetime(M{1}{s,3}{1}(1:10),'InputFormat','yyyy_MM_dd');
end



%Convert to 'operative' cells
for n=1:2
    M{n}=table2cell(M{n});
    for m=1:3;A{n}{m}=M{n}(:,m);end
    for m=1:3;B{n,m}=cat(1,A{n}{m}{:});end
    
end
M=B;

contRec=1;
while 1
    NB(1)=size(B{1,1},1);
    NB(2)=size(B{2,1},1);
    if contRec>NB(1);break;end
    fo=0;
    for n=1:NB(2)        
        if all(B{1,1}(contRec,:)==B{2,1}(n,:)) && B{1,2}(contRec,:)==B{2,2}(n,:) && B{1,3}(contRec,:)==B{2,3}(n,:)
            fo=1;
            for s=1:3             
                B{1,s}(contRec,:)=[];B{2,s}(n,:)=[];              
            end
            break;
        end
    end
    if ~fo;contRec=contRec+1;end
end

n=1;
while 1
    fo=0;
    NB(1)=size(B{1,1},1);
    if n>NB(1);break;end
    if B{1,3}(n,:)>=datetime('21-Feb-2019')
        fo=1;
        fprintf('Date of %s - %d - %s not included in RedCAP spreadsheet\n',B{1,1}(n,:),B{1,2}(n,:),B{1,3}(n,:));
        for s=1:3;B{1,s}(n,:)=[];end
    end
    if ~fo;n=n+1;end
end
fprintf('----------------------------\n');

NB(2)=size(B{2,1},1);
n=1;
while 1
    fo=0;
    NB(2)=size(B{2,1},1);
    if n>NB(2);break;end
    if B{2,2}(n)==0 || isnan(B{2,2}(n))
        fo=1;
        fprintf('Not valid scan id in RedCAP for %s - %d - %s\n',B{2,1}(n,:),B{2,2}(n,:),B{2,3}(n,:));
        for s=1:3;B{2,s}(n,:)=[];end
    end
    if ~fo;n=n+1;end
end
fprintf('----------------------------\n');

NB(2)=size(B{2,1},1);
n=1;
while 1
    fo=0;
    NB(2)=size(B{2,1},1);
    if n>NB(2);break;end
    c=0;
    for m=1:size(M{2,1},1)        
        if all(M{2,1}(m,:)==B{2,1}(n,:)) && M{2,2}(m,:)==B{2,2}(n,:) && M{2,3}(m,:)==B{2,3}(n,:)
            c=c+1;
            if c==2
                fo=1;
                fprintf('Duplicated RedCAP entry for %s - %d - %s\n',B{2,1}(n,:),B{2,2}(n,:),B{2,3}(n,:));
                for s=1:3;B{2,s}(n,:)=[];end
                break;
            end
        end
    end
    if ~fo;n=n+1;end
end
fprintf('----------------------------\n');

NB(2)=size(B{2,1},1);
n=1;
while 1
    fo=0;
    NB(2)=size(B{2,1},1);
    if n>NB(2);break;end
    if ~rawFolderDetection(strcat(datestr(B{2,3}(n,:),'yyyy_mm_dd'),'/',sprintf('%d',B{2,2}(n,:))),'/home/lcg13/Data/rawSource',1)       
        for m=1:size(B{1,1},1)
            if all(B{2,1}(n,:)==B{1,1}(m,:)) && B{2,3}(n,:)==B{1,3}(m,:)
                fo=1;
                fprintf('No raw data for RedCAP entry. Most likely data as for matched subject ID and date for REDCap %s - %d - %s and Recon %s - %d - %s. Change scan id from %d to %d would solve\n',B{2,1}(n,:),B{2,2}(n,:),B{2,3}(n,:),B{1,1}(m,:),B{1,2}(m,:),B{1,3}(m,:),B{2,2}(n,:),B{1,2}(m,:));
                for s=1:3;B{2,s}(n,:)=[];end
                for s=1:3;B{1,s}(m,:)=[];end
                break;
            end
        end
        
    end
    if ~fo;n=n+1;end
end
fprintf('----------------------------\n');


NB(2)=size(B{2,1},1);
n=1;
while 1
    fo=0;
    NB(2)=size(B{2,1},1);
    if n>NB(2);break;end
    if ~rawFolderDetection(strcat(datestr(B{2,3}(n,:),'yyyy_mm_dd'),'/',sprintf('%d',B{2,2}(n,:))),'/home/lcg13/Data/rawSource',1)
        fo=1;
        fprintf('No raw data for REDCap entry: %s - %d - %s\n',B{2,1}(n,:),B{2,2}(n,:),B{2,3}(n,:));
        for s=1:3;B{2,s}(n,:)=[];end
    end
    if ~fo;n=n+1;end
end
fprintf('----------------------------\n');


NB(2)=size(B{2,1},1);
n=1;
while 1
    fo=0;
    NB(2)=size(B{2,1},1);
    if n>NB(2);break;end
    for m=1:size(M{2,1},1)
        if any(M{2,1}(m,:)~=B{2,1}(n,:)) && M{2,2}(m,:)==B{2,2}(n,:) && M{2,3}(m,:)==B{2,3}(n,:)
            fo=1;
            fprintf('Two subjects assigned to the same scan in RedCAP: %s - %d - %s and  %s - %d - %s.',B{2,1}(n,:),B{2,2}(n,:),B{2,3}(n,:),M{2,1}(m,:),M{2,2}(m,:),M{2,3}(m,:));                
            forec=0;
            for l=1:size(M{1,1},1)
                if all(M{1,1}(l,:)==B{2,1}(n,:)) && M{1,3}(l,:)==M{2,3}(m,:)
                    forec=1;
                    fprintf(' Existing match %s - %d - %s in the reconstructed data.',M{1,1}(l,:),M{1,2}(l,:),M{1,3}(l,:));
                    for s=1:size(B{1,1},1)
                        if all(B{1,1}(s,:)==M{1,1}(l,:)) && B{1,2}(s,:)==M{1,2}(l,:) && B{1,3}(s,:)==M{1,3}(l,:)
                            for r=1:3;B{1,r}(s,:)=[];end
                            break
                        end
                    end
                    break
                end
            end                                
            if forec;fprintf(' Change scan id from %d to %d would solve\n',B{2,2}(n,:),M{1,2}(l,:));else fprintf(' Uncertain. Discussion required!\n');end
            for s=1:3;B{2,s}(n,:)=[];end
            break
        end
    end                                                       
    if ~fo;n=n+1;end
end
fprintf('----------------------------\n');


NB(1)=size(B{1,1},1);
for n=1:NB(1)
    fprintf('Case %s - %d - %s not found in RedCAP\n',B{1,1}(n,:),B{1,2}(n,:),B{1,3}(n,:));
end
fprintf('----------------------------\n');

NB(2)=size(B{2,1},1);
for n=1:NB(2)
    fprintf('Case %s - %d - %s not found in Recons\n',B{2,1}(n,:),B{2,2}(n,:),B{2,3}(n,:));
end

