function Par=feedHeader(Par,check)

%FEEDHEADER   Feeds the header of a given rec object with information
%directly read from the .par files
%   PAR=SORTDATA(PAR,NAMES,{CHECK})
%   * PAR are the parameters of a reconstruction structure
%   * {CHECK} is a flag to perform consistency checks. It defaults to 0
%   (replacing of information without consistency check)
%   * PAR are the parameters of a reconstruction structure
%

if nargin<2 || isempty(check);check=0;end

%DATA TYPES INFORMATION
RFIndexFields={'typ','mix','dyn','card','echo','loca', ...
               'chan','extr1','extr2','ky','kz','na',...
               'aver','sign','rf','grad','enc','rtop',...
               'rr','size','offset','random_phase','meas_phase','pda_index',...
               'pda_fac','dyn_time','coded_size','chan_grp','format','ky_label',...
                'kz_label'};
RFIndexTypes={'uint8','uint16','uint16','uint16','uint16','uint16',...
              'uint16','uint16','uint16','int32','int32','uint16',...
              'uint16','int8','uint16','uint16','int16','uint16',...
              'uint16','uint64','uint64','uint16','uint16','uint16',...
              'double','uint32','uint32','uint32','uint8','int32',...
              'int32'};
RFIndexMaps={[],'mix_nr','dynamic_scan_nr','cardiac_phase_nr','echo_nr','location', ...
             [],'row_nr','extra_attr_nr',[],[],'e3_profile_nr',...
             'measurement_nr',[],'rf_echo_nr','grad_echo_nr',[],'rtop_offset',...
             'rr_interval',[],[],'rand_phase',[],[],...
             [],[],'coded_data_size','channels_active',[],'e1_profile_nr',...
             'e2_profile_nr'};
ParFileFields={'data_size','coded_data_size','src_code','dst_code','seq_nr','label_type',...
               'control','monitoring_flag','measurement_phase','measurement_sign','gain_setting_index','spare_1',...
               'spare_2','progress_cnt','mix_nr','dynamic_scan_nr','cardiac_phase_nr','echo_nr',...
               'location','row_nr','extra_attr_nr','measurement_nr','e1_profile_nr','e2_profile_nr',...
               'e3_profile_nr','rf_echo_nr','grad_echo_nr','enc_time','rand_phase','rr_interval',...
               'rtop_offset','channels_active'};
ParFileTypes={'uint32','uint32','uint16','uint16','uint16','uint16',...
              'uint8','uint8','uint8','uint8','uint8','uint8',...
              'uint16','uint16','uint16','uint16','uint16','uint16',...
              'uint16','uint16','uint16','uint16','uint16','uint16',...
              'uint16','uint16','uint16','uint16','uint16','uint16',...
              'uint16','uint32'};
ParFileLocations={1:4,5:8,9:10,11:12,13:14,15:16,...
                  17,18,19,20,21,22,...
                  23:24,25:26,27:28,29:30,31:32,33:34,...
                  35:36,37:38,39:40,41:42,43:44,45:46,...
                  47:48,49:50,51:52,53:54,55:56,57:58,...
                  59:60,61:64};
                           
%READ LAB
s=dir(Par.Filename.Parameter);
the_size=s.bytes;assert(mod(the_size/64,1)==0,'Not forming blocks of 64 bytes');
fidlab=fopen(Par.Filename.Parameter,'r','ieee-le');
fileContent=fread(fidlab,inf,'uint8=>uint8');
fclose(fidlab);
fileContent=reshape(fileContent,64,[]);%64 bytes
fileContentS=sum(abs(fileContent),1);
endContent=find(fileContentS~=0,1,'last');
%fprintf('Number of entries before/after removing zeros: %d/%d\n',size(fileContent,2),endContent);
fileContent=fileContent(:,1:endContent);  

%ASSIGN LAB FIELDS
for t=1:length(ParFileFields)
    ParFile.(ParFileFields{t})=fileContent(ParFileLocations{t},:);
    ParFile.(ParFileFields{t})=typecast(ParFile.(ParFileFields{t})(:),ParFileTypes{t});
end

if ~check && ~any(ParFile.control==19);return;end%If there is no navigator data we do nothing

%GET THE CHANNELS (AS A BINARY ARRAY)
channels=single(de2bi(ParFile.channels_active,32));     

%CALCULATE THE OFFSETS
offsets=zeros(size(ParFile.data_size),'uint64');
offsets(1)=512;
offsets(2:end)=offsets(1)+uint64(cumsum(ParFile.data_size(1:end-1)));
%%This corresponds to the first half and second half of coded_data_size, all zeros in the examples we have seen
%%info.labels.LeadingDummies.vals   = bitshift (bitand(unparsed_labels(2,:), (2^16-1)),  -0);
%%info.labels.TrailingDummies.vals  = bitshift (bitand(unparsed_labels(2,:), (2^32-1)), -16);
%%for k=2:info.nLabels,
%%    info.fseek_offsets(k) = info.fseek_offsets(k-1)+ info.labels.DataSize.vals(k-1) - info.labels.TrailingDummies.vals(k-1)  - info.labels.LeadingDummies.vals(k-1);
%%end

%INDEXES OF INTERPRETABLE DATA SIZE
if check;indIndex=find(ParFile.control==0 | ParFile.control==3 | (ParFile.control==2 & ParFile.dst_code==9) | ParFile.control==16);%Respectively for data, EPI ghosting, Phase, and Noise
else indIndex=find(ParFile.control==0 | ParFile.control==3 | (ParFile.control==2 & ParFile.dst_code==9) | ParFile.control==19 | ParFile.control==16);%Respectively for data, EPI ghosting, Phase, Navigator, and Noise
end
channelsRF=channels(indIndex,:);%Channels of interest
NRF=sum(channelsRF(:));%Total number of profiles
NRep=sum(channelsRF,2);%Number of channel measures per profile
controlRF=ParFile.control(indIndex);%Types of interest   
addV=cell(1,length(NRep));%To add sizes to the offset
for r=1:length(NRep);addV{r}=(0:NRep(r)-1)';end
addV=uint64(cat(1,addV{:}));           
for t=1:length(RFIndexFields);IndexNew.(RFIndexFields{t})=zeros(NRF,1,RFIndexTypes{t});end%Booking memory    
controlRF=repelem(controlRF,NRep);%Replicated types    

%JUNK METADATA. READ THE RAW FOR KYMIN
indMeta=offsets(ParFile.control==2 & ParFile.dst_code==1);
sizMeta=ParFile.data_size(ParFile.control==2 & ParFile.dst_code==1);
metaData=[];
fidraw=fopen(Par.Filename.Data,'r','ieee-le');    
for l=1:length(indMeta)
    fseek(fidraw,indMeta(l),-1);    
    metaData=cat(1,metaData,fread(fidraw,sizMeta(l)/2,'single=>single'));        
end 
fclose(fidraw);
if isempty(metaData);return;end

%SOME INTERPRETABLE LABELS---NOT FULLY EXPLOITED YET
LabelsNew=interpretLabels(metaData,Par.Labels);

%MAPPING TO RF
Index=Par.Labels.Index;
for t=1:length(RFIndexMaps)
    if ~isempty(RFIndexMaps{t});IndexNew.(RFIndexFields{t})=repelem(ParFile.(RFIndexMaps{t})(indIndex),NRep);end
end
if check;IndexNew.typ(controlRF==0)=1;IndexNew.typ(controlRF==3)=3;IndexNew.typ(controlRF==2)=4;IndexNew.typ(controlRF==16)=5;    
else IndexNew.typ(controlRF==0)=1;IndexNew.typ(controlRF==3)=3;IndexNew.typ(controlRF==19)=4;IndexNew.typ(controlRF==2)=2;IndexNew.typ(controlRF==16)=5;
end
IndexNew.chan=mod(find(channelsRF')-1,32);
IndexNew.ky=repelem(int32(ParFile.e1_profile_nr(indIndex))+int32(LabelsNew.Kmin(ParFile.echo_nr(indIndex)+1,2)),NRep);
IndexNew.kz=repelem(int32(ParFile.e2_profile_nr(indIndex))+int32(LabelsNew.Kmin(ParFile.echo_nr(indIndex)+1,3)),NRep);
IndexNew.sign=repelem(1-2*int8(ParFile.measurement_sign(indIndex)),NRep);
IndexNew.enc=int16(IndexNew.ky);
IndexNew.enc(IndexNew.typ==3)=IndexNew.enc(IndexNew.typ==3)+int16(IndexNew.grad(IndexNew.typ==3))+int16(LabelsNew.Kmin(IndexNew.echo(IndexNew.typ==3)+1,2));
%IndexNew.enc=repelem(int16(ParFile.grad_echo_nr(indIndex)),NRep);
IndexNew.size=uint64(repelem(ParFile.data_size(indIndex)./uint32(NRep),NRep));
IndexNew.coded_size=uint32(int64(IndexNew.size)-int64(IndexNew.coded_size));
IndexNew.offset=repelem(offsets(indIndex),NRep)+addV.*IndexNew.size;
IndexNew.meas_phase=repelem(uint16(ParFile.measurement_phase(indIndex)),NRep);
IndexNew.pda_index=repelem(uint16(ParFile.gain_setting_index(indIndex)),NRep);
%PROBLEMS WITH PDA CHANGING FROM CHANNEL TO CHANNEL
gainChannel=Index.pda_fac(Index.typ==5)./Par.Labels.PDAFactors(Index.pda_index(Index.typ==5)+1).';
chanUn=unique(IndexNew.chan)+1;
NC=length(chanUn);
gainPerChannel=ones(1,max(chanUn));
for c=1:NC;gainPerChannel(chanUn(c))=mean(gainChannel(Index.chan(Index.typ==5)==chanUn(c)-1));end
IndexNew.pda_fac=(gainPerChannel(IndexNew.chan+1).*Par.Labels.PDAFactors(IndexNew.pda_index+1)).';%WE SHOULD READ THE RAW FOR PDA BUT IT IS VARIABLE LOCATION
IndexNew.dyn_time=repelem(uint32(ParFile.enc_time(indIndex)),NRep);
IndexNew.format(:)=0;%Not sure about this!!

%CHECKS FOR CONSISTENCY
if check    
    indCheck=ismember(Index.typ,[1 3 5]);
    indCheckNew=ismember(IndexNew.typ,[1 3 5]);
    for t=1:length(RFIndexFields)
        assert(all(size(IndexNew.(RFIndexFields{t})(indCheckNew))==size(Index.(RFIndexFields{t})(indCheck))),'Inconsistent sizes for field %s in file %s (new: %d / old: %d)',RFIndexFields{t},Par.Filename.Parameter,length(IndexNew.(RFIndexFields{t})(indCheckNew)),length(Index.(RFIndexFields{t})(indCheck)));
        assert(all(abs(single(IndexNew.(RFIndexFields{t})(indCheckNew))-single(Index.(RFIndexFields{t})(indCheck))))<1e-6,'Inconsistent results for field %s in file %s',RFIndexFields{t},Par.Filename.Parameter);
    end
    fprintf('Consistent feeding for %s\n',Par.Filename.Parameter)
else
    Par.Labels.Index=IndexNew;
end

%SOME INFORMATION FROM ISMRM DATA FORMAT
% 
% #ifndef PHILIPS_H
% #define PHILIPS_H
% /**********************************************************************************
%  Header file describing the layout of the Philips MRI raw data label files (*.lab)
%  The information in this header has been obtained by reverse engineering
%  using a hex editor and the output readPhilipsExports module and associated
%  ReadPhilips node in the GPI framework. Those tools can be obtained from:
%  https://github.com/gpilab/philips-data-reader
%  Michael S. Hansen (michael.hansen@nih.gov)
%  April 2015
% ***********************************************************************************/
% namespace philips
% {
% 
%   /*
%     
%     #Reverse engineered from GPI module readPhilipsExports
%     #In Python:
%     import readPhilipsExports as rp   
%     for l in range(100000):
%         if rp.LabelTypeEnum(l) != 'LABEL_TYPE_MAX':
%             print str(rp.LabelTypeEnum(l)) + " = " + str(hex(l))
%    */  
%   typedef enum
%   {
%     LABEL_TYPE_MIN = 0x7f00,
%     LABEL_TYPE_STANDARD = 0x7f01,
%     LABEL_TYPE_IGEO = 0x7f02,
%     LABEL_TYPE_IGEO_PMC = 0x7f03,
%     LABEL_TYPE_COIL_POS = 0x7f04,
%     LABEL_TYPE_NON_LIN = 0x7f05,
%     LABEL_TYPE_MAX
%   } label_types;
% 
% 
%   /*
%     #Reverse engineered from GPI module readPhilipsExports
%     #In Python:
%     import readPhilipsExports as rp   
%     for l in range(100000):
%         if rp.CtrlEnum(l) != 'CTRL_MAX':
%             print str(rp.CtrlEnum(l)) + " = " + str(hex(l)) + "," 
%    */
%   typedef enum
%   {
%     CTRL_NORMAL_DATA = 0x0,
%     CTRL_DC_OFFSET_DATA = 0x1,
%     CTRL_JUNK_DATA = 0x2,
%     CTRL_ECHO_PHASE_DATA = 0x3,
%     CTRL_NO_DATA = 0x4,
%     CTRL_NEXT_PHASE = 0x5,
%     CTRL_SUSPEND = 0x6,
%     CTRL_RESUME = 0x7,
%     CTRL_TOTAL_END = 0x8,
%     CTRL_INVALIDATION = 0x9,
%     CTRL_TYPE_NR_END = 0xa,
%     CTRL_VALIDATION = 0xb,
%     CTRL_NO_OPERATION = 0xc,
%     CTRL_DYN_SCAN_INFO = 0xd,
%     CTRL_SELECTIVE_END = 0xe,
%     CTRL_FRC_CH_DATA = 0xf,
%     CTRL_FRC_NOISE_DATA = 0x10,
%     CTRL_REFERENCE_DATA = 0x11,
%     CTRL_DC_FIXED_DATA = 0x12,
%     CTRL_NAVIGATOR_DATA = 0x13,
%     CTRL_FLUSH = 0x14,
%     CTRL_RECON_END = 0x15,
%     CTRL_IMAGE_STATUS = 0x16,
%     CTRL_TRACKING = 0x17,
%     CTRL_FLUOROSCOPY_TOGGLE = 0x18,
%     CTRL_REJECTED_DATA = 0x19,
%     CTRL_PROGRESS_INFO = 0x1a,
%     CTRL_END_PREP_PHASE = 0x1b,
%     CTRL_CHANNEL_DEFINITION = 0x1c,
%     CTRL_MAX
%   } ctrl_types;
% 
% 
%   /*
%     
%     The structure of the label (*.lab) files has been reverse engineered in the following way:
%     With access to two different versions of the files (an older and a newer) it was determined that there are two different formats. 
%     
%     The readPhilipsExports Python Module from GPI was used find a) the names of the header fields and b) the location of the header fields. 
%     The readPhilipsExports module has a function:
%     readLab(filename, isMira)
%        Parse a *.lab header accompanying the *.raw data file.
%     
%     The isMira argument is apparently used to indicate if this is a new header (True) or an old one (False).
%     Both new and old headers are 64 bytes long. 
%     By visual inspection of the result of reading new and old files using this function, the location of two header fields
%     "data_size" and "label_type" was found. These fields are located in the same bytes in the new and old header. 
%     A file with test labels was then generated using the following function:
%     int generate_test_lab(const char* filename)
%     {
%       typedef struct
%       {
%         uint32_t data_size;
% 	char test_bank1[10];
% 	uint16_t label_type;
% 	char test_bank2[48];
%       } test_label;
%   
%       test_label l;
%   
%       std::ofstream f(filename, std::ofstream::binary);
%       
%       if (!f.is_open()) {
%         std::cout << "Error opening file" << std::endl;
%         return -1;
%       }
%       for (size_t i = 0; i < 10; ++i) {
%         memset(&l,0,sizeof(test_label));
% 	l.label_type = philips::LABEL_TYPE_STANDARD;
% 	l.test_bank1[i] = 1;
% 	f.write(reinterpret_cast<char*>(&l),sizeof(test_label));
%       }
%       for (size_t i = 0; i < 48; ++i) {
%         memset(&l,0,sizeof(test_label));
% 	l.label_type = philips::LABEL_TYPE_STANDARD;
% 	l.test_bank2[i] = 1;
% 	f.write(reinterpret_cast<char*>(&l),sizeof(test_label));
%       }
%       f.close();
%       return 0;
%     }
%     The test label file was then read in python with code like:
%     #Python script for inspecting test labels
%     import readPhilipsExports as rp
%     lab = rp.rp.readLab('test.lab', True) #For new label layout
%     lab = rp.rp.readLab('test.lab', False) #For old label layout
%     #End of Python script for inspecting test labels
%     For example to figure out which bytes contribute to the echo_nr field, simply inspect:
%     
%     lab['echo_nr'] 
%     and find non-zero entries. 
%     The names of the header fields are listed below.
%     Old header field:
%     ['control', 'grad_echo_nr', 'rr_interval', 'echo_nr', 'cardiac_phase_nr', 'dst_code', 'leading_dummies_size', 'monitoring_flag', 'e2_profile_nr', 'row_nr', 'trailing_dummies_size', 'label_type', 'raw_format', 'src_code', 'rf_echo_nr', 'rtop_offset', 'dynamic_scan_nr', 'random_phase', 'location_nr', 'data_size', 'gain_setting_index', 'e3_profile_nr', 'progress_cnt', 'measurement_sign', 'channels_active', 'enc_time', 'measurement_phase', 'extra_attr_nr', 'e1_profile_nr', 'measurement_nr', 'mix_nr', 'seq_nr', 'spare_1']
%     New header fields:
%     ['control', 'grad_echo_nr', 'rr_interval', 'echo_nr', 'cardiac_phase_nr', 'measurement_nr', 'monitoring_flag', 'e2_profile_nr', 'row_nr', 'label_type', 'e1_profile_nr', 'rf_echo_nr', 'rtop_offset', 'dynamic_scan_nr', 'random_phase', 'location_nr', 'data_size', 'gain_setting_index', 'e3_profile_nr', 'progress_cnt', 'measurement_sign', 'coded_data_size', 'channels_active', 'enc_time', 'measurement_phase', 'extra_attr_nr', 'normalization_factor', 'raw_format', 'mix_nr', 'seq_nr', 'spare_1']
%    */
% 
%   typedef struct
%   {
%     uint32_t data_size;
%     uint32_t coded_data_size;
%     uint16_t src_code;
%     uint16_t dst_code;
%     uint16_t seq_nr;
%     uint16_t label_type;
%     char control;
%     char monitoring_flag;
%     char measurement_phase;
%     char measurement_sign;
%     char gain_setting_index;
%     char spare_1;
%     uint16_t spare_2;
%     uint16_t progress_cnt;
%     uint16_t mix_nr;
%     uint16_t dynamic_scan_nr;
%     uint16_t cardiac_phase_nr;
%     uint16_t echo_nr;
%     uint16_t location_nr;
%     uint16_t row_nr;
%     uint16_t extra_attr_nr;
%     uint16_t  measurement_nr;
%     uint16_t e1_profile_nr;
%     uint16_t e2_profile_nr;
%     uint16_t e3_profile_nr;
%     uint16_t rf_echo_nr;
%     uint16_t grad_echo_nr;
%     uint16_t enc_time;
%     uint16_t random_phase;
%     uint16_t rr_interval;
%     uint16_t rtop_offset;
%     uint32_t channels_active;
%   } old_label;
% 
%   typedef struct
%   {
%     uint32_t data_size;
%     uint32_t coded_data_size;
%     float normalization_factor;
%     uint16_t seq_nr;
%     uint16_t label_type;
%     char control;
%     char monitoring_flag;
%     char measurement_phase;
%     char measurement_sign;
%     char gain_setting_index;
%     char raw_format;
%     uint16_t spare_1;
%     uint16_t progress_cnt;
%     uint16_t mix_nr;
%     uint16_t dynamic_scan_nr;
%     uint16_t cardiac_phase_nr;
%     uint16_t echo_nr;
%     uint16_t location_nr;
%     uint16_t row_nr;
%     uint16_t extra_attr_nr;
%     uint16_t measurement_nr;
%     uint16_t e1_profile_nr;
%     uint16_t e2_profile_nr;
%     uint16_t e3_profile_nr;
%     uint16_t rf_echo_nr;
%     uint16_t grad_echo_nr;
%     uint16_t enc_time;
%     uint16_t random_phase;
%     uint16_t rr_interval;
%     uint16_t rtop_offset;
%     uint32_t channels_active;
%   } new_label; //Total size 64 bytes
% 
%   typedef struct
%   {
%     union
%     {
%       old_label old_;
%       new_label new_;
%     }; 
%   } label; //Total size 64 bytes
% 
% }
% 
% #endif //PHILIPS_H


%SOME FURTHER INFORMATION
% %% LOADLABRAW     Load a Philips LABRAW file
% %
% % [DATA,INFO] = LOADLABRAW(FILENAME)
% %
% %   FILENAME is a string containing a file prefix or name of the LAB
% %   hexadecimal label file or RAW data file, e.g. RAW_001 or RAW_001.LAB or RAW_001.RAW
% %
% %   DATA is an N-dimensional array holding the raw k-space data.
% %
% %   INFO is a structure containing details from the LAB hexadecimal label file
% %
% % [DATA,INFO] = LOADLABRAW([])
% %
% %   When the passed FILENAME is not provided or is an empty array or empty 
% %   string.  The user chooses a file using UIGETFILE.
% %
% % [DATA,INFO] = LOADLABRAW(FILENAME,'OptionName1',OptionValue1,...)
% %
% %   Options can be passed to LOADLABRAW to control the range/pattern of
% %   loaded data, verbose output, etc.  The list below shows the avialable 
% %   options.  Names are case-sensitive
% %
% %       OptionName          OptionValue       Description    
% %       ----------          -----------     ---------------
% %       'coil'              numeric         coils  
% %       'kx'                numeric         k-space kx samples
% %       'ky'                numeric         k-space ky rows (E1)
% %       'kz'                numeric         k-space kz rows (E2)
% %       'e3'                numeric         k-space 3rd encoding dim
% %       'loc'               numeric         locations
% %       'ec'                numeric         echoes          
% %       'dyn'               numeric         dynamics        
% %       'ph'                numeric         cardiac phases  
% %       'row'               numeric         rows  
% %       'mix'               numeric         mixes  
% %       'avg'               numeric         averages  
% %       'verbose'           logical         [ true |{false}]
% %       'savememory'        logical         [{true}| false ]
% %
% %       When 'savememory' is true, SINGLE precision is used instead of DOUBLE
% %
% %   Example:
% %       myfile = 'example.lab';
% %       [data,info] = loadLabRaw(myfile,'coil',[1 5],'verbose',true);
% %
% % [DATA,INFO] = LOADLABRAW(FILENAME,LOADOPTS)
% %
% %   LOADOPTS is a structure with fieldnames equal to any of the possible
% %   OptionNames.
% %
% %   Example:
% %       loadopts.coil = [1 5];
% %       loadopts.verbose = true;
% %       [data,info] = loadLabRaw(myfile,loadopts);
% %
% %   For any dimension, values may be repeated and appear in any order.
% %   Values that do not intersect with the available values for that
% %   dimension will be ignored.  If the intersection of the user-defined
% %   dimension values and the available dimension range has length zero, an
% %   error is generated.  The order of the user-defined pattern is preserved.
% %
% %   Example:
% %       % load a specific pattern of locations (-1 will be ignored)
% %       loadopts.loc = [1 1 2 1 1 2 -1];
% %       [data,info] = loadLabRaw(myfile,loadopts);
% %
% % INFO = LOADLABRAW(FILENAME)
% %
% %   If only one return argument is provided, the INFO structure will be
% %   returned.  DATA will not be loaded (fast execution).
% %
% % INFO structure contents
% %
% %   The INFO structure contains all the information from the LAB file in
% %   a filed names LABELS as well as other useful information to describe 
% %   and to work with the loaded DATA array.  The list below describes some
% %   of the additional fields found within INFO
% %
% %   FieldName              Description
% %   ---------              ------------------------------------------------
% %   FILENAME               filename of the loaded data
% %   LOADOPTS               structure containing the load options (see above)
% %   DIMS                   structure containing the DATA dimension names and values
% %   LABELS                 structure containing label names and values
% %   LABELS_ROW_INDEX_ARRAY (see below)
% %   LABEL_FIELDNAMES       names used for the labels
% %   IDX                    structure of arrays of index of different label types
% %   FSEEK_OFFSETS          byte offsets in the .RAW file for each data vector 
% %   NLABELS                # of total labels avialable in the LAB file
% %   NLOADEDLABELS          # of labels loaded from the LABRAW file
% %   NDATALABELS            # of labels in the returned data array (may contain repeats)
% %   DATASIZE               array showing the size of the returned DATA array
% %   FRC_NOISE_DATA         array of the FRC noise data
% %
% %   The INFO.LABELS_ROW_INDEX_ARRAY is a special array that is the same 
% %   size as the DATA array (minus the first two dimensions used to store 
% %   COIL and KX).  A given index for a given raw data vector in the DATA 
% %   array will return the label index number describing the details of that 
% %   raw data vector in the INFO.LABELS array when that same index is used 
% %   with the INFO.LABELS_ROW_INDEX_ARRAY.  This provides a quick way to 
% %   recall the details of any individual raw data vector contained within DATA.  
% %   If the INFO.TABLE_ROW_INDEX_ARRAY holds a ZERO for a given index, there 
% %   was no label from the LAB file that matched the dimension location in DATA.  
% %
% %  See also: LOADPARREC, LOADXMLREC, LOADDICOM
% %
% %  Dependencies: none
% %
% 
% %% Revision History
% % * 2008.11.07    initial version - brianwelch
% 
% %% Function definition
% function [data,info] = loadLabRaw(filename,varargin)
% 
% %% Start execution time clock and initialize DATA and INFO to empty arrays
% tic;
% data=[];
% info=[];
% 
% %% Initialize INFO structure
% % Serves to fix the display order
% info.filename = [];
% info.loadopts = [];
% info.dims = [];
% info.labels = [];
% info.labels_row_index_array = [];
% info.label_fieldnames = [];
% info.idx = [];
% info.fseek_offsets = [];
% info.nLabels = [];
% info.nLoadedLabels = [];
% info.nDataLabels = [];
% info.nNormalDataLabels = [];
% info.datasize = [];
% 
% %% Allow user to select a file if input FILENAME is not provided or is empty
% if nargin<1 | length(filename)==0,
%     [fn, pn] = uigetfile({'*.raw'},'Select a RAW file');
%     if fn~=0,
%         filename = sprintf('%s%s',pn,fn);
%     else
%         disp('LOADLABRAW cancelled');
%         return;
%     end
% end
% 
% %% Parse the filename.
% % It may be the LAB filename, RAW filename or just the filename prefix
% % Instead of REGEXP, use REGEXPI which igores case
% toks = regexpi(filename,'^(.*?)(\.lab|\.raw)?$','tokens');
% prefix = toks{1}{1};
% labname = sprintf('%s.lab',prefix);
% rawname = sprintf('%s.raw',prefix);
% info.filename = filename;
% 
% %% Open LAB file and read all hexadecimal labels
% labfid = fopen(labname,'r');
% if labfid==-1,
%     error( sprintf('Cannot open %s for reading', labname) );
% end
% 
% %% Read all hexadecimal labels
% [unparsed_labels, readsize] = fread (labfid,[16 Inf], 'uint32=>uint32');
% info.nLabels = size(unparsed_labels,2);
% fclose(labfid);
% 
% %% Parse hexadecimal labels
% % Inspired by Holger Eggers' readRaw.m.  Thanks Holger! 
% % See arsrcglo1.h for more details.
% info.labels.DataSize.vals         = unparsed_labels(1,:);
% 
% info.labels.LeadingDummies.vals   = bitshift (bitand(unparsed_labels(2,:), (2^16-1)),  -0);
% info.labels.TrailingDummies.vals  = bitshift (bitand(unparsed_labels(2,:), (2^32-1)), -16);
% 
% info.labels.SrcCode.vals          = bitshift (bitand(unparsed_labels(3,:), (2^16-1)),  -0);
% info.labels.DstCode.vals          = bitshift (bitand(unparsed_labels(3,:), (2^32-1)), -16);
% 
% info.labels.SeqNum.vals           = bitshift (bitand(unparsed_labels(4,:), (2^16-1)),  -0);
% info.labels.LabelType.vals        = bitshift (bitand(unparsed_labels(4,:), (2^32-1)), -16);
% 
% info.labels.ControlType.vals      = bitshift( bitand(unparsed_labels(5,:),  (2^8-1)),  -0);
% info.labels.MonitoringFlag.vals   = bitshift( bitand(unparsed_labels(5,:), (2^16-1)),  -8);
% info.labels.MeasurementPhase.vals = bitshift( bitand(unparsed_labels(5,:), (2^24-1)), -16);
% info.labels.MeasurementSign.vals  = bitshift( bitand(unparsed_labels(5,:), (2^32-1)), -24);
% 
% info.labels.GainSetting.vals      = bitshift( bitand(unparsed_labels(6,:),  (2^8-1)),  -0);
% info.labels.Spare1.vals           = bitshift( bitand(unparsed_labels(6,:), (2^16-1)),  -8);
% info.labels.Spare2.vals           = bitshift (bitand(unparsed_labels(6,:), (2^32-1)), -16);
% 
% info.labels.ProgressCnt.vals      = bitshift (bitand(unparsed_labels(7,:), (2^16-1)),  -0);
% info.labels.Mix.vals              = bitshift (bitand(unparsed_labels(7,:), (2^32-1)), -16);
% 
% info.labels.Dynamic.vals          = bitshift (bitand(unparsed_labels(8,:), (2^16-1)),  -0);
% info.labels.CardiacPhase.vals     = bitshift (bitand(unparsed_labels(8,:), (2^32-1)), -16);
% 
% info.labels.Echo.vals             = bitshift (bitand(unparsed_labels(9,:), (2^16-1)),  -0);
% info.labels.Location.vals         = bitshift (bitand(unparsed_labels(9,:), (2^32-1)), -16);
% 
% info.labels.Row.vals              = bitshift (bitand(unparsed_labels(10,:), (2^16-1)),  -0);
% info.labels.ExtraAtrr.vals        = bitshift (bitand(unparsed_labels(10,:), (2^32-1)), -16);
% 
% info.labels.Measurement.vals      = bitshift (bitand(unparsed_labels(11,:), (2^16-1)),  -0);
% info.labels.E1.vals               = bitshift (bitand(unparsed_labels(11,:), (2^32-1)), -16);
% 
% info.labels.E2.vals               = bitshift (bitand(unparsed_labels(12,:), (2^16-1)),  -0);
% info.labels.E3.vals               = bitshift (bitand(unparsed_labels(12,:), (2^32-1)), -16);
% 
% info.labels.RfEcho.vals           = bitshift (bitand(unparsed_labels(13,:), (2^16-1)),  -0);
% info.labels.GradEcho.vals         = bitshift (bitand(unparsed_labels(13,:), (2^32-1)), -16);
% 
% info.labels.EncTime.vals          = bitshift (bitand(unparsed_labels(14,:), (2^16-1)),  -0);
% info.labels.RandomPhase.vals      = bitshift (bitand(unparsed_labels(14,:), (2^32-1)), -16);
% 
% info.labels.RRInterval.vals       = bitshift (bitand(unparsed_labels(15,:), (2^16-1)),  -0);
% info.labels.RTopOffset.vals       = bitshift (bitand(unparsed_labels(15,:), (2^32-1)), -16);
% 
% info.labels.ChannelsActive.vals   = unparsed_labels(16,:);
% 
% clear unparsed_labels;
% 
% %% Find unique values of each label field
% info.label_fieldnames = fieldnames(info.labels);
% for k=1:length(info.label_fieldnames),
%     info.labels.(info.label_fieldnames{k}).uniq = unique( info.labels.(info.label_fieldnames{k}).vals ); 
% end
% 
% %% Calculate fseek offsets
% info.fseek_offsets = zeros(info.nLabels,1);
% info.fseek_offsets(1)=512; % add mysterious 512 byte offset to begin reading file
% for k=2:info.nLabels,
%     info.fseek_offsets(k) = info.fseek_offsets(k-1)+ info.labels.DataSize.vals(k-1) - info.labels.TrailingDummies.vals(k-1)  - info.labels.LeadingDummies.vals(k-1);
% end
% info.idx.no_data = find(info.labels.DataSize.vals==0);
% info.fseek_offsets(info.idx.no_data) = -1;
% 
% %% Find indices of different label control types
% % See arsrcglo1.h for more details.
% standard_labels = info.labels.LabelType.vals==32513;
% info.idx.NORMAL_DATA         = find(info.labels.ControlType.vals== 0 & standard_labels);
% info.idx.DC_OFFSET_DATA      = find(info.labels.ControlType.vals== 1 & standard_labels);
% info.idx.JUNK_DATA           = find(info.labels.ControlType.vals== 2 & standard_labels);
% info.idx.ECHO_PHASE_DATA     = find(info.labels.ControlType.vals== 3 & standard_labels);
% info.idx.NO_DATA             = find(info.labels.ControlType.vals== 4 & standard_labels);
% info.idx.NEXT_PHASE          = find(info.labels.ControlType.vals== 5 & standard_labels);
% info.idx.SUSPEND             = find(info.labels.ControlType.vals== 6 & standard_labels);
% info.idx.RESUME              = find(info.labels.ControlType.vals== 7 & standard_labels);
% info.idx.TOTAL_END           = find(info.labels.ControlType.vals== 8 & standard_labels);
% info.idx.INVALIDATION        = find(info.labels.ControlType.vals== 9 & standard_labels);
% info.idx.TYPE_NR_END         = find(info.labels.ControlType.vals==10 & standard_labels);
% info.idx.VALIDATION          = find(info.labels.ControlType.vals==11 & standard_labels);
% info.idx.NO_OPERATION        = find(info.labels.ControlType.vals==12 & standard_labels);
% info.idx.DYN_SCAN_INFO       = find(info.labels.ControlType.vals==13 & standard_labels);
% info.idx.SELECTIVE_END       = find(info.labels.ControlType.vals==14 & standard_labels);
% info.idx.FRC_CH_DATA         = find(info.labels.ControlType.vals==15 & standard_labels);
% info.idx.FRC_NOISE_DATA      = find(info.labels.ControlType.vals==16 & standard_labels);
% info.idx.REFERENCE_DATA      = find(info.labels.ControlType.vals==17 & standard_labels);
% info.idx.DC_FIXED_DATA       = find(info.labels.ControlType.vals==18 & standard_labels);
% info.idx.DNAVIGATOR_DATA     = find(info.labels.ControlType.vals==19 & standard_labels);
% info.idx.FLUSH               = find(info.labels.ControlType.vals==20 & standard_labels);
% info.idx.RECON_END           = find(info.labels.ControlType.vals==21 & standard_labels);
% info.idx.IMAGE_STATUS        = find(info.labels.ControlType.vals==22 & standard_labels);
% info.idx.TRACKING            = find(info.labels.ControlType.vals==23 & standard_labels);
% info.idx.FLUOROSCOPY_TOGGLE  = find(info.labels.ControlType.vals==24 & standard_labels);
% info.idx.REJECTED_DATA       = find(info.labels.ControlType.vals==25 & standard_labels);
% info.idx.UNKNOWN27           = find(info.labels.ControlType.vals==27 & standard_labels);
% info.idx.UNKNOWN28           = find(info.labels.ControlType.vals==28 & standard_labels);
% 
% %% Calculate number of standard, normal data labels
% info.nNormalDataLabels = length(info.idx.NORMAL_DATA);
% 
% %% Dimension names
% dimnames = {'coil','kx','ky','kz','E3','loc','ec','dyn','ph','row','mix','avg'};
% dimfields = {'N/A','N/A','E1','E2','E3','Location','Echo','Dynamic','CardiacPhase','Row','Mix','Measurement'};
% 
% %% Initialize dimension data to zero
% info.dims.nCoils         = 0;
% info.dims.nKx            = 0;
% info.dims.nKy            = 0;
% info.dims.nKz            = 0;
% info.dims.nE3            = 0;
% info.dims.nLocations     = 0;
% info.dims.nEchoes        = 0;
% info.dims.nDynamics      = 0;
% info.dims.nCardiacPhases = 0;
% info.dims.nRows          = 0;
% info.dims.nMixes         = 0;
% info.dims.nMeasurements  = 0;
% 
% %% Calculate max number of active coils
% maxChannelsActiveMask = 0;
% for k=1:length(info.labels.ChannelsActive.uniq),
%     maxChannelsActiveMask = bitor(maxChannelsActiveMask,info.labels.ChannelsActive.uniq(k));
% end
% while maxChannelsActiveMask > 0
%     if bitand(maxChannelsActiveMask, 1),
%         info.dims.nCoils = info.dims.nCoils + 1;
%     end
%     maxChannelsActiveMask = bitshift (maxChannelsActiveMask, -1);
% end
% %% VINCENT - hardcoded to 16 channels
% %info.dims.nCoils = 8;
% 
% %% Calculate dimensions of normal data
% info.dims.nKx            = max(info.labels.DataSize.vals(info.idx.NORMAL_DATA)) / info.dims.nCoils / 2 / 2;
% info.dims.nKy            = length(unique(info.labels.E1.vals(info.idx.NORMAL_DATA)));
% info.dims.nKz            = length(unique(info.labels.E2.vals(info.idx.NORMAL_DATA)));
% info.dims.nE3            = length(unique(info.labels.E3.vals(info.idx.NORMAL_DATA)));
% info.dims.nLocations     = length(unique(info.labels.Location.vals(info.idx.NORMAL_DATA)));
% info.dims.nEchoes        = length(unique(info.labels.Echo.vals(info.idx.NORMAL_DATA)));
% info.dims.nDynamics      = length(unique(info.labels.Dynamic.vals(info.idx.NORMAL_DATA)));
% info.dims.nCardiacPhases = length(unique(info.labels.CardiacPhase.vals(info.idx.NORMAL_DATA)));
% info.dims.nRows          = length(unique(info.labels.Row.vals(info.idx.NORMAL_DATA)));
% info.dims.nMixes         = length(unique(info.labels.Mix.vals(info.idx.NORMAL_DATA)));
% info.dims.nMeasurements  = length(unique(info.labels.Measurement.vals(info.idx.NORMAL_DATA)));
% 
% %% With known possible dimension names, the load options can now be parsed
% p = inputParser;
% p.StructExpand = true;
% p.CaseSensitive = true;
% p.KeepUnmatched = false; % throw an error for unmatched inputs
% p.addRequired('filename', @ischar);
% for k=1:length(dimnames),
%     p.addParamValue(dimnames{k}, [], @isnumeric);
% end
% p.addParamValue('verbose', false, @islogical);
% p.addParamValue('savememory', true, @islogical);
% p.parse(filename, varargin{:});
% 
% %% Return loadopts structure inside INFO structure
% % remove filename field - it is passed as the first required argument
% info.loadopts = rmfield(p.Results,'filename');
% 
% %% Find the unique set of values for each dimension name
% info.dims.coil = [1:info.dims.nCoils];
% info.dims.kx   = [1:info.dims.nKx];
% for k=3:length(dimnames), % skip coil and kx
%     info.dims.(dimnames{k}) = unique(info.labels.(dimfields{k}).vals(info.idx.NORMAL_DATA));
% end
% 
% %% Find intersection of available dimensions with LOADOPTS dimensions
% for k=1:length(dimnames),
%     if ~isempty(info.loadopts.(dimnames{k})),
%         info.dims.(dimnames{k}) = intersect_a_with_b(info.loadopts.(dimnames{k}),info.dims.(dimnames{k}));
%     end
% end
% 
% %% Calculate data size
% datasize = []; 
% for k=1:length(dimnames),
%     datasize = [datasize length(info.dims.(dimnames{k}))];
% end
% info.datasize = datasize;
% 
% % throw error if any dimension size is zero
% if any(info.datasize==0),
%     zero_length_str = sprintf(' ''%s'' ', dimnames{find(info.datasize==0)});
%     error('size of selected data to load has zero length along dimension(s): %s', zero_length_str);
% end
% 
% %% Skip data loading if only one output argument is provided, return INFO
% if nargout==1,
%     info.labels_row_index_array = [1:size(info.labels,1)];
%     data=info;
%     return;
% end
% 
% %% Create array to hold label row numbers for loaded data
% % skip the coil and kx dimensions
% info.labels_row_index_array = zeros(datasize(3:end));
% 
% %% Pre-allocate DATA array
% if info.loadopts.savememory==true,
%     data = zeros(info.datasize,'single');
% else
%     data = zeros(info.datasize);
% end
% 
% %% Read RAW data for selected dimension ranges
% fidraw = fopen(rawname,'r','ieee-le');
% if fidraw<0,
%     error(sprintf('cannot open RAW file: %s', rawname));
% end
% info.nLoadedLabels=0;
% 
% raw_data_fread_size = double(info.dims.nCoils * info.dims.nKx * 2);
% rawdata_2d = complex(zeros(info.dims.nCoils,info.dims.nKx),zeros(info.dims.nCoils,info.dims.nKx));
% 
% for n=1:length(info.idx.NORMAL_DATA),
%     
%     load_flag=1;
%     dim_assign_indices_full_array = [];
%         
%     label_idx = info.idx.NORMAL_DATA(n);
%     
%     for k=3:length(dimfields),
%         
%         dimval = info.labels.(dimfields{k}).vals(label_idx);
%         
%         % it is allowed that the dimval appears more than once 
%         % in the requested dimension ranges to be loaded
%         dim_assign_indices = find(dimval==info.dims.(dimnames{k}));
%         
%         if isempty(dim_assign_indices),
%             load_flag=0;
%             break;
%         else
%            
%             if k>3,
%                 
%                 dim_assign_indices_full_array_new = zeros( size(dim_assign_indices_full_array,1)*length(dim_assign_indices), size(dim_assign_indices_full_array,2)+1);
%                 
%                 mod_base_a = size(dim_assign_indices_full_array,1);
%                 mod_base_b = length(dim_assign_indices);
%                 
%                 for d=1:size(dim_assign_indices_full_array_new,1),
%                     dim_assign_indices_full_array_new(d,:) = [dim_assign_indices_full_array(mod(d,mod_base_a)+1,:) dim_assign_indices(mod(d,mod_base_b)+1)];
%                 end
%                 
%             else
%                 dim_assign_indices_full_array_new = dim_assign_indices(:);
%             end
%             
%             dim_assign_indices_full_array = dim_assign_indices_full_array_new;
%             
%         end
%     end
%     
%     if load_flag,
%         
%         info.nLoadedLabels = info.nLoadedLabels+1;
%         
%         byte_offset = info.fseek_offsets(label_idx);
%         status = fseek(fidraw, byte_offset, 'bof');
%         rawdata_1d = double(fread(fidraw, raw_data_fread_size, 'int16'));
%         
%         % Phase correction
%         RandomPhase = double(info.labels.RandomPhase.vals(label_idx));
%         MeasurementPhase = double(info.labels.MeasurementPhase.vals(label_idx));;
%         c = exp (- i * pi * (2 * RandomPhase / (2^16-1) + MeasurementPhase / 2));
% 
%         % parse into real and imaginary parts for each coil
%         for kx=1:info.dims.nKx,
%           for coil=1:info.dims.nCoils,
%             re_idx = 2*info.dims.nKx*(coil-1) + 2*(kx-1) + 1;
%             im_idx = re_idx + 1;
%             rawdata_2d(coil,kx) = c * (rawdata_1d(re_idx) + i*rawdata_1d(im_idx)); 
%           end
%         end
%         
%         % account for measurement sign
%         if(info.labels.MeasurementSign.vals(label_idx)),
%             rawdata_2d = fliplr(rawdata_2d);
%         end
%         
%         % select choosen coils
%         rawdata_2d = rawdata_2d(info.dims.coil,:);
%         
%         % select choosen kx
%         rawdata_2d = rawdata_2d(:,info.dims.kx);
%         
%         % insert rawdata_2d into proper locations of the data array
%         for d=1:size(dim_assign_indices_full_array,1),
%             
%             dim_assign_str = sprintf(',%d', dim_assign_indices_full_array(d,:) );
%             
%             % delete initial comma
%             dim_assign_str(1) = [];
%             
%             % assign index to table_index table
%             eval( sprintf('info.labels_row_index_array(%s)=%d;', dim_assign_str, label_idx) );
%         
%             % assign read image to correct location in data array
%             eval( sprintf('data(:,:,%s)=rawdata_2d;', dim_assign_str) );
%                     
%         end
%     
%     end
%     
% end
% 
% %% Read FRC noise data
% 
% % AR Include an if-statement to avoid error if no noise data is present
% 
% if (~isempty(info.idx.FRC_NOISE_DATA))
%     frc_noise_idx = info.idx.FRC_NOISE_DATA(1);
%     frc_noise_samples_per_coil = info.labels.DataSize.vals(frc_noise_idx) / 2 / 2 / info.dims.nCoils;
%     info.FRC_NOISE_DATA = zeros(info.dims.nCoils,frc_noise_samples_per_coil,'single');
%     byte_offset = info.fseek_offsets(frc_noise_idx);
%     status = fseek(fidraw, byte_offset, 'bof');
%     rawdata_1d = double(fread(fidraw, double(info.labels.DataSize.vals(frc_noise_idx)/2) , 'int16'));
%     for sample=1:frc_noise_samples_per_coil,
%       for coil=1:info.dims.nCoils,
%         re_idx = 2*frc_noise_samples_per_coil*(coil-1) + 2*(sample-1) + 1;
%         im_idx = re_idx + 1;
%         info.FRC_NOISE_DATA(coil,sample) = rawdata_1d(re_idx) + i*rawdata_1d(im_idx); 
%       end
%     end
% end
% 
% % Close RAW file
% fclose(fidraw);
% 
% %% Calculate total raw data blobs
% size_data = size(data);
% max_img_dims = size_data(3:end);
% info.nDataLabels = prod(max_img_dims);
% 
% %% If VERBOSE, display execution information
% if info.loadopts.verbose==true,
%     disp( sprintf('Loaded %d of %d available normal data labels', info.nLoadedLabels, info.nNormalDataLabels) );
%     tmpstr = '';
%     for k=1:length(dimnames),
%         tmpstr = sprintf('%s, # %s: %d', tmpstr, dimnames{k}, length(info.dims.(dimnames{k})) );
%     end
%     disp( sprintf('Data contains %d raw labels - %s', info.nDataLabels, tmpstr(3:end)) );
%     disp( sprintf('Total execution time = %.3f seconds', toc) );
% end
% 
% %% Find intersection of vector a with vector b without sorting 
% function c = intersect_a_with_b(a,b)
% c = a;
% % work backwards in order to use [] assignment
% for k=length(a):-1:1,
%     if isempty(find(a(k)==b)),
%         c(k)=[]; 
%     end
% end
% 
% % force c to be a row vector for easier display
% c = c(:).';
