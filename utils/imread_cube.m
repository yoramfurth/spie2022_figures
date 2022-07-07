function data = imread_cube(file_to_open)
%IMREAD_CUBE opens a datacube binary file.
%
%Description: 
%    This function opens datacube binary files, given in BSQ, BIP, BIL, IMG formats.
%    The data file is assumed to have an matching header-file with ".hdr" extension.
% 
%Inputs: 
%    file_to_open - The datacube name, including extension. 
%
%Outputs: 
% 	 data - The datacube, loaded from the binary file.
% 
%Example: 
%    data = imread_cube('mydir\self_test_rad.img');  % read a data-cube stored in "img" format.
% 
%See also IMREAD
%

%    Updated by Ori Raviv, 2003 as "cubeRead.m"
%    Dept. Electrical & Computer Engineering, BGU Israel.
%
%    Updated by Yoram Furth, 2017 and renamed to "imread_cube.m"
%    Dept. Electrical & Computer Engineering, BGU Israel.

% ensures that an ".hdr" file exists
if (strcmpi(file_to_open(length(file_to_open)-2:length(file_to_open)),'hdr')) %check last 3 chars
    error('Please use the data file name, including extension (if exists)');
end
    
% reads the header-file
if file_to_open(length(file_to_open)-3)=='.'  %if the file has an extension (other than hdr)
	temp_file=zeros([1 length(file_to_open)]);
	temp_file=char(temp_file);
	temp_file(1:length(temp_file))=file_to_open(1:length(temp_file));
	temp_file(length(temp_file)-2:length(temp_file))='hdr';
	[samples, lines, bands, data_type, interleave]=get_header_param(temp_file);
	clear temp_file;
else %if the file has no extension
	temp_file=zeros([1 length(file_to_open)+4]);
	temp_file=char(temp_file);
	temp_file(1:length(file_to_open))=file_to_open(1:length(file_to_open));
	temp_file(length(temp_file)-3:length(temp_file))='.hdr'; %add extension .hdr
	[samples, lines, bands, data_type, interleave] = get_header_param(temp_file);
	clear temp_file;
end

% determins the data-type
bands_2_read=(1:bands);
switch data_type
    case 1
        type='uint8'; % 8 bits
    case 2
        type='short'; %signed integer,  16 bits
    case 3
        type='ulong';
    case 4
        type='float';
end

% reads the data-cube
num_bytes=[1 2 4 4];
fid=fopen(file_to_open,'r');
switch interleave
    case 'bsq' % ************* Band Sequential (BSQ)
        w=waitbar(0,'Reading datacube, BSQ format');
        waitbar(0,w);
        j=1; % Index for data if not all datacube is read.
        for i=bands_2_read
            waitbar(i/max(bands_2_read),w)
            indx=samples*lines*(i-1)*num_bytes(data_type);
            [~]=fseek(fid,indx,'bof');
            [help_mat, ~]=fread(fid,[samples lines],type);
            data(:,:,j)=help_mat(:,:)'; %#ok
            j=j+1;
        end
        close(w)
    case 'bip' %Band Interleaved by Pixel
        j=1;% Index for data if not all datacube is read.
        for i=bands_2_read  %:bands
            fseek(fid,0,-1);
            [~,~]=fread(fid,i-1,type);
            bytes_to_skip=num_bytes(data_type)*(bands-1);
            [help_mat,~]=fread(fid,[samples lines],type,bytes_to_skip);
            data(:,:,j)=help_mat(:,:)'; %#ok
            j=j+1;
        end
    case 'bil' %Band Interleaved by Line
        fseek(fid,0,-1); %rewind
        for j=1:lines
            for i=1:bands
                [help_mat,~]=fread(fid,samples,type);
                data(j,:,i)=help_mat(:)'; %#ok
            end
        end
end
fclose(fid);


function [samples, lines, bands, data_type, interleave] = get_header_param(header_file)
% reads the header file
fid=fopen(header_file,'r');
t=fgetl(fid);
while t~=-1
    i=findstr('samples',t);%#ok
    if isempty(i)~=1
        samples=t(11:length(t));
        samples=str2double(samples);
    end
    i=findstr('lines',t);%#ok
    if isempty(i)~=1
        lines=t(11:length(t));
        lines=str2double(lines);
    end
    i=findstr('bands',t);%#ok
    if isempty(i)~=1
        bands=t(11:length(t));
        bands=str2double(bands);
    end
    i=findstr('data type',t);%#ok
    if isempty(i)~=1
        data_type=t(length(t));
        data_type=str2double(data_type);
    end
    i=findstr('interleave',t);%#ok
    if isempty(i)~=1
        interleave=t(14:length(t));
    end
    t=fgetl(fid);
    
end
fclose(fid);

