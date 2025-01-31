function J=mtfftc2(data,tapers,nfft,Fs,memFlag)
% Multi-taper fourier transform - continuous data
%
% Usage:
% J=mtfftc(data,tapers,nfft,Fs) - all arguments required
% Input: 
%       data (in form samples x channels/trials or a single vector) 
%       tapers (precalculated tapers from dpss) 
%       nfft (length of padded data)
%       Fs   (sampling frequency)
%                                   
% Output:
%       J (fft in form frequency index x taper index x channels/trials)
if ~exist('memFlag','var'); memFlag = false; end

if nargin < 4; error('Need all input arguments'); end;
data=change_row_to_column(data);
[NC,C]=size(data); % size of data
[NK,K]=size(tapers); % size of tapers
if NK~=NC; error('length of tapers is incompatible with length of data'); end;
if ~memFlag
    J = computeVectorized(tapers,data,NC,K,C,nfft,Fs);
else
    nBloc = round(NC*K*C/(1e9)); if nBloc==0; nBloc = 1; end
    disp('to large time X space X taper matrix')
    % disp('using single precision')
    % tapers = single(tapers);
    % data = single(data);
    disp(['spliting in ' num2str(nBloc) ' blocs'])
    disp('and averaging power across tapers (phase will be lost)')
    dist = repmat(round(C/nBloc),[1 nBloc]); dist(end) = dist(end) + (C - sum(dist));
    data = mat2cell(data,NC,dist);
    J = cell(size(data));
    for bloc = 1:nBloc
        tic
        disp([num2str(bloc) '/' num2str(nBloc) ': computing'])
        C = dist(bloc);
        J{bloc} = computeVectorized(tapers,data{bloc},NC,K,C,nfft,Fs);
        % average power to save space (phase information is lost)
        J{bloc} = sqrt(mean(conj(J{bloc}).*J{bloc},2));
        % delete data to save space
        data{bloc} = {};
        disp([num2str(bloc) '/' num2str(nBloc) ': done'])
        toc
    end
    J = cat(3,J{:});
end
% if NC*K*C < 1000000000
%     J = computeVectorized(tapers,data,NC,K,C,nfft,Fs);
% elseif NC*K*C < 10000000000/2
%     disp('to large time X space X taper matrix')
%     disp('using single precision')
%     tapers = single(tapers);
%     data = single(data);
%     J = computeVectorized(tapers,data,NC,K,C,nfft,Fs);
% elseif NC*K*C < 10000000000
%     disp('to large time X space X taper matrix')
%     disp('using single precision and parallelization')
%     tapers = single(tapers);
%     data = single(data);
%     J = computeParalellized(tapers,data,NC,K,C,nfft,Fs);
% else
%     nBloc = round(NC*K*C/(1e9));
%     disp('to large time X space X taper matrix')
%     % disp('using single precision')
%     % tapers = single(tapers);
%     % data = single(data);
%     disp('using parallelization')
%     disp([' and spliting in ' num2str(nBloc) ' blocs'])
%     disp('averaging power across tapers')
%     dist = repmat(round(C/nBloc),[1 nBloc]); dist(end) = dist(end) + (C - sum(dist));
%     data = mat2cell(data,NC,dist);
%     J = cell(size(data));
%     for bloc = 1:nBloc
%         tic
%         disp([num2str(bloc) '/' num2str(nBloc) ': computing'])
%         C = dist(bloc);
%         J{bloc} = computeVectorized(tapers,data{bloc},NC,K,C,nfft,Fs);
%         % average power to save space (phase information is lost)
%         J{bloc} = sqrt(mean(conj(J{bloc}).*J{bloc},2));
%         % delete data to save space
%         data{bloc} = {};
%         disp([num2str(bloc) '/' num2str(nBloc) ': done'])
%         toc
%     end
%     J = cat(3,J{:});
% end

function J = computeVectorized(tapers,data,NC,K,C,nfft,Fs)
tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers
data=data(:,:,ones(1,K)); % add taper indices to data
data=permute(data,[1 3 2]); % [time X tapers X vox] [NC K C] reshape data to get dimensions to match those of tapers
data=data.*tapers; clear tapers % [time X tapers X vox] [NC K C] product of data with tapers
% varList = {{'data'}}; for i = 1:length(varList); eval(['tmp = whos(''' strjoin(varList{i},''',''') ''');']); eval(['Gb.' strjoin(varList{i},'') ' = sum([tmp.bytes])./1e9;']); end
% tic
J=fft(data,nfft)/Fs; %clear data   % fft of projected data
% varList = {'J'}; eval(['tmp = whos(''' strjoin(varList,''',''') ''');']); eval(['Gb.' strjoin(varList,'') ' = sum([tmp.bytes])./1e9;']);
% disp(['cost of ' strjoin(varList,'') ' is ' num2str(Gb.(strjoin(varList,''))./Gb.data) ' the size of data'])
% toc
% 
% tic
% [u,s,v] = svd(permute(sum(data,1),[3 2 1]),"econ"); % [vox X tapers] [C K]
% s = diag(s).^2;
% c = s/sum(s);
% varList = {'u' 's' 'v' 'c'}; eval(['tmp = whos(''' strjoin(varList,''',''') ''');']); eval(['Gb.' strjoin(varList,'') ' = sum([tmp.bytes])./1e9;']);
% disp(['cost of ' strjoin(varList,'') ' is ' num2str(Gb.(strjoin(varList,''))./Gb.data) ' the size of data'])
% toc


function J = computeParalellized(tapers,data,NC,K,C,nfft,Fs)
J = complex(zeros(NC,K,C,'single'),zeros(NC,K,C,'single'));
parfor k = 1:K
    J(:,k,:)=fft(data.*tapers(:,k),nfft)/Fs;   % fft of projected data
end