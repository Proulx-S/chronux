function [sv,sp,fm] = spsvd(data,params,mdkp)
% 2022-03-07: modified by Sebastien Proulx (jsproulx@mgh.harvard.edu) to
% perform svd with the taper dimension augmented with frequency.
% 
% 
% Space frequency SVD of input data - continuous processes
% Usage: [sv,sp,fm] = spsvd(data,params,mdkp)
% Inputs:
% data       (data matrix in timexchannels form)-required
%       params      structure containing parameters - params has the
%       following fields: tapers, Fs, fpass, pad
%           tapers : precalculated tapers from dpss or in the one of the following
%                    forms:
%                   (1) A numeric vector [TW K] where TW is the
%                       time-bandwidth product and K is the number of
%                       tapers to be used (less than or equal to
%                       2TW-1).
%                   (2) A numeric vector [W T p] where W is the
%                       bandwidth, T is the duration of the data and p
%                       is an integer such that 2TW-p tapers are used. In
%                       this form there is no default i.e. to specify
%                       the bandwidth, you have to specify T and p as
%                       well. Note that the units of W and T have to be
%                       consistent: if W is in Hz, T must be in seconds
%                       and vice versa. Note that these units must also
%                       be consistent with the units of params.Fs: W can
%                       be in Hz if and only if params.Fs is in Hz.
%                       The default is to use form 1 with TW=3 and K=5
%
%	        Fs 	        (sampling frequency) -- optional. Defaults to 1.
%           fpass       (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
% mdkp       (number of dimensions to be kept)-optional. Default is the
%               maximum possible modes determined by taper parameters
%
% Outputs:
% sv sp fm  : singular values, space modes, frequency modes


if nargin < 1; error('Need data'); end;
if nargin < 2 || isempty(params); params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if isfield(params,'anaType')
    anaType = params.anaType;
end
clear err trialave params
[N,NCHAN]=size(data);
tapers=dpsschk(tapers,N,Fs);
nfft=max(2^(nextpow2(N)+pad),N);% number of points in fft
[N,K]=size(tapers);

tvec=1/Fs *(1:N)';
tvec=repmat(tvec,[1 K]);
tvec=tvec*2*pi*1i;
f=getfgrid(Fs,nfft,fpass);
nf=length(f);
if exist('anaType','var') && strcmp(anaType,'proulx')
    if nargin<3 || isempty(mdkp); mdkp=min(K*nf,NCHAN);
    elseif mdkp > min(K*nf,NCHAN); error('mdkp has to be less than both K and NCHAN');end;
    
%     sp=zeros(NCHAN,1,mdkp);
%     sp=sp+1i*sp;
%     fm=zeros(K,nf,mdkp);
%     fm=fm+1i*fm;
%     sv=zeros(1,mdkp);

    proj = tapers.*exp(permute(-f,[1 3 2]).*tvec); % (time X taper X freq)
    tmp=data'*proj(:,:); % projected data (chan X taper+freq)
    [sp,sv,fm]= svd(tmp,0); % svd

%     sv = diag(sv(1:mdkp,1:mdkp))';
    sv = diag(sv)';
    sp = permute(sp(:,1:mdkp),[1 3 2]);
    fm = permute( reshape( permute(fm(:,1:mdkp),[2 1]) ,[mdkp K nf]) ,[2 3 1]);
    
%     figure('WindowStyle','docked');
%     plot(s)
% 
%     % reduce
%     kf = 25;
%     u(:,kf+1:end) = [];
%     s(:,kf+1:end) = [];
%     v(:,kf+1:end) = [];
% 
%     u; % space X mode
%     s; % mode
%     v = permute(v,[2 1]);
%     v = reshape(v,[kf K nf]); % mode X taper X freq
%     v = squeeze(mean(v,2))'; % freq X mode
%     
% 
%     figure('WindowStyle','docked');
%     for kfInd = 1:size(v,1)
%         y = squeeze(mean(abs(v(kfInd,:,:)),2)) + 0.02*kfInd;
%         plot(f,y);
%         text(0.5,mean(y(end-10:end)),num2str(kfInd))
%         hold on
%     end
%     ax = gca;
%     ax.XScale = 'log';
% %     ax.YScale = 'log';
%     grid on
%     xlim([0.01 0.5])
    
else
    if nargin<3 || isempty(mdkp); mdkp=min(K,NCHAN);
    elseif mdkp > min(K,NCHAN); error('mdkp has to be less than both K and NCHAN');end;

    sp=zeros(NCHAN,nf,mdkp);
    sp=sp+1i*sp;
    fm=zeros(K,nf,mdkp);
    fm=fm+1i*fm;
    sv=zeros(nf,min([K,NCHAN]));
    % Original Mitra
    if nf>10
        parfor j=1:nf
            %     for k=1:K
            %       proj(:,k)=tapers(:,k).*exp(-f0*tvec');
            %     end
            proj=tapers.*exp(-f(j)*tvec);
            tmp=data'*proj; % projected data
            [u,s,v]= svd(tmp,0); % svd
            for mk=1:mdkp,
                sp(:,j,mk)=u(:,mk)';
                fm(:,j,mk)=v(:,mk)';
            end
            sv(j,:)=diag(s);
        end
    else
        for j=1:nf
            %     for k=1:K
            %       proj(:,k)=tapers(:,k).*exp(-f0*tvec');
            %     end
            proj=tapers.*exp(-f(j)*tvec);
            tmp=data'*proj; % projected data
            [u,s,v]= svd(tmp,0); % svd
            for mk=1:mdkp,
                sp(:,j,mk)=u(:,mk)';
                fm(:,j,mk)=v(:,mk)';
            end
            sv(j,:)=diag(s);
        end
    end
end
