function [sv,sp,fm,tapers,tvec,f] = spsvd_runAugmtd(data,params,mdkp,M)
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

if ~isfield(params,'scaleFlag')
    params.scaleFlag = 0;
end
if nargin < 1; error('Need data'); end;
if nargin < 2 || isempty(params); params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if isfield(params,'anaType')
    anaType = params.anaType;
end
scaleFlag = params.scaleFlag;
clear err trialave params
[N,NCHAN,NRUN]=size(data,[1 2 4]);
tapers=dpsschk(tapers,N,Fs);
nfft=max(2^(nextpow2(N)+pad),N);% number of points in fft
[N,K]=size(tapers);

tvec=1/Fs *(1:N)';
tvec=repmat(tvec,[1 K]);
tvec=tvec*2*pi*1i;
f=getfgrid(Fs,nfft,fpass);
nf=length(f);

if nargin<3 || isempty(mdkp); mdkp=min(K,NCHAN)*NRUN;
elseif mdkp > min(K,NCHAN); error('mdkp has to be less than both K and NCHAN');end;

sp=zeros(NCHAN,nf,mdkp);
sp=sp+1i*sp;
fm=zeros(K*NRUN,nf,mdkp);
fm=fm+1i*fm;
sv=zeros(nf,min([K,NCHAN])*NRUN);

if ~exist('M','var') || isempty(M); M=1; end
if nf<10; M=0; end
tic
parfor (j = 1:nf, M)
    a=data(:,:)'*(tapers.*exp(-f(j)*tvec));
    a = reshape(a',[K NCHAN NRUN]);
    if scaleFlag
        a = a ./ mean(abs(a),1);
    end
    a = permute(a,[2 1 3]);
    [u,s,v]= svd(a(:,:),0); a = [];
    for mk=1:mdkp
        sp(:,j,mk)=u(:,mk)';
        fm(:,j,mk)=v(:,mk)';
    end
    sv(j,:)=diag(s);
end
toc
