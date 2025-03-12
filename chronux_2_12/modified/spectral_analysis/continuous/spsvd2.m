function [u,s,v,f,bandV] = spsvd2(data,params,mdkp)
% Modified by Sebastien Proulx to accept multiple frequency range
% definitions and concatenate the resulting sets of tapers into a single
% svd.

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

funTs = data;
data = funTs.vec;
tr = funTs.tr/1000;
nFrame = size(funTs.vec,1);

if nargin < 1; error('Need data'); end;
if nargin < 2 || isempty(params); params=[]; end;
paramsOrig = params;

proj = cell(size(paramsOrig.fpass,1),1);
tapers = cell(size(paramsOrig.fpass,1),1);
f = cell(size(paramsOrig.fpass,1),1);
fpassReal = cell(size(paramsOrig.fpass,1),1);
for bandInd = 1:size(paramsOrig.fpass,1)
    param.Fs = paramsOrig.Fs;
    fpass = paramsOrig.fpass(bandInd,:);
    W = diff(fpass)/2;
    T = tr.*nFrame;
    TW = T*W;
    K = round(TW*2-1);
    TW = (K+1)/2;
    param.tapers = [TW K];
    [~,fx] = mtspectrumc(data(:,1), param);
    f0 = fpass(1)+W; [~,b] = min(abs(fx - f0)); f0 = fx(b);
    param.fpass = [f0 f0];
    mdkp = [];
    %%% Display actual frequency band used
    Wreal = W;
    fpassReal{bandInd} = f0+[-1 1].*(TW/T);
    
    [tapers{bandInd},pad,Fs,fpass,err,trialave,param]=getparams(param);
    clear err trialave param
    [N,NCHAN]=size(data);
    tapers{bandInd}=dpsschk(tapers{bandInd},N,Fs);
    nfft=max(2^(nextpow2(N)+pad),N);% number of points in fft
    [N,K]=size(tapers{bandInd});
    if nargin<3 || isempty(mdkp); mdkp=min(K,NCHAN);
    elseif mdkp > min(K,NCHAN); error('mdkp has to be less than both K and NCHAN');end;

    tvec=1/Fs *(1:N)';
    tvec=repmat(tvec,[1 K]);
    tvec=tvec*2*pi*i;
    f{bandInd}=getfgrid(Fs,nfft,fpass);
    nf=length(f{bandInd});

    proj{bandInd}=tapers{bandInd}.*exp(-f{bandInd}*tvec);
end

% Catenate tapers from all bands, first keeping track of which taper
% belongs to which frequency band
Kind = [];
for bandInd = 1:size(paramsOrig.fpass,1)
    K = size(tapers{bandInd},2);
    Kind = cat(1,Kind,ones(K,1).*bandInd);
end
proj = cat(2,proj{:});
tapers = cat(2,tapers{:});
f = cat(2,f{:})';

% Project data
tmp=data'*proj;
% Perform svd
[u,s,v]= svd(tmp,0); % svd

% Average v within each frequency band (should I the sqrt)
bandV = nan(1,size(paramsOrig.fpass,1));
for bandInd = 1:size(paramsOrig.fpass,1)
    bandV(bandInd) = mean(abs(v(Kind==bandInd)));
end
bandV = bandV';
