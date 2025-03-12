function [sv,sp,fm,faMTS] = spsvd3(data,params,mdkp)
% 2022-03-07: modified by Sebastien Proulx (jsproulx@mgh.harvard.edu) to
% perform svd with the taper dimension augmented with frequency.
% 2023-10-12: further modified by Sebastien Proulx
% (jsproulx@mgh.harvard.edu) to generalize the augmentation of the taper
% dimension with frequency

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

paramOrig = params;

%% Define a multitaper space (MTS) for each analysis band
MTS = cell(length(paramOrig.BW),1);
MTS_tpInd = cell(length(paramOrig.BW),1);
MTS_f = nan(length(paramOrig.BW),2);
MTS_bandInd = cell(length(paramOrig.BW),1);
[N,NCHAN]=size(data);
for bandInd = 1:length(paramOrig.BW)
    param = paramOrig;
    param.tapers = paramOrig.tapers(bandInd,:);
    param.fpass = paramOrig.fpass(bandInd,:);
    param.BW = paramOrig.BW(bandInd,:);

    [MTS{bandInd},MTS_f(bandInd,:)] = getMTS(param,N,NCHAN);
    MTS_tpInd{bandInd} = 1:size(MTS{bandInd},2);
    MTS_bandInd{bandInd} = ones(size(MTS_tpInd{bandInd})).*bandInd;
end
if nnz(diff(paramOrig.fpass,[],2)~=0); error('in param.fpass, first and second column should be the same (fpass is used as a single frequency instead of a band here)'); end
if nnz(paramOrig.fpass(:,1)~=MTS_f); error('fpass in param and out of getMTS should be the same'); end

%% Create frequency augmented multitaper space (faMTS) by concatenating MTS from each frequency band
faMTS.proj = cat(2,MTS{:}); clear MTS
faMTS.tpInd = cat(2,MTS_tpInd{:}); clear MTS_tpInd
faMTS.bandInd = cat(2,MTS_bandInd{:}); clear MTS_bandInd
faMTS.bandFreq = permute(MTS_f(:,1),[2 1]);
faMTS.info = 'time x taper/freq x vox';

%% Project data into the faMTS
data = data'*faMTS.proj;

%% Scale
faMTS.projScale = ones(size(data,1),length(faMTS.bandFreq));
if param.scale
    faMTS.projScale = mean(abs(data),2);
    data = data ./ faMTS.projScale;
end
faMTS.projScale = permute(faMTS.projScale,[3 2 1]);

%% Scale each frequency subspaces
faMTS.projScale = ones(size(data,1),length(faMTS.bandFreq));
if param.scale
    for i = 1:length(faMTS.bandFreq)
        faMTS.projScale(:,i) = mean(abs(data(:,faMTS.bandInd==i)),2);
        data(:,faMTS.bandInd==i) = data(:,faMTS.bandInd==i) ./ faMTS.projScale(:,i);
    end
end
faMTS.projScale = permute(faMTS.projScale,[3 2 1]);

%% Perform svd in frequency augmented multitaper space
[u,s,v]= svd(data,0);
if ~exist('mdkp','var') || isempty(mdkp); mdkp = size(u,2); end
sp = permute(u(:,1:mdkp)',[2 3 1]); % WARNING! The ' operator on complex number also outputs its complex conjugate. This is desired here as this is what the original Chronux toolbox is doing.
fm = permute(v(:,1:mdkp)',[2 3 1]); % WARNING! The ' operator on complex number also outputs its complex conjugate. This is desired here as this is what the original Chronux toolbox is doing.
sv=permute(diag(s(1:mdkp,1:mdkp))',[1 3 2]);
% whos sp fm sv

