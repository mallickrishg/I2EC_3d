function [bigK,bigKn,kfile] = computeKernels(rcv,name,kdir,nCheck)
% COMPUTEKERNELS loads cached traction kernels for a mesh or computes and
% caches them
%
%   [bigK,bigKn,kfile] = unicycle.greens.computeKernels(rcv,name,kdir,nCheck)
%
% INPUT:
%   rcv    - unicycle.geometry.triangleReceiver
%   name   - (optional) label for the cache file, default 'mesh'
%   kdir   - (optional) cache folder, default 'kernels'
%   nCheck - (optional) number of random source patches recomputed on load
%            to spot-check the cached kernels, default 2 (0 to skip)
%
% OUTPUT:
%   bigK  - [Kss,Kds;Ksd,Kdd] (2N x 2N) shear traction kernels
%   bigKn - [Ksn,Kdn] (N x 2N) normal traction kernels
%   kfile - cache file used
%
% the cache file is <kdir>/<name>_<hash>.mat, where hash is a SHA-256
% fingerprint of the mesh vertices, triangle connectivity, Green's function
% class, rigidity and Poisson's ratio. A different geometry or elastic
% model therefore maps to a different file. On load the stored fingerprint
% is compared in full and nCheck kernel columns are recomputed and compared.
%
% SEE ALSO: unicycle.geometry.triangleReceiver

if nargin < 2 || isempty(name), name = 'mesh'; end
if nargin < 3 || isempty(kdir), kdir = 'kernels'; end
if nargin < 4 || isempty(nCheck), nCheck = 2; end

fingerprint = meshFingerprint(rcv);
kfile = fullfile(kdir,sprintf('%s_%s.mat',name,fingerprint(1:12)));

if exist(kfile,'file') == 2
    stored = load(kfile,'fingerprint');
    assert(isfield(stored,'fingerprint') && strcmp(stored.fingerprint,fingerprint),...
        'computeKernels: %s does not match the mesh geometry',kfile)
    load(kfile,'Kss','Kds','Ksd','Kdd','Ksn','Kdn')
    fprintf('loaded kernels from %s\n',kfile)

    % recompute a few columns and compare
    if nCheck > 0
        idx = randperm(rcv.N,min(nCheck,rcv.N));
        src = struct('x',rcv.x,'vertices',rcv.vertices(idx,:),'slip',zeros(numel(idx),1));
        [cKss,cKds,cKsd,cKdd] = rcv.earthModel.tractionKernels(src,rcv);
        cached = [Kss(:,idx),Kds(:,idx),Ksd(:,idx),Kdd(:,idx)];
        fresh = [cKss,cKds,cKsd,cKdd];
        relerr = max(abs(cached(:)-fresh(:)))/max(abs(fresh(:)));
        assert(relerr < 1e-8,...
            'computeKernels: cached kernels in %s differ from recomputed columns (rel. error %.2e)',kfile,relerr)
    end
else
    fprintf('computing kernels for %d patches\n',rcv.N)
    [Kss,Kds,Ksd,Kdd,Ksn,Kdn] = rcv.tractionKernels(rcv);
    if ~exist(kdir,'dir'), mkdir(kdir), end
    save(kfile,'Kss','Kds','Ksd','Kdd','Ksn','Kdn','fingerprint','-v7.3')
    fprintf('saved kernels to %s\n',kfile)
end

bigK = [Kss,Kds;Ksd,Kdd];
bigKn = [Ksn,Kdn];

end

function h = meshFingerprint(rcv)
% SHA-256 of the geometry and elastic model as a hex string
desc = sprintf('%s|G=%.17g|nu=%.17g|N=%d',class(rcv.earthModel),...
    rcv.earthModel.G,rcv.earthModel.nu,rcv.N);
bytes = [uint8(desc),typecast(double(rcv.x(:))','uint8'),...
    typecast(double(rcv.vertices(:))','uint8')];
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes);
h = lower(reshape(dec2hex(typecast(md.digest(),'uint8'))',1,[]));
end
