% detrend data (1850 - 2100)

vars={'thetao' 'so' 'dissic' 'talk' 'o2' 'po4'};
% define model and variable
model='GFDL-CM4';
run='ssp245';
variant='r1i1p1f1';

% loop over vars
for l = 1:length(vars)
    
var=vars{l};%'thetao';
if strcmp(run,'historical')
    N0=1;
else
    N0=165; % 1 for historical, % (2015-1850)= 165 for future scenario runs
end

% define file name
fn0=[model,'/',model,'_piControl_',variant,'_2x2L33_ann.nc'];
fn1=[model,'/',model,'_',run,'_',variant,'_2x2L33_ann.nc'];
wn=[model,'/',var,'_',model,'_',run,'_',variant,'_2x2L33_ann_dft.mat'];
disp(['working on ',var]);

% load picontrol data
pi=ncread(fn0,var,[1 1 1 1],[180 90 33 250]);

% load data
data=ncread(fn1,var);
x=ncread(fn1,'lon');
y=ncread(fn1,'lat');
z=ncread(fn1,'lev');
time=ncread(fn1,'time');
Nt=length(time);

% calculate model drift
N=size(pi);
X=[1:N(4)]-0.5;
tmp=zeros(N);
for i=1:N(1)
%    disp(['   working on i = ',num2str(i),' ...']);
    for j=1:N(2)
        for k=1:N(3)
            if ~isnan(pi(i,j,k,1))
                Y=squeeze(pi(i,j,k,:));
                C=cov(X,Y);
                A(i,j,k)=C(1,2)/C(1,1);
                tmp(i,j,k,:)=A(i,j,k)*(X+N0-1);
            end
        end
    end
end
clear X Y C;

% drift corrected data
dataDC=data-tmp(:,:,:,N0:(N0+Nt-1));
eval([var,'=dataDC;']);

save(wn,'-v7.3',var);

end

disp('done!');


