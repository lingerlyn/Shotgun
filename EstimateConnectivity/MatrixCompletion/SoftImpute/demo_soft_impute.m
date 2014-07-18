%% test matrix completion / demo for using Soft-Impute.


%***** Install matlab-mex files:
% Type in matlab command line (root directory of SoftImpute)
% >> install_mex
% Creates mex files  PROPACK_utils/dbdqr.c    PROPACK_utils/reorth.c   Matlab_files/project_obs_UV.c

%%**** requires folders Matlab_files and PROPACK_utils in path

addpath('Matlab_files/')
addpath('PROPACK_utils/')


%% generate data
clear
randn('state',2009);
rand('state',2009);

% nrow=2000; ncol=2000; rnk=10; 
nrow=1000; ncol=1000; rnk=5;
density=0.05;
% nrow=30; ncol=30; rnk=3;
% density=.2;

R=sprand(nrow,ncol,density);
[i_row,j_col,temp_ones]=find(R); clear temp_ones  R*

%------- low-rank factorization
Ueff=randn(nrow,rnk);  Veff=randn(ncol,rnk); 
X_clean=Ueff*Veff';  % population outer-product 

temp_vec=Ueff(i_row,:);  temp_vecr=Veff(j_col,:);
TTemp=dot(temp_vec,temp_vecr,2); 
sig=0.05;
TTemp=TTemp+ sig*randn(length(TTemp),1); % this creates the observed matrix with NOISE

%------ sparse-observed matrix "GXobs" 
GXobs=sparse(i_row,j_col,TTemp,nrow,ncol);  clear TTemp
GPm=sparse(i_row,j_col,1,nrow,ncol);

%% declare fields of structure OPTS
OPTS=[];
OPTS.TOLERANCE=10^-4; OPTS.MAXITER=100; OPTS.SMALL_SCALE=0; 
lambda_max=spectral_norm(GXobs); % approximates the lambda value for which solution is zero

%%**************************************
%% create a path of solutions
%%************************************** 
path_length=100;
error=zeros(1,path_length);
lambda_seq=linspace(lambda_max*.9,lambda_max/5,path_length);

INIT=[];

%%
for i = 1:path_length
[a11,a22,a33,obj]=soft_impute(GXobs,lambda_seq(i),OPTS,INIT); % warm-start specified via INIT

%Do what you want to do with the output, for example compute the relative precition error
error(i)= norm((1-GPm).*(a11*a22*a33' - X_clean),'fro')/norm((1-GPm).*(X_clean),'fro');
disp(rank(a11*a22*a33'))
% keyboard

INIT=struct('U', a11 , 'D', a22,  'V', a33 ); % specify warm-starts for smaller lambda value

end



%%**************************************
%% solution at a single lambda
%%************************************** 

%% To obtain solution at a single lambda value, it is preferable to use the above strategy ie creating a path of solutions.
%% However, the following stand-alone call to soft_impute(...) also works 

lambda=lambda_seq(5);

%% STRATEGY 1: stand-alone call to soft_impute(..)
% call soft_impute() directly for a single value of lambda
% to keep a bound on the max possible rank change the field value

OPTS.MAX_RANK=50; 
% Warning: The above assignment assumes the optimal solution at lambda has rank < MAX_RANK, 

[a11,a22,a33,obj]=soft_impute(GXobs,lambda,OPTS);


%% STRATEGY 2: (preferable) create a path of solutions to obtain solution at lambda
% call soft_impute_path() directly for a single value of lambda
OPT=[];
OPTS.NUMBER_GRID=100; %10 grid values of lambda, smallest one at lambda
% [a11,a22,a33,obj]=soft_impute_path(GXobs,lambda,OPTS); 
%%
lambda=0;
OPTS.NUMBER_GRID=200;
OPTS.MAX_RANK=5;
[Glr_mat_u,Glr_mat_d,Glr_mat_v,allGlr_mat_u,allGlr_mat_d,allGlr_mat_v,obj]=...
    soft_impute_path(GXobs,lambda,OPTS);

ranks=zeros(OPTS.NUMBER_GRID,1);
imputeF=cell(OPTS.NUMBER_GRID,1);
errors=ranks;
for i=1:OPTS.NUMBER_GRID
    imputeF{i}=allGlr_mat_u{i}*allGlr_mat_d{i}*allGlr_mat_v{i}';
    ranks(i)=rank(imputeF{i});
%     errors(i)=norm(imputeF{i}-GXobs,'fro')/norm(GXobs,'fro');
    errors(i)=norm(imputeF{i}-X_clean,'fro')/norm(X_clean,'fro');
end

figure;
subplot(2,1,1); plot(ranks,'g'); title('ranks')
subplot(2,1,2); plot(errors,'r'); title('errors')















