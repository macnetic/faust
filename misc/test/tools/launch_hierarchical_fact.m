%% launch the mexfunction mexhierarchical_fact 
%% and makes time comparison between faust and data matrix multiplication with vector



function [mexlambda,mexfact,fc]=mexhierarchical_fact(params)


[mexlambda,mexfact]=mexHierarchical_fact(params);
disp(['MEX LAMBDA ' num2str(mexlambda)]);

%% launch the same algorithm with the Luc's function 
% [lambda,fact]=hierarchical_fact(params);
% disp([' LAMBDA ' num2str(lambda) ' MEX_LAMBDA ' num2str(mexlambda)]);

% sum_coeff = 0;
% threshold = 0.00001;
% disp('\n\n COMPARAISON ');
% for i=1:params.nfacts
%    comp=abs(mexfact{i} - fact{i})> threshold;
%    nb_coeff_diff = sum(comp(:));
%    disp([int2str(i) ' nb_coeff diff : ' int2str(nb_coeff_diff)]);
%    sum_coeff = nb_coeff_diff + sum_coeff; 
%     
% end
% disp([' nb coeff different ' int2str(sum_coeff)]);

% 
% mex_error = norm(params.data- lambda*dvp(fact));


mexfact{1}=mexlambda*mexfact{1};
fc=matlab_faust(mexfact);
mex_error = norm(params.data - get_product(fc));

disp(['relative error :  ' num2str(mex_error)]);



%% speed-up test for multiplication with a vector
disp('*** product data matrix-vector vs product faust-vector ***');
nbiter = 100;
dense_mat = params.data;
[nl,nc]=size(dense_mat);
y_faust=zeros(nl,1);
y_dense=zeros(nl,1);
tps_dense=zeros(nbiter,1);
tps_faust=zeros(nbiter,1);

for i=1:nbiter
    x=rand(nc,1);
   
    tic
        y_faust = fc*x;
    tps_faust=toc;
    
    tic
        y_dense = dense_mat*x;
    tps_dense=toc;    
end


tps_dense=mean(tps_dense);
tps_faust=mean(tps_faust);
speed_up = tps_dense/tps_faust;
disp(['dense time : ' num2str(tps_dense) ' faust time : ' num2str(tps_faust)]);
disp(['speed_up : ' num2str(speed_up)]);

end


