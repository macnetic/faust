%% test the mexhierarchical_fact function in the configuration described in 
% the file paramsfile
function hier_fact_test(paramsfile,expectedLambda, expectedLambdaPrecision)



% load the hierarchical_fact configuration
disp(['*** LOADING PARAMS FILE ***']);
disp([paramsfile]);
disp(' ');
disp(' ');
 
load(paramsfile);

%% factorisation (mexfile)
disp('*** MEX FACTORISATION ***');  
tic
[mexlambda,mexfact]=mexHierarchical_fact(params);
t=toc;

disp(['time factorisation (mex) : ' num2str(t)]);  



%% check if the result are ok
disp(['lambda value () : ' num2str(mexlambda)]);
if (abs(mexlambda - expectedLambda) > expectedLambdaPrecision)
    disp(' ');
    
    disp([ 'expected lamba value : ' int2str(expectedLambda) ' in the precision of ' int2str(expectedLambdaPrecision) ]);	
    error('invalid lambda value');
end


mexfact{1}=mexlambda*mexfact{1};
fc=Faust(mexfact);
mex_error = norm(params.data - get_product(fc));

disp(['relative error :  ' num2str(mex_error)]);




%% factorisation (matlab)
disp(' ');
disp(' ');
disp(['*** MATLAB FACTORISATION ***']); 
tic
[lambda,fact]=old_hierarchical_fact(params);
t=toc;

disp(['time factorisation (matlab) : ' num2str(t)]);  
disp(['lambda value (MATLAB) : ' num2str(lambda)]);





%% speed-up test for multiplication with a vector
disp(' ');
disp(' ');
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
