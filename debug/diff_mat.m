cd /home/tgautrai/faust2/debug

disp('Creation de la liste de fichiers a traiter')
file=[];
fileList=dir;
cmpt=1;
for k=1:length(fileList)
    if ~isempty(strfind(fileList(k).name,'_0_device'))
        file{cmpt}=fileList(k).name;
        cmpt = cmpt+1;
    end
end

%dev=cell(length(file),1);
%host=cell(length(file),1);
try
    for k=1:length(file)
        dev=load(file{k});
        host=load(strrep(file{k},'_device','_host'));
        if ~isempty(dev)
            err_abs = max(max(abs(dev-host)));
            err_rel = max(max(abs(dev-host)./host));
            if (~isnan(err_rel) && abs(err_rel)>0.01) || (err_abs>1e-6)
                fprintf('\nerr rel %s : %f%%\n',file{k},err_rel*100);
                fprintf('\nerr abs %s : %e\n\n',file{k},err_abs);
            else
                fprintf('.');
            end
        else
            fprintf('.');
        end
    end
catch err
    keyboard;
end
fprintf('\n');

a=0;

L1d=load('LorR1_avant_0_device.tmp');
L2d=load('LorR1_apres_0_device.tmp');
L2h=load('LorR1_apres_0_host.tmp');
L1h=load('LorR1_avant_0_host.tmp');
Sh=load('S1_4_0_host.tmp');
Sd=load('S1_4_0_device.tmp');
res=L1d*Sd;

max(max(abs((L1d-L1h)./L1h)))
max(max(abs((Sh-Sd)./Sh)))

max(max(abs((res-L2h))))
max(max(abs((res-L2d))))
