cd /home/tgautrai/faust2/debug

file=[];
fileList=dir;
cmpt=1;
for k=1:length(fileList)
    if ~isempty(strfind(fileList(k).name,'_device'))
        file{cmpt}=fileList(k).name;
        cmpt = cmpt+1;
    end
end


for k=1:length(file)
    dev=load(file{k});
    host=load(strrep(file{k},'_device','_host'));
    err_abs = max(max(abs(dev-host)));
    err_rel = max(max(abs(dev-host)./host));
    fprintf('err rel %s : %f\n',file{k},err_rel);
    fprintf('err abs %s : %f\n\n',file{k},err_abs);
end

a=0;
