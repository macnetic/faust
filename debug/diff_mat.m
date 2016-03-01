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

%dev=cell(length(file),1);
%host=cell(length(file),1);
for k=1:length(file)
    dev=load(file{k});
    host=load(strrep(file{k},'_device','_host'));
    err_abs = max(max(abs(dev-host)));
    err_rel = max(max(abs(dev-host)./host));
    if (~isnan(err_rel) && abs(err_rel)>0.01) && (err_abs>1e-6)
      fprintf('err rel %s : %f%%\n',file{k},err_rel*100);
      fprintf('err abs %s : %e\n\n',file{k},err_abs);
   end
end

a=0;
