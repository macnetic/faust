list_im={'hier_fact_tree','.png';
         'hier_fact_pseudo_code','.png'};


for i=1:size(list_im,1)
   name_file=list_im{i,1};
   ext_file=list_im{i,2};
   im=imread([name_file,ext_file]);
   imwrite(im,[name_file,'.jpg']);
   
    
end