function new_index = BitReversalPermutation(index)
n=length(index);



if (n == 1)
   new_index = index; 
else
    index1=BitReversalPermutation(index(1:2:end));
    index2=BitReversalPermutation(index(2:2:end));
    new_index=[index1,index2];
end


    

end

