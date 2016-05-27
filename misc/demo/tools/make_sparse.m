function new_facts = make_sparse(facts)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

for i=1:length(facts)
   new_facts{i}=sparse(facts{i}); 
end

end

