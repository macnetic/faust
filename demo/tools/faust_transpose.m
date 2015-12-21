function facts_transpose = faust_transpose(facts)

facts_transpose=cell(1,length(facts));
for i=1:length(facts);
%     length(facts)+1-i
    facts_transpose{i}=facts{length(facts)+1-i}';
end

end

