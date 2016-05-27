function y = f_mult(facts,x)
% y=ones(size(facts{1},2),size(x,2));
% disp('bla');
% for i=1:100
% end
y=x;
for i=length(facts):-1:1
    y=facts{i}*y;
end


end

