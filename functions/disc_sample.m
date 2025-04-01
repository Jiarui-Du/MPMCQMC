function x = disc_sample(n,p,r)
p = p(:);
sp = cumsum(p);
sp = [0;sp];
% x = 1;
x = [];
for k = 1:n
    if (sp(k) <= r) && (r < sp(k+1))
    x = k;
    break
    end
end
if isempty(x)
    disp('error\n')
    % disp(num2str(n))
    % disp(num2str(p))
    % disp(num2str(r))
end

        