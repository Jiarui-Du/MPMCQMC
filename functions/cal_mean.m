function result = cal_mean(xx,Ns,m,ml,R,starv)
if isempty(starv)
    flag = 0;
else
    flag = 1;
end
if ml ~= 1 && m == 1
    m = ml;
end
dd = size(xx,2);
rmean = zeros(Ns,dd,m);
rvar = zeros(Ns,dd,m);
rmse = zeros(Ns,dd,m);
rbias = zeros(Ns,dd,m);
for j = 1:m 
    ij = (j-1)*Ns*R;
    for i = 1:Ns 
        temp = xx(1+(i-1)*R+ij:i*R+ij,:);
        rmean(i,:,j) = mean(temp,1);
        rvar(i,:,j) = var(temp,0);
        if flag == 1
            rmse(i,:,j) = mean((temp-starv).^2);
            rbias(i,:,j) = (mean(temp)-starv).^2;
        end
    end
end
result_mean = zeros(m,Ns,dd);
result_var = zeros(m,Ns,dd);
result_bias = zeros(m,Ns,dd);
result_mse = zeros(m,Ns,dd);
for i = 1:m
    for  j = 1:Ns
        for k = 1:dd
            result_mean(i,j,k) = rmean(j,k,i);
            result_var(i,j,k) = rvar(j,k,i);
            if flag == 1
                result_mse(i,j,k) = rmse(j,k,i);
                result_bias(i,j,k) = rbias(j,k,i);
            end
        end
    end
end
result_msefactor = zeros(m,Ns,dd);
result_factor = zeros(m,Ns,dd);
for k = 1:dd
    for j = 1:Ns
        if flag == 1
            result_msefactor(:,j,k) = result_mse(:,1,k)./result_mse(:,j,k);
        else
            result_factor(:,j,k) = result_var(:,1,k)./result_var(:,j,k);
        end
    end
end
if flag == 1
    result = [result_mean;result_bias;result_var;result_mse;result_msefactor];
else
    result = [result_mean;result_var;result_factor];
end
end

