function [u,bi] = FELFSR(m,N)

gs = [115,291,172,267,332,388,283,514,698,706,1304,920,1336,1236,1511];
% g的取法与N-1互素，后续u的取值也刚好以N-1为周期
g = gs(m-9);
switch m-9
    case 1
        a = [0,3];
    case 2
        a = [0,2];
    case 3
        a = [0,1,4,6];
    case 4
        a = [0,1,3,4];
    case 5
        a = [0,1,3,5];
    case 6
        a = [0,1];
    case 7
        a = [0,2,3,5];
    case 8
        a = [0,3];
    case 9
        a = [0,7];
    case 10
        a = [0,1,2,5];
    case 11
        a = [0,3];
    case 12
        a = [0,2];
    case 13
        a = [0,1];
    case 14
        a = [0,5];
    case 15
        a = [0,1,3,4];
end
a = flip(m-a);
mm = N-1;
b = false(mm,1);
bi = false(m,mm);
% b的周期为N-1
b(1) = 1;
k = length(a);
for i = (m+1):mm
    temp = 0;
    for j = 1:k
        temp = temp+b(i-a(j));
    end
    b(i) = mod(temp,2);
end
b = [b;b(1:m-1)];
for i = 1:mm
    temp1 = mod((i-1)*g+1,mm);
    temp1 = temp1+(temp1==0)*mm;
    temp2 = temp1+m-1;
    bi(:,i) = b(temp1:temp2);
end
u = (2.^(-(1:m)))*bi;
end