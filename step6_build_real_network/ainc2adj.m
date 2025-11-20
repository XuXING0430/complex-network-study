function b  = ainc2adj( x )
%AINC2ADJ   x为邻接表（可在工作区中新建数据），输出b为邻接矩阵
%   此函数通过邻接表生成临接矩阵的
if min(x(:))==0;
   x=x+1;
end
 d=length(x);
 a=max(max(x));
 b=zeros(a,a);
 
for i=1:d
    if x(i,1)==x(i,2);
       b(x(i,1),x(i,2))=0;
    else
       b(x(i,1),x(i,2))=1;
        b(x(i,2),x(i,1))=1;
    end
end

