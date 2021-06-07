%caculate lncRNA functional similarity
function lncFunsim=test(MDA,DS)
[m,n]=size(MDA);
R=zeros(n,n);
for i=1:n
    for j=1:n
    sum1=0;n1=0;
    for q=1:m
        if MDA(q,i)==1
            max1=0;n1=n1+1;
            for p=1:m
                if MDA(p,j)==1&&DS(q,p)>max1
                    max1=DS(q,p);
                end
            end
                sum1=sum1+max1;
        end
    end
        sum2=0;n2=0;
        for p=1:m
            if MDA(p,j)==1
                max2=0;n2=n2+1;
                for q=1:m
                    if MDA(q,i)==1&&DS(q,p)>max2
                        max2=DS(q,p);
                    end
                end
                sum2=sum2+max2;
            end
        end
        R(i,j)=(sum1+sum2)/(n1+n2);
    end
end
for i=1:n
    for j=1:n
        if R(i,j)>0.1
            R(i,j)=R(i,j);
        else
            R(i,j)=0;
        end
    end
end
lncFunsim=R;
end
