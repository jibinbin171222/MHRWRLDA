function [kd,km] = gaussiansimilarity(interaction,nd,nm)
 gamad = nd/(norm(interaction,'fro')^2);
%calculate Gaussian kernel for the similarity between disease: kd
C=interaction;
kd=zeros(nd,nd);
D=C*C';
for i=1:nd
    for j=i:nd
        kd(i,j)=exp(-gamad*(D(i,i)+D(j,j)-2*D(i,j)));
    end
end
kd=kd+kd'-diag(diag(kd));
gamam = nm/(norm(interaction,'fro')^2);
%calculate Gaussian kernel for the similarity between miRNA: km
km=zeros(nm,nm);
E=C'*C;
for i=1:nm
    for j=i:nm
        km(i,j)=exp(-gamam*(E(i,i)+E(j,j)-2*E(i,j)));
    end
end
km=km+km'-diag(diag(km));
end

%extractVariablesFromFieldnames.m
function [Var_out] = extractVariablesFromFieldnames(   Var_in, fieldnameSet )  
if nargin<1
  Var_in.q=[1:5]';
  Var_in.ww='ddfsdfs';
  Var_in.qqq=1111;
  fieldnameSet = fieldnames(Var_in);  
  fieldnameSet = {'q','qqq','q'};
  Var_in 
end
if isa(Var_in,'struct')
    tableVar_in = struct2table(Var_in,  'AsArray',true)   ; 
    Var_out = tableVar_in(:, fieldnameSet ); 
    Var_out = table2struct(Var_out);    
elseif isa(Var_in,'table')
    Var_out = Var_in(:, fieldnameSet ); 
end
end
