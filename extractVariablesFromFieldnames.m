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
