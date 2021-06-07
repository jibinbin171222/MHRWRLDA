function [L_matrix,N_row, N_col, fieldnameset, SetType] = getMatrixSetInfo( MatrixSet  )
    L_matrix =0; N_row =[]; N_col =[]; fieldnameset =[]; SetType =[]; 
    if isa(MatrixSet,'struct')
        SetType      = 'struct'; 
        fieldnameset = fieldnames( MatrixSet); 
        L_matrix     =  length( fieldnameset ); 
        if L_matrix>0
            [N_row, N_col ] =  size(  MatrixSet.(fieldnameset{1}) ); 
        end       
        return; 
    elseif isa(MatrixSet,'table')
        SetType      = 'table'; 
        fieldnameset = MatrixSet.Properties.VariableNames; 
        L_matrix     =  length( fieldnameset ); 
        if L_matrix>0
            [N_row, N_col ] =  size(  MatrixSet.(fieldnameset{1}) ); 
        end
        return; 
    elseif isa(MatrixSet,'cell')
        SetType      = 'cell';         
        L_matrix     = length(  MatrixSet ); 
        fieldnameset = [];        
        if L_matrix>0
            [N_row, N_col ] =  size(  MatrixSet{1} );             
        end
        return; 
    else        
        error(  'There is no definition of this MatrixSet.' ); 
    end
end
