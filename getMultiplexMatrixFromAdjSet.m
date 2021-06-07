function[A_nLxnL,N_node,L_net]=getMultiplexMatrixFromAdjSet(AdjSet,pro_jump, NormalizationType, InterLayerMappingType, SetIsolatedNodeSelfLoop)
if ~exist('NormalizationType','var') 
    NormalizationType = 'col';
elseif isempty(NormalizationType) || strcmpi(NormalizationType,'none') 
    NormalizationType = 'None';
end
if ~exist('InterLayerMappingType','var') 
    InterLayerMappingType ='FullMapping';
end
if ~exist('SetIsolatedNodeSelfLoop','var') 
    SetIsolatedNodeSelfLoop = true;  
end
    delta =pro_jump;
    if isa(AdjSet,'struct') || isa(AdjSet,'table')
        if isa(AdjSet,'struct')
            fieldnameset = fieldnames( AdjSet); 
        elseif isa(AdjSet,'table')
            fieldnameset = AdjSet.Properties.VariableNames; 
        end
        L_net  =  length( fieldnameset ); 
        N_node =  length( AdjSet.(fieldnameset{1}));
        NxL = N_node*L_net ; 
        if L_net==1;  delta = 0; end %% no jumping for only single layer. 
        if strcmpi(InterLayerMappingType, 'PartialMapping' )
            k_colvec = []; k_rowvec = [];                
            for ii_net = 1: L_net
                    k_colvec = [k_colvec; sum(AdjSet.(fieldnameset{ii_net}),2) ] ;    
                    k_rowvec = [k_rowvec, sum(AdjSet.(fieldnameset{ii_net}),1) ] ;          
            end  
            A_nLxnL = repmat( speye( N_node, N_node), L_net,L_net); 
            Kcutinter = 0 ; 
            A_nLxnL(k_colvec<=Kcutinter,:  ) = 0;% k==Kcutinter 节点无层间连接  
            A_nLxnL(: , k_rowvec<=Kcutinter) = 0;  
            ind = sub2ind(  size(A_nLxnL) , [1:NxL],[1:NxL] );  
            A_nLxnL(ind) =0;% 对角清零 
            switch NormalizationType
case{'column','col','ProbabilityNormalizationColumn','ProbabilityNormalizationCol'}
                    k_rowvec_interlayer = sum( A_nLxnL,1);%列，等于连接的他层网络数  
                    if delta~=0; A_nLxnL = A_nLxnL.*( delta./sparse(k_rowvec_interlayer +eps) ) ;  end
                    ind_setselfloop = k_rowvec==0; 
                case { 'row' ,'ProbabilityNormalizationRow',    }
                    k_colvec_interlayer = sum( A_nLxnL ,2 );%行，等于连接的他层网络数  
                    size( A_nLxnL )
                    size( k_colvec_interlayer )
                    if delta~=0; A_nLxnL = A_nLxnL.*( delta./sparse(k_colvec_interlayer +eps) ) ;   end 
                    ind_setselfloop = k_colvec==0; 
                case 'LaplacianNormalization' % 
                    k_rowvec_interlayer = sum( A_nLxnL,1);  
                    k_colvec_interlayer = sum( A_nLxnL , 2 );      
                    if delta~=0; A_nLxnL = A_nLxnL./sparse( k_colvec_interlayer.^0.5 +eps )./sparse( k_rowvec_interlayer.^0.5 +eps );  end
                    ind_setselfloop = (k_colvec==0) & (k_rowvec'==0);  
                case 'None'    
                    if delta~=0; A_nLxnL = delta.*A_nLxnL; end
                    ind_setselfloop = (k_colvec==0) & (k_rowvec'==0);  
                otherwise
                    error(['There is no type of normalization: ',char( string(NormalizationType) )] );
            end   
            for ii_net = 1: L_net
                idx = N_node*(ii_net-1)+[1: N_node ] ; 
                if strcmpi(NormalizationType,'None') 
                    A_nLxnL(idx,idx)=  (1-delta).*AdjSet.(fieldnameset{ii_net}); 
                else
                    A_nLxnL(idx,idx)=(1-delta).*getNormalizedMatrix( AdjSet.(fieldnameset{ii_net}) , NormalizationType, false ); 
                end                   
            end
            if SetIsolatedNodeSelfLoop
                ii = find( ind_setselfloop ); 
                ind = sub2ind( size(A_nLxnL), ii,ii ); 
                A_nLxnL(ind) = 1;  
            end        
        elseif strcmpi(InterLayerMappingType, 'FullMapping' )
            A_nLxnL = repmat( (delta/(L_net-1+eps)).*speye( N_node, N_node), L_net,L_net);            
            for ii_net = 1: L_net
                idx = N_node*(ii_net-1)+[1: N_node ] ; 
                if strcmpi(NormalizationType,'None') 
                    A_nLxnL(idx,idx)=  (1-delta).*AdjSet.(fieldnameset{ii_net}); 
                else
                    A_nLxnL(idx,idx)=(1-delta).*getNormalizedMatrix(AdjSet.(fieldnameset{ii_net}),NormalizationType,SetIsolatedNodeSelfLoop ); 
                end                   
            end 
        else
            error(['There is no type of InterLayerType: ',char( string(InterLayerMappingType) )] );
        end

