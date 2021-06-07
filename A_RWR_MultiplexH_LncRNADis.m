function [P_G_set, P_G_m, P_D_set, P_D_m, TableScores, Mall] = A_RWR_MultiplexH_LncRNADis(AdjGfGSet, AdjGfD, AdjDfDSet, P0_G,P0_D, options_argin, UsedGlobalVar) 
istest = false ;    
    options.NormalizationType = options_argin.NormalizationType ; %  'col',  'row'  
    options.NormalizationMode=options_argin.NormalizationMode;   
    % default values of options 
    options.restart       = 0.9; 
    options.InterLayerMappingType='FullMapping';%InterLayerMappingType= 'PartialMapping';  
    options.SetIsolatedNodeSelfLoop= true;   
    options.pro_jump_nettype  ={ 0.1,   0.9; ...
                              0.9,   0.1    };
    options.iniProRatioVec_interlayers_net1 =[] ; % factor of assigning initial probabilities in same-type gene networks 
    options.iniProRatioVec_interlayers_net2 =[];    % factor of assigning initial probabilities in same-type disease networks
    options.iniProRatio_inter_net2 = 0.9 ;  % factor of assigning initial probabilities in different-type networks    
    options.isdebug                 = false;  
    [options, n_updated_fields, FieldNameSet_changed, ~ ] = UpdatingOptionsParameters( options,  options_argin,[],false   ); 
    NormalizationType = options.NormalizationType ; %  'col',  'row'  
    NormalizationMode = options.NormalizationMode ;    
    restart = options.restart;    
    InterLayerMappingType = options.InterLayerMappingType; 
    SetIsolatedNodeSelfLoop = options.SetIsolatedNodeSelfLoop; 
    pro_jump_interlayer_GeneNet     = options.pro_jump_nettype{1,1} ;
    pro_jump_interlayer_DiseaseNet  = options.pro_jump_nettype{2,2} ;
    pro_jump_inter_net_GtoD        = options.pro_jump_nettype{1,2} ; 
    pro_jump_inter_net_DtoG         = options.pro_jump_nettype{2,1} ;                     
    iniProRatioVec_interlayers_GeneNet=options.iniProRatioVec_interlayers_net1;
iniProRatioVec_interlayers_DiseaseNet =options.iniProRatioVec_interlayers_net2 ; eta_iniProRatio_DiseaseNet = options.iniProRatio_inter_net2 ;      
    isdebug     =options.isdebug ; 
    global   Global_Var_RWR_MultiplexH 
    if ~exist('UsedGlobalVar','var') || isempty (UsedGlobalVar)
        UsedGlobalVar = false; 
    end   
    if not( UsedGlobalVar )  || isempty( Global_Var_RWR_MultiplexH )           
        [Global_Var_RWR_MultiplexH.NormalizedMatrix,IsNormalized]=getNormalizedMatrix_MultiplexHeter(AdjGfGSet,AdjGfD,AdjDfDSet, NormalizationType,NormalizationMode,  options) ; 
        Global_Var_RWR_MultiplexH.IsNormalized = true;  
    end
    Mall = Global_Var_RWR_MultiplexH.NormalizedMatrix;
    IsNormalized = true;  
If any( strcmpi(NormalizationType,{'col','ProbabilityNormalizationColumn',
'NormalizationColumn', 'Column'}) )
        P0_G = P0_G./(sum(P0_G,1)+eps);  
        P0_D = P0_D./(sum(P0_D,1)+eps);  
    end 
    [L_genenet,N_genenode, ~, GeneNetNameSet] = getMatrixSetInfo( AdjGfGSet  ) ; 
    [L_disnet,N_disnode, ~, DisNetNameSet] = getMatrixSetInfo( AdjDfDSet  ) ;  
    if isempty(iniProRatioVec_interlayers_GeneNet)
        iniProRatioVec_interlayers_GeneNet = repmat(1/L_genenet,L_genenet,1); 
    end             
    P0_G_L_genenet = kron( reshape(iniProRatioVec_interlayers_GeneNet,[],1),  (P0_G ) ) ;  
    if isempty(iniProRatioVec_interlayers_DiseaseNet)
        iniProRatioVec_interlayers_DiseaseNet = repmat(1/L_disnet,L_disnet,1); 
    end          
    P0_D_L_disnet  =  kron( reshape(iniProRatioVec_interlayers_DiseaseNet,[],1), (P0_D ) ) ;  
    P0_all=[(1-eta_iniProRatio_DiseaseNet)*P0_G_L_genenet; eta_iniProRatio_DiseaseNet*P0_D_L_disnet] ; 
    Pt_all = A_RWRplus(Mall, restart, P0_all , [],[], IsNormalized);   
	N_gene = N_genenode;  N_disease =N_disnode ;   
    n_column = size(Pt_all,2); 
    P_G_set = Pt_all(1:N_gene*L_genenet, :) ;
    P_D_set = Pt_all(N_gene*L_genenet+1:end, :) ; 
    P_G_set_cell = mat2cell(P_G_set, repmat(N_gene,L_genenet,1), n_column  ); 
    P_G_set = zeros(N_gene,n_column , L_genenet );
    for i_cell = 1: L_genenet
       P_G_set(:,:,i_cell) = P_G_set_cell{i_cell};          
    end    
    P_D_set_cell = mat2cell(P_D_set, repmat(N_disease,L_disnet,1), n_column  ); 
    P_D_set = zeros( N_disease,n_column  , L_disnet  );
    for i_cell = 1: L_disnet
       P_D_set(:,:,i_cell) = P_D_set_cell{i_cell};          
    end        
    P_G_m.geomean = geomean(P_G_set,3) ; 
    P_G_m.mean    = mean(P_G_set,3) ; 
    P_D_m.geomean = geomean(P_D_set,3) ; 
    P_D_m.mean    = mean(P_D_set, 3 ) ; 
    if isa(AdjGfGSet,'struct')
        AdjSetFieldnameset = fieldnames( AdjGfGSet);  
    elseif isa(AdjGfGSet,'table')
        AdjSetFieldnameset = AdjGfGSet.Properties.VariableNames; 
    end  
    TableScores = table; 
    if L_genenet>1   
        TableScores.(['ScoreMean' ] ) = mean( P_G_set , 3 ) ;   
        TableScores.(['ScoreGeoMean' ] ) = geomean( P_G_set , 3  ) ; 
        Ranks = getRankingOfScoreList(   P_G_set, 'descend' ) ;
        TableScores.(['RankMean' ] ) = -mean( Ranks , 3  ) ; 
        TableScores.(['RankGeoMean'] ) = -geomean( Ranks , 3 ) ;     
        TableScores.(['RankHarmMean'] )= -harmmean( Ranks,3) ;   
    end
    for ii_L_net = 1: L_genenet
        TableScores.([AdjSetFieldnameset{ii_L_net}] ) = P_G_set(:,:, ii_L_net) ;    
    end 
    if isa(AdjDfDSet,'struct')
        AdjSetFieldnameset = fieldnames( AdjDfDSet);  
    elseif isa(AdjDfDSet,'table')
        AdjSetFieldnameset = AdjDfDSet.Properties.VariableNames; 
    end  
    TableDisScores = table; 
    if L_disnet>1   
        TableDisScores.(['ScoreMean' ] ) = mean( P_D_set , 3  ) ;   
        TableDisScores.(['ScoreGeoMean' ] ) = geomean( P_D_set ,3  ) ;    
        Ranks = getRankingOfScoreList( P_D_set, 'descend' ) ;
        TableDisScores.(['RankMean' ] ) = -mean( Ranks , 3  ) ;   
        TableDisScores.(['RankGeoMean'] ) = -geomean( Ranks , 3 ) ;  
        TableDisScores.(['RankHarmMean'] ) = -harmmean( Ranks , 3 ) ;     
    end
    for ii_L_net = 1: L_disnet
        TableDisScores.([AdjSetFieldnameset{ii_L_net}] ) = P_D_set(:, :, ii_L_net) ;    
    end 
    if istest
        sum(full( Mall ))  
        allP=sum(Pt_all)
        ss = sum(Mall(:))  
        size(Mall) 
        pg = sum( P_G_set(:) )
        pd = sum( P_D_set(:)   )
        pp =pg + pd
    end  
