function[combMatrix,IsNormalized,options_fieldnames]=getNormalizedMatrix_MultiplexHeter(AdjGfGSet,AdjGfD,AdjDfDSet,NormalizationType,NormalizationMode, options_argin) 
istesting = false; 
if nargin <1
    N_gene=30;
    N_disease = 10; 
    AdjGfGSet.Ag = rand(N_gene,N_gene)>0.2 ; 
    AdjGfGSet.Bg = rand(N_gene,N_gene)>0.2 ; 
    AdjGfD = rand(N_gene,N_disease)>0.2 ; 
    AdjDfDSet.Ad = rand(N_disease,N_disease)>0.2 ; 
    AdjDfDSet.Bd = rand(N_disease,N_disease)>0.2 ; 
    P0_G = zeros(N_gene,1); P0_G(1:3)=1; P0_G = P0_G./sum(P0_G);
    P0_D = zeros(N_disease,1); P0_D(1:2)=1; P0_D = P0_D./sum(P0_D);
    warning('TestTestTestTestTestTestTestTestTestTestTestTestTestTestTest'); 
    isdebug =true; 
    istesting =true; 
    NormalizationType =[];
    NormalizationMode =[];  
    options_argin.NormalizationType = 'ProbabilityNormalizationColumn';
    options_argin.SetIsolatedNodeSelfLoop = true ; 
    options_argin.NormalizationMode = 'IndividualNormalization';
end
    if strcmpi(AdjGfGSet,'getoptions') || ~exist('AdjGfGSet','var') || isempty(AdjGfGSet)  || strcmpi(options_argin,'getoptions') 
        options_fieldnames = fieldnames(  options  ) ; 
        combMatrix =[]; IsNormalized =[];
        return; 
    end
    if ~isempty( NormalizationType )
        options_argin.NormalizationType = NormalizationType; 
    end
    if ~isempty( NormalizationMode)
        options_argin.NormalizationMode = NormalizationMode; 
    end  
    [options,n_updated_fields,FieldNameSet_changed, ~ ] = UpdatingOptionsParameters( options,  options_argin, [], options_argin.isdebug  ); 
    NormalizationType = options.NormalizationType ;
    NormalizationMode = options.NormalizationMode ;
    InterLayerMappingType = options.InterLayerMappingType; 
    SetIsolatedNodeSelfLoop = options.SetIsolatedNodeSelfLoop; 
    pro_jump_interlayer_GeneNet     = options.pro_jump_nettype{1,1} ;
    pro_jump_interlayer_DiseaseNet  = options.pro_jump_nettype{2,2} ;
    pro_jump_inter_net_GtoD         = options.pro_jump_nettype{1,2} ;
    pro_jump_inter_net_DtoG         = options.pro_jump_nettype{2,1} ;
    isdebug     =options.isdebug ; 
    if isempty( AdjDfDSet )
       AdjDfDSet.DisMatrixEYE = speye(size( AdjGfD,2 ));  
       warning('AdjDfDSet is empty.');
    end  
    if any( strcmpi(NormalizationMode, {'IndividualNormalization', 'Individual'})) %逐一标准化
        IndividualAdjNormalizationType = NormalizationType;  
        IsUnifiedNormalization         = false;  
    elseif any( strcmpi(NormalizationMode, {'UnifiedNormalization', 'Unified'})) %统一标准化 
        % 该模式先用'None'方式组合各个矩阵，然后统一标准化 
        IndividualAdjNormalizationType = 'None'; 
        IsUnifiedNormalization         = true;  
    else
        error( [ NormalizationMode, ':NormalizationMode has No definition.']);
    end
[M_GfG,N_gene,L_genenet]=getMultiplexMatrixFromAdjSet(AdjGfGSet,pro_jump_interlayer_GeneNet,IndividualAdjNormalizationType, InterLayerMappingType, SetIsolatedNodeSelfLoop) ; 
[M_DfD,N_disease,L_disnet]=getMultiplexMatrixFromAdjSet(AdjDfDSet,pro_jump_interlayer_DiseaseNet,IndividualAdjNormalizationType,InterLayerMappingType,SetIsolatedNodeSelfLoop) ; 
    AdjGfD = repmat(AdjGfD, L_genenet, L_disnet);  
    pro_jump_GtoD = pro_jump_inter_net_GtoD   ;  
    pro_jump_DtoG = pro_jump_inter_net_DtoG  ;      
    IsNormalized = true; 
    switch lower( NormalizationType) % NormalizationTyp确定哪种组合方式以及标准化方式   
        case lower( {'None'} )
            combMatrix = [ M_GfG,    AdjGfD; ...
                           AdjGfD',  M_DfD      ] ;
            IsNormalized = false;
        case lower( {'Weight'} )
            pro_jump = (pro_jump_GtoD+pro_jump_DtoG)/2; 
            combMatrix = [ (1-pro_jump).*M_GfG,    pro_jump.*AdjGfD; ...
                           pro_jump.*AdjGfD',      (1-pro_jump).*M_DfD    ] ;
            IsNormalized = false;
        Case lower({'col','ProbabilityNormalizationColumn','NormalizationColumn', 'Column'}) 
            idxDis_WithDiseaseGene =  sum( AdjGfD, 1)~=0;
            idxGene_WithDisease    = (sum( AdjGfD, 2)~=0)';
            M_GfD = getNormalizedMatrix(AdjGfD   , IndividualAdjNormalizationType, false );  % probabilities from disease space to gene space 
            M_DfG = getNormalizedMatrix(AdjGfD'  , IndividualAdjNormalizationType, false );  % probabilities from gene space to disease space
            M_GfG(:,idxGene_WithDisease)=(1pro_jump_GtoD).*M_GfG(:,idxGene_WithDisease); 
            M_DfD(:,idxDis_WithDiseaseGene)=(1-pro_jump_DtoG).*M_DfD(:,idxDis_WithDiseaseGene ) ; 
            M_GfD                           = pro_jump_DtoG.*M_GfD;
            M_DfG                           = pro_jump_GtoD.*M_DfG;
            combMatrix = [ M_GfG, M_GfD; ...
                           M_DfG, M_DfD    ] ; 
        case lower( {'row','ProbabilityNormalizationRow','NormalizationRow'} )
            idxDis_WithDiseaseGene = (sum( AdjGfD, 1)~=0);  
            idxGene_WithDisease    = (sum( AdjGfD, 2)~=0);   
            M_GfD =getNormalizedMatrix(AdjGfD , IndividualAdjNormalizationType, false );  
            M_DfG = getNormalizedMatrix(AdjGfD', IndividualAdjNormalizationType, false );  
            M_GfG(idxGene_WithDisease,:)=(1-pro_jump_GtoD).*M_GfG(idxGene_WithDisease,:) ; 
            M_DfD(idxDis_WithDiseaseGene,:)   = (1-pro_jump_DtoG).*M_DfD(idxDis_WithDiseaseGene,: ) ; 
            M_GfD = pro_jump_GtoD.*M_GfD; 
            M_DfG = pro_jump_DtoG.*M_DfG;  
combMatrix = [ M_GfG, M_GfD; ...
                           M_DfG, M_DfD    ] ;            
        case lower( {'LaplacianNormalization'} )
            M_GfD=getNormalizedMatrix(AdjGfD,IndividualAdjNormalizationType, false );  % probabilities from disease space to gene space 
            M_DfG = M_GfD';  % probabilities from gene space to disease space
            pro_jump = (pro_jump_GtoD+pro_jump_DtoG)/2; 
            combMatrix = [ (1-pro_jump).*M_GfG, pro_jump.*M_GfD; pro_jump.*M_DfG , (1-pro_jump).*M_DfD    ] ;        
        otherwise
            error( [ IndividualAdjNormalizationType, ':IndividualAdjNormalizationType has No definition.']);
    end 
    if IsUnifiedNormalization && ~any( strcmpi( NormalizationType,  {'None','Weight'}  ) ) 
        combMatrix = getNormalizedMatrix(combMatrix   , NormalizationType, true ); 
        IsNormalized = true; 
    end 
if  istesting  
    combMatrix
     sumcol = sum(combMatrix,1)
     sumcol = sum(combMatrix,2)
     combMatrix =[];
end
if isdebug  
    disp('*******************************************************') 
end
end
