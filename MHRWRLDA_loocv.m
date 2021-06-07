clear;
clc;
AdjDfDSet.A=load('diseasesim_DOSEwang_matrix.mat');
AdjGfD1=load('lncRNA_Disease_Matrix.mat');
DS=load('diseasesim_DOSEwang_matrix.txt');
options_argin.NormalizationType = 'ProbabilityNormalizationColumn';   
options_argin.NormalizationMode = 'IndividualNormalization';
UsedGlobalVar = false;
[M_lncRNA, N_disease] =size(AdjGfD1);
AUClist =nan(M_lncRNA, N_disease);
AUPRlist =nan(M_lncRNA, N_disease);
for j_dis=1:N_disease
    label=AdjGfD1(:,j_dis);
    P0_D=zeros(N_disease,1); 
    P0_D(j_dis) = 1 ; 
    Vp=find(AdjGfD1(:,j_dis)==1); 
    Vn=find(AdjGfD1(:,j_dis)==0); 
    n_pos = numel(Vp); 
    n_neg = numel(Vn); 
    for i_pos = 1: n_pos
        id_pos_test = Vp( i_pos ); 
        ind_test = [Vn; id_pos_test];
        label=AdjGfD1(:,j_dis);
        label( id_pos_test )=0;
        P0_G=label;
        matDT=AdjGfD1;
        matDT(:,j_dis)=label;
        AdjGfGSet.A= test(matDT',DS);
        [kd,km] = gaussiansimilarity(matDT',N_disease,M_lncRNA);
        AdjDfDSet.B=kd;
        AdjGfGSet.B=km;
        [P_G_set, P_G_m, P_D_set, P_D_m, TableScores, Mall]= A_RWR_MultiplexH_LncRNADis(AdjGfGSet, matDT, AdjDfDSet, P0_G, P0_D, options_argin, UsedGlobalVar);    
        MatPredict = P_G_m.mean;
        [AUC,AUPR,Acc,Sen,Spe,Pre]=ROCcompute(MatPredict,AdjGfD1(:,j_dis),1);
        AUClist(i_pos,j_dis)=AUC;
        AUPRlist(i_pos,j_dis)=AUPR;
    end
end
Res =[];
Res.AUCmean_all  =mean(AUClist(:),'omitnan');
Res.AUPRmean_all =mean(AUPRlist(:),'omitnan');
Res.AUCmean_dis  =mean(AUClist ,1, 'omitnan');
Res.AUPRmean_dis =mean(AUPRlist ,1 ,  'omitnan');
Res
