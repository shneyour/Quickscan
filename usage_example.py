import QuickSCan as qs

# Using quickscan for comparing all the clusters
qs.analyze_sc(data_path='/Users/arielshneyour/Desktop/work/BBB/AD_files/HuVascAD.integrated_DF_withctx_sub_AD.h5ad',
              output_path='/Users/arielshneyour/Desktop/work/immune system/young_patients/quick scan test/',ann_column='cellID')

# 1 vs 1 comparsion
qs.analyze_sc(data_path='/Users/arielshneyour/Desktop/work/BBB/AD_files/HuVascAD.integrated_DF_withctx_sub_AD.h5ad',
              output_path='/Users/arielshneyour/Desktop/work/immune system/young_patients/quick scan test/',
              ann_column='cellID', general=False, cell1='Arterial', cell2='Capillary')
