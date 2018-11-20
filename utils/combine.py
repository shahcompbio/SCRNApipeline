from interface.singlecellexperiment import SingleCellExperiment

def combine(sces, filename):
    combined_sce = SingleCellExperiment.fromRData(sces[0])
    for sce in sces[1:]:
        combined_sce = combined_sce + sce
    combined_sce.save(combined_sce,filename)
