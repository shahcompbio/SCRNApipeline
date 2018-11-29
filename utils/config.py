reference = "/igo_large/reference/GRCh38"
build = reference.split("/")[-1]
#build = "hg19"
print("Running for Build {}".format(build))
rho_matrix = "/home/nceglia/codebase/refdata/rho.json"
copy_number_data = None
scviz_embedding = None
import os
if not os.path.exists(rho_matrix):
    rho_matrix = "/Users/ceglian/Codebase/refdata/rho.json"
perplexes = [5,10,20,30]
mdistances = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
components = [10,20,50]
chemistry = "SC5P-R2"
