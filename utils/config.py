reference = "/home/ceglian/codebase/refdata/GRCh38"
build = reference.split("/")[-1]
print("Running for Build {}".format(build))
rho_matrix = "/home/ceglian/codebase/refdata/rho_hgsc.json"
copy_number_data = None
scviz_embedding = None
import os
if not os.path.exists(rho_matrix):
    rho_matrix = "/Users/ceglian/Codebase/refdata/rho_hgsc.json"
