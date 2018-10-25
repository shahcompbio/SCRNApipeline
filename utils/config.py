reference = "/home/ceglian/codebase/refdata/GRCh38"
build = reference.split("/")[-1]
print("Running for Build {}".format(build))
rho_matrix = "/Users/ceglian/Codebase/refdata/rho_big.json"
copy_number_data = None
scviz_embedding = None
import os
