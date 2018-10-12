reference = "/juno/work/shah/ceglian/refdata/human/GRCh38"
build = reference.split("/")[-1]
print("Running for Build {}".format(build))
rho_matrix = "utils/rho_matrix.txt"
cnv_data = "utils/cnv.txt"
