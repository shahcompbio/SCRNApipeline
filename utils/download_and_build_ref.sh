wget ftp://ftp.ensembl.org/pub/grch37/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz


wget ftp://ftp.ensembl.org/pub/grch37/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh37.82.gtf.gz
gunzip Homo_sapiens.GRCh37.82.gtf.gz


wget ftp://ftp.ensembl.org/pub/release-84/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz


wget ftp://ftp.ensembl.org/pub/release-84/gtf/mus_musculus/Mus_musculus.GRCm38.84.gtf.gz
gunzip Mus_musculus.GRCm38.84.gtf.gz


cellranger mkgtf Homo_sapiens.GRCh37.82.gtf Homo_sapiens.GRCh37.82.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense


cellranger mkgtf Mus_musculus.GRCm38.84.gtf Mus_musculus.GRCm38.84.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense


cellranger mkref --genome=hg19 \
                 --fasta=Homo_sapiens.GRCh37.dna.primary_assembly.fa \
                 --genes=Homo_sapiens.GRCh37.82.filtered.gtf \
                 --genome=mm10 \
                 --fasta=Mus_musculus.GRCm38.dna.primary_assembly.fa \
                 --genes=Mus_musculus.GRCm38.84.filtered.gtf \
                 --ref-version=1.2.0
