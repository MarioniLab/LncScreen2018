# Annotation sources

Homo_sapiens.GRCh38.91.gtf.gz was obtained from http://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/ and processed as described below.

```sh
zcat Homo_sapiens.GRCh38.91.gtf.gz | \
    sed -r "s/^([0-9MXY])/chr\1/" | \
    sed "s/^chrMT/chrM/g" | \
    awk '$3 == "exon"' > hg38.gtf
```
