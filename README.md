# Sea4Blue


# short reads
#/media/shared1/atlantic_sea4blue/16S_Illumina/Lapo_18May/RAWREADS


```sh
for r1_file in *R1_001.fastq.gz; do
    # Estrae il nome base (es. "Lapo10_S10" da "Lapo10_S10_L001_R1_001.fastq.gz")
    sample="${r1_file%_L001_R1_001.fastq.gz}"
    
    # Costruisce il nome del file R2 corrispondente
    r2_file="${sample}_L001_R2_001.fastq.gz"
    
    echo "Processing: $sample"
    echo "R1 file: $r1_file"
    echo "R2 file: $r2_file"
    
    # Comando EMU (decommenta quando sei sicuro che funzioni)
     emu abundance --type sr "$r1_file" "$r2_file" \
         --output-dir EMU_db_emu \
        --db /media/shared1/EMU_DB \
         --output-basename "$sample" \
         --threads 10
done
```


# long reads original nanopore 
#/media/shared1/atlantic_sea4blue/16S_nanopore/nanopore_wholeGene/backup_rawreads

```

 for file in *.fastq; do emu abundance $file --output-dir /media/shared1/atlantic_sea4blue/EMU_short_long --db /media/shared1/EMU_DB  --threads 10 ; done
```
# long reads genomiphy 
#/media/shared1/atlantic_sea4blue/genomiphy_targeted/nuovo_basecalling/home/nanopore/Documents/lapo_19_dec_2022_16s_groel/no_sample/20221219_1656_MC-113519_FAS65719_e1c30dd4/basecalled/pass/ALL/16S/rawreads

```

 for file in *.fastq.gz; do emu abundance $file --output-dir /media/shared1/atlantic_sea4blue/EMU_short_long --db /media/shared1/EMU_DB  --threads 10 ; done
```

emu combine-outputs results species
