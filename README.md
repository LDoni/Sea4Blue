# Sea4Blue


## rifaccio prima fastp


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

## con fastp
 

for r1_file in *R1_001.fastq.gz; do
    # Estrae il nome base (es. "Lapo10_S10" da "Lapo10_S10_L001_R1_001.fastq.gz")
    sample="${r1_file%_L001_R1_001.fastq.gz}"

    # Costruisce il nome del file R2 corrispondente
    r2_file="${sample}_L001_R2_001.fastq.gz"

    echo "Processing: $sample"
    echo "R1 file: $r1_file"
    echo "R2 file: $r2_file"

   # fastp -i "$r1_file" -I "$r2_file" -o fastp/"${sample}_trimmed_R1.fastq.gz" -O fastp/"${sample}_trimmed_R2.fastq.gz"

    # Esegui EMU (decommenta la riga seguente se desideri eseguire il comando)
    emu abundance --type sr \
        "fastp/${sample}_trimmed_R1.fastq.gz" "fastp/${sample}_trimmed_R2.fastq.gz" \
        --output-dir /media/shared1/atlantic_sea4blue/EMU_short_long/fastp_denoised \
        --db /media/shared1/EMU_DB \
        --output-basename "$sample" \
        --threads 15

done






```


# long reads original nanopore 



#/media/shared1/atlantic_sea4blue/16S_nanopore/nanopore_wholeGene/backup_rawreads

```
mkdir -p fastp

for file in *.fastq; do
    # Ottieni il nome base senza estensione
    base=$(basename "$file" .fastq)
    
    # Definisci il file di output
    outfile="fastp/${base}_trimmed.fastq"

    echo "Processing $file → $outfile"
    
    fastplong -i "$file" -o "$outfile"
done

 for file in *.fastq; do emu abundance $file --output-dir /media/shared1/atlantic_sea4blue/EMU_short_long --db /media/shared1/EMU_DB  --threads 10 ; done
```
# long reads genomiphy 
#/media/shared1/atlantic_sea4blue/genomiphy_targeted/nuovo_basecalling/home/nanopore/Documents/lapo_19_dec_2022_16s_groel/no_sample/20221219_1656_MC-113519_FAS65719_e1c30dd4/basecalled/pass/ALL/16S/rawreads

```
mkdir -p fastp

for file in *.fastq.gz; do
    # Rimuove .fastq.gz per ottenere il nome base (es. barcode09)
    base=$(basename "$file" .fastq.gz)
    
    # Definisce il file di output
    outfile="fastp/${base}_trimmed.fastq.gz"


    echo "Processing $file → $outfile"
    
    fastplong -i "$file" -o "$outfile" --thread 10
done
 for file in *.fastq.gz; do emu abundance $file --output-dir /media/shared1/atlantic_sea4blue/EMU_short_long --db /media/shared1/EMU_DB  --threads 10 ; done
```

emu combine-outputs results species
