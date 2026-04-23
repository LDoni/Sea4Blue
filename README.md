# Sea4Blue Atlantic 16S Workspace

Repository di lavoro per le analisi 16S del transect atlantico.

## Main entry points

- `dada2.R`: preprocessing e denoising delle letture 16S.
- `scripts/run_final_analyses.R`: entry point del workflow finale ripulito.
- `scripts/`: script modulari usati per alpha/beta diversity, distance-decay, assembly, moduli, TINA/PINA e oceanografia.
- `original_workflow/`: workflow originale da cui sono stati derivati i moduli e le analisi di rete.
- `quarto/quarto_nb_sea4blue.qmd`: notebook Quarto usato come base narrativa/riproducibile.

## Repository structure

- `scripts/`: workflow R principale e funzioni di supporto.
- `original_workflow/`: script storici/originali del workflow modulare.
- `quarto/`: notebook Quarto.
- root: script esplorativi e helper ancora usati nel progetto.

Il repository e stato popolato con gli script effettivamente usati per le analisi e con il materiale Quarto rilevante. Dati grezzi, output pesanti e file temporanei restano esclusi dal versionamento.
