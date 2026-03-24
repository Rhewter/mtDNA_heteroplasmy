# mtDNA_heteroplasmy

Pipeline para análise de heteroplasmia em mtDNA a partir de dados de sequenciamento do SRA-NCBI.

## Descrição

Este repositório contém um script Bash (`run_sra_novoplasty_heteroplasmy.sh`) que automatiza o seguinte pipeline:

1. **Download de dados SRA**: Usando `prefetch` e `fasterq-dump` para baixar dados de sequenciamento.
2. **Controle de qualidade**: Trimmomatic com parâmetros configuráveis (padrão: `SLIDINGWINDOW:4:15`).
3. **Montagem mtDNA**: NOVOPlasty para assemblear o genoma mitocondrial.
4. **Análise de heteroplasmia**: Mapeamento com BWA, chamada de variantes com bcftools, e extração de heteroplasmias.

O script inclui checkpoints para evitar reprocessamento desnecessário e parâmetros configuráveis.

## Dependências

Certifique-se de que as seguintes ferramentas estão instaladas e disponíveis no PATH:

- `prefetch` e `fasterq-dump` (do pacote sra-tools)
- `trimmomatic` (jar ou comando)
- `perl` (para NOVOPlasty)
- `bwa`, `samtools`, `bcftools`
- NOVOPlasty.pl (script Perl)

## Como usar

### Preparação

1. Crie um arquivo de texto com os IDs SRA, um por linha (exemplo: `ids.txt`):
   ```
   SRR31925970
   SRR31925971
   ```

2. Execute o script:
   ```bash
   bash run_sra_novoplasty_heteroplasmy.sh -i ids.txt
   ```

### Parâmetros

- `-i FILE`: Arquivo com lista de IDs SRA (obrigatório).
- `-o DIR`: Diretório de saída (padrão: `./mtDNA_pipeline_out`).
- `-t N`: Número de threads (padrão: 8).
- `-n PATH`: Caminho para o jar do Trimmomatic (padrão: `trimmomatic` no PATH).
- `-a PATH`: Arquivo de adaptadores (padrão: `TruSeq3-PE.fa`).
- `-p PATH`: Caminho para o script NOVOPlasty.pl (padrão: `NOVOPlasty.pl` no PATH).
- `-r PATH`: Referência mtDNA para heteroplasmia (padrão: usa a assembleia do NOVOPlasty).
- `--trim-params STR`: Parâmetros de trimming (padrão: `SLIDINGWINDOW:4:15 MINLEN:36`).
- `--step NAME`: Executar a partir de um passo: `download`, `trim`, `novo`, `hetero`, `all` (padrão: `all`).
- `-h|--help`: Mostra ajuda.

### Exemplos

- Rodar apenas download:
  ```bash
  bash run_sra_novoplasty_heteroplasmy.sh -i ids.txt --step download
  ```

- Com parâmetros customizados:
  ```bash
  bash run_sra_novoplasty_heteroplasmy.sh -i ids.txt -o /path/to/output -t 16 --trim-params "SLIDINGWINDOW:4:20 MINLEN:50"
  ```

### Saídas

- `01_raw/`: Dados brutos baixados.
- `02_trimmed/`: Dados após trimming.
- `03_novoplasty/`: Assembleias mtDNA e logs.
- `04_heteroplasmy/`: Arquivos BAM, VCF e tabela de heteroplasmias (`*_heteroplasmy.tsv`).

## Observações

- O script assume dados paired-end.
- Para heteroplasmia, filtra variantes com MAF > 0.01 e DP > 100.
- Se precisar de referência externa, use `-r`.
- Ajustes no NOVOPlasty podem ser feitos editando o arquivo de config gerado por amostra.

