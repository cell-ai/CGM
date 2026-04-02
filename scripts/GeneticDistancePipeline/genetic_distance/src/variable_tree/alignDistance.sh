#!/bin/bash

# This script receives FASTA files as input and performs multiple sequence alignment with Mafft (--auto) and removes gaps and poorly aligned regions with TRIMAL

# Next, the search for the best evolutionary model is carried out using MODELFINDER from IQTREE and the calculation of genetic distances between pairs of sequences with RAxML using the evolutionary model indicated by MODELFINDER. 

# Esse script recebe como input arquivos FASTA e realiza o alinhamento múltiplo de sequências com Mafft (--auto) e remove gaps e regiões mal alinhadas com TRIMAL. 

# Em seguida, é feita a busca pelo melhor modelo evolutivo utilizando MODELFINDER de IQTREE e o cálculo de distâncias genéticas entre pares de sequências com RAxML utilizando o modelo evolutivo indicado por MODELFINDER. 
THREADS=32

align_sequences() {
    local INPUT_DIR=$1
    local OUTPUT_DIR=$2
    local LOG_FILE="$OUTPUT_DIR/error_log.txt"

    if [ ! -d "$INPUT_DIR" ]; then
        echo "Erro: diretório de entrada $INPUT_DIR não encontrado!"
        exit 1
    fi

    for FILE in $INPUT_DIR/*.fasta; do
        #if [[ $FILE =~ .*_filtered\.fasta$ ]] && ! [[ $FILE =~ ^\./\. ]]; then
        if [[ $FILE =~ .*\.fasta$ ]] && ! [[ $FILE =~ ^\./\. ]]; then
            #local BASENAME=$(basename -- "$FILE" .fasta | sed 's/_filtered//')
            local BASENAME=$(basename -- "$FILE" .fasta)
            
            if [[ $BASENAME == .* ]]; then
                continue
            fi
        
            #local MAFFT_OUTPUT="$OUTPUT_DIR/filtered_${BASENAME}-aligned.fasta"
            #local TRIMAL_OUTPUT="$OUTPUT_DIR/filtered_${BASENAME}-aligned-trim.fasta"
            local MAFFT_OUTPUT="$OUTPUT_DIR/${BASENAME}-aligned.fasta"
            local TRIMAL_OUTPUT="$OUTPUT_DIR/${BASENAME}-aligned-trim.fasta"

            mafft --auto --thread $THREADS $FILE > $MAFFT_OUTPUT
            if [ $? -ne 0 ]; then
                echo "Erro no MAFFT processando $FILE" >> $LOG_FILE
                continue
            fi

            trimal -in $MAFFT_OUTPUT -out $TRIMAL_OUTPUT -automated1
            if [ $? -ne 0 ]; then
                echo "Erro no TRIMAL processando $MAFFT_OUTPUT" >> $LOG_FILE
                continue
            fi

            if [ ! -f $TRIMAL_OUTPUT ]; then
                echo "Erro: arquivo $TRIMAL_OUTPUT não encontrado após TRIMAL" >> $LOG_FILE
                continue
            fi

            calculate_distances $TRIMAL_OUTPUT $BASENAME $OUTPUT_DIR $LOG_FILE
        fi
    done 
}

calculate_distances() {
    local INPUT=$1
    local BASENAME=$2
    local OUTPUT_DIR=$3
    local LOG_FILE=$4

    if [ ! -f "$INPUT" ]; then
        echo "Erro: O arquivo $INPUT não foi encontrado!" >> $LOG_FILE
        return
    fi

    raxmlHPC -p $(date +%s%N) -s $INPUT -m PROTGAMMAAUTO -f x -n ${BASENAME}_matrix -w $OUTPUT_DIR -T $THREADS
    if [ $? -ne 0 ]; then
        echo "Erro no RAxML processando $INPUT" >> $LOG_FILE
    fi
}

align_sequences "$1" "$2"
echo "Processo concluído."