#!/bin/bash
# filepath: /home/davi/CAD/atv21/test_scaling.sh

#SBATCH --partition gpu-4-a100
#SBATCH --gpus-per-node=1
#SBATCH --nodes=4               # Máximo de nós que você quer testar
#SBATCH --time=01:00:00         # Tempo total para todos os testes
#SBATCH --job-name=chapa_scaling
#SBATCH --output=scaling_results-%j.out

cd $SLURM_SUBMIT_DIR
ulimit -s unlimited

# Compilar o programa
echo "Compilando o programa..."
mpicc -O3 -march=native -fopenmp chapa.c -o chapa -lm
echo "Compilação concluída."
echo "================================================"

# Arrays com configurações para testar
NODES_ARRAY=(1 2 4)
PROCESSES_PER_NODE_ARRAY=(2 4 8)

# Arquivo para salvar resultados
RESULTS_FILE="scaling_results_$(date +%Y%m%d_%H%M%S).txt"

# Cabeçalho do arquivo de resultados
echo "Nodes,Processes_per_node,Total_processes,Execution_time,Speedup" > $RESULTS_FILE

# Configuração base para cálculo de speedup (1 nó, 1 processo)
BASE_TIME=""

echo "Iniciando testes de escalabilidade..."
echo "Resultados serão salvos em: $RESULTS_FILE"
echo "================================================"

# Loop para diferentes números de nós
for NODES in "${NODES_ARRAY[@]}"; do
    # Loop para diferentes números de processos por nó
    for PROCS_PER_NODE in "${PROCESSES_PER_NODE_ARRAY[@]}"; do
        TOTAL_PROCS=$((PROCS_PER_NODE))
        
        # Pular configurações que excedem recursos disponíveis
        # if [ $TOTAL_PROCS -gt 32 ]; then
        #     echo "Pulando configuração: $NODES nós, $PROCS_PER_NODE proc/nó (total: $TOTAL_PROCS - muito alto)"
        #     continue
        # fi
        
        echo "Testando: $NODES nós, $PROCS_PER_NODE processos/nó (total: $TOTAL_PROCS processos)"
        
        # # Configurar número de threads OpenMP baseado no número de processos por nó
        # if [ $PROCS_PER_NODE -eq 1 ]; then
        #     export OMP_NUM_THREADS=4
        # elif [ $PROCS_PER_NODE -eq 2 ]; then
        #     export OMP_NUM_THREADS=2
        # else
        #     export OMP_NUM_THREADS=1
        # fi
        
        # echo "  - OMP_NUM_THREADS: $OMP_NUM_THREADS"
        
        # Executar o programa e capturar o tempo
        START_TIME=$(date +%s.%N)
        
        if [ $NODES -eq 1 ]; then
            # Para 1 nó, usar mpirun simples
            mpirun -np $TOTAL_PROCS ./chapa > temp_output_${NODES}_${PROCS_PER_NODE}.txt 2>&1
        else
            # Para múltiplos nós, especificar distribuição
            mpirun --map-by ppr:${PROCS_PER_NODE}:node -np $TOTAL_PROCS ./chapa > temp_output_${NODES}_${PROCS_PER_NODE}.txt 2>&1
        fi
        
        END_TIME=$(date +%s.%N)
        EXECUTION_TIME=$(echo "$END_TIME - $START_TIME" | bc -l)
        
        # Extrair tempo do programa (se disponível na saída)
        PROGRAM_TIME=$(grep "Tempo de execução" temp_output_${NODES}_${PROCS_PER_NODE}.txt | awk '{print $NF}' | sed 's/segundos//')
        
        if [ -z "$PROGRAM_TIME" ]; then
            PROGRAM_TIME=$EXECUTION_TIME
        fi
        
        # Calcular speedup (baseado na primeira execução)
        if [ -z "$BASE_TIME" ]; then
            BASE_TIME=$PROGRAM_TIME
            SPEEDUP=1.00
        else
            SPEEDUP=$(echo "scale=2; $BASE_TIME / $PROGRAM_TIME" | bc -l)
        fi
        
        echo "  - Tempo de execução: ${PROGRAM_TIME} segundos"
        echo "  - Speedup: ${SPEEDUP}x"
        echo ""
        
        # Salvar resultados no arquivo
        echo "$NODES,$PROCS_PER_NODE,$TOTAL_PROCS,$PROGRAM_TIME,$SPEEDUP" >> $RESULTS_FILE
        
        # Limpar arquivo temporário
        rm -f temp_output_${NODES}_${PROCS_PER_NODE}.txt
    done
done

echo "================================================"
echo "Testes concluídos! Resultados em: $RESULTS_FILE"
echo ""
echo "Resumo dos resultados:"
echo "================================================"
column -t -s ',' $RESULTS_FILE

# Encontrar a melhor configuração
echo ""
echo "Melhor configuração:"
tail -n +2 $RESULTS_FILE | sort -t',' -k5 -nr | head -1 | while IFS=',' read nodes procs_per_node total_procs time speedup; do
    echo "Nós: $nodes, Processos/nó: $procs_per_node, Total: $total_procs, Speedup: ${speedup}x"
done

echo ""
echo "Script de teste finalizado."