# Relatório atividade 15
Aluno: Davi Diógenes Ferreira de Almeida
Disciplina: IMD1116 Computação de Alto Desempenho/2025.1

Para esta atividade, implementamos um programa em MPI que simula a propagação de calor ao longo de uma barra unidimensional. A barra deve ser dividida em 3 partes, uma para cada processo, e os processos devem se comunicar por meio de mensageria MPI para informar aos processos de partes adjacentes das barras o valor dos pontos nas extremidades das repartições que sejam necessários para mais de um processo.

Foram implementados 3 programas com diferentes métodos de mensagem MPI. O primeiro programa usou MPI_Send (bloqueante) e MPI_Recv. O segundo programa usa MPI_ISend e MPI_IRecv (não bloqueantes) com MPI_Wait para bloquear até que a informação seja recebida. A terceira implementação também usa MPI_ISend e MPI_IRecv para enviar e indicar intenção de receber, porém usa MPI_Test para receber as informações. 

O código da versão usando send e receive bloqueantes foi:
```c
if (rank % 2 == 0) { // Processos pares: Envia e depois recebe
    if (rank < mpi_size - 1) { // Se tem alguem a direita, envia
        MPI_Send(&u_curr[local_n], 1, MPI_DOUBLE, rank + 1, TAG_DATA_GOES_RIGHT, MPI_COMM_WORLD);
    }
    if (rank > 0) { // Se tem alguem a esquerda, recebe
        MPI_Recv(&u_curr[0], 1, MPI_DOUBLE, rank - 1, TAG_DATA_GOES_RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
} else { // Impares: Recebe e depois envia
    if (rank > 0) { // Se tem alguem a esquerda, recebe
        MPI_Recv(&u_curr[0], 1, MPI_DOUBLE, rank - 1, TAG_DATA_GOES_RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank < mpi_size - 1) { // Se tem alguem a direita, envia
        MPI_Send(&u_curr[local_n], 1, MPI_DOUBLE, rank + 1, TAG_DATA_GOES_RIGHT, MPI_COMM_WORLD);
    }
}
// Enviar dados para a esquerda e receber da direita
if (rank % 2 == 0) { // Pares: Envia e depois recebe
    if (rank > 0) { // Se tem alguem a esquerda, envia
        MPI_Send(&u_curr[1], 1, MPI_DOUBLE, rank - 1, TAG_DATA_GOES_LEFT, MPI_COMM_WORLD);
    }
    if (rank < mpi_size - 1) { // Se tem alguem a direita, recebe
        MPI_Recv(&u_curr[local_n + 1], 1, MPI_DOUBLE, rank + 1, TAG_DATA_GOES_LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
} else { // Impares: Recebe e depois envia
    if (rank < mpi_size - 1) { // Se tem alguem a direita, recebe
        MPI_Recv(&u_curr[local_n + 1], 1, MPI_DOUBLE, rank + 1, TAG_DATA_GOES_LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank > 0) { // Se tem alguem a esquerda, envia
        MPI_Send(&u_curr[1], 1, MPI_DOUBLE, rank - 1, TAG_DATA_GOES_LEFT, MPI_COMM_WORLD);
    }
}
// Computa novos valores de u_next
```

A segunda implementação, com isend e irecv, está descrita abaixo:

```c
// 1. Iniciar todos os MPI_Irecv
// Receber do vizinho da esquerda (se existir) -> preenche u_curr[0]
if (rank > 0) {
    MPI_Irecv(&u_curr[0], 1, MPI_DOUBLE, rank - 1, TAG_DATA_GOES_RIGHT,
                MPI_COMM_WORLD, &requests[num_reqs++]);
}
// Receber do vizinho da direita (se existir) -> preenche u_curr[local_n + 1]
if (rank < mpi_size - 1) {
    MPI_Irecv(&u_curr[local_n + 1], 1, MPI_DOUBLE, rank + 1, TAG_DATA_GOES_LEFT,
                MPI_COMM_WORLD, &requests[num_reqs++]);
}
// 2. Iniciar todos os MPI_Isend
// Enviar para o vizinho da direita (se existir) o ponto u_curr[local_n]
if (rank < mpi_size - 1) {
    MPI_Isend(&u_curr[local_n], 1, MPI_DOUBLE, rank + 1, TAG_DATA_GOES_RIGHT,
                MPI_COMM_WORLD, &requests[num_reqs++]);
}
// Enviar para o vizinho da esquerda (se existir) o ponto u_curr[1]
if (rank > 0) {
    MPI_Isend(&u_curr[1], 1, MPI_DOUBLE, rank - 1, TAG_DATA_GOES_LEFT,
                MPI_COMM_WORLD, &requests[num_reqs++]);
}
// 3. Esperar que todas as comunicações não bloqueantes completem usando MPI_Wait individualmente
if (num_reqs > 0) {
    for (int i = 0; i < num_reqs; i++) {
        // Para cada requisição no array 'requests', esperamos sua conclusão.
        MPI_Wait(&requests[i], MPI_STATUS_IGNORE); 
    }
}
// Computa novos valores de u_next
```

Agora com o Test no lugar do Wait

```c
// 1. Iniciar todas as comunicações não bloqueantes (envios e recebimentos)
// Receber do vizinho da esquerda (se existir) -> preenche u_curr[0]
if (rank > 0) {
    MPI_Irecv(&u_curr[0], 1, MPI_DOUBLE, rank - 1, TAG_DATA_GOES_RIGHT,
                MPI_COMM_WORLD, &requests[num_reqs++]);
}
// Receber do vizinho da direita (se existir) -> preenche u_curr[local_n + 1]
if (rank < mpi_size - 1) {
    MPI_Irecv(&u_curr[local_n + 1], 1, MPI_DOUBLE, rank + 1, TAG_DATA_GOES_LEFT,
                MPI_COMM_WORLD, &requests[num_reqs++]);
}
// Enviar para o vizinho da direita (se existir) o ponto u_curr[local_n]
if (rank < mpi_size - 1) {
    MPI_Isend(&u_curr[local_n], 1, MPI_DOUBLE, rank + 1, TAG_DATA_GOES_RIGHT,
                MPI_COMM_WORLD, &requests[num_reqs++]);
}
// Enviar para o vizinho da esquerda (se existir) o ponto u_curr[1]
if (rank > 0) {
    MPI_Isend(&u_curr[1], 1, MPI_DOUBLE, rank - 1, TAG_DATA_GOES_LEFT,
                MPI_COMM_WORLD, &requests[num_reqs++]);
}
// 2. Calcular pontos internos (que não dependem das células fantasmas desta iteração)
for (int i = 2; i <= local_n - 1; ++i) {
    u_next[i] = u_curr[i] + factor * (u_curr[i - 1] - 2.0 * u_curr[i] + u_curr[i + 1]);
}
// 3. Esperar/Testar a conclusão das comunicações das células fantasmas
if (num_reqs > 0) {
    all_communications_done = 0;
    while (!all_communications_done) {
        MPI_Testall(num_reqs, requests, &all_communications_done, MPI_STATUSES_IGNORE);
    }
}
```

Tempos de execução:
- Programa 1 (bloqueante): 11.833754s
- Programa 2 (mpi_wait): 10.794234s
- Programa 3 (mpi_test): 8.710295s

Os resultados mostram que o programa com Send e Receive Bloqueantes foi o com o pior resultado.
Em seguida, vem o o segundo programa, com ISend, IRecv e Wait, que foi ligeiramente melhor que o primeiro, uma vez que todos os MPI_Isend são iniciados ao mesmo tempo, melhorando o desempenho do programa. Isso ocorre devido à sobreposição de comunicação.
O terceiro programa foi o mais rápido de todos, uma vez que, além de inicializar todas as operações de ISend e IRecv ao mesmo tempo, ele faz o processamento dos pontos internos da barra antes de usar MPI_Testall para receber os pontos das extremidades, fazendo uma sobreposição de trabalho.