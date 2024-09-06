import networkx as nx
import time
from genetic_algorithm_l321 import Genetic_Algorithm_L321

# Função para criar o grafo a partir do arquivo.
def create_graph_from_file(filename, graph_is_0_indexed):
    G = nx.Graph()
    with open(filename, 'r') as file:
        for line in file:
            if line[0] == 'c' or line[0] == 'p':
                continue
            if line[0] == 'e':
                line = line.rstrip()
                e = line.split(' ')
                u = int(e[1])
                v = int(e[2])
                if graph_is_0_indexed == False:
                    u = u - 1
                    v = v - 1
                G.add_edge(u,v)
    return G

# Função para executar o algoritmo genético e salvar os resultados em um arquivo.
def run_experiment(graph_name, graph, population_size, generations, crossover_rate, mutation_rate, elitism_rate, file_name):

    ga_l321 = Genetic_Algorithm_L321(graph=graph,
                                     population_size=population_size,
                                     generations=generations,
                                     crossover_rate=crossover_rate,
                                     mutation_rate=mutation_rate,
                                     elitism_rate=elitism_rate)

    start_time = time.time()
    best_fitness = ga_l321.run()
    end_time = time.time()
    execution_time = end_time - start_time

    with open(file_name, 'a') as file:
        file.write(f"Grafo: {graph_name}\n")
        file.write(f"Melhor fitness: {best_fitness}\n")
        file.write(f"Tempo de execução: {execution_time:.4f} segundos\n")
        file.write("\n")

    print(f"Grafo: {graph_name}")
    print(f"Melhor fitness encontrado: {best_fitness}")
    print(f"Tempo de execução: {execution_time:.4f} segundos\n")

# Função main.
# Instancia os grafos e chama a função run_experiment.
# Salva os resultados em um arquivo.
# Executa a função main.
def main():

    file_name = "result.txt"
    open(file_name, 'w').close()

    population_size = 20
    generations = 100
    crossover_rate = 0.8
    mutation_rate = 0.2
    elitism_rate = 0.1

 # Grafos gerados pela biblioteca NetworkX
    graphs = [
        ("complete_graph(20)", nx.complete_graph(20)),
        ("complete_graph(70)", nx.complete_graph(70)),
        ("complete_bipartite_graph(20, 10)", nx.complete_bipartite_graph(20, 10)),
        ("complete_bipartite_graph(50, 50)", nx.complete_bipartite_graph(50, 50)),
        ("cycle_graph(10)", nx.cycle_graph(10)),
        ("cycle_graph(30)", nx.cycle_graph(30)),
        ("cycle_graph(50)", nx.cycle_graph(50)),
        ("path_graph(50)", nx.path_graph(50)),
        ("path_graph(70)", nx.path_graph(70)),
        ("path_graph(100)", nx.path_graph(100)),
        ("grid_2d_graph(10, 8)", nx.grid_2d_graph(10, 8)),
        ("grid_2d_graph(7, 6)", nx.grid_2d_graph(7, 6)),
        ("grid_2d_graph(3, 35)", nx.grid_2d_graph(3, 35)),
    ]

    # Executa o algoritmo genético para os grafos da biblioteca NetworkX
    for graph_name, graph in graphs:
        run_experiment(graph_name, graph, population_size, generations, crossover_rate, mutation_rate, elitism_rate, file_name)

    # # Grafos gerados a partir dos arquivos.
    # file_graphs = [
    #     ("Grafo do arquivo 1", create_graph_from_file("dsjc250.5.col", False)),
    #     ("Grafo do arquivo 2", create_graph_from_file("dsjc500.1.col", False)),
    #     ("Grafo do arquivo 3", create_graph_from_file("dsjc500.5.col", False)),
    #     ("Grafo do arquivo 4", create_graph_from_file("dsjc1000.1.col", False)),
    # ]

    # # # Executa o algoritmo genético para os grafos dos arquivos
    # for graph_name, graph in file_graphs:
    #     run_experiment(graph_name, graph, population_size, generations, crossover_rate, mutation_rate, elitism_rate, file_name)


# Executa a função main.
if __name__ == "__main__":
    main()