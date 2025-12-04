#include <iostream>
#include <random>
#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <limits>

// ----------------------------
// Funzione: Dijkstra
// ----------------------------
std::vector<int> dijkstra(const std::vector<std::vector<int>>& M, size_t source) {
    size_t N = M.size();
    const int INF = std::numeric_limits<int>::max();
    std::vector<int> dist(N, INF);
    std::vector<bool> visited(N, false);

    dist[source] = 0;

    for (size_t k = 0; k < N; ++k) {
        int u = -1;
        int minDist = INF;

        for (size_t i = 0; i < N; ++i) {
            if (!visited[i] && dist[i] < minDist) {
                minDist = dist[i];
                u = static_cast<int>(i);
            }
        }

        if (u == -1) break; // Nessun nodo raggiungibile rimasto
        visited[u] = true;

        for (size_t v = 0; v < N; ++v) {
            if (M[u][v] > 0 && dist[u] + M[u][v] < dist[v]) {
                dist[v] = dist[u] + M[u][v];
            }
        }
    }

    return dist;
}

// ----------------------------
// Funzione: genera grafo casuale con P lati
// ----------------------------
std::vector<std::vector<int>> generateRandomGraph(size_t N, int P) {
    std::vector<std::vector<int>> M(N, std::vector<int>(N, 0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, N - 1);

    int edgesAdded = 0;
    while (edgesAdded < P) {
        size_t i = dist(gen);
        size_t j = dist(gen);
        if (i != j && M[i][j] == 0) {
            M[i][j] = M[j][i] = 1;
            edgesAdded++;
        }
    }

    return M;
}

// ----------------------------
// Funzione: calcola distanza media tra tutti i nodi
// ----------------------------
double averageDistance(size_t N, int P) {
    std::vector<std::vector<int>> M;
    std::vector<std::vector<int>> allDistances;
    const int MAX_TRIES = 5040;
    int tries = 0;
    bool validGraph = false;

    while (!validGraph) {
        M = generateRandomGraph(N, P);
        allDistances.clear();
        validGraph = true;

        for (size_t source = 0; source < N; ++source) {
            std::vector<int> dist = dijkstra(M, source);
            for (size_t i = 0; i < N; ++i) {
                if (dist[i] == std::numeric_limits<int>::max()) {
                    validGraph = false;
                    tries++;
                    if (tries > MAX_TRIES) throw std::runtime_error("Troppi tentativi connessi non validi");
                    break;
                }
            }
            if (!validGraph) break;
            allDistances.push_back(dist);
        }
    }

    // Calcolo distanza media (solo una volta per ogni coppia)
    int totalDist = 0;
    int count = 0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            totalDist += allDistances[i][j];
            count++;
        }
    }

    return static_cast<double>(totalDist) / count;
}

// Funzione che ogni thread esegue
void worker(std::queue<int>& tasks, int Ntests, std::mutex& queue_mutex) {
    while (true) {
        int N;

        // Prendere un task in modo thread-safe
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            if (tasks.empty()) return; // Nessun task rimasto
            N = tasks.front();
            tasks.pop();
        }

        int P = static_cast<int>(2.0 * N + 0.5);
        double avgDistance = 0.0;

        for (int t = 0; t < Ntests; ++t) {
            avgDistance += averageDistance(N, P);
        }

        // TODO: Stampa thread-safe
        {
            // std::lock_guard<std::mutex> lock(cout_mutex);
            std::cout << "(" << N << ", " << avgDistance / Ntests << ")," << std::endl;
        }
    }
}

// ----------------------------
// Main
// ----------------------------
int main() {
    const int Ntests = 1080;
    const int N_min = 12;
    const int N_max = 480;
    const int num_threads = std::thread::hardware_concurrency();

    // Mutex per accesso sicuro alla coda
    std::mutex queue_mutex;

    std::queue<int> tasks;
    for (int N = N_min; N < N_max; ++N) {
        tasks.push(N);
    }

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(worker, std::ref(tasks), Ntests,
                             std::ref(queue_mutex));
    }

    for (auto& th : threads) {
        th.join();
    }

    return 0;
}