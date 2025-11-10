#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <filesystem>

constexpr int L = 25;  // Size of the lattice
constexpr int N_CHEMICALS = 10;  // Number of different chemicals (excluding empty=0)
constexpr int MOLECULES_PER_SITE = 100;  // Total molecules per site
constexpr int N_STEPS = 10000;  // Number of timesteps
constexpr int MAX_RULES = 1000;  // Maximum number of rules
constexpr int SUPPLY_RATE = 10;  // Molecules supplied per timestep at top
constexpr int D = 25;  // Diffusion coefficient: number of molecules to move per site per step

constexpr int RECORDING_INTERVAL = 1000; // Interval for recording lattice state

// Rule structure: A+B+C -> D+E+F (values 0-N_CHEMICALS, 0=empty)
struct Rule {
    int reactants[3];  // A, B, C
    int products[3];   // D, E, F
    bool active;       // Whether this rule is in use
};

__device__ int get_mass(int chemical) {
    return chemical;  // Mass equals the chemical number (1-N_CHEMICALS)
}

__device__ bool is_rule_legal(int r1, int r2, int r3, int p1, int p2, int p3) {
    int reactant_mass = get_mass(r1) + get_mass(r2) + get_mass(r3);
    int product_mass = get_mass(p1) + get_mass(p2) + get_mass(p3);
    return product_mass <= reactant_mass;
}

__global__ void initialize_lattice(int *lattice, float *volatility, curandState *states) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    int id = idy * L + idx;

    if (idx < L && idy < L) {
        curandState localState = states[id];
        
        // Initialize all sites with empty space (0)
        for (int c = 0; c <= N_CHEMICALS; ++c) {
            lattice[id * (N_CHEMICALS + 1) + c] = (c == 0) ? MOLECULES_PER_SITE : 0;
        }
        
        // Initialize volatilities (only once per block, shared across lattice)
        if (id == 0) {
            for (int c = 0; c <= N_CHEMICALS; ++c) {
                volatility[c] = curand_uniform(&localState);
            }
        }
        
        states[id] = localState;
    }
}

__global__ void supply_top(int *lattice) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < L) {
        int id = idx;  // Top row (idy = 0)
        int empty_idx = id * (N_CHEMICALS + 1) + 0;
        int chem1_idx = id * (N_CHEMICALS + 1) + 1;
        
        // Convert up to SUPPLY_RATE empty molecules to chemical 1
        int to_convert = min(SUPPLY_RATE, lattice[empty_idx]);
        lattice[empty_idx] -= to_convert;
        lattice[chem1_idx] += to_convert;
    }
}

__global__ void apply_reactions(int *lattice, Rule *rules, int n_rules, 
                                float *volatility, curandState *states) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    int id = idy * L + idx;

    if (idx < L && idy < L) {
        curandState localState = states[id];
        
        // Get counts at this site
        int counts[N_CHEMICALS + 1];
        for (int c = 0; c <= N_CHEMICALS; ++c) {
            counts[c] = lattice[id * (N_CHEMICALS + 1) + c];
        }
        
        // Find eligible rules and compute weights
        long long rule_weights[MAX_RULES];
        int eligible_rules[MAX_RULES];
        int num_eligible = 0;
        long long total_weight = 0;
        
        for (int r = 0; r < n_rules; ++r) {
            if (!rules[r].active) continue;
            
            // Check if we have enough reactants, accounting for duplicates
            int req[N_CHEMICALS + 1] = {0};
            for (int i = 0; i < 3; ++i) {
                req[rules[r].reactants[i]]++;
            }
            bool can_react = true;
            for (int c = 0; c <= N_CHEMICALS; ++c) {
                if (req[c] > counts[c]) {
                    can_react = false;
                    break;
                }
            }
            
            if (can_react) {
                // Compute weight: product of (count[c] * volatility[c]) for each reactant, considering multiplicity
                long long weight = 1;
                for (int i = 0; i < 3; ++i) {
                    int c = rules[r].reactants[i];
                    weight *= (long long)(counts[c] * volatility[c] * 1000);
                }
                rule_weights[num_eligible] = weight;
                eligible_rules[num_eligible] = r;
                total_weight += weight;
                num_eligible++;
            }
        }
        
        if (num_eligible > 0 && total_weight > 0) {
            // Select a rule proportionally to its weight
            long long rand_val = (long long)(curand_uniform(&localState) * total_weight);
            long long cumsum = 0;
            int selected_rule = -1;
            for (int i = 0; i < num_eligible; ++i) {
                cumsum += rule_weights[i];
                if (cumsum > rand_val) {
                    selected_rule = eligible_rules[i];
                    break;
                }
            }
            
            if (selected_rule != -1) {
                // Apply the selected rule
                for (int i = 0; i < 3; ++i) {
                    counts[rules[selected_rule].reactants[i]]--;
                }
                for (int i = 0; i < 3; ++i) {
                    counts[rules[selected_rule].products[i]]++;
                }
            }
        }
        
        // Write back to lattice
        for (int c = 0; c <= N_CHEMICALS; ++c) {
            lattice[id * (N_CHEMICALS + 1) + c] = counts[c];
        }
        
        states[id] = localState;
    }
}

__global__ void count_global_chemicals(int *lattice, int *global_counts) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    int id = idy * L + idx;

    if (idx < L && idy < L) {
        for (int c = 0; c <= N_CHEMICALS; ++c) {
            atomicAdd(&global_counts[c], lattice[id * (N_CHEMICALS + 1) + c]);
        }
    }
}

__global__ void add_new_rule(Rule *rules, int *n_rules, int *global_counts, 
                            float *volatility, curandState *states) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        curandState localState = states[0];
        
        // Calculate total weighted by volatility
        long long total_weighted = 0;
        for (int c = 0; c <= N_CHEMICALS; ++c) {
            total_weighted += (long long)(global_counts[c] * volatility[c] * 1000);
        }
        
        if (total_weighted == 0) return;
        
        // Pick 3 reactants based on volatility-weighted sampling
        int reactants[3];
        for (int i = 0; i < 3; ++i) {
            long long rand_val = (long long)(curand_uniform(&localState) * total_weighted);
            long long cumsum = 0;
            for (int c = 0; c <= N_CHEMICALS; ++c) {
                cumsum += (long long)(global_counts[c] * volatility[c] * 1000);
                if (cumsum > rand_val) {
                    reactants[i] = c;
                    break;
                }
            }
        }
        
        // Sort reactants for consistency
        for (int i = 0; i < 2; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                if (reactants[i] > reactants[j]) {
                    int temp = reactants[i];
                    reactants[i] = reactants[j];
                    reactants[j] = temp;
                }
            }
        }
        
        // Check if rule already exists (only reactants, since products not yet generated)
        bool exists = false;
        for (int r = 0; r < *n_rules; ++r) {
            if (rules[r].active &&
                rules[r].reactants[0] == reactants[0] &&
                rules[r].reactants[1] == reactants[1] &&
                rules[r].reactants[2] == reactants[2]) {
                exists = true;
                break;
            }
        }
        
        if (!exists && *n_rules < MAX_RULES) {
            // Generate legal products
            int reactant_mass = get_mass(reactants[0]) + get_mass(reactants[1]) + get_mass(reactants[2]);
            
            // Enumerate some legal product combinations
            int products[3];
            int attempts = 0;
            bool found_legal = false;
            
            while (!found_legal && attempts < 100) {
                // Random product generation
                for (int i = 0; i < 3; ++i) {
                    products[i] = curand(&localState) % (N_CHEMICALS + 1);
                }
                
                // Sort products
                for (int i = 0; i < 2; ++i) {
                    for (int j = i + 1; j < 3; ++j) {
                        if (products[i] > products[j]) {
                            int temp = products[i];
                            products[i] = products[j];
                            products[j] = temp;
                        }
                    }
                }
                
                if (is_rule_legal(reactants[0], reactants[1], reactants[2],
                                 products[0], products[1], products[2])) {
                    found_legal = true;
                }
                attempts++;
            }
            
            if (found_legal) {
                int rule_idx = *n_rules;
                rules[rule_idx].reactants[0] = reactants[0];
                rules[rule_idx].reactants[1] = reactants[1];
                rules[rule_idx].reactants[2] = reactants[2];
                rules[rule_idx].products[0] = products[0];
                rules[rule_idx].products[1] = products[1];
                rules[rule_idx].products[2] = products[2];
                rules[rule_idx].active = true;
                (*n_rules)++;
            }
        }
        
        states[0] = localState;
    }
}

__global__ void init_curand(curandState *states, unsigned long seed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    int id = idy * L + idx;

    if (idx < L && idy < L) {
        curand_init(seed, id, 0, &states[id]);
    }
}

void save_lattice(const std::vector<int> &lattice, const std::string &filename) {
    std::ofstream file(filename);
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            int site_idx = (i * L + j) * (N_CHEMICALS + 1);
            for (int c = 0; c <= N_CHEMICALS; ++c) {
                file << lattice[site_idx + c] << " ";
            }
            file << "\t";
        }
        file << "\n";
    }
    file.close();
}

void save_rules(const std::vector<Rule> &rules, int n_rules, const std::string &filename) {
    std::ofstream file(filename);
    file << "Active rules: " << n_rules << "\n";
    for (int i = 0; i < n_rules; ++i) {
        if (rules[i].active) {
            file << rules[i].reactants[0] << "+" << rules[i].reactants[1] << "+" 
                 << rules[i].reactants[2] << " -> "
                 << rules[i].products[0] << "+" << rules[i].products[1] << "+" 
                 << rules[i].products[2] << "\n";
        }
    }
    file.close();
}


__global__ void diffuse(int *lattice, curandState *states, int phase) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    int id = idy * L + idx;

    if (idx >= L || idy >= L) return;

    // Calculate phase according to the pattern shown in readme
    int computed_phase = (idx + idy * 2) % 6;
    if (computed_phase != phase) {
        return;
    }

    curandState localState = states[id];
    int stride = N_CHEMICALS + 1;
    
    // Define all possible neighbors (up, down, left, right)
    int neighbors[4][2];
    int num_neighbors = 0;
    
    if (idx > 0) { neighbors[num_neighbors][0] = idx - 1; neighbors[num_neighbors][1] = idy; num_neighbors++; }      // left
    if (idx + 1 < L) { neighbors[num_neighbors][0] = idx + 1; neighbors[num_neighbors][1] = idy; num_neighbors++; }  // right
    if (idy > 0) { neighbors[num_neighbors][0] = idx; neighbors[num_neighbors][1] = idy - 1; num_neighbors++; }      // up
    if (idy + 1 < L) { neighbors[num_neighbors][0] = idx; neighbors[num_neighbors][1] = idy + 1; num_neighbors++; }  // down
    
    // Perform D swaps with EACH neighbor
    for (int n = 0; n < num_neighbors; ++n) {
        int nidx = neighbors[n][0];
        int nidy = neighbors[n][1];
        int nid = nidy * L + nidx;
        
        // Do D swaps with this specific neighbor
        for (int i = 0; i < D; ++i) {
        // Do D swaps with this specific neighbor
        for (int i = 0; i < D; ++i) {
            // Both sites always have exactly MOLECULES_PER_SITE total molecules
            // Weighted random sampling from current site
            int rand_current = curand(&localState) % MOLECULES_PER_SITE;
            int c_a = 0;
            int cumsum = lattice[id * stride + 0];
            while (cumsum <= rand_current && c_a < stride - 1) {
                c_a++;
                cumsum += lattice[id * stride + c_a];
            }
            
            // Weighted random sampling from neighbor site
            int rand_neighbor = curand(&localState) % MOLECULES_PER_SITE;
            int c_b = 0;
            cumsum = lattice[nid * stride + 0];
            while (cumsum <= rand_neighbor && c_b < stride - 1) {
                c_b++;
                cumsum += lattice[nid * stride + c_b];
            }
            
            // Perform swap - guaranteed to succeed since we sampled from existing molecules
            lattice[id * stride + c_a] -= 1;
            lattice[nid * stride + c_a] += 1;
            lattice[nid * stride + c_b] -= 1;
            lattice[id * stride + c_b] += 1;
            }
        }
    
    states[id] = localState;
    }

}


int main() {
    int *d_lattice;
    float *d_volatility;
    Rule *d_rules;
    int *d_n_rules;
    int *d_global_counts;
    curandState *d_states;

    cudaMalloc(&d_lattice, L * L * (N_CHEMICALS + 1) * sizeof(int));
    cudaMalloc(&d_volatility, (N_CHEMICALS + 1) * sizeof(float));
    cudaMalloc(&d_rules, MAX_RULES * sizeof(Rule));
    cudaMalloc(&d_n_rules, sizeof(int));
    cudaMalloc(&d_global_counts, (N_CHEMICALS + 1) * sizeof(int));
    cudaMalloc(&d_states, L * L * sizeof(curandState));

    dim3 blockSize(16, 16);
    dim3 gridSize((L + blockSize.x - 1) / blockSize.x, (L + blockSize.y - 1) / blockSize.y);
    dim3 topBlockSize(256);
    dim3 topGridSize((L + topBlockSize.x - 1) / topBlockSize.x);

    init_curand<<<gridSize, blockSize>>>(d_states, time(0));
    cudaDeviceSynchronize();

    initialize_lattice<<<gridSize, blockSize>>>(d_lattice, d_volatility, d_states);
    cudaDeviceSynchronize();

    // Initialize with specified rules
    std::vector<Rule> h_rules(MAX_RULES);
    h_rules[0] = {{0, 0, 1}, {0, 0, 1}, true};  // 1+0+0 -> 1+0+0 (inert)
    h_rules[1] = {{0, 1, 1}, {0, 0, 2}, true};  // 1+1+0 -> 2+0+0
    h_rules[2] = {{0, 1, 2}, {0, 0, 3}, true};  // 2+1+0 -> 3+0+0
    h_rules[3] = {{0, 0, 2}, {0, 0, 2}, true};  // 2+0+0 -> 2+0+0
    int h_n_rules = 4;
        
    cudaMemcpy(d_rules, h_rules.data(), MAX_RULES * sizeof(Rule), cudaMemcpyHostToDevice);
    cudaMemcpy(d_n_rules, &h_n_rules, sizeof(int), cudaMemcpyHostToDevice);

    std::vector<int> h_lattice(L * L * (N_CHEMICALS + 1));

    // Dynamically create outputs/latticeTimeseries/L_<L>_D_<D> directory
    std::filesystem::path output_dir = std::filesystem::current_path() / "outputs" / "latticeTimeseries" / ("L_" + std::to_string(L) + "_D_" + std::to_string(D));
    std::filesystem::create_directories(output_dir);

    int prev_n_rules = 0;  // Track previous number of rules

    for (int t = 0; t < N_STEPS; ++t) {
        
        // Supply chemical 1 at top
        supply_top<<<topGridSize, topBlockSize>>>(d_lattice);
        cudaDeviceSynchronize();
        cudaError_t error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA Error after supply_top at step " << t << ": " << cudaGetErrorString(error) << std::endl;
        }

        // Apply reactions
        apply_reactions<<<gridSize, blockSize>>>(d_lattice, d_rules, h_n_rules, 
                                                 d_volatility, d_states);
        cudaDeviceSynchronize();
        error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA Error after apply_reactions at step " << t << ": " << cudaGetErrorString(error) << std::endl;
        }

        // Diffuse molecules in 6 phases
        for (int phase = 0; phase < 6; ++phase) {
            diffuse<<<gridSize, blockSize>>>(d_lattice, d_states, phase);
            cudaDeviceSynchronize();
            error = cudaGetLastError();
            if (error != cudaSuccess) {
                std::cerr << "CUDA Error after diffuse phase " << phase << " at step " << t << ": " << cudaGetErrorString(error) << std::endl;
            }
        }

        // Count global chemicals
        cudaMemset(d_global_counts, 0, (N_CHEMICALS + 1) * sizeof(int));
        count_global_chemicals<<<gridSize, blockSize>>>(d_lattice, d_global_counts);
        cudaDeviceSynchronize();
        error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA Error after count_global_chemicals at step " << t << ": " << cudaGetErrorString(error) << std::endl;
        }

        // Add new rule every timestep
        add_new_rule<<<1, 1>>>(d_rules, d_n_rules, d_global_counts, 
                              d_volatility, d_states);
        cudaDeviceSynchronize();
        error = cudaGetLastError();
        if (error != cudaSuccess) {
            std::cerr << "CUDA Error after add_new_rule at step " << t << ": " << cudaGetErrorString(error) << std::endl;
        }

        // Update n_rules
        cudaMemcpy(&h_n_rules, d_n_rules, sizeof(int), cudaMemcpyDeviceToHost);

        // Check if a new rule was added
        if (h_n_rules > prev_n_rules) {
            std::cout << "New rule added at step " << t << std::endl;
        }
        prev_n_rules = h_n_rules;

        // Save state periodically
        if (t % RECORDING_INTERVAL == 0) {
            cudaMemcpy(h_lattice.data(), d_lattice, 
                      L * L * (N_CHEMICALS + 1) * sizeof(int), 
                      cudaMemcpyDeviceToHost);
            save_lattice(h_lattice, (output_dir / ("lattice_" + std::to_string(t) + ".txt")).string());
            
            cudaMemcpy(h_rules.data(), d_rules, MAX_RULES * sizeof(Rule), 
                      cudaMemcpyDeviceToHost);
            save_rules(h_rules, h_n_rules, (output_dir / ("rules_" + std::to_string(t) + ".txt")).string());
        }

    }

    cudaFree(d_lattice);
    cudaFree(d_volatility);
    cudaFree(d_rules);
    cudaFree(d_n_rules);
    cudaFree(d_global_counts);
    cudaFree(d_states);

    return 0;
}