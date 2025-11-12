#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <random>
#include <algorithm>
#include <filesystem>
#include <array>

constexpr int N_CHEMICALS = 100;  // Number of different chemicals (excluding empty=0)
constexpr int TOTAL_MOLECULES = 100000;  // Total molecules
constexpr int N_STEPS = 100000000;  // Number of timesteps
constexpr int MAX_RULES = 1000;  // Maximum number of rules
constexpr int SUPPLY_RATE = 10;  // Molecules supplied per timestep

constexpr int RECORDING_INTERVAL = N_STEPS / 1000; // Interval for recording state

// Rule structure: A+B+C -> D+E+F (values 0-N_CHEMICALS, 0=empty)
struct Rule {
    int reactants[3];  // A, B, C
    int products[3];   // D, E, F
    bool active;       // Whether this rule is in use
};

int get_mass(int chemical) {
    return chemical;  // Mass equals the chemical number (1-N_CHEMICALS)
}

void update(std::vector<int> &counts, std::vector<Rule> &rules, int &n_rules, 
            const std::vector<float> &volatility, std::mt19937 &gen, int timestep,
            std::ofstream &counts_file, std::ofstream &rules_file) {
    // Supply: convert empty to chemical 1
    int to_convert = std::min(SUPPLY_RATE, counts[0]);
    counts[0] -= to_convert;
    counts[1] += to_convert;
    
    // Sample 3 reactants based on abundance and volatility
    double total_weighted = 0.0;
    for (int c = 0; c <= N_CHEMICALS; ++c) {
        total_weighted += counts[c] * volatility[c];
    }
    
    if (total_weighted == 0) return;
    
    std::uniform_real_distribution<> sample_dis(0.0, total_weighted);
    int reactants[3];
    for (int i = 0; i < 3; ++i) {
        double rand_val = sample_dis(gen);
        double cumsum = 0.0;
        for (int c = 0; c <= N_CHEMICALS; ++c) {
            cumsum += counts[c] * volatility[c];
            if (cumsum > rand_val) {
                reactants[i] = c;
                break;
            }
        }
    }
    
    // Sort reactants
    std::sort(reactants, reactants + 3);
    
    // Check if we have enough reactants
    std::vector<int> req(N_CHEMICALS + 1, 0);
    for (int i = 0; i < 3; ++i) {
        req[reactants[i]]++;
    }
    bool can_react = true;
    for (int c = 0; c <= N_CHEMICALS; ++c) {
        if (req[c] > counts[c]) {
            can_react = false;
            break;
        }
    }
    
    if (!can_react) return;
    
    // Find matching rule
    int matching_rule = -1;
    for (int r = 0; r < n_rules; ++r) {
        if (rules[r].active &&
            rules[r].reactants[0] == reactants[0] &&
            rules[r].reactants[1] == reactants[1] &&
            rules[r].reactants[2] == reactants[2]) {
            matching_rule = r;
            break;
        }
    }
    
    // If no matching rule, create one
    if (matching_rule == -1 && n_rules < MAX_RULES) {
        int reactant_mass = get_mass(reactants[0]) + get_mass(reactants[1]) + get_mass(reactants[2]);
        
        // Enumerate all valid product combinations (mass <= reactant_mass)
        std::vector<std::array<int, 3>> valid_products;
        for (int p1 = 0; p1 <= N_CHEMICALS && p1 <= reactant_mass; ++p1) {
            for (int p2 = p1; p2 <= N_CHEMICALS && p1 + p2 <= reactant_mass; ++p2) {
                for (int p3 = p2; p3 <= N_CHEMICALS; ++p3) {
                    if (p1 + p2 + p3 <= reactant_mass) {
                        valid_products.push_back(std::array<int, 3>{p1, p2, p3});
                    }
                }
            }
        }
        
        if (!valid_products.empty()) {
            // Randomly select one valid product combination
            std::uniform_int_distribution<> prod_dis(0, valid_products.size() - 1);
            std::array<int, 3> selected = valid_products[prod_dis(gen)];
            
            rules[n_rules].reactants[0] = reactants[0];
            rules[n_rules].reactants[1] = reactants[1];
            rules[n_rules].reactants[2] = reactants[2];
            rules[n_rules].products[0] = selected[0];
            rules[n_rules].products[1] = selected[1];
            rules[n_rules].products[2] = selected[2];
            rules[n_rules].active = true;
            matching_rule = n_rules;
            n_rules++;
            
            // Write new rule to file
            rules_file << timestep << "\t" 
                      << reactants[0] << "\t" << reactants[1] << "\t" << reactants[2] << "\t"
                      << selected[0] << "\t" << selected[1] << "\t" << selected[2] << "\n";
            rules_file.flush();
        }
    }
    
    // Apply the rule if found or created
    if (matching_rule != -1) {
        for (int i = 0; i < 3; ++i) {
            counts[rules[matching_rule].reactants[i]]--;
        }
        for (int i = 0; i < 3; ++i) {
            counts[rules[matching_rule].products[i]]++;
        }
    }
    
    // Write counts periodically
    if (timestep % RECORDING_INTERVAL == 0) {
        counts_file << timestep;
        for (int c = 0; c <= N_CHEMICALS; ++c) {
            counts_file << "\t" << counts[c];
        }
        counts_file << "\n";
        counts_file.flush();
    }
}

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Global molecule counts
    std::vector<int> counts(N_CHEMICALS + 1);
    counts[0] = TOTAL_MOLECULES;  // All empty initially
    
    // Volatilities
    std::vector<float> volatility(N_CHEMICALS + 1);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int c = 0; c <= N_CHEMICALS; ++c) {
        // volatility[c] = dis(gen);
        volatility[c] = 1.0f;
    }
    
    // Initialize rules
    std::vector<Rule> rules(MAX_RULES);
    rules[0] = {{0, 0, 1}, {0, 0, 1}, true};  // 1+0+0 -> 1+0+0 (inert)
    rules[1] = {{0, 1, 1}, {0, 0, 2}, true};  // 1+1+0 -> 2+0+0
    rules[2] = {{0, 1, 2}, {0, 0, 3}, true};  // 2+1+0 -> 3+0+0
    rules[3] = {{0, 0, 2}, {0, 0, 2}, true};  // 2+0+0 -> 2+0+0
    int n_rules = 4;
    
    // Output directory
    std::filesystem::path output_dir = std::filesystem::current_path() / "outputs" / "timeseries" / 
                                     ("N_" + std::to_string(TOTAL_MOLECULES) + "_S_" + std::to_string(SUPPLY_RATE));
    std::filesystem::create_directories(output_dir);
    
    // Open output files
    std::ofstream counts_file(output_dir / "counts.tsv");
    std::ofstream rules_file(output_dir / "rules.tsv");
    
    // Write headers
    counts_file << "step";
    for (int c = 0; c <= N_CHEMICALS; ++c) {
        counts_file << "\tN_" << c;
    }
    counts_file << "\n";
    
    // Write volatility as first row
    counts_file << "volatility";
    for (int c = 0; c <= N_CHEMICALS; ++c) {
        counts_file << "\t" << volatility[c];
    }
    counts_file << "\n";
    
    rules_file << "step_added\tR_1\tR_2\tR_3\tP_1\tP_2\tP_3\n";
    
    // Write initial rules
    for (int r = 0; r < n_rules; ++r) {
        rules_file << 0 << "\t" 
                  << rules[r].reactants[0] << "\t" << rules[r].reactants[1] << "\t" << rules[r].reactants[2] << "\t"
                  << rules[r].products[0] << "\t" << rules[r].products[1] << "\t" << rules[r].products[2] << "\n";
    }
    
    // Write initial counts
    counts_file << 0;
    for (int c = 0; c <= N_CHEMICALS; ++c) {
        counts_file << "\t" << counts[c];
    }
    counts_file << "\n";
    
    for (int t = 1; t < N_STEPS; ++t) {
        update(counts, rules, n_rules, volatility, gen, t, counts_file, rules_file);
        if (t % 1000 == 0) {
            std::cout << "Progress:" << (t * 100) / N_STEPS << "%\r";
        }
    }
    
    counts_file.close();
    rules_file.close();
    
    return 0;
}