#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>

extern "C" {
    std::map<std::string, std::vector<int>> featureInfo;
    std::map<int, std::vector<int>> star_neighbors;
    std::vector<std::vector<std::string>> size2_patterns;
    std::map<std::vector<std::string>, std::map<std::vector<int>, std::vector<int>>> instance_table;
    std::map<std::vector<std::string>, std::map<std::string, std::set<int>>> hashmap;

    // read the featureInfo from csv file and store
    void read_featureInfo() {
        std::ifstream featureInfo_file("required_files/featureInfo.csv");

        if (!featureInfo_file.is_open()) {
            std::cerr << "Error opening file" << std::endl;
        }

        std::string line;
        std::getline(featureInfo_file, line);

        while (std::getline(featureInfo_file, line)) {
            std::istringstream lineStream(line);
            std::string cell;
            std::string key;
            std::vector<int> values;

            if (std::getline(lineStream, cell, ',')) {
                key = cell;
            }

            while (std::getline(lineStream, cell, ',')) {
                values.push_back(std::stoi(cell));
            }
            featureInfo[key] = values;
        }

        featureInfo_file.close();
    }

    // read the star_neighbors csv file and store
    void read_star_neighbors() {
        
        std::ifstream star_neighbors_file("required_files/starNeighbors.csv");

        if (!star_neighbors_file.is_open()) {
            std::cerr << "Error opening file" << std::endl;
        }

        std::string line1;
        std::getline(star_neighbors_file, line1);

        while (std::getline(star_neighbors_file, line1)) {
            std::stringstream ss(line1);
            int key;
            int value;
            std::vector<int> values;
            
            ss >> std::ws;
            ss >> key;
            ss.ignore(); 

            while (ss >> value) {
                values.push_back(value);
                ss.ignore(); 
            }
            star_neighbors[key] = values;
        }

        star_neighbors_file.close();
    }
    
    // generate the size 2 candidate patterns
    std::vector<std::vector<std::string>> generate_size2_combos() {
        std::vector<std::string> features;
        for (const auto &entry : featureInfo) {
            features.push_back(entry.first);
        }

        std::vector<std::vector<std::string>> size2_candidatePatterns;
        for (size_t i = 0; i < features.size(); ++i) {
            for (size_t j = i + 1; j < features.size(); ++j) {
                size2_candidatePatterns.push_back({features[i], features[j]});
            }
        }
        return size2_candidatePatterns;
    }

    // log(n) function that finds values within a given range
    std::vector<int> findNeighborsInRange(const std::vector<int>& arr, int x, int y) {
        auto start_it = std::lower_bound(arr.begin(), arr.end(), x);
        auto end_it = std::upper_bound(arr.begin(), arr.end(), y);
        return std::vector<int>(start_it, end_it);
    }
    
    // find the prevalent size 2 patterns
    std::vector<std::vector<std::string>> degree2Processing(
        std::vector<std::vector<std::string>> candidatePatterns, int candidatePatterns_size, 
        double prevalence_threshold) {   
        std::vector<std::vector<std::string>> size2_patterns;

        // iterate over each candidate pattern
        for (const auto& coloc : candidatePatterns) {
            std::string first_feature = coloc[0];
            std::string second_feature = coloc[1];

            // determine first feature start and end
            auto feature1 = featureInfo.find(first_feature);
            std::vector<int> values1 = feature1->second;
            int first_feature_start = values1[1];
            int first_feature_end = values1[2];
            // determine second feature start and end
            auto feature2 = featureInfo.find(second_feature);
            std::vector<int> values2 = feature2->second;
            int second_feature_start = values2[1];
            int second_feature_end = values2[2];

            // create the key for the instance table and the hashmap and initialize
            std::vector<std::string> coloc_key = {first_feature, second_feature};
            instance_table[coloc_key] = {};
            hashmap[coloc_key] = {};
            hashmap[coloc_key][first_feature] = {};
            hashmap[coloc_key][second_feature] = {};

            for (int index = first_feature_start; index <= first_feature_end; index++) {
                auto star_neighbor_it = star_neighbors.find(index);
                std::vector<int> neighbors = findNeighborsInRange(star_neighbor_it->second, 
                                                                  second_feature_start,
                                                                  second_feature_end);
                if (!neighbors.empty()) {
                    std::vector<int> index_tuple = {index};
                    instance_table[coloc_key][index_tuple] = neighbors;
                    hashmap[coloc_key][first_feature].insert(index);
                    for (int neighbor : neighbors) {
                        hashmap[coloc_key][second_feature].insert(neighbor);
                    }
                }
            }

            // calculate the participation ratios
            double pr_first_feature = static_cast<double>(hashmap[coloc_key][first_feature].size()) /
                featureInfo[first_feature][0];
            double pr_second_feature = static_cast<double>(hashmap[coloc_key][second_feature].size()) /
                featureInfo[second_feature][0];
            
            // calculate the participation index
            double PI = 0.0;
            if (pr_first_feature < pr_second_feature) {
                PI = pr_first_feature;
            } else {
                PI = pr_second_feature;
            }

            // store the pattern if it is prevalent
            if (PI >= prevalence_threshold) {
                   size2_patterns.push_back(coloc_key);
            }
        }

        std::cout << "Degree 2 Prevalent Patterns:" << std::endl;
        for (auto i : size2_patterns) {
            std::cout << "(" << i[0] << ", " << i[1] << ")" << std::endl;
        }
        
        // remove the non-prevalent patterns from the instance table
        for (auto it = instance_table.begin(); it != instance_table.end(); ) {
            if (std::find(size2_patterns.begin(), size2_patterns.end(), it->first) == size2_patterns.end()) {
                it = instance_table.erase(it);
            } else {
                ++it;
            }
        }
        
        return size2_patterns;
    }
    
    // helper function: generates the combinations for the candidate patterns with degree > 2
    void generateCombinations(const std::vector<std::string>& features, int degree, 
                                  std::vector<std::vector<std::string>>& result, 
                                  std::vector<std::string>& current, int start) {
        if (current.size() == degree) {
            result.push_back(current);
            return;
        }
        for (int i = start; i < features.size(); ++i) {
            current.push_back(features[i]);
            generateCombinations(features, degree, result, current, i + 1);
            current.pop_back();
        }
    }

    // helper function: check if all (degree-1)-subpatterns of a pattern are in the prevalent patterns
    bool allSubpatternsInPrevalent(const std::vector<std::string>& pattern, 
                                   const std::set<std::vector<std::string>>& prevalentPatterns, 
                                   int degree) {
        std::vector<std::vector<std::string>> subpatterns;
        std::vector<std::string> current;
        generateCombinations(pattern, degree - 1, subpatterns, current, 0);

        for (const auto& subpattern : subpatterns) {
            if (prevalentPatterns.find(subpattern) == prevalentPatterns.end()) {
                return false;
            }
        }
        return true;
    }

    // generate candidate patterns of degree > 2 (calls helper functions)
    std::vector<std::vector<std::string>> getCandidatePatterns(
        const std::vector<std::vector<std::string>>& prevalentPattern, int degree) {
        std::set<std::vector<std::string>> prevalentPatterns(prevalentPattern.begin(), prevalentPattern.end());
        // extract features from the keys of featureInfo
        std::vector<std::string> features;
        for (const auto& pair : featureInfo) {
            features.push_back(pair.first);
        }

        std::vector<std::vector<std::string>> _candidatePatterns;
        std::vector<std::vector<std::string>> _patterns;
        std::vector<std::string> current;

        generateCombinations(features, degree, _patterns, current, 0);

        for (const auto& pattern : _patterns) {
            if (allSubpatternsInPrevalent(pattern, prevalentPatterns, degree)) {
                _candidatePatterns.push_back(pattern);
            }
        }

        return _candidatePatterns;
    }
    
    // find prevalent patterns of size > 2
    std::vector<std::vector<std::string>> colocationGeneral(
            std::vector<std::vector<std::string>> candidatePatterns, int candidatePatterns_size, 
                           double prevalence_threshold, int degree) {
            
        std::vector<std::vector<std::string>> prevalent;

        for (const auto& currentPattern : candidatePatterns) {
            std::vector<std::string> basePattern;
            for (int j = 0; j < degree - 1; j++) {
                basePattern.push_back(currentPattern[j]);
            }
            std::string lastFeature = currentPattern[degree - 1];
            // add  and initialize entries in instance_table and hashmap
            instance_table[currentPattern] = {};
            hashmap[currentPattern] = {};

            // initialize each element in currentPattern in hashmap
            for (const auto& f : currentPattern) {
                hashmap[currentPattern][f] = {};
            }
            auto colocTableIt = instance_table.find(basePattern);
            std::map<std::vector<int>, std::vector<int>>& colocTable = colocTableIt->second;

            for (const auto& entry : colocTable) {
                const std::vector<int>& key = entry.first;
                std::set<int> commonLastNeighbors;
                for (int instanceID : key) {
                    auto star_neighbor_it = star_neighbors.find(instanceID);
                    if (commonLastNeighbors.empty()) {
                        std::vector<int> temp_vector = findNeighborsInRange(star_neighbor_it->second,
                                                                     featureInfo[lastFeature][1],
                                                                     featureInfo[lastFeature][2]);
                        commonLastNeighbors.insert(temp_vector.begin(), temp_vector.end());
                    } else {
                        std::vector<int> temp_vector = findNeighborsInRange(star_neighbor_it->second, 
                                                                              featureInfo[lastFeature][1],
                                                                              featureInfo[lastFeature][2]);
                        std::set<int> temp_commonLastNeighbors(temp_vector.begin(), temp_vector.end());
                        std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                              temp_commonLastNeighbors.begin(), temp_commonLastNeighbors.end(),
                                              std::inserter(commonLastNeighbors, commonLastNeighbors.begin()));
                    }
                }
                for (int n : colocTable[key]) {
                    auto star_neighbor_it = star_neighbors.find(n);
                    std::vector<int> temp_vect = findNeighborsInRange(star_neighbor_it->second,
                                                                      featureInfo[lastFeature][1],
                                                                      featureInfo[lastFeature][2]);
                    std::set<int> temp_neighbors(temp_vect.begin(), temp_vect.end());
                    std::set<int> neighbors;
                    std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                          temp_neighbors.begin(), temp_neighbors.end(),
                                          std::inserter(neighbors, neighbors.begin()));

                    if (!neighbors.empty()) {
                        std::vector<int> new_key = key;
                        new_key.push_back(n);
                        std::vector<int> intersectionVec(neighbors.begin(), neighbors.end());
                        instance_table[currentPattern][new_key] = intersectionVec;
        
                        for (size_t k = 0; k < new_key.size(); ++k) {
                            hashmap[currentPattern][currentPattern[k]].insert(new_key[k]);
                        }
                        hashmap[currentPattern][lastFeature].insert(neighbors.begin(), neighbors.end());
                    }
                }
            }
            
            // calculate participation ratios
            std::vector<double> pr;
            for (int m = 0; m < degree; ++m) {
                std::string f = currentPattern[m];
                double ratio = static_cast<double>(hashmap[currentPattern][f].size()) 
                    / featureInfo[f][0];
                pr.push_back(ratio);
            }
            
            // calculate participation index
            double PI = *std::min_element(pr.begin(), pr.end());
            if (PI >= prevalence_threshold) {
                prevalent.push_back(currentPattern);
            }
        }

        std::cout << "Degree " << degree << " Prevalent Patterns:"<< std::endl;
        for (auto& patternVec : prevalent) {
            std::cout << "(";
            for (size_t i = 0; i < patternVec.size(); i++) {
                std::cout << patternVec[i];
                if (i < patternVec.size() - 1) {
                    std::cout << ", ";
                }
            }
            std::cout << ")" << std::endl;
        }
        
        // delete non-prevalent patterns from instance table
        for (auto it = instance_table.begin(); it != instance_table.end(); ) {
            if (std::find(prevalent.begin(), prevalent.end(), it->first) == prevalent.end()) {
                it = instance_table.erase(it);
            } else {
                ++it;
            }
        }
        
        return prevalent;
    }
    
    // main function for processing
    void colocation_main(double prevalence_threshold) {
        read_featureInfo();
        read_star_neighbors();
        std::vector<std::vector<std::string>> size2_candidatePatterns = 
            generate_size2_combos();
        std::vector<std::vector<std::string>> size2_patterns = 
            degree2Processing(size2_candidatePatterns,
                                            size2_candidatePatterns.size(), 
                                            prevalence_threshold);
        
        int degree = 3;
        std::vector<std::vector<std::string>> candidatePatterns =
            getCandidatePatterns(size2_patterns, degree);
        while (!candidatePatterns.empty()) {
            std::vector<std::vector<std::string>> prevalent_patterns = 
            colocationGeneral(candidatePatterns, candidatePatterns.size(), 
                       prevalence_threshold, degree);
            
            degree += 1;
            if (prevalent_patterns.size() == 0) {
                break;
            }
            candidatePatterns = getCandidatePatterns(prevalent_patterns, degree);
            
            if (degree > 3) {
                // Find keys to delete
                std::vector<std::vector<std::string>> to_delete;
                for (const auto& outer_pair : instance_table) {
                    const std::vector<std::string>& key = outer_pair.first;
                    if (key.size() == degree - 2) {
                        to_delete.push_back(key);
                    }
                }

                // Delete keys from instance_table
                for (const auto& key : to_delete) {
                    instance_table.erase(key);
                }
            }
        }
    }
}
