#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>

extern "C" {
    // Holds all information pertaining to the sub-regions
    class SubRegion {
    public:
        int subregion_id;  // unique identifier
        std::map<std::string, std::vector<int>> featureInfo;
        std::map<int, std::vector<int>> star_neighbors;
        std::vector<std::vector<std::string>> size2_patterns;
        std::map<std::vector<std::string>, std::map<std::vector<int>, std::vector<int>>> instance_table;
        std::map<std::vector<std::string>, std::map<std::string, std::set<int>>> hashmap;
        
        SubRegion(int number) {
            subregion_id = number;  // set the unique identifier
        }

        // get the featureInfo
        void read_featureInfo() {
            std::ifstream featureInfo_file("required_files/featureInfo/featureInfo" + std::to_string(this->subregion_id) + ".csv");
            if (!featureInfo_file.is_open()) {
                std::cerr << "ERROR OPENING SUB-REGION " << this->subregion_id << " FEATURE INFO FILE!" << std::endl;
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
                this->featureInfo[key] = values;
            }
            featureInfo_file.close();
        }

        // get the star_neighbors from csv file and store
        void read_star_neighbors() {
            std::ifstream star_neighbors_file("required_files/starNeighbors/starNeighbors" + 
                                              std::to_string(this->subregion_id) + ".csv");
            if (!star_neighbors_file.is_open()) {
                std::cerr << "ERROR OPENING SUB-REGION " << this->subregion_id << " STAR NEIGHBORS FILE!" << std::endl;
            }
            std::string line;
            std::getline(star_neighbors_file, line);
            while (std::getline(star_neighbors_file, line)) {
                std::stringstream ss(line);
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
                this->star_neighbors[key] = values;
            }
            star_neighbors_file.close();
        }

        // generate the size 2 candidate patterns
        std::vector<std::vector<std::string>> generate_size2_combos() {
            std::vector<std::string> features;
            for (const auto &entry : this->featureInfo) {
                features.push_back(entry.first);
            }

            std::vector<std::vector<std::string>> size2_candidatePatterns;
            for (size_t i = 0; i < features.size(); i++) {
                for (size_t j = i + 1; j < features.size(); j++) {
                    size2_candidatePatterns.push_back({features[i], features[j]});
                }
            }
            return size2_candidatePatterns;
        }

        // log(n) search that finds neighbors within a specific start and end range
        std::vector<int> findNeighborsInRange(const std::vector<int>& arr, int x, int y) {
            auto start_it = std::lower_bound(arr.begin(), arr.end(), x);
            auto end_it = std::upper_bound(arr.begin(), arr.end(), y);
            return std::vector<int>(start_it, end_it);
        }

        // find the degree 2 prevalent patterns
        std::vector<std::vector<std::string>> degree2Processing(
            std::vector<std::vector<std::string>> candidatePatterns, int candidatePatterns_size, 
            double prevalence_threshold) {   
            std::vector<std::vector<std::string>> size2_patterns;
            
            // iterate over each candidate pattern
            for (const auto& coloc : candidatePatterns) {
                std::string first_feature = coloc[0];
                std::string second_feature = coloc[1];
                
                // determine first feature start and end
                auto feature1 = this->featureInfo.find(first_feature);
                std::vector<int> values1 = feature1->second;
                int first_feature_start = values1[1];
                int first_feature_end = values1[2];
                // determine second feature start and end
                auto feature2 = this->featureInfo.find(second_feature);
                std::vector<int> values2 = feature2->second;
                int second_feature_start = values2[1];
                int second_feature_end = values2[2];
                
                // create the key and store in instance table and hashmap
                std::vector<std::string> coloc_key = {first_feature, second_feature};
                this->instance_table[coloc_key] = {};
                this->hashmap[coloc_key] = {};
                this->hashmap[coloc_key][first_feature] = {};
                this->hashmap[coloc_key][second_feature] = {};
                
                for (int index = first_feature_start; index <= first_feature_end; index++) {
                    auto star_neighbor_it = this->star_neighbors.find(index);
                    std::vector<int> neighbors = findNeighborsInRange(star_neighbor_it->second, 
                                                                      second_feature_start,
                                                                      second_feature_end);
                    if (!neighbors.empty()) {
                        std::vector<int> index_tuple = {index};
                        this->instance_table[coloc_key][index_tuple] = neighbors;
                        this->hashmap[coloc_key][first_feature].insert(index);
                        for (int neighbor : neighbors) {
                            this->hashmap[coloc_key][second_feature].insert(neighbor);
                        }
                    }
                }

                // calculate the participation rations
                double pr_first_feature = static_cast<double>(this->hashmap[coloc_key][first_feature].size()) /
                    this->featureInfo[first_feature][0];
                double pr_second_feature = static_cast<double>(this->hashmap[coloc_key][second_feature].size()) /
                    this->featureInfo[second_feature][0];
                
                double PI = 0.0;
                // calculate the participation index
                if (pr_first_feature < pr_second_feature) {
                    PI = pr_first_feature;
                } else {
                    PI = pr_second_feature;
                }

                // determine if the pattern is prevalent or not
                if (PI >= prevalence_threshold) {
                       size2_patterns.push_back(coloc_key);
                }
            }
            
            std::cout << "Degree 2 Prevalent Patterns for Sub-Region " << this->subregion_id << ":" << std::endl;
            std::ofstream size2_file("size2_patterns_subregion" + std::to_string(this->subregion_id) + ".txt");
            for (auto i : size2_patterns) {
                std::cout << "(" << i[0] << ", " << i[1] << ")" << std::endl;
                size2_file << "(" << i[0] << ", " << i[1] << ")" << std::endl;
            }
            size2_file.close();
            return size2_patterns;
        }

        // generates the combinations for the candidate patterns with degree > 2
        void generateCombinations(const std::vector<std::string>& features, int degree, 
                                  std::vector<std::vector<std::string>>& result, 
                                  std::vector<std::string>& current, int start) {
            if (current.size() == degree) {
                result.push_back(current);
                return;
            }
            for (int i = start; i < features.size(); i++) {
                current.push_back(features[i]);
                generateCombinations(features, degree, result, current, i + 1);
                current.pop_back();
            }
        }

        /* helper function for getCandidatePatterns: check if all (degree-1)-subpatterns of a pattern 
        are in the prevalent patterns */
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

        // generate candidate patterns of degree > 2
        std::vector<std::vector<std::string>> getCandidatePatterns(
            const std::vector<std::vector<std::string>>& prevalentPattern, int degree) {
            std::set<std::vector<std::string>> prevalentPatterns(prevalentPattern.begin(), prevalentPattern.end());
            std::vector<std::string> features;
            for (const auto& pair : this->featureInfo) {
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

        // find prevalent patterns for degree > 2
        std::vector<std::vector<std::string>> colocationGeneral(
            std::vector<std::vector<std::string>> candidatePatterns, int candidatePatterns_size, 
                           double prevalence_threshold, int degree) {
            
            std::vector<std::vector<std::string>> prevalent;
            
            // iterate over each candidate pattern
            for (const auto& currentPattern : candidatePatterns) {
                // extract the basePattern (ex: basePattern of (A,B,C) is (A,B))
                std::vector<std::string> basePattern;
                for (int j = 0; j < degree - 1; j++) {
                    basePattern.push_back(currentPattern[j]);
                }
                std::string lastFeature = currentPattern[degree - 1];
                this->instance_table[currentPattern] = {};
                this->hashmap[currentPattern] = {};
                for (const auto& f : currentPattern) {
                    this->hashmap[currentPattern][f] = {};
                }
                
                // find the basePattern in the instance table
                auto colocTableIt = this->instance_table.find(basePattern);
                std::map<std::vector<int>, std::vector<int>>& colocTable = colocTableIt->second;
                
                // iterate over each entry in the value of the colocTable map
                for (const auto& entry : colocTable) {
                    const std::vector<int>& key = entry.first;
                    std::set<int> commonLastNeighbors;
                    /* iterate over each entry in the key.  for example, if entry is (A1, B1) -> [C1, C2, C3],
                       then we are iterating over (A1, B1) */
                    for (int instanceID : key) {
                        // find the instanceID in the star_neighbors structure
                        auto star_neighbor_it = this->star_neighbors.find(instanceID);
                        if (commonLastNeighbors.empty()) {
                            // find the star neighbors of instanceID that are of type lastFeature
                            std::vector<int> temp_vector = findNeighborsInRange(star_neighbor_it->second,
                                                                         this->featureInfo[lastFeature][1],
                                                                         this->featureInfo[lastFeature][2]);
                            commonLastNeighbors.insert(temp_vector.begin(), temp_vector.end());
                        } else {
                            // find the star neighbors of instanceID that are of type lastFeature
                            std::vector<int> temp_vector = findNeighborsInRange(star_neighbor_it->second, 
                                                                                  this->featureInfo[lastFeature][1],
                                                                                  this->featureInfo[lastFeature][2]);
                            std::set<int> temp_commonLastNeighbors(temp_vector.begin(), temp_vector.end());
                            // perform intersection of commonLastNeighbors and temp_commonLastNeighbors
                            std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                                  temp_commonLastNeighbors.begin(), temp_commonLastNeighbors.end(),
                                                  std::inserter(commonLastNeighbors, commonLastNeighbors.begin()));
                        }
                    }
                    /* iterate over each entry in the value portion of the map.  for example, if key is (A1, B1), then
                       we are iterating over its value [C1, C2, C3] */
                    for (int n : colocTable[key]) {
                        // find n in the star_neighbors structure
                        auto star_neighbor_it = this->star_neighbors.find(n);
                        // find the star neighbors of n of type lastFeature
                        std::vector<int> temp_vect = findNeighborsInRange(star_neighbor_it->second,
                                                                          this->featureInfo[lastFeature][1],
                                                                          this->featureInfo[lastFeature][2]);
                        std::set<int> temp_neighbors(temp_vect.begin(), temp_vect.end());
                        std::set<int> neighbors;
                        std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                              temp_neighbors.begin(), temp_neighbors.end(),
                                              std::inserter(neighbors, neighbors.begin()));

                        if (!neighbors.empty()) {
                            // create new key for the instance_table
                            std::vector<int> new_key = key;
                            new_key.push_back(n);
                            std::vector<int> intersectionVec(neighbors.begin(), neighbors.end());
                            this->instance_table[currentPattern][new_key] = intersectionVec;

                            // update hashmap
                            for (size_t k = 0; k < new_key.size(); k++) {
                                this->hashmap[currentPattern][currentPattern[k]].insert(new_key[k]);
                            }
                            this->hashmap[currentPattern][lastFeature].insert(neighbors.begin(), neighbors.end());
                        }
                    }
                }
            
                // calculate participation ratios
                std::vector<double> pr;
                for (int m = 0; m < degree; m++) {
                    std::string f = currentPattern[m];
                    double ratio = static_cast<double>(this->hashmap[currentPattern][f].size()) 
                        / this->featureInfo[f][0];
                    pr.push_back(ratio);
                }
                
                // calculate participation index
                double PI = *std::min_element(pr.begin(), pr.end());
                
                // determine if pattern is prevalent or not
                if (PI >= prevalence_threshold) {
                    prevalent.push_back(currentPattern);
                }
            }

            std::cout << "Degree " << degree << " Prevalent Patterns for Sub-Region " << this->subregion_id << ":"<< std::endl;
            std::ofstream patterns_file("size" + std::to_string(degree) + "_patterns_subregion" 
                                        + std::to_string(this->subregion_id) + ".txt");
            
            for (auto& patternVec : prevalent) {
                std::cout << "(";
                patterns_file << "(";
                for (size_t i = 0; i < patternVec.size(); i++) {
                    std::cout << patternVec[i];
                    patterns_file << patternVec[i];
                    if (i < patternVec.size() - 1) {
                        std::cout << ", ";
                        patterns_file << ", ";
                    }
                }
                std::cout << ")" << std::endl;
                patterns_file << ")" << std::endl;
            }
            patterns_file.close();
            return prevalent;
        }
    };
    std::vector<SubRegion> subregions;

    // this class holds all information pertaining to the border region
    class Border {
    public:
        int border_id;  // unique identifier
        std::map<std::string, std::vector<int>> featureInfo;
        std::map<int, std::vector<int>> star_neighbors;
        std::vector<std::vector<std::string>> size2_patterns;
        std::map<std::vector<std::string>, std::map<std::vector<int>, std::vector<int>>> instance_table;
        std::map<std::vector<std::string>, std::map<std::string, std::set<int>>> hashmap;
        
        Border(int number) {
            border_id = number;  // assign unique identifier to each border region
        }
        
        // read featureInfo from csv file
        void read_featureInfo() {
            std::ifstream featureInfo_file("required_files/border_featureInfo/featureInfo.csv");
        
            if (!featureInfo_file.is_open()) {
                std::cerr << "ERROR OPENING BORDER " << this->border_id << " FEATURE INFO FILE!" << std::endl;
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
                this->featureInfo[key] = values;
            }
            featureInfo_file.close();
        }
        
        void read_star_neighbors() {
            // read the star_neighbors from a csv file
            std::ifstream star_neighbors_file("required_files/border_starNeighbors/starNeighbors.csv");
            if (!star_neighbors_file.is_open()) {
                std::cerr << "ERROR OPENING BORDER " << this->border_id << " STAR NEIGHBORS FILE!" << std::endl;
            }
            std::string line;
            std::getline(star_neighbors_file, line);
            while (std::getline(star_neighbors_file, line)) {
                std::stringstream ss(line);
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
                this->star_neighbors[key] = values;
            }
            star_neighbors_file.close();
        }
    
        // generate degree 2 candidate patterns
        std::vector<std::vector<std::string>> generate_size2_combos() {
            std::vector<std::string> features;
            for (const auto &entry : this->featureInfo) {
                features.push_back(entry.first);
            }
            std::vector<std::vector<std::string>> size2_candidatePatterns;
            for (size_t i = 0; i < features.size(); i++) {
                for (size_t j = i + 1; j < features.size(); j++) {
                    size2_candidatePatterns.push_back({features[i], features[j]});
                }
            }
            return size2_candidatePatterns;
        }
        
        // log(n) search that finds neighbors within a specific start and end range
        std::vector<int> findNeighborsInRange(const std::vector<int>& arr, int x, int y) {
            auto start_it = std::lower_bound(arr.begin(), arr.end(), x);
            auto end_it = std::upper_bound(arr.begin(), arr.end(), y);
            return std::vector<int>(start_it, end_it);
        }
        
        // calculate degree 2 prevalent patterns
        std::vector<std::vector<std::string>> degree2Processing(
            std::vector<std::vector<std::string>> candidatePatterns, int candidatePatterns_size, 
            double prevalence_threshold) {   
            std::vector<std::vector<std::string>> size2_patterns;
            
            for (const auto& coloc : candidatePatterns) {
                std::string first_feature = coloc[0];
                std::string second_feature = coloc[1];
                
                // determine first feature start and end
                auto feature1 = this->featureInfo.find(first_feature);
                std::vector<int> values1 = feature1->second;
                int first_feature_start = values1[1];
                int first_feature_end = values1[2];
                // determine second feature start and end
                auto feature2 = this->featureInfo.find(second_feature);
                std::vector<int> values2 = feature2->second;
                int second_feature_start = values2[1];
                int second_feature_end = values2[2];
                
                std::vector<std::string> coloc_key = {first_feature, second_feature};
                instance_table[coloc_key] = {};
                this->hashmap[coloc_key] = {};
                this->hashmap[coloc_key][first_feature] = {};
                this->hashmap[coloc_key][second_feature] = {};
                
                for (int index = first_feature_start; index <= first_feature_end; index++) {
                    auto star_neighbor_it = this->star_neighbors.find(index);
                    std::vector<int> neighbors = 
                        findNeighborsInRange(star_neighbor_it->second, second_feature_start,
                                            second_feature_end);
                    if (!neighbors.empty()) {
                        std::vector<int> index_tuple = {index};
                        this->instance_table[coloc_key][index_tuple] = neighbors;
                        this->hashmap[coloc_key][first_feature].insert(index);
                        for (int neighbor : neighbors) {
                            this->hashmap[coloc_key][second_feature].insert(neighbor);
                        }
                    }
                }
                
                // calculate participation ratios
                double pr_first_feature = 
                    static_cast<double>(this->hashmap[coloc_key][first_feature].size()) /
                    this->featureInfo[first_feature][0];
                double pr_second_feature = 
                    static_cast<double>(this->hashmap[coloc_key][second_feature].size()) /
                    this->featureInfo[second_feature][0];
                
                // calculate participation index
                double PI = 0.0;
                if (pr_first_feature < pr_second_feature) {
                    PI = pr_first_feature;
                } else {
                    PI = pr_second_feature;
                }

                // determine if pattern is prevalent or not
                if (PI >= prevalence_threshold) {
                       size2_patterns.push_back(coloc_key);
                }
            }
            
            std::cout << "Degree 2 Prevalent Patterns for Border " << this->border_id << ":" << std::endl;
            for (auto i : size2_patterns) {
                std::cout << "(" << i[0] << ", " << i[1] << ")" << std::endl;
            }
            return size2_patterns;
        }
    };
    std::vector<Border> borders;

    // this class holds all information pertaining to the entire region
    class Region {
    public:
        std::map<std::vector<std::string>, std::map<std::vector<int>, std::vector<int>>> instance_table;
        std::map<std::vector<std::string>, std::map<std::string, std::set<int>>> hashmap;
        
        // generate degree 2 candidate patterns
        std::vector<std::vector<std::string>> generate_size2_combos() {
            std::set<std::vector<std::string>> patterns;
            for (const auto& entry : this->hashmap) {
                patterns.insert(entry.first);
            }
            std::vector<std::vector<std::string>> size2_candidatePatterns(patterns.begin(), patterns.end());
            return size2_candidatePatterns;
        }
        
        // calculate degree 2 prevalent patterns
        std::vector<std::vector<std::string>> degree2Processing(std::vector<std::vector<std::string>> candidatePatterns, 
                                                                int candidatePatterns_size, double prevalence_threshold,
                                                                int number_subregions) {
            std::vector<std::vector<std::string>> prevalent;
            
            // iterate over each candidate pattern
            for (const auto& coloc : candidatePatterns) {
                std::string first_feature = coloc[0];
                std::string second_feature = coloc[1];
                int total_first_feature = 0;
                int total_second_feature = 0;
                
                // total up the number of features in the sub-regions
                for (int i = 0; i < number_subregions; i++) {
                    auto feature1 = subregions[i].featureInfo.find(first_feature);
                    std::vector<int> values1 = feature1->second;
                    total_first_feature += values1[0];
                    
                    auto feature2 = subregions[i].featureInfo.find(second_feature);
                    std::vector<int> values2 = feature2->second;
                    total_second_feature += values2[0];
                }
                
                // calculate the participation ratio
                double pr_first_feature = static_cast<double>(this->hashmap[coloc][first_feature].size()) / total_first_feature;
                double pr_second_feature = static_cast<double>(this->hashmap[coloc][second_feature].size()) / total_second_feature;
                
                // calculate the participation index
                double PI = 0.0;
                if (pr_first_feature < pr_second_feature) {
                    PI = pr_first_feature;
                } else {
                    PI = pr_second_feature;
                }

                // determine if patten in prevalent or not
                if (PI >= prevalence_threshold) {
                       prevalent.push_back(coloc);
                }
            }
            
            std::cout << "Degree 2 Prevalent Patterns for Entire Region:" << std::endl;
            std::ofstream size2_file("size2_patterns_entire_region.txt");
            for (auto i : prevalent) {
                std::cout << "(" << i[0] << ", " << i[1] << ")" << std::endl;
                size2_file << "(" << i[0] << ", " << i[1] << ")" << std::endl;
            }
            size2_file.close();
            return prevalent;
        }
        
        // helper function for getCandidatePatterns
        void generateCombinations(const std::vector<std::string>& features, int degree, 
                                  std::vector<std::vector<std::string>>& result, 
                                  std::vector<std::string>& current, int start) {
            if (current.size() == degree) {
                result.push_back(current);
                return;
            }
            for (int i = start; i < features.size(); i++) {
                current.push_back(features[i]);
                generateCombinations(features, degree, result, current, i + 1);
                current.pop_back();
            }
        }

        /* helper function for getCandidatePatterns: checks if all (degree-1)-subpatterns of
        a pattern are in the prevalent patterns */
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

        // generate candidate patterns of degree > 2
        std::vector<std::vector<std::string>> getCandidatePatterns(
            const std::vector<std::vector<std::string>>& prevalentPattern, int degree, std::vector<std::string> features) {
            std::set<std::vector<std::string>> prevalentPatterns(prevalentPattern.begin(), prevalentPattern.end());
            
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
        
        // log(n) search that finds neighbors within a specific start and end range
        std::vector<int> findNeighborsInRange(const std::vector<int>& arr, int x, int y) {
            auto start_it = std::lower_bound(arr.begin(), arr.end(), x);
            auto end_it = std::upper_bound(arr.begin(), arr.end(), y);
            return std::vector<int>(start_it, end_it);
        }
        
        // calculate prevalent patterns of degree > 2
        std::vector<std::vector<std::string>> colocationGeneral(std::vector<std::vector<std::string>> 
                                                                candidatePatterns, 
                                                                int candidatePatterns_size, 
                                                                double prevalence_threshold,
                                                                int degree, int number_subregions) {
            std::vector<std::vector<std::string>> prevalent;
            
            // iterate over each candidate pattern
            for (const auto& currentPattern : candidatePatterns) {
                // the basePattern of (A, B, C) would be (A, B)
                std::vector<std::string> basePattern;
                for (int j = 0; j < degree - 1; j++) {
                    basePattern.push_back(currentPattern[j]);
                }
                
                // the lastFeature of (A, B, C) would be C
                std::string lastFeature = currentPattern[degree - 1];
                this->instance_table[currentPattern] = {};
                this->hashmap[currentPattern] = {};
                
                for (const auto& f : currentPattern) {
                    this->hashmap[currentPattern][f] = {};
                }
                
                // find basePattern instance_table
                auto colocTableIt = this->instance_table.find(basePattern);
                std::map<std::vector<int>, std::vector<int>>& colocTable = colocTableIt->second;
                
                // iterate over each entry in colocTable
                for (const auto& entry : colocTable) {
                    // if entry is (A1, B1) -> [C1, C2, C3], then the key would be (A1, B1)
                    const std::vector<int>& key = entry.first;
                    std::set<int> commonLastNeighbors;
                    // itereate over each instanceID in the key
                    for (int instanceID : key) {
                        std::set<int> lastFeatureNeighbors;
                        // iterate over each subregion
                        for (auto subregion : subregions) {
                            // get the start and end indices of the lastFeature
                            int lastFeatureStart = subregion.featureInfo[lastFeature][1];
                            int lastFeatureEnd = subregion.featureInfo[lastFeature][2];
                            
                            // find the instanceID in the subregion's star_neighbors structure
                            auto subregion_star_neighbor_it = subregion.star_neighbors.find(instanceID);
                            if (subregion_star_neighbor_it != subregion.star_neighbors.end()) {
                                // make a copy of the instanceID's star_neighborhood
                                std::vector<int> subregion_star_neighbor_list(subregion_star_neighbor_it->second.begin(), 
                                                   subregion_star_neighbor_it->second.end());
                                
                                // find the star_neighbors of instanceID that are of type lastFeature in the subregion's star_neighbor structure
                                std::vector<int> neighborsInRange_vector = this->findNeighborsInRange(subregion_star_neighbor_list,
                                                                                                      lastFeatureStart, lastFeatureEnd);
                                // convert the vector into a set
                                std::set<int> neighborsInRange(neighborsInRange_vector.begin(), neighborsInRange_vector.end());
                                lastFeatureNeighbors.insert(neighborsInRange.begin(), neighborsInRange.end());
                            }
                            // note: this framework can only handle 1 border region, that is why it is hardcoded as borders[0]
                            // find instanceID in the border regions star_neighbors
                            auto border_star_neighbor_it = borders[0].star_neighbors.find(instanceID);
                            if (border_star_neighbor_it != borders[0].star_neighbors.end()) {
                                std::vector<int> border_star_neighbor_list(border_star_neighbor_it->second.begin(), 
                                                   border_star_neighbor_it->second.end());
                                // find star_neighbors of instanceID of type lastFeature in the border region's star_neighbors structure
                                std::vector<int> neighborsInRange_vector = this->findNeighborsInRange(border_star_neighbor_list, lastFeatureStart, 
                                                                                      lastFeatureEnd);
                                // convert the vector to a set
                                std::set<int> neighborsInRange(neighborsInRange_vector.begin(), neighborsInRange_vector.end());
                                lastFeatureNeighbors.insert(neighborsInRange.begin(), neighborsInRange.end());
                            }
                        }
                        if (commonLastNeighbors.empty()) {
                            commonLastNeighbors = lastFeatureNeighbors;
                        } else {
                            std::set<int> intersection;
                            std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                                  lastFeatureNeighbors.begin(), lastFeatureNeighbors.end(),
                                                  std::inserter(intersection, intersection.begin()));
                            commonLastNeighbors = intersection;
                        }
                    }
                    // if key is (A1, B1), then we are iterating over [C1, C2, C3]
                    for (int n : colocTable[key]) {
                        std::set<int> all_neighbors;
                        // iterate over each subregion
                        for (auto subregion : subregions) {
                            // get the start and end indices of lastFeature
                            int lastFeatureStart = subregion.featureInfo[lastFeature][1];
                            int lastFeatureEnd = subregion.featureInfo[lastFeature][2];
                            std::set<int> neighbors_subregion;
                            // find n in the star_neighbors structure
                            auto subregion_star_neighbor_it = subregion.star_neighbors.find(n);
                            if (subregion_star_neighbor_it != subregion.star_neighbors.end()) {
                                std::vector<int> subregion_star_neighbor_list(subregion_star_neighbor_it->second.begin(), 
                                                   subregion_star_neighbor_it->second.end());
                                // find the star_neighbors of n of type lastFeature in the subregion's star_neighbors structure
                                std::vector<int> neighborsInRange_vector = this->findNeighborsInRange(subregion_star_neighbor_list, lastFeatureStart, 
                                                                                      lastFeatureEnd);
                                std::set<int> neighborsInRange(neighborsInRange_vector.begin(), neighborsInRange_vector.end());
                                std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                                     neighborsInRange.begin(), neighborsInRange.end(),
                                                     std::inserter(neighbors_subregion, neighbors_subregion.begin()));
                            }
                            std::set<int> neighbors_border;
                            // find n in the border region's star_neighbors structure
                            auto border_star_neighbor_it = borders[0].star_neighbors.find(n);
                            if (border_star_neighbor_it != borders[0].star_neighbors.end()) {
                                std::vector<int> border_star_neighbor_list(border_star_neighbor_it->second.begin(), 
                                                   border_star_neighbor_it->second.end());
                                // find the star_neighbors of n in the border region's star_neighbors structure
                                std::vector<int> neighborsInRange_vector = this->findNeighborsInRange(border_star_neighbor_list, lastFeatureStart, 
                                                                                      lastFeatureEnd);
                                std::set<int> neighborsInRange(neighborsInRange_vector.begin(), neighborsInRange_vector.end());
                                std::set_intersection(commonLastNeighbors.begin(), commonLastNeighbors.end(),
                                                     neighborsInRange.begin(), neighborsInRange.end(),
                                                     std::inserter(neighbors_border, neighbors_border.begin()));
                            }
                            all_neighbors.insert(neighbors_subregion.begin(), neighbors_subregion.end());
                            all_neighbors.insert(neighbors_border.begin(), neighbors_border.end());
                        }
                        if (!all_neighbors.empty()) {
                            // create a new key for the instance_table
                            std::vector<int> new_key = key;
                            new_key.push_back(n);
                            std::vector<int> intersectionVec(all_neighbors.begin(), all_neighbors.end());
                            this->instance_table[currentPattern][new_key] = intersectionVec;

                            // update hashmap
                            for (size_t k = 0; k < new_key.size(); k++) {
                                this->hashmap[currentPattern][currentPattern[k]].insert(new_key[k]);
                            }
                            this->hashmap[currentPattern][lastFeature].insert(all_neighbors.begin(), all_neighbors.end());
                        }
                    }
                }
                
                // calculate participation ratio
                std::vector<double> pr;
                for (int m = 0; m < degree; m++) {
                    std::string f = currentPattern[m];
                    int total_count = 0;
                    for (auto subregion : subregions) {
                        total_count += subregion.featureInfo[f][0];
                    }
                    double ratio = static_cast<double>(this->hashmap[currentPattern][f].size()) 
                        / total_count;
                    pr.push_back(ratio);
                }
                
                // calculate participation index
                double PI = *std::min_element(pr.begin(), pr.end());
                
                // determine if pattern is prevalent or not
                if (PI >= prevalence_threshold) {
                    prevalent.push_back(currentPattern);
                }
            }
            
            std::cout << "Degree " << degree << " Prevalent Patterns for Entire Region:" << std::endl;
            std::ofstream patterns_file("size" + std::to_string(degree) + "_patterns_entire_region.txt");
            for (auto& patternVec : prevalent) {
                std::cout << "(";
                patterns_file << "(";
                for (size_t i = 0; i < patternVec.size(); i++) {
                    std::cout << patternVec[i];
                    patterns_file << patternVec[i];
                    if (i < patternVec.size() - 1) {
                        std::cout << ", ";
                        patterns_file << ", ";
                    }
                }
                std::cout << ")" << std::endl;
                patterns_file << ")" << std::endl;
            }
            patterns_file.close();
            return prevalent;
        }
    };
    Region region;

    // processing for all the sub-regions
    void subregion_main(int number_subregions, double prevalence_threshold) {
        for (int i = 0; i < number_subregions; i++) {
            subregions.push_back(SubRegion(i));
        }
        
        for (int i = 0; i < number_subregions; i++) {
            subregions[i].read_featureInfo();
            subregions[i].read_star_neighbors();
            std::vector<std::vector<std::string>> size2_candidatePatterns = 
                subregions[i].generate_size2_combos();
            subregions[i].size2_patterns = 
                subregions[i].degree2Processing(size2_candidatePatterns,
                                                size2_candidatePatterns.size(), 
                                                prevalence_threshold);
            
            int degree = 3;
            std::vector<std::vector<std::string>> candidatePatterns =
                subregions[i].getCandidatePatterns(subregions[i].size2_patterns, degree);
            while (!candidatePatterns.empty()) {
                std::vector<std::vector<std::string>> prevalent_patterns = 
                subregions[i].colocationGeneral(candidatePatterns, candidatePatterns.size(), 
                           prevalence_threshold, degree);
                degree += 1;
                if (prevalent_patterns.size() == 0) {
                    break;
                }
                candidatePatterns = subregions[i].getCandidatePatterns(prevalent_patterns,
                                                                       degree);
            }
        }
    }

    // processing for the border region
    void border_main(int number_borders, double prevalence_threshold) {
        for (int i = 0; i < number_borders; i++) {
            borders.push_back(Border(i));
        }
        
        for (int i = 0; i < number_borders; i++) {
            borders[i].read_featureInfo();
            borders[i].read_star_neighbors();
            std::vector<std::vector<std::string>> size2_candidatePatterns = 
                borders[i].generate_size2_combos();
            borders[i].size2_patterns = 
                borders[i].degree2Processing(size2_candidatePatterns, 
                                             size2_candidatePatterns.size(), 
                                             prevalence_threshold);
        } 
    } 

    void update_border_info(int* ids, int ids_size, int i) {
        // update the border hashmap so that it has the original IDS not the indicies
        std::map<std::vector<std::string>, std::map<std::string, std::set<int>>> temp_hash;
        for (const auto& outer_pair : borders[i].hashmap) {
            const std::vector<std::string>& outer_key = outer_pair.first;
            const auto& inner_map = outer_pair.second;
            std::map<std::string, std::set<int>> temp_inner_map;
            for (const auto& inner_pair : inner_map) {
                const std::string& inner_key = inner_pair.first;
                const std::set<int>& inner_set = inner_pair.second;
                std::set<int> temp_inner_set;
                for (int index : inner_set) {
                    if (index < ids_size) {
                        temp_inner_set.insert(ids[index]);
                    }
                }
                temp_inner_map[inner_key] = temp_inner_set;
            }
            temp_hash[outer_key] = temp_inner_map;
        }
        borders[i].hashmap = temp_hash;
    
        // update the border instance table so that it has the original IDS not indicies
        std::map<std::vector<std::string>, std::map<std::vector<int>, std::vector<int>>> temp_instance_table;

        for (const auto& outer_pair : borders[i].instance_table) {
            const std::vector<std::string>& outer_key = outer_pair.first;
            const auto& inner_map = outer_pair.second;
            std::map<std::vector<int>, std::vector<int>> temp_inner_map;
            for (const auto& inner_pair : inner_map) {
                const std::vector<int>& inner_key = inner_pair.first;
                const std::vector<int>& inner_value = inner_pair.second;
                std::vector<int> new_inner_key;
                for (int index : inner_key) {
                    if (index < ids_size) {
                        new_inner_key.push_back(ids[index]);
                    }
                }
                std::vector<int> new_inner_value;
                for (int index : inner_value) {
                    if (index < ids_size) {
                        new_inner_value.push_back(ids[index]);
                    }
                }
                temp_inner_map[new_inner_key] = new_inner_value;
            }
            temp_instance_table[outer_key] = temp_inner_map;
        }
        borders[i].instance_table = temp_instance_table;
        
        // update the border star neighbors so that is has the original IDS not indices
        std::map<int, std::vector<int>> temp_star_neighbors;
        
        for (const auto& entry : borders[i].star_neighbors) {
            int original_id = ids[entry.first]; // Get the original ID corresponding to the index
            const std::vector<int>& neighbors_indices = entry.second;
            std::vector<int> original_neighbors;
            for (int index : neighbors_indices) {
                if (index < ids_size) {
                    original_neighbors.push_back(ids[index]); // Get the original ID for each neighbor
                }
            }
            temp_star_neighbors[original_id] = original_neighbors;
        }
        borders[i].star_neighbors = temp_star_neighbors;        
    }

    // combine the hashmaps of the subregions and border region
    void combine_hashmaps(int number_subregions, int number_borders) {
        std::set<std::vector<std::string>> keys_needed;
        
        // determine which keys from the subregions and borders need to be combined to form the regional hashmap
        for (int i = 0; i < number_subregions; i++) {
            std::set<std::vector<std::string>> temp_vect(subregions[i].size2_patterns.begin(), 
                                                         subregions[i].size2_patterns.end());  
            if (i == 0) {
                keys_needed = temp_vect;
            } else {
                std::set<std::vector<std::string>> temp_result;
                std::set_union(keys_needed.begin(), keys_needed.end(),
                               temp_vect.begin(), temp_vect.end(),
                               std::inserter(temp_result, temp_result.begin()));
                keys_needed = temp_result;
            }
        }
        
        for (int i = 0; i < number_borders; i++) {
            std::set<std::vector<std::string>> temp_vect(borders[i].size2_patterns.begin(), 
                                                         borders[i].size2_patterns.end());
            std::set<std::vector<std::string>> temp_result;
            std::set_union(keys_needed.begin(), keys_needed.end(),
                           temp_vect.begin(), temp_vect.end(),
                           std::inserter(temp_result, temp_result.begin()));
            keys_needed = temp_result;
        }
        
        // combine the hashmaps
        for (const auto& key : keys_needed) {
            for (int i = 0; i < number_subregions; i++) {
                auto it = subregions[i].hashmap.find(key);
                if (it != subregions[i].hashmap.end()) {
                    for (const auto& inner_pair : it->second) {
                        const auto& inner_key = inner_pair.first;
                        const auto& inner_value = inner_pair.second;
                        auto& region_inner_map = region.hashmap[key][inner_key];
                        region_inner_map.insert(inner_value.begin(), inner_value.end());
                    }
                }
            }

            for (int i = 0; i < number_borders; i++) {
                auto it = borders[i].hashmap.find(key);
                if (it != borders[i].hashmap.end()) {
                    for (const auto& inner_pair : it->second) {
                        const auto& inner_key = inner_pair.first;
                        const auto& inner_value = inner_pair.second;
                        auto& region_inner_map = region.hashmap[key][inner_key];
                        region_inner_map.insert(inner_value.begin(), inner_value.end());
                    }
                }
            }
        }  
    }

    // combine the instance tables of the subregions and the border region
    void combine_instance_tables(int number_subregions, int number_borders) {
        std::set<std::vector<std::string>> keys_needed;
        
        // determine which keys from the subregions and borders need to be combined to form the regional hashmap
        for (int i = 0; i < number_subregions; i++) {
            std::set<std::vector<std::string>> temp_vect(subregions[i].size2_patterns.begin(), 
                                                         subregions[i].size2_patterns.end());  
            if (i == 0) {
                keys_needed = temp_vect;
            } else {
                std::set<std::vector<std::string>> temp_result;
                std::set_union(keys_needed.begin(), keys_needed.end(),
                               temp_vect.begin(), temp_vect.end(),
                               std::inserter(temp_result, temp_result.begin()));
                keys_needed = temp_result;
            }
        }
        
        for (int i = 0; i < number_borders; i++) {
            std::set<std::vector<std::string>> temp_vect(borders[i].size2_patterns.begin(), 
                                                         borders[i].size2_patterns.end());
            std::set<std::vector<std::string>> temp_result;
            std::set_union(keys_needed.begin(), keys_needed.end(),
                           temp_vect.begin(), temp_vect.end(),
                           std::inserter(temp_result, temp_result.begin()));
            keys_needed = temp_result;
        }
        
        // Combine the instance tables
        for (const auto& key : keys_needed) {
            for (int i = 0; i < number_subregions; i++) {
                auto it = subregions[i].instance_table.find(key);
                if (it != subregions[i].instance_table.end()) {
                    for (const auto& inner_pair : it->second) {
                        const auto& inner_key = inner_pair.first;
                        const auto& inner_value = inner_pair.second;
                        auto& region_inner_map = region.instance_table[key];
                        if (region_inner_map.find(inner_key) != region_inner_map.end()) {
                            region_inner_map[inner_key].insert(region_inner_map[inner_key].end(), inner_value.begin(),
                                                               inner_value.end());
                        } else {
                            region_inner_map[inner_key] = inner_value;
                        }
                    }
                }
            }

            for (int i = 0; i < number_borders; i++) {
                auto it = borders[i].instance_table.find(key);
                if (it != borders[i].instance_table.end()) {
                    for (const auto& inner_pair : it->second) {
                        const auto& inner_key = inner_pair.first;
                        const auto& inner_value = inner_pair.second;
                        auto& region_inner_map = region.instance_table[key];
                        if (region_inner_map.find(inner_key) != region_inner_map.end()) {
                            region_inner_map[inner_key].insert(region_inner_map[inner_key].end(), inner_value.begin(),
                                                               inner_value.end());
                        } else {
                            region_inner_map[inner_key] = inner_value;
                        }
                    }
                }
            }
        }
    }

    // processing for the entire region
    void region_main(int number_subregions, double prevalence_threshold, char** features_ptr, int features_size) {
        std::vector<std::string> features(features_ptr, features_ptr + features_size);
        std::vector<std::vector<std::string>> size2_candidatePatterns = region.generate_size2_combos();
        std::vector<std::vector<std::string>> prevalent_patterns = region.degree2Processing(size2_candidatePatterns, 
                                                                size2_candidatePatterns.size(), prevalence_threshold,
                                                                                           number_subregions);
        int degree = 3;
        std::vector<std::vector<std::string>> candidatePatterns = region.getCandidatePatterns(prevalent_patterns,
                                                                                              degree, features);
        while (!candidatePatterns.empty()) {
            std::vector<std::vector<std::string>> prevalent_patterns = region.colocationGeneral(candidatePatterns, 
                                                                    candidatePatterns.size(), prevalence_threshold, 
                                                                    degree, number_subregions);
            degree += 1;
            if (prevalent_patterns.size() == 0) {
                break;
            }
            candidatePatterns = region.getCandidatePatterns(prevalent_patterns, degree, features);
        }
    }
}
