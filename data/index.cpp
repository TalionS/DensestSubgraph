#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

int main() {
    std::ifstream input_file("/home/yy/DensestSubgraph/data/moreno_innovation/out.moreno_innovation_innovation");
    std::ofstream output_file("MI.txt");

    if (!input_file.is_open()) {
        std::cerr << "Failed to open input file." << std::endl;
        return 1;
    }

    if (!output_file.is_open()) {
        std::cerr << "Failed to open output file." << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(input_file, line)) {
        std::vector<int> numbers;
        std::istringstream iss(line);
        int num;
        while (iss >> num) {
            numbers.push_back(num - 1);
        }
        for (auto num: numbers)
            std::cout << num;
        std::cout << std::endl;
        for (size_t i = 0; i < numbers.size(); i++) {
            output_file << numbers[i];
            if (i < numbers.size() - 1) {
                output_file << " ";
            }
        }
        output_file << "\n";
    }

    input_file.close();
    output_file.close();

    return 0;
}