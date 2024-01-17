#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <unistd.h>
#include <limits.h> // 用于 PATH_MAX

int main() {
    std::ifstream input_file("data/1111.txt");
    std::ofstream output_file("data/MR.txt");

    char currentPath[PATH_MAX];
    if (getcwd(currentPath, sizeof(currentPath)) != NULL) {
        std::cout << "Current path is: " << currentPath << std::endl;
    } else {
        std::cerr << "Error getting current directory" << std::endl;
        return 1;
    }


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
        for (int i = 0; i < 2; i++){
            iss >> num;
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