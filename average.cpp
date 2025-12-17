#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <string>

using namespace std;

struct AverageData {
    vector<double> columnSums;
    int count;
};

int main() {
    // Configurable parameters
    string inputFilename = "results_4x4x4.dat";
    string outputFilename = "averaged_results_4x4x4.dat";

    // Open input file
    ifstream inputFile(inputFilename);
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open file " << inputFilename << endl;
        return 1;
    }

    // Read and save header
    string header;
    getline(inputFile, header);

    // Store data grouped by first column (Beta)
    map<double, AverageData> dataMap;

    string line;
    int totalRows = 0;
    int columnCount = 0;

    // Read data line by line
    while (getline(inputFile, line)) {
        if (line.empty()) continue;

        stringstream ss(line);
        vector<double> row;
        double value;

        // Read all values in the row
        while (ss >> value) {
            row.push_back(value);
        }

        // Determine column count from first data row
        if (totalRows == 0) {
            columnCount = row.size();
        }

        // Check column count consistency
        if (row.size() != columnCount) {
            cerr << "Warning: Column count mismatch in line "
                 << totalRows + 2 << ". Expected: " << columnCount
                 << ", Found: " << row.size() << endl;
            continue;
        }

        // First column is the key for grouping (Beta)
        double key = row[0];

        // Initialize if key is new
        if (dataMap.find(key) == dataMap.end()) {
            AverageData avgData;
            avgData.columnSums = vector<double>(columnCount, 0.0);
            avgData.count = 0;
            dataMap[key] = avgData;
        }

        // Sum values for each column
        for (int i = 0; i < columnCount; i++) {
            dataMap[key].columnSums[i] += row[i];
        }
        dataMap[key].count++;

        totalRows++;
    }

    inputFile.close();

    // Open output file
    ofstream outputFile(outputFilename);
    if (!outputFile.is_open()) {
        cerr << "Error: Could not create output file " << outputFilename << endl;
        return 1;
    }

    // Write header to output file
    outputFile << header << endl;

    // Write averaged data to output file
    outputFile << fixed << setprecision(10);

    cout << "Processing complete. Writing averaged data to " << outputFilename << endl;
    cout << "Number of data rows processed: " << totalRows << endl;
    cout << "Number of unique Beta values: " << dataMap.size() << endl;
    cout << "Number of columns: " << columnCount << endl;

    for (const auto& pair : dataMap) {
        double beta = pair.first;
        const AverageData& avgData = pair.second;

        // Write beta value
        outputFile << beta;

        // Write averaged values for other columns
        for (int i = 1; i < columnCount; i++) {
            double average = avgData.columnSums[i] / avgData.count;
            outputFile << "\t" << average;
        }
        outputFile << endl;
    }

    outputFile.close();

    // Optional: Display first few averaged values for verification
    cout << "\nFirst 3 averaged values for verification:" << endl;
    cout << header << endl;

    int count = 0;
    for (const auto& pair : dataMap) {
        if (count >= 15) break;

        double beta = pair.first;
        const AverageData& avgData = pair.second;

        cout << fixed << setprecision(6) << beta;
        for (int i = 1; i < columnCount; i++) {
            double average = avgData.columnSums[i] / avgData.count;
            cout << "\t" << average;
        }
        cout << endl;
        count++;
    }

    return 0;
}
