#include "csvExporter.h"
#include <fstream>
#include <sstream>

// Constructor
DataExporter::DataExporter() {}

// Method to export matrix to CSV file with column names
bool DataExporter::exportToCSV(const std::vector<cPoint>& matrix, const std::string& fileName, const std::vector<std::string>& columnNames) {
    std::ofstream file(fileName);
    if (!file.is_open()) {
        return false;
    }

    if (!columnNames.empty()) {
        file << columnNamesToCSV(columnNames) << "\n";
    }

    for (const auto& row : matrix) {
        file << rowToCSV(row) << "\n";
    }

    file.close();
    return true;
}

// Helper method to convert matrix row to CSV format
std::string DataExporter::rowToCSV(const cPoint& point) {
    std::ostringstream oss;
    std::vector<double> row = { static_cast<double>(point.getProperties().index),point.x,point.y,point.mu,point.getProperties().flow_angle, point.getProperties().mach,point.getProperties().pm_angle,static_cast<double>(point.getProperties().type) };
    for (size_t i = 0; i < row.size(); ++i) {
        if (i != 0) {
            oss << ",";
        }
        oss << row[i];
    }
    return oss.str();
}

// Helper method to convert column names to CSV format
std::string DataExporter::columnNamesToCSV(const std::vector<std::string>& columnNames) {
    std::ostringstream oss;
    for (size_t i = 0; i < columnNames.size(); ++i) {
        if (i != 0) {
            oss << ",";
        }
        oss << columnNames[i];
    }
    return oss.str();
}
