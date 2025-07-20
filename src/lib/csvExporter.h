#ifndef DATAEXPORTER_H
#define DATAEXPORTER_H

#include <vector>
#include <string>
#include "point.h"

class DataExporter {
public:
    DataExporter();

    // Method to export matrix to CSV file with column names
    bool exportToCSV(const std::vector<cPoint>& matrix, const std::string& fileName, const std::vector<std::string>& columnNames = {});

private:
    std::string rowToCSV(const cPoint& row);
    std::string columnNamesToCSV(const std::vector<std::string>& columnNames);
};

#endif // DATAEXPORTER_H
