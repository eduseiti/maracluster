#ifndef MARACLUSTER_MSEMBEDDINGS_H_
#define MARACLUSTER_MSEMBEDDINGS_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "MyException.h"
#include "PvalueVectors.h"

namespace maracluster {

class MSEmbeddings {

    public:
        static const unsigned int EMBEDDINGS_DIMENSIONS = 30;

        static constexpr const char *EMBEDDINGS_FILE_EXTENSION = ".bin";
        static constexpr const char *EMBEDDINGS_ORIGINAL_FILE_EXTENSION = ".txt";

        MSEmbeddings() {}
        ~MSEmbeddings();

        void readEmbeddings(std::string& embeddingsFilename);

        double calculateCosineDistance(PvalueVectorsDbRow& pvecRow,
                                       PvalueVectorsDbRow& queryPvecRow);

    protected:
        float **allEmbeddings = NULL;

        std::map<std::string, float**> spectraFileEmbeddingsMap_;
        std::vector<std::string> spectraFilenames_;
};

extern MSEmbeddings embeddedSpectra_;

}
#endif