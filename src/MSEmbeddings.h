#ifndef MARACLUSTER_MSEMBEDDINGS_H_
#define MARACLUSTER_MSEMBEDDINGS_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "MyException.h"
#include "PvalueVectors.h"

namespace maracluster {

class MSEmbeddings {

    public:
        static const unsigned int EMBEDDINGS_DIMENSIONS = 30;

        MSEmbeddings() {}
        ~MSEmbeddings();

        void readEmbeddings(std::string& embeddingsFilename);

        double calculateCosineDistance(PvalueVectorsDbRow& pvecRow,
                                       PvalueVectorsDbRow& queryPvecRow);

    protected:
        float **allEmbeddings = NULL;
        
};

#ifdef USE_EMBEDDINGS
extern MSEmbeddings embeddedSpectra_;
#endif

}
#endif