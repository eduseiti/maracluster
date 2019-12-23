#ifndef MARACLUSTER_MSEMBEDDINGS_H_
#define MARACLUSTER_MSEMBEDDINGS_H_

#include <iostream>
#include <fstream>
#include <vector>

#include "MyException.h"
#include "PvalueVectors.h"

namespace maracluster {

class MSEmbeddings {

    static:
        const unsigned int EMBEDDINGS_DIMENSIONS = 30;

    public:
        MSEmbeddings(char *embeddingsFilename);
        ~MSEmbeddings();

        double calculateCosineDistance(PvalueVectorsDbRow& pvecRow,
                                       PvalueVectorsDbRow& queryPvecRow);

    protected:
        float *allEmbeddings = NULL;
        
}

}
#endif