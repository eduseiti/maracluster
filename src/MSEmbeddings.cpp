#include "MSEmbeddings.h"

namespace maracluster {

MSEmbeddings::~MSEmbeddings() {
    if (allEmbeddings) {
        delete allEmbeddings;

        allEmbeddings = NULL;
    }
}

double MSEmbeddings::calculateCosineDistance(PvalueVectorsDbRow& pvecRow,
                                             PvalueVectorsDbRow& queryPvecRow) {

    float *embedding_A = allEmbeddings[pvecRow.scannr.getAbsoluteScanIndex()];
    float *embedding_B = allEmbeddings[queryPvecRow.scannr.getAbsoluteScanIndex()];

    double multi = 0.0;
    double norm_A = 0.0;
    double norm_B = 0.0;

    for (size_t i = 0; i < MSEmbeddings::EMBEDDINGS_DIMENSIONS; i++) {
        multi += embedding_A[i] * embedding_B[i];
        norm_A += embedding_A[i] * embedding_A[i];
        norm_B += embedding_B[i] * embedding_B[i];
    }

    return -100 * multi / (sqrt(norm_A) * sqrt(norm_B));
}

void MSEmbeddings::readEmbeddings(std::string& embeddingsFilename) {

    std::cerr << "Loading embeddings from file " << embeddingsFilename << std::endl;

    std::ifstream inputFile(embeddingsFilename, std::ios::in | std::ios::binary);

    inputFile.seekg(0, std::ios::end);
    size_t fileSize = inputFile.tellg();
    inputFile.close();

    uint32_t numOfEmbeddedSpectra = fileSize / (sizeof(float) * MSEmbeddings::EMBEDDINGS_DIMENSIONS);

    std::cerr << "Number of embedded spectra: " << numOfEmbeddedSpectra << std::endl;

    allEmbeddings = (float **) new float[numOfEmbeddedSpectra][MSEmbeddings::EMBEDDINGS_DIMENSIONS];

    inputFile.seekg(0, std::ios::beg);

    inputFile.read((char *) allEmbeddings, sizeof(allEmbeddings));

    inputFile.close();
}



}