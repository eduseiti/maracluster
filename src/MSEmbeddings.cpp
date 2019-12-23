#include "MSEmbeddings.h"

namespace maracluster {

MSEmbeddings::MSEmbeddings(char *embeddingsFilename) {

    std::cerr << "Loading embeddings from file " << embeddingsFilename << std::endl;

    ifstream inputFile(embeddingsFilename, ios::in | ios::binary);

    inputFile.seekg(0, ios::end);
    size_t fileSize = inputFile.tellg();
    inputFile.close();

    uint32_t numOfEmbeddedSpectra = fileSize / (sizeof(float) * MSEmbeddings.EMBEDDINGS_DIMENSIONS;

    std::cerr << "Number of embedded spectra: " << numOfEmbeddedSpectra << std::endl;

    allEmbeddings = new float[numOfEmbeddedSpectra];

    inputFile.seekg(0, ios::beg);

    inputFile.read((char *) allEmbeddings, sizeof(allEmbeddings));

    inputFile.close();
}

MSEmbeddings::~MSEmbeddings() {
    if (allEmbeddings) {
        delete allEmbeddings;

        allEmbeddings = NULL;
    }
}

}