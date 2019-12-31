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



    if (Globals::VERB > 5) {
        std::cerr << "calculateCosineDistance" << std::endl;

        std::cerr << "1. spectraFileEmbeddingsMap_.size()=" << spectraFileEmbeddingsMap_.size() << std::endl;
        std::cerr << "2. spectraFilenames_.size()=" << spectraFilenames_.size() << std::endl;
    
        std::cerr << "3. pvecRow.scannr.getScanIndex()=" << pvecRow.scannr.getScanIndex() << ", pvecRow.scannr.getFileIndex()=" << pvecRow.scannr.getFileIndex() << std::endl;
        std::cerr << "4. queryPvecRow.scannr.getScanIndex()=" << queryPvecRow.scannr.getScanIndex() << ", queryPvecRow.scannr.getFileIndex()=" << queryPvecRow.scannr.getFileIndex() << std::endl;

        std::cerr << "5. spectraFilenames_[pvecRow.scannr.getFileIndex()]=" << spectraFilenames_[pvecRow.scannr.getFileIndex()] << std::endl;
        std::cerr << "6. spectraFileEmbeddingsMap_[spectraFilenames_[pvecRow.scannr.getFileIndex()]]=" << spectraFileEmbeddingsMap_[spectraFilenames_[pvecRow.scannr.getFileIndex()]] << std::endl;

        std::cerr << "7. spectraFilenames_[queryPvecRow.scannr.getFileIndex()]=" << spectraFilenames_[queryPvecRow.scannr.getFileIndex()] << std::endl;
        std::cerr << "8. spectraFileEmbeddingsMap_[spectraFilenames_[queryPvecRow.scannr.getFileIndex()]]=" << spectraFileEmbeddingsMap_[spectraFilenames_[queryPvecRow.scannr.getFileIndex()]] << std::endl;
    }

    float *embedding_A = spectraFileEmbeddingsMap_[spectraFilenames_[pvecRow.scannr.getFileIndex()]][pvecRow.scannr.getScanIndex()].embedding;
    float *embedding_B = spectraFileEmbeddingsMap_[spectraFilenames_[queryPvecRow.scannr.getFileIndex()]][queryPvecRow.scannr.getScanIndex()].embedding;

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

    //
    // Read the embeddings entire array
    // 

    std::cerr << "Loading embeddings from file " << embeddingsFilename +  MSEmbeddings::EMBEDDINGS_FILE_EXTENSION << std::endl;

    std::ifstream inputFile(embeddingsFilename + MSEmbeddings::EMBEDDINGS_FILE_EXTENSION, std::ios::in | std::ios::binary);

    inputFile.seekg(0, std::ios::end);
    size_t fileSize = inputFile.tellg();

    uint32_t numOfEmbeddedSpectra = fileSize / (sizeof(float) * MSEmbeddings::EMBEDDINGS_DIMENSIONS);

    std::cerr << "Number of embedded spectra: " << numOfEmbeddedSpectra << std::endl;

    allEmbeddings = new embeddings[numOfEmbeddedSpectra];

    inputFile.seekg(0, std::ios::beg);

    inputFile.read((char *) allEmbeddings, fileSize);

    inputFile.close();

    //
    // Now read the embeddings original spectra filename
    //

    std::string line;
    std::string fileBeingHandled;

    inputFile.open(embeddingsFilename + MSEmbeddings::EMBEDDINGS_ORIGINAL_FILE_EXTENSION);

    uint32_t currentEmbeddingsIndex = 0;
    uint32_t currentFileStartingIndex = 0;


    while (getline(inputFile, line)) {
        std::map<std::string, embeddings*>::iterator it = spectraFileEmbeddingsMap_.find(line);

        if (it == spectraFileEmbeddingsMap_.end()) {
            if (currentEmbeddingsIndex > 0) {
                std::cerr << "-- Original spectra filename " << fileBeingHandled << " had total number of spectra of " << currentEmbeddingsIndex - currentFileStartingIndex << std::endl;
            }

            std::cerr << "Mapping filename " << line << " to embedding (absolute) index " << currentEmbeddingsIndex << std::endl;

            spectraFileEmbeddingsMap_[line] = &allEmbeddings[currentEmbeddingsIndex];

            currentFileStartingIndex = currentEmbeddingsIndex;
            fileBeingHandled = line;

            std::cerr << "Filename " << line << " stored as file index " << spectraFilenames_.size() << std::endl;

            spectraFilenames_.push_back(line);
        }

        currentEmbeddingsIndex++;
    }

    std::cerr << "Mapped " << currentEmbeddingsIndex << " spectra to filenames" << std::endl;   
}

}