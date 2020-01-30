/******************************************************************************  

 This module allows calculating the spectra distance using embeddings previouly
 calculated outside MaRaCluster.

 This module is included throught he USE_EMBEDDINGS compilation flag.

 In order to work properly, this module depends on the provided embedded list
 respect the exact order the spectra is read from the MaRaCluster input files, 
 according to the BATCH command: not only the spectra input files ― e.g. .mgf
 or .ms2 files ― have to follow the same order, but also the spectra order
 within each file shall be the same, as currently the embedding load at the
 distance calculation process ― performed by the "calculateCosineDistance()"
 method in this file ― uses the scan loading order (inside 
 "SpectrumFiles::getPeakCountsAndPrecursorMzs()" method) as the embeddings
 access index.

******************************************************************************/

#include "MSEmbeddings.h"

namespace maracluster {

double MSEmbeddings::calculateCosineDistance(PvalueVectorsDbRow& pvecRow,
                                             PvalueVectorsDbRow& queryPvecRow) {

    if (Globals::VERB > 5) {
        std::cerr << "calculateCosineDistance" << std::endl;

        std::cerr << "1. spectraFileEmbeddingsMap_.size()=" << spectraFileEmbeddingsMap_.size() << std::endl;
        std::cerr << "2. spectraFilenames_.size()=" << spectraFilenames_.size() << std::endl;
    
        std::cerr << "3. pvecRow.scannr.getScanIndex()=" << pvecRow.scannr.getScanIndex() << 
                     ", pvecRow.scannr.getFileIndex()=" << pvecRow.scannr.getFileIndex() << std::endl;

        std::cerr << "4. queryPvecRow.scannr.getScanIndex()=" << queryPvecRow.scannr.getScanIndex() << 
                     ", queryPvecRow.scannr.getFileIndex()=" << queryPvecRow.scannr.getFileIndex() << std::endl;

        std::cerr << "5. spectraFilenames_[pvecRow.scannr.getFileIndex()]=" << spectraFilenames_[pvecRow.scannr.getFileIndex()] << std::endl;

        std::cerr << "6. spectraFileEmbeddingsMap_[spectraFilenames_[pvecRow.scannr.getFileIndex()]]=" << spectraFileEmbeddingsMap_[spectraFilenames_[pvecRow.scannr.getFileIndex()]] << std::endl;

        std::cerr << "7. spectraFilenames_[queryPvecRow.scannr.getFileIndex()]=" << spectraFilenames_[queryPvecRow.scannr.getFileIndex()] << std::endl;

        std::cerr << "8. spectraFileEmbeddingsMap_[spectraFilenames_[queryPvecRow.scannr.getFileIndex()]]=" << 
                     spectraFileEmbeddingsMap_[spectraFilenames_[queryPvecRow.scannr.getFileIndex()]] << std::endl;
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

    double cosine_similarity = multi / (sqrt(norm_A) * sqrt(norm_B));

    if (saveComparisons_) {
        comparison_data thisComparison;

        thisComparison.spectrum_1_file_index = pvecRow.scannr.getFileIndex();
        thisComparison.spectrum_1_scannr     = pvecRow.scannr.getScanIndex();
        thisComparison.spectrum_2_file_index = queryPvecRow.scannr.getFileIndex();
        thisComparison.spectrum_2_scannr     = queryPvecRow.scannr.getScanIndex();
        thisComparison.cosine_similarity     = cosine_similarity;

#pragma omp critical (write_comparisons) 
{
        outputFile_.write((char *)&thisComparison, sizeof(comparison_data));
}

        comparisonsCount_++;
    }


    // -100 is the factor to take the cosine similarity until the p-value scale.

    return -100 * cosine_similarity;
}



//
// This method reads the embeddings to be used for the spectrum distance calculation. The embeddings info shall be provided 
// in 2 different files, described bellow; the "embeddingsFilename" parameter captures the filenames root part:
//
// <embeddingsFilename>.bin: Binary file with the each spectrum embeddings for all the N spectrum available.
//
// <embeddingsFilename>.txt: Text file with the original file name from which each particular spectrum was read; contains one
//                           line for each of the N spectrum available.
//

void MSEmbeddings::readEmbeddings(std::string& embeddingsFilename, bool saveComparisons) {

    //
    // Read the embeddings entire array
    // 

    std::cerr << "Loading embeddings from file " << embeddingsFilename +  MSEmbeddings::EMBEDDINGS_FILE_EXTENSION << std::endl;

    std::ifstream inputFile(embeddingsFilename + MSEmbeddings::EMBEDDINGS_FILE_EXTENSION, std::ios::in | std::ios::binary);

    inputFile.seekg(0, std::ios::end);
    size_t fileSize = inputFile.tellg();

    std::cerr << "Filesize: " << fileSize << std::endl;

    uint32_t numOfEmbeddedSpectra = fileSize / (sizeof(float) * MSEmbeddings::EMBEDDINGS_DIMENSIONS);

    std::cerr << "Number of embedded spectra: " << numOfEmbeddedSpectra << std::endl;

    allEmbeddings_ = new embeddings[numOfEmbeddedSpectra];

    inputFile.seekg(0, std::ios::beg);

    inputFile.read((char *) allEmbeddings_, fileSize);

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

            spectraFileEmbeddingsMap_[line] = &allEmbeddings_[currentEmbeddingsIndex];

            currentFileStartingIndex = currentEmbeddingsIndex;
            fileBeingHandled = line;

            std::cerr << "Filename " << line << " stored as file index " << spectraFilenames_.size() << std::endl;

            spectraFilenames_.push_back(line);
        }

        currentEmbeddingsIndex++;
    }

    std::cerr << "-- Original spectra filename " << fileBeingHandled << " had total number of spectra of " << currentEmbeddingsIndex - currentFileStartingIndex << std::endl;
    std::cerr << "Mapped " << currentEmbeddingsIndex << " spectra to filenames" << std::endl;


    //
    // Prepare to save the embeddings comparisons, if the flag has been set
    //

    if (saveComparisons) {
        saveComparisons_ = saveComparisons;

        outputFile_.open(embeddingsFilename + MSEmbeddings::EMBEDDINGS_COMPARISONS_FILE_EXTENSION, std::ios::out | std::ios::binary);
    }
}



void MSEmbeddings::finalizeEmbeddings() {
    outputFile_.close();

    std::cerr << "Stored " << comparisonsCount_ << " comparisons data." << std::endl;
}

}