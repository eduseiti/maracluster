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
        // static const unsigned int EMBEDDINGS_DIMENSIONS = 30;
        static const unsigned int EMBEDDINGS_DIMENSIONS = 40;

        static constexpr const double EPSILON = 0.00000001;

        static constexpr const char *EMBEDDINGS_FILE_EXTENSION = ".bin";
        static constexpr const char *EMBEDDINGS_ORIGINAL_FILE_EXTENSION = ".txt";
        static constexpr const char *EMBEDDINGS_COMPARISONS_FILE_EXTENSION = ".comparisons";

        // Describe each spectrum embeddings.

        typedef struct _embeddings {
            float embedding[EMBEDDINGS_DIMENSIONS];
        } embeddings;


        // Describe each spectra embeddings comparison

        typedef struct _comparison_data {
            unsigned char spectrum_1_file_index;
            unsigned int  spectrum_1_scannr;
            unsigned char spectrum_2_file_index;
            unsigned int  spectrum_2_scannr;
            double        cosine_similarity;
        } comparison_data;


        MSEmbeddings() : allEmbeddings_(NULL), comparisonsCount_(0u), saveComparisons_(false) {}

        void readEmbeddings(std::string& embeddingsFilename, bool saveComparisons);

        double calculateCosineDistance(PvalueVectorsDbRow& pvecRow,
                                       PvalueVectorsDbRow& queryPvecRow);

        void finalizeEmbeddings();

    protected:
        // Points to the entire embeddings array.
        // Spectrum order here MUST follow the same order of the "File with spectrum files to be processed in batch", either 
        // the spectrum files as the spectrum order within each file.

        embeddings *allEmbeddings_ = NULL;

        
        // Maps each original spectra filename (e.g. .mgf files) into the entire embeddings array, properly indicating the first
        // embeddings for the first spectrum within that file.

        std::map<std::string, embeddings*> spectraFileEmbeddingsMap_;


        // Array of original spectra filenames, in the order read from the embeddings files.

        std::vector<std::string> spectraFilenames_;


        // File handle for storing the embeddings comparisons

        std::ofstream outputFile_;

        unsigned int comparisonsCount_;

        bool saveComparisons_;
};

extern MSEmbeddings embeddedSpectra_;

}
#endif