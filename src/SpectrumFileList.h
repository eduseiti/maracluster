/******************************************************************************  
  Copyright 2015 Matthew The <matthew.the@scilifelab.se>
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
 ******************************************************************************/
 
#ifndef MARACLUSTER_SPECTRUMFILELIST_H_
#define MARACLUSTER_SPECTRUMFILELIST_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>

#include "Globals.h"
#include "MyException.h"
#include "ScanId.h"

namespace maracluster {

class SpectrumFileList {
#ifdef USE_EMBEDDINGS
  protected:
    /* A counter for the total number of spectra inside a given file; used to associated a spectrum index */
    std::vector<unsigned int> fileIndexScansCount_;

    uint32_t totalScansCount_;
#endif    

  public:
#ifdef USE_EMBEDDINGS  
    SpectrumFileList() : totalScansCount_(0) { }
#else
    SpectrumFileList() { }
#endif    
    
    ScanId getScanId(const std::string& filePath, const unsigned int scannr);
    
    unsigned int getMergeOffset() const { return size(); }
    
    inline unsigned int getScannr(const ScanId& scanId) const { 
      return scanId.scannr; 
    }

    inline std::string getFilePath(const ScanId& scanId) const {
      return getFilePath(scanId.fileIdx); 
    }
    
    inline std::string getFilePath(const unsigned int fileIdx) const {
      if (fileIdx < fileIndexVector_.size()) {
        return fileIndexVector_[fileIdx]; 
      } else {
        std::stringstream ss;
        ss << "(SpectrumFileList.h) FileIdx " << fileIdx << " out of range" << std::endl;
        throw MyException(ss);
      }
    }
    
    inline size_t size() const {
      return fileIndexVector_.size();
    }
    
    inline unsigned int getFileIdx(const ScanId scanId) const {
      return scanId.fileIdx; 
    }
    
    inline const std::vector<std::string>& getFilePaths() const {
      return fileIndexVector_;
    }
#ifdef USE_EMBEDDINGS
    unsigned int addFileEmbeddings(const std::string& filePath);
#endif
    void addFile(const std::string& filePath);

    void initFromFile(const std::string& fileListFN);
  protected:
    std::map<std::string, unsigned int> fileIndexMap_;
    std::vector<std::string> fileIndexVector_;
};

} /* namespace maracluster */

#endif /* MARACLUSTER_SPECTRUMFILELIST_H_ */
