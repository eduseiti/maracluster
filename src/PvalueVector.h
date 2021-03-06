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
 
#ifndef MARACLUSTER_PVALUEVECTOR_H_
#define MARACLUSTER_PVALUEVECTOR_H_

#include "PvalueCalculator.h"
#include "ScanId.h"

namespace maracluster {

struct PvalueVector {
  double precMass, precMz, retentionTime;
  double polyfit[PvalueCalculator::kPolyfitDegree + 1];
  short peakBins[PvalueCalculator::kMaxScoringPeaks];
  short peakScores[PvalueCalculator::kMaxScoringPeaks];
  int charge, queryCharge;
  ScanId scannr;
};

} /* namespace maracluster */

#endif /* MARACLUSTER_BATCHPVALUEVECTOR_H_ */
