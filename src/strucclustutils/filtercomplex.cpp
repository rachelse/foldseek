#include "DBWriter.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "MemoryMapped.h"
#include "Coordinate16.h"
#include "tmalign/basic_fun.h"
#include "createcomplexreport.h"
#include "LDDT.h"
#include "CalcProbTP.h"
#include <map>
#include <numeric>

#ifdef OPENMP
#include <omp.h>
#endif

struct ComplexFilterCriteria {
    ComplexFilterCriteria() {}
    ComplexFilterCriteria(unsigned int dbKey, unsigned int qTotalAlnLen, unsigned int tTotalAlnLen, double qTM, double tTM, double avglddt, double qChainTm, double tChainTm) :
                        dbKey(dbKey), qTotalAlnLen(qTotalAlnLen), tTotalAlnLen(tTotalAlnLen), qTM(qTM), tTM(tTM) {
                            alignedQChainTmScores.push_back(qChainTm);
                            alignedTChainTmScores.push_back(tChainTm);
                            chainlDDTs.push_back(avglddt);
                        }
    ~ComplexFilterCriteria() {
        alignedQChainTmScores.clear();
        alignedTChainTmScores.clear();
        chainlDDTs.clear();
    }

    bool hasChainnum(bool sameChainNum, int qChainNum, int tChainNum){
        switch (sameChainNum){
            case 1:
                if (qChainNum != tChainNum){
                    return false;
                }else{return true;}
            case 0:
                return true;
}
    }
    bool hasTM(float TMThr, int covMode){
        switch (covMode) {
            case Parameters::COV_MODE_BIDIRECTIONAL:
                return ((qTM>= TMThr) && (tTM >= TMThr));
            case Parameters::COV_MODE_TARGET:
                return (tTM >= TMThr);
            case Parameters::COV_MODE_QUERY:
                return (qTM >= TMThr);
            case Parameters::COV_MODE_LENGTH_QUERY :
            case Parameters::COV_MODE_LENGTH_TARGET :
            case Parameters::COV_MODE_LENGTH_SHORTER :
                return true;
            default:
                return true;
        }
    }
    // bool hasChainTm(float chainTMThr, int covMode, unsigned int qChainNum, unsigned int tChainNum) {
    //     if (chainTMThr > 0 ){
    //         switch (covMode) {
    //             case Parameters::COV_MODE_BIDIRECTIONAL:
    //                 if (alignedQChainTmScores.size()<std::min(qChainNum, tChainNum)){
    //                     return false;
    //                 }
    //                 for (size_t i = 0; i < alignedQChainTmScores.size(); i++) {
    //                     if (alignedQChainTmScores[i] < chainTMThr || alignedTChainTmScores[i] < chainTMThr) {
    //                         return false;
    //                     }
    //                 }
    //                 break;
    //             case Parameters::COV_MODE_TARGET:
    //                 if (alignedQChainTmScores.size()<std::min(qChainNum, tChainNum)){
    //                     return false;
    //                 }
    //                 for (size_t i = 0; i < alignedTChainTmScores.size(); i++) {
    //                     if (alignedTChainTmScores[i] < chainTMThr) {
    //                         return false;
    //                     }
    //                 }
    //                 break;
    //             case Parameters::COV_MODE_QUERY:
    //                 if (alignedQChainTmScores.size()<std::min(qChainNum, tChainNum)){
    //                     return false;
    //                 }
    //                 for (size_t i = 0; i < alignedQChainTmScores.size(); i++) {
    //                     if (alignedQChainTmScores[i] < chainTMThr) {
    //                         return false;
    //                     }
    //                 }
    //                 break;
    //             case Parameters::COV_MODE_LENGTH_QUERY :
    //             case Parameters::COV_MODE_LENGTH_TARGET :
    //             case Parameters::COV_MODE_LENGTH_SHORTER :
    //                 break;
    //         }
    //     }
    //     return true;
    // }
    
    bool hasChainlDDT(float chainlDDTthr, unsigned int qChainNum, unsigned int tChainNum ) {
        if (chainlDDTs.size()<std::min(qChainNum, tChainNum)){
            return false;
        }
        for (size_t i = 0; i < chainlDDTs.size(); i++) {
            if (chainlDDTs[i] < chainlDDTthr) {
                return false;
            }
        }
        return true;
    }

    bool satisfy(int covMode, float covThr, float TMThr, float chainTMThr, bool sameChainNum, int qChainNum, int tChainNum, float chainlDDTthr) {
        const bool covOK = Util::hasCoverage(covThr, covMode, qCov, tCov);
        const bool TMOK = hasTM(TMThr, covMode);
        // const bool chainTMOK = hasChainTm(chainTMThr, covMode, qChainNum, tChainNum);
        const bool chainlDDTOK = hasChainlDDT(chainlDDTthr,  qChainNum, tChainNum);
        const bool numOK = hasChainnum(sameChainNum, qChainNum, tChainNum);
        // return (covOK && TMOK && chainTMOK && numOK);
        return (covOK && TMOK && chainlDDTOK && numOK);
    }

    void update(unsigned int qTotalAlnLen, unsigned int tTotalAlnLen, double avglddt, double qChainTm, double tChainTm) {
        this->qTotalAlnLen += qTotalAlnLen;
        this->tTotalAlnLen += tTotalAlnLen;
        this->alignedQChainTmScores.push_back(qChainTm);
        this->alignedTChainTmScores.push_back(tChainTm);
        chainlDDTs.push_back(avglddt);
    }

    void calcCov(unsigned int qLen, unsigned int tLen) {
        qCov = static_cast<float>(qTotalAlnLen) / static_cast<float>(qLen);
        tCov = static_cast<float>(tTotalAlnLen) / static_cast<float>(tLen);
    }

    unsigned int dbKey;
    unsigned int qTotalAlnLen;
    unsigned int tTotalAlnLen;
    float qCov;
    float tCov;
    double qTM;
    double tTM;

    std::vector<double> alignedQChainTmScores;
    std::vector<double> alignedTChainTmScores;
    std::vector<double> chainlDDTs;
};

unsigned int adjustAlnLen(unsigned int qcov, unsigned int tcov, int covMode) {
    switch (covMode) {
        case Parameters::COV_MODE_BIDIRECTIONAL:
            return (qcov+tcov)/2;
        case Parameters::COV_MODE_TARGET:
            return qcov;
        case Parameters::COV_MODE_QUERY:
            return tcov;
        case Parameters::COV_MODE_LENGTH_QUERY :
        case Parameters::COV_MODE_LENGTH_TARGET :
        case Parameters::COV_MODE_LENGTH_SHORTER :
            return 0;
        default:
            return 0;
    }
}

void fillUArr(const std::string &uString, float (&u)[3][3]) {
    std::string tmp;
    int i = 0;
    int j=0;
    const int ulen = static_cast<int>(uString.size());
    for (int k=0; k < ulen; k++) {
        if (k==ulen-1) {
            u[i][j] = std::stof(tmp);
        } else if (uString[k] == ',') {
            u[i][j] = std::stof(tmp);
            tmp.clear();
            j++;
        } else {
            tmp.push_back(uString[k]);
        }
        if (j == 3) {
            i++;
            j = 0;
        }
    }
}

void fillTArr(const std::string &tString, float (&t)[3]) {
    std::string tmp;
    int i = 0;
    const int tlen = static_cast<int>(tString.size());
    for (int k=0; k<tlen; k++) {
        if (k ==tlen-1) {
            t[i] = std::stof(tmp);
        } else if (tString[k] == ',') {
            t[i] = std::stof(tmp);
            tmp.clear();
            i++;
        } else {
            tmp.push_back(tString[k]);
        }
    }
}

unsigned int cigarToAlignedLength(const std::string &cigar){
    std::string backtrace = Matcher::uncompressAlignment(cigar);
    unsigned int alni = 0;
    for (size_t btPos = 0; btPos < backtrace.size(); btPos++) {
        if (backtrace[btPos] == 'M') {
            alni++;
        }
    }
    return alni;
}

void fillMatchedCoord(float * qdata, float * tdata, 
                Coordinates &qm, Coordinates &tm,
                const std::string &cigar, int qStartPos, int tStartPos, int qLen, int tLen) {
    int qi = qStartPos;
    int ti = tStartPos;
    int mi = 0;

    std::string backtrace = Matcher::uncompressAlignment(cigar);
    for (size_t btPos = 0; btPos < backtrace.size(); btPos++) {
        if (backtrace[btPos] == 'M') {
            qm.x[mi] = qdata[qi];
            qm.y[mi] = qdata[qLen + qi];
            qm.z[mi] = qdata[2*qLen + qi];
            tm.x[mi] = tdata[ti];
            tm.y[mi] = tdata[tLen + ti];
            tm.z[mi] = tdata[2*tLen + ti];
            qi++;
            ti++;
            mi++;
        }
        else if (backtrace[btPos] == 'I') {
            qi++;
        }
        else {
            ti++;
        }
    }
}

double computeChainTmScore(Coordinates &qm, Coordinates &tm, float t[3], float u[3][3], unsigned int alnLen, int tLen) {
    double tmscore = 0;
    float d0 = 1.24*(cbrt(tLen-15)) -1.8;
    float d02 = d0*d0;
    Coordinates tmt(alnLen);
    BasicFunction::do_rotation(tm, tmt, alnLen, t, u);
    for (unsigned int k=0; k<alnLen; k++) {
        double di = BasicFunction::dist(qm.x[k], qm.y[k], qm.z[k], tmt.x[k], tmt.y[k], tmt.z[k]);
        tmscore += 1/(1+di/d02);
    }
    return tmscore;
}

unsigned int getComplexResidueLength( IndexReader *Dbr, std::vector<unsigned int> &ChainKeys) {
        unsigned int ResidueLen = 0;
        for (auto ChainKey: ChainKeys) {
            size_t id = Dbr->sequenceReader->getId(ChainKey);
            // Not accessible
            if (id == NOT_AVAILABLE_CHAIN_KEY)
                return 0;
            ResidueLen += Dbr->sequenceReader->getSeqLen(id);
        }
        return ResidueLen;
}

static void getlookupInfo(
        const std::string &file,
        std::map<unsigned int, std::string> &complexIdtoName,
        std::map<unsigned int, unsigned int> &chainKeyToComplexIdLookup,
        std::map<unsigned int, std::vector<unsigned int>> &complexIdToChainKeysLookup,
        std::vector<unsigned int> &complexIdVec
) {
    if (file.length() == 0) {
        return;
    }
    MemoryMapped lookupDB(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char *data = (char *) lookupDB.getData();
    char *end = data + lookupDB.mappedSize();
    const char *entry[255];
    int prevComplexId =  -1;
    while (data < end && *data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 3) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        auto chainKey = Util::fast_atoi<int>(entry[0]);
        std::string chainName(entry[1], (entry[2] - entry[1]) - 1);
        auto complexId = Util::fast_atoi<int>(entry[2]);
        chainKeyToComplexIdLookup.emplace(chainKey, complexId);
        
        size_t lastUnderscoreIndex = chainName.find_last_of('_');
        std::string complexName = chainName.substr(0, lastUnderscoreIndex);

        if (complexId != prevComplexId) {
            complexIdToChainKeysLookup.emplace(complexId, std::vector<unsigned int>());
            complexIdVec.emplace_back(complexId);
            complexIdtoName.emplace(complexId, complexName);
            prevComplexId = complexId;
        }
        complexIdToChainKeysLookup.at(complexId).emplace_back(chainKey);
        data = Util::skipLine(data);
    }
    lookupDB.close();
}
double calculateLDDT(Coordinates &queryCoords, Coordinates &targetCoords, unsigned int alnLen) {
    vector<double> lddt_scores(alnLen, 0.0);
    int count = 0;
    for (size_t i = 0; i < alnLen; ++i) {
        count = 0;
        for (size_t j = 0; j < alnLen; ++j) {
            if (i != j) {
                double queryDistance = BasicFunction::dist(queryCoords.x[i],queryCoords.y[i],queryCoords.z[i],queryCoords.x[j],queryCoords.y[j],queryCoords.z[j]);
                double targetDistance = BasicFunction::dist(targetCoords.x[i],targetCoords.y[i],targetCoords.z[i],targetCoords.x[j],targetCoords.y[j],targetCoords.z[j]);
                
                if (abs(queryDistance - targetDistance) <= 0.5) {
                    count++;
                }
                if (abs(queryDistance - targetDistance) <= 1.0) {
                    count++;
                }
                if (abs(queryDistance - targetDistance) <= 2.0) {
                    count++;
                }
                if (abs(queryDistance - targetDistance) <= 4.0) {
                    count++;
                }
            }
        }
        lddt_scores[i] = static_cast<double>(count) / ( 4 * (alnLen - 1) );
    }

    double lddt = accumulate(lddt_scores.begin(), lddt_scores.end(), 0.0) / alnLen;
    return lddt;
}

void getAverageCoord( float * qdata, float qavgCoords[3], int qlen, float * tdata, float tavgCoords[3], int tlen, float t[3], float u[3][3]){
    Coordinates tCoords(tlen);
    
    for (int i=0; i < tlen ; i++ ){
        tCoords.x[i] = tdata[i];
        tCoords.y[i] = tdata[tlen+i];
        tCoords.z[i] = tdata[2*tlen+i];
    }
    Coordinates Coordsrot(tlen);
    BasicFunction::do_rotation(tCoords, Coordsrot, tlen, t, u);
    float qxsum=0.0, qysum=0.0, qzsum=0.0;
    for (int i=0; i< qlen;i++){
        qxsum += qdata[i];
        qysum += qdata[qlen+i];
        qzsum += qdata[2*qlen+i];
    }
    qxsum /= qlen;
    qysum /= qlen;
    qzsum /= qlen;
    qavgCoords[0] = qxsum;
    qavgCoords[1] = qysum;
    qavgCoords[2] = qzsum;

    // qavgCoords[0] = qdata[qlen/3];
    // qavgCoords[1] = qdata[qlen + qlen/3];
    // qavgCoords[2] = qdata[qlen*2 + qlen/3];
    float txsum=0.0, tysum=0.0, tzsum=0.0;
    for (int i=0; i< qlen;i++){
        txsum += Coordsrot.x[i];
        tysum += Coordsrot.y[i];
        tzsum += Coordsrot.z[i];
    }
    txsum /= tlen;
    tysum /= tlen;
    tzsum /= tlen;
    tavgCoords[0] = txsum;
    tavgCoords[1] = tysum;
    tavgCoords[2] = tzsum;
    // tavgCoords[0] = Coordsrot.x[tlen/3];
    // tavgCoords[1] = Coordsrot.y[tlen/3];
    // tavgCoords[2] = Coordsrot.z[tlen/3];
}

int filtercomplex(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    int dbaccessMode = (DBReader<unsigned int>::USE_INDEX);

    //TOOD
    //maybe set LDDT threshold based on TMthreshold
    //or use LDDTthr 

    IndexReader* qDbr = NULL;
    qDbr = new IndexReader(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    // DBReader<unsigned int> qStructDbr((par.db1 + "_ca").c_str(), (par.db1 + "_ca.index").c_str(), 
    //                             par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    // qStructDbr.open(DBReader<unsigned int>::NOSORT);
    IndexReader *tStructDbr = NULL;
    IndexReader *qStructDbr = NULL;
    IndexReader* tDbr = NULL;
    // DBReader<unsigned int> *tStructDbr = NULL;
    qStructDbr = new IndexReader(
        par.db1,
            par.threads,
            IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
            touch ? IndexReader::PRELOAD_INDEX : 0,
            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
            "_ca"
    );
    if (sameDB) {
        tDbr = qDbr;
        tStructDbr = qStructDbr;
    }
    else{
        tDbr = new IndexReader(par.db2, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
        tStructDbr = new IndexReader(
                par.db2,
                par.threads,
                IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY_DB1),
                touch ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                "_ca");
    }
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX| DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    size_t localThreads = 1;

    // Debug(Debug::WARNING) << "Monomer will be treated as singleton\nMonomer chain key: \n";
#ifdef OPENMP
localThreads = std::max(std::min((size_t)par.threads, alnDbr.getSize()), (size_t)1);
#endif
    const bool shouldCompress = (par.compressed == true);
    const int db4Type = Parameters::DBTYPE_CLUSTER_RES;
    
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), 1, shouldCompress, db4Type);
    resultWriter.open();

    const int db5Type = Parameters::DBTYPE_GENERIC_DB;
    DBWriter resultWrite5(par.db5.c_str(), par.db5Index.c_str(), 1, shouldCompress, db5Type);
    resultWrite5.open();

    std::string qLookupFile = par.db1 + ".lookup";
    std::string tLookupFile = par.db2 + ".lookup";
    
    chainKeyToComplexId_t qChainKeyToComplexIdMap, tChainKeyToComplexIdMap;
    complexIdToChainKeys_t qComplexIdToChainKeyMap, tComplexIdToChainKeyMap;
    std::map<unsigned int, std::string> qcomplexIdToName, tcomplexIdToName;
    std::vector<unsigned int> qComplexIdVec, tComplexIdVec;
    getlookupInfo(qLookupFile, qcomplexIdToName,qChainKeyToComplexIdMap, qComplexIdToChainKeyMap, qComplexIdVec);
    getlookupInfo(tLookupFile, tcomplexIdToName, tChainKeyToComplexIdMap, tComplexIdToChainKeyMap, tComplexIdVec);
    qChainKeyToComplexIdMap.clear();
    Debug::Progress progress(qComplexIdVec.size());
    std::map<unsigned int, unsigned int> qComplexLength, tComplexLength;
    std::map<unsigned int, std::string> qComplexIdResult;

    for (size_t tComplexIdx = 0; tComplexIdx < tComplexIdVec.size(); tComplexIdx++) {
        unsigned int tComplexId = tComplexIdVec[tComplexIdx];
        std::vector<unsigned int> &tChainKeys = tComplexIdToChainKeyMap.at(tComplexId);
        if (tChainKeys.empty()) {
            continue;
        }
        unsigned int reslen = getComplexResidueLength(tDbr, tChainKeys);
        tComplexLength[tComplexId] =reslen;
    }
    for (size_t qComplexIdx = 0; qComplexIdx < qComplexIdVec.size(); qComplexIdx++) {
        unsigned int qComplexId = qComplexIdVec[qComplexIdx];
        std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeyMap.at(qComplexId);
        if (qChainKeys.empty()) {
            continue;
        }
        unsigned int reslen = getComplexResidueLength(qDbr, qChainKeys);
        qComplexLength[qComplexId] = reslen;
    }

    
    
#pragma omp parallel num_threads(localThreads) 
    {   
        resultToWrite_t result5;
        char buffer[32];
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::string result;
        std::map<unsigned int, std::string> tmpDBKEYut;
        std::map<unsigned int, ComplexFilterCriteria> localComplexMap;
        std::vector<unsigned int> assIdsToDelete;
        std::map<unsigned int, std::vector<unsigned int>> cmplIdToBestAssId; // cmplId : [assId, alnSum]
        std::vector<unsigned int> selectedAssIDs;
        Coordinate16 qcoords;
        Coordinate16 tcoords;
        
        Matcher::result_t res;
#pragma omp for schedule(dynamic, 1)    
        for (size_t queryComplexIdx = 0; queryComplexIdx < qComplexIdVec.size(); queryComplexIdx++) {
            progress.updateProgress();
            unsigned int qComplexId = qComplexIdVec[queryComplexIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeyMap.at(qComplexId);

            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++ ) {
                unsigned int qChainKey = qChainKeys[qChainIdx];
                unsigned int qChainAlnId = alnDbr.getId(qChainKey);
                unsigned int qChainDbId = qDbr->sequenceReader->getId(qChainKey);
                //handling monomer as singleton
                if (qChainAlnId == NOT_AVAILABLE_CHAIN_KEY){
                    char *outpos = Itoa::u32toa_sse2(qComplexId, buffer);
                    result.append(buffer, (outpos - buffer - 1));
                    result.push_back('\n');
                    result5.append(qcomplexIdToName.at(qComplexId) + "\t" + tcomplexIdToName.at(qComplexId) + "\t1.000000\t1.000000\t1.000000\t1.000000\n");
                    break;
                }
                char *data = alnDbr.getData(qChainAlnId, thread_idx);
                while (*data != '\0' ) {
                    ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                    if (!retComplex.isValid){
                        Debug(Debug::ERROR) << "No scorecomplex result provided";
                        EXIT(EXIT_FAILURE);
                    }

                    data = Util::skipLine(data);
                    unsigned int assId = retComplex.assId;
                    unsigned int tChainKey= res.dbKey;
                    unsigned int tChainAlnId = alnDbr.getId(tChainKey);
                    unsigned int tChainDbId = tDbr->sequenceReader->getId(tChainKey);
                    //if target is monomer, break to be singleton
                    if (tChainAlnId == NOT_AVAILABLE_CHAIN_KEY){
                        break;
                    }
                    float u[3][3];
                    float t[3];
                    fillUArr(retComplex.uString, u);
                    fillTArr(retComplex.tString, t);
                    tmpDBKEYut[assId]=retComplex.uString+","+retComplex.tString;

                    char *qcadata = qStructDbr->sequenceReader->getData(qChainDbId, thread_idx);
                    size_t qCaLength = qStructDbr->sequenceReader->getEntryLen(qChainDbId);
                    float* qdata = qcoords.read(qcadata, res.qLen, qCaLength);

                    char *tcadata = tStructDbr->sequenceReader->getData(tChainDbId, thread_idx);
                    size_t tCaLength = tStructDbr->sequenceReader->getEntryLen(tChainDbId);
                    float* tdata = tcoords.read(tcadata, res.dbLen, tCaLength);

                    unsigned int normlen = std::min(res.qLen, res.dbLen);
                    unsigned int alnLen = cigarToAlignedLength(res.backtrace);
                    Coordinates qm(alnLen), tm(alnLen);
                    fillMatchedCoord(qdata, tdata, qm, tm, res.backtrace, res.qStartPos, res.dbStartPos, res.qLen, res.dbLen);
                    double chainTm = computeChainTmScore(qm, tm, t, u, alnLen, res.dbLen);
                    double qChainTm = chainTm / res.qLen;
                    double tChainTm = chainTm/ res.dbLen;

                    float qAvgCoord[3], tAvgCoord[3];
                    getAverageCoord(qdata, qAvgCoord, res.qLen, tdata, tAvgCoord, res.dbLen, t, u) ;
                    // LDDTCalculator *lddtcalculator = NULL;
                    // lddtcalculator = new LDDTCalculator(qStructDbr.getMaxSeqLen() + 1,  tStructDbr->getMaxSeqLen() + 1);
                    // LDDTCalculator::LDDTScoreResult lddtres;
                    // lddtcalculator->initQuery(qCaLength, qm.x, qm.y, qm.z);
                    // lddtres = lddtcalculator->computeLDDTScore(res.dbLen, res.qStartPos, res.dbStartPos, res.backtrace, tm.x, tm.y, tm.z);
                    float lDDTScore = calculateLDDT (qm, tm, alnLen);   
  
                    unsigned int qtotalaln = (std::max(res.qStartPos, res.qEndPos) - std::min(res.qStartPos, res.qEndPos) + 1);
                    unsigned int ttotalaln = (std::max(res.dbStartPos, res.dbEndPos) - std::min(res.dbStartPos, res.dbEndPos) + 1);

                    if (localComplexMap.find(assId) == localComplexMap.end()) {
                        // ComplexFilterCriteria cmplfiltcrit = ComplexFilterCriteria(tChainKey, qtotalaln, ttotalaln, retComplex.qTmScore, retComplex.tTmScore, lddtres.avgLddtScore );
                        ComplexFilterCriteria cmplfiltcrit = ComplexFilterCriteria(tChainKey, qtotalaln, ttotalaln, retComplex.qTmScore, retComplex.tTmScore, lDDTScore, qChainTm, tChainTm);
                        localComplexMap[assId] = cmplfiltcrit;
                    } else {
                        // localComplexMap.at(assId).update(qtotalaln, ttotalaln, lddtres.avgLddtScore);
                        localComplexMap.at(assId).update(qtotalaln, ttotalaln, lDDTScore, qChainTm, tChainTm);
                    }
                    tcadata=nullptr;
                    qcadata=nullptr;
                    tdata=nullptr;
                    qdata=nullptr;
    
                } // while end
            }
            for (auto& assId_res : localComplexMap){
                unsigned int tComplexId = tChainKeyToComplexIdMap.at(assId_res.second.dbKey);
                assId_res.second.calcCov(qComplexLength.at(qComplexId), tComplexLength.at(tComplexId));
                std::vector<unsigned int> tChainKeys = tComplexIdToChainKeyMap.at(tComplexId);
                // if (!assId_res.second.satisfy(par.covMode, par.covThr, par.filtComplexTmThr, par.filtChainTmThr, par.sameChainNumber, qChainKeys.size(), tChainKeys.size())){
                if (!assId_res.second.satisfy(par.covMode, par.covThr, par.filtComplexTmThr, par.filtChainTmThr, par.sameChainNumber, qChainKeys.size(), tChainKeys.size(), par.filtChainLddtThr)){
                    assIdsToDelete.push_back(assId_res.first);
                }
            }

            for (const auto& key : assIdsToDelete) {
                localComplexMap.erase(key);
            }
            
            for (const auto& assId_res : localComplexMap){
                unsigned int tComplexId = tChainKeyToComplexIdMap.at(assId_res.second.dbKey);
                unsigned int alnlen = adjustAlnLen(assId_res.second.qTotalAlnLen, assId_res.second.tTotalAlnLen, par.covMode);
                if (cmplIdToBestAssId.find(tComplexId) == cmplIdToBestAssId.end()){
                    cmplIdToBestAssId[tComplexId] = {assId_res.first, alnlen};
                }
                else {
                    if (alnlen > cmplIdToBestAssId.at(tComplexId)[1]){
                        cmplIdToBestAssId[tComplexId] = {assId_res.first, alnlen};
                    }
                }
            }

            for (const auto& pair : cmplIdToBestAssId){
                selectedAssIDs.push_back(pair.second[0]);
            }

            for (unsigned int assIdidx = 0; assIdidx < selectedAssIDs.size(); assIdidx++){
                unsigned int assId = selectedAssIDs[assIdidx];
                unsigned int tComplexId = tChainKeyToComplexIdMap.at(localComplexMap.at(assId).dbKey);
                char *outpos = Itoa::u32toa_sse2(tComplexId, buffer);
                result.append(buffer, (outpos - buffer - 1));
                result.push_back('\n');
                result5.append(qcomplexIdToName.at(qComplexId) + "\t" + tcomplexIdToName.at(tComplexId) + "\t" + std::to_string(localComplexMap.at(assId).qCov) + "\t" + std::to_string(localComplexMap.at(assId).tCov) + "\t"+ std::to_string(localComplexMap.at(assId).qTM)+"\t"+ std::to_string(localComplexMap.at(assId).tTM)+ "\n");
            }
            #pragma omp critical
            {
                qComplexIdResult[qComplexId]= result;
            }
            result.clear();
            localComplexMap.clear();
            tmpDBKEYut.clear();
            assIdsToDelete.clear();
            cmplIdToBestAssId.clear();
            selectedAssIDs.clear();
        } // for end
        #pragma omp critical
        {
            resultWrite5.writeData(result5.c_str(), result5.length(), 0);
        }
        result5.clear();
    } // MP end
    for (auto &pair : qComplexIdResult){
        resultWriter.writeData(pair.second.c_str(), pair.second.length(), pair.first);
    }
    
    resultWriter.close(true);
    resultWrite5.close(par.dbOut == false);
    alnDbr.close();
    delete qStructDbr;
    delete qDbr;
    
    if (sameDB == false) {
        delete tDbr;
        delete tStructDbr;
    }
    qChainKeyToComplexIdMap.clear();
    tChainKeyToComplexIdMap.clear();
    qComplexIdToChainKeyMap.clear();
    tComplexIdToChainKeyMap.clear();
    qcomplexIdToName.clear();
    tcomplexIdToName.clear();
    qComplexIdVec.clear();
    tComplexIdVec.clear();
    qComplexLength.clear();
    tComplexLength.clear();
    
    return EXIT_SUCCESS;
}
