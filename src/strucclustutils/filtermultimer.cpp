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
#include "MultimerUtil.h"
#include "MultimerInterface.h"
#include "LDDT.h"
#include <map>
#include <set>

#ifdef OPENMP
#include <omp.h>
#endif

struct Complex {
    int complexId;
    unsigned int nChain;
    unsigned int complexLength;
    std::string complexName;
    std::vector<unsigned int> chainLengths;
    std::vector<unsigned int> chainKeys;

    // Coordinate16 Coords;

    Complex() : complexId(0), nChain(0), complexLength(0), complexName("") {}
    ~Complex() {
        chainKeys.clear();
    }
};

unsigned int adjustAlnLen(unsigned int qcov, unsigned int tcov, int covMode) {
    switch (covMode) {
        case Parameters::COV_MODE_BIDIRECTIONAL:
            return (qcov+tcov)/2;
        case Parameters::COV_MODE_TARGET:
            return qcov;
        case Parameters::COV_MODE_QUERY:
            return tcov;
        default:
            return 0;
    }
}

class ComplexFilterCriteria {
public:
    unsigned int targetComplexId;

    // per complex
    unsigned int qTotalAlnLen;
    unsigned int tTotalAlnLen;
    float qCov;
    float tCov;
    float interfaceLddt;
    double qTm;
    double tTm;
    float t[3];
    float u[3][3];

    // per chain : criteria for chainTmThr & lddtThr
    std::vector<unsigned int> qAlnChainKeys;
    std::vector<unsigned int> tAlnChainKeys;
    std::vector<double> qAlnChainTms;
    std::vector<double> tAlnChainTms;
    std::vector<std::vector<float>> qAlnChains;
    std::vector<std::vector<float>> tAlnChains;
    std::vector<float> talignedx, talignedy, talignedz, qalignedx, qalignedy, qalignedz;
    // std::map<unsigned int, std::vector<std::tuple<float, float, float>>> talignedInterface, qalignedInterface;

    ComplexFilterCriteria() {}
    ComplexFilterCriteria(unsigned int targetComplexId, double qTm, double tTm, float tstring[3], float ustring[3][3]) :
                            targetComplexId(targetComplexId), qTotalAlnLen(0), tTotalAlnLen(0), qTm(qTm), tTm(tTm) {
                                std::copy(tstring, tstring + 3, t);
                                for (int i = 0; i < 3; i++) {
                                    std::copy(ustring[i], ustring[i] + 3, u[i]);
                                }
                                talignedx.reserve(100);
                                talignedy.reserve(100);
                                talignedz.reserve(100);
                                qalignedx.reserve(100);
                                qalignedy.reserve(100);
                                qalignedz.reserve(100);
                            }
    ~ComplexFilterCriteria() {
        qAlnChainTms.clear();
        tAlnChainTms.clear();
        qAlnChainKeys.clear();
        tAlnChainKeys.clear();
        qAlnChains.clear();
        tAlnChains.clear();
    }

    bool hasTm(float TmThr, int covMode){
        switch (covMode) {
            case Parameters::COV_MODE_BIDIRECTIONAL:
                return ((qTm>= TmThr) && (tTm >= TmThr));
            case Parameters::COV_MODE_TARGET:   
                return (tTm >= TmThr);
            case Parameters::COV_MODE_QUERY:
                return (qTm >= TmThr);
            default:
                return true;
        }
    }

    // bool hasChainNum(int covMode, int filterMode, size_t qChainNum, size_t tChainNum ){
    //     switch (filterMode){
    //         case LocalParameters::FILTER_MODE_INTERFACE:
    //             switch (covMode) {
    //                 case Parameters::COV_MODE_BIDIRECTIONAL:
    //                     return (qAlnChainTms.size()==qChainNum && qChainNum==tChainNum);
    //                 case Parameters::COV_MODE_TARGET:
    //                     return (tAlnChainTms.size()==tChainNum);
    //                 case Parameters::COV_MODE_QUERY:
    //                     return (qAlnChainTms.size()==qChainNum);
    //                 default:
    //                     return true;
    //             }
    //         case LocalParameters::FILTER_MODE_CONFORMATION:
    //             switch (covMode) {
    //                 case Parameters::COV_MODE_BIDIRECTIONAL:
    //                     return (qChainNum==tChainNum);
    //                 case Parameters::COV_MODE_TARGET:
    //                     return (qChainNum>=tChainNum);
    //                 case Parameters::COV_MODE_QUERY:
    //                     return (qChainNum<=tChainNum);
    //                 default:
    //                     return true;
    //             }
    //         default:
    //             return true;
    //     }
    // } 

    bool hasChainTm(float chainTmThr, int covMode, unsigned int qChainNum, unsigned int tChainNum) {
        if (qAlnChainTms.size()<std::min(qChainNum, tChainNum)){
            return false;
        }
        switch (covMode) {
            case Parameters::COV_MODE_BIDIRECTIONAL:
                for (size_t i = 0; i < qAlnChainTms.size(); i++) {
                    if (qAlnChainTms[i] < chainTmThr || tAlnChainTms[i] < chainTmThr) {
                        return false;
                    }
                }
                break;
            case Parameters::COV_MODE_TARGET:
                for (size_t i = 0; i < tAlnChainTms.size(); i++) {
                    if (tAlnChainTms[i] < chainTmThr) {
                        return false;
                    }
                }
                break;
            case Parameters::COV_MODE_QUERY:
                for (size_t i = 0; i < qAlnChainTms.size(); i++) {
                    if (qAlnChainTms[i] < chainTmThr) {
                        return false;
                    }
                }
                break;
            default:
                return true;
        }
        return true;
    }

    bool satisfy(int covMode, float covThr, float TmThr, float chainTmThr, float iLddtThr, size_t qChainNum, size_t tChainNum ) {
        const bool covOK = covThr ? Util::hasCoverage(covThr, covMode, qCov, tCov) : true;
        const bool TmOK = TmThr ? hasTm(TmThr, covMode) : true;
        // const bool chainNumOK = hasChainNum(covMode, filterMode, qChainNum, tChainNum);
        const bool chainTmOK = chainTmThr ? hasChainTm(chainTmThr, covMode, qChainNum, tChainNum) : true;
        const bool lddtOK = iLddtThr ? (interfaceLddt >= iLddtThr) : true;
        // const bool conformationOK = isConformation(filterMode, chainTmThr);
        // return (covOK && TmOK && chainNumOK && chainTmOK);

        //TODO: If considering only cov, complextm, lddt -> should check chain number.
         return (covOK && TmOK && chainTmOK && lddtOK);
    }

    void updateAln(unsigned int qAlnLen, unsigned int tAlnLen) {
        qTotalAlnLen += qAlnLen;
        tTotalAlnLen += tAlnLen;
    }

    void updateChainTmAndInterface(unsigned int qChainKey, unsigned int tChainKey, unsigned int alnLen, float *qdata, float *tdata,
                            const std::string &cigar, int qStartPos, int tStartPos, int qLen, int tLen, float t[3], float u[3][3],
                            std::map<unsigned int, std::vector<unsigned int>> & qInterfaceIndex) {
        
        qAlnChainKeys.push_back(qChainKey);
        tAlnChainKeys.push_back(tChainKey);
        int qi = qStartPos;
        int ti = tStartPos;
        int mi = 0;
        std::string backtrace = Matcher::uncompressAlignment(cigar);
        Coordinates qm(alnLen), tm(alnLen);
        std::map<unsigned int, unsigned int> query_to_target;
        for (size_t btPos = 0; btPos < backtrace.size(); btPos++) {
            if (backtrace[btPos] == 'M') {
                qm.x[mi] = qdata[qi];
                qm.y[mi] = qdata[qLen + qi];
                qm.z[mi] = qdata[2*qLen + qi];
                tm.x[mi] = tdata[ti];
                tm.y[mi] = tdata[tLen + ti];
                tm.z[mi] = tdata[2*tLen + ti];
                query_to_target[qi] = mi;
                qi++;
                ti++;
                mi++;
            }
            else if (backtrace[btPos] == 'I') {
                query_to_target[qi] = -1;
                qi++;
            }
            else {
                ti++;
            }
        }
        double tmscore = 0;
        float d0 = 1.24*(cbrt(tLen-15)) -1.8;
        float d02 = d0*d0;
        Coordinates tmt(alnLen);
        BasicFunction::do_rotation(tm, tmt, alnLen, t, u);
        for (unsigned int k=0; k<alnLen; k++) {
            double di = BasicFunction::dist(qm.x[k], qm.y[k], qm.z[k], tmt.x[k], tmt.y[k], tmt.z[k]);
            tmscore += 1/(1+di/d02);
        }
        qAlnChainTms.push_back(tmscore/qLen);
        tAlnChainTms.push_back(tmscore/tLen);
        for (unsigned int qindex : qInterfaceIndex[qChainKey]){
            int aindex = query_to_target[qindex];
            if (query_to_target.find(qindex) != query_to_target.end() && aindex != -1){
                qalignedx.push_back(qm.x[aindex]);
                qalignedy.push_back(qm.y[aindex]);
                qalignedz.push_back(qm.z[aindex]);
                talignedx.push_back(tmt.x[aindex]);
                talignedy.push_back(tmt.y[aindex]);
                talignedz.push_back(tmt.z[aindex]);
            }
        }
    }
    void calcCov(unsigned int qLen, unsigned int tLen) {
        qCov = static_cast<float>(qTotalAlnLen) / static_cast<float>(qLen);
        tCov = static_cast<float>(tTotalAlnLen) / static_cast<float>(tLen);
    }

    void computeInterfaceLddt() {
        unsigned int alnLen = qalignedx.size();
        if (alnLen == 0){
            return;
        }
        std::string bt(alnLen, 'M');
        LDDTCalculator *lddtcalculator = NULL;
        lddtcalculator = new LDDTCalculator(alnLen+1, alnLen+1);
        lddtcalculator->initQuery(alnLen, qalignedx.data(), qalignedy.data(), qalignedz.data());
        LDDTCalculator::LDDTScoreResult lddtres = lddtcalculator->computeLDDTScore(alnLen, 0, 0, bt, talignedx.data(), talignedy.data(), talignedz.data());
        interfaceLddt = lddtres.avgLddtScore;
        delete lddtcalculator;
        qalignedx.clear();
        qalignedx.shrink_to_fit();
        qalignedy.clear();
        qalignedy.shrink_to_fit();
        qalignedz.clear();
        qalignedz.shrink_to_fit();
        talignedx.clear();
        talignedx.shrink_to_fit();
        talignedy.clear();
        talignedy.shrink_to_fit();
        talignedz.clear();
        talignedz.shrink_to_fit();
    }
};

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

void getComplexResidueLength( IndexReader *Dbr, std::vector<Complex> &complexes) {
    for (size_t complexIdx = 0; complexIdx < complexes.size(); complexIdx++) {
        Complex *complex = &complexes[complexIdx];
        std::vector<unsigned int> &chainKeys = complex->chainKeys;
        if (chainKeys.empty()) {
            continue;
        }
        unsigned int cmpllen = 0;
        for (auto chainKey: chainKeys) {
            size_t id = Dbr->sequenceReader->getId(chainKey);
            unsigned int reslen = Dbr->sequenceReader->getSeqLen(id);
            complex->chainLengths.push_back(reslen);
            cmpllen += Dbr->sequenceReader->getSeqLen(id);
        }
        complex->complexLength = cmpllen;
    }
}

static void getlookupInfo(
        const std::string &file,
        std::map<unsigned int, unsigned int> &chainKeyToComplexIdLookup,
        std::vector<Complex> &complexes,
        std::map<unsigned int, unsigned int> &complexIdtoIdx
) {
    if (file.length() == 0) {
        return;
    }
    MemoryMapped lookupDB(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
    char *data = (char *) lookupDB.getData();
    char *end = data + lookupDB.mappedSize();
    const char *entry[255];

    int prevComplexId =  -1;
    int nComplex = 0;
    while (data < end && *data != '\0') {
        const size_t columns = Util::getWordsOfLine(data, entry, 255);
        if (columns < 3) {
            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
            continue;
        }
        auto chainKey = Util::fast_atoi<int>(entry[0]);
        auto complexId = Util::fast_atoi<int>(entry[2]);
        chainKeyToComplexIdLookup.emplace(chainKey, complexId);
        
        std::string chainName(entry[1], (entry[2] - entry[1]) - 1);
        size_t lastUnderscoreIndex = chainName.find_last_of('_');
        std::string complexName = chainName.substr(0, lastUnderscoreIndex);

        if (complexId != prevComplexId) {
            
            Complex complex;
            complex.complexId = complexId;
            complex.complexName = complexName;
            complexIdtoIdx.emplace(complexId, nComplex);
            complexes.emplace_back(complex);

            prevComplexId = complexId;
            nComplex++;
        }
        complexes.back().chainKeys.emplace_back(chainKey);
        complexes.back().nChain++;
        data = Util::skipLine(data);
    }
    lookupDB.close();
}

int filtermultimer(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    int dbaccessMode = (DBReader<unsigned int>::USE_INDEX);

    IndexReader* qDbr = NULL;
    qDbr = new IndexReader(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    DBReader<unsigned int> qStructDbr((par.db1 + "_ca").c_str(), (par.db1 + "_ca.index").c_str(), 
                                par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    qStructDbr.open(DBReader<unsigned int>::NOSORT);

    IndexReader* tDbr = NULL;
    DBReader<unsigned int> *tStructDbr = NULL;
    if (sameDB) {
        tDbr = qDbr;
        tStructDbr = &qStructDbr;
    }
    else{
        tDbr = new IndexReader(par.db2, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
        tStructDbr = new DBReader<unsigned int>((par.db2 + "_ca").c_str(), (par.db2 + "_ca.index").c_str(),
                                           par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
        tStructDbr->open(DBReader<unsigned int>::NOSORT);
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
    
    std::vector<Complex> qComplexes, tComplexes;
    std::map<unsigned int, unsigned int> qComplexIdToIdx, tComplexIdToIdx;
    chainKeyToComplexId_t qChainKeyToComplexIdMap, tChainKeyToComplexIdMap;

    getlookupInfo(qLookupFile, qChainKeyToComplexIdMap, qComplexes, qComplexIdToIdx);
    getComplexResidueLength(qDbr, qComplexes);
    Debug::Progress progress(qComplexes.size());
    std::map<unsigned int, resultToWrite_t> qComplexIdResult;

    if (sameDB) {
        tChainKeyToComplexIdMap = qChainKeyToComplexIdMap;
        tComplexes = qComplexes;
        tComplexIdToIdx = qComplexIdToIdx;
    } else {
        getlookupInfo(tLookupFile, tChainKeyToComplexIdMap, tComplexes, tComplexIdToIdx);
        getComplexResidueLength(tDbr, tComplexes);
    }
    
#pragma omp parallel num_threads(localThreads) 
    {   
        resultToWrite_t result5;
        char buffer[32];
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        resultToWrite_t result;
        std::map<unsigned int, ComplexFilterCriteria> localComplexMap;
        std::map<unsigned int, std::vector<unsigned int>> cmplIdToBestAssId; // cmplId : [assId, alnSum]
        std::vector<unsigned int> selectedAssIDs;
        Coordinate16 qcoords;
        Coordinate16 tcoords;
        
        Matcher::result_t res;
#pragma omp for schedule(dynamic, 1)    
        for (size_t qComplexIdx = 0; qComplexIdx < qComplexes.size(); qComplexIdx++) {
            progress.updateProgress();
   
            Complex qComplex = qComplexes[qComplexIdx];
            unsigned int qComplexId = qComplex.complexId;
            std::vector<unsigned int> qChainKeys = qComplex.chainKeys;
            std::map<unsigned int, std::vector<unsigned int>> qInterfaceIndex; //chainKey: residueIdxVec
            if (qChainKeys.size() > 1){
                for (size_t qChainIdx1 = 0; qChainIdx1 < qChainKeys.size() -1 ; qChainIdx1++ ){
                    unsigned int qChainKey1 = qChainKeys[qChainIdx1];
                    unsigned int qChainDbId1 = qDbr->sequenceReader->getId(qChainKey1);
                    char *qcadata1 = qStructDbr.getData(qChainDbId1, thread_idx);
                    size_t qCaLength1 = qStructDbr.getEntryLen(qChainDbId1);
                    size_t qChainLen1 = qDbr->sequenceReader->getSeqLen(qChainDbId1);   
                    float* qdata1 = qcoords.read(qcadata1, qChainLen1, qCaLength1);
                    Interface *interface = NULL;
                    interface = new Interface(qChainLen1);  
                    interface->initQuery(qdata1, &qdata1[qChainLen1], &qdata1[qChainLen1 + qChainLen1], qChainKey1);
                    for (size_t qChainIdx2 = qChainIdx1+1; qChainIdx2 < qChainKeys.size(); qChainIdx2++ ){
                        unsigned int qChainKey2 = qChainKeys[qChainIdx2];
                        unsigned int qChainDbId2 = qDbr->sequenceReader->getId(qChainKey2);
                        char *qcadata2 = qStructDbr.getData(qChainDbId2, thread_idx);
                        size_t qCaLength2 = qStructDbr.getEntryLen(qChainDbId2);
                        size_t qChainLen2 = qDbr->sequenceReader->getSeqLen(qChainDbId2);   
                        float* qdata2 = qcoords.read(qcadata2, qChainLen2, qCaLength2);
                        interface->getinterface(qChainLen2, qdata2, &qdata2[qChainLen2], &qdata2[qChainLen2 + qChainLen2], qInterfaceIndex, qChainKey2);
                    }
                    delete interface;
                }
            }
            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++ ) {
                unsigned int qChainKey = qChainKeys[qChainIdx];
                unsigned int qChainAlnId = alnDbr.getId(qChainKey);
                unsigned int qChainDbId = qDbr->sequenceReader->getId(qChainKey);
                // Handling monomer as singleton
                if (qChainAlnId == NOT_AVAILABLE_CHAIN_KEY){
                    char *outpos = Itoa::u32toa_sse2(qComplexId, buffer);
                    result.append(buffer, (outpos - buffer - 1));
                    result.push_back('\n');
                    result5.append(qComplex.complexName + "\t" + tComplexes[qComplexIdx].complexName + "\t1.000000\t1.000000\t1.000000\t1.000000\n");
                    break;
                }

                char *qcadata = qStructDbr.getData(qChainDbId, thread_idx);
                size_t qCaLength = qStructDbr.getEntryLen(qChainDbId);
                size_t qChainLen = qDbr->sequenceReader->getSeqLen(qChainDbId);
                float* qdata = qcoords.read(qcadata, qChainLen, qCaLength);
                
                char *data = alnDbr.getData(qChainAlnId, thread_idx);
                while (*data != '\0' ) {
                    ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                
                    if (!retComplex.isValid){
                        Debug(Debug::ERROR) << "No scorecomplex result provided";
                        EXIT(EXIT_FAILURE);
                    }

                    data = Util::skipLine(data);
                    unsigned int assId = retComplex.assId;
                    unsigned int tChainKey = res.dbKey;
                    unsigned int tChainAlnId = alnDbr.getId(tChainKey);
                    //if target is monomer, break to be singleton
                    if (tChainAlnId == NOT_AVAILABLE_CHAIN_KEY){
                        break;
                    }
                    unsigned int tChainDbId = tDbr->sequenceReader->getId(tChainKey);
                    unsigned int tComplexId = tChainKeyToComplexIdMap.at(tChainKey);
                    unsigned int tComplexIdx = tComplexIdToIdx.at(tComplexId);
                    std::vector<unsigned int> tChainKeys = tComplexes[tComplexIdx].chainKeys;
                    
                    float u[3][3];
                    float t[3];
                    fillUArr(retComplex.uString, u);
                    fillTArr(retComplex.tString, t);
                    unsigned int qalnlen = (std::max(res.qStartPos, res.qEndPos) - std::min(res.qStartPos, res.qEndPos) + 1);
                    unsigned int talnlen = (std::max(res.dbStartPos, res.dbEndPos) - std::min(res.dbStartPos, res.dbEndPos) + 1);
                    if (localComplexMap.find(assId) == localComplexMap.end()) {
                        ComplexFilterCriteria cmplfiltcrit = ComplexFilterCriteria(tComplexId, retComplex.qTmScore, retComplex.tTmScore, t, u);
                        localComplexMap[assId] = cmplfiltcrit;
                    }
                    ComplexFilterCriteria &cmplfiltcrit = localComplexMap.at(assId);
                    cmplfiltcrit.updateAln(qalnlen, talnlen);
                    if (par.filtChainTmThr > 0.0f || par.filtInterfaceLddtThr > 0.0f) {
                        char *tcadata = tStructDbr->getData(tChainDbId, thread_idx);
                        size_t tCaLength = tStructDbr->getEntryLen(tChainDbId);
                        float* tdata = tcoords.read(tcadata, res.dbLen, tCaLength);
                        unsigned int alnLen = cigarToAlignedLength(res.backtrace);
                        cmplfiltcrit.updateChainTmAndInterface(qChainKey, tChainKey, alnLen, qdata, tdata, res.backtrace, res.qStartPos, res.dbStartPos, res.qLen, res.dbLen, t, u, qInterfaceIndex);
                    }
                } // while end
            }
            
            // Filter the target complexes and get the best alignment
            for (auto& assId_res : localComplexMap){
                unsigned int tComplexId  = assId_res.second.targetComplexId;
                unsigned int tComplexIdx = tComplexIdToIdx.at(tComplexId);
                Complex  tComplex = tComplexes[tComplexIdx];

                ComplexFilterCriteria &cmplfiltcrit = assId_res.second;
                cmplfiltcrit.calcCov(qComplex.complexLength, tComplex.complexLength);

                if (par.filtInterfaceLddtThr > 0.0) {
                    cmplfiltcrit.computeInterfaceLddt();
                }

                // Check if the criteria are met
                if (!(cmplfiltcrit.satisfy(par.covMode, par.covThr, par.filtMultimerTmThr, par.filtChainTmThr, par.filtInterfaceLddtThr, qComplex.nChain, tComplex.nChain))){\
                    continue;
                }

                unsigned int alnlen = adjustAlnLen(cmplfiltcrit.qTotalAlnLen, cmplfiltcrit.tTotalAlnLen, par.covMode);
                // Get the best alignement per each target complex
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
                ComplexFilterCriteria &cmplfiltcrit = localComplexMap.at(assId);
                unsigned int tComplexId = cmplfiltcrit.targetComplexId;
                unsigned int tComplexIdx = tComplexIdToIdx.at(tComplexId);
                Complex tComplex = tComplexes[tComplexIdx];
                
                char *outpos = Itoa::u32toa_sse2(tComplexId, buffer);
                result.append(buffer, (outpos - buffer - 1));
                result.push_back('\n');
                result5.append(qComplex.complexName + "\t" + tComplex.complexName + "\t" + std::to_string(cmplfiltcrit.qCov) + "\t" + std::to_string(cmplfiltcrit.tCov) + "\t"+ std::to_string(cmplfiltcrit.qTm)+"\t"+ std::to_string(cmplfiltcrit.tTm)+"\t"+ std::to_string(cmplfiltcrit.interfaceLddt)+"\n");
            }

            #pragma omp critical
            {
                qComplexIdResult[qComplexId]= result;
            }
            result.clear();
            localComplexMap.clear();
            cmplIdToBestAssId.clear();
            selectedAssIDs.clear();
            qInterfaceIndex.clear();
        } // for end
        #pragma omp critical
        {
            resultWrite5.writeData(result5.c_str(), result5.length(), 0);
            result5.clear();
        }
    } // MP end

    for (auto &pair : qComplexIdResult){
        resultWriter.writeData(pair.second.c_str(), pair.second.length(), pair.first);
    }
    
    resultWriter.close(true);
    resultWrite5.close(par.dbOut == false);
    qStructDbr.close();
    alnDbr.close();
    delete qDbr;
    if (sameDB == false) {
        delete tDbr;
        delete tStructDbr;
    }
    qChainKeyToComplexIdMap.clear();
    tChainKeyToComplexIdMap.clear();
    qComplexes.clear();
    tComplexes.clear();
    qComplexIdResult.clear();
    
    return EXIT_SUCCESS;
}