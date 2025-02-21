#include "DBWriter.h"
#include "Util.h"
#include "LocalParameters.h"
#include "DBReader.h"
#include "FileUtil.h"
#include "Coordinate16.h"
#include "tmalign/basic_fun.h"
#include "MultimerUtil.h"
#include "LDDT.h"
#include <map>
#include <set>

#ifdef OPENMP
#include <omp.h>
#endif
#define INTERFACE_THRESHOLD 8

unsigned int cigarToAlignedLength(const std::string &cigar) {
    std::string backtrace = Matcher::uncompressAlignment(cigar);
    unsigned int alni = 0;
    for (size_t btPos = 0; btPos < backtrace.size(); btPos++) {
        if (backtrace[btPos] == 'M') {
            alni++;
        }
    }
    return alni;
}

struct Complex {
    int complexId;
    unsigned int nChain;
    unsigned int complexLength;
    std::string complexName;
    // std::vector<unsigned int> chainLengths;
    std::vector<unsigned int> chainKeys;
    
    Complex() : complexId(0), nChain(0), complexLength(0), complexName("") {}
    ~Complex() {
        chainKeys.clear();
        // chainLengths.clear();
    }
};

typedef Coordinates AlignedCoordinate;

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

struct chainAlignment {
    unsigned int qKey;
    unsigned int tKey;
    unsigned int qLen;
    unsigned int tLen;
    unsigned int alnLen;
    unsigned int qStartPos;
    unsigned int tStartPos;
    std::string cigar;
    chainAlignment() : qKey(0), tKey(0), qLen(0), tLen(0), alnLen(0), qStartPos(0), tStartPos(0), cigar("") {}
    chainAlignment(unsigned int qKey, unsigned int tKey, unsigned int qLen, unsigned int tLen, unsigned int alnLen, unsigned int qStartPos, unsigned int tStartPos, const std::string &cigar) : 
        qKey(qKey), tKey(tKey), qLen(qLen), tLen(tLen), alnLen(alnLen), qStartPos(qStartPos), tStartPos(tStartPos), cigar(cigar) {}
    ~chainAlignment() {}
};

class ComplexFilterCriteria {
public:
    unsigned int targetComplexId;

    // per complex
    unsigned int qTotalAlnLen;
    unsigned int tTotalAlnLen;
    unsigned int interfaceAlnLen;
    float qCov;
    float tCov;
    float interfaceLddt;
    float qTm;
    float tTm;
    float avgTm;
    float t[3];
    float u[3][3];

    // per chain : criteria for chainTmThr & lddtThr
    std::vector<float> qAlnChainTms;
    std::vector<float> tAlnChainTms;
    std::vector<chainAlignment> alignedChains;

    ComplexFilterCriteria() {}
    ComplexFilterCriteria(
        unsigned int targetComplexId, float qTm, float tTm, float tstring[3], float ustring[3][3]
    ) :
        targetComplexId(targetComplexId), qTotalAlnLen(0), tTotalAlnLen(0), interfaceAlnLen(0),
        qCov(0), tCov(0), interfaceLddt(0), qTm(qTm), tTm(tTm), avgTm(0)
    {
        std::copy(tstring, tstring + 3, t);
        for (int i = 0; i < 3; i++) {
            std::copy(ustring[i], ustring[i] + 3, u[i]);
        }
    }

    ~ComplexFilterCriteria() {
        qAlnChainTms.clear();
        tAlnChainTms.clear();
        alignedChains.clear();
    }

    bool hasTm(float TmThr, int covMode) {
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

    bool hasChainTm(float chainTmThr, int covMode, unsigned int qChainNum, unsigned int tChainNum) {
        if (alignedChains.size()<std::min(qChainNum, tChainNum)) {
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

    bool hasChainNum(int covMode, unsigned int qChainNum, unsigned int tChainNum) {
        switch (covMode) {
            case Parameters::COV_MODE_BIDIRECTIONAL:
                if (qChainNum != tChainNum) {
                    return false;
                }
                break;
            default:
                return true;
        }
        return true;
    }

    // void calculateAvgTm(int covMode){
    //     switch (covMode) {
    //         case Parameters::COV_MODE_BIDIRECTIONAL:
    //             avgTm = ( qTm + tTm ) / 2 ;
    //             break;
    //         case Parameters::COV_MODE_TARGET:
    //             avgTm = tTm ;
    //             break;
    //         case Parameters::COV_MODE_QUERY:
    //             avgTm = qTm ;
    //             break;
    //         default :
    //             avgTm = ( qTm + tTm ) / 2 ;
    //     }
    // }

    bool hasInterfaceLDDT(float iLddtThr, unsigned int qChainNum, unsigned int tChainNum) {
        if (alignedChains.size()<std::min(qChainNum, tChainNum)) {
            return false;
        }
        return(interfaceLddt >= iLddtThr);
    }

    bool satisfy(int covMode, float covThr, float TmThr, float chainTmThr, float iLddtThr, size_t qChainNum, size_t tChainNum ) {
        const bool covOK = covThr ? Util::hasCoverage(covThr, covMode, qCov, tCov) : true;
        const bool TmOK = TmThr ? hasTm(TmThr, covMode) : true;
        const bool chainNumOK = hasChainNum(covMode, qChainNum, tChainNum);
        const bool chainTmOK = chainTmThr ? hasChainTm(chainTmThr, covMode, qChainNum, tChainNum) : true; 
        const bool lddtOK = iLddtThr ? hasInterfaceLDDT(iLddtThr, qChainNum, tChainNum) : true; 
        // calculateAvgTm(covMode);
        return (covOK && TmOK && chainNumOK && chainTmOK && lddtOK); 
    }

    bool satisfy_first(int covMode, float covThr, float TmThr, size_t qChainNum, size_t tChainNum ) {
        const bool covOK = covThr ? Util::hasCoverage(covThr, covMode, qCov, tCov) : true;
        const bool TmOK = TmThr ? hasTm(TmThr, covMode) : true;
        const bool chainNumOK = hasChainNum(covMode, qChainNum, tChainNum);
        return (covOK && TmOK && chainNumOK ); 
    }

    bool satisfy_second(int covMode, float chainTmThr, float iLddtThr, size_t qChainNum, size_t tChainNum ) {
        const bool chainTmOK = chainTmThr ? hasChainTm(chainTmThr, covMode, qChainNum, tChainNum) : true; 
        const bool lddtOK = iLddtThr ? hasInterfaceLDDT(iLddtThr, qChainNum, tChainNum) : true; 
        return (chainTmOK && lddtOK); 
    }

    void updateAln(unsigned int qAlnLen, unsigned int tAlnLen) {
        qTotalAlnLen += qAlnLen;
        tTotalAlnLen += tAlnLen;
    }

    void computeChainTmScore(AlignedCoordinate &qchain, AlignedCoordinate &tchain, unsigned int totalAlnLen) {
        AlignedCoordinate tmt(totalAlnLen);
        BasicFunction::do_rotation(tchain, tmt, totalAlnLen, t, u);

        unsigned int chainOffset = 0;
        for (unsigned int i=0; i<alignedChains.size(); i++) {
            chainAlignment &chainaln = alignedChains[i];
            unsigned int qLen = chainaln.qLen;
            unsigned int tLen = chainaln.tLen;
            unsigned int alnLen = chainaln.alnLen;
    
            float d0 = 1.24*(cbrt(tLen-15)) -1.8;
            float d02 = d0*d0;

            float tmScore = 0;
            for (unsigned int ci=chainOffset; ci<chainOffset+alnLen; ci++) {
                float xa_x = qchain.x[ci];
                float xa_y = qchain.y[ci];
                float xa_z = qchain.z[ci];
                float ya_x = tmt.x[ci];
                float ya_y = tmt.y[ci];
                float ya_z = tmt.z[ci];
                float di = BasicFunction::dist(xa_x, xa_y, xa_z, ya_x, ya_y, ya_z);
                float oneDividedDist = 1/(1+di/d02);
                tmScore += oneDividedDist;
            }

            float qtmscore = tmScore / qLen;
            float ttmscore = tmScore / tLen;
            updateChainTmScore(qtmscore, ttmscore);
            chainOffset += alnLen;
            
            // TODO: Implement in SIMD
            // simd_float vd02 = simdf32_set(d02);
            // simd_float one = simdf32_set(1.0);
            // simd_float acc = simdf32_set(0.0);
            // std::cout << "AlnKey: " << i << " AlnLen: " << alnLen << std::endl;
            // for (unsigned int ci=chainOffset; ci<chainOffset+alnLen-VECSIZE_FLOAT; ci+=VECSIZE_FLOAT) {
            //     std::cout << ci << std::endl;
            //     simd_float xa_x = simdf32_load(&qchain.x[ci]);
            //     // simd_float xa_y = simdf32_load(&qchain.y[ci]);
            //     // simd_float xa_z = simdf32_load(&qchain.z[ci]);
            //     // simd_float ya_x = simdf32_load(&tmt.x[ci]);
            //     // simd_float ya_y = simdf32_load(&tmt.y[ci]);
            //     // simd_float ya_z = simdf32_load(&tmt.z[ci]);
            //     // ya_x = simdf32_sub(xa_x, ya_x);
            //     // ya_y = simdf32_sub(xa_y, ya_y);
            //     // ya_z = simdf32_sub(xa_z, ya_z);
            //     // simd_float di = simdf32_add(simdf32_add(simdf32_mul(xa_x, xa_x), simdf32_mul(xa_y, xa_y)), simdf32_mul(xa_z, xa_z));
            //     // simd_float oneDividedDist = simdf32_div(one, simdf32_add(one, simdf32_div(di,vd02)));
            //     // acc = simdf32_add(acc, oneDividedDist);
            // }

            // // float sumArray[VECSIZE_FLOAT];
            // // float tmscore = 0;
            // // simdf32_store(sumArray, acc);
            // // for (size_t j = 0; j < VECSIZE_FLOAT; i++) {
            // //     tmscore += sumArray[i];
            // // }
            // // float qtmscore = tmscore / qLen;
            // // float ttmscore = tmscore / tLen;
            // // updateChainTmScore(qtmscore, ttmscore);
            // chainOffset += alnLen;
        }
    }

    void updateChainTmScore(float qChainTm, float tChainTm) {
        qAlnChainTms.push_back(qChainTm);
        tAlnChainTms.push_back(tChainTm);
    }

    void fillComplexAlignment(chainAlignment &alnchain, 
        // const std::string &cigar, int qStartPos, int tStartPos, int qLen, int tLen, 
        unsigned int &chainOffset, float *qdata, float *tdata, AlignedCoordinate &qAlnCoords, AlignedCoordinate &tAlnCoords) {
        int mi = chainOffset;
        int qi = alnchain.qStartPos;
        int ti = alnchain.tStartPos;
        int qLen = alnchain.qLen;
        int tLen = alnchain.tLen;
        std::string backtrace = Matcher::uncompressAlignment(alnchain.cigar);
                
        for (size_t btPos = 0; btPos < backtrace.size(); btPos++) {
            if (backtrace[btPos] == 'M') {
                qAlnCoords.x[mi] = qdata[qi];
                qAlnCoords.y[mi] = qdata[qLen + qi];
                qAlnCoords.z[mi] = qdata[2*qLen + qi];
                tAlnCoords.x[mi] = tdata[ti];
                tAlnCoords.y[mi] = tdata[tLen + ti];
                tAlnCoords.z[mi] = tdata[2*tLen + ti];
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
        chainOffset = mi;
    }

    void calcCov(unsigned int qLen, unsigned int tLen) {
        qCov = static_cast<float>(qTotalAlnLen) / static_cast<float>(qLen);
        tCov = static_cast<float>(tTotalAlnLen) / static_cast<float>(tLen);
    }

    void computeInterfaceLddt(AlignedCoordinate &qAlnCoords, AlignedCoordinate &tAlnCoords, float threshold = INTERFACE_THRESHOLD) {
        if (alignedChains.size() == 1) { // No interface if only one chain aligned
            interfaceLddt = 1;
            return;
        }
        std::vector<unsigned int> chainOffsets(alignedChains.size(), 0);
        unsigned int acc = 0;
        for (size_t i = 0; i < alignedChains.size(); i++) {
            chainOffsets[i] = acc;
            acc += alignedChains[i].alnLen;
        }
        
        float t2 = threshold * threshold;

        std::set<unsigned int> interfacePos;    
        unsigned int intLen = 0;

        // Find and save interface Coordinates
        for (size_t chainIdx = 0; chainIdx < chainOffsets.size(); chainIdx++) {
            unsigned int c1_start = chainOffsets[chainIdx];
            unsigned int c1_end = c1_start + alignedChains[chainIdx].alnLen;
            for (size_t resIdx1 = c1_start; resIdx1 < c1_end; resIdx1++) {
                for (size_t resIdx2 = c1_end; resIdx2 < acc; resIdx2++) { // Rest of the chainss
                    float dist = BasicFunction::dist(qAlnCoords.x[resIdx1], qAlnCoords.y[resIdx1], qAlnCoords.z[resIdx1],
                                                    qAlnCoords.x[resIdx2], qAlnCoords.y[resIdx2], qAlnCoords.z[resIdx2]);
                    if (dist < t2) {
                        if (interfacePos.find(resIdx1) == interfacePos.end()) {
                            interfacePos.insert(resIdx1);
                            intLen++;
                        }
                        if (interfacePos.find(resIdx2) == interfacePos.end()) {
                            interfacePos.insert(resIdx2);
                            intLen++;
                        }
                    }
                }
            }
        }

        interfaceAlnLen = intLen;

        if (intLen == 0) {
            return;
        }

        AlignedCoordinate qInterface(intLen);
        AlignedCoordinate tInterface(intLen);
        size_t idx = 0;
        //     // if (qInterfacePos[chainIdx].size() >= 4) { // TODO: Is it important? then change interfacePos into vector. But it can cause (intLen > idx) + downstream errors in lddt calculation
        for (size_t resIdx: interfacePos) {
            qInterface.x[idx] = qAlnCoords.x[resIdx];
            qInterface.y[idx] = qAlnCoords.y[resIdx];
            qInterface.z[idx] = qAlnCoords.z[resIdx];
            tInterface.x[idx] = tAlnCoords.x[resIdx];
            tInterface.y[idx] = tAlnCoords.y[resIdx];
            tInterface.z[idx] = tAlnCoords.z[resIdx];
            idx++;
        }
            // }    

        std::string bt(intLen, 'M');
        LDDTCalculator *lddtcalculator = NULL;
        lddtcalculator = new LDDTCalculator(intLen+1, intLen+1);
        lddtcalculator->initQuery(intLen, &qInterface.x[0], &qInterface.y[0], &qInterface.z[0]);
        LDDTCalculator::LDDTScoreResult lddtres = lddtcalculator->computeLDDTScore(intLen, 0, 0, bt, &tInterface.x[0], &tInterface.y[0], &tInterface.z[0]);
        interfaceLddt = lddtres.avgLddtScore;
        delete lddtcalculator;
    }
};

char* filterToBuffer(ComplexFilterCriteria cmplfiltcrit, char* tmpBuff){
    *(tmpBuff-1) = '\t';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.qCov, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.tCov, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.qTm, tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.tTm, tmpBuff);
    *(tmpBuff-1) = '\t';

    for (unsigned int i = 0; i < cmplfiltcrit.qAlnChainTms.size(); i++) {
        tmpBuff = fastfloatToBuffer(cmplfiltcrit.qAlnChainTms[i], tmpBuff);
        *(tmpBuff-1) = ',';
    }
    *(tmpBuff-1) = '\t';
    for (unsigned int i = 0; i < cmplfiltcrit.tAlnChainTms.size(); i++) {
        tmpBuff = fastfloatToBuffer(cmplfiltcrit.tAlnChainTms[i], tmpBuff);
        *(tmpBuff-1) = ',';
    }
    *(tmpBuff-1) = '\t';

    tmpBuff = fastfloatToBuffer(cmplfiltcrit.interfaceLddt, tmpBuff);    
    *(tmpBuff-1) = '\t';

    // tmpBuff = Itoa::u32toa_sse2(cmplfiltcrit.interfaceAlnLen, tmpBuff);
    // *(tmpBuff-1) = '\t';

    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[0][0], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[0][1], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[0][2], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[1][0], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[1][1], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[1][2], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[2][0], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[2][1], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.u[2][2], tmpBuff);
    *(tmpBuff-1) = '\t';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.t[0], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.t[1], tmpBuff);
    *(tmpBuff-1) = ',';
    tmpBuff = fastfloatToBuffer(cmplfiltcrit.t[2], tmpBuff);
    *(tmpBuff-1) = '\n';
    return tmpBuff;
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
            if (id == NOT_AVAILABLE_CHAIN_KEY) {
                break;
            }
            unsigned int reslen = Dbr->sequenceReader->getSeqLen(id);
            // complex->chainLengths.push_back(reslen);
            cmpllen += reslen;
        }
        complex->complexLength = cmpllen;
    }
}

static void getlookupInfo(
        IndexReader* dbr,
        const std::string &file,
        std::map<unsigned int, unsigned int> &chainKeyToComplexIdLookup,
        std::vector<Complex> &complexes,
        std::map<unsigned int, unsigned int> &complexIdtoIdx,
        std::map<unsigned int, std::string> &chainKeyToChainNameMap
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
        unsigned int chainDbId = dbr->sequenceReader->getId(chainKey);
        if (chainDbId != NOT_AVAILABLE_CHAIN_KEY) {
            auto complexId = Util::fast_atoi<int>(entry[2]);
            chainKeyToComplexIdLookup.emplace(chainKey, complexId);
            std::string chainName(entry[1], (entry[2] - entry[1]) - 1);
            size_t lastUnderscoreIndex = chainName.find_last_of('_');
            std::string complexName = chainName.substr(0, lastUnderscoreIndex);

            chainName = chainName.substr(lastUnderscoreIndex + 1, chainName.size()); // 7soy_1.pdb_A -> A
            chainKeyToChainNameMap.emplace(chainKey, chainName);

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
        }
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
#ifdef OPENMP
localThreads = std::max(std::min((size_t)par.threads, alnDbr.getSize()), (size_t)1);
#endif
    const bool shouldCompress = (par.compressed == true);
    const int db4Type = Parameters::DBTYPE_CLUSTER_RES;
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), par.threads, shouldCompress, db4Type);
    resultWriter.open();
    const int db5Type = Parameters::DBTYPE_GENERIC_DB;
    DBWriter resultWrite5((par.db4 + "_info").c_str(), (par.db4 + "_info.index").c_str(), par.threads, shouldCompress, db5Type);
    resultWrite5.open();

    std::string qLookupFile = par.db1 + ".lookup";
    std::string tLookupFile = par.db2 + ".lookup";
    
    std::vector<Complex> qComplexes, tComplexes;
    std::map<unsigned int, unsigned int> qComplexIdToIdx, tComplexIdToIdx;
    chainKeyToComplexId_t qChainKeyToComplexIdMap, tChainKeyToComplexIdMap;
    chainKeyToChainName_t qChainKeyToChainNameMap, tChainKeyToChainNameMap;

    getlookupInfo(qDbr, qLookupFile, qChainKeyToComplexIdMap, qComplexes, qComplexIdToIdx, qChainKeyToChainNameMap);
    getComplexResidueLength(qDbr, qComplexes);
    Debug::Progress progress(qComplexes.size());

    if (sameDB) {
        tChainKeyToComplexIdMap = qChainKeyToComplexIdMap;
        tChainKeyToChainNameMap = qChainKeyToChainNameMap;
        tComplexes = qComplexes;
        tComplexIdToIdx = qComplexIdToIdx;
    } else {
        getlookupInfo(tDbr, tLookupFile, tChainKeyToComplexIdMap, tComplexes, tComplexIdToIdx, tChainKeyToChainNameMap);
        getComplexResidueLength(tDbr, tComplexes);
    }
    
#pragma omp parallel num_threads(localThreads) 
    {   
        char buffer[32];
        char buffer2[4096];
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        resultToWrite_t result;
        std::map<unsigned int, ComplexFilterCriteria> localComplexMap;
        std::map<unsigned int, std::vector<unsigned int>> cmplIdToBestAssId;
        std::vector<unsigned int> selectedAssIDs;

        Matcher::result_t res;
#pragma omp for schedule(dynamic, 1) 
        for (size_t qComplexIdx = 0; qComplexIdx < qComplexes.size(); qComplexIdx++) {
            progress.updateProgress();
            Complex qComplex = qComplexes[qComplexIdx];
            unsigned int qComplexId = qComplex.complexId;
            std::vector<unsigned int> qChainKeys = qComplex.chainKeys;
            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++)
            {
                unsigned int qChainKey = qChainKeys[qChainIdx];
                unsigned int qChainAlnId = alnDbr.getId(qChainKey);
                // unsigned int qChainDbId = qDbr->sequenceReader->getId(qChainKey);
                // Handling monomer as singleton
                if (qChainAlnId == NOT_AVAILABLE_CHAIN_KEY) {
                    break;
                }
                
                char *data = alnDbr.getData(qChainAlnId, thread_idx);
                while (*data != '\0' ) {
                    ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                
                    if (!retComplex.isValid) {
                        Debug(Debug::ERROR) << "No scorecomplex result provided";
                        EXIT(EXIT_FAILURE);
                    }

                    data = Util::skipLine(data);
                    unsigned int assId = retComplex.assId;
                    unsigned int tChainKey = res.dbKey;
                    // unsigned int tChainDbId = tDbr->sequenceReader->getId(tChainKey); // RECOVER
                    unsigned int tComplexId = tChainKeyToComplexIdMap.at(tChainKey);

                    //if target is monomer, but user doesn't want, continue
                    unsigned int tChainAlnId = alnDbr.getId(tChainKey);
                    if (tChainAlnId == NOT_AVAILABLE_CHAIN_KEY) {
                        continue;
                    }
                    float u[3][3];
                    float t[3];
                    fillUArr(retComplex.uString, u);
                    fillTArr(retComplex.tString, t);
                    unsigned int qalnlen = (std::max(res.qStartPos, res.qEndPos) - std::min(res.qStartPos, res.qEndPos) + 1);
                    unsigned int talnlen = (std::max(res.dbStartPos, res.dbEndPos) - std::min(res.dbStartPos, res.dbEndPos) + 1);

                    if (localComplexMap.find(assId) == localComplexMap.end()) {
                        ComplexFilterCriteria cmplfiltcrit(tComplexId, retComplex.qTmScore, retComplex.tTmScore, t, u);
                        localComplexMap[assId] = cmplfiltcrit;
                    }
                    ComplexFilterCriteria &cmplfiltcrit = localComplexMap.at(assId);
                    cmplfiltcrit.updateAln(qalnlen, talnlen);
    
                    unsigned int matchLen = cigarToAlignedLength(res.backtrace);
                    chainAlignment chainaln = chainAlignment(qChainKey, tChainKey, res.qLen, res.dbLen, matchLen, res.qStartPos, res.dbStartPos, res.backtrace);
                    cmplfiltcrit.alignedChains.push_back(chainaln);
                } // while end
            }

            // Filter the target complexes and get the best alignment
            for (auto& assId_res : localComplexMap) {
                unsigned int tComplexId  = assId_res.second.targetComplexId;
                
                unsigned int tComplexIdx = tComplexIdToIdx.at(tComplexId);
                Complex  &tComplex = tComplexes[tComplexIdx];

                ComplexFilterCriteria &cmplfiltcrit = assId_res.second;
                cmplfiltcrit.calcCov(qComplex.complexLength, tComplex.complexLength);

                // Check the criteria first before calculating the rest (chain-tm-score, interface-lddt)
                if (!(cmplfiltcrit.satisfy_first(par.covMode, par.covThr, par.filtMultimerTmThr, qComplex.nChain, tComplex.nChain))) {
                    continue;
                }

                if (par.filtChainTmThr || par.filtInterfaceLddtThr) {
                    /* Fill aligned coords */
                    unsigned int totalAlnLen = 0;
                    for (size_t i = 0; i < cmplfiltcrit.alignedChains.size(); i++) {
                        totalAlnLen += cmplfiltcrit.alignedChains[i].alnLen;
                    }

                    AlignedCoordinate qAlnCoords = AlignedCoordinate(totalAlnLen);
                    AlignedCoordinate tAlnCoords = AlignedCoordinate(totalAlnLen);
                    Coordinate16 qcoords, tcoords;
                    unsigned int chainOffset = 0;
                    
                    std::set<chainToResidue> interface_chain2residue;
                    for (size_t chainIdx = 0; chainIdx < cmplfiltcrit.alignedChains.size(); chainIdx++) {
                        // Bring Coordinates from cadb
                        chainAlignment &alnchain = cmplfiltcrit.alignedChains[chainIdx];
                        unsigned int qChainKey = alnchain.qKey;
                        unsigned int qChainDbId = qDbr->sequenceReader->getId(qChainKey);
                        char *qcadata = qStructDbr.getData(qChainDbId, thread_idx);
                        size_t qCaLength = qStructDbr.getEntryLen(qChainDbId);
                        size_t qChainLen = qDbr->sequenceReader->getSeqLen(qChainDbId);
                        float* qdata = qcoords.read(qcadata, qChainLen, qCaLength);
                        
                        unsigned int tChainKey = alnchain.tKey;
                        unsigned int tChainDbId = tDbr->sequenceReader->getId(tChainKey);
                        char *tcadata = tStructDbr->getData(tChainDbId, thread_idx);
                        size_t tCaLength = tStructDbr->getEntryLen(tChainDbId);
                        size_t tChainLen = tDbr->sequenceReader->getSeqLen(tChainDbId);
                        float* tdata = tcoords.read(tcadata, tChainLen, tCaLength);

                        // Save each chain into Alignedcoords
                        cmplfiltcrit.fillComplexAlignment(alnchain, chainOffset, qdata, tdata, qAlnCoords, tAlnCoords);

                        // Retrieve interface residues from Query complex
                        if (par.filtInterfaceLddtThr) {
                            // std::vector<chainToResidue> local_chain2residue;
                            float d2 = INTERFACE_THRESHOLD * INTERFACE_THRESHOLD;
                            for (size_t otherIdx = chainIdx+1; otherIdx < cmplfiltcrit.alignedChains.size(); otherIdx++) {
                                chainAlignment &otherchain = cmplfiltcrit.alignedChains[otherIdx];
                                unsigned int qChainKey2 = otherchain.qKey;
                                unsigned int qChainDbId2 = qDbr->sequenceReader->getId(qChainKey2);
                                char *qcadata2 = qStructDbr.getData(qChainDbId2, thread_idx);
                                size_t qCaLength2 = qStructDbr.getEntryLen(qChainDbId2);
                                size_t qChainLen2 = qDbr->sequenceReader->getSeqLen(qChainDbId2);
                                float* qdata2 = qcoords.read(qcadata2, qChainLen2, qCaLength2);

                                for (size_t chainResIdx=0; chainResIdx < qChainLen; chainResIdx++) {
                                    bool isInterface = false;
                                    for (size_t otherResIdx=0; otherResIdx < qChainLen2; otherResIdx++) {
                                        float dist = BasicFunction::dist(qdata[chainResIdx], qdata[qChainLen + chainResIdx], qdata[2*qChainLen + chainResIdx],
                                                                        qdata2[otherResIdx], qdata2[qChainLen2 + otherResIdx], qdata2[2*qChainLen2 + otherResIdx]);
                                        if (dist < d2) {
                                            isInterface = true;
                                            // local_chain2residue.push_back({otherIdx, otherResIdx});
                                            interface_chain2residue.insert({otherIdx, otherResIdx});
                                        }
                                    }
                                    if (isInterface) {
                                        // local_chain2residue.push_back({chainIdx, chainResIdx});
                                        interface_chain2residue.insert({chainIdx, chainResIdx});
                                    }
                                }
                            }
                            // interface_chain2residue.insert(local_chain2residue.begin(), local_chain2residue.end());
                        }
                    }

                    if (par.filtChainTmThr > 0.0) {
                        cmplfiltcrit.computeChainTmScore(qAlnCoords, tAlnCoords, totalAlnLen);
                        if (!(cmplfiltcrit.hasChainTm(par.filtChainTmThr, par.covMode, qComplex.nChain, tComplex.nChain))) {
                            continue;
                        }
                    }

                    if (par.filtInterfaceLddtThr > 0.0) {
                        cmplfiltcrit.computeInterfaceLddt(qAlnCoords, tAlnCoords);
                        if (!(cmplfiltcrit.hasInterfaceLDDT(par.filtInterfaceLddtThr, qComplex.nChain, tComplex.nChain))) {
                            continue;
                        }
                        // if (cmplfiltcrit.alignedChains.size() > 1 && (cmplfiltcrit.interfaceAlnLen > interface_chain2residue.size())) {
                        //     std::cout << "QComplex: " << qComplex.complexName << " N aligned: " << cmplfiltcrit.alignedChains.size() <<" Total AlnLen: " << totalAlnLen << " Interface Len: " << interface_chain2residue.size() << " Interface AlnLen: " << cmplfiltcrit.interfaceAlnLen << " Interface LDDT: " << cmplfiltcrit.interfaceLddt << std::endl;
                        // }                        
                    }

                    // if (!(cmplfiltcrit.satisfy_second(par.covMode, par.filtChainTmThr, par.filtInterfaceLddtThr, qComplex.nChain, tComplex.nChain))) {
                    //     continue;
                    // }
                }

                // Check if the rest criteria are met
                // if (!(cmplfiltcrit.satisfy(par.covMode, par.covThr, par.filtMultimerTmThr, par.filtChainTmThr, par.filtInterfaceLddtThr, qComplex.nChain, tComplex.nChain))) {
                //     continue;
                // }

                unsigned int alnlen = adjustAlnLen(cmplfiltcrit.qTotalAlnLen, cmplfiltcrit.tTotalAlnLen, par.covMode);
                // Get the best alignement per each target complex   
                if (cmplIdToBestAssId.find(tComplexId) == cmplIdToBestAssId.end()) {
                    cmplIdToBestAssId[tComplexId] = {assId_res.first, alnlen};
                    // cmplIdToBestAssId[tComplexId] = {static_cast<float>(assId_res.first), cmplfiltcrit.avgTm};
                    // cmplIdToBestAssId[tComplexId] = {static_cast<float>(assId), cmplfiltcrit.avgTm};
                } else {
                    if (alnlen > cmplIdToBestAssId.at(tComplexId)[1]) {
                        cmplIdToBestAssId[tComplexId] = {assId_res.first, alnlen};
                    }

                }
                
            }

            for (const auto& pair : cmplIdToBestAssId) {
                selectedAssIDs.push_back(pair.second[0]);
            }
            if (selectedAssIDs.size() == 0) {
                float t[3];
                float u[3][3];
                for (int i=0; i < 3; i++) {
                    t[i] = 0.0;
                }
                for (int i=0; i < 3; i++) {
                    for (int j=0; j < 3; j++) {
                        u[i][j] = 0.0;
                    }
                }
                ComplexFilterCriteria cmplfiltcrit(qComplexId, 1.0, 1.0, t, u);
                cmplfiltcrit.qCov = 1.0;
                cmplfiltcrit.tCov = 1.0;
                cmplfiltcrit.interfaceLddt = 1.0;

                selectedAssIDs.push_back(0);
                localComplexMap.insert({0, cmplfiltcrit});
            }
            
            resultWrite5.writeStart(thread_idx);
            for (unsigned int assIdidx = 0; assIdidx < selectedAssIDs.size(); assIdidx++) {
                unsigned int assId = selectedAssIDs.at(assIdidx);
                ComplexFilterCriteria &cmplfiltcrit = localComplexMap.at(assId);

                unsigned int tComplexId = cmplfiltcrit.targetComplexId;
                unsigned int tComplexIdx = tComplexIdToIdx.at(tComplexId);
                Complex tComplex = tComplexes.at(tComplexIdx);

                std::string qComplexName = qComplex.complexName;
                std::string tComplexName = tComplex.complexName;
                std::string qChainNames = "";
                std::string tChainNames = "";

                // Add chain names if chain alignment is saved
                if (cmplfiltcrit.alignedChains.size() > 0) { // If chain alignment is saved : chainTmThr & lddtThr
                    for (size_t chainIdx = 0; chainIdx < cmplfiltcrit.alignedChains.size(); chainIdx++) {
                        chainAlignment &alnchain = cmplfiltcrit.alignedChains[chainIdx];
                        std::string qChainName = qChainKeyToChainNameMap.at(alnchain.qKey);
                        std::string tChainName = tChainKeyToChainNameMap.at(alnchain.tKey);
                        if (chainIdx == 0) {
                            qChainNames+= qChainName;
                            tChainNames+= tChainName;
                        } else {
                            qChainNames+= ","+qChainName;
                            tChainNames+= ","+tChainName;
                        }
                    }
                }
                
                char *outpos = Itoa::u32toa_sse2(tComplexId, buffer);
                result.append(buffer, (outpos - buffer - 1));
                result.push_back('\n');

                char * tmpBuff = buffer2 + sprintf(buffer2, "%d\t%s\t%s\t%s", qComplexId, tComplexName.c_str(), qChainNames.c_str(), tChainNames.c_str()) +1;
                tmpBuff = filterToBuffer(cmplfiltcrit, tmpBuff);
                resultWrite5.writeAdd(buffer2, tmpBuff - buffer2, thread_idx);
            }

            resultWriter.writeData(result.c_str(), result.length(), qComplexId, thread_idx);
            resultWrite5.writeEnd(qComplexId, thread_idx);
            result.clear();
            localComplexMap.clear();
            cmplIdToBestAssId.clear();
            selectedAssIDs.clear();
        } // for end
    } // MP end
    
    resultWriter.close(false);
    resultWrite5.close(false);
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
    return EXIT_SUCCESS;
}
