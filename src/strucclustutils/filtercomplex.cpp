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


#ifdef OPENMP
#include <omp.h>
#endif


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
    }
}

// bool hasChainnum(bool sameChainNum, int qChainNum, int tChainNum){
//     switch (sameChainNum){
//         case 1:
//             if (qChainNum != tChainNum){
//                 return false;
//             }else{return true;}
//         case 0:
//             return true;
//     }
// }

bool hasChainTm(float chainTMThr, int covMode, std::vector<double> &qChainTmScores, std::vector<double> &tChainTmScores, unsigned int qChainNum, unsigned int tChainNum) {
    if (chainTMThr > 0 ){
        switch (covMode) {
            case Parameters::COV_MODE_BIDIRECTIONAL:
                if (qChainTmScores.size()<std::min(qChainNum, tChainNum)){
                    return false;
                }
                for (size_t i = 0; i < qChainTmScores.size(); i++) {
                    if (qChainTmScores[i] < chainTMThr || tChainTmScores[i] < chainTMThr) {
                        return false;
                    }
                }
                break;
            case Parameters::COV_MODE_TARGET:
                if (qChainTmScores.size()<std::min(qChainNum, tChainNum)){
                    return false;
                }
                for (size_t i = 0; i < qChainTmScores.size(); i++) {
                    if (tChainTmScores[i] < chainTMThr) {
                        return false;
                    }
                }
                break;
            case Parameters::COV_MODE_QUERY:
                if (qChainTmScores.size()<std::min(qChainNum, tChainNum)){
                    return false;
                }
                for (size_t i = 0; i < qChainTmScores.size(); i++) {
                    if (qChainTmScores[i] < chainTMThr) {
                        return false;
                    }
                }
                break;
            case Parameters::COV_MODE_LENGTH_QUERY :
            case Parameters::COV_MODE_LENGTH_TARGET :
            case Parameters::COV_MODE_LENGTH_SHORTER :
                break;
        }
    }
    return true;
}

struct Calpha
{
    int pos;
    float x,y,z;
    Calpha(int idx, float xVal, float yVal, float zVal)
        : pos(idx), x(xVal), y(yVal), z(zVal) {}
};
struct ChainForInterface {
    int chainid;
    std::vector<Calpha> Calphas;
    ChainForInterface(int idx)
        : chainid(idx) {}
    void addCalpha(const Calpha& calpha){
        Calphas.push_back(calpha);
    }
    Calpha getCalpha(size_t idx) const {
        if (idx < Calphas.size()) {
            return Calphas[idx];
        }
    }
};
struct Interface{
    std::vector<ChainForInterface> chainsinInterface;
    void addChain(const ChainForInterface& chain) {
        chainsinInterface.push_back(chain);
    }
    ChainForInterface getChain(size_t idx) const {
        if (idx < chainsinInterface.size()) {
            return chainsinInterface[idx];
        }
    }
};

struct Complex {
    std::vector<ChainForInterface> chainsinComplex;
    std::vector<Interface> InterfaceVec;
    void addChain(const ChainForInterface& chain) {
        chainsinComplex.push_back(chain);
    }
    void findInterface(float distanceThreshold) {
        for (size_t i = 0; i < chainsinComplex.size(); ++i) {
            for (size_t j = i + 1; j < chainsinComplex.size(); ++j) {
                Interface interface;
                const ChainForInterface& chainA = chainsinComplex[i];
                const ChainForInterface& chainB = chainsinComplex[j];
                ChainForInterface chainInterfaceA(chainA.chainid);
                ChainForInterface chainInterfaceB(chainB.chainid);
                interface.addChain(chainInterfaceA);
                interface.addChain(chainInterfaceB); 
                std::vector<int> AvisitedCa, BvisitedCa;
                for (const auto& calphaA : chainA.Calphas) {
                    for (const auto& calphaB : chainB.Calphas) {
                        float distance = BasicFunction::dist(calphaA.x, calphaA.y, calphaA.z, calphaB.x, calphaB.y, calphaB.z);
                        if (distance < distanceThreshold) {
                            if (std::find(AvisitedCa.begin(), AvisitedCa.end(), calphaA.pos)== AvisitedCa.end()){
                                interface.chainsinInterface[0].addCalpha(calphaA);
                            }
                            if (std::find(BvisitedCa.begin(), BvisitedCa.end(), calphaB.pos)== BvisitedCa.end()){
                                interface.chainsinInterface[1].addCalpha(calphaB);
                            }
                            AvisitedCa.push_back(calphaA.pos);
                            BvisitedCa.push_back(calphaB.pos);
                        }
                    }
                }
                if (!AvisitedCa.empty()){
                    InterfaceVec.push_back(interface);
                }
            }
        }
    }
};

class ComplexFilterCriteria {
public:
    unsigned int dbKey;
    unsigned int qTotalAlnLen;
    unsigned int tTotalAlnLen;
    float qCov;
    float tCov;
    double qTM;
    double tTM;
    Complex qComplex;
    Complex tComplex;
    Complex qnewComplex;
    Complex tnewComplex;
    std::vector<double> alignedQChainTmScores;
    std::vector<double> alignedTChainTmScores;
    ComplexFilterCriteria() {}
    ComplexFilterCriteria(unsigned int dbKey, unsigned int qTotalAlnLen, unsigned int tTotalAlnLen, double qTM, double tTM, double qChainTm, double tChainTm, const Complex &qComplex, const Complex &tComplex, Complex &qnewComplex, Complex &tnewComplex) :
                        dbKey(dbKey), qTotalAlnLen(qTotalAlnLen), tTotalAlnLen(tTotalAlnLen), qTM(qTM), tTM(tTM), qComplex(qComplex), tComplex(tComplex),qnewComplex(qnewComplex), tnewComplex(tnewComplex) {
                            alignedQChainTmScores.push_back(qChainTm);
                            alignedTChainTmScores.push_back(tChainTm);
                        }
    ~ComplexFilterCriteria() {
        alignedQChainTmScores.clear();
        alignedTChainTmScores.clear();
    }

    bool hasTM(float TMThr, int covMode, int filterMode){
        switch (filterMode){
            case LocalParameters::FILTER_MODE_INTERFACE:
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
                }
        }
    }

    bool hasChainNum(int covMode, int filterMode, int qChainNum, int tChainNum ){
        switch (filterMode){
            case LocalParameters::FILTER_MODE_INTERFACE:
                switch (covMode) {
                    case Parameters::COV_MODE_BIDIRECTIONAL:
                        return (alignedQChainTmScores.size()==qChainNum && qChainNum==tChainNum);
                    case Parameters::COV_MODE_TARGET:
                        return (alignedTChainTmScores.size()==tChainNum);
                    case Parameters::COV_MODE_QUERY:
                        return (alignedQChainTmScores.size()==qChainNum);
                    case Parameters::COV_MODE_LENGTH_QUERY :
                    case Parameters::COV_MODE_LENGTH_TARGET :
                    case Parameters::COV_MODE_LENGTH_SHORTER :
                        return true;
                }
            case LocalParameters::FILTER_MODE_CONFORMATION:
                switch (covMode) {
                    case Parameters::COV_MODE_BIDIRECTIONAL:
                        return (qChainNum==tChainNum);
                    default:
                        return true;
                }
            case LocalParameters::FILTER_MODE_LOOSE:
                return true;

        }
    } 

    // bool hasMatchedCoord(int filterMode){
    //     switch (filterMode) {
    //         case LocalParameters::FILTER_MODE_INTERFACE:
    //             return (sameCoord);
    //         case LocalParameters::FILTER_MODE_CONFORMATION:
    //         case LocalParameters::FILTER_MODE_LOOSE:
    //             return true;
    //     }
    // }

    bool satisfy(int covMode, int filterMode, float covThr, float TMThr, float chainTMThr, int qChainNum, int tChainNum ) {
        //TODO
        // add filtermode in hasTM.
        // find other criteria for other filtermodes
        const bool covOK = Util::hasCoverage(covThr, covMode, qCov, tCov);
        const bool TMOK = hasTM(TMThr, covMode, filterMode);
        const bool chainNumOK = hasChainNum(covMode, filterMode, qChainNum, tChainNum);
        // const bool coordOK = hasMatchedCoord(filterMode);
        const bool chainTMOK = hasChainTm(chainTMThr, covMode, alignedQChainTmScores, alignedTChainTmScores, qChainNum, tChainNum);
        // const bool numOK = hasChainnum(sameChainNum, qChainNum, tChainNum);
        // return (covOK && TMOK && chainTMOK && numOK);
        return (covOK && TMOK && chainNumOK && chainTMOK);
    }

    void update(unsigned int qTotalAlnLen, unsigned int tTotalAlnLen, double qChainTm, double tChainTm) {
        this->qTotalAlnLen += qTotalAlnLen;
        this->tTotalAlnLen += tTotalAlnLen;
        this->alignedQChainTmScores.push_back(qChainTm);
        this->alignedTChainTmScores.push_back(tChainTm);
    }

    void calcCov(unsigned int qLen, unsigned int tLen) {
        qCov = static_cast<float>(qTotalAlnLen) / static_cast<float>(qLen);
        tCov = static_cast<float>(tTotalAlnLen) / static_cast<float>(tLen);
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

// static void getInterfaceIndex(float* &data1, int size1, float* &data2, int size2, int interfaceDistThr, std::vector<float> &x, std::vector<float> &y, std::vector<float> &z){
//     unsigned int interLen = 0;
//     for (int i =0; i < size1; i ++){
//         for (int j =0; j < size2; j ++){
//             float atomDist = BasicFunction::dist(data1[i], data1[i+size1], data1[i+size1*2], data2[j], data2[j+size2], data2[j+size2*2]);
//             if (atomDist<= interfaceDistThr){
//                 x.push_back(data1[i]);
//                 y.push_back(data1[i+size1]);
//                 z.push_back(data1[i+size1*2]);
//                 x.push_back(data2[j]);
//                 y.push_back(data2[j+size2]);
//                 z.push_back(data2[j+size2*2]);
//                 interLen ++;
//             }
//         }
//     }
//     Coordinates 
// }

// static void getInterface(std::vector<unsigned int> &ChainKeys, DBReader<unsigned int> &StructDbr, std::vector<Coordinates> &interfaceIndexVec, IndexReader* &DBr, unsigned int thread_idx, int interfaceDistThr){
//     for (size_t i = 0; i < ChainKeys.size(); ++i) {
//         unsigned int ChainDbId1 = DBr->sequenceReader->getId(ChainKeys[i]);
//         char *cadata1 = StructDbr.getData(ChainDbId1, thread_idx);
//         size_t CaLength1 = StructDbr.getEntryLen(ChainDbId1);
//         size_t seqLength1 = StructDbr.getSeqLen(ChainDbId1);
//         Coordinate16 coords1;
//         float* data1 = coords1.read(cadata1, seqLength1, CaLength1);
//         for (size_t j = i + 1; j < ChainKeys.size(); ++j) {
//             unsigned int ChainDbId2 = DBr->sequenceReader->getId(ChainKeys[j]);
//             char *cadata2 = StructDbr.getData(ChainDbId2, thread_idx);
//             size_t CaLength2 = StructDbr.getEntryLen(ChainDbId2);
//             size_t seqLength2 = StructDbr.getSeqLen(ChainDbId2);
//             Coordinate16 coords2;
//             float* data2 = coords2.read(cadata2, seqLength2, CaLength2);
//             std::vector<float> x;
//             std::vector<float> y;
//             std::vector<float> z;
//             getInterfaceIndex(data1, seqLength1, data2, seqLength2, interfaceDistThr, x, y, z);
//             interfaceIndexVec.push_back(interfaceCoord);
//         }
//     }
// }

int filtercomplex(int argc, const char **argv, const Command &command) {
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
    //TODO: remove resultWrite5 when done
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
        //TODO: set threshold for interface. now 15A
        int interfaceDistThr = 15;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::string result;
        std::map<unsigned int, ComplexFilterCriteria> localComplexMap;
        std::vector<unsigned int> assIdsToDelete;
        std::map<unsigned int, std::vector<unsigned int>> cmplIdToBestAssId; // cmplId : [assId, alnSum]
        std::vector<unsigned int> selectedAssIDs;
        Matcher::result_t res;
        
#pragma omp for schedule(dynamic, 1) 
        for (size_t queryComplexIdx = 0; queryComplexIdx < qComplexIdVec.size(); queryComplexIdx++) {
            Coordinate16 qcoords;
            Coordinate16 tcoords;   
            progress.updateProgress();
            unsigned int qComplexId = qComplexIdVec[queryComplexIdx];
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeyMap.at(qComplexId);
            std::map<unsigned int,std::map<unsigned int, unsigned int>> assIDtoqChainIdtotChainId;

            // std::vector<std::vector<int>> qInterfaceIndexVec;
            // getInterface(qChainKeys, qStructDbr, qInterfaceIndexVec, qDbr, thread_idx, interfaceDistThr);
            Complex qComplex;
            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++ ){
                unsigned int qChainKey = qChainKeys[qChainIdx];
                unsigned int qChainAlnId = alnDbr.getId(qChainKey);
                unsigned int qChainDbId = qDbr->sequenceReader->getId(qChainKey);
                ChainForInterface qchain(qChainDbId);

                //TODO: if monomer
                if (qChainAlnId == NOT_AVAILABLE_CHAIN_KEY){
                    break;
                }

                char *qcadata = qStructDbr.getData(qChainDbId, thread_idx);
                size_t qCaLength = qStructDbr.getEntryLen(qChainDbId);
                float* qdata = qcoords.read(qcadata, res.qLen, qCaLength);
                for (int pos = 0; pos <res.qLen; pos++){
                    Calpha calpha(pos, qdata[pos], qdata[res.qLen + pos], qdata[2*res.qLen + pos]);
                    qchain.Calphas.push_back(calpha);
                }
                qComplex.addChain(qchain);
            }
            qComplex.findInterface(interfaceDistThr);
            #pragma omp critical 
            {
                for (Interface interface : qComplex.InterfaceVec){
                    Debug(Debug::WARNING) << qcomplexIdToName.at(qComplexId) <<"\t" << interface.chainsinInterface.size() <<"\t"<<interface.chainsinInterface[0].chainid <<"\t" << interface.chainsinInterface[1].chainid <<"\n";
                }
            }
            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++ ) {
                unsigned int qChainKey = qChainKeys[qChainIdx];
                unsigned int qChainAlnId = alnDbr.getId(qChainKey);
                unsigned int qChainDbId = qDbr->sequenceReader->getId(qChainKey);
                char *data = alnDbr.getData(qChainAlnId, thread_idx);
                while (*data != '\0' ) {
                    ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                    data = Util::skipLine(data);
                    unsigned int assId = retComplex.assId;
                    unsigned int tChainKey= res.dbKey;
                    unsigned int tChainAlnId = alnDbr.getId(tChainKey);
                    unsigned int tChainDbId = tDbr->sequenceReader->getId(tChainKey);
                    assIDtoqChainIdtotChainId[assId][qChainDbId] = tChainDbId ;
                }
            }
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
                    char *qcadata = qStructDbr.getData(qChainDbId, thread_idx);
                    size_t qCaLength = qStructDbr.getEntryLen(qChainDbId);
                    float* qdata = qcoords.read(qcadata, res.qLen, qCaLength);
                
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
                    unsigned int tComplexId = tChainKeyToComplexIdMap.at(tChainKey);
                    std::vector<unsigned int> tChainKeys = tComplexIdToChainKeyMap.at(tComplexId);
                
                    float u[3][3];
                    float t[3];
                    fillUArr(retComplex.uString, u);
                    fillTArr(retComplex.tString, t);
                    char *tcadata = tStructDbr->getData(tChainDbId, thread_idx);
                    size_t tCaLength = tStructDbr->getEntryLen(tChainDbId);
                    float* tdata = tcoords.read(tcadata, res.dbLen, tCaLength);
                    unsigned int alnLen = cigarToAlignedLength(res.backtrace);
                    Coordinates qm(alnLen), tm(alnLen);
                    fillMatchedCoord(qdata, tdata, qm, tm, res.backtrace, res.qStartPos, res.dbStartPos, res.qLen, res.dbLen);
                    double chainTm = computeChainTmScore(qm, tm, t, u, alnLen, res.dbLen);
                    double qChainTm = chainTm / res.qLen;
                    double tChainTm = chainTm/ res.dbLen;
                    unsigned int qtotalaln = (std::max(res.qStartPos, res.qEndPos) - std::min(res.qStartPos, res.qEndPos) + 1);
                    unsigned int ttotalaln = (std::max(res.dbStartPos, res.dbEndPos) - std::min(res.dbStartPos, res.dbEndPos) + 1);
                    if (localComplexMap.find(assId) == localComplexMap.end()) {
                        Complex tComplex;
                        for (size_t tChainIdx = 0; tChainIdx < tChainKeys.size(); tChainIdx++ ){
                            unsigned int tChainKey = tChainKeys[tChainIdx];
                            unsigned int tChainAlnId = alnDbr.getId(tChainKey);
                            unsigned int tChainDbId = tDbr->sequenceReader->getId(tChainKey);
                            ChainForInterface tchain(tChainDbId);

                            //TODO: if monomer
                            if (tChainAlnId == NOT_AVAILABLE_CHAIN_KEY){
                                break;
                            }

                            char *tcadata = tStructDbr->getData(tChainDbId, thread_idx);
                            size_t tCaLength = tStructDbr->getEntryLen(tChainDbId);
                            float* tdata = tcoords.read(tcadata, res.dbLen, tCaLength);

                            for (int pos = 0; pos <res.dbLen; pos++){
                                Calpha calpha(pos, tdata[pos], tdata[res.dbLen + pos], tdata[2*res.dbLen + pos]);
                                tchain.Calphas.push_back(calpha);
                            }
                            tComplex.addChain(tchain);
                        }
                        tComplex.findInterface(interfaceDistThr);
                        Complex qnewComplex;
                        Complex tnewComplex;
                        unsigned int maxInterNum = max(tComplex.InterfaceVec.size(), qComplex.InterfaceVec.size()) ;
                        for(size_t intnum=0 ; intnum <= maxInterNum; intnum ++){
                            Interface inter;
                            qnewComplex.InterfaceVec.push_back(inter);
                            tnewComplex.InterfaceVec.push_back(inter);
                        }
                        std::map<unsigned int, std::map<unsigned int, unsigned int>> assIdtoqIntpostotIntpos;
                        for (size_t qIntpos = 0; qIntpos < qComplex.InterfaceVec.size(); qIntpos++){
                            for (size_t tIntpos = 0; tIntpos < tComplex.InterfaceVec.size(); tIntpos++){
                                if (tComplex.InterfaceVec[tIntpos].chainsinInterface[0].chainid == assIDtoqChainIdtotChainId[assId][qComplex.InterfaceVec[qIntpos].chainsinInterface[0].chainid] && tComplex.InterfaceVec[tIntpos].chainsinInterface[1].chainid == assIDtoqChainIdtotChainId[assId][qComplex.InterfaceVec[qIntpos].chainsinInterface[1].chainid]){
                                    assIdtoqIntpostotIntpos[assId][qIntpos] = tIntpos ;
                                }
                                else if (tComplex.InterfaceVec[tIntpos].chainsinInterface[1].chainid == assIDtoqChainIdtotChainId[assId][qComplex.InterfaceVec[qIntpos].chainsinInterface[0].chainid] && tComplex.InterfaceVec[tIntpos].chainsinInterface[0].chainid == assIDtoqChainIdtotChainId[assId][qComplex.InterfaceVec[qIntpos].chainsinInterface[1].chainid]) {
                                    assIdtoqIntpostotIntpos[assId][qIntpos] = tIntpos ;
                                }
                            }
                        }
                        ComplexFilterCriteria cmplfiltcrit = ComplexFilterCriteria(tChainKey, qtotalaln, ttotalaln, retComplex.qTmScore, retComplex.tTmScore, qChainTm, tChainTm, qComplex, tComplex, qnewComplex, tnewComplex);
                        localComplexMap[assId] = cmplfiltcrit;
                    } else {
                        localComplexMap.at(assId).update(qtotalaln, ttotalaln, qChainTm, tChainTm);
                    }
                } // while end
            }
            for (auto& assId_res : localComplexMap){
                unsigned int tComplexId = tChainKeyToComplexIdMap.at(assId_res.second.dbKey);
                assId_res.second.calcCov(qComplexLength.at(qComplexId), tComplexLength.at(tComplexId));
                std::vector<unsigned int> tChainKeys = tComplexIdToChainKeyMap.at(tComplexId);
                if (!assId_res.second.satisfy(par.covMode, par.filterMode, par.covThr, par.filtComplexTmThr, par.filtChainTmThr, qChainKeys.size(), tChainKeys.size())){
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
    qStructDbr.close();
    alnDbr.close();
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
