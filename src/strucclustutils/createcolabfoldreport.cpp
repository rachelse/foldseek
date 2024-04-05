/*
 * File: createcolabfoldreport.cpp
 * Project: strucclustutils
 * File Created: 5th Apr 2024
 * Author: Rachel Seongeun Kim (seamustard52@gmail.com)
 * -----
 * Copyright: Rachel Seongeun Kim
 */

#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "IndexReader.h"
#include "FileUtil.h"
#include "TranslateNucl.h"
#include "MemoryMapped.h"
#include "createcomplexreport.h"
#include "LDDT.h"
#include "CalcProbTP.h"
#include "Coordinate16.h"
#include "tmalign/Coordinates.h"
#include "tmalign/basic_fun.h"
#include <map>
#ifdef OPENMP
#include <omp.h>
#endif

double computeSimpleTm(Coordinates &qm, Coordinates &tm, float t[3], float u[3][3], unsigned int mlen, int normlen) {
    double tmscore = 0;
    float d0;
    // float score_d8 = 1.5*pow(normlen,0.3)+3.5;
    
    if (normlen<=19) {
        d0=0.168;
    }
    else {
        d0=1.24*pow((normlen-15),1.0/3)-1.8;
    }
    d0 += 0.8;

    Coordinates tmt(mlen);
    BasicFunction::do_rotation(tm, tmt, mlen, t, u);

    float d02 = d0*d0;
    // float score_d82 = score_d8*score_d8;
    for (unsigned int k=0; k<mlen; k++) {
        double di = BasicFunction::dist(qm.x[k], qm.y[k], qm.z[k], tmt.x[k], tmt.y[k], tmt.z[k]);
        // if (di < score_d82) {
        //     tmscore += 1/(1+di/d02);
        // }
        tmscore += 1/(1+di/d02);
    }
    return tmscore;
}

void fillArrU(const std::string &uString, float (&u)[3][3]) {
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

void fillArrT(const std::string &tString, float (&t)[3]) {
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

// void getComplexNameChainName(const chainName_t &chainName, compNameChainName_t &compAndChainName) {
//     size_t pos = chainName.rfind('_');
//     std::string comp = chainName.substr(0, pos);
//     std::string chain = chainName.substr(pos + 1);
//     compAndChainName = {comp, chain};
// }

// void getScoreComplexResults(
//         std::vector<ScoreComplexResult> &scoreComplexResults,
//         const std::vector<std::string> &qChainVector,
//         const std::vector<std::string> &tChainVector,
//         double qTMScore,
//         double tTMScore,
//         const std::string &u,
//         const std::string &t,
//         unsigned int assId
// ) {
//     char buffer[1024];
//     resultToWrite_t resultToWrite;
//     std::string qComplexName;
//     std::string tComplexName;
//     std::string qChainString;
//     std::string tChainString;
//     compNameChainName_t compAndChainName;
//     getComplexNameChainName(qChainVector[0], compAndChainName);
//     qComplexName = compAndChainName.first;
//     qChainString = compAndChainName.second;
//     getComplexNameChainName(tChainVector[0], compAndChainName);
//     tComplexName = compAndChainName.first;
//     tChainString = compAndChainName.second;
//     for (size_t qChainId = 1; qChainId < qChainVector.size(); qChainId++) {
//         getComplexNameChainName(qChainVector[qChainId], compAndChainName);
//         qChainString += ',' + compAndChainName.second;
//     }
//     for (size_t tChainId = 1; tChainId < tChainVector.size(); tChainId++) {
//         getComplexNameChainName(tChainVector[tChainId], compAndChainName);
//         tChainString += ',' + compAndChainName.second;
//     }
//     int count = snprintf(buffer,sizeof(buffer),"%s\t%s\t%s\t%s\t%1.5f\t%1.5f\t%s\t%s\t%d\n", qComplexName.c_str(), tComplexName.c_str(), qChainString.c_str(), tChainString.c_str(), qTMScore, tTMScore, u.c_str(), t.c_str(), assId);
//     resultToWrite.append(buffer, count);
//     scoreComplexResults.emplace_back(assId, resultToWrite);
// }

struct ColabfoldAlignment {
    ColabfoldAlignment(){};
    ColabfoldAlignment(chainName_t &qChainName, chainName_t &tChainName, double qTmScore, double tTmScore, std::string &u, std::string &t,  unsigned int assId) : qTMScore(qTmScore), tTMScore(tTmScore), u(u), t(t), assId(assId){
        qChainNames = {qChainName};
        tChainNames = {tChainName};
    };
    std::vector<chainName_t> qChainNames;
    std::vector<chainName_t> tChainNames;
    double qTMScore;
    double tTMScore;
    std::string u;
    std::string t;
    unsigned int assId;
    std::vector<float> qChainTms;
    std::vector<float> tChainTms;
};

int createcolabfoldreport(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    const bool sameDB = par.db1.compare(par.db2) == 0 ? true : false;
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    int dbaccessMode = (DBReader<unsigned int>::USE_INDEX);
    std::map<unsigned int, unsigned int> qKeyToSet;
    std::map<unsigned int, unsigned int> tKeyToSet;
    IndexReader qDbr(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    IndexReader qDbrHeader(par.db1, par.threads, IndexReader::SRC_HEADERS , (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    IndexReader *tDbrHeader;
    if (sameDB) {
        tDbrHeader = &qDbrHeader;
    } else {
        tDbrHeader = new IndexReader(par.db2, par.threads, IndexReader::SRC_HEADERS, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    size_t localThreads = 1;
#ifdef OPENMP
    localThreads = std::max(std::min((size_t)par.threads, alnDbr.getSize()), (size_t)1);
#endif

    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), 1, shouldCompress, dbType);
    resultWriter.open();
    const bool isDb = par.dbOut;
    std::string qLookupFile = par.db1 + ".lookup";
    TranslateNucl translateNucl(static_cast<TranslateNucl::GenCode>(par.translationTable));

    Matcher::result_t res;
    std::map<unsigned int, unsigned int> qChainKeyToComplexIdMap;
    std::map<unsigned int, std::vector<unsigned int>> qComplexIdToChainKeyMap;
    std::vector<unsigned int> qComplexIdVec;
    getKeyToIdMapIdToKeysMapIdVec(qLookupFile, qChainKeyToComplexIdMap, qComplexIdToChainKeyMap, qComplexIdVec);
    qChainKeyToComplexIdMap.clear();
    Debug::Progress progress(qComplexIdVec.size());

    std::vector<ScoreComplexResult> complexResults;
#pragma omp parallel num_threads(localThreads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Matcher::result_t res;
        std::vector<ScoreComplexResult> localComplexResults;
#pragma omp for schedule(dynamic, 10) nowait
        for (size_t queryComplexIdx = 0; queryComplexIdx < qComplexIdVec.size(); queryComplexIdx++) {
            progress.updateProgress();
            std::vector<unsigned int> assIdVec;
            std::vector<ColabfoldAlignment> compAlns;
            unsigned int qComplexId = qComplexIdVec[queryComplexIdx];

            Coordinate16 qcoords;
            Coordinate16 tcoords;
            std::vector<unsigned int> &qChainKeys = qComplexIdToChainKeyMap[qComplexId];
            for (size_t qChainIdx = 0; qChainIdx < qChainKeys.size(); qChainIdx++ ) {
                unsigned int qChainKey = qChainKeys[qChainIdx];
                unsigned int qChainDbKey = alnDbr.getId(qChainKey);
                if (qChainDbKey == NOT_AVAILABLE_CHAIN_KEY) {
                    continue;
                }
                size_t qHeaderId = qDbrHeader.sequenceReader->getId(qChainKey);
                const char *qHeader = qDbrHeader.sequenceReader->getData(qHeaderId, thread_idx);
                compNameChainName_t qCompAndChainName;
                chainName_t queryChainName = Util::parseFastaHeader(qHeader);
                int qChainLen = qDbr->sequenceReader->getSeqLen(qChainDbKey);
                char *qcadata = qStructDbr.getData(qChainDbKey, thread_idx);
                size_t qCaLength = qStructDbr.getEntryLen(qChainDbKey);
                float* qdata = qcoords.read(qcadata, qChainLen, qCaLength);
                // getComplexNameChainName(queryChainName, qCompAndChainName);
                char *data = alnDbr.getData(qChainDbKey, thread_idx);
                while (*data != '\0') {
                    ComplexDataHandler retComplex = parseScoreComplexResult(data, res);
                    if (!retComplex.isValid){
                        Debug(Debug::ERROR) << "No scorecomplex result provided";
                        EXIT(EXIT_FAILURE);
                    }
                    data = Util::skipLine(data);
                    size_t tHeaderId = tDbrHeader->sequenceReader->getId(res.dbKey);
                    const char *tHeader = tDbrHeader->sequenceReader->getData(tHeaderId, thread_idx);
                    chainName_t targetChainName = Util::parseFastaHeader(tHeader);
                    unsigned int assId = retComplex.assId;
                    unsigned int tChainDbKey = res.dbKey;
                    unsigned int compAlnIdx = std::find(assIdVec.begin(), assIdVec.end(), assId) - assIdVec.begin();

                    float u[3][3];
                    float t[3];
                    Coordinates qm(0), tm(0);
                    fillArrU(retComplex.uString, u);
                    fillArrT(retComplex.tString, t);

                    int tChainLen = tDbr->sequenceReader->getSeqLen(tChainDbKey);
                    char *tcadata = tStructDbr->getData(tChainDbKey, thread_idx);
                    size_t tCaLength = tStructDbr->getEntryLen(tChainDbKey);
                    float* tdata = tcoords.read(tcadata, tChainLen, tCaLength);
                    unsigned int normlen = std::min(res.qLen, res.dbLen);
                    double chainTm = computeSimpleTm(qm, tm, t, u, normlen, normlen);
                    double qTmScore = chainTm / qChainLen;
                    double tTmScore = chainTm / tChainLen;
                    if (compAlnIdx == compAlns.size()) {
                        assIdVec.emplace_back(assId);
                        compAlns.emplace_back(queryChainName, targetChainName, retComplex.qTmScore, retComplex.tTmScore, retComplex.uString, retComplex.tString, assId);
                        compAlns[compAlnIdx].qChainTms.emplace_back(qTmScore);
                        compAlns[compAlnIdx].tChainTms.emplace_back(tTmScore);
                    } else {
                        compAlns[compAlnIdx].qChainNames.emplace_back(queryChainName);
                        compAlns[compAlnIdx].tChainNames.emplace_back(targetChainName);
                        compAlns[compAlnIdx].qChainTms.emplace_back(qTmScore);
                        compAlns[compAlnIdx].tChainTms.emplace_back(tTmScore);
                    }

                } // while end
            }
            // for (size_t compAlnIdx = 0; compAlnIdx < compAlns.size(); compAlnIdx++) {
            //     const ColabfoldAlignment &aln = compAlns[compAlnIdx];
            //     getScoreComplexResults(localComplexResults, aln.qChainNames, aln.tChainNames, aln.qTMScore, aln.tTMScore, aln.u, aln.t, aln.assId);
            // }
        } // for end
#pragma omp critical
        {
            complexResults.insert(complexResults.end(), localComplexResults.begin(), localComplexResults.end());
        }
    } // MP end
    SORT_PARALLEL(complexResults.begin(), complexResults.end(), compareComplexResult);
    for (size_t complexResIdx = 0; complexResIdx < complexResults.size(); complexResIdx++) {
        const ScoreComplexResult& res = complexResults[complexResIdx];
        const resultToWrite_t& data = res.resultToWrite;
        resultWriter.writeData(data.c_str(), data.length(), res.assId, 0, isDb, isDb);
    }
    resultWriter.close(true);
    if (isDb == false) {
        FileUtil::remove(par.db4Index.c_str());
    }
    alnDbr.close();
    if (sameDB == false) {
        delete tDbrHeader;
    }
    return EXIT_SUCCESS;
}