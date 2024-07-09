#ifndef PFAAI_INTERFACE_HPP
#define PFAAI_INTERFACE_HPP

#include <string>
#include <vector>

enum PFAAI_ERROR_CODE {
    PFAAI_OK = 0,
    PFAAI_ERR_SQLITE_DB = 1,
    PFAAI_ERR_SQLITE_MEM_ALLOC = 2,
    PFAAI_ERR_CONSTRUCT = 3
};

template <typename IdType, typename IdPairType, typename IdMatType,
          typename ErrCodeType = int>
class DataBaseInterface {

  public:
    DataBaseInterface() {}
    //
    virtual inline bool isDBOpen() const = 0;
    virtual ErrCodeType getDBErrorCode() const = 0;
    virtual const char* getDBError() const = 0;
    virtual std::string getDBPath() const = 0;
    virtual PFAAI_ERROR_CODE validate() const = 0;
    //
    virtual int queryGenomeTetramers(const std::string protein,
                                     IdType tetramerStart, IdType tetramerEnd,
                                     std::vector<IdType>& Lc) = 0;
    virtual int
    queryProtienSetGPPairs(const std::vector<std::string>& proteinSet,
                           IdType tetramerStart, IdType tetramerEnd,
                           std::vector<IdType>& Lp,
                           std::vector<IdPairType>& F) = 0;
    virtual int queryMetaData(std::vector<std::string>& proteinSet,
                              std::vector<std::string>& genomeSet) = 0;
    virtual int
    queryProtienTetramerCounts(const std::vector<std::string>& proteinSet,
                               IdType proteinStart, IdType proteinEnd,
                               IdMatType& T) = 0;
    //
    virtual inline int closeDB() = 0;
    virtual ~DataBaseInterface() {}
};

template <typename IdType, typename IdPairType, typename IdMatrixType>
class DataStructInterface {
  public:
    DataStructInterface() {}
    //
    virtual const std::vector<IdType>& getLc() const = 0; 
    virtual const std::vector<IdType>& getLp() const = 0;
    virtual const IdMatrixType& getT() const = 0;
    virtual const std::vector<IdPairType>& getF() const  = 0;
    virtual float getSlackPercentage() const = 0;
    virtual IdType getTetramerCount() const = 0;
    virtual IdType getGenomeCount() const = 0;
    virtual IdType getGPCount() const = 0;
    virtual  IdType genomePairToJACIndex(IdType genomeA, IdType genomeB) const = 0;
    //
    virtual ~DataStructInterface() {}
};

#endif  // !PFAAI_INTERFACE_HPP
