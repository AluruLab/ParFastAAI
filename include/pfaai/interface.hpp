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
    virtual inline int closeDB() = 0;
    //
    virtual int
    queryGenomeTetramers(const std::string protein, IdType tetramerStart,
                         IdType tetramerEnd,
                         std::vector<IdType>& Lc) const = 0;  // NOLINT
    virtual int
    queryProtienSetGPPairs(const std::vector<std::string>& proteinSet,
                           IdType tetramerStart, IdType tetramerEnd,
                           std::vector<IdType>& Lp,                // NOLINT
                           std::vector<IdPairType>& F) const = 0;  // NOLINT
    virtual int
    queryMetaData(std::vector<std::string>& proteinSet,            // NOLINT
                  std::vector<std::string>& genomeSet) const = 0;  // NOLINT
    virtual int
    queryProtienTetramerCounts(const std::vector<std::string>& proteinSet,
                               IdType proteinStart, IdType proteinEnd,
                               IdMatType& T) const = 0;  // NOLINT
    //
    virtual ~DataBaseInterface() {}
};

template <typename IdType, typename IdPairType, typename IdMatrixType,
          typename JACType>
class DataStructInterface {
  public:
    // constants
    constexpr static IdType NTETRAMERS = (20 * 20 * 20 * 20);
    constexpr static float DEFAULT_SLACK_PCT = 0.0;
    //
    DataStructInterface() {}
    //
    virtual const std::vector<IdType>& getLc() const = 0;
    virtual const std::vector<IdType>& getLp() const = 0;
    virtual const IdMatrixType& getT() const = 0;
    virtual const std::vector<IdPairType>& getF() const = 0;
    virtual float getSlackPercentage() const = 0;
    virtual inline IdType getTetramerCount() const { return NTETRAMERS; }
    virtual IdType getQryGenomeCount() const = 0;
    virtual IdType getTgtGenomeCount() const = 0;
    virtual IdType getGPCount() const = 0;
    virtual IdType genomePairToJACIndex(IdType genomeA,
                                        IdType genomeB) const = 0;
    virtual void initJAC(std::vector<JACType>& jac_tuples) const = 0;  // NOLINT
    //
    virtual ~DataStructInterface() {}
};

#endif  // !PFAAI_INTERFACE_HPP
