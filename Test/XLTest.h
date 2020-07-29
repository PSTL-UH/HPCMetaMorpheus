#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include "tangible_filesystem.h"

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class CommonParameters; }
namespace TaskLayer { class XlSearchParameters; }

#include "Chemistry/Chemistry.h"
using namespace Chemistry;

#include "EngineLayer/EngineLayer.h"
using namespace EngineLayer;
using namespace EngineLayer::CrosslinkSearch;
using namespace EngineLayer::Indexing;

#include "MassSpectrometry/MassSpectrometry.h"
using namespace MassSpectrometry;

#include "MzLibUtil.h"
using namespace MzLibUtil;
//using namespace Nett;
//using namespace NUnit::Framework;

#include "Proteomics/Proteomics.h"
using namespace Proteomics;
using namespace Proteomics::AminoAcidPolymer;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

#include "TaskLayer/TaskLayer.h"
using namespace TaskLayer;

#include "UsefulProteomicsDatabases/UsefulProteomicsDatabases.h"
using namespace UsefulProteomicsDatabases;

namespace Test
{
    class XLTest final
    {
    private:
        static IndexingResults *privateindexResults;
        static CommonParameters *privatecommonParameters;
        static XlSearchParameters *privatexlSearchParameters;
        static std::vector<Protein*> privateproteinList;
        static std::vector<Modification*> privatevariableModifications;
        static std::vector<Modification*> privatefixedModifications;
        static Crosslinker *privatecrosslinker;
        static std::vector<PeptideWithSetModifications*> privatedigestedList;
        
        static IndexingResults *getindexResults();
        static void setindexResults(IndexingResults *value);
        static CommonParameters *getcommonParameters();
        static void setcommonParameters(CommonParameters *value);
        static XlSearchParameters *getxlSearchParameters();
        static void setxlSearchParameters(XlSearchParameters *value);
        static std::vector<Protein*> getproteinList();
        static void setproteinList(const std::vector<Protein*> &value);
        static std::vector<Modification*> getvariableModifications();
        static void setvariableModifications(const std::vector<Modification*> &value);
        static std::vector<Modification*> getfixedModifications();
        static void setfixedModifications(const std::vector<Modification*> &value);
        static Crosslinker *getcrosslinker();
        static void setcrosslinker(Crosslinker *value);
        static std::vector<PeptideWithSetModifications*> getdigestedList();
        static void setdigestedList(const std::vector<PeptideWithSetModifications*> &value);
        
    public:
        static void XlTestXlPosCal();
        
        static void XlTestGenerateIntensityRanks();
        
        static void XlTest_BSA_DSSO();
        
        static void XlTest_BSA_DSS_file();
        
        static void XlTest_GenerateUserDefinedCrosslinker();
        
        static void XlTest_DiffCrosslinkSites();
        
        /// <summary>
        /// Verifies that crosslinker is generated properly
        /// </summary>
        static void CrosslinkCreateTest();
        
        static void DeadendPeptideTest();
        
        /// <summary>
        /// Makes sure helper methods that generate indices function properly
        /// </summary>
        static void XLSearchWithGeneratedIndices();
        
        /// <summary>
        /// Tests that a crosslinker that links at proline ("P") generates the correct indices
        /// of potential crosslink sites in the sequence PEPTIDE. The indices should be positions 1 and 3
        /// </summary>
        static void TestGetPossibleCrosslinkerSites();
        
        /// <summary>
        /// Generate and test fragments for this loop-linked peptide:
        ///    _
        ///   | |
        /// PEPTIDE
        /// </summary>
        static void TestTheoreticalFragmentsLoop();
        
        /// <summary>
        /// Generate and test fragments for loop-linked peptides with a mod placed at several different locations
        /// </summary>
        static void TestTheoreticalLoopFragmentsWithMod();
        
        /// <summary>
        /// Generate and test fragments for this dead-end peptide (quenched with tris):
        ///   
        ///   |
        /// PEPTIDE[Methyl]
        /// </summary>
        static void TestDeadendTris();
        
        /// <summary>
        /// Generate and test fragments for this crosslinked peptide pair that was crosslinked
        /// with a non-cleavable crosslinker:
        ///  PRLTEIN
        ///   |
        /// PEPTIDE
        /// </summary>
        static void TestTheoreticalFragmentsNonCleavableCrosslink();
        
        /// <summary>
        /// Generate and test fragments for this crosslinked peptide pair that was crosslinked
        /// with a cleavable crosslinker:
        ///   PROTEIN
        ///   |
        /// PEPTIDE
        /// </summary>
        static void TestTheoreticalFragmentsCleavableCrosslink();
        
        static void TestWriteToPercolator();
        
        static void TestWriteNonSingleCross();
    };
    
    class XLTestDataFile : public MsDataFile
    {
        //Create DSSO crosslinked fake MS data. Include single, deadend, loop, inter, intra crosslinks ms2 data for match.
    public:
        XLTestDataFile();
        
        std::string getFilePath() const;
        
        std::string getName() const;
        
        void ReplaceFirstScanArrays(std::vector<double> &mz, std::vector<double> &intensities);
    };
    
    class XLTestDataFileDiffSite : public MsDataFile
    {
        //Create DSSO crosslinked fake MS data. Include single, deadend, loop, inter, intra crosslinks ms2 data for match.
    public:
        XLTestDataFileDiffSite();
        
        std::string getFilePath() const;
        
        std::string getName() const;
        
        void ReplaceFirstScanArrays(std::vector<double> &mz, std::vector<double> &intensities);
    };
}
