#include "MzIdentMLWriter.h"
#include "../../EngineLayer/PeptideSpectralMatch.h"
#include "../../EngineLayer/ProteinParsimony/ProteinGroup.h"
#include "../../EngineLayer/GlobalVariables.h"
#include "Chemistry/ClassExtensions.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace MzLibUtil;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace TaskLayer
{

    // Having a std::tuple as the key of a std::unordered_map is not directly allowed
    // in C++.
    // See https://stackoverflow.com/questions/11408934/using-a-stdtuple-as-key-for-stdunordered-map
    // for how to handle this.
    
    typedef std::tuple<std::string, int> MzIdentTuple;
    
    struct MzIdentTuple_hash: public std::unary_function<MzIdentTuple, std::size_t>{
        std::size_t operator() (const MzIdentTuple& k ) const
            {
                size_t h1= std::hash<std::string>{}(std::get<0>(k));
                size_t h2= std::hash<int>{}(std::get<1>(k));
                return h1 ^ (h2 << 1);
            }
    };
    
    struct MzIdentTuple_equal: public std::binary_function<MzIdentTuple, MzIdentTuple, bool>{
        bool operator() (const MzIdentTuple& lhs, const MzIdentTuple& rhs) const
            {
                return std::get<0>(lhs) == std::get<0>(rhs) &&
                    std::get<1>(lhs) == std::get<1>(rhs);
            }
    };
    
    
    void MzIdentMLWriter::WriteMzIdentMl(std::vector<PeptideSpectralMatch*> &psms,
                                         std::vector<EngineLayer::ProteinGroup*> &groups,
                                         std::vector<Modification*> &variableMods,
                                         std::vector<Modification*> &fixedMods,
                                         std::vector<Protease*> &proteases,
                                         double qValueFilter,
                                         Tolerance *productTolerance,
                                         Tolerance *parentTolerance,
                                         int missedCleavages,
                                         const std::string &outputPath)
    {
        
#ifdef ORIG
        psms = psms.Where([&] (std::any p)   {
                return p::FdrInfo::QValue <= qValueFilter && p::FdrInfo::QValueNotch <= qValueFilter;
            });
#endif
        std::vector<PeptideSpectralMatch*> tmp;
        for ( auto p : psms ) {
            if ( p->getFdrInfo()->getQValue() <= qValueFilter &&
                 p->getFdrInfo()->getQValueNotch() < qValueFilter ) {
                tmp.push_back(p);
            }
        }
        psms.clear();
        for ( auto p: tmp ) {
            psms.push_back(p);
        }
        
#ifdef ORIG
        std::vector<PeptideWithSetModifications*> peptides = psms.SelectMany([&] (std::any i) {
                i::BestMatchingPeptides->Select([&] (std::any v)  {
                        v::Peptide;
                    });
            }).Distinct().ToList();
#endif
        std::vector<PeptideWithSetModifications*> peptides;
        for ( auto i : psms ) {
            for ( auto v:  i->getBestMatchingPeptides() ) {
                bool found = false;
                for ( auto k: peptides ) {
                    if ( k == std::get<1>(v) ) {
                        found = true;
                        break;
                    }
                }
                if ( !found ) {
                    peptides.push_back(std::get<1>(v));
                }
            }
        }
        
#ifdef ORIG
        std::vector<Protein*> proteins = peptides.Select([&] (std::any p)	{
                p::Protein;
            }).Distinct().ToList();
#endif
        std::vector<Protein*> proteins;
        for ( auto p: peptides ) {
            bool found = false;
            for (auto k : proteins ) {
                if ( p->getProtein()  == k ) {
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                proteins.push_back(p->getProtein());
            }
        }
        
#ifdef ORIG
        std::vector<std::string> filenames = psms.Select([&] (std::any i)   {
                i::FullFilePath;
            }).Distinct().ToList();
#endif
        std::vector<std::string> filenames;
        for ( auto i: psms ) {
            bool found = false;
            for ( auto k : filenames ) {
                if ( k == i->getFullFilePath() ) {
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                filenames.push_back(i->getFullFilePath());
            }
        }
        
        std::unordered_map<std::string, std::string> database_reference;
#ifdef ORIG
        std::vector<std::string> databases = proteins.Select([&] (std::any p)   {
                p::DatabaseFilePath;
            }).Distinct().ToList();
#endif
        std::vector<std::string> databases;
        for ( auto i: proteins ) {
            bool found = false;
            for ( auto k : databases ) {
                if ( k == i->getDatabaseFilePath() ) {
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                filenames.push_back(i->getDatabaseFilePath());
            }
        }
        
        //UTF8Encoding *utf8EmitBOM = new UTF8Encoding(false);
        //XmlWriterSettings *settings = new XmlWriterSettings();
        //settings->NewLineChars = "\n";
        //settings->Indent = true;
        //settings->Encoding = utf8EmitBOM;
        //XmlSerializer *_indexedSerializer = new XmlSerializer(mzIdentML110::Generated::MzIdentMLType110::typeid);
        //auto _mzid = new mzIdentML110::MzIdentMLType110();
        auto _mzid = new mzIdentML110::MzIdentMLType();
        _mzid->version() = "1.1.0";
        _mzid->id() = "";
        
        _mzid->Provider() = (* new mzIdentML110::ProviderType());
        _mzid->Provider()->id() = "PROVIDER";
        _mzid->Provider()->ContactRole() = (* new mzIdentML110::ContactRoleType());
        _mzid->Provider()->ContactRole().get().contact_ref() = "UWMadisonSmithGroup";
        _mzid->Provider()->ContactRole().get().Role() = *(new mzIdentML110::RoleType());
        _mzid->Provider()->ContactRole().get().Role().cvParam() = *(new mzIdentML110::CVParamType());
        _mzid->Provider()->ContactRole().get().Role().cvParam().accession() = "MS:1001271";
        _mzid->Provider()->ContactRole().get().Role().cvParam().name() = "researcher";
        _mzid->Provider()->ContactRole().get().Role().cvParam().cvRef() = "PSI-MS";
        
        //_mzid->AuditCollection().get() = std::vector<mzIdentML110::AbstractContactType>(2);
        
        mzIdentML110::PersonType *tempVar = new mzIdentML110::PersonType();
        tempVar->id() = "UWMadisonSmithGroupPerson";

        mzIdentML110::CVParamType *tempVar2 = new mzIdentML110::CVParamType();
        tempVar2->accession() = "MS:1000589";
        tempVar2->name() = "contact email";
        tempVar2->cvRef() = "PSI-MS";
        tempVar2->value() = "mm_support@chem.wisc.edu";

        mzIdentML110::CVParamType *tempVar3 = new mzIdentML110::CVParamType();
        tempVar3->accession() = "MS:1000590";
        tempVar3->name() = "affiliation name";
        tempVar3->cvRef() = "PSI-MS";
        tempVar3->value() = "UWMadisonSmithGroup";
        tempVar->cvParam() = {tempVar2, tempVar3};
        auto t = new mzIdentML110::AuditCollectionType::Person_sequence();
        t->push_back(*tempVar);
        _mzid->AuditCollection().get().Person() = *t;
        
        mzIdentML110::OrganizationType *tempVar4 = new mzIdentML110::OrganizationType();
        tempVar4->id() = "UWMadisonSmithGroup";

        mzIdentML110::CVParamType *tempVar5 = new mzIdentML110::CVParamType();
        tempVar5->accession() = "MS:1000589";
        tempVar5->name() = "contact email";
        tempVar5->cvRef() = "PSI-MS";
        tempVar5->value() = "mm_support@chem.wisc.edu";

        mzIdentML110::CVParamType *tempVar6 = new mzIdentML110::CVParamType();
        tempVar6->accession() = "MS:1000590";
        tempVar6->name() = "affiliation name";
        tempVar6->cvRef() = "PSI-MS";
        tempVar6->value() = "UWMadisonSmithGroup";
        tempVar4->cvParam() = {tempVar5, tempVar6};
        auto t4 =  new mzIdentML110::AuditCollectionType::Organization_sequence();
        t4->push_back(*tempVar4);
        _mzid->AuditCollection().get().Organization() = *t4;
        
        //cvlist: URLs of controlled vocabularies used within the file.
        mzIdentML110::cvType *tempVar7 = new mzIdentML110::cvType();
        tempVar7->id() = "PSI-MS";
        tempVar7->fullName() = "Proteomics Standards Initiative Mass Spectrometry Vocabularies";
        tempVar7->uri() = "https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo";
        tempVar7->version() = "4.0.9";

        mzIdentML110::cvType *tempVar8 = new mzIdentML110::cvType();
        tempVar8->id() = "PSI-MOD";
        tempVar8->fullName() = "Proteomics Standards Initiative Modification Vocabularies";
        tempVar8->uri() = "http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/mod/data/PSI-MOD.obo";
        tempVar8->version() = "1.2";

        mzIdentML110::cvType *tempVar9 = new mzIdentML110::cvType();
        tempVar9->id() = "UNIMOD";
        tempVar9->fullName() = "UNIT-ONTOLOGY";
        tempVar9->uri() = "http://www.unimod.org/obo/unimod.obo";

        mzIdentML110::cvType *tempVar10 = new mzIdentML110::cvType();
        tempVar10->id() = "UO";
        tempVar10->fullName() = "UNIT-ONTOLOGY";
        tempVar10->uri() = "http://www.unimod.org/obo/unimod.obo";
        auto t10 = new mzIdentML110::CVListType::cv_sequence();            
        t10->push_back(*tempVar7);
        t10->push_back(*tempVar8);
        t10->push_back(*tempVar9);
        t10->push_back(*tempVar10);
        _mzid->cvList().cv(*t10);
        
        mzIdentML110::AnalysisSoftwareType *tempVar11 = new mzIdentML110::AnalysisSoftwareType();
        tempVar11->id() = "AS_MetaMorpheus";
        tempVar11->name() = "MetaMorpheus";
        tempVar11->version() = GlobalVariables::getMetaMorpheusVersion();
        tempVar11->uri() = "https://github.com/smith-chem-wisc/MetaMorpheus";
        tempVar11->SoftwareName() = ( *new mzIdentML110::ParamType());
        tempVar11->SoftwareName().cvParam() = (*new mzIdentML110::CVParamType());
        tempVar11->SoftwareName().cvParam().get().accession() = "MS:1002826";
        tempVar11->SoftwareName().cvParam().get().name() = "MetaMorpheus";
        tempVar11->SoftwareName().cvParam().get().cvRef() = "PSI-MS";
        tempVar11->ContactRole() = (* new mzIdentML110::ContactRoleType());
        tempVar11->ContactRole().get().contact_ref() = "UWMadisonSmithGroup";
        tempVar11->ContactRole().get().Role() = (* new mzIdentML110::RoleType());
        tempVar11->ContactRole().get().Role().cvParam() = *new mzIdentML110::CVParamType();
        tempVar11->ContactRole().get().Role().cvParam().accession() = "MS:1001267";
        tempVar11->ContactRole().get().Role().cvParam().name() = "software vendor";
        tempVar11->ContactRole().get().Role().cvParam().cvRef() = "PSI-MS";
        auto t11 = new mzIdentML110::AnalysisSoftwareListType::AnalysisSoftware_sequence ();
        t11->push_back(*tempVar11);
        _mzid->AnalysisSoftwareList().get().AnalysisSoftware() = *t11;
        
        _mzid->DataCollection() = * new mzIdentML110::DataCollectionType();
        _mzid->DataCollection().AnalysisData() = *new mzIdentML110::AnalysisDataType();

        mzIdentML110::SpectrumIdentificationListType *tempVar12 = new mzIdentML110::SpectrumIdentificationListType();
        tempVar12->id() = "SI";
        tempVar12->SpectrumIdentificationResult() = mzIdentML110::SpectrumIdentificationListType::SpectrumIdentificationResult_sequence(psms.size());
        auto t12 = new mzIdentML110::AnalysisDataType::SpectrumIdentificationList_sequence();
        t12->push_back(*tempVar12);
        _mzid->DataCollection().AnalysisData().SpectrumIdentificationList() = *t12;
        
        _mzid->DataCollection().Inputs() = *new mzIdentML110::InputsType();
        _mzid->DataCollection().Inputs().SearchDatabase() = mzIdentML110::InputsType::SearchDatabase_sequence(databases.size());
        _mzid->DataCollection().Inputs().SpectraData() = mzIdentML110::InputsType::SpectraData_sequence(filenames.size());
        _mzid->SequenceCollection() = * new mzIdentML110::SequenceCollectionType();
        _mzid->SequenceCollection().get().Peptide() = mzIdentML110::SequenceCollectionType::Peptide_sequence (peptides.size() );
        _mzid->SequenceCollection().get().DBSequence() = mzIdentML110::SequenceCollectionType::DBSequence_sequence(proteins.size());
        _mzid->SequenceCollection().get().PeptideEvidence() = mzIdentML110::SequenceCollectionType::PeptideEvidence_sequence(peptides.size() );
        
        _mzid->AnalysisCollection() = *new mzIdentML110::AnalysisCollectionType();

        mzIdentML110::SpectrumIdentificationType *tempVar13 = new mzIdentML110::SpectrumIdentificationType();
        tempVar13->id() = "SI";
        tempVar13->spectrumIdentificationList_ref() = "SI";
        tempVar13->spectrumIdentificationProtocol_ref() = "SIP";
        tempVar13->InputSpectra() = mzIdentML110::SpectrumIdentificationType::InputSpectra_sequence(filenames.size());
        tempVar13->SearchDatabaseRef() = mzIdentML110::SpectrumIdentificationType::SearchDatabaseRef_sequence(databases.size());
        auto t13 = new mzIdentML110::AnalysisCollectionType::SpectrumIdentification_sequence();
        t13->push_back(*tempVar13);
        _mzid->AnalysisCollection().SpectrumIdentification() = *t13;
        
        int database_index = 0;
        for (auto database : databases)
        {
            mzIdentML110::SearchDatabaseType *tempVar14 = new mzIdentML110::SearchDatabaseType();
            tempVar14->id() = "SDB_" + std::to_string(database_index);
            tempVar14->location() = database;
            tempVar14->DatabaseName() = *new mzIdentML110::ParamType();
            tempVar14->DatabaseName().cvParam() = *new mzIdentML110::CVParamType();
            tempVar14->DatabaseName().cvParam().get().accession() = "MS:1001073";
            tempVar14->DatabaseName().cvParam().get().name() = "database type amino acid";
            tempVar14->DatabaseName().cvParam().get().cvRef() = "PSI-MS";
            _mzid->DataCollection().Inputs().SearchDatabase()[database_index] = *tempVar14;            
            database_reference.emplace(database, "SDB_" + std::to_string(database_index));
            
            mzIdentML110::SearchDatabaseRefType *tempVar15 = new mzIdentML110::SearchDatabaseRefType();
            tempVar15->searchDatabase_ref() = "SDB_" + std::to_string(database_index);
            _mzid->AnalysisCollection().SpectrumIdentification()[0].SearchDatabaseRef()[database_index] = *tempVar15;
            database_index++;
        }
            
        int protein_index = 0;
        for (auto protein : proteins)
        {
            mzIdentML110::DBSequenceType *tempVar16 = new mzIdentML110::DBSequenceType();
            tempVar16->id() = "DBS_" + protein->getAccession();
            //tempVar16->lengthSpecified() = true;
            //tempVar16->length() = protein->getLength();
            tempVar16->length().set(protein->getLength());
            tempVar16->searchDatabase_ref() = database_reference[protein->getDatabaseFilePath()];
            tempVar16->accession() = protein->getAccession();
            tempVar16->Seq() = protein->getBaseSequence();

            mzIdentML110::CVParamType *tempVar17 = new mzIdentML110::CVParamType();
            tempVar17->accession() = "MS:1001088";
            tempVar17->name() = "protein description";
            tempVar17->cvRef() = "PSI-MS";
            tempVar17->value() = protein->getFullDescription();
            auto t17 = new mzIdentML110::DBSequenceType::cvParam_sequence();
            t17->push_back(*tempVar17);
            tempVar16->cvParam() = *t17;
            
            tempVar16->name() = protein->getName();
            _mzid->SequenceCollection().get().DBSequence()[protein_index] = *tempVar16;
            protein_index++;
        }
        
        std::unordered_map<std::string, int> spectral_ids; //key is datafile, value is datafile's id
        int spectra_data_id = 0;
        for (auto data_filepath : filenames)
        {
            std::string extension = data_filepath.substr(data_filepath.find_last_of("."));
            bool thermoRawFile = extension == ".raw";
            std::string spectral_data_id = "SD_" + std::to_string(spectra_data_id);
            spectral_ids.emplace(data_filepath, spectra_data_id);

            mzIdentML110::InputSpectraType *tempVar18 = new mzIdentML110::InputSpectraType();
            tempVar18->spectraData_ref() = spectral_data_id;
            _mzid->AnalysisCollection().SpectrumIdentification()[0].InputSpectra()[spectra_data_id] = *tempVar18;

            mzIdentML110::SpectraDataType *tempVar19 = new mzIdentML110::SpectraDataType();
            tempVar19->id() = spectral_data_id;
            tempVar19->name() = data_filepath.substr(0, data_filepath.find_last_of("."));
            tempVar19->location() = data_filepath;
            tempVar19->FileFormat() = *new mzIdentML110::FileFormatType();
            tempVar19->FileFormat().get().cvParam() = *new mzIdentML110::CVParamType();
            tempVar19->FileFormat().get().cvParam().accession() = thermoRawFile ? "MS:1000563" : "MS:1000584";
            tempVar19->FileFormat().get().cvParam().name() = thermoRawFile ? "Thermo RAW format" : "mzML format";
            tempVar19->FileFormat().get().cvParam().cvRef() = "PSI-MS";
            tempVar19->SpectrumIDFormat() = *new mzIdentML110::SpectrumIDFormatType();
            tempVar19->SpectrumIDFormat().cvParam() =* new mzIdentML110::CVParamType();
            tempVar19->SpectrumIDFormat().cvParam().accession() = thermoRawFile ? "MS:1000768" : "MS:1001530";
            tempVar19->SpectrumIDFormat().cvParam().name() = thermoRawFile ? "Thermo nativeID format" : "mzML unique identifier";
            tempVar19->SpectrumIDFormat().cvParam().cvRef() = "PSI-MS";
            _mzid->DataCollection().Inputs().SpectraData()[spectra_data_id] = *tempVar19;
            spectra_data_id++;
        }
            
        int sir_id = 0;
        int pe_index = 0;
        int p_index = 0;
        std::unordered_map<PeptideWithSetModifications*, int> peptide_evidence_ids;
        //key is peptide sequence, value is <peptide id for that peptide, peptide evidences>, list of spectra id's
        std::unordered_map<std::string, std::tuple<int, std::unordered_set<std::string>>> peptide_ids; 
        //key is <filename, scan numer> value is <scan result id, scan item id #'s (could be more than one ID per scan)>
        //std::unordered_map<std::tuple<std::string, int>, std::tuple<int, int>> psm_per_scan; 
        std::unordered_map<MzIdentTuple, std::tuple<int, int>, MzIdentTuple_hash, MzIdentTuple_equal> psm_per_scan; 
        
#ifdef ORIG
        auto unambiguousPsms = psms.Where([&] (std::any psm)   {
                delete _mzid;
                delete _indexedSerializer;
                delete settings;
                return psm::FullSequence != nullptr;
            });
#endif
        std::vector<PeptideSpectralMatch*> unambiguousPsms;
        for ( auto psm : psms ) {
            if ( !psm->getFullSequence().empty() ) {
                unambiguousPsms.push_back(psm);
            }
        }
        
        for (auto psm : unambiguousPsms)
        {
#ifdef ORIG
            for (PeptideWithSetModifications *peptide : psm->BestMatchingPeptides->Select([&] (std::any p)  {
                        p::Peptide;
                    }).Distinct());
#endif
            std::vector<PeptideWithSetModifications *> bestMatchingPeptides;
            for ( auto p: psm->getBestMatchingPeptides() ) {
                bool found=false;
                for ( auto q: bestMatchingPeptides ){
                    if ( q == std::get<1>(p) ) {
                        found = true;
                        break;
                    }
                }
                if ( !found ) {
                    bestMatchingPeptides.push_back(std::get<1>(p));
                }
            }
            
            for ( PeptideWithSetModifications *peptide : bestMatchingPeptides )
            {
                
                //if first peptide on list hasn't been added, add peptide and peptide evidence
                std::tuple<int, std::unordered_set<std::string>> peptide_id;
                std::unordered_map<std::string, std::tuple<int, std::unordered_set<std::string>>>::const_iterator peptide_ids_iterator = peptide_ids.find(peptide->getFullSequence());
                if (peptide_ids_iterator == peptide_ids.end())
                {
                    peptide_id = std::tuple<int, std::unordered_set<std::string>>(p_index, std::unordered_set<std::string>());
                    p_index++;
                    mzIdentML110::PeptideType *tempVar20 = new mzIdentML110::PeptideType();
                    tempVar20->PeptideSequence() = peptide->getBaseSequence();
                    tempVar20->id() = "P_" + std::get<0>(peptide_id); 
                    tempVar20->Modification() = *new mzIdentML110::PeptideType::Modification_sequence(peptide->getNumMods());
                    //_mzid->SequenceCollection().Peptide[peptide_id->Item1] = tempVar20;
                    _mzid->SequenceCollection().get().Peptide()[std::get<0>(peptide_id)] = *tempVar20;
                    int mod_id = 0;
                    for (auto mod : peptide->getAllModsOneIsNterminus())
                    {
                        mzIdentML110::ModificationType *tempVar21 = new mzIdentML110::ModificationType();
                        tempVar21->location().set(std::get<0>(mod) - 1);
                        //tempVar21->locationSpecified() = true;
                        tempVar21->monoisotopicMassDelta().set(std::get<1>(mod)->getMonoisotopicMass().value());
                        //tempVar21->monoisotopicMassDeltaSpecified() = true;
                        auto s21 = std::to_string(peptide->getBaseSequence()[std::min(std::max(0, mod.first - 2), peptide->getLength() - 1)]);
                        mzIdentML110::listOfChars t21a (1, s21.c_str());
                        tempVar21->residues().set(t21a);
                        auto t21 = new mzIdentML110::ModificationType::cvParam_sequence();
                        t21->push_back(*GetUnimodCvParam(std::get<1>(mod)));
                        tempVar21->cvParam() = *t21;
                        _mzid->SequenceCollection().get().Peptide()[std::get<0>(peptide_id)].Modification()[mod_id] = *tempVar21;
                        mod_id++;
                    }
                    peptide_ids.emplace(peptide->getFullSequence(), peptide_id);
                }
                else  {
                    peptide_id = peptide_ids_iterator->second;
                }
                
                if (peptide_evidence_ids.find(peptide) == peptide_evidence_ids.end())
                {
                    mzIdentML110::PeptideEvidenceType *tempVar22 = new mzIdentML110::PeptideEvidenceType();
                    tempVar22->id() = "PE_" + std::to_string(pe_index);
                    tempVar22->peptide_ref() = "P_" + std::get<0>(peptide_id);
                    tempVar22->dBSequence_ref() = "DBS_" + peptide->getProtein()->getAccession();
                    tempVar22->isDecoy() = peptide->getProtein()->getIsDecoy();
                    //tempVar22->startSpecified() = true;
                    tempVar22->start().set(peptide->getOneBasedStartResidueInProtein());
                    //tempVar22->endSpecified() = true;
                    tempVar22->end().set(peptide->getOneBasedEndResidueInProtein());
                    auto t22  = peptide->getPreviousAminoAcid(); //.ToString();
                    tempVar22->pre() = &t22;
                    //tempVar22->post() = (peptide->getOneBasedEndResidueInProtein() < (int)peptide->getProtein()->getBaseSequence().size() ) ? peptide->getProtein()[peptide->getOneBasedEndResidueInProtein()] : "-";
                    if (peptide->getOneBasedEndResidueInProtein() < (int)peptide->getProtein()->getBaseSequence().size() ) {
                        tempVar22->post().set(peptide->getProtein()[peptide->getOneBasedEndResidueInProtein()].getName().c_str());
                    }
                    else {
                        tempVar22->post().set("-");
                    }
                    _mzid->SequenceCollection().get().PeptideEvidence()[pe_index] = *tempVar22;
                    peptide_evidence_ids.emplace(peptide, pe_index);
                    pe_index++;
                }
            }
            
            std::tuple<int, int> scan_result_scan_item;
            MzIdentTuple this_tuple = std::make_tuple(psm->getFullFilePath(), psm->getScanNumber());
            //std::unordered_map<std::tuple<std::string, int>, std::tuple<int, int>>::const_iterator psm_per_scan_iterator = psm_per_scan.find(this_tuple);
            std::unordered_map<MzIdentTuple, std::tuple<int, int>, MzIdentTuple_hash, MzIdentTuple_equal>::const_iterator psm_per_scan_iterator = psm_per_scan.find(this_tuple);
            if (psm_per_scan_iterator == psm_per_scan.end()) //check to see if scan has already been added
            {
                //scan_result_scan_item = psm_per_scan_iterator->second;
                scan_result_scan_item = std::make_tuple(sir_id, 0);
                mzIdentML110::SpectrumIdentificationResultType *tempVar23 = new mzIdentML110::SpectrumIdentificationResultType();
                tempVar23->id() = "SIR_" + std::get<0>(scan_result_scan_item);
                tempVar23->spectraData_ref() = "SD_" + std::to_string(spectral_ids[psm->getFullFilePath()]);
                tempVar23->spectrumID() = "scan=" + std::to_string(psm->getScanNumber());
                auto t23 = new mzIdentML110::SpectrumIdentificationResultType::SpectrumIdentificationItem_sequence(500);
                tempVar23->SpectrumIdentificationItem() = *t23;

                mzIdentML110::CVParamType *tempVar24 = new mzIdentML110::CVParamType();
                tempVar24->name() = "scan start time";
                tempVar24->cvRef() = "PSI-MS";
                tempVar24->accession() = "MS:1000016";
                tempVar24->value() = std::to_string(psm->getScanRetentionTime());
                auto t24 = new mzIdentML110::SpectrumIdentificationResultType::cvParam_sequence();
                t24->push_back(*tempVar24);
                tempVar23->cvParam() = *t24;
                _mzid->DataCollection().AnalysisData().SpectrumIdentificationList()[0].SpectrumIdentificationResult()[std::get<0>(scan_result_scan_item)] = *tempVar23;
                psm_per_scan.emplace(std::make_tuple(psm->getFullFilePath(), psm->getScanNumber()), scan_result_scan_item);
                sir_id++;
            }
            else
            {
                scan_result_scan_item = psm_per_scan_iterator->second;
                psm_per_scan[std::make_tuple(psm->getFullFilePath(), psm->getScanNumber())] =
                    std::make_tuple(std::get<0>(scan_result_scan_item), std::get<1>(scan_result_scan_item) + 1);
                scan_result_scan_item = psm_per_scan[std::make_tuple(psm->getFullFilePath(), psm->getScanNumber())];
                }

#ifdef ORIG
            for (PeptideWithSetModifications *p : psm->BestMatchingPeptides->Select([&] (std::any p) {
                        p->Peptide;
                    }).Distinct());
#endif
            for (PeptideWithSetModifications *p : bestMatchingPeptides )    {
                std::get<1>(peptide_ids[p->getFullSequence()]).insert("SII_" +
                                                                   std::to_string(std::get<0>(scan_result_scan_item))  + "_" +
                                                                   std::to_string(std::get<1>(scan_result_scan_item)) );
            }
            mzIdentML110::CVParamType *tempVar25 = new mzIdentML110::CVParamType();
            tempVar25->name() = "MetaMorpheus:score";
            tempVar25->cvRef() = "PSI-MS";
            tempVar25->accession() = "MS:1002827";
            tempVar25->value() = std::to_string(psm->getScore());

            mzIdentML110::CVParamType *tempVar26 = new mzIdentML110::CVParamType();
            tempVar26->accession() = "MS:1002354";
            tempVar26->name() = "PSM-level q-value";
            tempVar26->cvRef() = "PSI-MS";
            tempVar26->value() = std::to_string(psm->getFdrInfo()->getQValue());

            auto t26 = new mzIdentML110::SpectrumIdentificationItemType();
            t26->rank() = 1;
            t26->chargeState() = psm->getScanPrecursorCharge();
            t26->id() = "SII_" + std::to_string(std::get<0>(scan_result_scan_item)) + "_"
                + std::to_string(std::get<1>(scan_result_scan_item));
            t26->experimentalMassToCharge() = std::round(psm->getScanPrecursorMonoisotopicPeakMz()*std::pow(10,5))/std::pow(10,5);
            t26->passThreshold() = psm->getFdrInfo()->getQValue() <= 0.01;
            t26->peptide_ref() = "P_" + std::get<0>(peptide_ids[psm->getFullSequence()]);
#ifdef ORIG
            t26->PeptideEvidenceRef() = std::vector<mzIdentML110::PeptideEvidenceRefType*>(psm->getBestMatchingPeptides()->Select([&] (std::any p) {
                        p::Peptide;
                    }).Distinct()->Count()), cvParam = {tempVar25, tempVar26};
            
#endif
            auto t26a = new mzIdentML110::SpectrumIdentificationItemType::PeptideEvidenceRef_sequence(bestMatchingPeptides.size());
            t26->PeptideEvidenceRef() = *t26a;
            auto t26b = new mzIdentML110::SpectrumIdentificationItemType::cvParam_sequence();
            t26b->push_back(*tempVar25);
            t26b->push_back(*tempVar26);
            t26->cvParam() = *t26b;
                
            _mzid->DataCollection().AnalysisData().SpectrumIdentificationList()[0].SpectrumIdentificationResult()[std::get<0>(scan_result_scan_item)].SpectrumIdentificationItem()[std::get<1>(scan_result_scan_item)] = *t26;
            
            if (psm->getPeptideMonisotopicMass().has_value())
            {
                auto tt =std::round(Chemistry::ClassExtensions::ToMz( psm->getPeptideMonisotopicMass().value(),
                                                                      psm->getScanPrecursorCharge() )
                                    * std::pow(10, 5)) / std::pow(10, 5);
                _mzid->DataCollection().AnalysisData().SpectrumIdentificationList()[0].SpectrumIdentificationResult()[std::get<0>(scan_result_scan_item)].SpectrumIdentificationItem()[std::get<1>(scan_result_scan_item)].calculatedMassToCharge().set(tt);
                
                //_mzid->DataCollection().AnalysisData().SpectrumIdentificationList()[0].SpectrumIdentificationResult()[std::get<0>(scan_result_scan_item)].SpectrumIdentificationItem()[std::get<1>(scan_result_scan_item)]->calculatedMassToChargeSpecified() = true;
            }
            
            int pe = 0;
#ifdef ORIG
            for (PeptideWithSetModifications *p : psm->BestMatchingPeptides->Select([&] (std::any p) {
                        p->Peptide;
                    }).Distinct());
#endif
            for (PeptideWithSetModifications *p : bestMatchingPeptides )
            {
                mzIdentML110::PeptideEvidenceRefType *tempVar27 = new mzIdentML110::PeptideEvidenceRefType();
                tempVar27->peptideEvidence_ref() = "PE_" + std::to_string(peptide_evidence_ids[p]);
                _mzid->DataCollection().AnalysisData().SpectrumIdentificationList()[0].SpectrumIdentificationResult()[std::get<0>(scan_result_scan_item)].SpectrumIdentificationItem()[std::get<1>(scan_result_scan_item)].PeptideEvidenceRef()[pe] = *tempVar27;
                pe++;
            }
        }
        
        _mzid->AnalysisProtocolCollection() = *new mzIdentML110::AnalysisProtocolCollectionType();
        mzIdentML110::SpectrumIdentificationProtocolType *tempVar28 = new mzIdentML110::SpectrumIdentificationProtocolType();
        tempVar28->id() = "SIP";
        tempVar28->analysisSoftware_ref() = "AS_MetaMorpheus";
        tempVar28->SearchType() = *new mzIdentML110::ParamType();
        tempVar28->SearchType().cvParam() = *new mzIdentML110::CVParamType();
        tempVar28->SearchType().cvParam().get().accession() = "MS:1001083";
        tempVar28->SearchType().cvParam().get().name() = "ms-ms search";
        tempVar28->SearchType().cvParam().get().cvRef() = "PSI-MS";
        tempVar28->AdditionalSearchParams().get() = *new mzIdentML110::ParamListType();
        
        mzIdentML110::CVParamType *tempVar29 = new mzIdentML110::CVParamType();
        tempVar29->accession() = "MS:1001211";
        tempVar29->cvRef() = "PSI-MS";
        tempVar29->name() = "parent mass type mono";
        
        mzIdentML110::CVParamType *tempVar30 = new mzIdentML110::CVParamType();
        tempVar30->accession() = "MS:1001255";
        tempVar30->name() = "fragment mass type mono";
        tempVar30->cvRef() = "PSI-MS";
        auto t30 = new mzIdentML110::ParamListType::cvParam_sequence();
        t30->push_back(*tempVar29);
        t30->push_back(*tempVar30);
        tempVar28->AdditionalSearchParams().get().cvParam() = *t30;
        tempVar28->ModificationParams().get().SearchModification() = mzIdentML110::ModificationParamsType::SearchModification_sequence(fixedMods.size() + variableMods.size());
        tempVar28->Enzymes() = * new mzIdentML110::EnzymesType();
        tempVar28->Enzymes().get().Enzyme() = mzIdentML110::EnzymesType::Enzyme_sequence(proteases.size());
        
        mzIdentML110::CVParamType *tempVar31 = new mzIdentML110::CVParamType();
        tempVar31->accession() = "MS:1001412";
        tempVar31->name() = "search tolerance plus value";
        tempVar31->value() = std::to_string(productTolerance->getValue());
        tempVar31->cvRef() = "PSI-MS";
        tempVar31->unitAccession() = dynamic_cast<PpmTolerance*>(productTolerance) != nullptr? "UO:0000169": "UO:0000221";
        tempVar31->unitName() = dynamic_cast<PpmTolerance*>(productTolerance) != nullptr? "parts per million" : "dalton";
        tempVar31->unitCvRef() = "UO";
        
        mzIdentML110::CVParamType *tempVar32 = new mzIdentML110::CVParamType();
        tempVar32->accession() = "MS:1001413";
        tempVar32->name() = "search tolerance minus value";
        tempVar32->value() = std::to_string(productTolerance->getValue());
        tempVar32->cvRef() = "PSI-MS";
        tempVar32->unitAccession() = dynamic_cast<PpmTolerance*>(productTolerance) != nullptr? "UO:0000169": "UO:0000221";
        tempVar32->unitName() = dynamic_cast<PpmTolerance*>(productTolerance) != nullptr? "parts per million" : "dalton";
        tempVar32->unitCvRef() = "UO";
        auto t32 = new mzIdentML110::ToleranceType::cvParam_sequence();
        t32->push_back(*tempVar31);
        t32->push_back(*tempVar32);
        tempVar28->FragmentTolerance().get().cvParam() = *t32;
        
        mzIdentML110::CVParamType *tempVar33 = new mzIdentML110::CVParamType();
        tempVar33->accession() = "MS:1001412";
        tempVar33->name() = "search tolerance plus value";
        tempVar33->value() = std::to_string(parentTolerance->getValue());
        tempVar33->cvRef() = "PSI-MS";
        tempVar33->unitAccession() = dynamic_cast<PpmTolerance*>(parentTolerance) != nullptr? "UO:0000169": "UO:0000221";
        tempVar33->unitName() = dynamic_cast<PpmTolerance*>(parentTolerance) != nullptr? "parts per million" : "dalton";
        tempVar33->unitCvRef() = "UO";
        
        mzIdentML110::CVParamType *tempVar34 = new mzIdentML110::CVParamType();
        tempVar34->accession() = "MS:1001413";
        tempVar34->name() = "search tolerance minus value";
        tempVar34->value() = std::to_string(parentTolerance->getValue());
        tempVar34->cvRef() = "PSI-MS";
        tempVar34->unitAccession() = dynamic_cast<PpmTolerance*>(parentTolerance) != nullptr? "UO:0000169": "UO:0000221";
        tempVar34->unitName() = dynamic_cast<PpmTolerance*>(parentTolerance) != nullptr? "parts per million" : "dalton";
        tempVar34->unitCvRef() = "UO";
        auto t33 = new mzIdentML110::ToleranceType::cvParam_sequence();
        t33->push_back(*tempVar33);
        t33->push_back(*tempVar34);
        tempVar28->ParentTolerance().get().cvParam() = *t33;
        tempVar28->Threshold() = *new mzIdentML110::ParamListType();
        
        mzIdentML110::CVParamType *tempVar35 = new mzIdentML110::CVParamType();
        tempVar35->accession() = "MS:1001448";
        tempVar35->name() = "pep:FDR threshold";
        tempVar35->cvRef() = "PSI-MS";
        tempVar35->value() = "0.01";
        //tempVar28->Threshold()->Items = {tempVar35};
        auto t35a = new mzIdentML110::ParamListType::cvParam_sequence();
        t35a->push_back(*tempVar35);
        tempVar28->Threshold().cvParam() = *t35a;
        auto t35 = new mzIdentML110::AnalysisProtocolCollectionType::SpectrumIdentificationProtocol_sequence();
        t35->push_back(*tempVar28);
        _mzid->AnalysisProtocolCollection().SpectrumIdentificationProtocol() = *t35;
        
        int protease_index = 0;
        for (auto protease : proteases)
        {
            mzIdentML110::EnzymeType *tempVar36 = new mzIdentML110::EnzymeType();
            tempVar36->id() = "E_" + std::to_string(protease_index);
            tempVar36->name() = protease->getName();
            tempVar36->semiSpecific() = protease->getCleavageSpecificity() == CleavageSpecificity::Semi;
            //tempVar36->missedCleavagesSpecified() = true;
            tempVar36->missedCleavages().set( missedCleavages);
            tempVar36->EnzymeName() = *new mzIdentML110::ParamListType();
            
            mzIdentML110::CVParamType *tempVar37 = new mzIdentML110::CVParamType();
            tempVar37->accession() = protease->getPsiMsAccessionNumber();
            tempVar37->name() = protease->getPsiMsName();
            tempVar37->cvRef() = "PSI-MS";
            auto t37 = new mzIdentML110::ParamListType::cvParam_sequence();
            t37->push_back(*tempVar37);
            tempVar36->EnzymeName().get().cvParam() = *t37;
            _mzid->AnalysisProtocolCollection().SpectrumIdentificationProtocol()[0].Enzymes().get().Enzyme()[protease_index] = *tempVar36;
            protease_index++;
        }
        
        int mod_index = 0;
        for (auto mod : fixedMods)
        {
            mzIdentML110::SearchModificationType *tempVar38 = new mzIdentML110::SearchModificationType();
            tempVar38->fixedMod() = true;
            tempVar38->massDelta() = static_cast<float>(mod->getMonoisotopicMass().value());
            tempVar38->residues() = mod->getTarget()->ToString();
            //tempVar38->cvParam() = {GetUnimodCvParam(mod)};
            auto t38 = new mzIdentML110::SearchModificationType::cvParam_sequence();
            t38->push_back(*GetUnimodCvParam(mod));
            tempVar38->cvParam() = *t38;
            _mzid->AnalysisProtocolCollection().SpectrumIdentificationProtocol()[0].ModificationParams().get().SearchModification()[mod_index] = *tempVar38;
            mod_index++;
        }
        
        for (auto mod : variableMods)
        {
            mzIdentML110::SearchModificationType *tempVar39 = new mzIdentML110::SearchModificationType();
            tempVar39->fixedMod() = false;
            tempVar39->massDelta() = static_cast<float>(mod->getMonoisotopicMass().value());
            tempVar39->residues() = mod->getTarget()->ToString();
            auto t39 = new mzIdentML110::SearchModificationType::cvParam_sequence();
            t39->push_back(*GetUnimodCvParam(mod));
            tempVar39->cvParam() = *t39;
            _mzid->AnalysisProtocolCollection().SpectrumIdentificationProtocol()[0].ModificationParams().get().SearchModification()[mod_index] = *tempVar39;
            mod_index++;
        }
        
        _mzid->AnalysisProtocolCollection().ProteinDetectionProtocol() = * new mzIdentML110::ProteinDetectionProtocolType();
        _mzid->AnalysisProtocolCollection().ProteinDetectionProtocol().get().id() = "PDP";
        _mzid->AnalysisProtocolCollection().ProteinDetectionProtocol().get().analysisSoftware_ref() = "AS_MetaMorpheus";

        _mzid->AnalysisProtocolCollection().ProteinDetectionProtocol().get().Threshold() = *new mzIdentML110::ParamListType();
        
        mzIdentML110::CVParamType *tempVar40 = new mzIdentML110::CVParamType();
        tempVar40->accession() = "MS:1001447";
        tempVar40->name() = "prot:FDR threshold";
        tempVar40->cvRef() = "PSI-MS";
        tempVar40->value() = "0.01";
        auto t40 = new mzIdentML110::ParamListType::cvParam_sequence();
        t40->push_back(*tempVar40);
        _mzid->AnalysisProtocolCollection().ProteinDetectionProtocol().get().Threshold().cvParam() = *t40;
        
        if (groups.size() > 0)
        {
            _mzid->DataCollection().AnalysisData().ProteinDetectionList().get() = *new mzIdentML110::ProteinDetectionListType();
            _mzid->DataCollection().AnalysisData().ProteinDetectionList().get().id() = "PD";
            _mzid->DataCollection().AnalysisData().ProteinDetectionList().get().ProteinAmbiguityGroup() = mzIdentML110::ProteinDetectionListType::ProteinAmbiguityGroup_sequence(groups.size());
            
            int group_id = 0;
            int protein_id = 0;
            for (auto proteinGroup : groups)
            {
                mzIdentML110::ProteinAmbiguityGroupType *tempVar41 = new mzIdentML110::ProteinAmbiguityGroupType();
                tempVar41->id() = "PAG_" + std::to_string(group_id);
                tempVar41->ProteinDetectionHypothesis() = mzIdentML110::ProteinAmbiguityGroupType::ProteinDetectionHypothesis_sequence(proteinGroup->getProteins().size());
                _mzid->DataCollection().AnalysisData().ProteinDetectionList().get().ProteinAmbiguityGroup()[group_id] = *tempVar41;
                int pag_protein_index = 0;
                for (auto protein : proteinGroup->getProteins())
                {
                    mzIdentML110::ProteinDetectionHypothesisType *tempVar42 = new mzIdentML110::ProteinDetectionHypothesisType();
                    tempVar42->id() = "PDH_" + std::to_string(protein_id);
                    tempVar42->dBSequence_ref() = "DBS_" + protein->getAccession();
                    tempVar42->passThreshold() = proteinGroup->getQValue() <= 0.01;
                    tempVar42->PeptideHypothesis() = mzIdentML110::ProteinDetectionHypothesisType::PeptideHypothesis_sequence(proteinGroup->getAllPeptides().size());
                    mzIdentML110::CVParamType *tempVar43 = new mzIdentML110::CVParamType();
                    tempVar43->accession() = "MS:1002828";
                    tempVar43->name() = "MetaMorpheus:protein score";
                    tempVar43->cvRef() = "PSI-MS";
                    tempVar43->value() = std::to_string(proteinGroup->getProteinGroupScore());
                    
                    mzIdentML110::CVParamType *tempVar44 = new mzIdentML110::CVParamType();
                    tempVar44->accession() = "MS:1002373";
                    tempVar44->name() = "protein group-level q-value";
                    tempVar44->cvRef() = "PSI-MS";
                    tempVar44->value() = std::to_string(proteinGroup->getQValue());
                    
                    mzIdentML110::CVParamType *tempVar45 = new mzIdentML110::CVParamType();
                    tempVar45->accession() = "MS:1001093";
                    tempVar45->name() = "sequence coverage";
                    tempVar45->cvRef() = "PSI-MS";
                    tempVar45->value() = std::to_string(proteinGroup->getSequenceCoveragePercent().front());
                    
                    mzIdentML110::CVParamType *tempVar46 = new mzIdentML110::CVParamType();
                    tempVar46->accession() = "MS:1001097";
                    tempVar46->name() = "distinct peptide sequences";
                    tempVar46->cvRef() = "PSI-MS";
                    tempVar46->value() = std::to_string(proteinGroup->getUniquePeptides().size());
                    auto t46 = new mzIdentML110::ProteinDetectionHypothesisType::cvParam_sequence();
                    t46->push_back(*tempVar43);
                    t46->push_back(*tempVar44);
                    t46->push_back(*tempVar45);
                    t46->push_back(*tempVar46);
                    tempVar42->cvParam() = *t46;
                    _mzid->DataCollection().AnalysisData().ProteinDetectionList().get().ProteinAmbiguityGroup()[group_id].ProteinDetectionHypothesis()[pag_protein_index] = *tempVar42;
                    
                    int peptide_id = 0;
                    for (auto peptide : proteinGroup->getAllPeptides())
                    {
                        if (peptide_evidence_ids.find(peptide) != peptide_evidence_ids.end())
                        {
                            if (peptide->getProtein() == protein)
                            {
                                mzIdentML110::PeptideHypothesisType *tempVar47 = new mzIdentML110::PeptideHypothesisType();
                                tempVar47->peptideEvidence_ref() = "PE_" + std::to_string(peptide_evidence_ids[peptide]);
                                tempVar47->SpectrumIdentificationItemRef() = mzIdentML110::PeptideHypothesisType::SpectrumIdentificationItemRef_sequence(std::get<1>(peptide_ids[peptide->getFullSequence()]).size());
                                _mzid->DataCollection().AnalysisData().ProteinDetectionList().get().ProteinAmbiguityGroup()[group_id].ProteinDetectionHypothesis()[pag_protein_index].PeptideHypothesis()[peptide_id] = *tempVar47;
                                
                                int i = 0;
                                for (std::string sii : std::get<1>(peptide_ids[peptide->getFullSequence()]))
                                {
                                    mzIdentML110::SpectrumIdentificationItemRefType *tempVar48 = new mzIdentML110::SpectrumIdentificationItemRefType();
                                    tempVar48->spectrumIdentificationItem_ref() = sii;
                                    _mzid->DataCollection().AnalysisData().ProteinDetectionList().get().ProteinAmbiguityGroup()[group_id].ProteinDetectionHypothesis()[pag_protein_index].PeptideHypothesis()[peptide_id].SpectrumIdentificationItemRef()[i] = *tempVar48;
                                    i++;
                                }
                                peptide_id++;
                            }
                        }
                    }
                    pag_protein_index++;
                    protein_id++;
                }
                group_id++;
            }
        }
        //XmlWriter *writer = XmlWriter::Create(outputPath, settings);
        //_indexedSerializer->Serialize(writer, _mzid);
        //writer->Close();
        //delete _indexedSerializer;

        // Serialize the object model to XML.
        //
        xml_schema::namespace_infomap map;
        map[""].name = "";
        map[""].schema = "/home/gabriel/XLMS/mzlib-master/MzIdentML/mzIdentML110.xsd";

        // Serialize to a file.
        try{
            std::ofstream ofs (outputPath);
            mzIdentML110::MzIdentML (ofs, *_mzid, map);
            ofs.close();
        }

        catch (const xml_schema::exception& e)
        {
            std::cerr << e << std::endl;
        }

        
        delete _mzid;
    }
    
    mzIdentML110::CVParamType *MzIdentMLWriter::GetUnimodCvParam(Modification *mod)
    {
        if ( !mod->getDatabaseReference().empty() &&
             mod->getDatabaseReference().find("Unimod") != mod->getDatabaseReference().end())
        {
            mzIdentML110::CVParamType *tempVar = new mzIdentML110::CVParamType();
            tempVar->accession() = "UNIMOD:" + mod->getDatabaseReference()["Unimod"].front();
            tempVar->name() = mod->getIdWithMotif();
            tempVar->cvRef() = "PSI-MS";
            return tempVar;
        }
        else
        {
            mzIdentML110::CVParamType *tempVar2 = new mzIdentML110::CVParamType();
            tempVar2->accession() = "MS:1001460";
            tempVar2->name() = "unknown modification";
            tempVar2->cvRef() = "UNIMOD";
            tempVar2->value() = mod->getIdWithMotif();
            return tempVar2;
        }
    }
}
