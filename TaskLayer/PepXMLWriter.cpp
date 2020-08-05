#include "time.h"
#include "PepXMLWriter.h"
#include "../EngineLayer/CommonParameters.h"
#include "../EngineLayer/PeptideSpectralMatch.h"
#include "DbForTask.h"
#include "../EngineLayer/GlobalVariables.h"

#include "pepXML/pepXML_v120.h"
#include <sstream>

using namespace EngineLayer;
using namespace Proteomics;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace TaskLayer
{

	void PepXMLWriter::WritePepXml(std::vector<PeptideSpectralMatch*> &psms,
                                       std::vector<DbForTask*> &database,
                                       std::vector<Modification*> &variableModifications,
                                       std::vector<Modification*> &fixedModifications,
                                       CommonParameters *CommonParameters,
                                       const std::string &outputPath,
                                       double qValueFilter)
	{

#ifdef ORIG
            // TODO: needs a unit test
            psms = psms.Where([&] (std::any p) {
                    return p::FdrInfo::QValue <= qValueFilter && p::FdrInfo::QValueNotch < qValueFilter;
                }).ToList();
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
            if ( psms.empty()) {
                return;
            }

            //XmlSerializer *_indexedSerializer = new XmlSerializer(pepXML::Generated::msms_pipeline_analysis::typeid);
            //auto _pepxml = new pepXML::Generated::msms_pipeline_analysis();
            auto _pepxml = new pepXML::msms_pipeline_analysis();

            time_t timer;
            time(&timer);
            struct tm *tmi = localtime(&timer);
            
            short zone_hours=0;
            short zone_minutes=0;

            ::xml_schema::date_time *dt = new ::xml_schema::date_time(tmi->tm_year, tmi->tm_mon, tmi->tm_mday, tmi->tm_hour,
                                 tmi->tm_min, (double)tmi->tm_sec, zone_hours, zone_minutes);
            _pepxml->date(*dt);
            _pepxml->summary_xml(psms[0]->getFullFilePath() + ".pep.XML");

#ifdef ORIG
            std::string proteaseNC = std::string::Join("", CommonParameters->getDigestionParams()->Protease.DigestionMotifs->Select([&] (std::any m) {
                        m::InducingCleavage;
                    }));
#endif
            std::string proteaseNC;
            for ( auto m : CommonParameters->getDigestionParams()->getProtease()->getDigestionMotifs() ) {
                proteaseNC += m->InducingCleavage;
            }

#ifdef ORIG
            std::string proteaseC = std::string::Join("", CommonParameters->getDigestionParams()->Protease.DigestionMotifs->Select([&] (std::any m){
			m::InducingCleavage;
                    }));
#endif
            std::string proteaseC = proteaseNC; 

            //std::string fileNameNoExtension = Path::GetFileNameWithoutExtension(psms[0]->getFullFilePath());
            std::string tmps = psms[0]->getFullFilePath();
            std::string fileNameNoExtension = tmps.substr(0, tmps.find_last_of("."));
            //std::string filePathNoExtension = Path::ChangeExtension(psms[0]->getFullFilePath(), "");
            std::string filePathNoExtension = tmps.substr(0, tmps.find_last_of("/"));

            //auto para = std::vector<pepXML::nameValueType*>();
            auto para = new pepXML::search_summary::parameter_sequence();
            {
                pepXML::nameValueType *tempVar = new pepXML::nameValueType();
                tempVar->name("threads");
                tempVar->value(std::to_string(CommonParameters->getMaxThreadsToUsePerFile()));
                para->push_back(*tempVar);
                delete tempVar;
                
                pepXML::nameValueType *tempVar2 = new pepXML::nameValueType();
                tempVar2->name("database");
                tempVar2->value(database.front()->getFilePath());
                para->push_back(*tempVar2);
                delete tempVar2;
                
                pepXML::nameValueType *tempVar3 = new pepXML::nameValueType();
                tempVar3->name("MS_data_file");
                tempVar3->value(psms[0]->getFullFilePath());
                para->push_back(*tempVar3);
                delete tempVar3;
                
                pepXML::nameValueType *tempVar4 = new pepXML::nameValueType();
                tempVar4->name("MaxMissed Cleavages");                
                tempVar4->value(std::to_string(CommonParameters->getDigestionParams()->getMaxMissedCleavages()));
                para->push_back(*tempVar4);
                delete tempVar4;
                
                pepXML::nameValueType *tempVar5 = new pepXML::nameValueType();
                tempVar5->name("Protease");
                tempVar5->value(CommonParameters->getDigestionParams()->getProtease()->getName());
                para->push_back(*tempVar5);
                delete tempVar5;
                
                pepXML::nameValueType *tempVar6 = new pepXML::nameValueType();
                tempVar6->name("Initiator Methionine");
                tempVar6->value(std::to_string((int)(CommonParameters->getDigestionParams()->getInitiatorMethionineBehavior())));
                para->push_back(*tempVar6);
                delete tempVar6;
                
                pepXML::nameValueType *tempVar7 = new pepXML::nameValueType();
                tempVar7->name("Max Modification Isoforms");
                tempVar7->value(std::to_string(CommonParameters->getDigestionParams()->getMaxModificationIsoforms()));
                para->push_back(*tempVar7);
                delete tempVar7;
                
                pepXML::nameValueType *tempVar8 = new pepXML::nameValueType();
                tempVar8->name("Min Peptide Len");
                tempVar8->value(std::to_string(CommonParameters->getDigestionParams()->getMinPeptideLength()));
                para->push_back(*tempVar8);
                delete tempVar8;
                
                pepXML::nameValueType *tempVar9 = new pepXML::nameValueType();
                tempVar9->name("Max Peptide Len");
                tempVar9->value(std::to_string(CommonParameters->getDigestionParams()->getMaxPeptideLength()));
                para->push_back(*tempVar9);
                delete tempVar9;
                
                pepXML::nameValueType *tempVar10 = new pepXML::nameValueType();
                tempVar10->name("Product Mass Tolerance");
                tempVar10->value(std::to_string(CommonParameters->getProductMassTolerance()->getValue()));
                para->push_back(*tempVar10);
                delete tempVar10;

                // TODO: check this
                pepXML::nameValueType *tempVar11 = new pepXML::nameValueType();
                tempVar11->name("Ions to search");
                // tempVar11->value() = std::string::Join(", ", DissociationTypeCollection::ProductsFromDissociationType[CommonParameters->getDissociationType()]);
                std:: string tmps2;
                for ( auto i: DissociationTypeCollection::ProductsFromDissociationType[CommonParameters->getDissociationType()] ) {
                    tmps2 += std::to_string((int)(i)); 
                }
                tempVar11->value(tmps2);
                para->push_back(*tempVar11);
                delete tempVar11;
                
                pepXML::nameValueType *tempVar12 = new pepXML::nameValueType();
                tempVar12->name("Q-value Filter");
                tempVar12->value(std::to_string(CommonParameters->getQValueOutputFilter()));
                para->push_back(*tempVar12);
                delete tempVar12;               

                for (auto item : fixedModifications)
                {
                    pepXML::nameValueType *tempVar13 = new pepXML::nameValueType();
                    tempVar13->name("Fixed Modifications: " + item->getIdWithMotif());                    
                    tempVar13->value(std::to_string(item->getMonoisotopicMass().value()));
                    para->push_back(*tempVar13);
                    delete tempVar13;
                }

                for (auto item : variableModifications)
                {
                    pepXML::nameValueType *tempVar14 = new pepXML::nameValueType();
                    tempVar14->name("Variable Modifications: " + item->getIdWithMotif());                    
                    tempVar14->value(std::to_string(item->getMonoisotopicMass().value()));
                    para->push_back(*tempVar14);
                    delete tempVar14;
                }
                
                pepXML::nameValueType *tempVar15 = new pepXML::nameValueType();
                tempVar15->name("Localize All Modifications");
                tempVar15->value("true");
                para->push_back(*tempVar15);
                delete tempVar15;
            }
            
            pepXML::msms_run_summary *tempVar16 = new pepXML::msms_run_summary();
            tempVar16->base_name(filePathNoExtension);
            tempVar16->raw_data_type("raw");
            tempVar16->raw_data(".mzM");

            auto topt = new pepXML::sample_enzyme();
            tempVar16->sample_enzyme(*topt);
            tempVar16->sample_enzyme()->name(CommonParameters->getDigestionParams()->getProtease()->getName());
            pepXML::specificity *tempVar17 = new pepXML::specificity();
            tempVar17->cut(proteaseC);
            tempVar17->no_cut(proteaseNC);

            auto t17 = new pepXML::sample_enzyme::specificity_sequence();
            t17->push_back(*tempVar17);
            delete tempVar17;
            
            tempVar16->sample_enzyme()->specificity(*t17);
            delete t17;            

            pepXML::search_summary *tempVar18 = new pepXML::search_summary();
            tempVar18->base_name(filePathNoExtension);
            tempVar18->search_engine_version(GlobalVariables::getMetaMorpheusVersion());
            tempVar18->precursor_mass_type(pepXML::massType::monoisotopic);
            tempVar18->fragment_mass_type(pepXML::massType::monoisotopic);
            tempVar18->search_id(1);

            tempVar18->search_database(*(new pepXML::search_database())); 
            tempVar18->search_database()->local_path(database.front()->getFilePath());

            tempVar18->search_database()->type(pepXML::type::value::AA);
            tempVar18->enzymatic_search_constraint(*(new pepXML::enzymatic_search_constraint()));
            tempVar18->enzymatic_search_constraint()->enzyme(CommonParameters->getDigestionParams()->getProtease()->getName());
            
            tempVar18->enzymatic_search_constraint()->max_num_internal_cleavages(CommonParameters->getDigestionParams()->getMaxMissedCleavages());
            tempVar18->parameter(*para);
            auto t18 = new pepXML::msms_run_summary::search_summary_sequence();
            t18->push_back(*tempVar18);
            delete tempVar18;
            
            tempVar16->search_summary(*t18);
            delete t18;

            auto tt18 = new pepXML::msms_pipeline_analysis::msms_run_summary_sequence();
            tt18->push_back(*tempVar16);
            delete tempVar16;
            
            _pepxml->msms_run_summary(*tt18);
            delete tt18;
            
            auto ttt18 = new pepXML::msms_run_summary::spectrum_query_sequence(psms.size());
            _pepxml->msms_run_summary()[0].spectrum_query(*ttt18);
            delete ttt18;
            
            auto searchHits = std::vector<pepXML::search_hit*>();
            
            for (auto psm : psms)
            {
                PeptideWithSetModifications *peptide = std::get<1>(psm->getBestMatchingPeptides().front());
                
                //auto mods = std::vector<pepXML::mod_aminoacid_mass*>();
                auto mods = new pepXML::modInfoDataType::mod_aminoacid_mass_sequence();
                for (auto mod : peptide->getAllModsOneIsNterminus())
                {
                    auto pepXmlMod = new pepXML::mod_aminoacid_mass();
                    pepXmlMod->mass(static_cast<double>(std::get<1>(mod)->getMonoisotopicMass().value()));
                    
                    pepXmlMod->position(std::get<0>(mod) - 1);
                    mods->push_back(*pepXmlMod);
                    delete pepXmlMod;
                }
                
#ifdef ORIG
                auto proteinAccessions = psm->BestMatchingPeptides->Select([&] (std::any p)  {
                        p::Peptide::Protein::Accession;
                    }).Distinct();
#endif
                std::vector<std::string> proteinAccessions;
                std::string proteinAccessionsString;
                for ( auto p : psm->getBestMatchingPeptides() ) {
                    if ( proteinAccessions.size() == 0 ) {
                        proteinAccessions.push_back(std::get<1>(p)->getProtein()->getAccession() );
                        proteinAccessionsString += std::get<1>(p)->getProtein()->getAccession();
                        continue;
                    }
                    bool found = false;
                    for ( auto w: proteinAccessions ) {
                        if (std::get<1>(p)->getProtein()->getAccession() == w ) {
                            found = true;
                            break;
                        }
                    }
                    if ( !found ) {
                        proteinAccessions.push_back(std::get<1>(p)->getProtein()->getAccession() );
                        proteinAccessionsString += "|" + std::get<1>(p)->getProtein()->getAccession();                    }
                }
                
                auto searchHit = new pepXML::search_hit();
                searchHit->hit_rank(1);
                searchHit->peptide(((psm->getBaseSequence() != "") ? psm->getBaseSequence() : "Ambiguous"));

                std::stringstream ss;
                ss << peptide->getPreviousAminoAcid();
                searchHit->peptide_prev_aa(ss.str() );
                ss.str("");
                ss << peptide->getNextAminoAcid();
                searchHit->peptide_next_aa(ss.str() );
                searchHit->protein(((peptide->getProtein()->getAccession() != "") ? peptide->getProtein()->getAccession() : proteinAccessionsString));
                searchHit->num_tot_proteins(static_cast<unsigned int>(proteinAccessions.size()));
                searchHit->calc_neutral_pep_mass(static_cast<float>((psm->getPeptideMonisotopicMass().has_value()) ?
                                                                    psm->getPeptideMonisotopicMass().value() : NAN));
#ifdef ORIG
                searchHit->massdiff() = ((psm->getPeptideMonisotopicMass().has_value()) ?
                                         std::to_string(psm->getScanPrecursorMass() - psm->getPeptideMonisotopicMass().value())
                                         : "Ambiguous");
#endif
                // massdiff() seems to be a double in this version, not a string.
                if (psm->getPeptideMonisotopicMass().has_value() ) {
                    searchHit->massdiff(psm->getScanPrecursorMass() - psm->getPeptideMonisotopicMass().value());
                }
                else {
                    searchHit->massdiff(-1);
                }
                
                pepXML::modInfoDataType *tempVar19 = new pepXML::modInfoDataType();
                tempVar19->mod_aminoacid_mass(*mods);
                if ( mods->empty () ) {
                    searchHit->modification_info(*tempVar19);
                }

                pepXML::nameValueType *tempVar20 = new pepXML::nameValueType();
                tempVar20->name("Score");
                tempVar20->value(std::to_string(psm->getScore()));

                pepXML::nameValueType *tempVar21 = new pepXML::nameValueType();
                tempVar21->name("Qvalue");
                tempVar21->value(std::to_string(psm->getFdrInfo()->getQValue()));

                //delete tempVar20;
                //delete tempVar21;
                searchHit->search_score() = {tempVar20, tempVar21};
                searchHits.push_back(searchHit);
                
                //C# TO C++ CONVERTER TODO TASK: A 'delete searchHit' statement was not added since searchHit
                // was passed to a method or constructor. Handle memory management manually.
            }
            
            for (int i = 0; i < (int)psms.size(); i++)
            {
                pepXML::spectrum_query *tempVar22 = new pepXML::spectrum_query();
                
                tempVar22->spectrum(fileNameNoExtension + "." + std::to_string(psms[i]->getScanNumber()));
                tempVar22->start_scan(static_cast<unsigned int>(psms[i]->getScanNumber()));
                tempVar22->end_scan(static_cast<unsigned int>(psms[i]->getScanNumber()));
                tempVar22->precursor_neutral_mass(static_cast<float>(psms[i]->getScanPrecursorMass()));
                
                tempVar22->assumed_charge(psms[i]->getScanPrecursorCharge());
                tempVar22->index((unsigned int)(i + 1));
                tempVar22->retention_time_sec(static_cast<float>(psms[i]->getScanRetentionTime() * 60));
                pepXML::search_result *tempVar23 = new pepXML::search_result();
                auto tt23  = new pepXML::search_result::search_hit_sequence();
                tt23->push_back(*searchHits[i]);
                delete searchHits[i];
                
                tempVar23->search_hit(*tt23);
                delete tt23;
                
                auto t23 = new pepXML::spectrum_query::search_result_sequence();
                t23->push_back(*tempVar23);
                delete tempVar23;
                
                tempVar22->search_result(*t23);
                delete t23;
                
                _pepxml->msms_run_summary()[0].spectrum_query()[i] = *tempVar22;
                delete tempVar22;
            }
            // Serialize the object model to XML.
            //
            xml_schema::namespace_infomap map;
            map[""].name = "";
            map[""].schema = "/home/gabriel/XLMS/mzlib-master/pepXML/pepXML_v120.xsd";

            // Serialize to a file.
            try{
                std::ofstream ofs (outputPath);
                pepXML::msms_pipeline_analysis_ (ofs, *_pepxml, map);
                ofs.close();
            }

            catch (const xml_schema::exception& e)
            {
                std::cerr << e << std::endl;
            }

            delete _pepxml;
	}
}
