#include "Bin.h"
#include "../PeptideSpectralMatch.h"
#include "../GlobalVariables.h"

#include "stringhelper.h"

using namespace Proteomics;
using namespace Proteomics::AminoAcidPolymer;
namespace EngineLayer
{
    namespace HistogramAnalysis
    {
        
        Bin::Bin(double massShift)
        {
            this->privateMassShift = massShift;
            //UniquePSMs = std::unordered_map<std::string, std::tuple<std::string, std::string, PeptideSpectralMatch*>>();
        }
        
        int Bin::getPepNlocCount() const
        {
            return privatePepNlocCount;
        }
        
        void Bin::setPepNlocCount(int value)
        {
            privatePepNlocCount = value;
        }
        
        int Bin::getPepClocCount() const
        {
            return privatePepClocCount;
        }
        
        void Bin::setPepClocCount(int value)
        {
            privatePepClocCount = value;
        }
        
        int Bin::getProtNlocCount() const
        {
            return privateProtNlocCount;
        }
        
        void Bin::setProtNlocCount(int value)
        {
            privateProtNlocCount = value;
        }
        
        int Bin::getProtClocCount() const
        {
            return privateProtClocCount;
        }
        
        void Bin::setProtClocCount(int value)
        {
            privateProtClocCount = value;
        }
        
        std::string Bin::getCombos() const
        {
            return privateCombos;
        }
        
        void Bin::setCombos(const std::string &value)
        {
            privateCombos = value;
        }
        
        std::string Bin::getUnimodDiffs() const
        {
            return privateUnimodDiffs;
        }
        
        void Bin::setUnimodDiffs(const std::string &value)
        {
            privateUnimodDiffs = value;
        }
        
        std::string Bin::getUniprotID() const
        {
            return privateUniprotID;
        }
        
        void Bin::setUniprotID(const std::string &value)
        {
            privateUniprotID = value;
        }
        
        std::string Bin::getUnimodFormulas() const
        {
            return privateUnimodFormulas;
        }

        void Bin::setUnimodFormulas(const std::string &value)
        {
            privateUnimodFormulas = value;
        }
        
        std::string Bin::getUnimodId() const
        {
            return privateUnimodId;
        }
        
        void Bin::setUnimodId(const std::string &value)
        {
            privateUnimodId = value;
        }
        
        double Bin::getMassShift() const
        {
            return privateMassShift;
        }
        
        int Bin::getCount() const
        {
            return UniquePSMs.size();
        }
        
        int Bin::getCountDecoy() const
        {
#ifdef ORIG
            return UniquePSMs.Values->Count([&] (std::any b)  {
                    b::Item3->IsDecoy;
                });
#endif
            int count = 0;
            for ( auto b:  UniquePSMs ) {
                auto q =  std::get<1>(b);
                if ( std::get<2>(q)->getIsDecoy() ){
                    count++;
                }
            }
            return count;
        }

        int Bin::getCountTarget() const
        {
            return getCount() - getCountDecoy();
        }
        
        int Bin::getLocalizeableTarget() const
        {
#ifdef ORIG
            return UniquePSMs.Values->Where([&] (std::any b)  {
                    return b::Item3->LocalizedScores != nullptr;
                })->Count([&] (std::any b)   {
                        return !b::Item3->IsDecoy && b::Item3->LocalizedScores.Max() >= b::Item3->Score + 1;
                    });
#endif
            int count=0;
            for ( auto b: UniquePSMs ) {
                auto q = std::get<1>(b);
                if ( std::get<2>(q)->getLocalizedScores().size() != 0  ) {
                    double max = std::get<2>(q)->getLocalizedScores()[0];
                    for (auto r : std::get<2>(q)->getLocalizedScores() ) {
                        if ( r > max ) max = r;
                    }
                    if ( std::get<2>(q)->getIsDecoy() &&
                         max >= (std::get<2>(q)->getScore()+1) ) {
                        count++;
                    }
                }
            }
            return count;
        }
        
        std::string Bin::getMine() const
        {
            return privateMine;
        }
        
        void Bin::setMine(const std::string &value)
        {
            privateMine = value;
        }
        
        std::unordered_map<char, int> Bin::getAAsInCommon() const
        {
            return privateAAsInCommon;
        }
        
        void Bin::setAAsInCommon(const std::unordered_map<char, int> &value)
        {
            privateAAsInCommon = value;
        }
        
        int Bin::getOverlapping() const
        {
            return privateOverlapping;
        }
        
        void Bin::setOverlapping(int value)
        {
            privateOverlapping = value;
        }
        
        double Bin::getFracWithSingle() const
        {
            return privateFracWithSingle;
        }
        
        void Bin::setFracWithSingle(double value)
        {
            privateFracWithSingle = value;
        }
        
        double Bin::getMedianLength() const
        {
            return privateMedianLength;
        }
        
        void Bin::setMedianLength(double value)
        {
            privateMedianLength = value;
        }
        
        void Bin::IdentifyResidues()
        {
            std::unordered_map<char, int> ResidueCount;

#ifdef ORIG
            //for (auto hehe : UniquePSMs.Values->Where([&] (std::any b)  {
            //            return b::Item3->LocalizedScores != nullptr;
            //        }))
#endif
            for ( auto h :  UniquePSMs ) 
            {
                auto hehe = std::get<1>(h);
                if ( std::get<2>(hehe)->getLocalizedScores().size() == 0 ) {
                    continue;
                }
                
#ifdef ORIG
                double bestScore = hehe::Item3->LocalizedScores.Max();
#endif
                double bestScore = std::get<2>(hehe)->getLocalizedScores()[0];
                for ( auto q: std::get<2>(hehe)->getLocalizedScores() ) {
                    if ( q > bestScore ) bestScore = q;                    
                }
                
                if ( bestScore >= std::get<2>(hehe)->getScore() + 1 &&
                     !std::get<2>(hehe)->getIsDecoy() )
                {
                    for (int i = 0; i < (int)std::get<0>(hehe).size(); i++)
                    {
                        if (bestScore - std::get<2>(hehe)->getLocalizedScores()[i] < 0.5)
                        {
                            if (ResidueCount.find(std::get<0>(hehe)[i]) != ResidueCount.end())
                            {
                                ResidueCount[std::get<0>(hehe)[i]]++;
                            }
                            else
                            {
                                ResidueCount.emplace(std::get<0>(hehe)[i], 1);
                            }
                        }
                    }
                    if ( bestScore - std::get<2>(hehe)->getLocalizedScores()[0] < 0.5)
                    {
                        setPepNlocCount(getPepNlocCount() + 1);
                        if (std::get<2>(hehe)->getOneBasedStartResidueInProtein().has_value() &&
                            std::get<2>(hehe)->getOneBasedStartResidueInProtein().value() <= 2)
                        {
                            setProtNlocCount(getProtNlocCount() + 1);
                        }
                    }
                    if ( bestScore - std::get<2>(hehe)->getLocalizedScores().back() < 0.5)
                    {
                        setPepClocCount(getPepClocCount() + 1);
                        if (std::get<2>(hehe)->getOneBasedEndResidueInProtein().has_value() &&
                            std::get<2>(hehe)->getProteinLength().has_value() &&
                            std::get<2>(hehe)->getOneBasedEndResidueInProtein().value() == std::get<2>(hehe)->getProteinLength().value() )
                        {
                            setProtClocCount(getProtClocCount() + 1);
                        }
                    }
                }
            }
        }
        
        void Bin::IdentifyCombos(double v, BinTuple_set &ok)
        {
            std::vector<std::string> okk;
            for (auto hm = ok.begin(); hm != ok.end(); hm++ )
            {
                if (std::abs(std::get<0>(*hm) + std::get<1>(*hm) - getMassShift()) <= v &&
                    getCountTarget() < std::get<2>(*hm))
                {
                    std::string s = "Combo " + std::to_string(std::min(std::get<0>(*hm), std::get<1>(*hm))) + " and " +
                        std::to_string(std::max(std::get<0>(*hm), std::get<1>(*hm)));
                    okk.push_back(s);
                }
            }
            std::string del = "|";
            setCombos(StringHelper::join(okk, del) );
        }
        
        double Bin::ComputeZ(double v)
        {
            return std::sqrt(getCount()) * (static_cast<double>(getCountDecoy()) / getCount() - v) / (v * (1 - v));
        }
        
        void Bin::IdentifyUniprotBins(double v)
        {
            std::unordered_set<std::string> modIdOptions;
            for (auto mod : GlobalVariables::getUniprotDeseralized())
            {
                if (mod->getMonoisotopicMass().has_value() &&
                    std::abs(mod->getMonoisotopicMass().value() - getMassShift()) <= v)
                {
                    modIdOptions.insert(mod->getIdWithMotif());
                }
            }
            std::string del = "|";
            setUniprotID(StringHelper::join(modIdOptions, del));
        }
        
        void Bin::IdentifyAA(double v)
        {
            std::unordered_set<std::string> ok;
            for (char c = 'A'; c <= 'Z'; c++)
            {
                Residue *residue;
                if (Residue::TryGetResidue(c, &residue))
                {
                    if (std::abs(residue->getMonoisotopicMass() - getMassShift()) <= v)
                    {
                        ok.insert("Add " + residue->getName());
                    }
                    if (std::abs(residue->getMonoisotopicMass() + getMassShift()) <= v)
                    {
                        ok.insert("Remove " + residue->getName());
                    }
                    for (char cc = 'A'; cc <= 'Z'; cc++)
                    {
                        Residue *residueCC;
                        if (Residue::TryGetResidue(cc, &residueCC))
                        {
                            if (std::abs(residueCC->getMonoisotopicMass() + residue->getMonoisotopicMass() - getMassShift()) <= v)
                            {
                                ok.insert("Add (" + residue->getName() + "+" + residueCC->getName() + ")");
                            }
                            if (std::abs(residueCC->getMonoisotopicMass() + residue->getMonoisotopicMass() + getMassShift()) <= v)
                            {
                                ok.insert("Remove (" + residue->getName() + "+" + residueCC->getName() + ")");
                            }
                        }
                    }
                }
            }
            std::string del = "|";
            AA = StringHelper::join(ok, del);
        }
        
        void Bin::IdentifyUnimodBins(double v)
        {
            std::unordered_set<std::string> ok;
            std::unordered_set<std::string> okformula;
            std::unordered_set<double> okDiff;
            for (auto hm : GlobalVariables::getUnimodDeserialized())
            {
                auto theMod = dynamic_cast<Modification*>(hm);
                if (std::abs(theMod->getMonoisotopicMass().value() - getMassShift()) <= v)
                {
                    ok.insert(hm->getIdWithMotif());
                    okformula.insert(theMod->getChemicalFormula()->getFormula());
                    okDiff.insert(theMod->getMonoisotopicMass().value() - getMassShift());
                }
            }
            std::string del="|";
            setUnimodId(StringHelper::join(ok, del));
            setUnimodFormulas(StringHelper::join(okformula, del));
            std::vector<std::string>tmpvec;
            for ( auto p =okDiff.begin(); p!=okDiff.end(); p++ ) {
                tmpvec.push_back(std::to_string(*p));
            }
            setUnimodDiffs(StringHelper::join(tmpvec, del));
        }
        
        void Bin::Add(PeptideSpectralMatch *ok)
        {
            if (ok->getFullSequence() != "")
            {
                if (UniquePSMs.find(ok->getFullSequence()) != UniquePSMs.end())
                {
                    auto current = UniquePSMs[ok->getFullSequence()];
                    if (std::get<2>(current)->getScore() < ok->getScore())
                    {
                        UniquePSMs[ok->getFullSequence()] = std::make_tuple(ok->getBaseSequence(),
                                                                            ok->getFullSequence(), ok);
                    }
                }
                else
                {
                    UniquePSMs.emplace(ok->getFullSequence(), std::make_tuple(ok->getBaseSequence(),
                                                                              ok->getFullSequence(), ok));
                }
            }
        }
    }
}
