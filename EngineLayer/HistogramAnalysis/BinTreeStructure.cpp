#include "BinTreeStructure.h"
#include "Bin.h"
#include "../PeptideSpectralMatch.h"
#include "OkBin.h"
#include "MyInfo.h"

#include "Math.h"

using namespace MassSpectrometry;
using namespace Proteomics::Fragmentation;
namespace EngineLayer
{
    namespace HistogramAnalysis
    {
        
        std::vector<Bin*> BinTreeStructure::getFinalBins() const
        {
            return privateFinalBins;
        }
        
        void BinTreeStructure::setFinalBins(const std::vector<Bin*> &value)
        {
            privateFinalBins = value;
        }
        
        void BinTreeStructure::GenerateBins(std::vector<PeptideSpectralMatch*> &targetAndDecoyMatches, double dc)
        {
#ifdef ORIG
            std::vector<double> listOfMassShifts = targetAndDecoyMatches.Where([&] (std::any b)   {
                    b::PeptideMonisotopicMass.HasValue;
                })->Select([&] (std::any b)  {
                        return b::ScanPrecursorMass - b::PeptideMonisotopicMass->Value;
                    }).OrderBy([&] (std::any b)  {
                            return b;
                        }).ToList();
#endif
            std::vector<double> listOfMassShifts;
            for ( auto b: targetAndDecoyMatches ) {
                if ( b->getPeptideMonisotopicMass().has_value() ) {
                    listOfMassShifts.push_back(b->getScanPrecursorMass() - b->getPeptideMonisotopicMass().value());
                }
            }
            std::sort(listOfMassShifts.begin(), listOfMassShifts.end() );
            
            double minMassShift = *(std::min(listOfMassShifts.begin(), listOfMassShifts.end()));
            double maxMassShift = *(std::max(listOfMassShifts.begin(), listOfMassShifts.end()));
            
            std::vector<int> p(listOfMassShifts.size());
            
            int firstIndex = 0;
            int lastIndex = 0;
            for (int i = 0; i < (int)listOfMassShifts.size(); i++)
            {
                auto thisMassShift = listOfMassShifts[i];
                
                while (thisMassShift - listOfMassShifts[firstIndex] > dc)
                {
                    firstIndex++;
                }
                while (lastIndex + 1 < (int)listOfMassShifts.size() &&
                       listOfMassShifts[lastIndex + 1] - thisMassShift <= dc)
                {
                    lastIndex++;
                }
                
                p[i] = lastIndex - firstIndex;
            }
            
            int maxP = *(std::max(p.begin(), p.end()));
            std::vector<double> sigma(listOfMassShifts.size());
            
            for (int i = 0; i < (int)listOfMassShifts.size(); i++)
            {
                auto thisMassShift = listOfMassShifts[i];
                auto thisP = p[i];
                if (thisP == maxP)
                {
                    sigma[i] = std::max(maxMassShift - thisMassShift, thisMassShift - minMassShift);
                }
                else
                {
                    // SIGMA IS THE DISTANCE TO THE CLOSEST MASS SHIFT THAT HAS A HIGHER P VALUE THAN ITSELF
                    
                    sigma[i] = GetSigma(thisMassShift, thisP, i, listOfMassShifts, p);
                }
            }
            
            std::vector<OkBin*> listokbin;
            for (int i = 0; i < (int)sigma.size(); i++)
            {
                OkBin tempVar(listOfMassShifts[i], sigma[i], p[i]);
                listokbin.push_back(&tempVar);
            }
            
            std::unordered_set<double> prelimBins;
#ifdef ORIG
            //for (OkBin *okbin : listokbin.OrderByDescending([&] (std::any b)    {
            //            b::P;
            //        }))
#endif
            std::sort(listokbin.begin(), listokbin.end(), [&] (OkBin* l, OkBin *r) {
                    return l->P > r->P;
                });
            for (OkBin *okbin : listokbin )
            {
                if (okbin->Sigma < dc || okbin->P < MinAdditionalPsmsInBin)
                {
                    continue;
                }
                bool add = true;
                for (auto a : prelimBins)
                {
                    if (std::abs(okbin->MassShift - a) <= dc)
                    {
                        add = false;
                        break;
                    }
                }
                if (add)
                {
                    prelimBins.insert(okbin->MassShift);
                }
            }
            
            auto forFinalBins = std::unordered_map<double, std::vector<double>>();
            for (auto ok : prelimBins)
            {
                forFinalBins.emplace(ok, std::vector<double>());
            }
            for (auto a : listOfMassShifts)
            {
                for (auto b : prelimBins)
                {
                    if (std::abs(a - b) <= dc)
                    {
                        forFinalBins[b].push_back(a);
                    }
                }
            }
            
#ifdef ORIG
            setFinalBins(forFinalBins.Select([&] (std::any b) {
                        new Bin(b->Value->Average());
                    }).ToList());
#endif
            std::vector<Bin*> newBins;
            for ( auto b: forFinalBins ) {
                double sum=0.0;
                for ( auto p: b.second) {
                    sum += p;
                }
                newBins.push_back(new Bin(sum/b.second.size()));
            }
            setFinalBins(newBins);
            
            for (int i = 0; i < (int)targetAndDecoyMatches.size(); i++)
            {
                for (auto bin : getFinalBins())
                {
                    if (targetAndDecoyMatches[i]->getPeptideMonisotopicMass().has_value() &&
                        std::abs(targetAndDecoyMatches[i]->getScanPrecursorMass() -
                                 targetAndDecoyMatches[i]->getPeptideMonisotopicMass().value() -
                                 bin->getMassShift()) <= dc)
                    {
                        bin->Add(targetAndDecoyMatches[i]);
                    }
                }
            }
            
#ifdef ORIG
            setFinalBins(getFinalBins().Where([&] (std::any b)   {
                        return b->Count > 1;
                    }).ToList());
#endif
            newBins.clear();
            for ( auto b : getFinalBins() ) {
                if ( b->getCount() > 1 ) {
                    newBins.push_back(b);
                }
            }
            setFinalBins(newBins);
            
            IdentifyUnimodBins(dc);
            IdentifyUniprotBins(dc);
            IdentifyAA(dc);
            
            IdentifyCombos(dc);
            
            IdentifyResidues();
            
            IdentifyMods();
            
            IdentifyAAsInCommon();
            
            IdentifyMine(dc);
            
            IdentifyFracWithSingle();
            IdentifyMedianLength();
        }
        
        double BinTreeStructure::GetSigma(double thisMassShift, int thisP, int i, std::vector<double> &listOfMassShifts,
                                          std::vector<int> &p)
        {
            int currentDown = i - 1;
            int currentUp = i + 1;
            int lookingAtP;
            double distDown, distUp;
            while (true)
            {
                distDown = currentDown == -1 ? std::numeric_limits<double>::max() :
                    thisMassShift - listOfMassShifts[currentDown];
                distUp = currentUp == (int)listOfMassShifts.size() ? std::numeric_limits<double>::max() :
                    listOfMassShifts[currentUp] - thisMassShift;
                if (distDown < distUp)
                {
                    lookingAtP = p[currentDown];
                    if (lookingAtP > thisP)
                    {
                        return distDown;
                    }
                    currentDown--;
                }
                else
                {
                    lookingAtP = p[currentUp];
                    if (lookingAtP > thisP)
                    {
                        return distUp;
                    }
                    currentUp++;
                }
            }
        }
        
        void BinTreeStructure::IdentifyFracWithSingle()
        {
            for (auto bin : getFinalBins())
            {
#ifdef ORIG
                auto numTarget = bin->UniquePSMs.Values->Count([&] (std::any b)  {
                        !b::Item3->IsDecoy;
                    });
#endif
                int numTarget = 0;
                for ( auto b: bin->UniquePSMs ) {
                    if ( !std::get<2>(b.second)->getIsDecoy() ) {
                        numTarget++;
                    }
                }
                if (numTarget > 0)
                {
#ifdef ORIG
                    bin->setFracWithSingle(static_cast<double>(bin->UniquePSMs.Values->Count([&] (std::any b) {
                                    return !b::Item3->IsDecoy && b::Item3->NumDifferentMatchingPeptides == 1;
                                })) / numTarget);
#endif
                    int count=0;
                    for ( auto b: bin->UniquePSMs) {
                        if ( !std::get<2>(b.second)->getIsDecoy() &&
                             std::get<2>(b.second)->getNumDifferentMatchingPeptides() == 1 ) {
                            count++;
                        }
                    }
                    bin->setFracWithSingle(static_cast<double>(count/numTarget));
                }
            }
        }
        
        void BinTreeStructure::IdentifyMedianLength()
        {
            for (auto bin : getFinalBins())
            {
#ifdef ORIG
                auto numTarget = bin->UniquePSMs.Values->Count([&] (std::any b) {
                        !b::Item3->IsDecoy;
                    });
#endif
                int numTarget = 0;
                for ( auto b: bin->UniquePSMs ) {
                    if ( !std::get<2>(b.second)->getIsDecoy() ) {
                        numTarget++;
                    }
                }
                
                if (numTarget > 0)
                {
#ifdef ORIG
                    bin->setMedianLength(Statistics::Median(bin->UniquePSMs.Values->Where([&] (std::any b) {
                                    !b::Item3->IsDecoy;
                                }).Where([&] (std::any b)  {
                                        b::Item3->PeptideLength.HasValue;
                                    })->Select([&] (std::any b)  {
                                            (double)b.Item3.PeptideLength.Value;
                                        })));
#endif
                    std::vector<double>tmpvec;
                    for ( auto b: bin->UniquePSMs ) {
                        auto q = b.second;
                        if ( std::get<2>(q)->getPeptideLength().has_value() ) {
                            tmpvec.push_back(std::get<2>(q)->getPeptideLength().value() );
                        }
                    }
                    bin->setMedianLength(Math::Median(tmpvec));                    
                }
            }
        }
        
        void BinTreeStructure::IdentifyAAsInCommon()
        {
            for (auto bin : getFinalBins())
            {
                bin->setAAsInCommon(std::unordered_map<char, int>());
#ifdef ORIG
                //for (auto hehe : bin->UniquePSMs.Values->Where([&] (std::any b)    {
                //            !b::Item3->IsDecoy;
                // }))
#endif
                for ( auto h : bin->UniquePSMs )                    
                {
                    auto hehe = h.second;
                    if ( std::get<2>(hehe)->getIsDecoy() ) {
                        continue;
                    }
                    
                    std::unordered_set<char> chars;
                    for (int i = 0; i < (int)std::get<0>(hehe).size(); i++)
                    {
                        chars.insert(std::get<0>(hehe)[i]);
                    }
                    for (auto ch : chars)
                    {
                        if (bin->getAAsInCommon().find(ch) != bin->getAAsInCommon().end())
                        {
                            bin->getAAsInCommon()[ch]++;
                        }
                        else
                        {
                            bin->getAAsInCommon().emplace(ch, 1);
                        }
                    }
                }
            }
        }
        
        void BinTreeStructure::IdentifyMods()
        {
            for (auto bin : getFinalBins())
            {
                bin->ModsInCommon = std::unordered_map<std::string, int>();
#ifdef ORIG
                //for (auto hehe : bin->UniquePSMs.Values->Where([&] (std::any b) {
                //            !b::Item3->IsDecoy;
                //        }))
#endif
                for ( auto h : bin->UniquePSMs ) 
                {
                    auto hehe = h.second;
                    if ( std::get<2>(hehe)->getIsDecoy() ) {
                        continue;
                    }
                    
                    int inModLevel = 0;
                    StringBuilder *currentMod = new StringBuilder();
                    std::unordered_set<std::string> modsHere;
                    for (int i = 0; i < (int)std::get<1>(hehe).size(); i++)
                    {
                        char ye = std::get<1>(hehe)[i];
                        if (ye == '[')
                        {
                            inModLevel++;
                            if (inModLevel == 1)
                            {
                                continue;
                            }
                        }
                        else if (ye == ']' )
                        {
                            inModLevel--;
                            if (inModLevel == 0)
                            {
                                if (!StringHelper::startsWith(currentMod->toString(),"Common Fixed:") )
                                {
                                    modsHere.insert(currentMod->toString());
                                }
                                currentMod->clear();
                            }
                            continue;
                        }
                        if (inModLevel > 0)
                        {
                            currentMod->append(ye);
                        }
                    }
                    for (auto modInHS : modsHere)
                    {
                        if (bin->ModsInCommon.find(modInHS) != bin->ModsInCommon.end())
                        {
                            bin->ModsInCommon[modInHS]++;
                        }
                        else
                        {
                            bin->ModsInCommon.emplace(modInHS, 1);
                        }
                    }
                    
                    delete currentMod;
                }
            }
        }
        
        void BinTreeStructure::IdentifyResidues()
        {
            for (auto bin : getFinalBins())
            {
                bin->IdentifyResidues();
            }
        }
        
        void BinTreeStructure::IdentifyUnimodBins(double v)
        {
            for (auto bin : getFinalBins())
            {
                bin->IdentifyUnimodBins(v);
            }
        }
        
        void BinTreeStructure::IdentifyUniprotBins(double v)
        {
            for (auto bin : getFinalBins())
            {
                bin->IdentifyUniprotBins(v);
            }
        }
        
        void BinTreeStructure::IdentifyCombos(double v)
        {
#ifdef ORIG
            double totalTargetCount = getFinalBins().Select([&] (std::any b) {
                    b::CountTarget;
                }).Sum();
#endif
            double totalTargetCount =0;
            for ( auto b:  getFinalBins() ) {
                totalTargetCount += b->getCountTarget();
            }
            //auto ok = std::unordered_set<std::tuple<double, double, double>>();
            BinTuple_set ok;
            
            // For every non-zero bin
#ifdef ORIG
            //for (auto bin : getFinalBins().Where([&] (std::any b) {
            //            return std::abs(b->getMassShift()) > v;
            //        }))
#endif
            for (auto bin : getFinalBins() )
            {
                if ( std::abs(bin->getMassShift()) <= v ) {
                    continue;
                }
                        
#ifdef ORIG
                //for (auto bin2 : getFinalBins().Where([&] (std::any b) {
                //            return std::abs(b->getMassShift()) > v;
                //        }))
#endif
                for ( auto bin2: getFinalBins() ) 
                {
                    if ( std::abs(bin2->getMassShift()) <= v ) {
                        continue;
                    }
                    if (bin->getCountTarget() * bin2->getCountTarget() >= totalTargetCount)
                    {
                        ok.insert(std::make_tuple(bin->getMassShift(), bin2->getMassShift(),
                                                  std::min(bin->getCountTarget(), bin2->getCountTarget())));
                    }
                }
            }
            
            for (auto bin : getFinalBins())
            {
                bin->IdentifyCombos(v, ok);
            }
        }
        
        void BinTreeStructure::IdentifyAA(double v)
        {
            for (auto bin : getFinalBins())
            {
                bin->IdentifyAA(v);
            }
        }
        
        void BinTreeStructure::IdentifyMine(double v)
        {
            auto myInfos = std::vector<MyInfo*>
                {
                    new MyInfo(0, "Exact match!"),
                    new MyInfo(-48.128629, "Phosphorylation-Lysine: Probably reverse is the correct match"),
                    new MyInfo(-76.134779, "Phosphorylation-Arginine: Probably reverse is the correct match"),
                    new MyInfo(1.0029, "1 MM"),
                    new MyInfo(2.0052, "2 MM"),
                    new MyInfo(3.0077, "3 MM"),
                    new MyInfo(173.051055, "Acetylation + Methionine: Usually on protein N terminus"),
                    new MyInfo(-91.009185, "neg Carbamidomethylation - H2S: Usually on cysteine."),
                    new MyInfo(-32.008456, "oxidation and then loss of oxidized M side chain"),
                    new MyInfo(-79.966331, "neg Phosphorylation."),
                    new MyInfo(189.045969, "Carboxymethylated + Methionine. Usually on protein N terminus"),
                    new MyInfo(356.20596, "Lysine+V+E or Lysine+L+D"),
                    new MyInfo(239.126988, "Lysine+H(5) C(5) N O(2), possibly Nmethylmaleimide"),
                    new MyInfo(-105.02484, "Methionine loss then acetaldehyde"),
                    new MyInfo(52.911464, "Fe[III]"),
                    new MyInfo(71.000729, "H C2 N O2"),
                    new MyInfo(50.000394, "H2 O3")
                };
            for (auto bin : getFinalBins())
            {
                bin->setMine("");
                for (auto myInfo : myInfos)
                {
                    if (std::abs(myInfo->getMassShift() - bin->getMassShift()) <= v)
                    {
                        bin->setMine(myInfo->infostring);
                    }
                }
            }
        }
    }
}
