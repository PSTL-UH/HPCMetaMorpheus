#include "MzIdentMLWriter.h"
#include "../../EngineLayer/PeptideSpectralMatch.h"
#include "../../EngineLayer/ProteinParsimony/ProteinGroup.h"
#include "../../EngineLayer/GlobalVariables.h"

using namespace Chemistry;
using namespace EngineLayer;
using namespace MzLibUtil;
using namespace Proteomics;
using namespace Proteomics::ProteolyticDigestion;

namespace TaskLayer
{

	void MzIdentMLWriter::WriteMzIdentMl(std::vector<PeptideSpectralMatch*> &psms, std::vector<EngineLayer::ProteinGroup*> &groups, std::vector<Modification*> &variableMods, std::vector<Modification*> &fixedMods, std::vector<Protease*> &proteases, double qValueFilter, Tolerance *productTolerance, Tolerance *parentTolerance, int missedCleavages, const std::wstring &outputPath)
	{
		psms = psms.Where([&] (std::any p)
		{
			return p::FdrInfo::QValue <= qValueFilter && p::FdrInfo::QValueNotch <= qValueFilter;
		});

		std::vector<PeptideWithSetModifications*> peptides = psms.SelectMany([&] (std::any i)
		{
			i::BestMatchingPeptides->Select([&] (std::any v)
			{
				v::Peptide;
			});
		}).Distinct().ToList();
		std::vector<Protein*> proteins = peptides.Select([&] (std::any p)
		{
			p::Protein;
		}).Distinct().ToList();
		std::vector<std::wstring> filenames = psms.Select([&] (std::any i)
		{
			i::FullFilePath;
		}).Distinct().ToList();
		std::unordered_map<std::wstring, std::wstring> database_reference;
		std::vector<std::wstring> databases = proteins.Select([&] (std::any p)
		{
			p::DatabaseFilePath;
		}).Distinct().ToList();

		UTF8Encoding *utf8EmitBOM = new UTF8Encoding(false);
		XmlWriterSettings *settings = new XmlWriterSettings();
		settings->NewLineChars = L"\n";
		settings->Indent = true;
		settings->Encoding = utf8EmitBOM;
		XmlSerializer *_indexedSerializer = new XmlSerializer(mzIdentML110::Generated::MzIdentMLType110::typeid);
		auto _mzid = new mzIdentML110::Generated::MzIdentMLType110();
		_mzid->version = L"1.1.0";
		_mzid->id = L"";

		_mzid->Provider = new mzIdentML110::Generated::ProviderType();
		_mzid->Provider->id = L"PROVIDER";
		_mzid->Provider->ContactRole = new mzIdentML110::Generated::ContactRoleType();
		_mzid->Provider->ContactRole->contact_ref = L"UWMadisonSmithGroup";
		_mzid->Provider->ContactRole->Role = new mzIdentML110::Generated::RoleType();
		_mzid->Provider->ContactRole.Role->cvParam = new mzIdentML110::Generated::CVParamType();
		_mzid->Provider->ContactRole.Role->cvParam->accession = L"MS:1001271";
		_mzid->Provider->ContactRole.Role->cvParam->name = L"researcher";
		_mzid->Provider->ContactRole.Role->cvParam->cvRef = L"PSI-MS";

		_mzid->AuditCollection = std::vector<mzIdentML110::Generated::AbstractContactType*>(2);

		mzIdentML110::Generated::PersonType *tempVar = new mzIdentML110::Generated::PersonType();
		tempVar->id = L"UWMadisonSmithGroupPerson";
		mzIdentML110::Generated::CVParamType *tempVar2 = new mzIdentML110::Generated::CVParamType();
		tempVar2->accession = L"MS:1000589";
		tempVar2->name = L"contact email";
		tempVar2->cvRef = L"PSI-MS";
		tempVar2->value = L"mm_support@chem.wisc.edu";
		mzIdentML110::Generated::CVParamType *tempVar3 = new mzIdentML110::Generated::CVParamType();
		tempVar3->accession = L"MS:1000590";
		tempVar3->name = L"affiliation name";
		tempVar3->cvRef = L"PSI-MS";
		tempVar3->value = L"UWMadisonSmithGroup";
		tempVar->cvParam = {tempVar2, tempVar3};
		_mzid->AuditCollection[0] = tempVar;

		mzIdentML110::Generated::OrganizationType *tempVar4 = new mzIdentML110::Generated::OrganizationType();
		tempVar4->id = L"UWMadisonSmithGroup";
		mzIdentML110::Generated::CVParamType *tempVar5 = new mzIdentML110::Generated::CVParamType();
		tempVar5->accession = L"MS:1000589";
		tempVar5->name = L"contact email";
		tempVar5->cvRef = L"PSI-MS";
		tempVar5->value = L"mm_support@chem.wisc.edu";
		mzIdentML110::Generated::CVParamType *tempVar6 = new mzIdentML110::Generated::CVParamType();
		tempVar6->accession = L"MS:1000590";
		tempVar6->name = L"affiliation name";
		tempVar6->cvRef = L"PSI-MS";
		tempVar6->value = L"UWMadisonSmithGroup";
		tempVar4->cvParam = {tempVar5, tempVar6};
		_mzid->AuditCollection[1] = tempVar4;

		//cvlist: URLs of controlled vocabularies used within the file.
		mzIdentML110::Generated::cvType *tempVar7 = new mzIdentML110::Generated::cvType();
		tempVar7->id = L"PSI-MS";
		tempVar7->fullName = L"Proteomics Standards Initiative Mass Spectrometry Vocabularies";
		tempVar7->uri = L"https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo";
		tempVar7->version = L"4.0.9";
		mzIdentML110::Generated::cvType *tempVar8 = new mzIdentML110::Generated::cvType();
		tempVar8->id = L"PSI-MOD";
		tempVar8->fullName = L"Proteomics Standards Initiative Modification Vocabularies";
		tempVar8->uri = L"http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/mod/data/PSI-MOD.obo";
		tempVar8->version = L"1.2";
		mzIdentML110::Generated::cvType *tempVar9 = new mzIdentML110::Generated::cvType();
		tempVar9->id = L"UNIMOD";
		tempVar9->fullName = L"UNIT-ONTOLOGY";
		tempVar9->uri = L"http://www.unimod.org/obo/unimod.obo";
		mzIdentML110::Generated::cvType *tempVar10 = new mzIdentML110::Generated::cvType();
		tempVar10->id = L"UO";
		tempVar10->fullName = L"UNIT-ONTOLOGY";
		tempVar10->uri = L"http://www.unimod.org/obo/unimod.obo";
		_mzid->cvList = {tempVar7, tempVar8, tempVar9, tempVar10};

		mzIdentML110::Generated::AnalysisSoftwareType *tempVar11 = new mzIdentML110::Generated::AnalysisSoftwareType();
		tempVar11->id = L"AS_MetaMorpheus";
		tempVar11->name = L"MetaMorpheus";
		tempVar11->version = GlobalVariables::getMetaMorpheusVersion();
		tempVar11->uri = L"https://github.com/smith-chem-wisc/MetaMorpheus";
		tempVar11->SoftwareName = new mzIdentML110::Generated::ParamType();
		tempVar11->SoftwareName->Item = new mzIdentML110::Generated::CVParamType();
		tempVar11->SoftwareName.Item->accession = L"MS:1002826";
		tempVar11->SoftwareName.Item->name = L"MetaMorpheus";
		tempVar11->SoftwareName.Item->cvRef = L"PSI-MS";
		tempVar11->ContactRole = new mzIdentML110::Generated::ContactRoleType();
		tempVar11->ContactRole->contact_ref = L"UWMadisonSmithGroup";
		tempVar11->ContactRole->Role = new mzIdentML110::Generated::RoleType();
		tempVar11->ContactRole.Role->cvParam = new mzIdentML110::Generated::CVParamType();
		tempVar11->ContactRole.Role->cvParam->accession = L"MS:1001267";
		tempVar11->ContactRole.Role->cvParam->name = L"software vendor";
		tempVar11->ContactRole.Role->cvParam->cvRef = L"PSI-MS";
		_mzid->AnalysisSoftwareList = {tempVar11};
		_mzid->DataCollection = new mzIdentML110::Generated::DataCollectionType();
		_mzid->DataCollection->AnalysisData = new mzIdentML110::Generated::AnalysisDataType();
		mzIdentML110::Generated::SpectrumIdentificationListType *tempVar12 = new mzIdentML110::Generated::SpectrumIdentificationListType();
		tempVar12->id = L"SIL";
		tempVar12->SpectrumIdentificationResult = std::vector<mzIdentML110::Generated::SpectrumIdentificationResultType*>(psms.size()());
		_mzid->DataCollection.AnalysisData->SpectrumIdentificationList = {tempVar12};
		_mzid->DataCollection->Inputs = new mzIdentML110::Generated::InputsType();
		_mzid->DataCollection.Inputs->SearchDatabase = std::vector<mzIdentML110::Generated::SearchDatabaseType*>(databases.size()());
		_mzid->DataCollection.Inputs->SpectraData = std::vector<mzIdentML110::Generated::SpectraDataType*>(filenames.size());

		_mzid->SequenceCollection = new mzIdentML110::Generated::SequenceCollectionType();
		_mzid->SequenceCollection->Peptide = std::vector<mzIdentML110::Generated::PeptideType*>(peptides.size());
		_mzid->SequenceCollection->DBSequence = std::vector<mzIdentML110::Generated::DBSequenceType*>(proteins.size());
		_mzid->SequenceCollection->PeptideEvidence = std::vector<mzIdentML110::Generated::PeptideEvidenceType*>(peptides.size());

		_mzid->AnalysisCollection = new mzIdentML110::Generated::AnalysisCollectionType();
		mzIdentML110::Generated::SpectrumIdentificationType *tempVar13 = new mzIdentML110::Generated::SpectrumIdentificationType();
		tempVar13->id = L"SI";
		tempVar13->spectrumIdentificationList_ref = L"SIL";
		tempVar13->spectrumIdentificationProtocol_ref = L"SIP";
		tempVar13->InputSpectra = std::vector<mzIdentML110::Generated::InputSpectraType*>(filenames.size());
		tempVar13->SearchDatabaseRef = std::vector<mzIdentML110::Generated::SearchDatabaseRefType*>(databases.size());
		_mzid->AnalysisCollection->SpectrumIdentification = {tempVar13};
		int database_index = 0;
		for (auto database : databases)
		{
			mzIdentML110::Generated::SearchDatabaseType *tempVar14 = new mzIdentML110::Generated::SearchDatabaseType();
			tempVar14->id = L"SDB_" + std::to_wstring(database_index);
			tempVar14->location = database;
			tempVar14->DatabaseName = new mzIdentML110::Generated::ParamType();
			tempVar14->DatabaseName->Item = new mzIdentML110::Generated::CVParamType();
			tempVar14->DatabaseName.Item->accession = L"MS:1001073";
			tempVar14->DatabaseName.Item->name = L"database type amino acid";
			tempVar14->DatabaseName.Item->cvRef = L"PSI-MS";
			_mzid->DataCollection.Inputs.SearchDatabase[database_index] = tempVar14;
			database_reference.emplace(database, L"SDB_" + std::to_wstring(database_index));
			mzIdentML110::Generated::SearchDatabaseRefType *tempVar15 = new mzIdentML110::Generated::SearchDatabaseRefType();
			tempVar15->searchDatabase_ref = L"SDB_" + std::to_wstring(database_index);
			_mzid->AnalysisCollection.SpectrumIdentification[0].SearchDatabaseRef[database_index] = tempVar15;
			database_index++;
		}

		int protein_index = 0;
		for (auto protein : proteins)
		{
			mzIdentML110::Generated::DBSequenceType *tempVar16 = new mzIdentML110::Generated::DBSequenceType();
			tempVar16->id = L"DBS_" + protein->Accession;
			tempVar16->lengthSpecified = true;
			tempVar16->length = protein->Length;
			tempVar16->searchDatabase_ref = database_reference[protein->DatabaseFilePath];
			tempVar16->accession = protein->Accession;
			tempVar16->Seq = protein->BaseSequence;
			mzIdentML110::Generated::CVParamType *tempVar17 = new mzIdentML110::Generated::CVParamType();
			tempVar17->accession = L"MS:1001088";
			tempVar17->name = L"protein description";
			tempVar17->cvRef = L"PSI-MS";
			tempVar17->value = protein->FullDescription;
			tempVar16->cvParam = {tempVar17};
			tempVar16->name = protein->Name;
			_mzid->SequenceCollection.DBSequence[protein_index] = tempVar16;
			protein_index++;
		}

		std::unordered_map<std::wstring, int> spectral_ids; //key is datafile, value is datafile's id
		int spectra_data_id = 0;
		for (auto data_filepath : filenames)
		{
			bool thermoRawFile = Path::GetExtension(data_filepath) == L".raw";
			std::wstring spectral_data_id = L"SD_" + std::to_wstring(spectra_data_id);
			spectral_ids.emplace(data_filepath, spectra_data_id);
			mzIdentML110::Generated::InputSpectraType *tempVar18 = new mzIdentML110::Generated::InputSpectraType();
			tempVar18->spectraData_ref = spectral_data_id;
			_mzid->AnalysisCollection.SpectrumIdentification[0].InputSpectra[spectra_data_id] = tempVar18;
			mzIdentML110::Generated::SpectraDataType *tempVar19 = new mzIdentML110::Generated::SpectraDataType();
			tempVar19->id = spectral_data_id;
			tempVar19->name = Path::GetFileNameWithoutExtension(data_filepath);
			tempVar19->location = data_filepath;
			tempVar19->FileFormat = new mzIdentML110::Generated::FileFormatType();
			tempVar19->FileFormat->cvParam = new mzIdentML110::Generated::CVParamType();
			tempVar19->FileFormat.cvParam->accession = thermoRawFile ? L"MS:1000563" : L"MS:1000584";
			tempVar19->FileFormat.cvParam->name = thermoRawFile ? L"Thermo RAW format" : L"mzML format";
			tempVar19->FileFormat.cvParam->cvRef = L"PSI-MS";
			tempVar19->SpectrumIDFormat = new mzIdentML110::Generated::SpectrumIDFormatType();
			tempVar19->SpectrumIDFormat->cvParam = new mzIdentML110::Generated::CVParamType();
			tempVar19->SpectrumIDFormat.cvParam->accession = thermoRawFile ? L"MS:1000768" : L"MS:1001530";
			tempVar19->SpectrumIDFormat.cvParam->name = thermoRawFile ? L"Thermo nativeID format" : L"mzML unique identifier";
			tempVar19->SpectrumIDFormat.cvParam->cvRef = L"PSI-MS";
			_mzid->DataCollection.Inputs.SpectraData[spectra_data_id] = tempVar19;
			spectra_data_id++;
		}

		int sir_id = 0;
		int pe_index = 0;
		int p_index = 0;
		std::unordered_map<PeptideWithSetModifications*, int> peptide_evidence_ids;
		std::unordered_map<std::wstring, std::tuple<int, std::unordered_set<std::wstring>>> peptide_ids; //key is peptide sequence, value is <peptide id for that peptide, peptide evidences>, list of spectra id's
		std::unordered_map<std::tuple<std::wstring, int>, std::tuple<int, int>> psm_per_scan; //key is <filename, scan numer> value is <scan result id, scan item id #'s (could be more than one ID per scan)>

		auto unambiguousPsms = psms.Where([&] (std::any psm)
		{
		delete _mzid;
		delete _indexedSerializer;
		delete settings;
//C# TO C++ CONVERTER TODO TASK: A 'delete utf8EmitBOM' statement was not added since utf8EmitBOM was assigned to another object. Handle memory management manually.
			return psm::FullSequence != nullptr;
		});

		for (auto psm : unambiguousPsms)
		{
			for (PeptideWithSetModifications *peptide : psm->BestMatchingPeptides->Select([&] (std::any p)
			{
				p::Peptide;
			}).Distinct())
			{
				//if first peptide on list hasn't been added, add peptide and peptide evidence
				std::tuple<int, std::unordered_set<std::wstring>> peptide_id;
				std::unordered_map<std::wstring, std::tuple<int, std::unordered_set<std::wstring>>>::const_iterator peptide_ids_iterator = peptide_ids.find(peptide.FullSequence);
				if (peptide_ids_iterator == peptide_ids.end())
				{
					peptide_id = peptide_ids_iterator->second;
					peptide_id = std::tuple<int, std::unordered_set<std::wstring>>(p_index, std::unordered_set<std::wstring>());
					p_index++;
					mzIdentML110::Generated::PeptideType *tempVar20 = new mzIdentML110::Generated::PeptideType();
					tempVar20->PeptideSequence = peptide::BaseSequence;
					tempVar20->id = L"P_" + peptide_id::Item1;
					tempVar20->Modification = std::vector<mzIdentML110::Generated::ModificationType*>(peptide::NumMods);
					_mzid->SequenceCollection.Peptide[peptide_id::Item1] = tempVar20;
					int mod_id = 0;
					for (auto mod : peptide::AllModsOneIsNterminus)
					{
						mzIdentML110::Generated::ModificationType *tempVar21 = new mzIdentML110::Generated::ModificationType();
						tempVar21->location = mod.first - 1;
						tempVar21->locationSpecified = true;
						tempVar21->monoisotopicMassDelta = mod.second::MonoisotopicMass->Value;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
						tempVar21->residues = {peptide::BaseSequence[std::min(std::max(0, mod.first - 2), peptide->Length - 1)].ToString()};
						tempVar21->monoisotopicMassDeltaSpecified = true;
						tempVar21->cvParam = {GetUnimodCvParam(mod.second)};
						_mzid->SequenceCollection.Peptide[peptide_id::Item1].Modification[mod_id] = tempVar21;
						mod_id++;
					}
					peptide_ids.emplace(peptide::FullSequence, peptide_id);
				}
				else
				{
					peptide_id = peptide_ids_iterator->second;
				}

				if (peptide_evidence_ids.find(peptide) == peptide_evidence_ids.end())
				{
					mzIdentML110::Generated::PeptideEvidenceType *tempVar22 = new mzIdentML110::Generated::PeptideEvidenceType();
					tempVar22->id = L"PE_" + std::to_wstring(pe_index);
					tempVar22->peptide_ref = L"P_" + peptide_id::Item1;
					tempVar22->dBSequence_ref = L"DBS_" + peptide::Protein::Accession;
					tempVar22->isDecoy = peptide::Protein::IsDecoy;
					tempVar22->startSpecified = true;
					tempVar22->start = peptide::OneBasedStartResidueInProtein;
					tempVar22->endSpecified = true;
					tempVar22->end = peptide::OneBasedEndResidueInProtein;
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					tempVar22->pre = peptide::PreviousAminoAcid.ToString();
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					tempVar22->post = (peptide::OneBasedEndResidueInProtein < peptide::Protein::BaseSequence->Length) ? peptide::Protein[peptide::OneBasedEndResidueInProtein].ToString() : L"-";
					_mzid->SequenceCollection.PeptideEvidence[pe_index] = tempVar22;
					peptide_evidence_ids.emplace(peptide, pe_index);
					pe_index++;
				}
			}

			std::tuple<int, int> scan_result_scan_item;
			std::unordered_map<std::tuple<std::wstring, int>, std::tuple<int, int>>::const_iterator psm_per_scan_iterator = psm_per_scan.find(new std::tuple<std::wstring, int>(psm.FullFilePath, psm.ScanNumber));
			if (psm_per_scan_iterator == psm_per_scan.end()) //check to see if scan has already been added
			{
				scan_result_scan_item = psm_per_scan_iterator->second;
				scan_result_scan_item = std::tuple<int, int>(sir_id, 0);
				mzIdentML110::Generated::SpectrumIdentificationResultType *tempVar23 = new mzIdentML110::Generated::SpectrumIdentificationResultType();
				tempVar23->id = L"SIR_" + scan_result_scan_item::Item1;
				tempVar23->spectraData_ref = L"SD_" + std::to_wstring(spectral_ids[psm->getFullFilePath()]);
				tempVar23->spectrumID = L"scan=" + std::to_wstring(psm->getScanNumber());
				tempVar23->SpectrumIdentificationItem = std::vector<mzIdentML110::Generated::SpectrumIdentificationItemType*>(500);
				mzIdentML110::Generated::CVParamType *tempVar24 = new mzIdentML110::Generated::CVParamType();
				tempVar24->name = L"scan start time";
				tempVar24->cvRef = L"PSI-MS";
				tempVar24->accession = L"MS:1000016";
				tempVar24->value = std::to_wstring(psm->getScanRetentionTime());
				tempVar23->cvParam = {tempVar24};
				_mzid->DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item::Item1] = tempVar23;
				psm_per_scan.emplace(std::tuple<std::wstring, int>(psm->getFullFilePath(), psm->getScanNumber()), scan_result_scan_item);
				sir_id++;
			}
			else
			{
				scan_result_scan_item = psm_per_scan_iterator->second;
				psm_per_scan[std::tuple<std::wstring, int>(psm->getFullFilePath(), psm->getScanNumber())] = std::tuple<int, int>(scan_result_scan_item::Item1, scan_result_scan_item::Item2 + 1);
				scan_result_scan_item = psm_per_scan[std::tuple<std::wstring, int>(psm->getFullFilePath(), psm->getScanNumber())];
			}
			for (PeptideWithSetModifications *p : psm->BestMatchingPeptides->Select([&] (std::any p)
			{
				p->Peptide;
			}).Distinct())
			{
				peptide_ids[p::FullSequence].Item2->Add(L"SII_" + scan_result_scan_item::Item1 + L"_" + scan_result_scan_item::Item2);
			}
			mzIdentML110::Generated::CVParamType *tempVar25 = new mzIdentML110::Generated::CVParamType();
			tempVar25->name = L"MetaMorpheus:score";
			tempVar25->cvRef = L"PSI-MS";
			tempVar25->accession = L"MS:1002827";
			tempVar25->value = std::to_wstring(psm->getScore());
			mzIdentML110::Generated::CVParamType *tempVar26 = new mzIdentML110::Generated::CVParamType();
			tempVar26->accession = L"MS:1002354";
			tempVar26->name = L"PSM-level q-value";
			tempVar26->cvRef = L"PSI-MS";
			tempVar26->value = std::to_wstring(psm->getFdrInfo()->getQValue());
			_mzid->DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item::Item1].SpectrumIdentificationItem[scan_result_scan_item::Item2] = new mzIdentML110::Generated::SpectrumIdentificationItemType() { rank = 1, chargeState = psm->getScanPrecursorCharge(), id = L"SII_" + scan_result_scan_item::Item1 + L"_" + scan_result_scan_item::Item2, experimentalMassToCharge = std::round(psm->getScanPrecursorMonoisotopicPeakMz() * std::pow(10, 5)) / std::pow(10, 5), passThreshold = psm->getFdrInfo()->QValue <= 0.01, peptide_ref = L"P_" + std::get<0>(peptide_ids[psm->getFullSequence()]), PeptideEvidenceRef = std::vector<mzIdentML110::Generated::PeptideEvidenceRefType*>(psm->BestMatchingPeptides->Select([&] (std::any p)
			{
				p::Peptide;
			}).Distinct()->Count()), cvParam = {tempVar25, tempVar26}
			};
			if (psm->getPeptideMonisotopicMass().HasValue)
			{
				_mzid->DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item::Item1].SpectrumIdentificationItem[scan_result_scan_item::Item2]->calculatedMassToCharge = std::round(psm->getPeptideMonisotopicMass().Value.ToMz(psm->getScanPrecursorCharge()) * std::pow(10, 5)) / std::pow(10, 5);
				_mzid->DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item::Item1].SpectrumIdentificationItem[scan_result_scan_item::Item2]->calculatedMassToChargeSpecified = true;
			}

			int pe = 0;
			for (PeptideWithSetModifications *p : psm->BestMatchingPeptides->Select([&] (std::any p)
			{
				p->Peptide;
			}).Distinct())
			{
				mzIdentML110::Generated::PeptideEvidenceRefType *tempVar27 = new mzIdentML110::Generated::PeptideEvidenceRefType();
				tempVar27->peptideEvidence_ref = L"PE_" + std::to_wstring(peptide_evidence_ids[p]);
				_mzid->DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[scan_result_scan_item::Item1].SpectrumIdentificationItem[scan_result_scan_item::Item2].PeptideEvidenceRef[pe] = tempVar27;
				pe++;
			}
		}

		_mzid->AnalysisProtocolCollection = new mzIdentML110::Generated::AnalysisProtocolCollectionType();
		mzIdentML110::Generated::SpectrumIdentificationProtocolType *tempVar28 = new mzIdentML110::Generated::SpectrumIdentificationProtocolType();
		tempVar28->id = L"SIP";
		tempVar28->analysisSoftware_ref = L"AS_MetaMorpheus";
		tempVar28->SearchType = new mzIdentML110::Generated::ParamType();
		tempVar28->SearchType->Item = new mzIdentML110::Generated::CVParamType();
		tempVar28->SearchType.Item->accession = L"MS:1001083";
		tempVar28->SearchType.Item->name = L"ms-ms search";
		tempVar28->SearchType.Item->cvRef = L"PSI-MS";
		tempVar28->AdditionalSearchParams = new mzIdentML110::Generated::ParamListType();
		mzIdentML110::Generated::CVParamType *tempVar29 = new mzIdentML110::Generated::CVParamType();
		tempVar29->accession = L"MS:1001211";
		tempVar29->cvRef = L"PSI-MS";
		tempVar29->name = L"parent mass type mono";
		mzIdentML110::Generated::CVParamType *tempVar30 = new mzIdentML110::Generated::CVParamType();
		tempVar30->accession = L"MS:1001255";
		tempVar30->name = L"fragment mass type mono";
		tempVar30->cvRef = L"PSI-MS";
		tempVar28->AdditionalSearchParams->Items = {tempVar29, tempVar30};
		tempVar28->ModificationParams = std::vector<mzIdentML110::Generated::SearchModificationType*>(fixedMods.size() + variableMods.size());
		tempVar28->Enzymes = new mzIdentML110::Generated::EnzymesType();
		tempVar28->Enzymes->Enzyme = std::vector<mzIdentML110::Generated::EnzymeType*>(proteases.size());
		mzIdentML110::Generated::CVParamType *tempVar31 = new mzIdentML110::Generated::CVParamType();
		tempVar31->accession = L"MS:1001412";
		tempVar31->name = L"search tolerance plus value";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		tempVar31->value = productTolerance->Value->ToString();
		tempVar31->cvRef = L"PSI-MS";
		tempVar31->unitAccession = dynamic_cast<PpmTolerance*>(productTolerance) != nullptr? L"UO:0000169": L"UO:0000221";
		tempVar31->unitName = dynamic_cast<PpmTolerance*>(productTolerance) != nullptr? L"parts per million" : L"dalton";
		tempVar31->unitCvRef = L"UO";
		mzIdentML110::Generated::CVParamType *tempVar32 = new mzIdentML110::Generated::CVParamType();
		tempVar32->accession = L"MS:1001413";
		tempVar32->name = L"search tolerance minus value";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		tempVar32->value = productTolerance->Value->ToString();
		tempVar32->cvRef = L"PSI-MS";
		tempVar32->unitAccession = dynamic_cast<PpmTolerance*>(productTolerance) != nullptr? L"UO:0000169": L"UO:0000221";
		tempVar32->unitName = dynamic_cast<PpmTolerance*>(productTolerance) != nullptr? L"parts per million" : L"dalton";
		tempVar32->unitCvRef = L"UO";
		tempVar28->FragmentTolerance = {tempVar31, tempVar32};
		mzIdentML110::Generated::CVParamType *tempVar33 = new mzIdentML110::Generated::CVParamType();
		tempVar33->accession = L"MS:1001412";
		tempVar33->name = L"search tolerance plus value";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		tempVar33->value = parentTolerance->Value->ToString();
		tempVar33->cvRef = L"PSI-MS";
		tempVar33->unitAccession = dynamic_cast<PpmTolerance*>(parentTolerance) != nullptr? L"UO:0000169": L"UO:0000221";
		tempVar33->unitName = dynamic_cast<PpmTolerance*>(parentTolerance) != nullptr? L"parts per million" : L"dalton";
		tempVar33->unitCvRef = L"UO";
		mzIdentML110::Generated::CVParamType *tempVar34 = new mzIdentML110::Generated::CVParamType();
		tempVar34->accession = L"MS:1001413";
		tempVar34->name = L"search tolerance minus value";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
		tempVar34->value = parentTolerance->Value->ToString();
		tempVar34->cvRef = L"PSI-MS";
		tempVar34->unitAccession = dynamic_cast<PpmTolerance*>(parentTolerance) != nullptr? L"UO:0000169": L"UO:0000221";
		tempVar34->unitName = dynamic_cast<PpmTolerance*>(parentTolerance) != nullptr? L"parts per million" : L"dalton";
		tempVar34->unitCvRef = L"UO";
		tempVar28->ParentTolerance = {tempVar33, tempVar34};
		tempVar28->Threshold = new mzIdentML110::Generated::ParamListType();
		mzIdentML110::Generated::CVParamType *tempVar35 = new mzIdentML110::Generated::CVParamType();
		tempVar35->accession = L"MS:1001448";
		tempVar35->name = L"pep:FDR threshold";
		tempVar35->cvRef = L"PSI-MS";
		tempVar35->value = L"0.01";
		tempVar28->Threshold->Items = {tempVar35};
		_mzid->AnalysisProtocolCollection->SpectrumIdentificationProtocol = {tempVar28};

		int protease_index = 0;
		for (auto protease : proteases)
		{
			mzIdentML110::Generated::EnzymeType *tempVar36 = new mzIdentML110::Generated::EnzymeType();
			tempVar36->id = L"E_" + std::to_wstring(protease_index);
			tempVar36->name = protease->Name;
			tempVar36->semiSpecific = protease->CleavageSpecificity == CleavageSpecificity::Semi;
			tempVar36->missedCleavagesSpecified = true;
			tempVar36->missedCleavages = missedCleavages;
			tempVar36->EnzymeName = new mzIdentML110::Generated::ParamListType();
			mzIdentML110::Generated::CVParamType *tempVar37 = new mzIdentML110::Generated::CVParamType();
			tempVar37->accession = protease->PsiMsAccessionNumber;
			tempVar37->name = protease->PsiMsName;
			tempVar37->cvRef = L"PSI-MS";
			tempVar36->EnzymeName->Items = {tempVar37};
			_mzid->AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].Enzymes.Enzyme[protease_index] = tempVar36;
			protease_index++;
		}

		int mod_index = 0;
		for (auto mod : fixedMods)
		{
			mzIdentML110::Generated::SearchModificationType *tempVar38 = new mzIdentML110::Generated::SearchModificationType();
			tempVar38->fixedMod = true;
			tempVar38->massDelta = static_cast<float>(mod->MonoisotopicMass);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			tempVar38->residues = mod->Target->ToString();
			tempVar38->cvParam = {GetUnimodCvParam(mod)};
			_mzid->AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ModificationParams[mod_index] = tempVar38;
			mod_index++;
		}

		for (auto mod : variableMods)
		{
			mzIdentML110::Generated::SearchModificationType *tempVar39 = new mzIdentML110::Generated::SearchModificationType();
			tempVar39->fixedMod = false;
			tempVar39->massDelta = static_cast<float>(mod->MonoisotopicMass);
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
			tempVar39->residues = mod->Target->ToString();
			tempVar39->cvParam = {GetUnimodCvParam(mod)};
			_mzid->AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ModificationParams[mod_index] = tempVar39;
			mod_index++;
		}

		_mzid->AnalysisProtocolCollection->ProteinDetectionProtocol = new mzIdentML110::Generated::ProteinDetectionProtocolType();
		_mzid->AnalysisProtocolCollection.ProteinDetectionProtocol->id = L"PDP";
		_mzid->AnalysisProtocolCollection.ProteinDetectionProtocol->analysisSoftware_ref = L"AS_MetaMorpheus";
		_mzid->AnalysisProtocolCollection.ProteinDetectionProtocol->Threshold = new mzIdentML110::Generated::ParamListType();
		mzIdentML110::Generated::CVParamType *tempVar40 = new mzIdentML110::Generated::CVParamType();
		tempVar40->accession = L"MS:1001447";
		tempVar40->name = L"prot:FDR threshold";
		tempVar40->cvRef = L"PSI-MS";
		tempVar40->value = L"0.01";
		_mzid->AnalysisProtocolCollection.ProteinDetectionProtocol.Threshold->Items = {tempVar40};

		if (groups.size() > 0)
		{
			_mzid->DataCollection.AnalysisData->ProteinDetectionList = new mzIdentML110::Generated::ProteinDetectionListType();
			_mzid->DataCollection.AnalysisData.ProteinDetectionList->id = L"PDL";
			_mzid->DataCollection.AnalysisData.ProteinDetectionList->ProteinAmbiguityGroup = std::vector<mzIdentML110::Generated::ProteinAmbiguityGroupType*>(groups.size());

			int group_id = 0;
			int protein_id = 0;
			for (auto proteinGroup : groups)
			{
				mzIdentML110::Generated::ProteinAmbiguityGroupType *tempVar41 = new mzIdentML110::Generated::ProteinAmbiguityGroupType();
				tempVar41->id = L"PAG_" + std::to_wstring(group_id);
				tempVar41->ProteinDetectionHypothesis = std::vector<mzIdentML110::Generated::ProteinDetectionHypothesisType*>(proteinGroup->getProteins().size());
				_mzid->DataCollection.AnalysisData.ProteinDetectionList.ProteinAmbiguityGroup[group_id] = tempVar41;
				int pag_protein_index = 0;
				for (auto protein : proteinGroup->getProteins())
				{
					mzIdentML110::Generated::ProteinDetectionHypothesisType *tempVar42 = new mzIdentML110::Generated::ProteinDetectionHypothesisType();
					tempVar42->id = L"PDH_" + std::to_wstring(protein_id);
					tempVar42->dBSequence_ref = L"DBS_" + protein->Accession;
					tempVar42->passThreshold = proteinGroup->getQValue() <= 0.01;
					tempVar42->PeptideHypothesis = std::vector<mzIdentML110::Generated::PeptideHypothesisType*>(proteinGroup->getAllPeptides().size());
					mzIdentML110::Generated::CVParamType *tempVar43 = new mzIdentML110::Generated::CVParamType();
					tempVar43->accession = L"MS:1002828";
					tempVar43->name = L"MetaMorpheus:protein score";
					tempVar43->cvRef = L"PSI-MS";
					tempVar43->value = std::to_wstring(proteinGroup->getProteinGroupScore());
					mzIdentML110::Generated::CVParamType *tempVar44 = new mzIdentML110::Generated::CVParamType();
					tempVar44->accession = L"MS:1002373";
					tempVar44->name = L"protein group-level q-value";
					tempVar44->cvRef = L"PSI-MS";
					tempVar44->value = std::to_wstring(proteinGroup->getQValue());
					mzIdentML110::Generated::CVParamType *tempVar45 = new mzIdentML110::Generated::CVParamType();
					tempVar45->accession = L"MS:1001093";
					tempVar45->name = L"sequence coverage";
					tempVar45->cvRef = L"PSI-MS";
//C# TO C++ CONVERTER TODO TASK: There is no native C++ equivalent to 'ToString':
					tempVar45->value = proteinGroup->getSequenceCoveragePercent().front().ToString();
					mzIdentML110::Generated::CVParamType *tempVar46 = new mzIdentML110::Generated::CVParamType();
					tempVar46->accession = L"MS:1001097";
					tempVar46->name = L"distinct peptide sequences";
					tempVar46->cvRef = L"PSI-MS";
					tempVar46->value = std::to_wstring(proteinGroup->getUniquePeptides().size());
					tempVar42->cvParam = {tempVar43, tempVar44, tempVar45, tempVar46};
					_mzid->DataCollection.AnalysisData.ProteinDetectionList.ProteinAmbiguityGroup[group_id].ProteinDetectionHypothesis[pag_protein_index] = tempVar42;
					int peptide_id = 0;
					for (auto peptide : proteinGroup->getAllPeptides())
					{
						if (peptide_evidence_ids.find(peptide) != peptide_evidence_ids.end())
						{
							if (peptide->Protein == protein)
							{
								mzIdentML110::Generated::PeptideHypothesisType *tempVar47 = new mzIdentML110::Generated::PeptideHypothesisType();
								tempVar47->peptideEvidence_ref = L"PE_" + std::to_wstring(peptide_evidence_ids[peptide]);
								tempVar47->SpectrumIdentificationItemRef = std::vector<mzIdentML110::Generated::SpectrumIdentificationItemRefType*>(peptide_ids[peptide->FullSequence].Item2->Count);
								_mzid->DataCollection.AnalysisData.ProteinDetectionList.ProteinAmbiguityGroup[group_id].ProteinDetectionHypothesis[pag_protein_index].PeptideHypothesis[peptide_id] = tempVar47;

								int i = 0;
								for (std::wstring sii : std::get<1>(peptide_ids[peptide->FullSequence]))
								{
									mzIdentML110::Generated::SpectrumIdentificationItemRefType *tempVar48 = new mzIdentML110::Generated::SpectrumIdentificationItemRefType();
									tempVar48->spectrumIdentificationItem_ref = sii;
									_mzid->DataCollection.AnalysisData.ProteinDetectionList.ProteinAmbiguityGroup[group_id].ProteinDetectionHypothesis[pag_protein_index].PeptideHypothesis[peptide_id].SpectrumIdentificationItemRef[i] = tempVar48;
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
		XmlWriter *writer = XmlWriter::Create(outputPath, settings);
		_indexedSerializer->Serialize(writer, _mzid);
		writer->Close();

//C# TO C++ CONVERTER TODO TASK: A 'delete _mzid' statement was not added since _mzid was passed to a method or constructor. Handle memory management manually.
		delete _indexedSerializer;
//C# TO C++ CONVERTER TODO TASK: A 'delete settings' statement was not added since settings was passed to a method or constructor. Handle memory management manually.
//C# TO C++ CONVERTER TODO TASK: A 'delete utf8EmitBOM' statement was not added since utf8EmitBOM was assigned to another object. Handle memory management manually.
	}

	mzIdentML110::Generated::CVParamType *MzIdentMLWriter::GetUnimodCvParam(Modification *mod)
	{
		if (mod->DatabaseReference != nullptr && mod->DatabaseReference->ContainsKey(L"Unimod"))
		{
			mzIdentML110::Generated::CVParamType *tempVar = new mzIdentML110::Generated::CVParamType();
			tempVar->accession = L"UNIMOD:" + mod->DatabaseReference[L"Unimod"].First();
			tempVar->name = mod->IdWithMotif;
			tempVar->cvRef = L"PSI-MS";
			return tempVar;
		}
		else
		{
			mzIdentML110::Generated::CVParamType *tempVar2 = new mzIdentML110::Generated::CVParamType();
			tempVar2->accession = L"MS:1001460";
			tempVar2->name = L"unknown modification";
			tempVar2->cvRef = L"UNIMOD";
			tempVar2->value = mod->IdWithMotif;
			return tempVar2;
		}
	}
}
