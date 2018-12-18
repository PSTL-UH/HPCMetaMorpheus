﻿#pragma once

#include "../MetaMorpheusEngine.h"
#include <string>
#include <vector>

//C# TO C++ CONVERTER NOTE: Forward class declarations:
namespace EngineLayer { class PeptideSpectralMatch; }
namespace EngineLayer { class CommonParameters; }
namespace EngineLayer { class MetaMorpheusEngineResults; }

using namespace MassSpectrometry;
using namespace Proteomics::Fragmentation;
using namespace Proteomics::ProteolyticDigestion;

namespace EngineLayer
{
	namespace Localization
	{
		/// <summary>
		/// The mass difference between an experimental precursor and the theoretical mass of the assigned peptide is determined. The LocalizationEngine attempts
		/// to localize this mass to one of the residues. It does this by adding the mass difference to each theoretical ion mass and looking for additional matches
		/// in the experimental spectrum. This engine should only be run for open, notch or custom searches. It should not be run for exact mass or missed
		/// monoisopic searches.
		/// </summary>
		class LocalizationEngine : public MetaMorpheusEngine
		{
		private:
			const std::vector<PeptideSpectralMatch*> AllResultingIdentifications;
			MsDataFile *const MyMsDataFile;

		public:
			virtual ~LocalizationEngine()
			{
				delete MyMsDataFile;
			}

			LocalizationEngine(std::vector<PeptideSpectralMatch*> &allResultingIdentifications, MsDataFile *myMsDataFile, CommonParameters *commonParameters, std::vector<std::wstring> &nestedIds);

		protected:
			MetaMorpheusEngineResults *RunSpecific() override;
		};
	}
}